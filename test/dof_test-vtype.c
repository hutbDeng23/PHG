#include "phg.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef PHG_TO_P4EST
# include "phg-to-p4est.h"
#endif

int dim = Dim;

static void
func(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    int dim0 = dim * DOF_DEFAULT->dim;

    if (dim0 == 1) {
	*values = 0.5 * (x * x + y * y + z * z);
	return;
    }
    *(values++) = 0.5 * x * x;
    *(values++) = 0.5 * y * y;
    if (dim0 >= 3) {
	*(values++) = 0.5 * z * z;
	if (dim0 > 3)
	    bzero(values, (dim0 - 3) * sizeof(*values));
    }
}

typedef enum {BasBas = 0, GradGrad = 1, CurlCurl = 2, DivDiv = 3} MTYPE_t;
static const char *mtype_name[] = {"BasBas", "GradGrad", "CurlCurl", "DivDiv"};
static int mtype_dim[] = {0, 0, Dim, Dim};	/* req. dimension, or <=0 */

static MAT *
get_matrix(DOF *u, MTYPE_t mt)
/* returns the mass matrix of the DOF "u" */
{
    MAP *map;
    MAT *A;
    GRID *g = u->g;
    ELEMENT *e = NULL;
    QCACHE *qc;
    int dorder = 0, fid = -1;

    phgDofSetDirichletBoundaryMask(u, 0);
    map = phgMapCreate(u, NULL);
    A = phgMapCreateMat(map, map);
    phgMapDestroy(&map);

    switch (mt) {
	case BasBas:	dorder =  0; fid = Q_BAS;  break;
	case GradGrad:	dorder = -2; fid = Q_GRAD; break;
	case CurlCurl:	dorder = -2; fid = Q_CURL; break;
	case DivDiv:	dorder = -2; fid = Q_DIV;  break;
    }

    qc = phgQCNew(QD_DEFAULT, u);
    ForAllElements(g, e) {
	int i, j, N = DofNBas(u, e);
	INT I, J, eno = e->index;
	FLOAT val, *rule;
	rule = phgQuadGetRule3D(g, e, DofTypeOrder(u,e) + dorder);
	phgQCSetRule(qc, rule, -1.);
	for (i = 0; i < N; i++) {
	    I = phgMapE2G(A->rmap, 0, e, i);
	    for (j = 0; j <= i; j++) {
		J = phgMapE2G(A->cmap, 0, e, j);
		val = phgQCIntegrate(qc, eno, fid, j, qc, eno, fid, i);
		phgMatAddGlobalEntry(A, I, J, val); 
		if (i != j)
		    phgMatAddGlobalEntry(A, J, I, val); 
	    }
	}
	phgFree(rule);
    }
    phgQCFree(&qc);
    phgMatAssemble(A);

    return A;
}

int
main(int argc, char **argv)
{
    GRID *g;
    DOF *u, *tmp;
    DOF_TYPE *vtype;
    char *fn = "../p4est/cube.mesh";
    INT refine0 = 0, refine = 0;
    int i, mt;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterInt("-dim", "Vector dimension (1 or 3)", &dim);
    phgOptionsRegisterInt("-refine0", "Initial refinement levels", &refine0);
    phgOptionsRegisterInt("-refine", "Repeated refinement levels", &refine);

    //phgOptionsPreset("-dof_type P2 -refine 0 -mesh_file tetra.dat");
    phgInit(&argc, &argv);

    vtype = phgDofTypeVector(DOF_DEFAULT, dim, NULL);

    g = phgNewGrid(-1);
    phgImport(g, fn, TRUE);

    for (i = 0; i < refine0; i++) {
	phgRefineAllElements(g, 1);
	phgBalanceGrid(g, 1.2, 0, NULL, 1.0);
	phgPrintf("Initial refinement: %"dFMT" elements.\n", g->nelem_global);
    }
    

    u = phgDofNew(g, vtype, 1, "u", DofInterpolation);
    phgDofSetDataByFunction(u, func);
    //phgDofDump(u);
    tmp = phgDofGradient(u, NULL, NULL, "grad");
    phgDofSetFunction(tmp, DofInterpolation);
#if 0
    ELEMENT *e;
    assert(dim == Dim);
    ForAllElements(g, e) {
# ifdef PHG_TO_P4EST
	FLOAT xref[] = {0.5,0.5,0.5};
# else	/* defined(PHG_TO_P4EST) */
	FLOAT xref[] = {0.25,0.25,0.25,0.25};
# endif	/* defined(PHG_TO_P4EST) */
	FLOAT x[Dim], val[Dim * Dim];
	phgGeomLambda2XYZ(g, e, xref, x, x + 1, x + 2);
	phgDofEval(u, e, xref, val);
	phgInfo(-1, "u(%g,%g,%g) = (%g,%g,%g)\n",
			x[0], x[1], x[2], val[0], val[1], val[2]);
	phgDofEval(tmp, e, xref, val);
	phgInfo(-1, "grad(%g,%g,%g) = ((%g,%g,%g),(%g,%g,%g),(%g,%g,%g))\n",
			x[0], x[1], x[2], val[0], val[1], val[2],
					  val[3], val[4], val[5],
					  val[6], val[7], val[8]);
    }
#endif

    i = 0;
    while (i < refine) {
	i += 1;
	phgRefineAllElements(g, 1);
	phgBalanceGrid(g, 1.2, 0, NULL, 1.0);
	phgPrintf("Repeated refinement: %"dFMT" elements.\n", g->nelem_global);
    }
    phgPrintf("|u|_2 = %g, |grad(u)|_2 = %g\n",
		(double)phgDofNormL2(u), (double)phgDofNormL2(tmp));

    phgDofFree(&tmp);

    tmp = phgDofNew(g, vtype->base_type, 1, "u0", DofNoData);
    for (mt = 0; mt < sizeof(mtype_name) / sizeof(mtype_name[0]); mt++) {
	MAT *A, *a;
	INT j;
	phgPrintf("%s matrix: ", mtype_name[mt]);
	if (mtype_dim[mt] > 0 && vtype->dim != mtype_dim[mt]) {
	    phgPrintf("skipped, wrong dimension (%d)\n", vtype->dim);
	    continue;
	}
	A = get_matrix(u, mt);
	phgPrintf("norm = %g, ", (double)phgMatNormInfty(A));
	if (mtype_dim[mt] > 0 && DOF_DEFAULT->dim != mtype_dim[mt]) {
	    phgPrintf("not verified, wrong base dimension (%d)\n",
							DOF_DEFAULT->dim);
	    phgMatDestroy(&A);
	    continue;
	}
	/* check that "A" == P("a"), P being the 1D->3D expansion operator */
	a = get_matrix(tmp, mt);	/* matrix for the base type */
	for (j = 0; j < a->rmap->nlocal; j++) {
	    const MAT_ROW *row = phgMatGetRow(a, j);
	    for (i = 0; i < row->ncols; i++) {
		/* Note: column indices in row are global */
		INT k = row->cols[i];
		FLOAT d0 = row->data[i];
		int ii;
		for (ii = 0; ii < dim; ii++) {
		    FLOAT d = phgMatGetLGEntry(A, j * dim + ii, k * dim + ii);
		    if (Fabs(d0 - d) > 100.0 * FLOAT_EPSILON * Fabs(d0)) {
			j += a->rmap->partition[g->rank]; /* -> global index */
			phgError(1, "check failed: A(%d,%d)=%g, a(%d,%d)=%g\n",
				    j * dim + ii, k * dim + ii, d, j, k, d0); 
		    }
		}
	    }
	}
	phgPrintf("verified.\n");
	phgMatDestroy(&a);
	phgMatDestroy(&A);
    }
    phgDofFree(&tmp);

    phgDofFree(&u);
    phgDofTypeFree(&vtype);
    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
