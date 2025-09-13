#include <phg.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*---------------------------------------------------------------------------
 * This program solves the following interface problem using XFEM:
 *	-\div (a \grad u) = f,		in \Omega^- \cup \Omega^+,
 *			u = g_D,	on \partial\Omega_D,
 *		 a\grad u = g_N\cdot n,	on \partial\Omega_N,
 *		      [u] = j_D,	on \Gamma,
 *	       [a\grad u] = j_N\cdot n,	on \Gamma,
 * Where
 * 	\Omega^- := \Omega \cap \{x: L(x)<0\},
 * 	\Omega^+ := \Omega \cap \{x: L(x)>0\},
 * 	\Gamma	 := \Omega \cap \{x: L(x)=0\},
 * 	a := COEFF1 for x\in\Omega^- and COEFF2 for x\in\Omega^+
 * L(x) is the level-set function defined by ls_func()
 *
 * The analytic test solution is given by func_u1() in \Omega^- and func_u2()
 * in \Omega^+.
 *
 * $Id: interface.c,v 1.292 2022/09/20 11:57:03 zlb Exp $
 *---------------------------------------------------------------------------*/

/* Note: the following commands should produce exactly the same results:
 *	./ipdg -test_case=1 -rel -refine0 1
 *	./interface -xfem_dtol=0 -xfem_vfrac=0 -coeff1=1 -r 2
 * 	./interface -coeff1=1 -coeff2=1 -ls_order=1 -c=0.5 -no_jump
 */

/* The coefficients */
static FLOAT coeff1 = 0.5, coeff2 = 2.0;

/* The IP parameters */
static FLOAT beta = 1.0, gamma0 = 10.0, gamma1 = 0.1;

static FLOAT xc = 0.5, yc = 0.5, zc = 0.5, c = 0.3;
static INT ls_order = 2;

/* extra quad. orders for interface elements */
static INT added_order = 3;

/* control parameters */
BOOLEAN no_jump = FALSE;
static INT sol_order = -1; /* order of analytic solution (< 0 ==> non poly.) */

static BOOLEAN dump_solver = FALSE;

/* analytic solution */
static void
func_u1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if (sol_order < 0)
	*value = Exp(x * y * z);
    else
	*value = 1.0 + Pow(x, sol_order) +
		 2.0 * Pow(y, sol_order) +
		 3.0 * Pow(z, sol_order);
}

static void
func_u2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if (sol_order < 0 && !no_jump)
	*value = Sin(x + y + z);
    else
	func_u1(x, y, z, value);
}

static void
func_jD(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT u2;

    func_u1(x, y, z, value);
    func_u2(x, y, z, &u2);
    *value -= u2;
}

static void
func_grad_u1(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    if (sol_order < 0) {
	values[0] = y * z * Exp(x * y * z);
	values[1] = x * z * Exp(x * y * z);
	values[2] = x * y * Exp(x * y * z);
    }
    else {
	values[0] = 1.0 * sol_order * Pow(x, sol_order - 1.0);
	values[1] = 2.0 * sol_order * Pow(y, sol_order - 1.0);
	values[2] = 3.0 * sol_order * Pow(z, sol_order - 1.0);
    }
}

static void
func_grad_u2(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    if (sol_order < 0 && !no_jump) {
	values[0] = Cos(x + y + z);
	values[1] = Cos(x + y + z);
	values[2] = Cos(x + y + z);
    }
    else {
	func_grad_u1(x, y, z, values);
    }
}

static void
func_jN(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    FLOAT un2[Dim];

    func_grad_u1(x, y, z, values);
    func_grad_u2(x, y, z, un2);
    values[0] = coeff1 * values[0] - coeff2 * un2[0];
    values[1] = coeff1 * values[1] - coeff2 * un2[1];
    values[2] = coeff1 * values[2] - coeff2 * un2[2];
}

/* right hand side */
static void
func_f1(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if (sol_order < 0) {
	*value = -(x * x * y * y + x * x * z * z + y * y * z * z)
						* Exp(x * y * z) * coeff1;
    }
    else {
	if (sol_order <= 1)
	    *value = -0. * coeff1;
	else
	    *value = -sol_order * (sol_order - 1.) * coeff1 *
		(1.0 * Pow(x, sol_order - 2.) +
		 2.0 * Pow(y, sol_order - 2.) +
		 3.0 * Pow(z, sol_order - 2.));
    }
}

static void
func_f2(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if (sol_order < 0 && !no_jump) {
	*value = 3. * Sin(x + y + z) * coeff2;
    }
    else {
	func_f1(x, y, z, value);
	*value *= coeff2 / coeff1;
    }
}

/* level set function */

static void
ls_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    assert(ls_order == 1 || ls_order == 2);
    if (ls_order == 2) {
	/* Sphere */
	x -= xc; y -= yc; z -= zc;
	*value = x * x + y * y + z * z - c * c;
    }
    else {
	/* plane */
	*value = x * xc + y * yc + z * zc - c;
    }
}

static void
ls_grad_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *grad)
{
    assert(ls_order == 1 || ls_order == 2);
    if (ls_order == 2) {
	/* Sphere */
	x -= xc; y -= yc; z -= zc;
	grad[0] = x + x; grad[1] = y + y; grad[2] = z + z;
    }
    else {
	/* plane */
	grad[0] = xc; grad[1] = yc; grad[2] = zc;
    }
}

#if 0
/* FIXME: linking against the static library libphg.a, or enabling the line
 * below, makes the code noticebly faster (>2x) */
# include "../src/quad-cache.c"
#endif

/* function numbers in the QCACHE objects */
static int Q_f, Q_gD, Q_gN, Q_jD, Q_jN;

static void
do_face(SOLVER *solver, XFEM_INFO *xi,
	QCACHE *qc,  ELEMENT *e,  int face,  int dof_id,  FLOAT coe,
	QCACHE *qc1, ELEMENT *e1, int face1, int dof_id1, FLOAT coe1)
/* This function computes integrations on a surface area, which can be:
 *
 *   1) an element face, in this case:
 * 	e != e1 && dof_id == dof_id1, e and e1 are both non interface elements
 *
 *   2) part of an element face, in this case:
 *   	e != e1 && dof_id == dof_id1, e and e1 are both interface elements
 *
 *   3) \Gamma \cap e, in this case:
 *	e == e1, dof_id = 0, dof_id1 = 1, and e is an interface element
 *
 *   4) an element face which is entirely contained in \Gamma, in this case:
 *	e1 != NULL, e1 != e, mark(e) < 0, mark(e1) > 0, dof_id = 0, dof_id1 = 1
 *
 * For cases 1) and 2), the normal vector stored in qc is the outward normal
 * vector of e, and e1 can be NULL (boundary face).
 *
 * For the cases 3) and 4) the normal vectors stored in qc points from
 * \Omega^- to \Omega^+.
 */
{
    GRID *g;
    DOF *u_h, *u1_h;
    int n, n1, p;   /* p = polynomial order */
    INT I, J;
    int i, j;
    FLOAT val;
    FLOAT G0, G1, h, a; /* G0 = coeff*gamma0*p^2/h, G1 = coeff*gamma1*h/p^2 */
    BTYPE bdry_flag;

    if (phgQCGetNP(qc) == 0)
	return;

    u_h = qc->fe;
    u1_h = qc1->fe;
    g = u_h->g;

    n = qc->qd->n_bas(qc, e->index);
    p = DofTypeOrder(u_h, e);

    h = e == e1 ? phgGeomGetDiameter(g, e) :
		  phgGeomGetFaceDiameter(g, e, face);

    if (e1 == NULL) {
	/* boundary face */
	bdry_flag = GetFaceBTYPE(g, e, face);
	if (bdry_flag != NEUMANN)
	    bdry_flag = DIRICHLET;
	n1 = 0;
	
	G0 = coe * gamma0 * p * p / h;
	G1 = coe * gamma1 * h / (p * (FLOAT)p);
	/* RHS */
	for (i = 0; i < n; i++) {
	    I = phgSolverMapE2G(solver, dof_id, e, i);

	    if (bdry_flag == DIRICHLET) {	/* Dirichlet boundary */
		/* -\beta\int_\Gamma_D g_D (A\grad v).n */
		a = phgQCIntegrateFace(
				qc, e->index, face, Q_gD,   PROJ_NONE, 0,
				qc, e->index, face, Q_GRAD, PROJ_DOT,  i)
				* coe;
		val = -beta * a;

		/* G0 \int_\Gamma_D g_D v */
		a = phgQCIntegrateFace(
				qc, e->index, face, Q_gD, PROJ_NONE, 0,
				qc, e->index, face, Q_BAS, PROJ_NONE, i);
		val += G0 * a;
	    }
	    else {				/* Neumann boundary */
		/* \int_\Gamma_N A g_N v */
		val = phgQCIntegrateFace(
				qc, e->index, face, Q_gN,  PROJ_DOT,  0,
				qc, e->index, face, Q_BAS, PROJ_NONE, i)
				* coe;
		/* \int_\Gamma_N G1 A g_N (A\grad v).n */
		a = phgQCIntegrateFace(
				qc, e->index, face, Q_gN,   PROJ_DOT, 0,
				qc, e->index, face, Q_GRAD, PROJ_DOT, i)
				* coe * coe;
		val += G1 * a;
	    }

	    phgSolverAddGlobalRHSEntry(solver, I, val);
	}
    }
    else {
	bdry_flag = INTERIOR;
	n1 = qc->qd->n_bas(qc1, e1->index);
	i = DofTypeOrder(u1_h, e1);
	if (p < i)
	    p = i;
	G0 = 0.5 * (coe + coe1) * gamma0 * p * (FLOAT)p / h;
	G1 = 0.5 * (coe + coe1) * gamma1 * h / (p * (FLOAT)p);
    }

    /* The following macros depend on the assertion */

#define Sel(i, o, o1)	(i < n ? o : o1)
#define Qc(i)	Sel(i, qc, qc1)
#define Ele(i)	Sel(i, e, e1)
#define Fac(i)	Sel(i, face, face1)
#define Bas(i)	Sel(i, i, i - n)
#define Dof(i)	Sel(i, dof_id, dof_id1)
#define Coe(i)	Sel(i, coe, coe1)
#define Eid(i)	Ele(i)->index
#define Quad(fid1, proj1, i1, fid2, proj2, i2) phgQCIntegrateFace( \
		Qc(i1), Eid(i1), Fac(i1), fid1, proj1, Bas(i1), \
		Qc(i2), Eid(i2), Fac(i2), fid2, proj2, Bas(i2))

    /* loop on {basis funcs in e} \cup {basis funcs in e1} */
    for (i = 0; i < n + n1; i++) {
	I = phgSolverMapE2G(solver, Dof(i), Ele(i), Bas(i));
	/* loop on {basis funcs in e} \cup {basis funcs in e1} */
	for (j = 0; j < n + n1; j++) {
	    J = phgSolverMapE2G(solver, Dof(j), Ele(j), Bas(j));

	    /* skip jumps for interior face and continuous element */
	    if (DofFESpace(u_h) == FE_H1 && e != e1 && e1 != NULL &&
		xi->info[e->index].mark * xi->info[e1->index].mark >= 0) {
		assert(bdry_flag == INTERIOR && dof_id == dof_id1);
		continue;
	    }

	    /*-----------------------------------------------------------------
	     * Note if the normal vector is reversed, the sign of [.] changes,
	     * while the sign of {.} is unaffected.
	     *----------------------------------------------------------------*/

	    val = 0.0;

	    if (bdry_flag != NEUMANN) {
		/* -\int {A\grad u}.n [v] (i<=>v, j<=>u, n<=>e) */
		a = Quad(Q_GRAD,  PROJ_DOT, j, Q_BAS, PROJ_NONE, i) * Coe(j);
		if (bdry_flag == INTERIOR)
		    a *= (i < n ? 0.5 : -0.5);
		val = -a;

		/* -\int \beta [u]{A\grad v}.n (note: i<=>v, j<=>u, n<=>e) */
		a = Quad(Q_BAS, PROJ_NONE, j, Q_GRAD,  PROJ_DOT, i) * Coe(i);
		if (bdry_flag == INTERIOR)
		    a *= (j < n ? 0.5 : -0.5);
		val += -beta * a;

		/* \int G0 [u][v] (i<=>v, j<=>u, n<=>e) */
		a = Quad(Q_BAS, PROJ_NONE, j, Q_BAS, PROJ_NONE, i);
		if (bdry_flag == INTERIOR && (i < n) != (j < n))
		    a = -a;
		val += G0 * a;
	    }

	    if (bdry_flag != DIRICHLET) {
		/* \int G1 [A\grad u].n [A\grad v].n */
		a = Quad(Q_GRAD,  PROJ_DOT, j, Q_GRAD,  PROJ_DOT, i)
			* Coe(j) * Coe(i);
		if (bdry_flag == INTERIOR && (i < n) != (j < n))
		    a = -a;
		val += G1 * a;
	    }

	    phgSolverAddGlobalMatrixEntry(solver, I, J, val);
	}

	if (e1 == NULL || dof_id == dof_id1)
	    continue;		/* non interface face */

	/* this face is part of the interface, apply jump conditions */
	assert(dof_id == 0 && dof_id1 == 1);
	assert((e == e1 && xi->info[e->index].mark == 0) ||
	       (xi->info[e->index].mark < 0 && xi->info[e1->index].mark > 0));

	/* \int [A*gN.n] {v} */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jN,  PROJ_DOT,  0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
	val = a * 0.5;

	/* G1 \int [A*gN.n] [(A\grad v).n] */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jN,   PROJ_DOT, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
	val += Coe(i) * G1 * (i < n ? a : -a);

	/* -beta \int [gD] {(A\grad v).n} (func[1] := gD) */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jD,   PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_GRAD, PROJ_DOT, Bas(i));
	val += -beta * 0.5 * Coe(i) * a;

	/* G0 \int [gD] [v] */
	a = phgQCIntegrateFace(Qc(i), Eid(i), Fac(i), Q_jD,  PROJ_NONE, 0,
			       Qc(i), Eid(i), Fac(i), Q_BAS, PROJ_NONE, Bas(i));
	val += G0 * (i < n ? a : -a);

	phgSolverAddGlobalRHSEntry(solver, I, val);
    }
#undef Sel
#undef Qc
#undef Ele
#undef Fac
#undef Bas
#undef Dof
#undef Coe
#undef Quad
}

static void
build_linear_system(XFEM_INFO *xi, SOLVER *solver, DOF *u1_h, DOF *u2_h,
		    DOF *f1_h, DOF *f2_h)
{
    GRID *g = xi->ls->g;
    QCACHE *qm, *qp;
    ELEMENT *e;
    
    assert(u1_h->dim == 1 && u2_h->dim == 1);

    qm = phgXFEMNewQC(xi, u1_h, 0);
    Q_f = phgQCAddFEFunction(qm, f1_h);			/* f */
    Q_gD = phgQCAddXYZFunction(qm, func_u1, 1);		/* gD */
    Q_gN = phgQCAddXYZFunction(qm, func_grad_u1, 3);	/* gN */
    Q_jD = phgQCAddXYZFunction(qm, func_jD, 1);		/* jD */
    Q_jN = phgQCAddXYZFunction(qm, func_jN, 3);		/* jN */

    qp = phgXFEMNewQC(xi, u2_h, 1);
    /* The following must be in exactly the same order as above */
    Q_f = phgQCAddFEFunction(qp, f2_h);			/* f  */
    Q_gD = phgQCAddXYZFunction(qp, func_u2, 1);		/* gD */
    Q_gN = phgQCAddXYZFunction(qp, func_grad_u2, 3);	/* gN */
    Q_jD = phgQCAddXYZFunction(qp, func_jD, 1);		/* jD */
    Q_jN = phgQCAddXYZFunction(qp, func_jN, 3);		/* jN */

#if USE_OMP
#pragma omp parallel for private(e)
#endif	/* USE_OMP */
    ForAllElementsBegin(g, e) {
	FLOAT *rule_m = NULL, *rule_0 = NULL, *rule_p = NULL;
	int N = qm->qd->n_bas(qm, e->index);
	INT I, J, eno = e->index;
	int i, j, k, face, order;
	FLOAT val;

	/* construct quadrature rules for the subdomains and the interface */
	i = DofTypeOrder(u1_h, e);
	j = DofTypeOrder(u2_h, e);
	order = 2 * (i > j ? i : j);
	if (xi->info[e->index].mark == 0)
	    phgQuadInterface(xi->ls, xi->ls_grad, e, order + added_order,
				&rule_m, &rule_0, &rule_p);
	else if (xi->info[e->index].mark < 0)
	    phgQuadInterface(NULL, xi->ls_grad, e, order /*- 2*/,
			     &rule_m, NULL, NULL);
	else if (xi->info[e->index].mark > 0)
	    phgQuadInterface(NULL, xi->ls_grad, e, order/* - 2*/,
			     NULL, NULL, &rule_p);

	/*==================== volume integrals ====================*/

	phgQCSetRule(qm, rule_m, -1.);
	phgQCSetRule(qp, rule_p, -1.);

	for (k = 0; k < 2; k++) {
	    QCACHE *q;
	    if ((k == 0 && xi->info[e->index].mark > 0) ||
		(k == 1 && xi->info[e->index].mark < 0))
		continue;
	    q = k == 0 ? qm : qp;
	    for (i = 0; i < N; i++) {
		I = phgSolverMapE2G(solver, k, e, i);
		for (j = 0; j <= i; j++) {
		    /* \int_T A\grad u1 . \grad v */
		    J = phgSolverMapE2G(solver, k, e, j);
		    val = phgQCIntegrate(q, eno, Q_GRAD, j, q, eno, Q_GRAD, i)
			  * (k == 0 ? coeff1 : coeff2);
		    phgSolverAddGlobalMatrixEntry(solver, I, J, val); 
		    if (i != j)
			phgSolverAddGlobalMatrixEntry(solver, J, I, val); 
		}

		/* \int_T f1 v */
		val = phgQCIntegrate(q, eno, Q_BAS, i, q, eno, Q_f, 0);
		phgSolverAddGlobalRHSEntry(solver, I, val);
	    }
	}
	phgFree(rule_m);
	phgFree(rule_p);
	rule_m = rule_p = NULL;

	/*==================== face integrals ====================*/

	if (xi->info[e->index].mark == 0) {	/* the interface */
	    phgQCSetRule(qm, rule_0, -1.);
	    phgQCSetRule(qp, rule_0, -1.);
	    do_face(solver, xi, qm, e, -1, /*dof_id*/0,  coeff1,
				qp, e, -1, /*dof_id1*/1, coeff2);
	}
	phgFree(rule_0);
	rule_0 = NULL;

	/* the faces of the element */
	for (face = 0; face < NFace; face++) {
	    FLOAT nv[Dim];
	    ELEMENT *e1 = phgGetNeighbour(g, e, face);
	    int face1 = -1;
	    if (e1 != NULL) {
		/* a face is only processed by the smaller of the two
		 * neighbouring elements, and by the one with smaller
		 * global index if the two elements are of the same size,
		 * to avoid double counting or redundant computation */
		if (e->generation < e1->generation)
		    continue;
		if (e->generation == e1->generation &&
		    GlobalElement(g, e->index) > GlobalElement(g, e1->index))
		    continue; /* process each interior face just once */
		face1 = phgOppositeFace(g, e, face, e1);
	    }
	    phgGeomGetFaceOutNormal(g, e, face, nv);
	    if (xi->info[e->index].mark == 0)
		/* the face might be cut into 2 parts by the interface */
		phgQuadInterfaceFace(xi->ls, xi->ls_grad, e, face,
				     order + added_order,
				     &rule_m, NULL, &rule_p);
	    else if (xi->info[e->index].mark < 0)
		phgQuadInterfaceFace(NULL, xi->ls_grad, e, face, order,
				     &rule_m, NULL, NULL);
	    else
		phgQuadInterfaceFace(NULL, xi->ls_grad, e, face, order,
				     NULL, NULL, &rule_p);

	    /* Note: e1 is either NULL or valid (since gen(e) >= gen(e1)) */
	    if (e1 != NULL &&
		xi->info[e->index].mark * xi->info[e1->index].mark < 0) {
		/* The face is flat and is part of the interface
		 * (e and e1 in different subdomains) */
		if (xi->info[e->index].mark < 0) {
		    phgQCSetRule(qm, rule_m, -1.);
		    phgQCSetConstantNormal(qm, nv);
		    phgQCSetRule(qp, rule_m, -1.);
		    phgQCSetConstantNormal(qp, nv);
		    do_face(solver, xi, qm, e,  face,  0, coeff1,
					qp, e1, face1, 1, coeff2);
		}
		else {
		    nv[0] = -nv[0]; nv[1] = -nv[1]; nv[2] = -nv[2];
		    phgQCSetRule(qp, rule_p, -1.);
		    phgQCSetConstantNormal(qp, nv);
		    phgQCSetRule(qm, rule_p, -1.);
		    phgQCSetConstantNormal(qm, nv);
		    do_face(solver, xi, qm, e1, face1, 0, coeff1,
					qp, e,  face,  1, coeff2);
		}
	    }
	    else {
		/* e and e1 in a same subdomain */
		if (xi->info[e->index].mark <= 0) {
		    phgQCSetRule(qm, rule_m, -1.);
		    phgQCSetConstantNormal(qm, nv);
		    do_face(solver, xi, qm, e,  face,  0, coeff1,
					qm, e1, face1, 0, coeff1);
		}
		if (xi->info[e->index].mark >= 0) {
		    phgQCSetRule(qp, rule_p, -1.);
		    phgQCSetConstantNormal(qp, nv);
		    do_face(solver, xi, qp, e,  face,  1, coeff2,
					qp, e1, face1, 1, coeff2);
		}
	    }
	    phgFree(rule_m);
	    phgFree(rule_p);
	    rule_m = rule_p = NULL;
	}
    } ForAllElementsEnd

    phgQCFree(&qm);
    phgQCFree(&qp);

    /* set the diagonal entry of empty rows to 1 */
    phgXFEMProcessEmptyRows(solver);

    if (dump_solver) {
	phgSolverAssemble(solver);
	phgSolverUpdateRHS(solver);
	phgSolverDumpMATLAB_(solver, "A", "b", NULL, TRUE);
    }
}

int
main(int argc, char *argv[])
{
    char *fn = "../p4est/cube.mesh";
    int i, level, total_level = 0;
    INT refine = 0;
#ifdef PHG_TO_P4EST
    INT refine0 = 1, refine_step = 1;
#else	/* defined(PHG_TO_P4EST) */
    INT refine0 = 0, refine_step = 3;
#endif	/* defined(PHG_TO_P4EST) */
    GRID *g;
    XFEM_INFO *xi;
    DOF *u1_h, *u2_h, *f1_h, *f2_h, *error1, *gerror1, *error2, *gerror2, *gu_h;
    DOF *u1_old, *u2_old, *ls, *ls_grad;
    SOLVER *solver;
    size_t mem_peak;
    double t;
    FLOAT L2norm, H1norm, L2err, H1err, d;
    BOOLEAN vtk = FALSE;
    ELEMENT *e;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", &fn);
    phgOptionsRegisterInt("-refine0", "Initial refinement levels", &refine0);
    phgOptionsRegisterInt("-refine", "Repeated refinement levels", &refine);
    phgOptionsRegisterInt("-refine_step", "Refinement step", &refine_step);
    phgOptionsRegisterFloat("-c",  "Radius or offset of the interface", &c);
    phgOptionsRegisterFloat("-xc", "Center or normal of the interface: x", &xc);
    phgOptionsRegisterFloat("-yc", "Center or normal of the interface: y", &yc);
    phgOptionsRegisterFloat("-zc", "Center or normal of the interface: z", &zc);
    phgOptionsRegisterFloat("-coeff1", "The coefficient coeff1", &coeff1);
    phgOptionsRegisterFloat("-coeff2", "The coefficient coeff2", &coeff2);
    phgOptionsRegisterFloat("-beta", "The parameter beta", &beta);
    phgOptionsRegisterFloat("-gamma0", "The parameter gamma0", &gamma0);
    phgOptionsRegisterFloat("-gamma1", "The parameter gamma1", &gamma1);
    phgOptionsRegisterInt("-added_order", "Extra quadrature orders for "
					"interface elements", &added_order);
    phgOptionsRegisterNoArg("-vtk", "Create \"interface.vtk\"", &vtk);

    phgOptionsRegisterNoArg("-dump_solver", "Output matrix/RHS as .m files",
					&dump_solver);

    phgOptionsRegisterInt("-ls_order", "Polynomial order of the levelset"
			    "function (1: plane, 2: sphere)", &ls_order);
    phgOptionsRegisterInt("-sol_order", "Polynomial order of the solution "
			    "(<0 => no polynomial)", &sol_order);
    phgOptionsRegisterNoArg("-no_jump", "No jump in the solution across the "
			    "interface", &no_jump);

    phgOptionsPreset("-dof_type DG1 -solver gmres");
#if USE_HYPRE
    /* BoomerAMG seems to work pretty well */
    phgOptionsPreset("-solver hypre -hypre_solver gmres -hypre_pc boomeramg");
#endif

    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);

    /* print data for the sphere or the plane to help debugging with ParaView */
    assert(ls_order == 1 || ls_order == 2);
    if (ls_order == 1) {
	/* find a parallelogram on the plane enclosing the mesh for ParaView */
#if DEBUG && defined(PHG_TO_P4EST)
	/* local coordinate system on the plan: (origin, v0, v1) */
	double d, bbox[2][2] = {{0.,0.},{0.,0.}}, origin[Dim], v0[Dim], v1[Dim];
	double nv[] = {xc, yc, zc};
	int flag = 0;
	/* normalize nv */
	d = 1.0 / sqrt(nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]);
	nv[0] *= d;
	nv[1] *= d;
	nv[2] *= d;
	/* origin := projection of the center */
	origin[0] = (g->bbox[0][0] + g->bbox[1][0]) * 0.5;
	origin[1] = (g->bbox[0][1] + g->bbox[1][1]) * 0.5;
	origin[2] = (g->bbox[0][2] + g->bbox[1][2]) * 0.5;
	/* compute d such that (origin + d*nv).(xc,yc,zc) == c */
	d *= c - origin[0] * xc - origin[1] * yc - origin[2] * zc;
	origin[0] += d * nv[0];
	origin[1] += d * nv[1];
	origin[2] += d * nv[2];
	/* i := the smallest component of (xc,yc,zc) */
	d = Fabs(xc); i = 0;
	if (d > Fabs(yc)) {d = Fabs(yc); i = 1;}
	if (d > Fabs(zc)) {d = Fabs(zc); i = 2;}
	switch (i) {
	    case 0: v0[0] = 0.; v0[1] = nv[2]; v0[2] = -nv[1]; break;
	    case 1: v0[1] = 0.; v0[0] = nv[2]; v0[2] = -nv[0]; break;
	    case 2: v0[2] = 0.; v0[0] = nv[1]; v0[1] = -nv[0]; break;
	}
	d = 1.0 / sqrt(v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]);
	v0[0] *= d;
	v0[1] *= d;
	v0[2] *= d;
	/* v1 := nv cross v0 */
	v1[0] = nv[1] * v0[2] - nv[2] * v0[1];
	v1[1] = nv[2] * v0[0] - nv[0] * v0[2];
	v1[2] = nv[0] * v0[1] - nv[1] * v0[0];
	for (i = 0; i < NEdge; i++) {
	    /* find the cut point with the i-th edge */
	    int i0 = GetEdgeVertex(i, 0), i1 = GetEdgeVertex(i, 1);
	    double c0[] = {g->bbox[i0 & 1][0],
			   g->bbox[(i0>>1) & 1][1],
			   g->bbox[(i0>>2) & 1][2]},
		   c1[] = {g->bbox[i1 & 1][0] - c0[0],
			   g->bbox[(i1>>1) & 1][1] - c0[1],
			   g->bbox[(i1>>2) & 1][2] - c0[2]};
	    /* solve c0[].{xc,yc,zc} + t*c1[].{xc,yc,zc} == c for t */
	    double a = c1[0] * xc + c1[1] * yc + c1[2] * zc,
		   b = c - c0[0] * xc - c0[1] * yc - c0[2] * zc;
	    if (a == 0.0)
		continue;
	    if ((d = b / a) < -FLOAT_EPSILON || d > 1 + FLOAT_EPSILON)
		continue;
	    c0[0] += d * c1[0];
	    c0[1] += d * c1[1];
	    c0[2] += d * c1[2];
	    c0[0] -= origin[0];
	    c0[1] -= origin[1];
	    c0[2] -= origin[2];
	    /* (a,b) := coordinates w.r.t. v0 and v1 */
	    a = c0[0] * v0[0] + c0[1] * v0[1] + c0[2] * v0[2];
	    b = c0[0] * v1[0] + c0[1] * v1[1] + c0[2] * v1[2];
	    if (flag == 0) {
		flag = 1;
		bbox[0][0] = bbox[1][0] = a;
		bbox[0][1] = bbox[1][1] = b;
		continue;
	    }
	    if (a < bbox[0][0])
		bbox[0][0] = a;
	    else if (a > bbox[1][0])
		bbox[1][0] = a;
	    if (b < bbox[0][1])
		bbox[0][1] = b;
	    else if (b > bbox[1][1])
		bbox[1][1] = b;
	}
	/* enlarge the bbox by 10% */
	d = (bbox[1][0] - bbox[0][0]) * 0.1;
	bbox[0][0] -= d;
	bbox[1][0] += d;
	d = (bbox[1][1] - bbox[0][1]) * 0.1;
	bbox[0][1] -= d;
	bbox[1][1] += d;
	phgPrintf("The plane:\n    origin: (%g %g %g)\n    "
		  "point 1: (%g %g %g)\n    point 2: (%g %g %g)\n",
			/* origin: (bbox[0][0], bbox[0][1]) */
			origin[0] + bbox[0][0] * v0[0] + bbox[0][1] * v1[0],
			origin[1] + bbox[0][0] * v0[1] + bbox[0][1] * v1[1],
			origin[2] + bbox[0][0] * v0[2] + bbox[0][1] * v1[2],
			/* point1: (bbox[1][0], bbox[0][1]) */
			origin[0] + bbox[1][0] * v0[0] + bbox[0][1] * v1[0],
			origin[1] + bbox[1][0] * v0[1] + bbox[0][1] * v1[1],
			origin[2] + bbox[1][0] * v0[2] + bbox[0][1] * v1[2],
			/* point2: (bbox[0][0], bbox[1][1]) */
			origin[0] + bbox[0][0] * v0[0] + bbox[1][1] * v1[0],
			origin[1] + bbox[0][0] * v0[1] + bbox[1][1] * v1[1],
			origin[2] + bbox[0][0] * v0[2] + bbox[1][1] * v1[2]);
#endif	/* DEBUG && defined(PHG_TO_P4EST) */
    }
    else {
	phgPrintf("Sphere: center=(%g %g %g), radius=%g\n",
				(double)xc, (double)yc, (double)zc, (double)c);
    }

    /* refine the mesh */
    phgPrintf("Initial refinement (refine0=%d): ", refine0);
    while (refine0 > 0) {
	level = refine0 > refine_step ? refine_step : refine0;
#if 1
	/* uniform refinement */
	phgRefineAllElements(g, level);
#else
	/* refine the north pole */
	ELEMENT *e;
	FLOAT pole[] = {xc, yc, zc + c};
	ForAllElements(g, e)
	    if (e->corners[0][0] <= pole[0] && e->corners[1][0] >= pole[0] &&
		e->corners[0][1] <= pole[1] && e->corners[1][1] >= pole[1] &&
		e->corners[0][2] <= pole[2] && e->corners[1][2] >= pole[2])
		e->mark = level;
	    else
		e->mark = 0;
	phgRefineMarkedElements(g);
#endif
    	phgBalanceGrid(g, 1.1, 1, NULL, 1.0);
	refine0 -= level;
	total_level += level;
    }
    phgPrintf("%"dFMT" element%s.\n", g->nelem_global,
					g->nelem_global > 1 ? "s" : "");

#if 1
    /* use analytic function */
    ls = phgDofNew(g, DOF_ANALYTIC, 1, "ls", ls_func);
    /*phgDofSetPolyOrder(ls, ls_order);*/
    ls_grad = phgDofNew(g, DOF_ANALYTIC, 3, "ls_grad", ls_grad_func);
    /*phgDofSetPolyOrder(ls_grad, ls_order <= 0 ? ls_order : ls_order - 1);*/
#else
    /* project the levelset function to a FE space */
    assert(ls_order >= 1);
    ls = phgDofNew(g, DOF_DGn[ls_order], 1, "ls", ls_func);
    ls_grad = phgDofNew(g, DOF_DGn[ls_order - 1], 3, "ls_grad", ls_grad_func);
    /*ls->userfunc = DofInterpolation;
    ls_grad->userfunc = DofInterpolation;*/
#endif

    u1_h = phgDofNew(g, DOF_DEFAULT, 1, "u1_h", DofInterpolation);
    u2_h = phgDofNew(g, DOF_DEFAULT, 1, "u2_h", DofInterpolation);
    
    t = phgGetTime(NULL);
    while (TRUE) {
	double h_min = 1e10, h_max = 0.0, nnz, nnz_d, nnz_o;	/* statistics */
	phgPrintf("\n********** Level %d, %d proc%s, %"
		  dFMT" elem%s, LIF %0.2lf, refine time: %0.4lg\n",
		  total_level, g->nprocs, g->nprocs > 1 ? "s" : "",
		  g->nleaf_global, g->nleaf_global > 1 ? "s" : "",
		  (double)g->lif, phgGetTime(NULL) - t);

	phgPrintf("Resolving/merging interface elements: ");
	t = phgGetTime(NULL);
	xi = phgXFEMInit(ls, ls_grad, ls_order,
					2 * u1_h->type->order + added_order);
	phgPrintf("%"dFMT" element%s, time %0.4lg\n",
			g->nelem_global, g->nelem_global > 1 ? "s" : "",
			phgGetTime(NULL) - t);
	phgPrintf("    Elements by region: -:%"dFMT", 0:%"dFMT", +:%"dFMT"\n",
			xi->ecnts[0], xi->ecnts[1], xi->ecnts[2]);

	phgSetupHalo(g, HALO_FACE);	/* HALO_FACE is needed for DG */

	if (beta == 1.0 && !phgOptionsIfUsed("-solver_symmetry"))
	    phgOptionsSetOptions("-solver_symmetry=spd");
	solver = phgSolverCreate(SOLVER_DEFAULT, u1_h, u2_h, NULL);
	solver->mat->handle_bdry_eqns = FALSE;
	phgPrintf("Building linear equations ...\n");
	t = phgGetTime(NULL);
	f1_h = phgDofNew(g, DOF_DEFAULT, 1, "f1_h",  func_f1);
	f2_h = phgDofNew(g, DOF_DEFAULT, 1, "f2_h",  func_f2);
	build_linear_system(xi, solver, u1_h, u2_h, f1_h, f2_h);
	phgDofFree(&f1_h);
	phgDofFree(&f2_h);
	/* get h_max and h_min */
	ForAllElements(g, e) {
	    FLOAT (*corners)[Dim] = phgGeomGetCorners(g, e, NULL);
	    for (i = 0; i < Dim; i++) {
		d = corners[1][i] - corners[0][i];
		if (h_max < d)
		    h_max = d;
		if (h_min > d)
		    h_min = d;
	    }
	}
	nnz_d = solver->mat->nnz_d;
	nnz_o = solver->mat->nnz_o;
	nnz = nnz_d + nnz_o;
#if USE_MPI
	if (g->nprocs > 1) {
	    double tmp0[4] = {-h_min, h_max, nnz_d, nnz_o}, tmp[4];
	    MPI_Reduce(tmp0, tmp, 4, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	    h_min = -tmp[0];
	    h_max = tmp[1];
	    nnz_d = tmp[2];
	    nnz_o = tmp[3];
	    tmp[0] = nnz;
	    MPI_Reduce(tmp, &nnz, 1, MPI_DOUBLE, MPI_SUM, 0, g->comm);
	}
#endif	/* USE_MPI */
	phgMemoryUsage(g, &mem_peak);
	phgPrintf("    %"dFMT" DOFs, h_max/h_min: %0.2le/%0.2le, "
		  "time: %0.2lg, mem: %0.2lfGB\n", solver->rhs->map->nglobal,
		  h_max, h_min, phgGetTime(NULL) - t,
		  (double)mem_peak / (1024. * 1024. * 1024.));
	phgPrintf("    # nonzeros: max inproc/offproc %0.0lf/%0.0lf, "
		  "total %0.0lf\n", nnz_d, nnz_o, nnz);
	phgPrintf("Solving linear equations ...\n");
	t = phgGetTime(NULL);
	u1_old = phgDofCopy(u1_h, NULL, NULL, "u1_old");    
	u2_old = phgDofCopy(u2_h, NULL, NULL, "u2_old");    

	phgXFEMSetInitialSolution(xi, u1_h, u2_h);
	phgSolverSolve(solver, TRUE, u1_h, u2_h, NULL);
	phgXFEMUpdateSolution(xi, u1_h, u2_h);

	phgPrintf("    solver %s; nits: %d; residual: %lg; time: %0.4lg\n",
	    solver->oem_solver->name, solver->nits, (double)solver->residual,
	    phgGetTime(NULL) - t);
	if (solver->cond > 0.0)
	    phgPrintf("    Condition number: %0.2le\n", (double)solver->cond);
	phgSolverDestroy(&solver);
    
	/* compute L2 and H1 errors */
	t = phgGetTime(NULL);

	/* error := u, projected to the FE space */
	error1 = phgDofNew(g, u1_h->type, 1, "error1", func_u1);
	error2 = phgDofNew(g, u2_h->type, 1, "error2", func_u2);
	L2norm = Sqrt(phgXFEMDot(xi, error1, error2, error1, error2));
    
	/* gerror := Grad(u), projected to the FE space */
	gerror1 = phgDofGradient(error1, NULL, NULL, NULL);
	gerror2 = phgDofGradient(error2, NULL, NULL, NULL);
	H1norm = Sqrt(phgXFEMDot(xi, gerror1, gerror2, gerror1, gerror2));
	H1norm += L2norm;

	d = Sqrt(phgXFEMDot(xi, u1_h, u2_h, u1_h, u2_h));
	if (L2norm < d)
	    L2norm = d;
	phgDofFree(&gerror1);
	phgDofFree(&gerror2);
	gerror1 = phgDofGradient(error1, NULL, NULL, NULL);
	gerror2 = phgDofGradient(error2, NULL, NULL, NULL);
	d = L2norm + Sqrt(phgXFEMDot(xi, gerror1, gerror2, gerror1, gerror2));
	if (H1norm < d)
	    H1norm = d;
    
	/* error := u_h - error */
	phgDofAXPY(-1.0, u1_h, &error1);
	phgDofAXPY(-1.0, u2_h, &error2);
	/* L2err := |error| */
	L2err = Sqrt(phgXFEMDot(xi, error1, error2, error1, error2));
    
	/* gerror := Grad(u_h) - gerror */
	gu_h = phgDofGradient(u1_h, NULL, NULL, NULL);
	phgDofAXPY(-1.0, gu_h, &gerror1);
	phgDofFree(&gu_h);
	gu_h = phgDofGradient(u2_h, NULL, NULL, NULL);
	phgDofAXPY(-1.0, gu_h, &gerror2);
	phgDofFree(&gu_h);
	/* H1err := L2err + |gerror| */
	H1err = Sqrt(phgXFEMDot(xi, gerror1, gerror2, gerror1, gerror2));
	H1err += L2err;

	phgMemoryUsage(g, &mem_peak);
	phgPrintf("L2err: %0.10le; H1err: %0.10le; "
		  "mem: %0.2lfGB; time: %0.2lg\n",
		    (double)L2err / (L2norm == 0. ? 1.0 : (double)L2norm),
		    (double)H1err / (H1norm == 0. ? 1.0 : (double)H1norm),
		    (double)mem_peak / (1024. * 1024. * 1024.),
		    phgGetTime(NULL) - t);
    
	phgDofAXPY(-1.0, u1_h, &u1_old);
	phgDofAXPY(-1.0, u2_h, &u2_old);
	phgPrintf("    |u_h-u_H| = %0.5le\n", (double)
			Sqrt(phgXFEMDot(xi, u1_old, u2_old, u1_old, u2_old)));

	if (vtk) {
	    char name[128];
	    ForAllElements(g, e)
		e->region_mark = xi->info[e->index].mark;
	    sprintf(name, "interface-%02d.vtk", total_level);
	    phgExportVTK(g, name, u1_h, u2_h, error1, error2, ls, NULL);
	    phgPrintf("\"%s\" created.\n", name);
	}

	phgXFEMFree(&xi);

	phgDofFree(&u1_old);
	phgDofFree(&u2_old);
	phgDofFree(&error1);
	phgDofFree(&gerror1);
	phgDofFree(&error2);
	phgDofFree(&gerror2);

	if (refine <= 0)
	    break;

	level = refine > refine_step ? refine_step : refine;
	phgRefineAllElements(g, level);
	phgBalanceGrid(g, 1.2, 0, NULL, 1.0);
	refine -= level;
	total_level += level;
    }

    phgDofFree(&u1_h);
    phgDofFree(&u2_h);
    phgDofFree(&ls);
    phgDofFree(&ls_grad);
    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
