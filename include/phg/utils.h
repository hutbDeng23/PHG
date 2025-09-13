/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: utils.h,v 1.47 2022/09/19 01:42:11 zlb Exp $ */

#ifndef PHG_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __WIN__
/* check for Windows specific macros */
# if defined(__WIN32) || defined(__WIN32__) || defined(_WIN32) || \
     defined(__WIN64) || defined(__WIN64__) || defined(_WIN64) || \
     defined(__WINNT) || defined(__WINNT__)
#  define __WIN__
# endif
#endif	/* !defined(__WIN__) */

#ifdef __WIN__
# define bzero(s, n)	memset(s, 0, n)
# define srand48	srand
# define drand48()	((double)rand()/(double)RAND_MAX)
/*# define sleep(i)	usleep((i) * 1000000)*/
#endif	/* defined(__WIN__) */

#if USE_MPI
# ifndef MPIAPI
#  define MPIAPI	/* This macro is needed for MicroSoft MSMPI */
# endif	/* defined(MPIAPI) */
#endif	/* USE_MPI */

/* Macro for avoiding 'unused variable ...' warnings */
#define Unused(v)	(void)(v)

/* Name of the current function */
#if defined __STDC_VERSION__
# if __STDC_VERSION__ < 199901L
#  if __GNUC__ >= 2
#   define __func__ __FUNCTION__
#  else
#   define __func__ "none"
#  endif
# endif
#endif

/* #define phgDebug(level, fmt, ...) phgInfo(level, fmt, ##__VA_ARGS__) */
#if DEBUG
# define phgDebug(s) phgInfo s
#else
# define phgDebug(s)
#endif

/* The macro CB_ARGS is used to control whether 'g' should be
 * in the argument list of the callback functions. */
#if 0
/* carry g in the argument list of callback functions */
# define CB_ARGS(e)	(GRID *g, ELEMENT *e)
# define CB_ARGS0(g, e)	(g, e)
#else
/* don't carry g in the argument list of callback functions */
# define CB_ARGS(e)	(ELEMENT *e)
# define CB_ARGS0(g, e)	(e)
#endif

#define IsLeaf(e) ((e)->children[0] == NULL && (e)->children[1] == NULL)
#define HasLocalNeighbour(e, v) \
    ((((ELEMENT *)(e))->bound_type[v] & (INTERIOR | REMOTE)) == INTERIOR)
#define HasRemoteNeighbour(e, v) \
    ((((ELEMENT *)(e))->bound_type[v] & (INTERIOR | REMOTE)) == \
	(INTERIOR | REMOTE))

#define GetNeighbour(e, v) \
    ((ELEMENT *)(HasLocalNeighbour(e, v) ? (e)->neighbours[v] : NULL))

#define GetRNeighbour(g, e, v) \
    ((HasRemoteNeighbour(e, v) && (g)->neighbours.count > 0) \
	? (g)->neighbours.list + (size_t)((e)->neighbours[v]) : \
	  (RNEIGHBOUR *)NULL)

/**< Magic numbers */
#define MagicNumber(a, b, c, d)	(((((((INT)d) << 8) | c) << 8 | b) << 8) | a)
#define MAGIC_DOF 	MagicNumber('D', 'O', 'F', ' ')
#define MAGIC_MAT 	MagicNumber('M', 'A', 'T', ' ')
#define MAGIC_VEC 	MagicNumber('V', 'E', 'C', ' ')

#define MagicCheck(m, var) {						\
    if ((var) != NULL && (var)->magic != MAGIC ## _ ## m) {		\
	int magic = MAGIC ## _ ## m;					\
	phgError(1, "%s (%s:%d): invalid or uninitialized %c%c%c object.\n", \
		__func__, __FILE__, __LINE__,				\
		magic & 255, (magic >> 8) & 255, (magic >> 16) & 255);	\
    }									\
}

/* Macro for warning an untested code segment */
#define WARN_UNTESTED	{						\
  static BOOLEAN warned = FALSE;					\
  if (!warned) {							\
    warned = TRUE;							\
    phgWarning("%s:%d: untested code segment!\n", __FILE__, __LINE__);	\
  }									\
}

#if USE_OMP
# if DEBUG
#  define CheckThread							\
    if (phgThreadId > 0)					\
	phgError(1, "The function %s (%s:%d) is not thread aware.\n",	\
			__func__, __FILE__, __LINE__);
# else	/* DEBUG */
#  define CheckThread
# endif	/* DEBUG */
#else	/* USE_OMP */
# define CheckThread
#endif	/* USE_OMP */

/**< Protocol for user functions mapping bdry types.
 *
 * The function should return a combination of PHG's bdry type, or -1 if
 * bc_type is unknown */
typedef int (*BDRY_MAP_FUNC)(int bc_type);

/* global variables */
extern const char *phgExeFn;
extern const char *phgCurrentFunctionName_;

/* Note: obj/build.sh defines DONT_WANT_phgThreadId as a work around for the
 * following error with certain compilers:
 * 	/usr/bin/ld: phgThreadId: TLS definition in utils.o section .tbss mismatches non-TLS reference in utils.o.
 */
#ifdef DONT_WANT_OMP			/*-----------------------------------*/
extern int phgMaxThreads_dummy, phgThreadId_dummy;
# define phgMaxThreads	phgMaxThreads_dummy
# define phgThreadId	phgThreadId_dummy
#else	/* defined(DONT_WANT_OMP) */	/*-----------------------------------*/
extern int phgMaxThreads, phgThreadId;
# if USE_OMP
# pragma omp threadprivate(phgThreadId)
# endif	/* USE_OMP */
#endif	/* defined(DONT_WANT_OMP) */	/*-----------------------------------*/

extern INT	_phg_debug_i;
extern FLOAT	_phg_debug_f;
extern void *	_phg_debug_p;
# if USE_OMP
# pragma omp threadprivate(_phg_debug_i, _phg_debug_f, _phg_debug_p)
# endif	/* USE_OMP */

typedef enum {FOR_NONE, FOR_OWNED, FOR_UNOWNED, FOR_ALL} FORALL_TYPE;

/* The macro ForAllElements_() loops on the elements in the submesh */

#define ForAllElements_(g, e, for_what) \
    /*if ((g)->i >= 0) phgError(1, "(%s:%d) nested use of ForAllElements detected.\n", __FILE__, __LINE__);*/ \
    for ((g)->i = 0; \
	 (e = _phg_step_element(g, &(g)->i, for_what)) != NULL; \
	 (g)->i++)

#define ForAllElements(g, e) ForAllElements_(g, e, FOR_OWNED)

/* The macro ForAllElements2() uses an user-provided loop index,
 * so it can be safely nested */
#define ForAllElements2(g, e, i) \
    for (i = 0; (e = _phg_step_element(g, &i, FOR_OWNED)) != NULL; i++)

/* The macro pair ForAllElementsBegin and ForAllElementsEnd should be used
 * with OpenMP (which only allows the following for-loop in the simplest form).
 * For example:
 *
 * 	#pragma omp parallel for private(e, ...)
 * 	ForAllElementsBegin(g, e)
 * 	    // the loop body
 * 	ForAllElementsEnd
 *
 * (note: 'e' must be included in the private list) */
extern INT _phg_omp_i;
#define ForAllElementsBegin_(g, e, for_what) \
    for (_phg_omp_i = \
	    (for_what == FOR_UNOWNED ? \
		((g)->halo == NULL ? (g)->nelem : (g)->halo->nelem_bak) : 0);\
	 _phg_omp_i < \
	    (for_what == FOR_NONE ? 0 : \
	     for_what != FOR_OWNED ? (g)->nelem : \
		((g)->halo == NULL ? (g)->nelem : (g)->halo->nelem_bak)); \
	 _phg_omp_i++) { \
	if ((e = (g)->elems[_phg_omp_i]) == NULL) \
	    continue;
#define ForAllElementsBegin(g, e) ForAllElementsBegin_(g, e, FOR_OWNED)
#define ForAllElementsEnd	} _phg_omp_i = -1;

#define ForAllFaces(g, e, face) \
    if ((g)->nface > 0 && (g)->Fmap_elems == NULL) { \
	VEF_MAP *vef = phgDofSetupVEFMap_(g, NULL, FACE_FLAG, FOR_OWNED); \
	phgDofFreeVEFMap(&vef); \
    } \
    for ((g)->i = 0; \
	 ((g)->i = _phg_step_face(g,(g)->i, &e, &face)) >= 0; \
	 (g)->i++)

#define ForAllFacesBegin(g, e, face) \
    for (_phg_omp_i = 0; \
	 _phg_omp_i < ((g)->halo == NULL ? (g)->nface : (g)->halo->nface_bak); \
	 _phg_omp_i++) { \
	if (!(g->types_face[_phg_omp_i] & OWNER)) \
	    continue; \
	e = g->Fmap_elems[_phg_omp_i]; \
	face = g->Fmap_faces[_phg_omp_i];
#define ForAllFacesEnd	} _phg_omp_i = -1;

extern FLOAT _phg_round_tmp;
extern INT phgVerbosity;
extern int phgTraverseDepth;
extern INT phgTraverseIndex;
extern ELEMENT *_phg_step_element(GRID *g, INT *i, FORALL_TYPE for_what);
extern INT _phg_step_face(GRID *g, INT i, ELEMENT **e, BYTE *face);
extern BOOLEAN phgInitialized;
extern BOOLEAN phgMasterSlave;
extern int _phg_bcmap(GRID *g, int input_bc);
extern BOOLEAN _phg_pause;

/* common auxiliary functions */
int phgCompFLOAT(const void *p0, const void *p1);	/* compare 2 FLOAT's */
int phgCompINT(const void *p0, const void *p1);		/* compares 2 'INT's */
int phgCompint(const void *p0, const void *p1);		/* compares 2 'int's */
INT phgBinarySearch(INT n, const void *table, const INT *ordering,
	const void *key, size_t size, int (*comp)(const void *, const void *));
INT phgBinarySearchINT(INT n, const INT *table, const INT *ordering, INT key);
/*#define phgBinarySearchINT(n, table, ordering, key) \
	phgBinarySearch(n, table, ordering, &(key), sizeof(INT), phgCompINT)*/

#define FunctionEntry_head						\
    double phg_t0[3];							\
    const char *phg_fname_save = phgCurrentFunctionName_;

#define FunctionEntry_tail						\
    phgCurrentFunctionName_ = __func__;					\
    phgGetTime(phg_t0);

#define FunctionEntry	FunctionEntry_head FunctionEntry_tail

/* Note: the dummy 'if (...)' here allows not including Return
 * in curly braces */
#define Return if ((phgCurrentFunctionName_ = phg_fname_save) || TRUE) return

#define PrintTime(vl) 							\
    if (vl >= 0 && vl <= phgVerbosity && IS_RANK0) {			\
	double phg_t1[3];						\
	phgGetTime(phg_t1);						\
	phgInfo(vl, "user=%f, sys=%f, wall=%f\n",			\
			phg_t1[0] - phg_t0[0], phg_t1[1] - phg_t0[1],	\
			phg_t1[2] - phg_t0[2]);				\
	phg_t0[0] = phg_t1[0]; phg_t0[1] = phg_t1[1]; phg_t0[2] = phg_t1[2]; \
    }

#if DEBUG
# define FreeAtExit_(ptr, level) {			\
    static BOOLEAN registered = FALSE;			\
    assert(phgThreadId == 0);				\
    if (!registered) {					\
	phgFreeAtExit((void **)(void *)&ptr, level);	\
	registered = TRUE;				\
    }							\
}
# define FreeAtExitNoCheck_(ptr, level) \
    phgFreeAtExit((void **)(void *)&ptr, level)
#else
# define FreeAtExit_(ptr, level)
# define FreeAtExitNoCheck_(ptr, level)
#endif
#define FreeAtExit(ptr)	FreeAtExit_(ptr, 0)
#define FreeAtExitNoCheck(ptr)	FreeAtExitNoCheck_(ptr, 0)

void phgTrapSignals_(BOOLEAN init);
void phgInit(int *argc, char ***argv);
void phgSetVerbosity_(GRID *g);
int  phgSetVerbosity(int verbosity);
void phgPause(int seconds);
void phgFinalize(void);
void phgInitElements(ELEMENT *s, INT count);
ELEMENT *phgNewElements(INT count);
ELEMENT *phgReallocElements(ELEMENT *s, INT oldcount, INT newcount);
COORD *phgNewVertices(INT count);
void phgNewGrid_(GRID *g);
GRID *phgNewGrid(int flags);
int phgGridId(const GRID *g);
#define phgGridIsValid(g)	(phgGridId(g) >= 0)
void *phgAlloc(size_t size);
void *phgCalloc(size_t nmemb, size_t size);

void *phgRealloc0(void *ptr, size_t size);
void *phgRealloc__(void *ptr, size_t size, size_t oldsize,
					const char *file, int line);

/* The macro phgRealloc_ is the function actually used in PHG's lib */
#if REALLOC_HACK
 /* For MVAPICH2 on DeepComp7000: use phgRealloc__ which falls back to malloc()
  *	+ memcpy() if realloc() fails (realloc() may fail when np > 64) */
# define phgRealloc_(ptr, size, oldsize) \
	 phgRealloc__(ptr, size, oldsize, __FILE__, __LINE__)
#else	/* REALLOC_HACK */
 /* Default: use phgRealloc0 which relies on realloc() ('oldsize' not used). */
# define phgRealloc_(ptr, size, oldsize) phgRealloc0(ptr, size)
#endif	/* REALLOC_HACK */

void phgFree(void *ptr);
int  phgMemcmp(const void *v0, const void *v1, size_t n);
void phgFreeAtExit(void **ptr, int level);
void phgWarning(const char *fmt, ...);
void phgAbort(int code);
void phgError(int code, const char *fmt, ...);
void phgInfo(int verbose_level, const char *fmt, ...);
int phgPrintf(const char *fmt, ...);
int phgPrintf2(GRID *g, const char *fmt, ...);
FILE *phgFOpen(GRID *g, const char *path, const char *mode);
int phgFPrintf(GRID *g, FILE *fp, char *fmt, ...);
int phgFClose(GRID *g, FILE *fp);
double phgGetTime(double tarray[]);

BYTE phgOppositeFace(GRID *g, ELEMENT *e, BYTE v, ELEMENT *e_op);
void phgFreeElement(ELEMENT **eptr);
void phgFreeGrid_(GRID *g);
void phgFreeGrid(GRID **gptr);
void phgTraverseElements(GRID *g, BOOLEAN(*cb)CB_ARGS(e));
void phgTraverseAllElements(GRID *g, BOOLEAN(*cb) CB_ARGS(e));
void phgCheckConformity(GRID *g);
void phgRemoveDegeneratedElements(GRID *g);
void phgUpdateNeighbours(GRID *g);
void phgUpdateEdges(GRID *g);
#if Dim == 3
void phgUpdateFaces(GRID *g);
#endif
void phgUpdateElementIndices(GRID *g, INT nleaf_old, INT nleaf_global_old);
int  phgGetPatch(GRID *g, ELEMENT *e, BYTE v0, BYTE v1, ELEMENT ***p);
void phgDumpElement(GRID *g, ELEMENT *e);
void phgDumpPatch(GRID *g, ELEMENT *e, int v0, int v1);
void phgDumpGridInfo(GRID *g);
void phgDumpGrid(GRID *g);
void phgDumpTree(GRID *g);
void phgDumpBranch(GRID *g, ELEMENT *e);

/* struct mapping input bc in the interval [min,max] to PHG type b */
typedef struct {
    int		min, max;
    BTYPE	btype;
} BCMAP;
int phgLoadBCmap(const char *mesh_fn);
int phgMapBC(int input_bc);
const char *phgBTypesString(BTYPE btype);
int phgBTypesValue(const char *btype_string);

void phgImportSetBdryMapFunc(BDRY_MAP_FUNC func);
BTYPE phgImportSetDefaultBdryType(BTYPE type);
BOOLEAN phgImport_(GRID *g, const char *filename, MPI_Comm comm_world);
#define phgImport(g, filename, distr) \
	phgImport_(g, filename, (distr) ? phgComm : MPI_COMM_SELF)

#if ALLOW_CURVED_BOUNDARY
void phgInitBoundaryFunctions(GRID *g);
#endif
FLOAT phgGetBoundaryError(GRID *g, ELEMENT *e);

GRID *phgDupGrid(GRID *g, BOOLEAN flatten);

#if USE_MPI
INT  phgCountRNeighbours(int nprocs, RNEIGHBOUR *list, INT count,
				int *cnts, int *dspls);
#endif	/* USE_MPI */

/* Macro to get the local index of a given vertex of a given edge,
 * 'edge' is the local edge no, 'vertex' is 0 or 1 */
#if PHG_GRID_TYPE == GRID_TYPE_TET
#define GetEdgeVertex(edge, vertex) ((int)phgTetEdgeVertex[edge][vertex])
extern char phgTetEdgeVertex[NEdge][2];
#elif PHG_GRID_TYPE == GRID_TYPE_HEX
#define GetEdgeVertex(edge, vertex) ((int)phgHexEdgeVertex[edge][vertex])
extern char phgHexEdgeVertex[NEdge][2];
#elif PHG_GRID_TYPE == GRID_TYPE_MIX
#define GetEdgeVertex(edge, vertex) (phgGetEdgeVertex(edge, vertex))
BYTE phgGetEdgeVertex(int edge, int vertex);
#endif

#define GetFaceBTYPE(g, e, face)	((e)->bound_type[face] & (BDRY_MASK))

/* Macro to get the local index of a given vertex of a given face,
 * 'face' is the local face no, 'vertex' is 0, 1, or 2 */
#if PHG_GRID_TYPE == GRID_TYPE_TET
#define GetFaceVertex(face, vertex) ((int)phgTetFaceVertex[face][vertex])
extern char phgTetFaceVertex[NFace][NVertFace];
#elif PHG_GRID_TYPE == GRID_TYPE_HEX
#define GetFaceVertex(face, vertex) ((int)phgHexFaceVertex[face][vertex])
extern char phgHexFaceVertex[NFace][NVertFace];
#elif PHG_GRID_TYPE == GRID_TYPE_MIX
#define GetFaceVertex(face, vertex) (phgGetFaceVertex(face, vertex))
BYTE phgGetFaceVertex(int face, int vert);
#endif

/* Macro to get the local index of the edge on the two given vertices */
#if PHG_GRID_TYPE == GRID_TYPE_TET
#define GetEdgeNo(v0, v1) ((int)phgTetEdgeNo[v0][v1])
extern char phgTetEdgeNo[NVert][NVert];
#elif PHG_GRID_TYPE == GRID_TYPE_HEX
#define GetEdgeNo(v0, v1) ((int)phgHexEdgeNo[v0][v1])
extern char phgHexEdgeNo[NVert][NVert];
#elif PHG_GRID_TYPE == GRID_TYPE_MIX
#define GetEdgeNo(v0, v1) (phgGetEdgeNo(v0, v1))
BYTE phgGetEdgeNo(int vert, int vert);
#endif

/* Macro to get the local index of the edge on the given faces */
#if PHG_GRID_TYPE == GRID_TYPE_TET
#define GetFaceEdge(face, edge) ((int)phgTetFaceEdge[face][edge])
extern char phgTetFaceEdge[NFace][NEdgeFace];
#elif PHG_GRID_TYPE == GRID_TYPE_HEX
#define GetFaceEdge(face, edge) ((int)phgHexFaceEdge[face][edge])
extern char phgHexFaceEdge[NFace][NEdgeFace];
#elif PHG_GRID_TYPE == GRID_TYPE_MIX
#define GetFaceEdge(face, edge) (phgGetFaceEdge(face, edge))
BYTE phgGetFaceEdge(int face, int edge);
#endif

/* utility function */
void phgUpdateBoundaryTypes(GRID *g);

#ifdef __cplusplus
}
#endif

/* Some utility macros for FLOAT */

#include <float.h>

#define IsZero(a) (Fabs(a) < 10.0 * FLOAT_MIN)
#define IsOne(a)  (Fabs(a - 1.0) < 10.0 * FLOAT_EPSILON)

/* Macros rounding FLOAT to various integer types.
 * Note: C99 built-in function not used because it can only round to int */
#define Round_(a, type)	\
    ((_phg_round_tmp = (a)) < 0. ? -(type)(-_phg_round_tmp + .5) \
				 :  (type)(_phg_round_tmp + .5))
#define Round(a)	Round_(a, int)
#define RoundINT(a)	Round_(a, INT)
#define Roundsize_t(a)	Round_(a, size_t)

#define PHG_UTILS_H
#endif
