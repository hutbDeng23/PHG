%{
/*
 * Runtime evaluation of analytical functions (expressions).
 *
 * $Id: expression.y,v 1.9 2014/06/17 02:39:36 zlb Exp $
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phg/expression.h"
#if 0
/* Hacky: rename symbols which are multiply defined in MVAPICH2 */
#define yy_scan_string          yy_scan_string_dummy
#define yy_delete_buffer        yy_delete_buffer_dummy
#define yy_scan_buffer          yy_scan_buffer_dummy
#define yy_scan_bytes           yy_scan_bytes_dummy
#define yy_switch_to_buffer     yy_switch_to_buffer_dummy
#define yy_load_buffer_state    yy_load_buffer_state_dummy
#define yy_flush_buffer         yy_flush_buffer_dummy
#define yy_create_buffer        yy_create_buffer_dummy
#define yyrestart               yyrestart_dummy
#define yylex                   yylex_dummy
#define yy_init_buffer          yy_init_buffer_dummy
#define yyparse                 yyparse_dummy
#define yyerror                 yyerror_dummy
#define yyin                    yyin_dummy
#define yyout                   yyout_dummy
#endif	/* 0 */

#define pi	3.1415926535897932384626L

static int tmp;
#define P(a, b) ( fabs(b - (tmp = (int)b)) > 1e-10 ? pow(a,b) : \
    (tmp == 0 ? 1.0 : (tmp == 1 ? a : (tmp == 2) ? a * a : (tmp == 3 ? a*a*a : \
	(tmp == 4 ? a*a*a*a : (tmp == 5 ? a*a*a*a*a : (tmp == 6 ? a*a*a*a*a*a :\
	(tmp == 7 ? a*a*a*a*a*a*a : (tmp == 8 ? a*a*a*a*a*a*a*a : \
	pow(a,b))))))))) )

static double x_value = 0, y_value = 0, z_value = 0, t_value = 0;

static EXPR *tree = NULL;
static const char *yy_expression = NULL;

%}

%union {EXPR *exp;}

%left '+' '-'
%left '*' '/'
%left '^'

%nonassoc UMINUS

%token NUMBER SIN COS TAN PI EXP LOG SQRT POW X Y Z T
%start line
%type <exp> expr NUMBER

%%

line	: 
	| expr			{tree = /*dup_expr*/($1);}
	;
expr	: NUMBER		{tree = $$ = new_expr(NUMBER, $1, NULL);}
	| PI			{tree = $$ = new_expr(PI, NULL, NULL);}
	| SIN '(' expr ')'	{tree = $$ = new_expr(SIN, $3, NULL);}
	| COS '(' expr ')'	{tree = $$ = new_expr(COS, $3, NULL);}
	| TAN '(' expr ')'	{tree = $$ = new_expr(TAN, $3, NULL);}
	| EXP '(' expr ')'	{tree = $$ = new_expr(EXP, $3, NULL);}
	| LOG '(' expr ')'	{tree = $$ = new_expr(LOG, $3, NULL);}
	| SQRT '(' expr ')'	{tree = $$ = new_expr(SQRT, $3, NULL);}
	| POW '(' expr ',' expr ')' {tree = $$ = new_expr(POW, $3, $5);}
	| X			{tree = $$ = new_expr(X, NULL, NULL);}
	| Y			{tree = $$ = new_expr(Y, NULL, NULL);}
	| Z			{tree = $$ = new_expr(Z, NULL, NULL);}
	| T			{tree = $$ = new_expr(T, NULL, NULL);}
	| '(' expr ')'		{tree = $$ = $2;}
	| expr '+' expr		{tree = $$ = new_expr('+', $1, $3);}
	| expr '-' expr		{tree = $$ = new_expr('-', $1, $3);}
	| expr '*' expr		{tree = $$ = new_expr('*', $1, $3);}
	| expr '/' expr		{tree = $$ = new_expr('/', $1, $3);}
	| expr '^' expr		{tree = $$ = new_expr('^', $1, $3);}
	| '-' expr %prec UMINUS {tree = $$ = new_expr(UMINUS, $2, NULL);}
	;
%%

int yyparse(void);
int yylex(void);
int _phg_yywrap(void) {return 1;}
void yyerror(char *s) {fprintf(stderr, "%s\n", s);}

static EXPR *
dup_expr(EXPR *e0)
{
    EXPR *e;

    if (e0 == NULL)
	return NULL;

    e = malloc(sizeof(EXPR));
    memcpy(e, e0, sizeof(EXPR));

    switch (e->op) {
	case NUMBER:
	    e->arg1 = malloc(sizeof(double));
	    *(double *)(e->arg1) = *(double *)(e0->arg1);
	    break;
	case PI:
	case X:
	case Y:
	case Z:
	case T:
	    break;
	case SIN:
	case COS:
	case TAN:
	case EXP:
	case LOG:
	case SQRT:
	case UMINUS :
	    e->arg1 = dup_expr(e0->arg1);
	    break;
	case '+':
	case '-':
	case '*':
	case '/':
	case '^':
	case POW:
	    e->arg1 = dup_expr(e0->arg1);
	    e->arg2 = dup_expr(e0->arg2);
	    break;
	default:
	    fprintf(stderr, "dup_expr: unhandled operation %d\n", e->op);
	    exit(1);
    }

    return e;
}

static void
free_expr(EXPR *e)
{
    if (e == NULL)
	return;

    switch (e->op) {
	case NUMBER:
	    free(e->arg1);
	    break;
	case PI:
	case X:
	case Y:
	case Z:
	case T:
	    break;
	case SIN:
	case COS:
	case TAN:
	case EXP:
	case LOG:
	case SQRT:
	case UMINUS :
	    free_expr(e->arg1);
	    break;
	case '+':
	case '-':
	case '*':
	case '/':
	case '^':
	case POW:
	    free_expr(e->arg1);
	    free_expr(e->arg2);
	    break;
	default:
	    fprintf(stderr, "free_expr: unhandled operation %d\n", e->op);
	    exit(1);
    }

    free(e);
}

static EXPR *
new_expr(int op, void *arg1, void *arg2)
{
    EXPR *e = NULL, *a1 = arg1, *a2 = arg2;
#define is_NUMBER(op) (op == NUMBER)

    switch (op) {
	case NUMBER:
	    e = arg1;
	    break;
	case PI:
	    e = malloc(sizeof(EXPR));
	    e->arg1 = malloc(sizeof(double));
	    *(double *)(e->arg1) = pi;
	    e->op = NUMBER;
	    break;
	case X:
	case Y:
	case Z:
	case T:
	    e = malloc(sizeof(EXPR));
	    e->op = op;
	    break;
	case SIN:
	    if (!is_NUMBER(a1->op)) goto new1;
	    *(double *)(a1->arg1) = sin(*(double *)(a1->arg1));
	    e = arg1;
	    break;
	case COS:
	    if (!is_NUMBER(a1->op)) goto new1;
	    *(double *)(a1->arg1) = cos(*(double *)(a1->arg1));
	    e = arg1;
	    break;
	case TAN:
	    if (!is_NUMBER(a1->op)) goto new1;
	    *(double *)(a1->arg1) = tan(*(double *)(a1->arg1));
	    e = arg1;
	    break;
	case EXP:
	    if (!is_NUMBER(a1->op)) goto new1;
	    *(double *)(a1->arg1) = exp(*(double *)(a1->arg1));
	    e = arg1;
	    break;
	case LOG:
	    if (!is_NUMBER(a1->op)) goto new1;
	    *(double *)(a1->arg1) = log(*(double *)(a1->arg1));
	    e = arg1;
	    break;
	case SQRT:
	    if (!is_NUMBER(a1->op)) goto new1;
	    *(double *)(a1->arg1) = sqrt(*(double *)(a1->arg1));
	    e = arg1;
	    break;
	case UMINUS:
	    if (!is_NUMBER(a1->op)) goto new1;
	    *(double *)(a1->arg1) = -(*(double *)(a1->arg1));
	    e = arg1;
	    break;
    new1:
	    e = malloc(sizeof(EXPR));
	    e->arg1 = arg1;
	    e->op = op;
	    break;
	case '+':
	    if (!is_NUMBER(a1->op) || !is_NUMBER(a2->op)) goto new2;
	    *(double *)(a1->arg1) = (*(double *)(a1->arg1)) +
				    (*(double *)(a2->arg1));
	    free_expr(arg2);
	    e = arg1;
	    break;
	case '-':
	    if (!is_NUMBER(a1->op) || !is_NUMBER(a2->op)) goto new2;
	    *(double *)(a1->arg1) = (*(double *)(a1->arg1)) -
				    (*(double *)(a2->arg1));
	    free_expr(arg2);
	    e = arg1;
	    break;
	case '*':
	    if (!is_NUMBER(a1->op) || !is_NUMBER(a2->op)) goto new2;
	    *(double *)(a1->arg1) = (*(double *)(a1->arg1)) *
				    (*(double *)(a2->arg1));
	    free_expr(arg2);
	    e = arg1;
	    break;
	case '/':
	    if (!is_NUMBER(a1->op) || !is_NUMBER(a2->op)) goto new2;
	    *(double *)(a1->arg1) = (*(double *)(a1->arg1)) /
				    (*(double *)(a2->arg1));
	    free_expr(arg2);
	    e = arg1;
	    break;
	case '^':
	case POW:
	    if (!is_NUMBER(a1->op) || !is_NUMBER(a2->op)) goto new2;
	    *(double *)(a1->arg1) = P(*(double *)(a1->arg1),
				      *(double *)(a2->arg1));
	    free_expr(arg2);
	    e = arg1;
	    break;
    new2:
	    e = malloc(sizeof(EXPR));
	    e->arg1 = arg1;
	    e->arg2 = arg2;
	    e->op = op;
	    break;
	default:
	    fprintf(stderr, "new_expr: unhandled operation %d\n", op);
	    exit(1);
    }

    return e;
}

static double
eval_expr(EXPR *e)
{
    switch (e->op) {
	case NUMBER:
	    return *(double *)(e->arg1);
	case PI:
	    return pi;
	case X:
	    return x_value;
	case Y:
	    return y_value;
	case Z:
	    return z_value;
	case T:
	    return t_value;
	case SIN:
	    return sin(eval_expr(e->arg1));
	case COS:
	    return cos(eval_expr(e->arg1));
	case TAN:
	    return tan(eval_expr(e->arg1));
	case EXP:
	    return exp(eval_expr(e->arg1));
	case LOG:
	    return log(eval_expr(e->arg1));
	case SQRT:
	    return sqrt(eval_expr(e->arg1));
	case UMINUS :
	    return -eval_expr(e->arg1);
	case '+':
	    return eval_expr(e->arg1) + eval_expr(e->arg2);
	case '-':
	    return eval_expr(e->arg1) - eval_expr(e->arg2);
	case '*':
	    return eval_expr(e->arg1) * eval_expr(e->arg2);
	case '/':
	    return eval_expr(e->arg1) / eval_expr(e->arg2);
	case '^':
	case POW:
	    return P(eval_expr(e->arg1), eval_expr(e->arg2));
	default:
	    fprintf(stderr, "eval_expr: unhandled operation %d\n", e->op);
	    exit(1);
    }

    return 0;	/* to make gcc happy */
}

char *
phgDump3DFunction(EXPR *e)
{
    static int recur_depth = 0;
    static char work[4096];
    char *p1, *p2;
    double v;

    if (!recur_depth && e == NULL)
	return "";

    recur_depth++;
    switch (e->op) {
	case NUMBER:
	    if ((v = *(double *)(e->arg1)) >= 0)
		sprintf(work, "%0.16lg", v);
	    else
		sprintf(work, "(%0.16lg)", v);
	    break;
	case PI:
	    sprintf(work, "pi");
	    break;
	case X:
	    sprintf(work, "x");
	    break;
	case Y:
	    sprintf(work, "y");
	    break;
	case Z:
	    sprintf(work, "z");
	    break;
	case T:
	    sprintf(work, "t");
	    break;
	case SIN:
	    sprintf(work, "sin(%s)", p1 = phgDump3DFunction(e->arg1));
	    free(p1);
	    break;
	case COS:
	    sprintf(work, "cos(%s)", p1 = phgDump3DFunction(e->arg1));
	    free(p1);
	    break;
	case TAN:
	    sprintf(work, "tan(%s)", p1 = phgDump3DFunction(e->arg1));
	    free(p1);
	    break;
	case EXP:
	    sprintf(work, "exp(%s)", p1 = phgDump3DFunction(e->arg1));
	    free(p1);
	    break;
	case LOG:
	    sprintf(work, "log(%s)", p1 = phgDump3DFunction(e->arg1));
	    free(p1);
	    break;
	case SQRT:
	    sprintf(work, "sqrt(%s)", p1 = phgDump3DFunction(e->arg1));
	    free(p1);
	    break;
	case UMINUS :
	    sprintf(work, "-(%s)", p1 = phgDump3DFunction(e->arg1));
	    free(p1);
	    break;
	case '+':
	    sprintf(work, "(%s+%s)", p1 = phgDump3DFunction(e->arg1),
				     p2 = phgDump3DFunction(e->arg2));
	    free(p1);
	    free(p2);
	    break;
	case '-':
	    sprintf(work, "(%s-%s)", p1 = phgDump3DFunction(e->arg1),
				     p2 = phgDump3DFunction(e->arg2));
	    free(p1);
	    free(p2);
	    break;
	case '*':
	    sprintf(work, "(%s*%s)", p1 = phgDump3DFunction(e->arg1),
				     p2 = phgDump3DFunction(e->arg2));
	    free(p1);
	    free(p2);
	    break;
	case '/':
	    sprintf(work, "(%s/%s)", p1 = phgDump3DFunction(e->arg1),
				     p2 = phgDump3DFunction(e->arg2));
	    free(p1);
	    free(p2);
	    break;
	case '^':
	case POW:
	    sprintf(work, "%s^%s", p1 = phgDump3DFunction(e->arg1),
				   p2 = phgDump3DFunction(e->arg2));
	    free(p1);
	    free(p2);
	    break;
	default:
	    fprintf(stderr, "phgDump3DFunction: unhandled operation %d\n",
			e->op);
	    exit(1);
    }
    recur_depth--;

    if (recur_depth > 0)
	return strdup(work);

    if (work[0] == '(') {
	/* remove toplevel '(' ')' */
	work[strlen(work) - 1] = '\0';
	return work + 1;
    }
    else {
	return work;
    }
}

/*----------------------- User functions -------------------------*/

EXPR *
phgDefine3DFunction(const char *function_expression)
{
    EXPR *e;

    free_expr(tree);
    tree = NULL;
    yy_expression = function_expression;
    if (yyparse() != 0) {
	free_expr(tree);
	tree = NULL;
    }
    e = tree;
    tree = NULL;
    return e;
}

double
phgEvaluate3DFunction(EXPR *tree, double x, double y, double z, double t)
{
    if (tree == NULL) {
	fprintf(stderr, "phgEvaluate3DFunction: function undefined\n");
	return 0.0;
    }
    x_value = x; y_value = y; z_value = z; t_value = t;
    return eval_expr(tree);
}

EXPR *
phgDup3DFunction(EXPR *tree)
{
    return dup_expr(tree);
}

void
phgFree3DFunction(EXPR *tree)
{
    free_expr(tree);
}

#ifdef TEST

// #define FUNC sin(x) + sin(y) + sin(z)
#define FUNC sin(x) + (y*y + z*z) + 12*sin(pi/4)
#define N	100

#ifdef __STRING
# define STRING(f) __STRING(f)
#else
# define STRING(f) "sin(x) + (y*y + z*z)"
#endif

static double
F(double x, double y, double z)
{
    return FUNC;
}

#include <sys/resource.h>
int
main(int argc, char *argv[])
{
    struct rusage RU;
    double sum, time0, time1, h = 1.0/N, x, y, z;
    int i, j, k;
    EXPR *func[4];

    for (i = 0; i < sizeof(func)/sizeof(func[0]); i++) {
	func[i] = phgDefine3DFunction(STRING(FUNC));
    }
    printf("Function: %s\n", phgDump3DFunction(func[0]));

    getrusage(RUSAGE_SELF, &RU);
    time0 = RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec*1e-6 +
	    RU.ru_stime.tv_sec + (double)RU.ru_stime.tv_usec*1e-6;
    sum = 0;
    for (i = 0; i <= N; i++) {
	x = i * h;
	for (j = 0; j <= N; j++) {
	    y = j * h;
	    for (k = 0; k <= N; k++) {
		z = k * h;
		sum += phgEvaluate3DFunction(func[0], x, y, z, 0.0);
	    }
	}
    }
    getrusage(RUSAGE_SELF, &RU);
    time1 = RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec*1e-6 +
	    RU.ru_stime.tv_sec + (double)RU.ru_stime.tv_usec*1e-6;
    printf("Dynamic:  %0.16lg, %0.16lg secs\n", sum, time1 - time0);

    getrusage(RUSAGE_SELF, &RU);
    time0 = RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec*1e-6 +
	    RU.ru_stime.tv_sec + (double)RU.ru_stime.tv_usec*1e-6;
    sum = 0;
    for (i = 0; i <= N; i++) {
	x = i * h;
	for (j = 0; j <= N; j++) {
	    y = j * h;
	    for (k = 0; k <= N; k++) {
		z = k * h;
		sum += F(x, y, z);
	    }
	}
    }
    getrusage(RUSAGE_SELF, &RU);
    time1 = RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec*1e-6 +
	    RU.ru_stime.tv_sec + (double)RU.ru_stime.tv_usec*1e-6;
    printf("Compiled: %0.16lg, %0.16lg secs\n", sum, time1 - time0);

    for (i = 0; i < sizeof(func)/sizeof(func[0]); i++)
	phgFree3DFunction(func[i]);

    return 0;
}
#endif

/* include expression.l */
#include "expression.inc"
