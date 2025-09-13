#ifndef lint
static const char yysccsid[] = "@(#)yaccpar	1.9 (Berkeley) 02/21/93";
#endif

#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define YYPATCH 20120115

#define YYEMPTY        (-1)
#define yyclearin      (yychar = YYEMPTY)
#define yyerrok        (yyerrflag = 0)
#define YYRECOVERING() (yyerrflag != 0)

#define YYPREFIX "yy"

#define YYPURE 0

#line 2 "expression.y"
/*
 * Runtime evaluation of analytical functions (expressions).
 *
 * $Id: expression.c,v 1.14 2014/06/17 02:39:48 zlb Exp $
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

#line 47 "expression.y"
#ifdef YYSTYPE
#undef  YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
#endif
#ifndef YYSTYPE_IS_DECLARED
#define YYSTYPE_IS_DECLARED 1
typedef union {EXPR *exp;} YYSTYPE;
#endif /* !YYSTYPE_IS_DECLARED */
#line 72 "y.tab.c"

/* compatibility with bison */
#ifdef YYPARSE_PARAM
/* compatibility with FreeBSD */
# ifdef YYPARSE_PARAM_TYPE
#  define YYPARSE_DECL() yyparse(YYPARSE_PARAM_TYPE YYPARSE_PARAM)
# else
#  define YYPARSE_DECL() yyparse(void *YYPARSE_PARAM)
# endif
#else
# define YYPARSE_DECL() yyparse(void)
#endif

/* Parameters sent to lex. */
#ifdef YYLEX_PARAM
# define YYLEX_DECL() yylex(void *YYLEX_PARAM)
# define YYLEX yylex(YYLEX_PARAM)
#else
# define YYLEX_DECL() yylex(void)
# define YYLEX yylex()
#endif

/* Parameters sent to yyerror. */
#ifndef YYERROR_DECL
#define YYERROR_DECL() yyerror(const char *s)
#endif
#ifndef YYERROR_CALL
#define YYERROR_CALL(msg) yyerror(msg)
#endif

extern int YYPARSE_DECL();

#define UMINUS 257
#define NUMBER 258
#define SIN 259
#define COS 260
#define TAN 261
#define PI 262
#define EXP 263
#define LOG 264
#define SQRT 265
#define POW 266
#define X 267
#define Y 268
#define Z 269
#define T 270
#define YYERRCODE 256
static const short yylhs[] = {                           -1,
    0,    0,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
    1,    1,
};
static const short yylen[] = {                            2,
    0,    1,    1,    1,    4,    4,    4,    4,    4,    4,
    6,    1,    1,    1,    1,    3,    3,    3,    3,    3,
    3,    2,
};
static const short yydefred[] = {                         0,
    0,    3,    0,    0,    0,    4,    0,    0,    0,    0,
   12,   13,   14,   15,    0,    0,    0,   22,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,   16,    0,
    0,    0,    0,   21,    5,    6,    7,    8,    9,   10,
    0,    0,   11,
};
static const short yydgoto[] = {                         16,
   17,
};
static const short yysindex[] = {                       -40,
  -40,    0,  -26,  -21,  -19,    0,  -14,  -12,   -7,   -5,
    0,    0,    0,    0,  -40,    0,   28,    0,  -40,  -40,
  -40,  -40,  -40,  -40,  -40,  -39,  -40,  -40,  -40,  -40,
  -40,  -32,  -25,  -18,  -11,   -4,    3,   21,    0,  -35,
  -35,  -93,  -93,    0,    0,    0,    0,    0,    0,    0,
  -40,   11,    0,
};
static const short yyrindex[] = {                        40,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,   42,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,   96,
  111,   82,   91,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,
};
static const short yygindex[] = {                         0,
   79,
};
#define YYTABLESIZE 230
static const short yytable[] = {                         15,
   31,   39,   29,   27,    1,   28,   29,   30,   45,   29,
   27,   30,   28,   19,   30,   46,   29,   27,   20,   28,
   21,   30,   47,   29,   27,   22,   28,   23,   30,   48,
   29,   27,   24,   28,   25,   30,   49,   29,   27,    1,
   28,    2,   30,   50,   29,   27,    0,   28,    0,   30,
    0,   53,   29,   27,   31,   28,    0,   30,   31,    0,
    0,   31,   29,   27,   51,   28,    0,   30,   31,   29,
   27,    0,   28,    0,   30,   31,    0,    0,    0,   18,
    0,   19,   31,    0,    0,    0,    0,    0,    0,   31,
   20,    0,    0,   26,    0,   17,   31,   32,   33,   34,
   35,   36,   37,   38,   31,   40,   41,   42,   43,   44,
   18,    0,    0,    0,   31,    0,    0,    0,    0,    0,
    0,   31,   19,   19,   19,   19,   19,    0,   19,   52,
    0,   20,   20,   20,   20,   20,   17,   20,   17,   17,
   17,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,   18,    0,   18,   18,   18,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    2,    3,    4,
    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,
};
static const short yycheck[] = {                         40,
   94,   41,   42,   43,   45,   45,   42,   47,   41,   42,
   43,   47,   45,   40,   47,   41,   42,   43,   40,   45,
   40,   47,   41,   42,   43,   40,   45,   40,   47,   41,
   42,   43,   40,   45,   40,   47,   41,   42,   43,    0,
   45,    0,   47,   41,   42,   43,   -1,   45,   -1,   47,
   -1,   41,   42,   43,   94,   45,   -1,   47,   94,   -1,
   -1,   94,   42,   43,   44,   45,   -1,   47,   94,   42,
   43,   -1,   45,   -1,   47,   94,   -1,   -1,   -1,    1,
   -1,    0,   94,   -1,   -1,   -1,   -1,   -1,   -1,   94,
    0,   -1,   -1,   15,   -1,    0,   94,   19,   20,   21,
   22,   23,   24,   25,   94,   27,   28,   29,   30,   31,
    0,   -1,   -1,   -1,   94,   -1,   -1,   -1,   -1,   -1,
   -1,   94,   41,   42,   43,   44,   45,   -1,   47,   51,
   -1,   41,   42,   43,   44,   45,   41,   47,   43,   44,
   45,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   41,   -1,   43,   44,   45,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,  258,  259,  260,
  261,  262,  263,  264,  265,  266,  267,  268,  269,  270,
};
#define YYFINAL 16
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 270
#if YYDEBUG
static const char *yyname[] = {

"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,"'('","')'","'*'","'+'","','","'-'",0,"'/'",0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"'^'",0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
"UMINUS","NUMBER","SIN","COS","TAN","PI","EXP","LOG","SQRT","POW","X","Y","Z",
"T",
};
static const char *yyrule[] = {
"$accept : line",
"line :",
"line : expr",
"expr : NUMBER",
"expr : PI",
"expr : SIN '(' expr ')'",
"expr : COS '(' expr ')'",
"expr : TAN '(' expr ')'",
"expr : EXP '(' expr ')'",
"expr : LOG '(' expr ')'",
"expr : SQRT '(' expr ')'",
"expr : POW '(' expr ',' expr ')'",
"expr : X",
"expr : Y",
"expr : Z",
"expr : T",
"expr : '(' expr ')'",
"expr : expr '+' expr",
"expr : expr '-' expr",
"expr : expr '*' expr",
"expr : expr '/' expr",
"expr : expr '^' expr",
"expr : '-' expr",

};
#endif

int      yydebug;
int      yynerrs;

int      yyerrflag;
int      yychar;
YYSTYPE  yyval;
YYSTYPE  yylval;

/* define the initial stack-sizes */
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH  YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 10000
#define YYMAXDEPTH  500
#endif
#endif

#define YYINITSTACKSIZE 500

typedef struct {
    unsigned stacksize;
    short    *s_base;
    short    *s_mark;
    short    *s_last;
    YYSTYPE  *l_base;
    YYSTYPE  *l_mark;
} YYSTACKDATA;
/* variables for the parser stack */
static YYSTACKDATA yystack;
#line 86 "expression.y"

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
#line 783 "y.tab.c"

#if YYDEBUG
#include <stdio.h>		/* needed for printf */
#endif

#include <stdlib.h>	/* needed for malloc, etc */
#include <string.h>	/* needed for memset */

/* allocate initial stack or double stack size, up to YYMAXDEPTH */
static int yygrowstack(YYSTACKDATA *data)
{
    int i;
    unsigned newsize;
    short *newss;
    YYSTYPE *newvs;

    if ((newsize = data->stacksize) == 0)
        newsize = YYINITSTACKSIZE;
    else if (newsize >= YYMAXDEPTH)
        return -1;
    else if ((newsize *= 2) > YYMAXDEPTH)
        newsize = YYMAXDEPTH;

    i = data->s_mark - data->s_base;
    newss = (short *)realloc(data->s_base, newsize * sizeof(*newss));
    if (newss == 0)
        return -1;

    data->s_base = newss;
    data->s_mark = newss + i;

    newvs = (YYSTYPE *)realloc(data->l_base, newsize * sizeof(*newvs));
    if (newvs == 0)
        return -1;

    data->l_base = newvs;
    data->l_mark = newvs + i;

    data->stacksize = newsize;
    data->s_last = data->s_base + newsize - 1;
    return 0;
}

#if YYPURE || defined(YY_NO_LEAKS)
static void yyfreestack(YYSTACKDATA *data)
{
    free(data->s_base);
    free(data->l_base);
    memset(data, 0, sizeof(*data));
}
#else
#define yyfreestack(data) /* nothing */
#endif

#define YYABORT  goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR  goto yyerrlab

int
YYPARSE_DECL()
{
    int yym, yyn, yystate;
#if YYDEBUG
    const char *yys;

    if ((yys = getenv("YYDEBUG")) != 0)
    {
        yyn = *yys;
        if (yyn >= '0' && yyn <= '9')
            yydebug = yyn - '0';
    }
#endif

    yynerrs = 0;
    yyerrflag = 0;
    yychar = YYEMPTY;
    yystate = 0;

#if YYPURE
    memset(&yystack, 0, sizeof(yystack));
#endif

    if (yystack.s_base == NULL && yygrowstack(&yystack)) goto yyoverflow;
    yystack.s_mark = yystack.s_base;
    yystack.l_mark = yystack.l_base;
    yystate = 0;
    *yystack.s_mark = 0;

yyloop:
    if ((yyn = yydefred[yystate]) != 0) goto yyreduce;
    if (yychar < 0)
    {
        if ((yychar = YYLEX) < 0) yychar = 0;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
    }
    if ((yyn = yysindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, yystate, yytable[yyn]);
#endif
        if (yystack.s_mark >= yystack.s_last && yygrowstack(&yystack))
        {
            goto yyoverflow;
        }
        yystate = yytable[yyn];
        *++yystack.s_mark = yytable[yyn];
        *++yystack.l_mark = yylval;
        yychar = YYEMPTY;
        if (yyerrflag > 0)  --yyerrflag;
        goto yyloop;
    }
    if ((yyn = yyrindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
        yyn = yytable[yyn];
        goto yyreduce;
    }
    if (yyerrflag) goto yyinrecovery;

    yyerror("syntax error");

    goto yyerrlab;

yyerrlab:
    ++yynerrs;

yyinrecovery:
    if (yyerrflag < 3)
    {
        yyerrflag = 3;
        for (;;)
        {
            if ((yyn = yysindex[*yystack.s_mark]) && (yyn += YYERRCODE) >= 0 &&
                    yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yystack.s_mark, yytable[yyn]);
#endif
                if (yystack.s_mark >= yystack.s_last && yygrowstack(&yystack))
                {
                    goto yyoverflow;
                }
                yystate = yytable[yyn];
                *++yystack.s_mark = yytable[yyn];
                *++yystack.l_mark = yylval;
                goto yyloop;
            }
            else
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *yystack.s_mark);
#endif
                if (yystack.s_mark <= yystack.s_base) goto yyabort;
                --yystack.s_mark;
                --yystack.l_mark;
            }
        }
    }
    else
    {
        if (yychar == 0) goto yyabort;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
        yychar = YYEMPTY;
        goto yyloop;
    }

yyreduce:
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
    yym = yylen[yyn];
    if (yym)
        yyval = yystack.l_mark[1-yym];
    else
        memset(&yyval, 0, sizeof yyval);
    switch (yyn)
    {
case 2:
#line 62 "expression.y"
	{tree = /*dup_expr*/(yystack.l_mark[0].exp);}
break;
case 3:
#line 64 "expression.y"
	{tree = yyval.exp = new_expr(NUMBER, yystack.l_mark[0].exp, NULL);}
break;
case 4:
#line 65 "expression.y"
	{tree = yyval.exp = new_expr(PI, NULL, NULL);}
break;
case 5:
#line 66 "expression.y"
	{tree = yyval.exp = new_expr(SIN, yystack.l_mark[-1].exp, NULL);}
break;
case 6:
#line 67 "expression.y"
	{tree = yyval.exp = new_expr(COS, yystack.l_mark[-1].exp, NULL);}
break;
case 7:
#line 68 "expression.y"
	{tree = yyval.exp = new_expr(TAN, yystack.l_mark[-1].exp, NULL);}
break;
case 8:
#line 69 "expression.y"
	{tree = yyval.exp = new_expr(EXP, yystack.l_mark[-1].exp, NULL);}
break;
case 9:
#line 70 "expression.y"
	{tree = yyval.exp = new_expr(LOG, yystack.l_mark[-1].exp, NULL);}
break;
case 10:
#line 71 "expression.y"
	{tree = yyval.exp = new_expr(SQRT, yystack.l_mark[-1].exp, NULL);}
break;
case 11:
#line 72 "expression.y"
	{tree = yyval.exp = new_expr(POW, yystack.l_mark[-3].exp, yystack.l_mark[-1].exp);}
break;
case 12:
#line 73 "expression.y"
	{tree = yyval.exp = new_expr(X, NULL, NULL);}
break;
case 13:
#line 74 "expression.y"
	{tree = yyval.exp = new_expr(Y, NULL, NULL);}
break;
case 14:
#line 75 "expression.y"
	{tree = yyval.exp = new_expr(Z, NULL, NULL);}
break;
case 15:
#line 76 "expression.y"
	{tree = yyval.exp = new_expr(T, NULL, NULL);}
break;
case 16:
#line 77 "expression.y"
	{tree = yyval.exp = yystack.l_mark[-1].exp;}
break;
case 17:
#line 78 "expression.y"
	{tree = yyval.exp = new_expr('+', yystack.l_mark[-2].exp, yystack.l_mark[0].exp);}
break;
case 18:
#line 79 "expression.y"
	{tree = yyval.exp = new_expr('-', yystack.l_mark[-2].exp, yystack.l_mark[0].exp);}
break;
case 19:
#line 80 "expression.y"
	{tree = yyval.exp = new_expr('*', yystack.l_mark[-2].exp, yystack.l_mark[0].exp);}
break;
case 20:
#line 81 "expression.y"
	{tree = yyval.exp = new_expr('/', yystack.l_mark[-2].exp, yystack.l_mark[0].exp);}
break;
case 21:
#line 82 "expression.y"
	{tree = yyval.exp = new_expr('^', yystack.l_mark[-2].exp, yystack.l_mark[0].exp);}
break;
case 22:
#line 83 "expression.y"
	{tree = yyval.exp = new_expr(UMINUS, yystack.l_mark[0].exp, NULL);}
break;
#line 1073 "y.tab.c"
    }
    yystack.s_mark -= yym;
    yystate = *yystack.s_mark;
    yystack.l_mark -= yym;
    yym = yylhs[yyn];
    if (yystate == 0 && yym == 0)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        yystate = YYFINAL;
        *++yystack.s_mark = YYFINAL;
        *++yystack.l_mark = yyval;
        if (yychar < 0)
        {
            if ((yychar = YYLEX) < 0) yychar = 0;
#if YYDEBUG
            if (yydebug)
            {
                yys = 0;
                if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
                if (!yys) yys = "illegal-symbol";
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, yychar, yys);
            }
#endif
        }
        if (yychar == 0) goto yyaccept;
        goto yyloop;
    }
    if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
        yystate = yytable[yyn];
    else
        yystate = yydgoto[yym];
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yystack.s_mark, yystate);
#endif
    if (yystack.s_mark >= yystack.s_last && yygrowstack(&yystack))
    {
        goto yyoverflow;
    }
    *++yystack.s_mark = (short) yystate;
    *++yystack.l_mark = yyval;
    goto yyloop;

yyoverflow:
    yyerror("yacc stack overflow");

yyabort:
    yyfreestack(&yystack);
    return (1);

yyaccept:
    yyfreestack(&yystack);
    return (0);
}
