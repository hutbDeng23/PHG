%{
/* $Id: expression.l,v 1.1 2006/02/11 04:12:37 zlb Exp $ */
#define YY_NO_UNPUT	/* to suppress warning "yyunput unused" */
#if 1
#ifdef YY_INPUT
# undef YY_INPUT
#endif
#define YY_INPUT(buf, result, max_size)				\
    {								\
	int n = 0;						\
	if (yy_expression != NULL) {				\
	    while (max_size > n && *yy_expression != '\0') {	\
		buf[n++] = *(yy_expression++);			\
	    }							\
	    result = n;						\
	}							\
	else {							\
	    result = fread(buf, 1, max_size, yyin);		\
	}							\
    }
#else
static int getc_bak(FILE *fin) {return getc(fin);}
#ifdef getc
# undef getc
#endif
#define getc(fin) ((yy_expression == NULL) ? getc_bak(fin) : \
	(*yy_expression == '\0' ? EOF : *(yy_expression++)))
#endif
%}

%%

\.[0-9]+([eE][+-]?[0-9]+)? |
[0-9]+(\.[0-9]*)?([eE][+-]?[0-9]+)? {
			yylval.exp = malloc(sizeof(EXPR));
			yylval.exp->op = NUMBER;
			yylval.exp->arg1 = malloc(sizeof(double));
			*(double *)(yylval.exp->arg1) = atof(yytext);
			return NUMBER;}
"**"			{return '^';}
[-+*/^(),]		{return *yytext;}
pi			{return PI;}
sin			{return SIN;}
cos			{return COS;}
tan			{return TAN;}
exp			{return EXP;}
log			{return LOG;}
sqrt			{return SQRT;}
pow			{return POW;}
x			{return X;}
y			{return Y;}
z			{return Z;}
t			{return T;}
[ \t\n]			/* do nothing */
.			{fprintf(stderr, "Unknown character (%c)\n", *yytext);}

%%
