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

/* $Id: option.h,v 1.3 2021/05/05 11:55:08 zlb Exp $ */

#ifndef PHG_OPTION_H

/* All options are stored in a list 'options'. When parsing a cmdline option,
 * the corresponding variable, *o->var, is updated with the given value.
 * o->var points to the location of the variable with one of the following
 * types:
 *
 * - For VT_INT, o->var is intepreted as (INT *)o->var. The operators
 *   '+=', '-=', '*=', '/=', '%=', '|=', '&=' and '^=' are allowed.
 *
 * - For VT_FLOAT, o->var is intepreted as (FLOAT *)o->var. The operators
 *   '+=', '-=', '*=' and '/=' are allowed.
 *
 * - For options without argument (VT_NONE), o->var is intepreted as
 *   (BOOLEAN *)o->var, the cmdline option '-name' sets it to TRUE,
 *   and '+name' sets it to FALSE.
 *
 * - For VT_STRING and VT_FILENAME, o->var is intepreted as
 *   (const char **)o->var, i.e., it points to the address of the string
 *   pointer, calling the option will modify the pointer to point to the new
 *   string which is allocated and freed by phgOptionXXXXX functions (thus user
 *   code should never free the string after calling phgOptionRegister).
 *
 *   Note: VT_STRING allows the '+=' operator while VT_FILENAME does not.
 *
 * - For VT_KEYWORD, o->keys is a NULL terminated list of valid keywords
 *   for this option, and o->var o->var is intepreted as (int *)o->var
 *   which is an index to the list of keywords.
 *
 * - For VT_HANDLER, (OPTION_HANDLER)o->var is the handler function for the
 *   option (which can be NULL), and the argument is stored in o->keys[0]
 *   through o->keys[n-1] (with o->keys[n] == NULL), depending on the value
 *   of o->append:
 *   	FALSE: n==1 and each call overrides the previous one
 *	TRUE:  all values are stored and each call increments n by 1
 *
 *   Note: after parsing all options, the value of the option '-oem_options'
 *   will be used as the cmdline arguments for external packages.
 *
 *   TODO: remove the member o->append, honor the '+=' operator in VT_HANDLER.
 */

typedef enum {
    VT_NONE, VT_INT, VT_FLOAT, VT_STRING, VT_FILENAME, VT_KEYWORD, VT_HANDLER,
    VT_INCLUDE, VT_TITLE, VT_INIT
} VTYPE;

typedef struct {
    char	*name;		/* option name without leading dash */
    char	*help;		/* help text for this option */
    char	**keys;		/* list of key words if type is VT_KEYWORD */
    void 	*var;		/* address of the variable to assign value to
				   It's assumed to have the following type:
					- VT_HANDLER	(OPTION_HANDLER)var
					- VT_FILENAME	(const char **)var
					- VT_STRING	(const char **)var
					- VT_KEYWORD	(int *)var
					- VT_INT	(INT *)var
					- VT_FLOAT	(FLOAT *)var
					- VT_NONE	(BOOLEAN *)var */
    VTYPE	type;		/* type of the variable */
    BOOLEAN	used;		/* whether the option is specified in cmdline */
    BOOLEAN	append;		/* for VT_HANDLER: append or replace */
} OPTION;

/* Option handler function protocol.
 *	'o' is the name of the option,
 *	'arg' is the argument with the following conventions:
 *		arg == NULL ==> print help
 *		arg[0] == '\0' ==> keep (and return) the current value
 *	'prefix' is the char prefixing '=' (e.g., '-myoption+=string'), or '\0'.
 * Return the string represeting the current value of the option (as a string
 * which can be used as the argument for the option), or NULL if error. */
typedef const char *(*OPTION_HANDLER)(OPTION *o, const char *arg, int prefix);

#ifdef __cplusplus
extern "C" {
#endif

#define phgOptionsRegisterInit_() \
    phgOptionsRegister(NULL, NULL, NULL, NULL, VT_INIT, FALSE)
#define phgOptionsRegisterNoArg(name, help, var) \
    phgOptionsRegister(name, help, NULL, var, VT_NONE, FALSE)
#define phgOptionsRegisterInt(name, help, var) \
    phgOptionsRegister(name, help, NULL, var, VT_INT, FALSE)
#define phgOptionsRegisterFloat(name, help, var) \
    phgOptionsRegister(name, help, NULL, var, VT_FLOAT, FALSE)
#define phgOptionsRegisterString(name, help, var) \
    phgOptionsRegister(name, help, NULL, var, VT_STRING, FALSE)
#define phgOptionsRegisterFilename(name, help, var) \
    phgOptionsRegister(name, help, NULL, var, VT_FILENAME, FALSE)
#define phgOptionsRegisterKeyword(name, help, keys, var) \
    phgOptionsRegister(name, help, keys, var, VT_KEYWORD, FALSE)
#define phgOptionsRegisterHandler(name, help, func, append) \
    phgOptionsRegister(name, help, NULL, func, VT_HANDLER, append)
#define phgOptionsRegisterTitle(str, help, category) \
    /* Note: category will be stored in '->var' */ \
    phgOptionsRegister(str, help, NULL, (void *)category, VT_TITLE, FALSE)

void phgOptionsRegister(const char *name, const char *help, const char **keys,
			void *var, VTYPE type, BOOLEAN append);
void phgOptionsModify(const char *name, const char *help, const char **keys,
		      void *var, VTYPE type, BOOLEAN append);
void phgOptionsReset(void);
void phgOptionsShowCmdline(void);
void phgOptionsShowUsed(void);
BOOLEAN phgOptionsIfUsed(const char *op_name);
void phgOptionsPrintfHelp(OPTION *o, const char *help);
void phgOptionsHelp(void);
void phgOptionsParseCmdline(int *argc, char ***argv);
void phgOptionsPreset(const char *str);

void phgOptionsSetOptions(const char *str);
BOOLEAN phgOptionsSetNoArg(const char *op_name, BOOLEAN value);
BOOLEAN phgOptionsSetInt(const char *op_name, INT value);
BOOLEAN phgOptionsSetFloat(const char *op_name, FLOAT value);
BOOLEAN phgOptionsSetKeyword(const char *op_name, const char *value);
BOOLEAN phgOptionsSetString(const char *op_name, const char *value);
BOOLEAN phgOptionsSetFilename(const char *op_name, const char *value);
BOOLEAN phgOptionsSetHandler(const char *op_name, const char *value);

void phgOptionsPush(void);
void phgOptionsPop(void);

BOOLEAN phgOptionsGetNoArg(const char *op_name);
INT phgOptionsGetInt(const char *op_name);
FLOAT phgOptionsGetFloat(const char *op_name);
const char *phgOptionsGetKeyword(const char *op_name);
const char *phgOptionsGetString(const char *op_name);
const char *phgOptionsGetFilename(const char *op_name);
const char *phgOptionsGetHandler(const char *op_name);

#ifdef __cplusplus
}
#endif

#define PHG_OPTION_H
#endif
