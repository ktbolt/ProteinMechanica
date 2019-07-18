
/* -------------------------------------------------------------------------- *
 *                              Protein Mechanica                             *
 * -------------------------------------------------------------------------- *
 * This is part of the Protein Mechanica coarse-grained molecular motor       *
 * modeling application originating from Simbios, the NIH National Center for *
 * Physics-Based Simulation of Biological Structures at Stanford, funded      *
 * under the NIH Roadmap for Medical Research, grant U54 GM072970.            *
 * See https://simtk.org.                                                     *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: David Parker                                                      *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

//*============================================================*
//* cmd:                 c o m m a n d                         *
//*============================================================*
// pm command processing.

#include "pm/pm.h"

const int CMD_MAX_LINE = 3000;
const int KEY_ESCAPE   =   27;

typedef struct CmdState {
  char c, line[CMD_MAX_LINE];
  int n, len;
  FILE *fp;
  } CmdState;

typedef struct CmdMacroList {
  char *str;
  struct CmdMacroList *next; 
  } CmdMacroList;

typedef struct CmdMacro {
  char *name;
  CmdMacroList *args;
  CmdMacroList *body;
  struct CmdMacro *next;
  } CmdMacro; 

#define CmdDataNull(dv) (!dv || dv[0] == '\0') 

typedef struct CmdKeyboardProcData {
  char cmd[100];
  char key[30];
  void (*func)(char*);
  } CmdKeyboardProcData; 


typedef struct CmdData {
  char name[80], val[3000];
  } CmdData;

typedef struct CmdDataList {
  int num;
  int count;
  CmdData list[20];
  } CmdDataList;

typedef struct Var {
  char *name;
  float val;
  } Var;

typedef struct StrVar {
  char *name;
  char val[3000];
  //char val[400];
  } StrVar;


typedef struct Cmd {
  char *name;
  void (*func)(char *cmd);
  } Cmd;

typedef struct Abrev {
  char *name;
  char *full_name;
  void (*func)(char *cmd);
  } Abrev;

void
pm_CmdDomain   (char *cmd),
pm_CmdDomains  (char *cmd),
pm_CmdMolecule (char *cmd),
pm_CmdProc     (char *cmd),
pm_CmdVar      (char *cmd),
pm_CmdHelp     (char *cmd);

static Cmd cmd_table[] = {
  {"domain",   pm_CmdDomain},
  {"domains",  pm_CmdDomains},
  {"molecule", pm_CmdMolecule},
  {"var",      pm_CmdVar},
  {"help",     pm_CmdHelp},
  {"?",        pm_CmdHelp},
  {NULL,       NULL}};

/*
static Abrev abrev_table[] = {
  {"da",     "dynamics motion animate",       pm_CmdDyn},
  {"dc",     "dynamics contact",              pm_CmdDyn},
  {"di",     "dynamics run incremental",      pm_CmdDyn},
  {"dm",     "dynamics motion",               pm_CmdDyn},
  {"dr",     "dynamics reset",                pm_CmdDyn},
  {"ds",     "dynamics state",                pm_CmdDyn},
  {"dv",     "dynamics variable",             pm_CmdDyn},
  {"dom",    "domain",                        pm_CmdDomain},
  {"gc",     "graphics center",               pm_CmdGraphics},
  {"gcp",    "graphics center pick",          pm_CmdGraphics},
  {"gp",     "graphics pick",                 pm_CmdGraphics},
  {"gr",     "graphics record",               pm_CmdGraphics},
  {"grot",   "graphics rotate",               pm_CmdGraphics},
  {"gw",     "graphics write",                pm_CmdGraphics},
  {"ma",     "molecule atom",                 pm_CmdMolecule},
  {"mc",     "molecule chain",                pm_CmdMolecule},
  {"mol",    "molecule",                      pm_CmdMolecule},
  {"msd",    "measure distance",              pm_CmdMeasure},
  {"n",      "node",                          pm_CmdNode},
  {"nm",     "modes",                         pm_CmdModes},
  {"nma",    "modes animate",                 pm_CmdModes},
  {"r",      "read",                          pm_CmdProc},
  {"s",      "dynamics increment",            pm_CmdDyn},
  {"reset",  "graphics reset",                pm_CmdGraphics},
  {NULL,       NULL}};
*/

static Abrev abrev_table[] = {
  {NULL,       NULL}};


char last_cmd[1000];

static Var g_vars[100];
static int g_num_vars = 0;

static StrVar g_str_vars[100];
static int g_num_str_vars = 0;

static int g_num_keyboard_proc = 0;
static CmdKeyboardProcData g_key_board_proc[20]; 

static CmdMacro *g_macros = NULL; 

static FILE *g_fp = NULL; 



//////////////////////////////////////////////////////////////////
//           u t i l i t y   p r o t o t y p e s               //
////////////////////////////////////////////////////////////////

int
pm_CmdBoolProc (char *s);

void
pm_CmdContact (CmdDataList *dlist);

void
pm_CmdContinueProc();

int
pm_CmdDataFloatGet (CmdDataList *dlist, char *name, float *val);

int
pm_CmdDataFvecGet (CmdDataList *dlist, char *name, float val[3]);

void
pm_CmdDataFvecParse (char *s,  float val[3]);

int
pm_CmdDataNext (CmdDataList *dlist, char **name, char **val);

void
pm_CmdDataParse (char *cmd, CmdDataList *dlist);

int
pm_CmdDataValGet (CmdDataList *dlist, char *name, char **val);

void
pm_CmdFloatEval (char *s, float *val);

void
pm_CmdIntEval (char *s, int *val);

void
pm_CmdInterruptProc ();

void
pm_CmdKeyboardProc (char *cmd);

void
pm_CmdMacroDef (char *cmd);

void
pm_CmdMacroProc (CmdMacro *macro, char *cmd);

void
pm_CmdStrVarGet (char *name, StrVar **var);

void
pm_CmdStrVarValGet (char *name, char **val);

void
pm_CmdVarGet (char *name, Var **var);

