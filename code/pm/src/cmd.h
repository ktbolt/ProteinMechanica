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
#include "model.h"
#include "grid.h"
#include "part.h"
#include "gc.h"
#include "bexp.h"
#include <sstream>
#include <fstream>
#include <string.h>

namespace ProteinMechanica {

const char NC = '\0';
const char COMMENTCHAR = '#';
const char VARIABLECHAR = '$';
const char CMDPR[] = "    >>> ";
const char CMDSP[] = "        ";

typedef enum PmCmdSymbolType {
  PM_SYMBOL_UNKNOWN,
  PM_SYMBOL_ALPHA,
  PM_SYMBOL_AND,
  PM_SYMBOL_BACKSLASH,
  PM_SYMBOL_BLANK,
  PM_SYMBOL_BRACE,
  PM_SYMBOL_BRACKET,
  PM_SYMBOL_COLON,
  PM_SYMBOL_COMMA,
  PM_SYMBOL_DOLLAR,
  PM_SYMBOL_DIGIT,
  PM_SYMBOL_EQUAL,
  PM_SYMBOL_GT,
  PM_SYMBOL_LT,
  PM_SYMBOL_OPERATOR,
  PM_SYMBOL_OR,
  PM_SYMBOL_PARENTHESES,
  PM_SYMBOL_PERIOD,
  PM_SYMBOL_QUOTE,
  PM_SYMBOL_SIGN,
  PM_SYMBOL_SLASH,
  PM_SYMBOL_STAR 
  } PmCmdSymbolType;

typedef enum PmTokenType {
  PM_TOKEN_UNKNOWN,
  PM_TOKEN_BRACE,
  PM_TOKEN_BRACKET,
  PM_TOKEN_COMMENT,
  PM_TOKEN_CONTINUE,
  PM_TOKEN_EQUAL,
  PM_TOKEN_NUMBER,
  PM_TOKEN_PARENTHESES,
  PM_TOKEN_SEPARATOR,
  PM_TOKEN_STRING,
  PM_TOKEN_VARIABLE
  } PmTokenType;

typedef enum PmDataType {
  PM_DATA_UNKNOWN,
  PM_DATA_NONE,
  PM_DATA_NUMBER,
  PM_DATA_NUMBER_LIST,
  PM_DATA_STRING,
  PM_DATA_STRING_LIST,
  PM_DATA_VARIABLE
  } PmDataType;

const int CMD_MAX_LINE = 3000;
const int KEY_ESCAPE   =   27;
const char *delim = " ";

typedef struct CmdState {
  char c, line[CMD_MAX_LINE];
  int n, len;
  FILE *fp;
  } CmdState;

typedef struct CmdMacro {
  string name;
  vector <string> args;
  vector <string> body;
  } CmdMacro;

static vector<CmdMacro*> g_macros;

//-------------------------------------------------------------
//                    PmBooleanExpression                     |
//-------------------------------------------------------------

static vector<PmBooleanExpression*> g_boolean_expressions;

//-------------------------------------------------------------
//                    CmdKeyboardProcData                     |
//-------------------------------------------------------------

typedef struct CmdKeyboardProcData {
  string cmd;
  string key;
  void (*func)(char*);
  } CmdKeyboardProcData;

char last_cmd[1000];
static vector<CmdKeyboardProcData> g_key_board_proc;

//-------------------------------------------------------------
//                            CmdToken                        |
//-------------------------------------------------------------

class CmdToken {
  public:
    string val;
    PmTokenType type;
    bool closed;
    bool error;
    bool eol;
    CmdToken() { init(); }
    void init() { val.clear(); type = PM_TOKEN_UNKNOWN;
                  closed = false; error = false; eol = false; }
  };

//-------------------------------------------------------------
//                            PmCmdData                       |
//-------------------------------------------------------------

class PmCmdData {
  public:
    string name;
    PmDataType type;
    bool read; 
    vector<float> numbers;
    vector<string> strings;
    PmCmdData(){init();};
    void init() { 
      type = PM_DATA_UNKNOWN; 
      name.clear(); 
      read = false;
      numbers.clear(); 
      strings.clear(); 
      }
    bool getBoolean();
    float getFloat();
    int getInt();
    void getString(string& s);
    bool getVector(PmVector3& vec);
    bool getColor(PmVector3& vec);
    void getStringList(vector<string>& slist);
  };

//-------------------------------------------------------------
//                            PmCmdDataList                   |
//-------------------------------------------------------------

class PmCmdDataList {
  public:
    vector<PmCmdData> list;
    bool error;
    int count; 
    PmCmdDataList() { init(); }; 
    void init() { error = false; list.clear(); count = 0; }
    bool check();
    bool getBoolean(const string name, bool& val);
    bool getNext(PmCmdData& data, string& val);  
    bool getNext(PmCmdData& data);  
    bool getData(const string name, PmCmdData& data);
    bool getFloat(const string name, float& val);
    bool getInt(const string name, int& val);
    bool getString(const string name, string& s);
    bool getStringList(const string name, vector<string>& slist);
    bool getVector(const string name, PmVector3& vec);
    bool getColor(const string name, PmVector3& vec);
  };

//-------------------------------------------------------------
//                            PmVariable                      |
//-------------------------------------------------------------

typedef struct PmVariable {
  string name;
  string value;
  PmCmdData data;
  PmBooleanExpression *bexp;
  bool persistent;
  } PmVariable;

static vector<PmVariable> g_vars;


//-------------------------------------------------------------
//                            PmCmd                           |
//-------------------------------------------------------------

class PM_EXPORT PmCmd {
  public:
    PmCmd(); 
    void setLine (string line);
    void getDataList (PmCmdDataList& dlist) ;
    void printDataList (PmCmdDataList& dlist);
    bool getChar();
    bool getToken(CmdToken& token);
    bool isDelimiter(char c);
    PmCmdSymbolType classifyChar(char c);
    bool checkChar();
    void callMacro(FILE *fp, CmdMacro *macro, PmCmdDataList& dlist);
    void procArgList(CmdToken& token, PmCmdDataList& dlist);
    char getCurrChar( ) { return curr_c; };
    int getCurrCnum( ) { return cnum; };
    static void testScript (const string file_name);
    static void processScript (const string file_name);

  private:
     string line;
     int line_size;
     char curr_c, next_c;
     int cnum;
     bool in_brace;
     vector<PmVariable> variables;
  };


//-------------------------------------------------------------
//           command function tables                          |
//-------------------------------------------------------------

typedef struct Cmd {
  char *name;
  void (*func)(PmCmdDataList& dlist);
  } Cmd;

typedef struct Abbrev {
  char *name;
  char *full_name;
  void (*func)(PmCmdDataList& dlist);
  } Abbrev;

void
pm_CmdBindSimulation (PmCmdDataList& dlist),
pm_CmdBody           (PmCmdDataList& dlist),
pm_CmdBodies         (PmCmdDataList& dlist),
pm_CmdCurve          (PmCmdDataList& dlist),
pm_CmdDatabase       (PmCmdDataList& dlist),
pm_CmdDomain         (PmCmdDataList& dlist),
pm_CmdDomains        (PmCmdDataList& dlist),
pm_CmdForce          (PmCmdDataList& dlist),
pm_CmdGraphics       (PmCmdDataList& dlist),
pm_CmdGrid           (PmCmdDataList& dlist),
pm_CmdJoint          (PmCmdDataList& dlist),
pm_CmdModel          (PmCmdDataList& dlist),
pm_CmdMeasure        (PmCmdDataList& dlist),
pm_CmdMeasurement    (PmCmdDataList& dlist),
pm_CmdMolecule       (PmCmdDataList& dlist),
pm_CmdMotor          (PmCmdDataList& dlist),
pm_CmdMultiBody      (PmCmdDataList& dlist),
pm_CmdParticle       (PmCmdDataList& dlist),
pm_CmdPotential      (PmCmdDataList& dlist),
pm_CmdSimulation     (PmCmdDataList& dlist),
pm_CmdSimulations    (PmCmdDataList& dlist),
pm_CmdSolid          (PmCmdDataList& dlist),
pm_CmdSurface        (PmCmdDataList& dlist),
pm_CmdSystem         (PmCmdDataList& dlist),
pm_CmdTrace          (PmCmdDataList& dlist),
pm_CmdProc           (PmCmdDataList& dlist),
pm_CmdUnits          (PmCmdDataList& dlist),
pm_CmdVector         (PmCmdDataList& dlist),
pm_CmdHelp           (PmCmdDataList& dlist);

// command table //

static Cmd cmd_table[] = {
  {"binding_simulation",   pm_CmdBindSimulation},
  {"body",                 pm_CmdBody},
  {"bodies",               pm_CmdBodies},
  {"curve",                pm_CmdCurve},
  {"database",             pm_CmdDatabase},
  {"domain",               pm_CmdDomain},
  {"domains",              pm_CmdDomains},
  {"force",                pm_CmdForce},
  {"graphics",             pm_CmdGraphics},
  {"grid",                 pm_CmdGrid},
  {"joint",                pm_CmdJoint},
  {"measure",              pm_CmdMeasure},
  {"measurement",          pm_CmdMeasurement},
  {"model",                pm_CmdModel},
  {"molecule",             pm_CmdMolecule},
  {"motor",                pm_CmdMotor},
  {"multibody",            pm_CmdMultiBody},
  {"particle",             pm_CmdParticle},
  {"potential",            pm_CmdPotential},
  {"simulation",           pm_CmdSimulation},
  {"simulations",          pm_CmdSimulations},
  {"solid",                pm_CmdSolid},
  {"surface",              pm_CmdSurface},
  {"system",               pm_CmdSystem},
  {"trace",                pm_CmdTrace},
  {"units",                pm_CmdUnits},
  {"vector",               pm_CmdVector},
  {"help",                 pm_CmdHelp},
  {"?",                    pm_CmdHelp},
  {NULL,                   NULL}};

// command abreviation table //

static Abbrev abbrev_table[] = {
  {"gcp",      "graphics center pick", pm_CmdGraphics},
  {"gr",       "graphics record",      pm_CmdGraphics},
  {"gw",       "graphics write",       pm_CmdGraphics},
  {"md",       "measure distance pick=true", pm_CmdMeasure},
  {"s",        "simulation step",      pm_CmdSimulation},
  {"sa",       "simulations step",     pm_CmdSimulations},
  {"step",     "simulation step",      pm_CmdSimulation},
  {"sr",       "simulation replay",    pm_CmdSimulation},
  {NULL,       NULL, NULL}};

// system variables //

static const string system_vars[] = {
   "potential_energy",
   "step",
   ""
   };

void pm_CmdQuit();
void pm_CmdAddVariable(const string line);
void pm_CmdAddVariable(const string name, const string value);
bool pm_CmdGetVariable(const string name, PmVariable **var);
bool pm_CmdGetVariableData(const string name, PmCmdData& data);
bool pm_CmdGetVariableValue(PmCmd *cmd, int cnum, const string sym, int& n, 
                             string& value);
bool pm_CmdGetVariableValue(const string name, std::string& value);
bool pm_CmdHasVariable(const string name);
bool pm_CmdSetVariableValue(const string name, const string value); 

bool pm_CmdGetBooleanExpression(const string name, PmBooleanExpression **bexp);

void pm_CmdLineProcess (FILE *fp, const string line);
void PM_EXPORT pm_CmdProcess (FILE *fp, const string line);
bool pm_CmdSubstVariable (FILE *fp, string& line);

}
