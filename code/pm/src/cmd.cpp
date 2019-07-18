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

#include "bsim.h"
#include "cmd.h"
#include "db.h"
#include "msr.h"
#include "pm/mth.h"
#include "part.h"
#include "pot.h"
#include "rbsim.h"
#include "rest.h"
#include "sim.h"
#include "solid.h"
#include "trace.h"
#include "geom.h"

namespace ProteinMechanica {

#define ndbg_getToken
#define ndbg_getDataList

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmCmd::PmCmd() 
  {
  curr_c = NC;
  next_c = NC;
  cnum = 0;
  in_brace = false;
  }

//*============================================================*
//*==========              convToFloat               ==========*
//*============================================================*
// convert a string into a float.

bool
convToFloat (string s, float& val)
  {
  //fprintf (stderr,"convToFloat  %s  \n", s.c_str());

  if (s == "inf") {
    val = std::numeric_limits<float>::max();
    return true;
    }

  std::istringstream iss(s);
  char c;

  if (!(iss >> val) || iss.get(c)) {
    fprintf (stderr,"**** Error: convToFloat  of \"%s\" failed. \n", s.c_str());
    return false;
    }

  return true;
  }

//*============================================================*
//*==========              setLine                   ==========*
//*============================================================*
// set the line for command processing.

void
PmCmd::setLine (string line) 
  {
  this->line = line;
  line_size = line.size();
  curr_c = NC;
  next_c = NC;
  cnum = 0;
  }

//*============================================================*
//*==========              isDelimiter               ==========*
//*============================================================*
// check if a character is a delimiter.

bool 
PmCmd::isDelimiter(char c) 
  {
  if ((c == '=') ||
      (c == '[') ||
      (c == ']') ||
      (c == '{') ||
      (c == '}') ||
      (c == ',') ||
      (c == ' ') ||
      (c == '\"')) {
    return true;     
    }

  return false;
  }

//*============================================================*
//*==========              getChar                   ==========*
//*============================================================*
// get the next character from line[].

bool
PmCmd::getChar() 
  {
  int n;
  static string var_value;
  static int var_n = 0;
  static int var_size = 0;

  if (var_size) {
    curr_c = var_value[var_n];

    if (var_n < var_size-1) {
      next_c = var_value[var_n+1];
      }
    else if (cnum < line_size) {
      next_c = line[cnum];
      // NOTE: changed this Fri Nov  6 14:55:02 PST 2009
      //cnum += 1;
      }
    else {
      next_c = NC;
      }

    var_n += 1;

    if (var_n == var_size) {
      var_size = 0;
      var_n = 0;
      var_value.clear();
      }

    //fprintf (stderr, ">>>> var: curr %c  next %c  \n", curr_c, next_c);
    return true;
    }

  if (cnum == line_size) {
    curr_c = NC;
    return false;
    }

  curr_c = line[cnum];
  //fprintf (stderr, ">>>> curr_c = %c  \n", curr_c);

  #define nproc_var_PmCmd_getChar
  #ifdef proc_var_PmCmd_getChar

  if (curr_c == '$') {
    if (!pm_CmdGetVarValue(this, cnum, line, n, var_value)) { 
      //fprintf (stderr, "**** Error: bad variable name at %d. \n", cnum);
      curr_c = NC;
      next_c = NC;
      cnum = line_size;
      var_value.clear();
      return false;
      }

    cnum += n;
    var_n = 0;
    var_size = var_value.size();
    curr_c = var_value[var_n];

    fprintf (stderr, ">>>>>> PmCmd::getChar: proc variable \n");
    fprintf (stderr, ">>> var_value=\"%s\" n=%d cnum=%d \n", var_value.c_str(), n, cnum);
    //fprintf (stderr, ">>>> var size = %d  \n", var_size);
    //fprintf (stderr, ">>>> var sym size = %d  \n", n);
    //fprintf (stderr, ">>>> var line[%d] = %c  \n", cnum, line[cnum]);

    if (var_size > 1) {
      next_c = var_value[var_n+1];
      var_n += 1;
      return true;
      }
    else {
      var_size = 0;
      var_n = 0;
      var_value.clear();
      cnum -= 1;
      }
    } 
  #endif

  //fprintf (stderr, ">>>> line[%d] = %c  \n", cnum, line[cnum]);

  if (cnum < line_size-1) {
    next_c = line[cnum+1];
    }
  else {
    next_c = NC;
    }

  //fprintf (stderr, ">>>> next_c = %c  \n", next_c);
  cnum += 1;
  return true;
  }

//*============================================================*
//*==========              classifyChar              ==========*
//*============================================================*
// classify a character. only classify those characters used
// in pm. all others will be classified as unknown. 

PmCmdSymbolType
PmCmd::classifyChar(char c)
  {
  PmCmdSymbolType sym = PM_SYMBOL_UNKNOWN;

  if (isalpha(c) || (c == '_')) {
    sym = PM_SYMBOL_ALPHA;
    }
  else if (isdigit(c)) {
    sym = PM_SYMBOL_DIGIT;
    }
  else if ((c == '-') || (c == '+')) {
    sym = PM_SYMBOL_SIGN;
    }
  else if ((c == '[') || (c == ']')) { 
    sym = PM_SYMBOL_BRACKET;
    }
  else if ((c == '{') || (c == '}')) {
    sym = PM_SYMBOL_BRACE;
    }
  else if (c == ',') {
    sym = PM_SYMBOL_COMMA;
    }
  else if (c == '.') {
    sym = PM_SYMBOL_PERIOD;
    }
  else if ((c == '(') || (c == ')')) {
    sym = PM_SYMBOL_PARENTHESES;
    }
  else if (c == '=') {
    sym = PM_SYMBOL_EQUAL;
    }
  else if (c == ':') {
    sym = PM_SYMBOL_COLON;
    }
  else if ((c == ' ') || (c == '\0')) {
    sym = PM_SYMBOL_BLANK;
    }
  else if (c == '/') {
    sym = PM_SYMBOL_SLASH;
    }
  else if (c == '\\') {
    sym = PM_SYMBOL_BACKSLASH;
    }
  else if ((c == '\'') || (c == '\"')) {
    sym = PM_SYMBOL_QUOTE;
    }
  else if (c == '$') {
    sym = PM_SYMBOL_DOLLAR;
    }
  else if (c == '<') {
    sym = PM_SYMBOL_LT;
    }
  else if (c == '>') {
    sym = PM_SYMBOL_GT;
    }
  else if (c == '&') {
    sym = PM_SYMBOL_AND;
    }
  else if (c == '|') {
    sym = PM_SYMBOL_OR;
    }
  else if (c == '*') {
    sym = PM_SYMBOL_STAR;
    }

  return sym;
  }

//*============================================================*
//*==========              checkChar                 ==========*
//*============================================================*
// check for valid combinations of characters.

bool
PmCmd::checkChar()
  {
  bool status;
  PmCmdSymbolType csym = classifyChar(curr_c);
  PmCmdSymbolType nsym = classifyChar(next_c);
  //fprintf (stderr,"[%c:%d]  [%c:%d]", curr_c, curr_c, next_c, next_c);
  //fprintf (stderr,"[%c:%d_%c:%d]\n", curr_c, csym, next_c, nsym);
  if (nsym == PM_SYMBOL_BLANK) return true;

  if ((csym == PM_SYMBOL_UNKNOWN) || (nsym == PM_SYMBOL_UNKNOWN)) {
    return false;
    }

  else if (csym == PM_SYMBOL_SIGN) {
    if ( (nsym != PM_SYMBOL_DIGIT) && 
         (nsym != PM_SYMBOL_ALPHA)  &&
         (nsym != PM_SYMBOL_PERIOD) ) { 
      fprintf (stderr,"\n>>>>>> PmCmd::checkChar failed 1\n");
      return false;
      }
    }

  // certain combinations are legal only within a brace //
  // where a domain descriptor may be used.             //

  else if ((csym == PM_SYMBOL_BRACKET) || (csym == PM_SYMBOL_BRACE)) {
    if ((curr_c == ']') && in_brace) {
      return true;
      }

    if ((curr_c == '}') && in_brace) {
      return true;
      }

    if ( (nsym != PM_SYMBOL_DIGIT)        && 
         (nsym != PM_SYMBOL_PERIOD)       &&
         (nsym != PM_SYMBOL_ALPHA)        &&
         (nsym != PM_SYMBOL_COLON)        &&
         (nsym != PM_SYMBOL_COMMA)        &&
         (nsym != PM_SYMBOL_PARENTHESES)  &&
         (nsym != PM_SYMBOL_DOLLAR)       && 
         (nsym != PM_SYMBOL_BRACKET)      &&
         (nsym != PM_SYMBOL_BRACE)        &&
         (nsym != PM_SYMBOL_SIGN) ) { 
      fprintf (stderr,"\n>>>>>> PmCmd::checkChar failed 2\n");
      return false;
      }
    }

  else if (csym == PM_SYMBOL_EQUAL) {
    if ( (nsym != PM_SYMBOL_DIGIT)     &&
         (nsym != PM_SYMBOL_PERIOD)    &&
         (nsym != PM_SYMBOL_ALPHA)     &&
         (nsym != PM_SYMBOL_BRACKET)   &&
         (nsym != PM_SYMBOL_BRACE)     &&
         (nsym != PM_SYMBOL_DOLLAR)    &&
         (nsym != PM_SYMBOL_SLASH)     &&
         (nsym != PM_SYMBOL_BACKSLASH) &&
         (nsym != PM_SYMBOL_SIGN) ) {
      fprintf (stderr,"\n>>>>>> PmCmd::checkChar failed 3\n");
      return false;
      }
    }

  else if (csym == PM_SYMBOL_COMMA) {
    if ( (nsym != PM_SYMBOL_DIGIT)   &&
         (nsym != PM_SYMBOL_PERIOD)  &&
         (nsym != PM_SYMBOL_ALPHA)   &&
         (nsym != PM_SYMBOL_SIGN) ) {
      fprintf (stderr,"\n>>>>>> PmCmd::checkChar failed 4\n");
      return false;
      }
    }

  else if (csym == PM_SYMBOL_DOLLAR) {
    //if ((nsym != PM_SYMBOL_ALPHA) && (nsym != PM_SYMBOL_BRACKET)) { 
    if (nsym != PM_SYMBOL_BRACE) { 
      fprintf (stderr,"\n>>>>>> PmCmd::checkChar failed 5\n");
      return false;
      }
    }

  else if ((csym == PM_SYMBOL_SLASH) || (csym == PM_SYMBOL_BACKSLASH)) {
    if ((nsym != PM_SYMBOL_ALPHA) && (nsym != PM_SYMBOL_PERIOD) && 
        (nsym != PM_SYMBOL_DIGIT)) { 
      fprintf (stderr,"\n>>>>>> PmCmd::checkChar failed 6\n");
      fprintf (stderr,">>> curr_c='%c'  next_c='%c' \n", curr_c, next_c);
      return false;
      }
    }

  return true;
  }

//*============================================================*
//*==========              getToken                  ==========*
//*============================================================*
// get the next token from line[].

bool
PmCmd::getToken(CmdToken& token)
  {
  token.val.clear();
  token.type = PM_TOKEN_UNKNOWN;
  token.closed = false;
  bool first = true;
  int num_q = 0;
  #ifdef dbg_getToken
  fprintf (stderr,"\n token: ");
  #endif

  //fprintf (stderr,"\n------ get token ------  \n");

  while (getChar()) {
    //fprintf (stderr,">>> curr=%c next=%c type=%d \n", curr_c, next_c, token.type); 

    if (!checkChar()) {
      token.error = true; 
      //fprintf (stderr,"################ \n");
      //fprintf (stderr,"curr %c  next %c \n", curr_c, next_c); 
      return false;
      }

    #ifdef dbg_getToken
    fprintf (stderr,"_%c_[%d]", curr_c, token.type);
    #endif

    if (curr_c == COMMENTCHAR) { 
      token.type = PM_TOKEN_COMMENT;
      token.val.push_back(curr_c);
      break;
      }

    else if (curr_c == '\\') { 
      token.type = PM_TOKEN_CONTINUE;
      token.val.push_back(curr_c);
      break;
      }

    else if (curr_c == '(') {
      token.type = PM_TOKEN_PARENTHESES;
      break;
      }

    else if (curr_c == '$') {
      if (next_c == '{') {
        token.type = PM_TOKEN_VARIABLE;
        token.val.push_back(curr_c);
        in_brace = true;
        //fprintf (stderr,">>> PmCmd::getToken start token variable \n");
        }
      }

    else if (token.type == PM_TOKEN_VARIABLE) {
      token.val.push_back(curr_c);

      if (curr_c == '}') {
        token.closed = true;
        in_brace = false;
        //fprintf (stderr,">>> PmCmd::getToken define token variable=%s \n", token.val.c_str());
        }
      }

    // most symbols will go into a string if a  //
    // a string is being processed.             //

    else if (token.type == PM_TOKEN_STRING) {
      token.val.push_back(curr_c);
      //fprintf (stderr,"_%c_[%d]", curr_c, token.type);
      }

    else if ((curr_c == '{') || (curr_c == '}')) {
      token.type = PM_TOKEN_BRACE;
      token.val.push_back(curr_c);

      if (curr_c == '}') {
        token.closed = true;
        in_brace = false;
        }
      else {
        in_brace = true;
        }
      break;
      }

    // brackets are used for vectors and specifying  //
    // a domain so check to see if we are processing //
    // a string.                                     //

    else if ((curr_c == ']') || (curr_c == '[')) { 
      if (token.type != PM_TOKEN_STRING) {
        token.type = PM_TOKEN_BRACKET;
        token.val.push_back(curr_c);
        if (curr_c == ']') token.closed = true;
        break;
        }
      else {
        token.val.push_back(curr_c);
        }
      }

    else if (curr_c == ',') {
      token.val.push_back(curr_c);

      if (token.type != PM_TOKEN_STRING) {
        token.type = PM_TOKEN_SEPARATOR;
        break;
        }
      }

    else if (curr_c == NC) {
      token.eol = true;
      break;
      }

    // for certain commands we want to include  //
    // '=' in a string.                         //

    else if ((curr_c == '=') && !in_brace) {
      token.val.push_back(curr_c);
      token.type = PM_TOKEN_EQUAL;
      break;
      }

    else if (curr_c == '\"') {
      if (first) {
        first = false;
        token.type = PM_TOKEN_STRING;
        }

      num_q += 1;

      if (num_q == 2) {
        //fprintf (stderr, ">>>> token=[%s] \n", token.val.c_str()); 
        return true;
        }
      }

    else if (curr_c == '/') {
      token.type = PM_TOKEN_STRING;
      token.val.push_back(curr_c);
      }

    else if (curr_c == '.') {
      if (next_c == '.') {
        first = false;
        token.type = PM_TOKEN_STRING;
        token.val.push_back(curr_c);
        }
      else {
        token.type = PM_TOKEN_NUMBER;
        token.val.push_back(curr_c);
        }
      }

    else if (isalpha(curr_c)) {
      token.val.push_back(curr_c);

      if (first && (next_c != '-') && (next_c != '+')) {
        first = false;
        token.type = PM_TOKEN_STRING;
        }
      }

    else if (isdigit(curr_c)) {
      if ((next_c == '[') || in_brace) {
        token.type = PM_TOKEN_STRING;
        token.val.push_back(curr_c);
        }
      else {
        token.type = PM_TOKEN_NUMBER;
        token.val.push_back(curr_c);
        }
      }

    else if ((curr_c == '-') || (curr_c == '+')) {
      if (isdigit(next_c) || (next_c == '.')) {
        token.val.push_back(curr_c);
        token.type = PM_TOKEN_NUMBER;
        }
      }

    // if a delimiter is encountered then finish //

    if (isDelimiter(next_c)) {
      if ((next_c == '[') || (next_c == ']') && (token.type == PM_TOKEN_STRING)) {
        continue;
        }

      if ((next_c == '=') && (token.type == PM_TOKEN_STRING) && in_brace) {
        continue;
        }

      if ((next_c == ',') && (token.type == PM_TOKEN_STRING)) {
        continue;
        }

      if (token.type == PM_TOKEN_VARIABLE) {
        continue;
        }

      if (curr_c != ' ') break;
      } 
    }

  if (curr_c == NC) token.eol = true;
  return true;
  }

//*============================================================*
//*==========              getDataList               ==========*
//*============================================================*
// get a list of data items from line[].

void
PmCmd::getDataList (PmCmdDataList& dlist) 
  {
  //fprintf (stderr, "\n>>>>>> PmCmd::getDataList line=%s \n", line.c_str());
  CmdToken prev_token, token;
  PmCmdData data; 
  string val;
  dlist.init();
  dlist.error = true;

  if (!getToken(prev_token)) return;
  #ifdef dbg_getDataList
  fprintf (stderr, "prev_token[%s] type[%d]\n", prev_token.val.c_str(), prev_token.type);
  #endif

  // process argument list. //
  // NOTE [20Nov2009] This is in two places depending on form of the macro call. //
  // Need to figure out this some day.                                           //

  if (prev_token.type == PM_TOKEN_PARENTHESES) {
    procArgList(prev_token, dlist);
    //fprintf (stderr, "\n####### here 1 \n");
    dlist.error = false;
    return;
    }

  while (1) {
    if (!getToken(token)) {
      if (dlist.error) {
        fprintf (stderr, "**** Error parsing data list. \n");
        fprintf (stderr, "**** Error at [%s] \n", token.val.c_str());
        fprintf (stderr, "     token=\"%s\" type=%d eol=%d\n", 
                 token.val.c_str(), token.type,token.eol);
        fprintf (stderr, "     previous token=\"%s\" \n", 
                 prev_token.val.c_str());
        }
      return;
      }

    #ifdef dbg_getDataList
    fprintf (stderr, "token[%s] type[%d] eol[%d]\n", token.val.c_str(), 
             token.type,token.eol);
    #endif

    //===== process value definition:  <name> = <value> =====//

    if (token.type == PM_TOKEN_EQUAL) {
      data.name = prev_token.val; 

      if (!getToken(token)) {
        if (dlist.error) {
          fprintf (stderr, "**** Error parsing equal symbol. \n");
          fprintf (stderr, "**** Error at [%s] \n", token.val.c_str());
          }
        return;
        }

      #ifdef dbg_getDataList
      fprintf (stderr, ">>> token=%s type=%d\n", token.val.c_str(), token.type);
      #endif

      if (token.type == PM_TOKEN_STRING) {
        data.type = PM_DATA_STRING;
        data.strings.push_back(token.val);  
        dlist.list.push_back(data);
        data.init(); 
        }
      else if (token.type == PM_TOKEN_VARIABLE) {
        data.type = PM_DATA_STRING;
        data.strings.push_back(token.val);  
        dlist.list.push_back(data);
        data.init(); 
        }
      else if (token.type == PM_TOKEN_NUMBER) {
        data.type = PM_DATA_NUMBER;
        float fval; 

        if (!convToFloat(token.val, fval)) {
          fprintf (stderr, "**** Error parsing float. \n");
          fprintf (stderr, "**** Error at [%s] \n", token.val.c_str());
          return;
          }

        data.numbers.push_back(fval);  
        dlist.list.push_back(data);
        data.init(); 
        }

      // process a list:  <name> = { <val1> <val2> ... <valn> }  //

      else if (token.type == PM_TOKEN_BRACE) {
        while (1) {
          if (!getToken(token)) {
            if (dlist.error) {
              fprintf (stderr, "**** Error in \"{\" at [%s] \n", token.val.c_str());
              }
            return;
            }

          #ifdef dbg_getDataList
          fprintf (stderr, "token[%s] type[%d] \n", token.val.c_str(), token.type);
          #endif

          if ((token.type == PM_TOKEN_STRING) || (token.type == PM_TOKEN_NUMBER)) {
            data.type = PM_DATA_STRING_LIST;
            data.strings.push_back(token.val);
            }

          else if (token.type == PM_TOKEN_BRACE) {
            dlist.list.push_back(data);
            data.init();
            prev_token.init();
            break;
            }
          else if (token.val[0] != ',') {
            fprintf (stderr, "**** Error parsing brace. \n");
            fprintf (stderr, "**** Error at [%s] \n", token.val.c_str());
            return;
            }
          }
        }

      // process a vector:  <name> = [ <val1> <val2> <val3> ]  //

      else if (token.type == PM_TOKEN_BRACKET) {
        int num = 0;

        while (1) {
          if (!getToken(token)) return;
          #ifdef dbg_getDataList
          fprintf (stderr, "token[%s] type[%d] \n", token.val.c_str(), token.type);
          #endif

          if (token.type == PM_TOKEN_NUMBER) {
            num += 1;
            data.type = PM_DATA_NUMBER_LIST;
            float fval; 

            if (num > 3) {
              fprintf (stderr, "**** Error: to many values for a vector.\n");
              return;
              }

            if (!convToFloat(token.val, fval)) {
              fprintf (stderr, "**** Error parsing number. \n");
              fprintf (stderr, "**** Error at [%s] \n", token.val.c_str());
              return;
              }

            data.numbers.push_back(fval);  
            }

          else if (token.type == PM_TOKEN_BRACKET) {
            dlist.list.push_back(data);
            data.init(); 
            prev_token.init();

            if (num != 3) {
              fprintf (stderr, "**** Error: need three values for a vector.\n");
              return;
              }

            break;
            }
          else if (token.val[0] != ',') {
            fprintf (stderr, "**** Error parsing bracket. \n");
            fprintf (stderr, "**** Error at [%s] \n", token.val.c_str());
            return;
            }
          }
        }

      prev_token.init();
      }

    // lone string //

    else if ((token.type == PM_TOKEN_STRING) || 
             (token.type == PM_TOKEN_NUMBER) || token.eol) {
      if (!prev_token.val.empty()) {
        data.name = prev_token.val; 
        data.type = PM_DATA_NONE;
        //data.strings.push_back(prev_token.val);  
        dlist.list.push_back(data);
        data.init(); 
        #ifdef dbg_getDataList
        fprintf (stderr, "add lone string[%s]\n", prev_token.val.c_str());
        #endif
        }

      if (token.eol && !token.val.empty()) {
        data.name = token.val; 
        data.type = PM_DATA_NONE;
        // NOTE: not sure about this //
        data.strings.push_back(token.val);  
        dlist.list.push_back(data);
        }

      prev_token = token;
      }

    // process an argument list //

    else if (token.type == PM_TOKEN_PARENTHESES) {
      if (prev_token.type == PM_TOKEN_EQUAL) {
        break;
        }
      else {
        procArgList(prev_token, dlist);
        //fprintf (stderr, "\n####### here 2 \n");
        break;
        }
      }

    if (token.eol) break;
    }

  dlist.error = false;
  }

//*============================================================*
//*==========              procArgList               ==========*
//*============================================================*
// process an argument list. 

void
PmCmd::procArgList (CmdToken& token, PmCmdDataList& dlist)
  {
  #define ndbg_PmCmd_procArgList 
  #ifdef dbg_PmCmd_procArgList 
  fprintf (stderr, "\n--------- proc arg list ----------- \n");
  fprintf (stderr, ">>> token val = %s \n", token.val.c_str());
  #endif

  string arg;
  PmCmdData data; 
  int n;
  data.name = token.val; 
  data.type = PM_DATA_NONE;
  dlist.list.push_back(data);
  n = 0;

  while (getChar()) {
    if ((curr_c == ',') || (curr_c == ')')) {
      n = arg.size(); 
      data.name.clear(); 

      for (int i = n-1; i >= 0; i--) { 
        if (arg[i] != ' ') {
          n = i + 1;
          break;
          }
        }

      for (int i = 0; i < n; i++) { 
        data.name.push_back(arg[i]); 
        }

      data.type = PM_DATA_NONE;
      dlist.list.push_back(data);
      #ifdef dbg_PmCmd_procArgList 
      fprintf (stderr, ">>> arg=%s  data.name=%s \n", arg.c_str(), data.name.c_str());
      #endif
      arg.clear();
      n = 0;
      }

    else { 
      if ((curr_c == ' ') && (n == 1)) {
        arg.push_back(curr_c);
        }
      else if (curr_c != ' ') { 
        arg.push_back(curr_c);
        n = 1;
        }
      }
    }
  }

//*============================================================*
//*==========              printDataList             ==========*
//*============================================================*
// print data list items.

void
PmCmd::printDataList (PmCmdDataList& dlist) 
  {
  fprintf (stderr, "\n--------- data list ----------- \n");
  int n = dlist.list.size();
  fprintf (stderr, ">>>>>> n [%d] \n", n);

  for (int i = 0; i < n; i++) {
    if (dlist.list[i].type == PM_DATA_NONE) {
      fprintf (stderr, "%d: name [%s] \n", i+1, dlist.list[i].name.c_str());
      }
    else {
      fprintf (stderr, "%d: name [%s] ", i+1, dlist.list[i].name.c_str());

      if (dlist.list[i].type == PM_DATA_STRING) {
        fprintf (stderr, " string val = %s \n", dlist.list[i].strings[0].c_str());
        }

      else if (dlist.list[i].type == PM_DATA_NUMBER) {
        fprintf (stderr, " number val = %f \n", dlist.list[i].numbers[0]);
        }

      else if (dlist.list[i].type == PM_DATA_NUMBER_LIST) {
        fprintf (stderr, " number list = ");
        for (unsigned int j = 0; j < dlist.list[i].numbers.size(); j++) {
          fprintf (stderr, " %f ", dlist.list[i].numbers[j]);
          }
        fprintf (stderr, "\n");
        }

      else if (dlist.list[i].type == PM_DATA_STRING_LIST) {
        fprintf (stderr, " string list = ");
        for (unsigned int j = 0; j < dlist.list[i].strings.size(); j++) {
          fprintf (stderr, " %s ", dlist.list[i].strings[j].c_str());
          }
        fprintf (stderr, "\n");
        }

      else { 
        fprintf (stderr, " unknown type = %d \n", dlist.list[i].type);
        }
      }
    }
  }

//*============================================================*
//*==========              processScript             ==========*
//*============================================================*
// process a script.                                         

void
PmCmd::processScript (const string file_name)
  {
  char c, last_c; 
  string line;
  int cont;
  FILE *fp;
  bool in_comment, start; 

  //static CmdState state;
  cont = 0;

  if ((fp = pmSystem.getCmdFilePtr()) != NULL) {  
    //fprintf (stderr, " ---- wait  get fp = %x \n", fp); 
    }
  else if (file_name.empty()) {
    return;
    }
  else {
    fp = fopen (file_name.c_str(), "r");
    //fprintf (stderr, " fp = %x \n", fp); 
    }

  if (!fp) {
    if (file_name != ".pm") {
      pm_ErrorReport (PM, "can't open script named \"%s\".", "*", file_name.c_str());
      }
    return;
    }

  //===== process input line =====//

  // remove leading spaces //

  last_c = '\0';
  in_comment = false; 
  start = true;

  while ((c = fgetc(fp)) != EOF) {
    //fprintf (stderr, " c = %c  last_c = %c  inc = %d \n", c, last_c, in_comment); 

    if (start && (c == ' ')) {
      continue;
      }

    start = false;

    if ((c < 32) && (c != '\n')) {
      c = ' ';
      }

    if (c == '\n') {
      if (line != "") {
        if (pmSystem.getCmdEcho()) {
          fprintf (stderr, "pm> %s\n", line.c_str());
          }

        pm_CmdProcess (fp, line);
        line.clear();
        in_comment = false; 
        }

      if (pmSystem.getCmdWait()) {  
        pmSystem.setCmdFilePtr (fp);
        //fprintf (stderr, " ++++ wait  set fp = %x \n", fp); 
        return;
        }
      }
    else if (c == '\\') {
      while (((c = fgetc(fp)) != EOF) && (c != '\n'));
      while (((c = fgetc(fp)) != EOF) && (c == ' '));
      line.push_back(' '); 

      if (c > 31) { 
        line.push_back(c); 
        }
      }
    else if (c == '#') { 
      in_comment = true; 
      line.push_back(c); 
      }

    else if ((c == '{') && (last_c == '$') && (!in_comment)) { 
      if (!pm_CmdSubstVariable(fp, line)) {
        return;
        }
      }
    else {
      if ((c != ' ') || (last_c != ' ') || (c != '\r') || (c != '\t')) {
        line.push_back(c); 
        }
      }

    last_c = c;
    }

  fclose (fp);
  }

//*============================================================*
//*==========              testScript                ==========*
//*============================================================*
// test script processing.

void
PmCmd::testScript (const string file_name)
  {

  PmCmdDataList dlist;
  PmCmd cmd; 
  string line;
  ifstream cmd_file; 

  cmd_file.open (file_name.c_str());

  while (!cmd_file.eof()) {
    getline (cmd_file, line);

    if (line == "") {
      continue;
      }

    cmd.setLine (line);
    cmd.getDataList(dlist);
    cmd.printDataList(dlist);

    if (!dlist.list.size()) {
      continue;
      }

    PmCmdData data;
    dlist.getNext (data);
    fprintf (stderr, "\n------ cmd [%s] ------ \n", data.name.c_str()); 

    while (dlist.getNext (data)) {
      fprintf (stderr, "data [%s] \n", data.name.c_str()); 

      if (data.name == "name") {
        string s;
        data.getString(s);
        fprintf (stderr, "    name data [%s] \n",s.c_str()); 
        }

      else if (data.name == "domain") {
        string s;
        //data.getString(s);
        //fprintf (stderr, "    domain data [%s] \n",s.c_str()); 

        vector<string> slist;
        data.getStringList(slist);
        fprintf (stderr, "    domain data:  "); 

        for (unsigned int i = 0; i < slist.size(); i++) {
          fprintf (stderr, " %s ", slist[i].c_str()); 
          }

        fprintf (stderr, "\n"); 
        }

      else if (data.name == "size") {
        float size = data.getFloat();
        fprintf (stderr, "    size data [%f] \n", size); 
        }

      else if (data.name == "pos") {
        PmVector3 vec;
        data.getVector(vec);
        fprintf (stderr, "    pos data (%f %f %f) \n", vec[0], vec[1], vec[2]); 
        }
      }

    fprintf (stderr, "-------------------------------- \n"); 

    string s;
    dlist.getString ("name", s);
    fprintf (stderr, ">>>> name string [%s] \n",s.c_str()); 

    dlist.getString ("type", s);
    fprintf (stderr, ">>>> type string [%s] \n",s.c_str()); 

    float size;
    dlist.getFloat ("size", size);
    fprintf (stderr, ">>>> size number [%g] \n", size); 
    }
  }

//*============================================================*
//*==========              callMacro                 ==========*
//*============================================================*
// process a macro call.         

void
PmCmd::callMacro (FILE *fp, CmdMacro *macro, PmCmdDataList& dlist)
  {
  PmCmdData data;
  vector<string> args;
  string line, sline, var;
  PmCmd cmd;
  bool found_var;
  bool echo = pmSystem.getCmdEcho();

  #define ndbg_callMacro 
  #ifdef dbg_callMacro 
  fprintf (stderr, "\n>>>>>> PmCmd::callMacro \n");
  fprintf (stderr, ">>> name=%s\n", macro->name.c_str()); 
  fprintf (stderr, ">>> args names=");

  for (int i = 0; i <  macro->args.size(); i++) {
    fprintf (stderr, "%s ", macro->args[i].c_str());
    }

  fprintf (stderr, "\n");
  fprintf (stderr, ">>> args values=");
  #endif

  while (dlist.getNext (data)) {
    #ifdef dbg_callMacro 
    fprintf (stderr, "%s ", data.name.c_str()); 
    #endif
    args.push_back(data.name);
    }

  #ifdef dbg_callMacro 
  fprintf (stderr, "\n");
  #endif

  if (args.size() != macro->args.size()) {
    pm_ErrorReport (PM, "number of arguments don't agree.", "*");
    return;
    }

  #ifdef old_PmCmd_callMacro 
 
  for (unsigned int i = 0; i < args.size(); i++) {
    cmd.addVariable(macro->args[i], args[i]);
    }

  // perform argument substitution //

  for (unsigned int i = 0; i < macro->body.size(); i++) {
    line = macro->body[i];
    cmd.setLine(line);
    found_var = false;
    #ifdef dbg_callMacro 
    fprintf (stderr, ">>> line=%s\n", line.c_str()); 
    #endif

    while (cmd.getChar()) {
      if (cmd.curr_c < 32) cmd.curr_c = 32;
      sline.push_back(cmd.curr_c);
      //fprintf (stderr, " %c(%d) \n", cmd.curr_c, cmd.curr_c); 
      }

    #ifdef dbg_callMacro 
    fprintf (stderr, "    sline=%s\n", sline.c_str()); 
    #endif
    pm_CmdProcess (fp, sline);
    sline.clear();
    }
  #else

  for (unsigned int i = 0; i < macro->body.size(); i++) {
    line = macro->body[i];
    #ifdef dbg_callMacro 
    fprintf (stderr, "\n>>> %d line=%s\n", i, line.c_str()); 
    #endif
    sline.clear();
    string name, value;
    bool found;

    for (unsigned int j = 0; j < line.size();) {
      char c = line[j++];

      if (c == '$') {
        name.clear();
        value.clear();

        while (j < line.size()) {
          c = line[j++];

          if (c == '}') {
            break;
            }
          else if (c != '{') {
            name.push_back(c);
            }
          }

        #ifdef dbg_callMacro 
        fprintf (stderr, ">>> name=%s\n", name.c_str()); 
        #endif
        found = false;

        // search macro arguments //

        for (int k = 0; k <  macro->args.size(); k++) {
          if (name == macro->args[k]) {
            value = args[k];
            found = true;
            break;
            }
          }

        // if not found then search variables //

        if (!found) {
          found = pm_CmdGetVariableValue(name, value);
          }

        if (found) {
          #ifdef dbg_callMacro 
          fprintf (stderr, ">>> subst name=%s with value=%s \n", name.c_str(), 
                   value.c_str()); 
          #endif

          for (int k = 0; k < value.size(); k++) {
            sline.push_back(value[k]);
            }
          }

        else {
          sline.push_back('$');
          sline.push_back('{');
          for (int k = 0; k <  name.size(); k++) {
            sline.push_back(name[k]);
            }
          sline.push_back('}');
          }
        }
      else {
        sline.push_back(c);
        }
      }

    #ifdef dbg_callMacro 
    fprintf (stderr, ">>> sline=%s\n", sline.c_str()); 
    #endif
    if (echo) {
      fprintf (stderr, "%s\n", sline.c_str()); 
      }
    pm_CmdLineProcess (fp, sline);
    sline.clear();
    }

  #endif
  }

//*============================================================*
//*==========              setVariableValue          ==========*
//*============================================================*
// set the value of a variable.

bool
pm_CmdSetVariableValue(const string name, const string line) 
  {
  //fprintf (stderr, "######## setVariableValue:  name %s\n", name.c_str()); 
  string vname, value;
  int n;
  bool equal_found=false;
  bool plus_found=false;
  bool minus_found=false;

  for (n=0; n < line.size(); n++) {
    if (line[n] == '=') {
      n += 1;
      equal_found = true;
      break;
      }
    else if (line[n] != ' ') {
      vname.push_back(line[n]);
      }
    }

  if (!equal_found) {
    pm_ErrorReport (PM, "bad variable set value", "*");
    return false;
    }

  for (; n < line.size(); n++) {
    if (line[n] != ' ') {
      break;
      }
    }

  for (; n < line.size(); n++) {
    if (line[n] == '+') {
      plus_found = true;
     }
    else if (line[n] == '-') {
      minus_found = true;
     }
    else {
       value.push_back(line[n]);
     }
  }

  //fprintf (stderr, "######## setVariableValue:  value %s\n", value.c_str()); 

  PmVariable *var;
  pm_CmdGetVariable(name, &var);

  if (plus_found) { 
    int ival1 = atoi(var->value.c_str());
    int ival2 = atoi(value.c_str());
    int sum = ival1 + ival2;
    ostringstream convert;   
    convert << sum;  
    value = convert.str(); 
  }
  else if (minus_found) { 
    int ival1 = atoi(var->value.c_str());
    int ival2 = atoi(value.c_str());
    int sum = ival1 - ival2;
    ostringstream convert;   
    convert << sum;  
    value = convert.str(); 
  }

  var->value = value;
  return true;
  }

//*============================================================*
//*==========         pm_CmdSubstVariable            ==========*
//*============================================================*
// substitute a variable.

bool 
pm_CmdSubstVariable (FILE *fp, string& line)
  {
  //fprintf (stderr, "\n>>>>>> pm_CmdSubstVariable: \n");
  //fprintf (stderr, ">>> begin line=%s \n", line.c_str());
  string name, value;
  bool bad_name, ifdef_cmd;
  char c;
  int n;
  string::iterator it;

  //===== process variable name =====//

  bad_name = false;

  while ((c = fgetc(fp)) != EOF) {
    if (c == '}') {
      break;
      }
    else if (isalpha(c) || isdigit(c) || (c == '_')) {
      name.push_back(c);
      }
    else {
      bad_name = true;
      break;
      }
    }

  if (bad_name) {
    pm_ErrorReport (PM, "bad variable name \"%s\" ", "*", name.c_str());
    return false;
    } 

  if (line.compare(0, 6, "ifdef ") == 0) { 
    line.push_back('{');

    for (unsigned int i = 0; i < name.size(); i++) {
      line.push_back(name[i]);
      }

    line.push_back('}');
    return true;
    }

  //fprintf (stderr, "   >>> name = %s \n", name.c_str());

  //===== seach for variable name =====//

  PmVariable *var;

  if (pm_CmdGetVariable(name, &var)) {
    value = var->value;
    }

  if (value.empty()) {
    pm_ErrorReport (PM, "unknown variable named \"%s\" ", "*", name.c_str());
    return false;
    }

  if (var->persistent) {
    line.push_back('{');

    for (unsigned int i = 0; i < name.size(); i++) {
      line.push_back(name[i]);
      }

    line.push_back('}');
    return true;
    }

  // remove '$' in line[] //

  n = line.size();
  it = line.begin() + n-1;
  line.erase (it);

  // if the variable is a boolean expression then evaluate it //

  /*
  string bval;

  if (var.bexp) {
    var.bexp->evaluate(bval);
    for (unsigned int i = 0; i < bval.size(); i++) {
      line.push_back(bval[i]);
      }
    }
  */

  // add variable value to input line //

  for (unsigned int i = 0; i < value.size(); i++) {
    line.push_back(value[i]);
    }

  return true;
  }

////////////////////////////////////////////////////////////////
//                    PmCmdDataList                          //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              check                     ==========*
//*============================================================*
// check to see if all data has been read.

bool
PmCmdDataList::check()
  {
  int n = list.size();
  PmCmdData data;
  int num_unread = 0;

  if (pmSystem.getCmdError()) {
    return true;
    }

  for (int i = 1; i < n; i++) {
    data = list[i];

    if (!data.read) {
      pm_ErrorReport (PM, "unknown argument \"%s\" ", "*", data.name.c_str());
      num_unread += 1;
      }
    }

  if (num_unread != 0) { 
    pmSystem.setCmdError(true);
    return false;
    }

  return true;
  }

//*============================================================*
//*==========              getBoolean                ==========*
//*============================================================*
// get a boolean from a data list.

bool
PmCmdDataList::getBoolean(const string name, bool& val)
  {
  PmCmdData data;
  val = false;

  if (getData(name, data)) {
    val = data.getBoolean();
    return true;
    }

  return false;
  }

//*============================================================*
//*==========              getNext                   ==========*
//*============================================================*
// get the next data item from the list that has not been read. 

bool 
PmCmdDataList::getNext(PmCmdData& data)
  {
  data.init();
  int n = list.size();

  if (count == n) {
    return false;
    } 

  while (count < n) {
    data = list[count];
    count++;

    if (!data.read) {
      list[count-1].read = true;
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========              getNext                   ==========*
//*============================================================*
// get the next data item from the list.

bool
PmCmdDataList::getNext(PmCmdData& data, string& val)
  {
  val.clear();

  if (getNext(data)) {
    data.getString(val);
    return true;
    }

  return false;
  }

//*============================================================*
//*==========              getData                   ==========*
//*============================================================*
// get a data item from the list.

bool
PmCmdDataList::getData (const string name, PmCmdData& data)
  {
  data.init();

  for (unsigned int i = 0; i < list.size(); i++) {
    if (name == list[i].name) {
      data = list[i];
      list[i].read = true;
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========              getFloat                  ==========*
//*============================================================*
// get a float from a data list.

bool 
PmCmdDataList::getFloat (const string name, float& val)
  {
  PmCmdData data;

  if (getData(name, data)) {
    val = data.getFloat();
    return true; 
    }

  return false;
  }

//*============================================================*
//*==========              getInt                    ==========*
//*============================================================*
// get an integer from a data list.

bool
PmCmdDataList::getInt (const string name, int& val)
  {
  PmCmdData data;

  if (getData(name, data)) {
    float fval = data.getFloat();
    val = (int)fval;
    return true;
    }

  return false;
  }

//*============================================================*
//*==========              getString                 ==========*
//*============================================================*
// get a string from a data list.

bool 
PmCmdDataList::getString(const string name, string& s)
  {
  s.clear();
  PmCmdData data;

  if (getData(name, data)) {
    data.getString(s);
    return true;
    }

  return false;
  }

//*============================================================*
//*==========              getStringList             ==========*
//*============================================================*
// get a list string data value.

bool
PmCmdDataList::getStringList(const string name, vector<string>& slist)
  {
  PmCmdData data;

  if (getData(name, data)) {
    data.getStringList(slist);
    return true;
    }

  return false;
  }

//*============================================================*
//*==========              getVector                 ==========*
//*============================================================*
// get a vector from a data list.

bool
PmCmdDataList::getVector(const string name, PmVector3& vec)
  {
  PmCmdData data;

  if (getData(name, data)) {
    return data.getVector(vec);
    }

  return false;
  }

//*============================================================*
//*==========              getColor                  ==========*
//*============================================================*
// get a color vector from a data list.

bool
PmCmdDataList::getColor(const string name, PmVector3& vec)
  {
  PmCmdData data;

  if (getData(name, data)) {
    return data.getColor(vec);
    }

  return false;
  }

////////////////////////////////////////////////////////////////
//                    PmCmdData                              //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              getBoolean                ==========*
//*============================================================*
// evaluate a boolean data value.        

bool
PmCmdData::getBoolean ()
  { 
  if (strings.empty()) {
    return true;
    }

  if ((strings[0] == "true") || 
      (strings[0] == "on")) {
    return true;
    }

  return false;
  }

//*============================================================*
//*==========              getFloat                  ==========*
//*============================================================*
// get a float data value.

float
PmCmdData::getFloat()
  {
  if (numbers.empty()) {
    if (strings.empty()) {
      return 0.0;
      }

    float val;

    if (convToFloat(strings[0], val)) {
      return val;
      }

    return 0.0;
    }

  return numbers[0];
  }

//*============================================================*
//*==========              getInt                    ==========*
//*============================================================*
// get an integer data value.

int
PmCmdData::getInt()
  {
  float fval = getFloat();
  return (int)fval;
  }

//*============================================================*
//*==========              getString                 ==========*
//*============================================================*
// get a string data value.

void 
PmCmdData::getString(string& s)
  {
  s.clear();

  if (!numbers.empty()) {
    int val = (int)numbers[0];
    stringstream out;
    out << val;
    s = out.str();
    return;
    }

  if (!strings.empty()) {
    s = strings[0];
    }
  }

//*============================================================*
//*==========              getStringList             ==========*
//*============================================================*
// get a list string data value.

void
PmCmdData::getStringList(vector<string>& slist)
  {
  slist.clear();

  if (strings.empty()) {
    return;
    }

  for (unsigned int i = 0; i < strings.size(); i++) {
    slist.push_back (strings[i]);
    }
  }

//*============================================================*
//*==========              getVector                 ==========*
//*============================================================*
// get a float vector data value.

bool  
PmCmdData::getVector(PmVector3& vec)
  {
  vec.set(0,0,0);

  if (numbers.empty()) {
    return false;
    }

  vec[0] = numbers[0];
  vec[1] = numbers[1];
  vec[2] = numbers[2];
  return true;
  }

//*============================================================*
//*==========              getColor                  ==========*
//*============================================================*
// get a color float vector data value.

bool
PmCmdData::getColor(PmVector3& vec)
  {
  string s;
  vec.set(0,0,0);

  if (!numbers.empty()) {
    vec[0] = numbers[0];
    vec[1] = numbers[1];
    vec[2] = numbers[2];
    return true ;
    }

  if (strings.empty()) {
    return false; 
    }

  s = strings[0];

  if (s == "red") {
    vec.set(1,0,0);
    }
  else if (s == "green") {
    vec.set(0,1,0);
    }
  else if (s == "blue") {
    vec.set(0,0,1);
    }
  else {
    pm_ErrorReport (PM, "unknown color \"%s\" ", "*", s.c_str());
    pmSystem.setCmdError(true);
    return false; 
    }

  return true;
  }

//*============================================================*
//*==========            pm_CmdInterruptProc         ==========*
//*============================================================*
// process an interrupt.

void
pm_CmdContinueProc()
  {
  //fprintf (stderr, "\n\n**** continue proc **** \n\n");
  PmCmd::processScript("");
  }

//*============================================================*
//*==========            pm_CmdInterruptProc         ==========*
//*============================================================*
// process an interrupt.

void
pm_CmdInterruptProc()
  {
  fprintf (stderr, "\n\n**** interrupt **** \n\n");
  pmSystem.procInterrupt();
  }

//*============================================================*
//*==========              pm_CmdKeyboardProc        ==========*
//*============================================================*

void
pm_CmdAddKeyboardEvent (const string cmd, const string key)
  {
  /*
  fprintf (stderr, "\n>>>>>> pm_CmdAddKeyboardEvent \n");
  fprintf (stderr, ">>> cmd=%s \n", cmd.c_str());
  fprintf (stderr, ">>> key=%s \n", key.c_str());
  */
  CmdKeyboardProcData data;
  data.cmd = cmd;
  data.key = key;
  g_key_board_proc.push_back(data);
  }

//*============================================================*
//*==========              pm_CmdKeyboardProc        ==========*
//*============================================================*

void
pm_CmdKeyboardProc (const string line)
  {
  /*
  fprintf (stderr, "\n>>>>>> pm_CmdKeyboardProc \n");
  fprintf (stderr, ">>> line=%s \n", line.c_str());
  */
  PmCmd cmd;
  PmCmdData data;
  PmCmdDataList dlist;
  string dn, dv;

  cmd.setLine (line);
  cmd.getDataList(dlist);

  if (dlist.error) {
    pm_ErrorReport (PM, "bad keyboard syntax.", "*");
    return;
    }

  dlist.getNext(data);
  dlist.getNext(data);
  //fprintf (stderr, ">>> data name=%s \n", data.name.c_str());

  for (unsigned int i = 0; i < g_key_board_proc.size(); i++) {
    if (data.name == g_key_board_proc[i].key) {
      pm_CmdProcess (0, g_key_board_proc[i].cmd);
      //fprintf (stderr, ">>> found cmd=%s \n", g_key_board_proc[i].cmd.c_str()); 
      }
    }
  }

//*============================================================*
//*==========              pm_CmdForLoop             ==========*
//*============================================================*
// process a for loop.   

void
pm_CmdForLoop (FILE *fp, const string cmd)
  {
  #ifdef dbg_pm_CmdForLoop 
  fprintf (stderr, ">>>>>> pm_CmdForLoop: \n");
  fprintf (stderr, ">>> cmd=%s \n", cmd.c_str());
  #endif
  vector <string> tokens;
  string token, loop_var, loop_range;
  char c, next_c='\0';
  bool in_token=false;
  int lmin, lmax, linc; 
  bool echo = pmSystem.getCmdEcho();

  for (unsigned int i = 0; i < cmd.size(); i++) {
    c = cmd[i];

    if ((c == ' ') && (next_c == ' ')) {
      continue;
      }

    if ((c == ' ') || (c == '=') || (c == ':')) { 
      if (in_token) {
        tokens.push_back(token);
        token.clear();
        in_token = false;
        }
      }
    else {
      in_token = true;
      token.push_back(c);
      }
    }

  tokens.push_back(token);

  #ifdef dbg_pm_CmdForLoop 
  fprintf (stderr, ">>> tokens= ");
  for (unsigned int i = 0; i < tokens.size(); i++) {
    fprintf (stderr, "%s,", tokens[i].c_str());
    }
  fprintf (stderr, "\n");
  #endif

  if (tokens.size() < 4) {
    pm_ErrorReport (PM, "bad for loop syntax.", "*");
    return;
    }

  loop_var = tokens[1];
  lmin = atoi(tokens[2].c_str());
  lmax = atoi(tokens[3].c_str());

  if (tokens.size() == 5) {
    linc = atoi(tokens[4].c_str());
    }
  else {
    linc = 1; 
    }

  #ifdef dbg_pm_CmdForLoop 
  fprintf (stderr, ">>> loop var=%s\n", loop_var.c_str());
  fprintf (stderr, ">>> loop min=%d\n", lmin); 
  fprintf (stderr, ">>> loop max=%d\n", lmax); 
  fprintf (stderr, ">>> loop inc=%d\n", linc); 
  #endif

  //==== parse for body =====//

  int i, n, num_blines;
  char *s, last_c;
  string line;
  vector<string> body_lines;

  last_c = '\0';
  num_blines = 0;

  while ((c = fgetc(fp)) != EOF) {
    if (c == '\n') {
      if ((line.find("end ") != string::npos) && (line.find(" for") != string::npos))  {
        break;
        }

      num_blines += 1;
      body_lines.push_back(line);
      #ifdef dbg_pm_CmdForLoop 
      fprintf (stderr, ">>> line=%s\n", line.c_str()); 
      #endif
      line.clear();
      }

    else if (c == '\\') {
      while (((c = fgetc(fp)) != EOF) && (c != '\n'));
      while (((c = fgetc(fp)) != EOF) && (c == ' '));
      line.push_back(' ');
      line.push_back(c);
      }
    else {
      if ((c != ' ') || (last_c != ' ')) {
        line.push_back(c);
        }
      }

    last_c = c;
    }

  #ifdef dbg_pm_CmdForLoop 
  fprintf (stderr, ">>> number of body lines=%d\n", body_lines.size()); 
  #endif

  //===== process for body =====//

  string sline, name, value;
  stringstream dss;
  i = lmin;

  while (true) {
    for (unsigned int j = 0; j < body_lines.size(); j++) {
      line = body_lines[j];

      for (unsigned int k = 0; k < line.size();) {
        c = line[k++];

        if (c == '$') {
          name.clear();

          while (k < line.size()) {
            c = line[k++];

            if (c == '}') {
              break;
              }
            else if (c != '{') {
              name.push_back(c);
              }
            }

          #ifdef dbg_pm_CmdForLoop 
          fprintf (stderr, ">>> name=%s\n", name.c_str()); 
          #endif

          if (name == loop_var) {
            dss << i;
            value = dss.str();
            dss.str(std::string());

            for (int m = 0; m < value.size(); m++) {
              sline.push_back(value[m]);
              }
            }

          else {
            if (!pm_CmdGetVariableValue(name, value)) {
              pm_ErrorReport (PM, "unknown variable \"%s\".", "*", name.c_str());
              return;
              }

            for (int m = 0; m < value.size(); m++) {
              sline.push_back(value[m]);
              }
            }
          }
        else {
          sline.push_back(c);
          }
        }

      if (echo) {
        fprintf (stderr, "%s\n", sline.c_str()); 
        }

      pm_CmdProcess (fp, sline);
      sline.clear();
      }

    if (lmin <= lmax) {
      i += linc;
      if (i > lmax) break;
      }
    else {
      i += linc;
      if (i < lmax) break;
      }
    }
  }

//*============================================================*
//*==========              pm_CmdIfdef               ==========*
//*============================================================*
// process an ifdef.

void
pm_CmdIfdef (FILE *fp, const string cmd)
  {
  int n, num_blines;
  char c, next_c, *s, last_c, str[1000];
  string name, line; 
  bool eval, in_name, in_else, start, in_comment;

  #define ndbg_pm_CmdIfdef      
  #ifdef dbg_pm_CmdIfdef      
  fprintf (stderr, "\n------ ifdef ------ \n");
  #endif
  in_name = false;

  for (unsigned int i = 0; i < cmd.size()-1; i++) {
    c = cmd[i];
    next_c = cmd[i+1];

    if (next_c == '}') {
      name.push_back(c);
      break;
      }

    if (c == '$') {
      if (next_c != '{') {
        pm_ErrorReport (PM, "bad variable name.", "*");
        return;
        } 
      in_name = true;
      i++;
      }
    else if (in_name) {
      name.push_back(c);
      }
    }

  if (next_c != '}') {
    pm_ErrorReport (PM, "bad variable name.", "*");
    return;
    }

  #ifdef dbg_pm_CmdIfdef      
  fprintf (stderr, "variable=%s \n", name.c_str());
  #endif

  //===== look for variable =====//

  bool var_found = pm_CmdHasVariable(name);

  //===== process ifdef body =====//

  last_c = '\0';
  num_blines = 0;
  in_else = false;
  start = true;
  in_comment = false; 

  while ((c = fgetc(fp)) != EOF) {
    if (start && (c == ' ')) {
      continue;
      }

    start = false;

    if (c == '\n') {
      if ((line.find("end ") != string::npos) && (line.find(" ifdef") != string::npos))  {
        #ifdef dbg_pm_CmdIfdef      
        fprintf (stderr, ">>> found end \n"); 
        #endif
        break;
        }

      if ((line.find("else") != string::npos))  {
        in_else = true;
        #ifdef dbg_pm_CmdIfdef      
        fprintf (stderr, ">>> found else \n"); 
        #endif
        }

      else if (!in_comment) {
        if (var_found && !in_else) {
          #ifdef dbg_pm_CmdIfdef      
          fprintf (stderr, ">>> process line=%s \n", line.c_str()); 
          #endif
          pm_CmdProcess (fp, line);
          }
        else if (in_else && !var_found) {
          pm_CmdProcess (fp, line);
          }
        }

      in_comment = false;
      num_blines += 1;
      line.clear();
      }

    else if (c == '#') {
      in_comment = true;
      }

    else if (c == '\\') {
      while (((c = fgetc(fp)) != EOF) && (c != '\n'));
      while (((c = fgetc(fp)) != EOF) && (c == ' '));
      line.push_back(' ');
      line.push_back(c);
      }

    else if ((c == '{') && (last_c == '$') && (!in_comment) && (var_found)) {
      if (!pm_CmdSubstVariable(fp, line)) {
        return;
        }
      }

    else {
      if ((c != ' ') || (last_c != ' ')) {
        line.push_back(c);
        }
      }

    last_c = c;
    }
  }

//*============================================================*
//*==========              pm_CmdConditional         ==========*
//*============================================================*
// process a conditional.        

void
pm_CmdConditional (FILE *fp, const string cmd)
  {
  int n, num_blines;
  char c, *s, last_c;
  char str[1000];
  string line, lhs, rhs, bline;
  bool eval, in_else;

  #define ndbg_pm_CmdConditional
  #ifdef dbg_pm_CmdConditional
  fprintf (stderr, "\n------ conditional ------ \n");
  #endif

  for (unsigned int i = 0; i < cmd.size(); i++) {
    str[i] = cmd[i];
    }

  str[cmd.size()] = '\0';
  #ifdef dbg_pm_CmdConditional
  fprintf (stderr, ">>> str=\"%s\" \n", str);
  #endif

  //===== process if clause =====//

  s = &str[0];

  while (*s != '\0') {
    if ((*s != ' ') || (*s == '\n')) {
      break;
      }

    s++;
    }

  if (*s == '\n') {
    return;
    }

  s += 2;
  n = 0;

  while (*s != '\0') {
    if ((*s == '=') || ((*s == ' ') && n)) {
      s++;
      break;
      }

    else if (*s != ' ') {
      lhs.push_back(*s);
      n += 1;
      }

    s++;
    }

  while (*s != '\0') {
    if ((*s != ' ') && (*s != '=')) {
      rhs.push_back(*s);
      }

    s++;
    }

  #ifdef dbg_pm_CmdConditional
  fprintf (stderr, ">>> lhs=\"%s\" \n", lhs.c_str());
  fprintf (stderr, ">>> rhs=\"%s\" \n", rhs.c_str());
  #endif

  if (rhs.empty() && (lhs == "true")) {
    eval = true;
    }
  else if (lhs == rhs) {
    eval = true;
    }
  else {
    eval = false;
    }

  #ifdef dbg_pm_CmdConditional
  fprintf (stderr, ">>> eval=%d \n", eval); 
  #endif

  //===== process if body =====//

  last_c = '\0';
  num_blines = 0;
  in_else = false;

  while ((c = fgetc(fp)) != EOF) {
    if ((c == ' ') && line.empty()) {
      continue;
      }

    if (c == '\n') {
      if ((line.find("end ") != string::npos) && (line.find(" if") != string::npos))  {
        #ifdef dbg_pm_CmdConditional
        fprintf (stderr, ">>> found end \n"); 
        #endif
        break;
        }

      if ((line.find("else") != string::npos))  {
        in_else = true;
        #ifdef dbg_pm_CmdConditional
        fprintf (stderr, ">>> found else \n"); 
        #endif
        fprintf (stderr, ">>> found else \n"); 
        }

      else {
        if (eval && !in_else) {
          pm_CmdLineProcess (fp, line);
          //pm_CmdProcess (fp, line);
          #ifdef dbg_pm_CmdConditional
          fprintf (stderr, ">>> process if line=\"%s\"\n", line.c_str()); 
          #endif
          }
        else if (!eval && in_else) {
          pm_CmdLineProcess (fp, line);
          //pm_CmdProcess (fp, line);
          #ifdef dbg_pm_CmdConditional
          fprintf (stderr, ">>> process in_else line=%s \n", line.c_str()); 
          #endif
          }
        }

      num_blines += 1;
      line.clear();
      }

    else if (c == '\\') {
      while (((c = fgetc(fp)) != EOF) && (c != '\n'));
      while (((c = fgetc(fp)) != EOF) && (c == ' '));
      line.push_back(' ');
      line.push_back(c);
      }
    else {
      if ((c != ' ') || (last_c != ' ')) {
        line.push_back(c);
        }
      }

    last_c = c;
    }
  }

//*============================================================*
//*==========              pm_CmdMacroDefine         ==========*
//*============================================================*
// process a macro definition.

void
pm_CmdMacroDefine (FILE *fp, const string cmd)
  {
  int n, num_blines;
  char c, *s, last_c; 
  char str[1000]; 
  string line, name, arg;
  CmdMacro *macro;

  #define ndbg_pm_CmdMacroDefine 
  #ifdef dbg_pm_CmdMacroDefine 
  fprintf (stderr, "\n------ define macro ------ \n");
  #endif

  for (unsigned int i = 0; i < cmd.size(); i++) {
    str[i] = cmd[i];
    }

  // process macro name //

  s = &str[0];

  while (*s != '\0') {
    if ((*s != ' ') || (*s == '\n')) {
      break;
      }

    s++;
    }

  if (*s == '\n') {
    return;
    }

  s += 5;
  n = 0;

  while (*s != '\0') {
    if ((*s == '(') || ((*s == ' ') && n)) { 
      break;
      }

    else if (*s != ' ') { 
      name.push_back(*s);
      n += 1;
      }

    s++;
    }

  #ifdef dbg_pm_CmdMacroDefine 
  fprintf (stderr, ">>> name = \"%s\" \n", name.c_str());
  #endif

  // create macro object //

  macro = new CmdMacro;
  macro->name = name;

  // process argument list //

  while (*s != '\0') {
    if ((*s == ')') || (*s == ','))  {
      #ifdef dbg_pm_CmdMacroDefine 
      fprintf (stderr, ">>> arg = \"%s\" \n", arg.c_str());
      #endif
      macro->args.push_back(arg);
      arg.clear();

      if (*s == ')') {
        break;
        }
      }

    else if (isalpha(*s) || isdigit(*s) || (*s == '_')) {
      arg.push_back(*s);
      }

    s++;
    }

  //===== process body =====//

  #ifdef dbg_pm_CmdMacroDefine 
  fprintf (stderr, ">>> body -----> \n");
  #endif

  last_c = '\0';
  num_blines = 0;
  int nc = 0;

  while ((c = fgetc(fp)) != EOF) {
    //fprintf (stderr, ">>> %d:  c=%c (%d)\n", nc++, c, c);

    if ((c < 32) && (c != '\n')) {
      c = ' ';
      }

    if (c == '\n') {
      if (line != "") {
        //fprintf (stderr, ">>> line = %s\n", line.c_str());
        }

      if ((line.find("end") != string::npos) && (line.find("macro") != string::npos)) { 
        break;
        }

      macro->body.push_back(line);
      num_blines += 1;
      line.clear();
      }

    else if (c == '\\') {
      while (((c = fgetc(fp)) != EOF) && (c != '\n'));
      while (((c = fgetc(fp)) != EOF) && (c == ' '));
      line.push_back(' ');
      line.push_back(c);
      }
    else {
      if ((c != ' ') || (last_c != ' ')) {
        line.push_back(c);
        }
      }

    last_c = c;
    }

  if (num_blines == 0) { 
    pm_ErrorReport (PM, "macro has no commands.", "*");
    return;
    }

  #ifdef dbg_pm_CmdMacroDefine 
  fprintf (stderr, "<------- body \n");
  #endif
  g_macros.push_back(macro);
  #ifdef dbg_pm_CmdMacroDefine 
  fprintf (stderr, "------ finished define macro ------ \n");
  #endif
  }

//*============================================================*
//*==========              pm_CmdLineProcess         ==========*
//*============================================================*
// process a line.   

void
pm_CmdLineProcess (FILE *fp, const string line)
  {
  string name, value, pline;
  unsigned int i;

  for (i = 0; i < line.size()-1; i++) {
    if (line[i] != ' ') {
      break;
      }
    }

  for (; i < line.size();) {
    char c = line[i++];

    if (c == '$') {
      name.clear();
      value.clear();

      while (i < line.size()) {
        c = line[i++];

        if (c == '}') {
          break;
          }
        else if (c != '{') {
          name.push_back(c);
          }
        }

      if (pm_CmdGetVariableValue(name, value)) {
        for (int k = 0; k < value.size(); k++) {
          pline.push_back(value[k]);
          }
        }
      else {
        pline.push_back('$');
        pline.push_back('{');

        for (int k = 0; k <  name.size(); k++) {
          pline.push_back(name[k]);
          }

        pline.push_back('}');
        }
      }
    else {
      pline.push_back(c);
      }
    }

  //fprintf (stderr, ">>> pline=\"%s\" \n", pline.c_str());
  pm_CmdProcess (fp, pline);
  }

//*============================================================*
//*==========              pm_CmdProcess             ==========*
//*============================================================*
// process a command.

void
pm_CmdProcess (FILE *fp, const string line)
  {
  #define ndbg_pm_CmdProcess
  #ifdef dbg_pm_CmdProcess
  fprintf (stderr, "\n>>>>>> pm_CmdProcess: \n");
  fprintf (stderr, ">>> line=%s\n", line.c_str());
  #endif
  unsigned int i, n;

  void (*func)(PmCmdDataList& dlist);
  n = line.size();

  if (n == 0) { 
    return;
    }

  /*
  fprintf (stderr, ">>> line=%s \n", line.c_str()); 
  for (int i = 0; i < n; i++) {
    fprintf (stderr, "%d: %c (%d) \n", i, line[i], line[i]);
    }
  */

  char c = line[0];

  if (c == COMMENTCHAR) {
    return;
    }

  if ((line == "quit") || (c == 'q')) {
    pm_CmdQuit();
    }

  if (c == KEY_ESCAPE) {
    pm_CmdInterruptProc();
    return;
    }

  /*
  if (pmSystem.getCmdWait()) { 
    return;
    }
  */
  pmSystem.setCmdError(false);

  // just print the line //

  if ((n = line.find("print")) != string::npos) { 
    if (n <= 1) {
      fprintf (stderr, "%s\n", line.c_str()+5+n);
      return;
      }
    }

  //===== add a variable definition =====//

  if (line.compare(0, 8, "variable") == 0) { 
    pm_CmdAddVariable(line);
    return;
    }

  //===== process a macro definition =====//

  if (line.compare(0, 5, "macro") == 0) { 
    pm_CmdMacroDefine(fp, line);
    return;
    }

  //===== process a for-loop =====//

  if (line.compare(0, 4, "for ") == 0) {
    pm_CmdForLoop(fp, line);
    return;
    }

  //===== process a conditional =====//

  if (line.compare(0, 3, "if ") == 0) {
    pm_CmdConditional(fp, line);
    return;
    }

  //===== process an ifdef =====//

  if (line.compare(0, 6, "ifdef ") == 0) {
    pm_CmdIfdef (fp, line);
    return;
    }

  //===== parse line[] breaking it into data items (dlist), =====//
  //===== usually of the form  <name> = <value>.            ===== //

  PmCmd cmd;
  PmCmdData data;
  PmCmdDataList dlist;
  cmd.setLine (line);
  cmd.getDataList(dlist);

  #ifdef pm_CmdProcess
  cmd.printDataList (dlist) ;
  #endif

  if (dlist.error) {
    char curr_c;
    int cnum;
    curr_c = cmd.getCurrChar();
    cnum = cmd.getCurrCnum();
    pm_ErrorReport (PM, "bad command syntax at %d ", "*", cnum);
    cmd.printDataList (dlist) ;
    return;
    }

  // get the first data item on the line.
  // this should be a command name.

  dlist.getNext(data);

  if (data.name == "") {
    return;
    }

  //=====  read in one or more scripts  =====//

  if (data.name == "read") {
    while (dlist.getNext(data)) {
      //fprintf (stderr, "   >>> data[%s] \n", data.name.c_str());
      PmCmd::processScript(data.name);
      }
    return;
    }

  //===== process special function key =====//

  if (data.name == "keyboard") {
    pm_CmdKeyboardProc(line);
    return;
    }

  else if (data.name == "wait") {
    /*
    if (pmSystem.useGraphics()) {
      pmSystem.setCmdWait(true);
      }
    */
    pmSystem.setCmdWait(true);
    return;
    }

  // search for cmd name in the cmd table  //

  for (i = 0; cmd_table[i].name != NULL; i++) {
    if (data.name == cmd_table[i].name) {
      func = cmd_table[i].func;
      //fprintf (stderr, "   >>> dlist[%s] \n", dlist.list[0].name.c_str());
      (*func)(dlist);
      dlist.check();

      if (pmSystem.getCmdError()) {
        pm_ErrorReport (PM, "line = \"%s\" ", "*", line.c_str());
        }

      // Note [29dec09] will this have performance issues? //
      pmSystem.forceUpdateGraphics();
      return;
      }
    }

  // search for cmd name in the abreviated cmd table  //

  for (i = 0; abbrev_table[i].name != NULL; i++) {
    if (data.name == abbrev_table[i].name) {
      char str[200], *s;
      PmCmdData adata;
      strcpy(str, abbrev_table[i].full_name);
      //fprintf (stderr, ">>> abbrev str=\"%s\" \n", str); 

      if (abbrev_table[i].name == "md") {
        pm_CmdProcess(NULL, str);
        return;
        }

      s = strtok(str, " ");
      vector<PmCmdData>::iterator it; 
      int n = 0;

      while (s) {
        adata.init();
        adata.name = s;
        adata.type = PM_DATA_STRING;
        adata.strings.push_back(s);

        if (n) {
          it = dlist.list.begin();
          dlist.list.insert(it+n, adata);
          }
        else {
          dlist.list[0] = adata;
          }

        n += 1;
        s = strtok(NULL, " ");
        }

      //cmd.printDataList (dlist) ;

      func = abbrev_table[i].func;
      (*func)(dlist);
      dlist.check();

      if (pmSystem.getCmdError()) {
        pm_ErrorReport (PM, "line = \"%s\" ", "*", line.c_str());
        }

      return;
      }
    }

  //===== search for macro =====//

  for (unsigned int i = 0; i < g_macros.size(); i++) {
    if (g_macros[i]->name == data.name) {
      cmd.callMacro(fp, g_macros[i], dlist);
      return;
      }
    }

  //===== search for variable =====//

  if (pm_CmdHasVariable(data.name)) {
    if (pm_CmdSetVariableValue(data.name, line)) {
      return;
      }
    }

  pm_ErrorReport (PM, "unknown command \"%s\" ", "*", data.name.c_str());
  }

//*============================================================*
//*==========              pm_CmdHelp                ==========*
//*============================================================*
// process cmd help command.

void
pm_CmdHelp (PmCmdDataList& dlist)
  {

  PmCmdData data;

  string str;

  int i;

  dlist.getNext(data); 

  if (data.name != "") {
    pm_CmdProcess (0, data.name);
    return;
    }

  fprintf (stderr, "\n---------- commands ---------- \n");

  for (i = 0; cmd_table[i].name != NULL; i++) {
    fprintf (stderr, "%s \n", cmd_table[i].name);
    }

  fprintf (stderr, "\n---------- abbreviations ----------\n");

  for (i = 0; abbrev_table[i].name != NULL; i++) {
    fprintf (stderr, "%s --> %s \n", abbrev_table[i].name, abbrev_table[i].full_name);
    }
  }

//*============================================================*
//*==========         pm_CmdCheckSystemVar           ==========*
//*============================================================*
// check for a system variable.

bool 
pm_CmdCheckSystemVar(const string name) 
  {

  for (int i = 0; system_vars[i].size() != 0; i++) {
    if (system_vars[i] == name) {
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========      pm_CmdAddBooleanExpression        ==========*
//*============================================================*
// add a boolean expression.                   

void
pm_CmdAddBooleanExpression(const string name, const string value, 
                           PmBooleanExpression **pbexp)
  {
  #ifdef dbg_pm_CmdAddBooleanExpression
  fprintf (stderr, "\n>>>>>> pm_CmdAddBooleanExpression \n");
  #endif
  vector<string> vars(2);
  PmBooleanRelOpType relop;
  PmBooleanOpType op;
  vector<PmBooleanOpType> ops;
  vector<PmBooleanCondition> conds;
  char c, nc;
  bool first_alnum;
  unsigned int n=0; 
  int num_vars, num_cond=0, num_paren;
  op = PM_BOOLEAN_OP_UNKNOWN;
  num_paren = 0;
  *pbexp = NULL;

  while (n < value.size()) {
    c = value[n++];
    if (n < value.size()-1) {
      nc = value[n+1];
      }
    else {
      nc = '\0'; 
      }

    if (c == '(') {
      num_vars=0;
      num_paren += 1;
      vars[0].clear();
      vars[1].clear();
      vars.clear();

      while (n < value.size()) {
        c = value[n++];
        //fprintf (stderr, ">>> c=%c\n", c);

        if (c == ' ') {
          continue;
          }

        if (c == ')') {
          num_cond += 1;
          num_paren += 1;
          break;
          }

        if (isalnum(c) || c == '_' || c == '-' || c == '.') {
          vars[num_vars].push_back(c);
          }

        else {
          num_vars += 1;
          #ifdef dbg_pm_CmdAddBooleanExpression
          fprintf (stderr, ">>> num vars=%d\n", num_vars);
          #endif

          if (num_vars > 2) {
            return;
            }
          }

        if (c == '<') {
          if (nc == '=') {
            relop = PM_BOOLEAN_RELOP_LE;
            n += 1;
            }
          else {
            relop = PM_BOOLEAN_RELOP_LT;
            }
          } 

        else if (c == '>') {
          if (nc == '=') {
            relop = PM_BOOLEAN_RELOP_GE;
            n += 1;
            }
          else {
            relop = PM_BOOLEAN_RELOP_GT;
            }
          }

        else if (c == '=') {
          if (nc == '=') {
            relop = PM_BOOLEAN_RELOP_EQ;
            n += 1;
            }
          }

        else if (c == '!') {
          if (nc == '=') {
            relop = PM_BOOLEAN_RELOP_NE;
            n += 1;
            }
          }
        }

      #ifdef dbg_pm_CmdAddBooleanExpression
      fprintf (stderr, ">>> num_paren=%d num_vars=%d \n", num_paren, num_vars);
      #endif

      if ((num_paren != 2) || (num_vars != 1)) {
        return;
        }

      #ifdef dbg_pm_CmdAddBooleanExpression
      fprintf (stderr, ">>> cond %d = %s %d %s \n", num_cond, vars[0].c_str(), relop, 
               vars[1].c_str());
      #endif
      num_paren = 0;

      for (int i = 0; i < 2; i++) {
        string var = vars[i];

        if (isalpha(var[0])) {
          if (!pm_CmdCheckSystemVar(var)) { 
            pm_ErrorReport (PM, "unknown system variable \"%s\" ", "*", var.c_str());
            return;
            }
          }
        }

      PmBooleanCondition cond;
      cond.var1 = vars[0];
      cond.var2 = vars[1];
      cond.op = relop;
      conds.push_back(cond);
      }

    // process operator //

    else if (c == '&') {
      op = PM_BOOLEAN_OP_AND;
      ops.push_back(op);
      }
    else if (c == '|') {
      op = PM_BOOLEAN_OP_OR;
      ops.push_back(op);
      }
    else if (c == '!') {
      op = PM_BOOLEAN_OP_NOT;
      ops.push_back(op);
      }
    }

  PmBooleanExpression *bexp = new PmBooleanExpression;
  bexp->name = name;
  bexp->ops = ops;
  bexp->conditions = conds;
  g_boolean_expressions.push_back(bexp);
  *pbexp = bexp;
  }

//*============================================================*
//*==========      pm_CmdGetBooleanExpression        ==========*
//*============================================================*
// get a boolean expression.

bool 
pm_CmdGetBooleanExpression(const string name, PmBooleanExpression **bexp)
  {
  //fprintf (stderr, ">>>>>> pm_CmdGetBooleanExpression name=%s \n", name.c_str());
  string str;
  char c;
  *bexp = NULL;

  for (unsigned int i = 0; i < name.size(); i++) {
    c = name[i];

    if ((c != '$') && (c != '}') && (c != '{') && (c != ' ')) {
      str.push_back(name[i]);
      } 
    }

  for (unsigned int i = 0; i < g_boolean_expressions.size(); i++) {
    if (g_boolean_expressions[i]->name == str) {
      *bexp = g_boolean_expressions[i];
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========              pm_CmdAddVariable         ==========*
//*============================================================*
// process variable statements of the form: 
//                                  
//       variable <name> = <value> 

void
pm_CmdAddVariable (string line)
  {
  #ifdef dbg_pm_CmdAddVar 
  fprintf (stderr, ">>>>>> pm_CmdAddVar \n");
  fprintf (stderr, "line=\"%s\" \n", line.c_str());
  #endif

  PmCmdData data;
  PmVariable var;
  string name, value;
  unsigned int n;

  n = line.find("variable") + 8;

  for (; n < line.size(); n++) {
    if (line[n] == '=') {
      n += 1;
      break;
      }
    else if (line[n] != ' ') {
      name.push_back(line[n]);
      }
    }

  for (; n < line.size(); n++) {
    if (line[n] != ' ') {
      break;
      }
    }

 for (; n < line.size(); n++) {
   value.push_back(line[n]);
   }

  #ifdef dbg_pm_CmdAddVar 
  fprintf (stderr, "\n---------- define variable ----------\n");
  fprintf (stderr, ">>> name  = %s \n", name.c_str());
  fprintf (stderr, ">>> value = \"%s\" \n", value.c_str());
  #endif

  if (pm_CmdHasVariable(name)) {
    pm_ErrorReport (PM, "variable \"%s\" is already defined.", "*", name.c_str());
    return;
    }

  //===== process substitution for body names =====//

  size_t pos;

  if ((pos = value.find("body:")) != string::npos) {
    string str (value.begin()+5, value.end()); 
    string base_str, bname;
    char c;
    vector<PmBody*> bodies;
    //fprintf (stderr, ">>> str = \"%s\" \n", str.c_str());

    for (unsigned int i = 0; i < str.size(); i++) {
      c = str[i];

      if (c == '*') {
        break;
        }

      base_str.push_back(c);
      }

    //fprintf (stderr, ">>> base_str = \"%s\" \n", base_str.c_str());
    pmSystem.getBodies(bodies);
    value.clear();

    for (unsigned int i = 0; i < bodies.size(); i++) {
      bodies[i]->getName(bname);
      pos = bname.find(base_str);

      if (pos == 0) {
        value.append(bname);
        value.append(" ");
        }
      }

    //fprintf (stderr, ">>> value = \"%s\" \n", value.c_str());
    }

  //===== process substitution for domain names =====//

  if ((pos = value.find("domain:")) != string::npos) {
    string str (value.begin()+5, value.end());
    string base_str, bname;
    char c;
    vector<PmMolecule*> domains;

    for (unsigned int i = 0; i < str.size(); i++) {
      c = str[i];

      if (c == '*') {
        break;
        }

      base_str.push_back(c);
      }

    //fprintf (stderr, ">>> base_str = \"%s\" \n", base_str.c_str());
    pmSystem.getDomains(domains);
    value.clear();

    for (unsigned int i = 0; i < domains.size(); i++) {
      domains[i]->getName(bname);
      pos = bname.find(base_str);

      if (pos == 0) {
        value.append(bname);
        value.append(" ");
        }
      }
    }

  //===== check for a boolean expression =====//

  int num_paren = 0;
  int num_ops = 0;

  for (unsigned int i = 0; i < value.size(); i++) {
    if (value[i] == '(') {
      num_paren += 1;
      }
    }

  // check for a boolean expressions //

  PmBooleanExpression *bexp = NULL;

  if (num_paren != 0) { 
    pm_CmdAddBooleanExpression(name, value, &bexp);
    var.persistent = true;
    }
  else {
    var.persistent = false;
    }

  var.name = name;
  var.value = value;
  var.bexp = bexp;
  g_vars.push_back(var);
  }

//*============================================================*
//*==========          pm_CmdAddVariable             ==========*
//*============================================================*
// add a commnad-line variable. this is only called in pmsys.cpp.

void 
pm_CmdAddVariable(const string name, const string value)
  {
  //fprintf (stderr, ">>> add variable name=%s \n", name.c_str());
  PmVariable var;
  var.name = name;
  var.value = value;
  var.persistent = false;
  g_vars.push_back(var);
  }

//*============================================================*
//*==========          pm_CmdGetVariable             ==========*
//*============================================================*
// get a variable.                                                      

bool
pm_CmdGetVariable(const string name, PmVariable **var) 
  {
  for (unsigned int i = 0; i < g_vars.size(); i++) {
    if (g_vars[i].name == name) {
      *var = &g_vars[i];
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========          pm_CmdGetVariableValue        ==========*
//*============================================================*
// get a variable value

bool
pm_CmdGetVariableValue (const string name, string& value)
  {
  PmVariable *var;

  if (pm_CmdGetVariable(name, &var)) {
    value = var->value; 
    return true;
    }
   
  return false;
  }

//*============================================================*
//*==========          pm_CmdGetVariableValue        ==========*
//*============================================================*
// get a variable value

bool
pm_CmdGetVariableValue (PmCmd *cmd, int cnum, string line, int& n, string& value) 
  {
  string name, sym;
  //fprintf (stderr, ">>>>>> pm_CmdGetVarValue \n");

  for (unsigned int i = cnum; i < line.size(); i++) {
    sym.push_back(line[i]);

    if (line[i] == '}') {
      break;
      }
    }

  //fprintf (stderr, ">>> sym = %s \n", sym.c_str());
  n = sym.size();
  value.clear();

  if ((n < 4) || (sym[1] != '{') || (sym[n-1] != '}')) {
    fprintf (stderr, "**** Error: bad variable name \"%s\". \n", sym.c_str());
    return false;
    }

  for (unsigned int i = 2; i < sym.size()-1; i++) {
    name.push_back(sym[i]);
    }

  //fprintf (stderr, ">>> name = \"%s\" \n", name.c_str());

  if (pm_CmdGetVariableValue(name, value)) {
    return true;
    }

  fprintf (stderr, "**** Error: can't find variable name \"%s\". \n", sym.c_str());
  return false;
  }

//*============================================================*
//*==========              pm_CmdGetVariableData     ==========*
//*============================================================*
// get a variable.                            

bool
pm_CmdGetVariableData (const string name, PmCmdData& data)
  {
  for (unsigned int i = 0; i < g_vars.size(); i++) {
    }
  return true;
  }

//*============================================================*
//*==========              pm_CmdHasVariable         ==========*
//*============================================================*

bool
pm_CmdHasVariable (const string name)
  {
  PmVariable *var;
  return pm_CmdGetVariable(name, &var);
  }

//*============================================================*
//*==========         pm_CmdProcTimeInterval         ==========*
//*============================================================*
// process a time interval statement.

void
pm_CmdProcTimeInterval(PmCmdData& data, PmTimeInterval& time) 
  {
  if (!data.type == PM_DATA_STRING_LIST) {
    pm_ErrorReport (PM, "time interval must be a list.", "*");
    return;
    }

  vector<string> tlist;
  float bt, et; 
  data.getStringList(tlist);

  if (tlist.size() < 2) {
    pm_ErrorReport (PM, "time interval must have at least two items.", "*");
    return;
    }

  if (tlist.size() % 2) {
    pm_ErrorReport (PM, "time interval must have pairs of (begin,end) times.", "*");
    return;
    }

  for (unsigned int i = 0; i < tlist.size(); i += 2) {
    if (!convToFloat(tlist[i], bt)) {
      pm_ErrorReport (PM, "bad time value \"%s\".", "*", tlist[2*i].c_str());
      return;
      }

    if (!convToFloat(tlist[i+1], et)) {
      pm_ErrorReport (PM, "bad time value \"%s\".", "*", tlist[2*i+1].c_str());
      return;
      }
    
    time.addInterval(bt, et);
    }
  }

//*============================================================*
//*==========         pm_CmdProcPosition             ==========*
//*============================================================*
// process a position, point or location statement.

bool
pm_CmdProcPosition(PmCmdData& data, PmVector3& position) 
  {
  //fprintf (stderr, "\n---------- pm_CmdProcPosition ---------- \n");

  // data is a vector [ <x> <y> <z> ] //

  if (data.getVector(position)) {
    return true;
    }

  // data is a string of the form  <domain> : <residues> : <atoms> | mc | sc // 

  string str, dname, desc, atoms;
  vector<string> names(3);
  char c;
  PmAtomFilter filter;
  PmMolecule *domain;
  vector<PmVector3> coords;
  bool helix_center = false;
  int n = 0;

  data.getString(str);
  //fprintf (stderr, "   >>> str[%s]  \n", str.c_str()); 

  for (unsigned int i = 0; i < str.size(); i++) {
    c = str[i];

    if (c == ':') {
      n += 1;
      }
    else {
      names[n].push_back(c);
      }
    }

  dname = names[0];
  //fprintf (stderr, "   >>> domain[%s]  \n", dname.c_str()); 
  pmSystem.getDomain (dname, &domain);

  if (!domain) {
    pmSystem.getMolecule(dname, &domain);

    if (!domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
      return false;
      }
    }

  if (n > 0) {
    desc = names[1];
    }

  if (n == 2) {
    atoms = names[2];

    if (atoms == "mc") {
      filter.mainchain = true;
      domain->getBackboneAtomNames(filter.names);
      }
    else if (atoms == "sc") {
      domain->getBackboneAtomNames(filter.names);
      filter.sidechain = true;
      filter.exclude = true;
      }
    else if (atoms == "hc") {
      helix_center = true;
      }
    else { 
      filter.names.push_back(atoms);
      }
    }

  if (helix_center) {
    PmHelixProps helix_props;
    domain->getHelixProps(desc, helix_props);

    if (helix_props.origin.size() != 0) { 
      coords.push_back(helix_props.origin[0]);
      }
    }
  else {
    domain->getAtomCoords(desc, filter, coords);
    }

  if (!coords.size()) {
    pm_ErrorReport (PM, "residue/atoms \"%s\" not found.", "*", desc.c_str());
    return false;
    }

  position.set(0,0,0);

  for (unsigned int i = 0; i < coords.size(); i++) {
    position = position + coords[i];
    }

  position = (1.0/coords.size())*position;
  return true;
  }

//*============================================================*
//*==========             pm_CmdQuit                 ==========*
//*============================================================*
// process the quit command.                             

void 
pm_CmdQuit()
  {
  pmSystem.procQuit();
  }

////////////////////////////////////////////////////////////////
//         p r o c e s s   p m   c o m m a n d s             //
//////////////////////////////////////////////////////////////

// body commands
#include "cmd_body.cpp"

// binding simulation commands
#include "cmd_bsim.cpp"

// curve commands
#include "cmd_curve.cpp"

// database commands
#include "cmd_db.cpp"

// domain commands
#include "cmd_dom.cpp"

// force commands
#include "cmd_force.cpp"

// graphics commands
#include "cmd_gr.cpp"

// grid commands
#include "cmd_grid.cpp"

// joint commands
#include "cmd_joint.cpp"

// measurement commands
#include "cmd_msr.cpp"

// mulitbody commands
#include "cmd_mbody.cpp"

// model commands
#include "cmd_model.cpp"

// molecule commands
#include "cmd_mol.cpp"

// motor commands
#include "cmd_motor.cpp"

// particle commands
#include "cmd_part.cpp"

// potential commands
#include "cmd_pot.cpp"

// simulation commands
#include "cmd_sim.cpp"

// solid commands
#include "cmd_solid.cpp"

// surface commands
#include "cmd_surf.cpp"

// system commands
#include "cmd_sys.cpp"

// trace commands
#include "cmd_trace.cpp"

// units commands
#include "cmd_units.cpp"

// vector commands
#include "cmd_vector.cpp"


}

