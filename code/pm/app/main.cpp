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
//* app:    p r o t e i n   m e c h a n i c a   a p p          *
//*============================================================*

#include "pm/pm.h"

using namespace ProteinMechanica;

namespace ProteinMechanica {
  void pm_CmdProcess (FILE *fp, const string line);
  }

bool silent = false;

//*============================================================*
//*==========               procArgs                 ==========*
//*============================================================*

void
procArgs (int num_args, char **args, char **script)
  {
  char *str, c;
  string name, value;
  int n, i;
  *script = NULL;
  vector<string> cmd_files; 

  //===== process arguments =====//

  //fprintf (stderr,"##### proc args ##### \n");

  for (int j = 1; j < num_args; j++) {
    str = args[j];

    for (i = 0; i < strlen(str); i++) {
      if ((str[i] < 31) && (str[i] != '\n')) {
        str[i] = ' ';
        }
      }

    if (!strstr(str, "=")) {
      cmd_files.push_back(str);
      }

    // process variable of the form <var name> = <value> //

    else {
      n = strlen(str);
      name.clear();
      value.clear();

      for (i = 0; i < n; i++) {
        if (str[i] == '=') {
          i += 1;
          break;
          }
        else if (str[i] != ' ') {
          name.push_back(str[i]);
          }
        }

      for (; i < n; i++) {
        if (str[i] != ' ') {
          break;
          }
        }

      for (; i < n; i++) {
        value.push_back(str[i]);
        }

      pmSystem.addCmdVariable(name, value);

      if (name == "silent") {
        if (value == "true") {
          pmSystem.setSilent(true);
          }
        }
      }
    }

  // create script file with cmd file names //

  n = 0;

  for (unsigned int i = 0; i < cmd_files.size(); i++) {
    //fprintf (stderr,"##### cmd file=%s \n", cmd_files[i].c_str());
    n += cmd_files[i].size() + 1;
    }

  *script = new char[n+1];
  n = 0;

  for (unsigned int i = 0; i < cmd_files.size(); i++) {
    for (unsigned int j = 0; j < cmd_files[i].size(); j++) {
      (*script)[n++] = cmd_files[i][j];
      }

    (*script)[n++] = ' '; 
    }

  (*script)[n] = '\0'; 
  }

//*============================================================*
//*==========               procInitFile             ==========*
//*============================================================*

void
procInitFile() 
  {
  FILE *fp;
  char *fname = ".pm";
  char name[80], *str;
  bool found;

  //===== look in current directory =====// 

  strcpy (name, fname);

  if ((fp = fopen(name, "rb")) != NULL) {
    found = true;
    fclose(fp);
    }

  //===== look in home dir =====//

  else {
    if ((str = getenv("HOME")) != NULL) {
      for (int i = 0; i < strlen(str); i++) {
        if (str[i] == '\\') str[i] = '/';
        }

      sprintf (name, "%s/%s", str, fname);

      if ((fp = fopen(name, "rb")) != NULL) {
        found = true;
        fclose(fp);
        }
      }
    }

  //==== read .pm file as a command script =====//

  if (found) {
    string cmd = "read ";
    cmd.append(name);
    if (!silent) fprintf (stderr, ">>> %s \n", cmd.c_str());
    pm_CmdProcess (NULL, cmd);
    }
  }

//*============================================================*
//*==========                 main                   ==========*
//*============================================================*

int
main (int num_args, char **args) 
  {
  char *script;

  // initialize protein modeler //

  for (int j = 1; j < num_args; j++) {
    if (!strncmp(args[j], "silent", 6)) {
      silent = true;
      }
    }

  pmSystem.init(silent);

  // process .pm init file //

  procInitFile();

  // process arguments //

  procArgs(num_args, args, &script);

  // enter event processing loop //

  //fprintf (stderr,"##### script=%s \n", script);

  pmSystem.procEvents (script);

  return (0);
  }

