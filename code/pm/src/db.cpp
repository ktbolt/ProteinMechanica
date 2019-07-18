
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
// * db:                 d a t a b a s e                       *
//*============================================================*

#include "db.h"

namespace ProteinMechanica {

static int pm_DbLineGet (const FILE *fp, char& line);
static void pm_DbLineStrip (char *s);

#define DB_MAX_LINE 3000

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

//PmDbInterface::PmDbInterface(const string name) {} 

/*
PmDbInterface::PmDbInterface(const string name) { 
  this->name = name;
  }
*/

PmDbInterface::~PmDbInterface() { }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// creates interface object defined by <type>.

PmDbInterface *
PmDbInterfaceSelect::create (const string name, const PmDbType type)
  {
  PmDbInterface *ip;

  switch (type) { 
    case PM_DB_CERFACS:
      ip = new PmDbCerfacsInterface(name);
    break;

    case PM_DB_GROMACS:
      ip = new PmDbGroInterface(name);
    break;

    case PM_DB_MRC:
      ip = new PmDbMrcInterface(name);
    break;

    case PM_DB_PDB:
      ip = new PmDbPdbInterface(name);
    break;

    case PM_DB_PM:
      ip = new PmDbPmInterface(name);
    break;

    default:
      ip = NULL; 
    }

  return (ip);
  }

//////////////////////////////////////////////////////////////////
//                    p u b l i c                              //
////////////////////////////////////////////////////////////////

//*============================================================*
//*==========              pm_DbFormatConv           ==========*
//*============================================================*
// convert format string to type. 

PmDbType
pm_DbFormatConv (const string& format) 
  {
  PmDbType type = PM_DB_UNKNOWN;

  if (format == "pdb") {
    type = PM_DB_PDB;
    }
  else if (format == "cerfacs") {
    type = PM_DB_CERFACS;
    }
  else if (format == "mrc") {
    type = PM_DB_MRC;
    }
  else if (format == "gromacs") {
    type = PM_DB_GROMACS;
    }
  else if (format == "pm") {
    type = PM_DB_PM;
    }
  else if (format == "text") {
    type = PM_DB_TEXT;
    }
  else if (format == "coord") {
    type = PM_DB_COORD;
    }

  return (type);
  }

//*============================================================*
//*==========              convModeType              ==========*
//*============================================================*
// convert mode type string to type.

PmDbModeType
PmDbInterface::convModeType(const string mode) 
  {
  PmDbModeType type;

  if (mode == "read") {
    type = PM_DB_MODE_READ;
    }
  else if (mode == "write") {
    type = PM_DB_MODE_WRITE;
    }
  else if (mode == "append") {
    type = PM_DB_MODE_APPEND;
    }
  else {
    type = PM_DB_MODE_UNKNOWN;
    }

  return (type);
  }

//*============================================================*
//*==========              pm_DbLineGet              ==========*
//*============================================================*
// pm_DbLineGet 

static int
pm_DbLineGet (FILE *fp, char line[]) {

  while (1) {
    if (fgets(line, 2000, fp) == NULL) {
      return (0);
      }
    else {
      pm_DbLineStrip (line);

      if (line[0] != '\n') {
        return (1);
        }
      }
    }
  }

//*============================================================*
//*==========              pm_DbLineStrip            ==========*
//*============================================================*
// remove blanks from a line.

static void
pm_DbLineStrip (char *s)
  {

  unsigned int i, n;

  for (i = 0, n = 0; i < strlen(s); i++) {
    if (n != 0) {
      s[n++] = s[i];
      }
    else if (s[i] != ' ') {
      s[n++] = s[i];
      }
    }

  s[n] = '\0';
  }

//*============================================================*
//*==========              pm_DbStrParse             ==========*
//*============================================================*

/*
static void
pm_DbStrParse (char *line, char *val) 
  {
  for (unsigned int i = 0; i < strlen(line); ++i) {
    if (line[i] == '"') {
      break;
      }
    }
  }
*/

//*============================================================*
//*==========              pm_DbRecordStrGet         ==========*
//*============================================================*

static void
pm_DbRecordStrGet (char *line, int start, int end, char *str) {

  int i, j;

  start--;
  end--;

  for (j = 0, i = start; i <= end; i++, j++) {
    str[j] = line[i];
    }

  if (start != end) {
    str[j] = '\0';
    }
  }

//*============================================================*
//*==========              pm_DbRecordIntGet         ==========*
//*============================================================*

static void
pm_DbRecordIntGet (char *line, int start, int end, int& val) {
  char str[80];
  pm_DbRecordStrGet (line, start, end, str);
  val = atoi (str);
  }

//*============================================================*
//*==========              pm_DbRecordIntGet         ==========*
//*============================================================*

static void
pm_DbRecordRealGet (char *line, int start, int end, float& val) {
  char str[80];
  pm_DbRecordStrGet (line, start, end, str);
  val = atof (str);
  }


////////////////////////////////////////////////////////////////
//                    include db interfaces                  //
//////////////////////////////////////////////////////////////

#include "db_cerfacs.cpp"
#include "db_gro.cpp"
#include "db_mrc.cpp"
#include "db_pdb.cpp"
#include "db_pm.cpp"

}
