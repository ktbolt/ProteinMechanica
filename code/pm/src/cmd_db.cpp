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
//*==========              pm_CmdDatabase            ==========*
//*============================================================*
// process database command.                                 

void
pm_CmdDatabase (PmCmdDataList& dlist)
  {
  string dv, fname, name, format_str, data_type;
  PmDbInterfaceSelect db_sel;
  PmExtent extent;
  PmCmdData data;
  PmDbType format;

  while (dlist.getNext(data,dv)) {
    if (data.name == "open") {
      if (!dlist.getString("name", name)) { 
        pm_ErrorReport (PM, "no database name given.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      if (!dlist.getString("format", format_str)) {
        format = PM_DB_PDB;
        }
      else {
        format = pm_DbFormatConv (format_str);

        if (format == PM_DB_UNKNOWN) {
          pm_ErrorReport (PM, "uknown database format \"%s\".", "*", format_str.c_str());
          pmSystem.setCmdError(true);
          return;
          }
        }

      if (!dlist.getString ("file", fname)) {
        pm_ErrorReport (PM, "no file name given.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      dlist.getString ("data", data_type);
      PmDbInterface *db = db_sel.create(name, format);
      db->open(fname, PM_DB_MODE_READ, data_type);
      pmSystem.addDatabase(db);
      }
    }
  }

