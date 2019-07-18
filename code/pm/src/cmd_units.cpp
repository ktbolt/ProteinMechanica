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
//*==========              pm_CmdUnits               ==========*
//*============================================================*
// process units command.

void PM_EXPORT
pm_CmdUnits (PmCmdDataList& dlist)
  {

  string dv, uname, name;

  PmUnitType type;

  float scale;

  PmCmdData data;

  dlist.getNext (data, dv);

  if (data.name == "") {
    pmSystem.printUnits();
    return;
    }

  //fprintf (stderr, ">>> pm_CmdUnits:  data.name[%s] \n", data.name.c_str());

  if (data.name == "scale") {
    if (!dlist.getString("name", name)) {
      return;
      }

    pmSystem.convUnitType(name, type);

    if (type == PM_UNIT_UNKNOWN) {
      pm_ErrorReport (PM, "unknown unit name \"%s\".", "*", name.c_str());
      return;
      }

    if (!dlist.getString("uname", uname)) {
      return;
      }

    if (!dlist.getFloat("value", scale)) {
      return;
      }

    pmSystem.setUnitsScale (type, uname, scale);
    }
  }
