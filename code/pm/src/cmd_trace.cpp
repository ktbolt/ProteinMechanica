
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

////////////////////////////////////////////////////////////////
//              t r a c e   c o m m a n d s                  //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========             pm_CmdTrace                ==========*
//*============================================================*
// process trace command.

void
pm_CmdTrace (PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, dname; 
  PmTrace *trace;
  PmExtent extent;
  PmVector3 color, dir;
  int ipt;

  dlist.getNext (data);

  // create a new force and add it to the system //

  if (data.name == "read") {
    string fname, format_str;
    PmDbInterfaceSelect db_sel;
    PmDbType format;

    dlist.getString ("name", name);
    dlist.getString ("format", format_str);
    dlist.getString ("file", fname);
    format = pm_DbFormatConv (format_str);
    PmDbInterface *db = db_sel.create ("tmp", format);

    if (!db) {
      pm_ErrorReport ("pm> ", "unknown format \"%s\".", "*", format_str.c_str());
      return;
      }

    db->open(fname, PM_DB_MODE_READ, "");
    db->read(name);
    db->close();
    trace = db->getTrace();

    if (!trace) {
      return;
      }

    pmSystem.addTrace(trace);

    // set extent //
    trace->getExtent(extent);
    pmSystem.setExtent(extent);
    return;
    }
  else {
    pmSystem.getTrace(data.name, &trace);
    }

  if (!trace) {
    pm_ErrorReport (PM, "no trace named \"%s\".", "*", data.name.c_str());
    return;
    }

  bool show_set = true;
  bool show = true;

  while (dlist.getNext(data)) {
    data.getString(dv);

    // color = [ <r:0-1> <g:0-1> <b:0-1> ] //

    if (data.name == "color") {
      data.getVector(color);
      trace->setColor(color);
      }

    // find an interval within a trace //

    else if (data.name == "interval") {
      float t, tol, min_dl, max_dl;
      tol = data.getFloat();
      dlist.getBoolean("show", show);
      trace->findInterval(tol, show, ipt, t, min_dl, max_dl);
      fprintf (stderr, "%s ipt = %d  min_dl = %f   max_dl = %f \n", CMDPR, ipt, min_dl,
               max_dl);
      }

    else if (data.name == "length") {
      dlist.getVector("direction", dir);

      if (!dlist.getInt("point", ipt)) {
        ipt = -1;
        }

      float length = trace->compLength(dir, ipt);
      fprintf (stderr, "%s length = %f \n", CMDPR, length);
      }

    }

  if (show_set) {
    trace->display(show);
    }
  }

