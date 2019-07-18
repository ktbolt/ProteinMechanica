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
//              c u r v e     c o m m a n d s                //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========             pm_CmdCurve                ==========*
//*============================================================*
// process curve command.

void
pm_CmdCurve (PmCmdDataList& dlist)
  {
  PmCmdData data, pdata;
  string dv, name, str; 
  vector<string> spts; 
  float width, fval; 
  PmGeometry *geom;
  PmGeometryCurve *cgeom;
  PmBody *body;
  PmVector3 color, point;
  vector<PmVector3> points;
  int num_pts = 0;

  dlist.getNext (data);

  // create a new curve and add it to the system //

  if (data.name == "create") {
    dlist.getString("name", name);
    dlist.getStringList("points", spts); 
    int j = 0;

    for (int i = 0; i < spts.size(); i++) {
      str = spts[i];

      if (convToFloat(str, fval)) {
        point[j++] = fval;

        if (j == 3) {
          points.push_back(point);
          j = 0;
          }
        }
      }

    cgeom = new PmGeometryCurve(name, points);
    geom = dynamic_cast<PmGeometry*>(cgeom);
    pmSystem.addGeometry(geom);
    }

  else {
    pmSystem.getGeometry(data.name, &geom);
    cgeom = dynamic_cast<PmGeometryCurve*>(geom);
    }

  if (!cgeom) {
    pm_ErrorReport (PM, "no curve named \"%s\".", "*", data.name.c_str());
    return;
    }

  //===== process curve parameters =====//

  bool show = false, show_set = false;

  while (dlist.getNext(data)) {
    data.getString(dv);

    if (data.name == "body") {
      pmSystem.getBody (dv, &body);

      if (!body) {
        pm_ErrorReport (PM, "no body named \"%s\". ", "*", dv.c_str());
        return;
        }

      cgeom->setBody(body);
      //body->setGeometry(geom);
      }

    // color = [ <r:0-1> <g:0-1> <b:0-1> ] //

    else if (data.name== "color") {
      data.getVector(color);
      cgeom->setColor(color);
      }

    else if (data.name == "show") {
      show = data.getBoolean();
      show_set = true;
      }

    else if (data.name== "width") {
      width = data.getFloat();
      cgeom->setWidth(width);
      }
    }

  if (show_set) {
    cgeom->display(show);
    }
  }

