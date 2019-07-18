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
//              v e c t o r   c o m m a n d s                //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========             pm_CmdVector               ==========*
//*============================================================*
// process vector command.

void
pm_CmdVector (PmCmdDataList& dlist)
  {
  PmCmdData data, pdata;
  string dv, name, bname, str; 
  PmVector3 dir, origin, pt1, pt2, pt3;
  float scale, width, angle;
  PmGeometry *geom;
  PmGeometryVector *vgeom;
  PmBody *body;
  PmVector3 v, v1, v2, v3, color;
  bool use_pca, negate, cross, proj_plane;
  int num_pts = 0;

  dlist.getNext (data);

  // create a new vector and add it to the system //

  if (data.name == "create") {
    dlist.getString("name", name);
    dlist.getBoolean("use_pca", use_pca); 
    dlist.getBoolean("negate", negate); 

    // use pca axes from a domain //

    if (use_pca) {
      PmMolecule *domain;
      string desc, dname;
      PmAtomFilter filter;
      PmPcaResults pca;
      PmVector3 axis, axis1, axis2, axis3;
      PmVector3 pt; 
      vector<PmVector3> coords;
      int pca_axis;
      float val, wd[3];
      PmExtent proj;

      if (dlist.getString("domain", dname)) { 
        pmSystem.getDomain (dname, &domain);

        if (!domain) {
          pm_ErrorReport (PM, "no pca domain named \"%s\".", "*", dname.c_str());
          return;
          }
        }
      else {
        pm_ErrorReport (PM, "no pca domain given.", "*");
        return;
        }

      if (!dlist.getString("pca_res", desc)) { 
        pm_ErrorReport (PM, "no pca residues given.", "*");
        return;
        }

      dlist.getStringList("atom_names", filter.names);

      domain->getAtomCoords(desc, filter, coords);

      if (coords.size() == 0) {
        pm_ErrorReport (PM, "no coordinates for pca.", "*");
        return;
        }

      pm_MathPrincipalComponents(coords, pca);
      pm_MathPrincipalComponentsProj(coords, pca, proj, wd);

      if (!dlist.getInt("pca_axis", pca_axis)) { 
        pca_axis = 1;
        }

      if (pca_axis == 1) {
        dir = pca.axis1;
        }
      else if (pca_axis == 2) {
        dir = pca.axis2;
        }
      else if (pca_axis == 3) {
        dir = pca.axis3;
        }
      else {
        pm_ErrorReport (PM, "unknown axis %d specified for pca.", "*", pca_axis);
        return;
        }

      if (negate) {
        dir = -dir;
        }

      if (dlist.getFloat("origin", val)) {
        origin = pca.com - 0.5*dir + val*dir;
        }
      else if (dlist.getVector("origin", v)) {
        origin = v;
        }
      else if (dlist.getData("origin", pdata)) {
        if (!pm_CmdProcPosition(pdata, origin)) {
          pm_ErrorReport (PM, "couldn't process vector origin.", "*");
          return;
          }
        }
      else {
        origin = pca.com - wd[pca_axis-1]*dir;
        pt = pca.com + wd[pca_axis-1]*dir;
        dir = pt - origin;
        }
      }

    // specify points directly //

    else {

      if (dlist.getVector("direction", v)) {
        dir = v;
        }

      if (dlist.getVector("origin", v)) {
        origin = v;
        }
      else if (dlist.getData("origin", pdata)) {
        if (!pm_CmdProcPosition(pdata, origin)) {
          pm_ErrorReport (PM, "couldn't process vector origin.", "*");
          return;
          }
        }

      if (dlist.getVector("point1", v)) {
        pt1 = v;
        num_pts = 1;
        }
      else if (dlist.getData("point1", pdata)) {
        if (!pm_CmdProcPosition(pdata, pt1)) {
          pm_ErrorReport (PM, "couldn't process vector point1.", "*");
          return;
          }

        num_pts = 1;
        }

      if (dlist.getVector("point2", v)) {
        pt2 = v;
        num_pts += 1;
        }
      else if (dlist.getData("point2", pdata)) {
        if (!pm_CmdProcPosition(pdata, pt2)) {
          pm_ErrorReport (PM, "couldn't process vector point2.", "*");
          return;
          }

        num_pts += 1;
        }

      if (dlist.getVector("point3", v)) {
        pt3 = v;
        num_pts += 1;
        }
      else if (dlist.getData("point3", pdata)) {
        if (!pm_CmdProcPosition(pdata, pt3)) {
          pm_ErrorReport (PM, "couldn't process vector point 3.", "*");
          return;
          }

        num_pts += 1;
        }

      //===== define 2nd point as projection onto a vector =====//

      if (dlist.getString("project", str)) {
        PmGeometryVector *vgeom;
        PmVector3 line[2], proj_pt;
        float dist; 
        vector<PmVector3> pts;
        pmSystem.getGeometry(str, &geom);
        vgeom = dynamic_cast<PmGeometryVector*>(geom);

        if (!vgeom) {
          pm_ErrorReport (PM, "no vector named \"%s\".", "*", dv.c_str());
          return;
          }

        if (num_pts != 1) {
          pm_ErrorReport (PM, "no point1 given.", "*");
          return;
          }

        vgeom->getCurrVector(v1);
        vgeom->getCurrPoints(pts);
        line[0] = pts[0];
        line[1] = pts[1];
        pm_GcLinePointProj (line, pt1, dist, pt2);
        num_pts += 1;
        fprintf (stderr, ">>> vector %s point1=(%f %f %f) point2=(%f %f %f) \n", 
                 name.c_str(), pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2]); 
        }

      //===== define 2nd point as projection onto a plane =====//

      else if (dlist.getBoolean("project_plane", proj_plane)) {
        PmGeometryVector *vgeom;
        PmVector3 plane_pt, plane_normal, line[2], proj_pt;
        float pd, dist;
        vector<PmVector3> pts;

        if (num_pts != 1) {
          pm_ErrorReport (PM, "no point1 given.", "*");
          return;
          }

        if (!dlist.getVector("plane_point", plane_pt)) {
          if (dlist.getData("plane_point", pdata)) {
            if (!pm_CmdProcPosition(pdata, plane_pt)) {
              pm_ErrorReport (PM, "couldn't process vector plane point.", "*");
              return;
              }
            }
          }

        if (!dlist.getVector("plane_normal", plane_normal)) {
          pm_ErrorReport (PM, "no plane normal given.", "*");
          return;
          }

        pd = plane_normal*plane_pt;
        dist = pt1*plane_normal - pd;
        pt2 = pt1 - plane_normal*dist;
        num_pts += 1;
        fprintf (stderr, ">>> vector point2=(%f %f %f) \n", pt2[0], pt2[1], pt2[2]); 
        }

      if (num_pts == 2) { 
        origin = pt1;
        dir = pt2 - pt1;
        }

      else if (num_pts == 3) { 
        origin = pt1;
        v1 = pt2 - pt1;
        v2 = pt3 - pt1;
        dir = v1.cross(v2);
        float mag = dir.length();

        if (mag == 0.0) {
          pm_ErrorReport (PM, "points are colinear.", "*");
          return;
          }

        dir = dir / mag;
        }

      if (negate) {
        dir = -dir;
        }
      }

    pm_PrintMsg (CMDPR, "vector \"%s\" added.", name.c_str());
    pm_PrintMsg (CMDPR, "   origin = %f %f %f", origin[0], origin[1], origin[2]); 
    pm_PrintMsg (CMDPR, "   dir    = %f %f %f", dir[0], dir[1], dir[2]); 

    vgeom = new PmGeometryVector(name, origin, dir);
    geom = dynamic_cast<PmGeometry*>(vgeom);
    pmSystem.addGeometry(geom);
    }

  else {
    pmSystem.getGeometry(data.name, &geom);
    vgeom = dynamic_cast<PmGeometryVector*>(geom);
    }

  if (!vgeom) {
    pm_ErrorReport (PM, "no vector named \"%s\".", "*", data.name.c_str());
    return;
    }

  // process vector parameters //

  bool show = false, show_set = false;

  while (dlist.getNext(data)) {
    data.getString(dv);

    if (data.name == "body") {
      pmSystem.getBody (dv, &body);

      if (!body) {
        pm_ErrorReport (PM, "no body named \"%s\". ", "*", dv.c_str());
        return;
        }

      vgeom->setBody(body);
      //body->setGeometry(geom);
      }

    // color = [ <r:0-1> <g:0-1> <b:0-1> ] //

    else if (data.name== "color") {
      data.getVector(color);
      vgeom->setColor(color);
      }

    else if (data.name== "angle") {
      PmVector3 v1, v2;
      PmGeometryVector *vgeom1;
      string name, name1;

      vgeom->getCurrVector(v1);
      pmSystem.getGeometry(dv, &geom);
      vgeom1 = dynamic_cast<PmGeometryVector*>(geom);

      if (!vgeom1) {
        pm_ErrorReport (PM, "no vector named \"%s\".", "*", dv.c_str());
        return;
        }

      // compute angle using cross product //

      dlist.getBoolean("cross", cross); 

      vgeom1->getCurrVector(v2);
      v1.normalize();
      v2.normalize();
 
      if (cross) {
        PmVector3 v;
        v = v1.cross(v2);
        angle = 180.0*asin(v.length()) / M_PI;
        }
      else {
        angle = 180.0*acos(v1*v2) / M_PI;
        }

      vgeom->getName(name);
      vgeom1->getName(name1);
      fprintf (stderr, "\n>>> angle between vectors \"%s\" and \"%s\" = %f degs \n", 
               name.c_str(), name1.c_str(), angle);
      }

    else if (data.name == "print") {
      PmVector3 v;
      string name;
      vgeom->getName(name);
      vgeom->getCurrVector(v);
      fprintf (stderr, "\n>>> vector \"%s\" value=%f %f %f length=%g \n", name.c_str(),
               v[0], v[1], v[2], v.length());
      }

    else if (data.name== "scale") {
      scale = data.getFloat();
      vgeom->setScale(scale);
      }

    else if (data.name == "show") {
      show = data.getBoolean();
      show_set = true;
      }

    else if (data.name== "width") {
      width = data.getFloat();
      vgeom->setWidth(width);
      }

    else if ((data.name == "write") || (data.name == "append")) {
      vector<PmVector3> pts;
      string fname, format, color;
      FILE *fp;

      vgeom->getCurrPoints(pts);
      dlist.getString("file", fname);
      dlist.getString("format", format);

      if (data.name == "write") {
        fp = fopen(fname.c_str(), "w");
        }
      else {
        fp = fopen(fname.c_str(), "a");
        }

      if (format == "vmd") {
        if (dlist.getString("vmd_color", color)) {
          fprintf (fp, "draw color %s \n", color.c_str());
          }

        fprintf (fp, "draw cylinder \"%f %f %f\"  \"%f %f %f\" radius 1.0 \n",
                 10.0*pts[0][0], 10.0*pts[0][1], 10.0*pts[0][2], 
                 10.0*pts[1][0], 10.0*pts[1][1], 10.0*pts[1][2]); 
        }

      fflush(fp);
      fclose (fp);
      }
    }

  if (show_set) {
    vgeom->display(show);
    }
  }

