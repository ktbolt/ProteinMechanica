
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
//              s u r f a c e   c o m m a n d s              //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========         pm_CmdSurfaceCreate            ==========*
//*============================================================*
// process surface create command.

void
pm_CmdSurfaceCreate (string name, PmCmdDataList& dlist, PmSurface **p_surf)
  {
  string desc, dv, type;
  PmCmdData data;
  PmSurface *surf;
  PmExtent extent;
  bool plate, plates;
  float thickness;

  *p_surf = NULL;

  dlist.getString ("name", name);
  dlist.getString ("type", type);
  dlist.getBoolean("plate", plate);
  dlist.getBoolean("plates", plates);

  if (!dlist.getFloat("thickness", thickness)) {
    thickness = 1.0;
    }

  if (plate) {
    PmVector3 point;
    vector<PmVector3> coords;
    float val;
    vector<string> points;
    string spt; 
    int i, j;

    if (dlist.getStringList("points", points)) {
      j = 0;

      for (int i = 0; i < points.size(); i++) {
        spt = points[i];

        if (convToFloat(spt, val)) {
          point[j++] = val;

          if (j == 3) {
            coords.push_back(point);
            j = 0;
            }
          }
        }
      }

    else {
      if (dlist.getVector("point1", point)) {
        coords.push_back(point);
        }

      if (dlist.getVector("point2", point)) {
        coords.push_back(point);
        }

      if (dlist.getVector("point3", point)) {
        coords.push_back(point);
        }
      }

    surf = new PmSurface(name, coords, thickness);
    }

  else if (plates) {
    vector<PmSurface*> surfs;
    PmSurface::genPlates(name, type, surfs);

    for (int i = 0; i < surfs.size(); i++) {
      surf = surfs[i];
      pmSystem.addSurface (surf);
      surf->getExtent(extent);
      pmSystem.setExtent (extent);
      surf->display(true);
      }

    return;
    }
  else {
    surf = new PmSurface(name, type);
    }

  if (!surf) {
    *p_surf = NULL;
    return;
    }

  // add surface to system //

  pmSystem.addSurface (surf);

  // set extent //

  surf->getExtent(extent);
  pmSystem.setExtent (extent);
  *p_surf = surf;
  }

//*============================================================*
//*==========         pm_CmdSurfaceDefRegion         ==========*
//*============================================================*
// process surface region command.

void
pm_CmdSurfaceDefRegion (PmSurface *surf, string name, PmCmdDataList& dlist)
  {
  string desc, dv, str;
  vector<string> strlist;
  PmCmdData data;
  PmSurfRegionParameters params;
  PmVector3 point;
  float distance, fval;
  bool use_spheres, show;
  float charge_range[2];

  if (!dlist.getString("indices", desc)) {
    if (dlist.getVector("point", point)) {
      if (!dlist.getFloat("distance", distance)) {
        pm_ErrorReport (PM, "distance not given", "*");
        return;
        }

      params.use_distance = true;
      params.distance = distance;
      params.point = point;
      }
    }

  if (dlist.getFloat("radius", fval)) {
    params.radius = fval;
    params.radius_set = true;
    }

  dlist.getBoolean("use_face", params.use_face);
  dlist.getString("data", params.data_name);

  // process charge range //

  if (dlist.getStringList("charge_range", strlist)) {
    if (strlist.size() != 2) {
      pm_ErrorReport (PM, "charge range should be two values.", "*");
      return;
      }

    for (unsigned int i = 0; i < strlist.size(); i++) {
      str = strlist[i];

      if (convToFloat(str, fval)) {
        charge_range[i] = fval;
        }
      else {
        pm_ErrorReport (PM, "improper value \"%s\" for charge range.", "*", 
                        str.c_str());
        return;
        }
      }

    params.charge_min = fmin(charge_range[0], charge_range[1]);
    params.charge_max = fmax(charge_range[0], charge_range[1]);
    }

  // create a region for the surface //

  surf->defineRegion(name, desc, params);

  if (dlist.getBoolean("show", show) && show) {
    PmVector3 color;

    if (!dlist.getVector("color", color)) {
      color.set(1,1,1);
      }

    dlist.getBoolean("use_spheres", use_spheres);
    surf->displayRegion(name, color, use_spheres);
    }
  }

//*============================================================*
//*==========         pm_CmdSurfaceXform             ==========*
//*============================================================*
// process surface xform command.

void
pm_CmdSurfaceXform (PmSurface *surface, PmCmdDataList& dlist)
  {
  string dv;
  PmCmdData data;
  PmVector3 translation, rotation, axis;
  PmXform xform;
  PmMassProperties props;
  float angle = 0.0;
  bool rot_about_axis = false;
  PmSurface *surf_copy = NULL;

  // set center of xform  //
  surface->getMassProps(props);
  xform.center = props.com;

  while (dlist.getNext (data, dv)) {
    if (data.name == "angle") {
      angle = data.getFloat();
      }

    else if (data.name == "axis") {
      data.getVector(axis);
      rot_about_axis = true;
      }

    else if (data.name == "copy") {
      //domain->copy(dv, &dom_copy);
      //pmSystem.addDomain(dom_copy);
      }

    else if (data.name == "translation") {
      data.getVector(xform.translation);
      }

    else if (data.name == "rotation") {
      PmVector3 angles;
      data.getVector(angles);
      xform.setRotation(angles[0], angles[1], angles[2]);
      }
    }

  // apply the transformation //

  if (rot_about_axis) {
    xform.setRotation(angle, axis);
    }

  if (surf_copy) {
    surf_copy->xformCoordinates(xform);
    }
  else {
    surface->xformCoordinates(xform);
    }
  }

//*============================================================*
//*==========              pm_CmdSurface             ==========*
//*============================================================*
// process surface command.                                 

void
pm_CmdSurface (PmCmdDataList& dlist)
  {
  string dn, dv, name;
  vector<string> strlist; 
  PmSurface *surf;
  PmExtent extent;
  PmVector3 color;
  PmCmdData data;

  dlist.getNext (data);

  if (data.name== "read") {
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
    surf = db->getSurface();

    if (!surf) {
      return;
      }

    // add surface to system //

    pmSystem.addSurface (surf);

    // set extent //

    surf->getExtent(extent);
    pmSystem.setExtent (extent);
    }

  else if (data.name== "create") {
    pm_CmdSurfaceCreate (dv, dlist, &surf);
    if (!surf) return;
    }

  else {
    pmSystem.getSurface (data.name, &surf);
    }

  if (!surf) {
    pm_ErrorReport (PM, "no surface named \"%s\".", "*", data.name.c_str());
    return;
    }

  pmSystem.setCurrentSurface (surf);
  bool show_set = true;
  bool show = true;

  while (dlist.getNext(data, dv)) { 

    // color = [ <r:0-1> <g:0-1> <b:0-1> ] //

    if (data.name == "color") {
      if (dv == "charge") {
        surf->setMapCharge(true);
        }
      else {
        data.getVector(color);
        surf->setColor(color);
        }
      }

    else if (data.name == "charge_range") {
      data.getStringList(strlist);
      string str;
      float fval, charge_range[2], charge_min, charge_max;

      if (strlist.size() != 2) {
        pm_ErrorReport (PM, "charge range should be two values.", "*");
        return;
        }

      for (unsigned int i = 0; i < strlist.size(); i++) {
        str = strlist[i];

        if (convToFloat(str, fval)) {
          charge_range[i] = fval;
          }
        else {
          pm_ErrorReport (PM, "improper value \"%s\" for charge range.", "*",
                          str.c_str());
          return;
          }
        }

      charge_min = fmin(charge_range[0], charge_range[1]);
      charge_max = fmax(charge_range[0], charge_range[1]);
      surf->setChargeMap(charge_min, charge_max);
      }

    // display = point | line | solid // 

    else if (data.name == "display") {
      PmGeometryDisplayType dtype;
      PmGraphicsGeometry::convDisplayType (dv, dtype);

      if (dtype == PM_GEOMETRY_DISPLAY_UNKNOWN) {
        pm_ErrorReport (PM, "unknown display type \"%s\".", "*", dv.c_str());
        }
      else {
        surf->setDisplayType(dtype);
        }
      }

    else if (data.name == "lighting") {
      if (dv == "two-sided") {
        surf->setLighting(true);
        }
      }

    // shading = none | flat | color | normal //

    else if (data.name== "shading") {
      PmGeometryShadingType stype;
      PmGraphicsGeometry::convShadingType (dv, stype);

      if (stype == PM_GEOMETRY_SHADING_UNKNOWN) {
        pm_ErrorReport (PM, "unknown shading type \"%s\".", "*", dv.c_str());
        }
      else {
        surf->setShadingType(stype);
        }
      }

    // define //

    else if (data.name== "define") {
      dlist.getNext(data, dv);

      if (data.name == "region") {
        pm_CmdSurfaceDefRegion(surf, dv, dlist);
        }

      return;
      }

    // transform the surface //

    else if (data.name== "xform") {
      pm_CmdSurfaceXform (surf, dlist);
      }
    }

  if (show_set) {
    surf->display(show);
    }
  }

