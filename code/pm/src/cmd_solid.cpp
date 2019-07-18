
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
//              s o l i d       c o m m a n d s              //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========         pm_CmdSolidGetAxes             ==========*
//*============================================================*

bool 
pm_CmdSolidGetAxes (PmCmdDataList& dlist, string dname, PmPcaResults& pca,
                    PmExtent& proj, float wd[3])
  { 
  PmAtomFilter filter;
  string region; 
  vector<PmVector3> coords;
  PmMolecule *domain;

  pmSystem.getDomain (dname, &domain);

  if (!domain) {
    pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
    return false;
    }

  if (!dlist.getString("region", region)) {
    pm_ErrorReport (PM, "no region specified.", "*");
    return false;
    }

  // get atoms names for computation.  //

  dlist.getStringList("atoms", filter.names);
  domain->getAtomCoords(region, filter, coords);
  pm_MathPrincipalComponents(coords, pca);
  pm_MathPrincipalComponentsProj(coords, pca, proj, wd);
  return true;
  }

//*============================================================*
//*==========         pm_CmdSolidDefRegion           ==========*
//*============================================================*
// process solid define region command.

void
pm_CmdSolidDefRegion (PmSolid *solid, string name, PmCmdDataList& dlist)
  {
  string desc, dv;
  PmCmdData data;
  PmSolidRegionParameters params;
  PmVector3 point;
  float distance;
  bool use_spheres, show;

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

  solid->defineRegion(name, desc, params);
  pm_PrintMsg (CMDPR, "region \"%s\" defined. ", name.c_str());

  if (dlist.getBoolean("show", show) && show) {
    PmVector3 color;

    if (!dlist.getVector("color", color)) {
      color.set(1,1,1);
      }

    dlist.getBoolean("use_spheres", use_spheres);
    solid->displayRegion(name, color, use_spheres);
    }
  }

//*============================================================*
//*==========         pm_CmdSolidXform               ==========*
//*============================================================*
// process solid xform command.

void
pm_CmdSolidXform (PmSolid *solid, PmCmdDataList& dlist)
  {
  string dv;
  PmCmdData data;
  PmVector3 translation, rotation, axis;
  PmXform xform;
  PmMassProperties props;
  float angle = 0.0;
  bool rot_about_axis = false;
  PmSolid *copy = NULL;

  // set center of xform  //
  solid->getMassProps(props);
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

  if (copy) {
    //copy->xformCoordinates(xform);
    }
  else {
    //solid->xformCoordinates(xform);
    }
  }

//*============================================================*
//*==========              pm_CmdSolidCreatEllipsoid ==========*
//*============================================================*
// process solid creat ellipsoid command.

void
pm_CmdSolidCreatEllipsoid (string name, PmCmdDataList& dlist)
  {
  string rstr, astr, pstr;
  PmCmdData pdata;
  float r;
  PmVector3 vec, pt1, pt2;
  vector<PmVector3> axes(3);
  vector<float> radii(3);
  PmVector3  center;
  PmEllipsoid *ellipsoid;

  // create an ellipsoid using two points //

  if (dlist.getString("point1", pstr)) {
    dlist.getData("point1", pdata);

    if (!pm_CmdProcPosition(pdata, pt1)) {
      pm_ErrorReport (PM, "couldn't process point1 data \"%s\".", "*", pstr.c_str());
      return;
      }

    if (dlist.getString("point2", pstr)) {
      dlist.getData("point2", pdata);

      if (!pm_CmdProcPosition(pdata, pt2)) {
        pm_ErrorReport (PM, "couldn't process point2 data \"%s\".", "*",
                        pstr.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no point2 specified.", "*");
      return;
      }

    if (!dlist.getFloat("radius", r)) {
      pm_ErrorReport (PM, "no radius specified.", "*");
      return;
      }

    center = 0.5*(pt1 + pt2);
    vec = pt2 - pt1;
    radii[0] = vec.length() / 2,0;
    radii[1] = r; 
    radii[2] = r;
    vec.normalize(); 

    axes[0] = vec;
    pm_MathBasisCompute (vec, axes[1], axes[2]);

    PmSolid::getEllipsoid(name, &ellipsoid);
    ellipsoid->setCenter(center);
    ellipsoid->setAxes(axes);
    ellipsoid->setRadii(radii);
    }

  // create using a point and an axis //

  else if (dlist.getString("point", pstr)) {
    dlist.getData("point", pdata);

    if (!pm_CmdProcPosition(pdata, pt1)) {
      pm_ErrorReport (PM, "couldn't process point data \"%s\".", "*", pstr.c_str());
      return;
      }

    if (!dlist.getVector("axis", vec)) {
      pm_ErrorReport (PM, "no axis given", "*");
      return;
      }

    if (dlist.getFloat("radius1", r)) {
      radii[0] = r;
      }
    else {
      pm_ErrorReport (PM, "no radius1 specified.", "*");
      return;
      }

    if (dlist.getFloat("radius2", r)) {
      radii[1] = r;
      radii[2] = r;
      }
    else {
      pm_ErrorReport (PM, "no radius2 specified.", "*");
      return;
      }

    vec.normalize();
    axes[0] = vec;
    pm_MathBasisCompute (vec, axes[1], axes[2]);
    pt2 = pt1 + 2.0*radii[0]*axes[0];
    center = 0.5*(pt1 + pt2);

    PmSolid::getEllipsoid(name, &ellipsoid);
    ellipsoid->setCenter(center);
    ellipsoid->setAxes(axes);
    ellipsoid->setRadii(radii);

    pm_PrintMsg (CMDPR, "pt1 = %g %g %g ", pt1[0], pt1[1], pt1[2]); 
    pm_PrintMsg (CMDPR, "pt2 = %g %g %g ", pt2[0], pt2[1], pt2[2]); 
    pm_PrintMsg (CMDPR, "center = %g %g %g ", center[0], center[1], center[2]); 
    }

  // uses axes and radii //

  else {
    astr = "axis1";
    rstr = "radius1";

    for (int i = 0; i < 3; i++) {
      astr[4] = '1' + i;
      rstr[6] = '1' + i;

      if (dlist.getVector(astr, vec)) {
        axes[i] = vec;
        }
      else {
        pm_ErrorReport (PM, "no %s given.", "*", astr.c_str());
        return;
        }

      if (dlist.getFloat(rstr, r)) {
        radii[i] = r;
        }
      else {
        pm_ErrorReport (PM, "no radius1 given.", "*");
        return;
        }
      }

    if (dlist.getString("center", pstr)) {
      dlist.getData("center", pdata);

      if (!pm_CmdProcPosition(pdata, center)) {
        pm_ErrorReport (PM, "couldn't process center data \"%s\".", "*",
                        pstr.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no center given.", "*");
      return;
      }

    PmSolid::getEllipsoid(name, &ellipsoid);
    ellipsoid->setCenter(center);
    ellipsoid->setAxes(axes);
    ellipsoid->setRadii(radii);
    }
  }

//*============================================================*
//*==========              pm_CmdSolidCreatSphere    ==========*
//*============================================================*
// process solid creat sphere command.

void
pm_CmdSolidCreatSphere (string name, PmCmdDataList& dlist)
  {
  float radius;
  PmCmdData pdata;
  string pstr, dname;
  PmVector3 vec, pt, center, pt1, pt2;
  PmSphere *sphere;
  PmSolid::getSphere(name, &sphere);

  //===== create sphere using two points =====//

  if (dlist.getString("point1", pstr)) {
    dlist.getData("point1", pdata);

    if (!pm_CmdProcPosition(pdata, pt1)) {
      pm_ErrorReport (PM, "couldn't process point1 data \"%s\".", "*", pstr.c_str());
      return;
      }

    if (dlist.getString("point2", pstr)) {
      dlist.getData("point2", pdata);

      if (!pm_CmdProcPosition(pdata, pt2)) {
        pm_ErrorReport (PM, "couldn't process point2 data \"%s\".", "*", 
                        pstr.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no point2 specified.", "*");
      return;
      }

    center = 0.5*(pt1 + pt2);
    vec = pt2 - pt1;
    radius = vec.length() / 2,0;
    sphere->setCenter(center);
    sphere->setRadius(radius);

    pm_PrintMsg (CMDPR, "pt1 = %g %g %g ", pt1[0], pt1[1], pt1[2]); 
    pm_PrintMsg (CMDPR, "pt2 = %g %g %g ", pt2[0], pt2[1], pt2[2]); 
    pm_PrintMsg (CMDPR, "center = %g %g %g ", center[0], center[1], center[2]); 
    }

  //===== create sphere from a domain =====//

  else if (dlist.getString("domain", dname)) {
    PmMolecule *domain;
    PmPcaResults pca;
    PmAtomFilter filter;
    string region, allign, allign_point;
    vector<PmVector3> coords;
    PmExtent proj;
    float wd[3], dist;

    if (!pm_CmdSolidGetAxes(dlist, dname, pca, proj, wd)) {
      return;
      }

    center = pca.com;
    radius = pca.s1;
    sphere->setCenter(center);
    sphere->setRadius(radius);
    }

  //===== specify parameters =====//

  else {
    if (dlist.getFloat("radius", radius)) {
      sphere->setRadius(radius);
      }
    else {
      pm_ErrorReport (PM, "no radius specified.", "*");
      return;
      }

    if (dlist.getVector("axis", vec)) {
      if (!dlist.getVector("point", pt)) {
        pm_ErrorReport (PM, "no point specified.", "*");
        return;
        }

      vec.normalize();
      center = pt + radius*vec;
      sphere->setCenter(center);
      }

    else if (dlist.getString("center", pstr)) {
      dlist.getData("center", pdata);

      if (!pm_CmdProcPosition(pdata, center)) {
        pm_ErrorReport (PM, "couldn't process center data \"%s\".", "*",
                        pstr.c_str());
        return;
        }

      sphere->setCenter(center);
      }
    else {
      pm_ErrorReport (PM, "no center given.", "*");
      return;
      }
    }

  pm_PrintMsg (CMDPR, "create sphere \"%s\" ", name.c_str());
  pm_PrintMsg (CMDSP, "center = %g %g %g ", center[0], center[1], center[2]); 
  pm_PrintMsg (CMDSP, "radius = %g ", radius); 
  }

//*============================================================*
//*==========              pm_CmdSolidCreatCyl       ==========*
//*============================================================*
// process solid create cylinder command.

void
pm_CmdSolidCreatCyl (string name, PmCmdDataList& dlist)
  {
  string dname, pstr;
  PmSolid *solid;
  PmVector3 center, axis, pt1, pt2; 
  float radius, length, val;
  PmVector3 vec; 
  PmCylinder *cyl;
  PmCmdData pdata;

  // create cyl using two points //

  if (dlist.getString("point1", pstr)) {
    PmCmdData pdata;
    dlist.getData("point1", pdata);

    if (!pm_CmdProcPosition(pdata, pt1)) {
      pm_ErrorReport (PM, "couldn't process point1 data \"%s\".", "*", pstr.c_str());
      return;
      }

    if (dlist.getString("point2", pstr)) {
      dlist.getData("point2", pdata);

      if (!pm_CmdProcPosition(pdata, pt2)) {
        pm_ErrorReport (PM, "couldn't process point2 data \"%s\".", "*", 
                        pstr.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no point2 given", "*");
      return;
      }

    if (!dlist.getFloat("radius", radius)) {
      pm_ErrorReport (PM, "no radius given", "*");
      return;
      }

    center = 0.5*(pt2 + pt1);
    axis = pt2 - pt1;
    length = axis.length();
    axis.normalize();

    PmSolid::getCylinder(name, &cyl);
    cyl->setAxis(axis);
    cyl->setCenter(center);
    cyl->setLength(length);
    cyl->setRadius(radius);
    }

  else if (dlist.getString("domain", dname)) {
    PmMolecule *domain;
    PmPcaResults pca;
    PmAtomFilter filter;
    string region, allign, allign_point;
    vector<PmVector3> coords;
    PmExtent proj;
    float wd[3], dist;
    PmCylinder *acyl;

    pmSystem.getDomain (dname, &domain);

    if (!domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
      return;
      }

    if (!dlist.getString("region", region)) {
      pm_ErrorReport (PM, "no region specified.", "*");
      return;
      }

    // get atoms names for computation.  //

    dlist.getStringList("atoms", filter.names);
    domain->getAtomCoords(region, filter, coords);
    pm_MathPrincipalComponents(coords, pca);
    pm_MathPrincipalComponentsProj(coords, pca, proj, wd);
    //domain->displayPca(pca, true);

    PmSolid::getCylinder(name, &cyl);

    if (dlist.getString("allign", allign)) {
      PmSolid::getCylinder(allign, &acyl);
      acyl->getParameters(center, radius, length, axis);

      if (!dlist.getFloat("allign_point", dist)) {
        dist = 1.0;
        }

      pt1 = center + 0.5*length*axis*dist;
      pt2 = pca.com + wd[0]*pca.axis1*dist;
      center = 0.5*(pt2 + pt1);
      axis = pt2 - pt1;
      length = axis.length();
      axis.normalize();
      cyl->setAxis(axis);
      cyl->setCenter(center);
      cyl->setLength(length);
      }
    else {
      cyl->setCenter(pca.com);
      cyl->setAxis(pca.axis1);
      length = 2.0*wd[0];
      cyl->setLength(length);
      }

    if (!dlist.getFloat("radius", radius)) {
      radius = (wd[1] + wd[2]) / 2.0;
      }

    cyl->setRadius(radius);
    }

  else {
    PmSolid::getCylinder(name, &cyl);

    if (dlist.getFloat("radius", val)) {
      cyl->setRadius(val);
      }
    else {
      pm_ErrorReport (PM, "no radius given.", "*");
      return;
      }

    if (dlist.getFloat("length", val)) {
      cyl->setLength(val);
      }

    if (dlist.getVector("axis", vec)) {
      cyl->setAxis(vec);
      }

    if (dlist.getString("center", pstr)) {
      dlist.getData("center", pdata);

      if (!pm_CmdProcPosition(pdata, center)) {
        pm_ErrorReport (PM, "couldn't process center data \"%s\".", "*",
                        pstr.c_str());
        return;
        }

      cyl->setCenter(center);
      }
    else {
      pm_ErrorReport (PM, "no center given.", "*");
      return;
      }
    }

  cyl->getParameters(center, radius, length, axis);
  pm_PrintMsg (CMDPR, "create cylinder \"%s\" ", name.c_str());
  pm_PrintMsg (CMDSP, "center = %g %g %g ", center[0], center[1], center[2]); 
  pm_PrintMsg (CMDSP, "radius = %g ", radius); 
  pm_PrintMsg (CMDSP, "length = %g ", length); 
  pm_PrintMsg (CMDSP, "axis = %g %g %g ", axis[0], axis[1], axis[2]); 

  pt1 = center - 0.5*length*axis;
  pt2 = center + 0.5*length*axis;
  pm_PrintMsg (CMDPR, "end pt1 = %g %g %g ", pt1[0], pt1[1], pt1[2]); 
  pm_PrintMsg (CMDPR, "end pt2 = %g %g %g ", pt2[0], pt2[1], pt2[2]); 
  }

//*============================================================*
//*==========              pm_CmdSolid               ==========*
//*============================================================*
// process solid command.                                 

void
pm_CmdSolid (PmCmdDataList& dlist)
  {

  string dn, dv, name, type_str, dname, pstr;
  PmSolid *solid;
  PmSolidType type;
  PmVector3 center, axis, pt1, pt2; 
  float radius, length; 
  PmExtent extent;
  PmVector3 color;
  PmCmdData data;

  // get the next data item //

  dlist.getNext (data);

  // create a solid //

  if (data.name== "create") {
    dlist.getString ("name", name);
    dlist.getString ("type", type_str);
    PmSolid::convSolidType(type_str, type);

    if (type == PM_SOLID_UNKNOWN) {
      pm_ErrorReport ("pm> ", "unknown solid type \"%s\".", "*", type_str.c_str());
      return;
      }

    solid = PmSolid::create(name, type);
    pmSystem.addSolid(solid);

    // get parameters for particular solid type //

    if (type == PM_SOLID_CYLINDER) {
      pm_CmdSolidCreatCyl(name, dlist);
      }

    else if (type == PM_SOLID_SPHERE) {
      pm_CmdSolidCreatSphere (name, dlist);
      }

    else if (type == PM_SOLID_ELLIPSOID) {
      pm_CmdSolidCreatEllipsoid (name, dlist);
      }

    // set extent //

    solid->getExtent(extent);

    pm_PrintMsg (CMDPR, "extent  min (%g %g %g)  max (%g %g %g) ",
                 extent.min[0], extent.min[1], extent.min[2],
                 extent.max[0], extent.max[1], extent.max[2]);

    //pmSystem.setExtent (extent);
    pmSystem.updateExtent (extent);

    PmMassProperties props;
    solid->getMassProps(props);
    }

  else {
    pmSystem.getSolid(data.name, &solid);
    }

  if (!solid) {
    pm_ErrorReport (PM, "no solid named \"%s\".", "*", data.name.c_str());
    return;
    }

  //pmSystem.setCurrentSurface (surf);
  bool show_set = true;
  bool show = true;

  while (dlist.getNext(data)) { 
    data.getString(dv);

    // color = [ <r:0-1> <g:0-1> <b:0-1> ] //

    if (data.name== "color") {
      data.getVector(color);
      solid->setColor(color);
      }

    // display = point | line | solid 

    else if (data.name== "display") {
      PmGeometryDisplayType dtype;
      PmGraphicsGeometry::convDisplayType (dv, dtype);

      if (dtype == PM_GEOMETRY_DISPLAY_UNKNOWN) {
        pm_ErrorReport (PM, "unknown display type \"%s\".", "*", dv.c_str());
        }
      else {
        solid->setDisplayType(dtype);
        }
      }

    // shading = none | flat | color | normal

    else if (data.name== "shading") {
      PmGeometryShadingType stype;
      PmGraphicsGeometry::convShadingType (dv, stype);

      if (stype == PM_GEOMETRY_SHADING_UNKNOWN) {
        pm_ErrorReport (PM, "unknown shading type \"%s\".", "*", dv.c_str());
        }
      else {
        solid->setShadingType(stype);
        }
      }

    else if (data.name == "show") {
      show = data.getBoolean();
      }

     // ===== define ===== //

    else if (data.name == "define") {
      dlist.getNext(data, dv);

      if (data.name == "region") {
        pm_CmdSolidDefRegion(solid, dv, dlist);
        }

      return;
      }

    // transform the solid //

    else if (data.name== "xform") {
      pm_CmdSolidXform (solid, dlist);
      }
    }

  if (show_set) {
    solid->display(show);
    }
  }

