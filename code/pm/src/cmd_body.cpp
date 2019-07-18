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
//              b o d y     c o m m a n d s                  //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========           pm_CmdBodyContact            ==========*
//*============================================================*
// process body contact command.

void
pm_CmdBodyContact(PmBody *body, PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, bname;
  PmVector3 color1, color2;
  bool flag, show, show_set;
  vector<string> body_names;
  PmBody *bdy;
  vector<PmBody*> bodies;
  float tol = 0.0;

  show_set = true;
  show = true;

  if (!dlist.getStringList("list", body_names)) {
    pm_ErrorReport (PM, "no body names given.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  for (unsigned int i = 0; i < body_names.size(); i++) {
    bname = body_names[i]; 
    pmSystem.getBody(bname, &bdy);

    if (!bdy) {
      pm_ErrorReport (PM, "no body named \"%s\".", bname.c_str());
      pmSystem.setCmdError(true);
      return;
      }
    else {
      bodies.push_back(bdy);
      }
    }

  if (dlist.getBoolean("show", flag)) {
    show_set = true;
    show = flag;
    }

  if (!dlist.getFloat("tolerance", tol)) {
    tol = 0.0;
    }

  if (!dlist.getVector("color1", color1)) {
    color1.set(1,0,0);
    }

  if (!dlist.getVector("color2", color2)) {
    color2.set(0,1,0);
    }

  for (unsigned int i = 0; i < bodies.size(); i++) {
    bdy = bodies[i];
    body->checkContact(bdy, tol, color1, color2, show);
    }
  }

//*============================================================*
//*==========     pm_CmdBodyAddPotentialGeometry     ==========*
//*============================================================*
// process body add potential geometry command.

void
pm_CmdBodyAddPotentialGeometry (PmBody *body, string name, PmCmdDataList& dlist)
  {
  PmCmdData data, pdata;
  string str, rgn_name, bname; 
  PmPotentialType ptype;
  PmPotentialGeomType pgeom_type;
  bool show, sidechain, display_spheres;
  PmPotentialParameters params;
  PmGeometryDisplayType dtype;
  PmVector3 color;
  float width, radius, charge;

  //===== process potentail geometry type =====//

  if (dlist.getString("geometry", str)) {
    PmPotentialGeom::convType(str, pgeom_type);

    if (pgeom_type == PM_POTENTIAL_GEOM_UNKNOWN) {
      pm_ErrorReport (PM, "unknown potential geometry type \"%s\".", "*", str.c_str());
      pmSystem.setCmdError(true);
      return;
      }
    }
  else {
    pm_ErrorReport (PM, "no potential geometry given.", "*");
    pmSystem.setCmdError(true);
    return;
    }

 //===== process body region =====//

 if (dlist.getString("region", str)) {
   if (body->hasRegion(str)) {
     rgn_name = str;
     }
   else { 
     pm_ErrorReport (PM, "no region \"%s\" associated with body.", "*", str.c_str());
     pmSystem.setCmdError(true);
     return;
     }
   }
 else if (pgeom_type == PM_POTENTIAL_GEOM_PLANE) {
   PmVector3 point, normal, scales;

   if (!dlist.getVector("point", point)) {
     pm_ErrorReport (PM, "no plane point given.", "*");
     pmSystem.setCmdError(true);
     return;
     }

   if (!dlist.getVector("normal", normal)) {
     pm_ErrorReport (PM, "no plane normal given.", "*");
     pmSystem.setCmdError(true);
     return;
     }

   params.setPlaneData(point, normal);

   if (dlist.getVector("scales", scales)) {
     params.setPlaneScales(scales);
     }
   }
 else { 
   pm_ErrorReport (PM, "no region specified.", "*");
   pmSystem.setCmdError(true);
   return;
   }

  if (dlist.getBoolean("sidechains", sidechain)) {
    params.setUseSideChains(sidechain);
    }

  PmPotentialGeom *pgeom = PmPotentialGeom::create(name, pgeom_type, body, rgn_name, 
                                                   params);
  body->addPotentialGeometry(pgeom);
  body->getName(bname);
  pm_PrintMsg (CMDPR, "potential geometry \"%s\" added to body \"%s\" ", name.c_str(),
               bname.c_str()); 
  pmSystem.addPotentialGeometry(pgeom);

  //===== set geometry parameters =====//

  if (pgeom_type == PM_POTENTIAL_GEOM_SPHERE) {
    float radius;

    if (dlist.getFloat("radius", radius)) {
      pm_PrintMsg (CMDPR, "sphere radius[%g] \n", radius); 
      PmPotentialSphere *sgeom = dynamic_cast<PmPotentialSphere*>(pgeom);
      sgeom->setSphereRadius(radius);
      }
    }

  else if (pgeom_type == PM_POTENTIAL_GEOM_ELLIPSOID) {
    float scale;

    if (dlist.getFloat("scale", scale)) {
      PmPotentialEllipsoid *egeom = dynamic_cast<PmPotentialEllipsoid*>(pgeom);
      egeom->setEllipsoidScale(scale);
      }
    }

  // lines //

  else if (pgeom_type == PM_POTENTIAL_GEOM_LINES) {
    vector<string> axes_names, fclist;
    vector<int> axes_ids;
    string str;
    int id;
    PmPotentialLines *lgeom = dynamic_cast<PmPotentialLines*>(pgeom);

    if (dlist.getStringList("axes", axes_names)) {
      if (axes_names.size() > 2) {
        pm_ErrorReport (PM, "only two axes can be given.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      for (unsigned int i = 0; i < axes_names.size(); i++) {
        id = atoi(axes_names[i].c_str());

        if ((id < 1) || (id > 3)) {
          pm_ErrorReport (PM, "axis id must be 1, 2 or 3.", "*");
          pmSystem.setCmdError(true);
          return;
          }

        axes_ids.push_back(id);
        }

      lgeom->setAxes(axes_ids);
      }
    }

  else if (pgeom_type == PM_POTENTIAL_GEOM_POINTS) {
    PmPotentialPoints *sgeom = dynamic_cast<PmPotentialPoints*>(pgeom);

    if (dlist.getBoolean("display_spheres", display_spheres)) {
      pgeom->setDisplaySpheres(display_spheres);
      }

    if (dlist.getFloat("radius", radius)) {
      sgeom->setSphereRadius(radius);
      }
    }

  //===== set properties =====//

  if (dlist.getFloat("charge", charge)) {
    pgeom->setCharge(charge);
    }

  //===== set graphics parameters =====//

  if (dlist.getString("display", str)) {
    PmGraphicsGeometry::convDisplayType(str, dtype);
    pgeom->setDisplayType(dtype);
    //fprintf (stderr, ">>>> display \n");
    }

  if (dlist.getVector("color", color)) {
    pgeom->setColor(color);
    }

  if (dlist.getFloat("width", width)) {
    pgeom->setLineWidth(width);
    }

  if (dlist.getBoolean("show", show)) {
    if (show) {
      // note: say what? //
      //body->displayPotentialGeom(name, show);
      pgeom->display(show);
      }
    }

  //pmSystem.addPotential(pot);
  }

//*============================================================*
//*==========             pm_CmdBody                 ==========*
//*============================================================*
// process body command.

void
pm_CmdBody (PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, bname, dname, res, type_str, use_type, body_name;
  string domain_prefix;
  vector<string> body_names;
  vector<string> domain_names;
  vector<string> xform_names;

  PmBody *body;
  PmBodyType btype;
  PmVector3 bpos;
  PmMolecule *domain;
  PmSurface *surface;
  PmParticle *particle;
  PmSolid *solid;
  vector<PmBody*> bodies;

  // get the next data item //

  dlist.getNext (data);

  //===== create a new body and add it to the system =====//

  if (data.name == "create") {
    dlist.getString("name", name);
    dlist.getString("type", type_str);
    body_name = name;

    if (type_str != "") {
      PmBody::convBodyType(type_str, btype);

      if (btype == PM_BODY_UNKNOWN) {
        pm_ErrorReport (PM, "unknown body type \"%s\".", "*", type_str.c_str());
        pmSystem.setCmdError(true);
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no body type specified.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    // create a body of a particular type and add it to the system //

    if (dlist.getStringList("names", body_names)) {
      if (dlist.getStringList("domains", domain_names)) {
        if (domain_names.size() != body_names.size()) {
          pm_ErrorReport (PM, "number of body name != number of domain names", "*");
          pmSystem.setCmdError(true);
          return;
          }

        for (unsigned i = 0; i < body_names.size(); i++) {
          bname = body_names[i];
          body = new PmBody(bname, btype);
          pmSystem.addBody (body);

          dname = domain_names[i];
          pmSystem.getDomain (dname, &domain);

          if (!domain) {
            pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
            pmSystem.setCmdError(true);
            return;
            }

          body->setPhysicalObj(domain);
          }
        }
      }

    //===== domains are given but not names =====//

    else if (dlist.getStringList("domains", domain_names)) {
      stringstream dss;

      for (unsigned i = 0; i < domain_names.size(); i++) {
        dname = domain_names[i];
        pmSystem.getDomain (dname, &domain);

        if (!domain) {
          pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
          pmSystem.setCmdError(true);
          return;
          }

        dss << dname << "Body";
        bname = dss.str();
        dss.str(std::string());

        body = new PmBody(bname, btype);
        pmSystem.addBody (body);
        body->setPhysicalObj(domain);
        //fprintf (stderr, ">>>> add body \"%s\" \n", bname.c_str());
        }
      }

    //===== domain prefix is given =====//

    else if (dlist.getString("domain_prefix", domain_prefix)) {
      vector<PmMolecule*> domains;
      pmSystem.getDomains(domains);
      size_t pos;
      stringstream dss;

      for (unsigned int i = 0; i < domains.size(); i++) {
        domains[i]->getName(dname);
        pos = dname.find(domain_prefix);

        if (pos == 0) {
          dss << dname << "Body";
          bname = dss.str();
          dss.str(std::string());
          pmSystem.getBody (bname, &body);

          if (!body) {
            body = new PmBody(bname, btype);
            pmSystem.addBody (body);
            body->setPhysicalObj(domains[i]);
            pm_PrintMsg (CMDPR, "add body \"%s\" ", bname.c_str());
            }
          }
        }
      }
    else {
      body = new PmBody(name, btype);
      pmSystem.addBody (body);
      body_name = name;
      pm_PrintMsg (CMDPR, "create body \"%s\" ", name.c_str());

      if (dlist.getString("domain", dname)) {
        pmSystem.getDomain (dname, &domain);

        if (!domain) {
          pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
          pmSystem.setCmdError(true);
          return;
          }

        body->setPhysicalObj(domain);
        }

      else if (dlist.getString("particle", dname)) {
        pmSystem.getParticle(dname, &particle);

        if (!particle) {
          pm_ErrorReport (PM, "no particle named \"%s\".", "*", dname.c_str());
          pmSystem.setCmdError(true);
          return;
          }

        body->setPhysicalObj(particle);
        }

      else if (dlist.getString("surface", dname)) {
        pmSystem.getSurface(dname, &surface);

        if (!surface) {
          pm_ErrorReport (PM, "no surface named \"%s\".", "*", dname.c_str());
          pmSystem.setCmdError(true);
          return;
          }

        body->setPhysicalObj(surface);
        }

      else if (dlist.getString("solid", dname)) {
        pmSystem.getSolid(dname, &solid);

        if (!solid) {
          pm_ErrorReport (PM, "no solid named \"%s\".", "*", dname.c_str());
          pmSystem.setCmdError(true);
          return;
          }

        body->setPhysicalObj(solid);
        }
      }

    // add xform objects //

    if (dlist.getStringList("xform_objects", xform_names)) {
      PmPhysicalObj *pobj;
      string pname;

      for (unsigned int i = 0; i < xform_names.size(); i++) {
        pname = xform_names[i];

        if (!PmPhysicalObj::get(pname, &pobj)) {
          pm_ErrorReport (PM, "no physical object named \"%s\".", "*", pname.c_str());
          pmSystem.setCmdError(true);
          return;
          }

        body->addXformPhysicalObj(pobj);
        }
      }
    }

  else {
    pmSystem.getBody (data.name, &body);
    body_name = data.name;
    }

  if (!body) {
    pm_ErrorReport (PM, "no body named \"%s\".", "*", data.name.c_str());
    pmSystem.setCmdError(true);
    return;
    }

  // process body parameters //

  bool show_set = true;

  while (dlist.getNext(data)) {
    if (data.name == "add") {
      dlist.getNext(data, dv);

      if (data.name == "potential_geometry") {
        pm_CmdBodyAddPotentialGeometry (body, dv, dlist);
        }
      }

    // ===== contact ===== //

    else if (data.name == "contact") {
      pm_CmdBodyContact(body, dlist);
      return;
      }

    // ===== damping ===== //

    else if (data.name == "damping") {
      float val = data.getFloat();
      body->setDampingFactor(val);
      }

    // ===== print ===== //

    else if (data.name == "print") {
      string str;

      if (dlist.getString("properties", str)) {
        if (str == "mass") {
          PmMassProperties mass_props;
          body->getMassProps(mass_props);
          fprintf (stderr, "\n------ body \"%s\" mass properties ------ \n", 
                   body_name.c_str());
          fprintf (stderr, "mass = %f \n", mass_props.mass);
          fprintf (stderr, "center of mass = (%f, %f, %f) \n", mass_props.com[0],
                   mass_props.com[1], mass_props.com[2]);
          fprintf (stderr, "inertia tensor = [ %f  %f  %f  \n", mass_props.inertia(0,0), 
                                                                mass_props.inertia(0,1),
                                                                mass_props.inertia(0,2));
          fprintf (stderr, "                   %f  %f  %f  \n", mass_props.inertia(1,0), 
                                                                mass_props.inertia(1,1),
                                                                mass_props.inertia(1,2));
          fprintf (stderr, "                   %f  %f  %f ]\n", mass_props.inertia(2,0), 
                                                                mass_props.inertia(2,1),
                                                                mass_props.inertia(2,2));
          }
        }
      }

    // add xform objects //

    else if (data.name == "xform_objects") {
      PmPhysicalObj *pobj;
      string pname;
      data.getStringList(xform_names);

      for (unsigned int i = 0; i < xform_names.size(); i++) {
        pname = xform_names[i];

        if (!PmPhysicalObj::get(pname, &pobj)) {
          pm_ErrorReport (PM, "no physical object named \"%s\".", "*", pname.c_str());
          pmSystem.setCmdError(true);
          return;
          }

        body->addXformPhysicalObj(pobj);
        }
      }
    }

  if (show_set) {
    //body->display (show);
    }
  }

//*============================================================*
//*==========             pm_CmdBodies               ==========*
//*============================================================*
// process bodies command.

void
pm_CmdBodies (PmCmdDataList& dlist)
  {
  PmCmdData data;
  PmVector3 color;
  bool flag, show = false;
  float val;
  vector<PmBody*> bodies;
  PmBody* body;
  vector<string> bnames;
  string prefix;

  // get body names //

  // use a list of names //

  if (dlist.getStringList("list", bnames)) {
    for (unsigned int i = 0; i < bnames.size(); i++) {
      pmSystem.getBody(bnames[i], &body);

      if (!body) {
        pm_ErrorReport (PM, "no body named \"%s\".", "*", bnames[i].c_str());
        pmSystem.setCmdError(true);
        return;
        }

      bodies.push_back(body);
      }
    }

  // use a prefix //

  else if (dlist.getString("prefix", prefix)) {
    vector<PmBody*> search_bodies;
    string bname;
    pmSystem.getBodies(search_bodies);

    for (unsigned int i = 0; i < search_bodies.size(); i++) {
      search_bodies[i]->getName(bname);

      if (bname.find(prefix) == 0) {
        bodies.push_back(search_bodies[i]);
        }
      }

    if (bodies.empty()) { 
      pm_ErrorWarnReport (PM, "no bodies found with prefix = \"%s\".", "*", 
                          prefix.c_str());
      }
    }

  // use all the bodies defined in the system //

  else {
    pmSystem.getBodies(bodies);
    }


  // process parameters //

  while (dlist.getNext(data)) {

    // ===== contact ===== //

    if (data.name == "contact") {
      color.set(1,1,1);

      if (dlist.getBoolean("show", flag)) {
        show = flag;
        }

      dlist.getVector("color", color);
      PmBody::checkContact(color, show);
      return;
      }

    // ===== damping ===== //

    else if (data.name == "damping") {
      float val = data.getFloat();

      for (int i = 0; i < bodies.size(); i++) {
        bodies[i]->setDampingFactor(val);
        }
      }
    }
  }

