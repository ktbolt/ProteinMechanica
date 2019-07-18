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
//     r i g i d   s i m u l a t i o n   c o m m a n d s     //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              pm_CmdSimulations         ==========*
//*============================================================*
// process simualtions command.

void
pm_CmdSimulations (PmCmdDataList& dlist)
  {
  string dv, name, dname, res, type_str, sim_name;
  vector<PmSimulation*> simulations;
  PmSimulation *sim;
  PmSimulationType stype;
  PmCmdData data;
  dlist.getNext (data, dv);

  pmSystem.getSimulations (simulations);
  }

//*============================================================*
//*==========          pm_CmdRigidSimAddTrace        ==========*
//*============================================================*
// process rigid simualtion trace command.

void
pm_CmdRigidSimAddTrace(string tname, PmRigidSimulation *sim, PmCmdDataList& dlist)
  {
  string dv, bname, name;
  PmCmdData data;
  PmVector3 point, dir;
  PmCmdData pdata;
  PmBody *body;
  PmTrace *trace;
  PmVector3 color;
  bool global_frame;
  float width;

  if (!dlist.getString("name", name)) { 
    pm_ErrorReport (PM, "no trace name specified.", "*");
    return;
    }

  if (!dlist.getData("point", pdata)) {
    pm_ErrorReport (PM, "trace point not specified.", "*");
    return;
    }

  if (!pm_CmdProcPosition(pdata, point)) {
    pm_ErrorReport (PM, "couldn't process trace point.", "*");
    return;
    }

  if (!dlist.getString("body", bname)) {
    pm_ErrorReport (PM, "trace body not specified.", "*");
    return;
    }

  if (!sim->getBody(bname, &body)) {
    pm_ErrorReport (PM, "no body \"%s\" in simulation.", "*", bname.c_str());
    return;
    }

  trace = new PmTrace(name, body, point);

  if (dlist.getBoolean("global_frame", global_frame)) {
    trace->setGlobalFrame(global_frame);
    }

  if (dlist.getVector("color", color)) {
    trace->setColor(color);
    }

  if (dlist.getFloat("width", width)) {
    trace->setLineWidth(width);
    }

  if (dlist.getString("display", dv)) {
    PmGeometryDisplayType dtype;
    PmGraphicsGeometry::convDisplayType (dv, dtype);
    trace->setDisplayType(dtype);
    }

  // add trace to simulation //

  if (sim->addTrace(trace)) {
    pm_PrintMsg (CMDPR, "trace \"%s\" added to simulation \"%s\" ", name.c_str(),
                 name.c_str());
    }
  }

//*============================================================*
//*==========     pm_CmdRigidSimAddBodyInteraction   ==========*
//*============================================================*
// process rigid simulation add interaction command with body list.

void
pm_CmdRigidSimAddBodyInteraction (string name, PmSimulation *sim, PmPotentialType ptype,
                                  vector<string> body_names, string association,
                                  PmCmdDataList& dlist)
  {
  fprintf (stderr, "\n #######  pm_CmdRigidSimAddBodyInteraction  \n");
  string dv, bname1, bname2, rgn_name, pgeom_name, pot_name, int_name; 
  stringstream dss;
  PmCmdData data, tdata;
  vector<string> pot_names;
  vector<PmBody*> bodies1, bodies2;
  PmBody *body, *body1, *body2;
  vector<PmPotential*> potentials1, potentials2;
  vector<PmPotentialGeom*> glist; 
  PmPhysicalObj *pobj1, *pobj2;
  PmPotentialGeom *pgeom1, *pgeom2; 
  vector<PmRegion*> rgns;
  PmTimeInterval time;
  string bname, pname, tstr, name1, name2;

  //===== check body names =====//

  for (unsigned int i = 0; i < body_names.size(); i++) {
    bname = body_names[i];
    pmSystem.getBody(bname, &body);
    fprintf (stderr, ">>> add body %s \n", bname.c_str());

    if (!body) {
      pm_ErrorReport (PM, "no body named \"%s\".", "*", bname.c_str());
      return;
      }
    }

  if (dlist.getData("time", tdata)) {
    pm_CmdProcTimeInterval(tdata, time);
    }

  //===== create spring potential between bodies =====//

  if (ptype == PM_POTENTIAL_SPRING) {
    float cutoff;
    bool hbonds, use_sidechains, use_ca, show_geom, no_mainchain;
    string chain; 
    PmPotentialGeomType pgeom_type = PM_POTENTIAL_GEOM_POINTS;
    PmPotentialParameters params;
    PmInteraction *inter;
    vector<PmPotential*> plist; 
    PmPotential *pot; 
    PmAtomFilter filter;
    PmMolRegionParameters rgn_params;
    PmObjectType ptype1, ptype2;
    vector<string> map_scale;
    float map_min, map_max, max_strain;
    vector<string> atom_names;
    int num_bodies = body_names.size();

    // get parameters for spring potential //

    if (!dlist.getFloat("cutoff", cutoff)) { 
      cutoff = 1.0;
      }

    if (!dlist.getFloat("max_strain", max_strain)) {  
      max_strain = 0.0;
      }

    if (!dlist.getBoolean("hbonds", hbonds)) {
      hbonds = false;
      }

    dlist.getBoolean("show_geom", show_geom);

    if (dlist.getBoolean("use_sidechains", use_sidechains)) {
      filter.exclude = true;
      filter.sidechain = true;
      bool bval;
      dlist.getBoolean("ljspring", bval);

      if (bval) {
        pm_ErrorReport (PM, "can't use LJ springs with use_sidechain option.", "*");
        return;
        }
      }
    else if (dlist.getStringList("atom_names", atom_names)) {
      if (atom_names.size() == 1) {
        if (atom_names[0] == "mc") {
          filter.mainchain = true;
          }
        else if (atom_names[0] == "sc") {
          filter.sidechain = true;
          }
        else if (atom_names[0] != "all") {
          filter.names = atom_names;
          }
        }
      else if (atom_names.size() != 0) {
        filter.names = atom_names;
        }
      }

    if (dlist.getBoolean("use_ca", use_ca)) {
      if (use_ca) filter.names.push_back("CA");
      }

    if (hbonds) { 
      if (!dlist.getString("chain", chain)) {
        pm_ErrorReport (PM, "no chain id given.", "*");
        return;
        }

      params.setHydrogenBonds(true);
      }

    // create regions and potential geometries //

    pm_PrintMsg (CMDPR, ">>> create regions and potential geometries");
    pm_PrintMsg (CMDPR, "num bodies=%d", num_bodies);
    pm_PrintMsg (CMDPR, "cutoff=%f", cutoff);
    pm_PrintMsg (CMDPR, "association=%s", association.c_str());
    vector<string> blist1, blist2;

    if (association == "full") {
      for (int i = 0; i < num_bodies-1; i++) {
        bname1 = body_names[i]; 

        for (int j = i+1; j < num_bodies; j++) {
          bname2 = body_names[j]; 
          blist1.push_back(bname1);
          blist2.push_back(bname2);
          }
        }
      }
    else { 
      for (int i = 0; i < num_bodies-1; i++) {
        bname1 = body_names[i]; 
        bname2 = body_names[i+1]; 
        blist1.push_back(bname1);
        blist2.push_back(bname2);
        }
      }

    pm_PrintMsg (CMDPR, "num assoc=%d", blist1.size()); 

    for (int i = 0; i < blist1.size(); i++) {
      bname1 = blist1[i]; 
      pmSystem.getBody(bname1, &body1);
      bname2 = blist2[i]; 
      pmSystem.getBody(bname2, &body2);
      body1->getPhysicalObj(&pobj1);
      body2->getPhysicalObj(&pobj2);
      ptype1 = pobj1->getObjectType();
      ptype2 = pobj2->getObjectType();

      // create regions for body1, body2 //

      if (hbonds) {
        pobj1->createHbondRegions(name, pobj2, cutoff, chain, rgns);
        }
      else {
        pobj1->createRegions(name, pobj2, cutoff, use_sidechains, filter, rgns);
        }

      if (rgns.size() == 0) {
        //pm_ErrorWarnReport (PM, "can't create regions for potential.", "*");
        continue;
        }

      pm_PrintMsg (CMDPR, ">>> body1=%s  body2=%s", bname1.c_str(), bname2.c_str());

      // create potential geometries //

      for (int k = 0; k < rgns.size(); k++) {
        rgn_name = rgns[k]->name;

        if (k == 0) {
          body = body1;
          }
        else {
          body = body2;
          }

        body->getName(name1);
        dss << rgn_name << name1;
        pgeom_name = dss.str();
        dss.str(std::string());

        PmPotentialGeom *pgeom = PmPotentialGeom::create(pgeom_name, pgeom_type, body, 
                                                         rgn_name, params);
        body->addPotentialGeometry(pgeom);
        pmSystem.addPotentialGeometry(pgeom);
        pgeom->display(show_geom);
        glist.push_back(pgeom);
        pm_PrintMsg (CMDPR, ">>> potential geometry=%s", pgeom_name.c_str());
        }
      }

    // create potentials //

    pm_PrintMsg (CMDPR, "------ create potentials ------");

    if (glist.size() == 0) {
      pm_ErrorWarnReport (PM, "can't create potentials.", "*");
      return;
      }

    for (int i = 0; i < glist.size()-1; i += 1) {
      pgeom1 = glist[i];
      pgeom2 = glist[i+1];
      pgeom1->getName(name1);
      pgeom2->getName(name2);
      dss << name1 << '-' << name2 << "Potential";
      pot_name = dss.str();
      dss.str(std::string());

      pot = PmPotential::create(pot_name, ptype, pgeom1, pgeom2);
      pmSystem.addPotential(pot);
      PmSpringPotential *spot = dynamic_cast<PmSpringPotential*>(pot);
      spot->setCutoff(cutoff);
      spot->setMaxStrain(max_strain);
      plist.push_back(pot);
      pm_PrintMsg (CMDPR, ">>> potential=%s", pot_name.c_str());
      //pot->displayGeometry(true);
      }

    // create interaction //

    inter = new PmInteraction(name, plist);
    inter->setTimeInterval(time);
    sim->addInteraction(inter);

    // set parameters //

    pm_CmdPotentialParams (plist, ptype, dlist);
    }

  else {
    pm_ErrorReport (PM, "interaction type \"%s\" not supported for body list.", "*", 
                    tstr.c_str());
    }
  }

//*============================================================*
//*==========       pm_CmdRigidSimAddInteraction     ==========*
//*============================================================*
// process rigid simulation add interaction command.

void
pm_CmdRigidSimAddInteraction (string name, PmSimulation *sim, PmCmdDataList& dlist)
  {
  string dv;
  PmCmdData data, tdata;
  PmBody *body;
  PmPotential *pot;
  PmPotentialGeom *pgeom, *pgeom1, *pgeom2;
  vector<PmPotentialGeom*> glist1, glist2; 
  PmInteraction *inter;
  string bname, pname, str, association;
  vector<string> body_names, pot_names, list1, list2, list;
  vector<PmPotential*> potentials;
  PmPotentialType ptype;
  PmTimeInterval time; 
  PmVector3 color;
  bool verbose, pairwise, report_energy;
  fprintf (stderr, "\n #######  pm_CmdRigidSimAddInteraction \n");

  //===== get time interval =====// 

  if (dlist.getData("time", tdata)) {
    pm_CmdProcTimeInterval(tdata, time);
    }

  if (!dlist.getBoolean("report_energy", report_energy)) {
    report_energy = true;
    }

  //====== check for list of potential names =====//

  if (dlist.getStringList("potentials", pot_names)) { 
    for (unsigned int i = 0; i < pot_names.size(); i++) {
      pname = pot_names[i];
      pmSystem.getPotential(pname, &pot);

      if (!pot) {
        pm_ErrorReport (PM, "no potential named \"%s\".", "*", pname.c_str());
        return;
        }

      potentials.push_back(pot);
      }

    inter = new PmInteraction(name, potentials); 
    inter->setTimeInterval(time);
    sim->addInteraction(inter);
    pm_PrintMsg (CMDPR, "interaction \"%s\" added \n", name.c_str());
    return;
    }

  //===== get potential type =====//

  if (dlist.getString("type", str)) { 
    if (str != "") {
      PmPotential::convType(str, ptype);

      if (ptype == PM_POTENTIAL_UNKNOWN) {
        pm_ErrorReport (PM, "unknown potential type \"%s\".", "*", str.c_str());
        return;
        }
      }
    }
  else {
    pm_ErrorReport (PM, "no potential type specified.", "*");
    return;
    }

  //===== get association =====//

  if (!dlist.getString("association", association)) {
    association = "pairwise";
    }

  if (association == "pairwise") { 
    pairwise = true;
    }
  else if (association == "full") { 
    pairwise = false;
    }
  else {
    pm_ErrorReport (PM, "unknown association type \"%s\".", "*", association.c_str());
    return;
    }

  //====== check for two lists of potential geometries =====//

  if (dlist.getStringList("list1", list1)) { 
    vector<PmPotentialGeom*> pgeom1, pgeom2;
    string pname;

    if (!dlist.getStringList("list2", list2)) { 
      pm_ErrorReport (PM, "no list2 specified.", "*");
      return;
      }

    if (pairwise && (list1.size() != list2.size())) {
      pm_ErrorReport (PM, "list1 and list2 are not the same size. ", "*");
      return;
      }

    // check potential geometry names //

    vector<PmPotentialGeom*> tlist1, tlist2; 

    for (unsigned int i = 0; i < list1.size(); i++) {
      pname = list1[i];
      pmSystem.getPotentialGeometry(pname, &pgeom);

      if (!pgeom) {
        pm_ErrorReport (PM, "no potential geometry named \"%s\".", "*", pname.c_str());
        return;
        }

      tlist1.push_back(pgeom);
      }

    for (unsigned int i = 0; i < list2.size(); i++) {
      pname = list2[i];
      pmSystem.getPotentialGeometry(pname, &pgeom);

      if (!pgeom) {
        pm_ErrorReport (PM, "no potential geometry named \"%s\".", "*", pname.c_str());
        return;
        }

      tlist2.push_back(pgeom);
      }

    if (pairwise) {
      for (unsigned int i = 0; i < list1.size(); i++) {
        glist1.push_back(tlist1[i]);
        glist2.push_back(tlist2[i]);
        }
      }
    else {
      for (unsigned int i = 0; i < list1.size(); i++) {
        pgeom = tlist1[i];

        for (unsigned int j = 0; j < list2.size(); j++) {
          glist1.push_back(pgeom);
          glist2.push_back(tlist2[j]);
          }
        }
      }
    }

  //====== check for single list of potential geometries =====//

  else if (dlist.getStringList("list", list)) { 
    vector<PmPotentialGeom*> pgeoms;
    string pname;

    // check potential geometry names //

    for (unsigned int i = 0; i < list.size(); i++) {
      pname = list[i];
      pmSystem.getPotentialGeometry(pname, &pgeom);

      if (!pgeom) {
        pm_ErrorReport (PM, "no potential geometry named \"%s\".", "*", pname.c_str());
        return;
        }

      pgeoms.push_back(pgeom);
      }

    if (pairwise) { 
      for (unsigned int i = 0; i < pgeoms.size()-1; i++) {
        glist1.push_back(pgeoms[i]);
        glist2.push_back(pgeoms[i+1]);
        }
      }
    else {
      for (unsigned int i = 0; i < pgeoms.size()-1; i++) {
        for (unsigned int j = i+1; j < pgeoms.size(); j++) {
          glist1.push_back(pgeoms[i]);
          glist2.push_back(pgeoms[j]);
          }
        }
      }
    }

  //====== check for list of body names =====//

  else if (dlist.getStringList("bodies", body_names)) { 
    fprintf (stderr, "\n #######  bodies \n");
    pm_CmdRigidSimAddBodyInteraction (name, sim, ptype, body_names, association, dlist);
    return;
    }

  //===== create potentials =====// 

  stringstream dss;
  string gname1, gname2, pot_name;

  for (unsigned int i = 0; i < glist1.size(); i++) {
    pgeom1 = glist1[i];
    pgeom1->getName(gname1);
    pgeom2 = glist2[i];
    pgeom2->getName(gname2);
    dss << name << '_' << gname1 << '-' << gname2;
    pot_name = dss.str();
    dss.str(std::string());
    pot = PmPotential::create(pot_name, ptype, pgeom1, pgeom2);
    pmSystem.addPotential(pot);
    potentials.push_back(pot);
    pm_PrintMsg (CMDPR, "potential \"%s\" added", pot_name.c_str());
    pgeom1->addPotential(pot);
    pgeom2->addPotential(pot);
    pot->setReportEnergy(report_energy);
    //fprintf (stderr, "\n #######  add pot to pgeom1  %s \n", gname1.c_str());
    //fprintf (stderr, "          add pot to pgeom2  %s \n", gname2.c_str());
    }

  //===== set potential parameters =====// 

  pm_CmdPotentialParams (potentials, ptype, dlist);

  //===== create interaction =====// 

  inter = new PmInteraction(name, potentials);
  inter->setTimeInterval(time);
  sim->addInteraction(inter);

  if (dlist.getBoolean("verbose", verbose)) {
    inter->setVerbose(verbose);
    }
  }

//*============================================================*
//*==========       pm_CmdRigidSimAddRestraint       ==========*
//*============================================================*
// process rigid simualtion add restraint command.

void
pm_CmdRigidSimAddRestraint (string name, PmSimulation *sim, PmCmdDataList& dlist)
  {
  PmCmdData data, tdata;
  PmBody *body1 = NULL, *body2 = NULL;
  string str, bname, rgn_name, region1, region2;
  PmRestraint *res;
  PmTimeInterval tint; 
  float force_const, distance, power_params[2], ramp;
  float width;
  PmVector3 color;
  bool show, abs_dist, use_pca, comp_energy;
  vector<string> strl; 
  PmRestraintType type;

  if (dlist.getString("body1", bname)) {
    pmSystem.getBody(bname, &body1);

    if (!body1) {
      pm_ErrorReport (PM, "no body named \"%s\".", "*", bname.c_str());
      return;
      }

    if (dlist.getString("region1", rgn_name)) {
      if (body1->hasRegion(rgn_name)) {
        region1 = rgn_name;
        }
      else {
        pm_ErrorReport (PM, "no region \"%s\" associated with body1.", "*",
                       rgn_name.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no region1 specified.", "*");
      return;
      }
    }
  else {
    pm_ErrorReport (PM, "no body1.", "*");
    return;
    }

  if (dlist.getString("body2", bname)) {
    pmSystem.getBody(bname, &body2);

    if (!body2) {
      pm_ErrorReport (PM, "no body named \"%s\".", "*", bname.c_str());
      return;
      }

    if (dlist.getString("region2", rgn_name)) {
      if (body2->hasRegion(rgn_name)) {
        region2 = rgn_name;
        }
      else {
        pm_ErrorReport (PM, "no region \"%s\" associated with body2.", "*",
                       rgn_name.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no region2 specified.", "*");
      return;
      }
    }
  else {
    pm_ErrorReport (PM, "no body2.", "*");
    return;
    }

  if (!dlist.getFloat("force_const", force_const)) {
    force_const = 1.0;
    }

  if (!dlist.getFloat("distance", distance)) {
    distance = 1.0;
    }

  if (dlist.getData("time", tdata)) {
    pm_CmdProcTimeInterval(tdata, tint);
    }

  /*
  if (dlist.getBoolean("use_pca", use_pca) && use_pca) {
    fprintf (stderr, "\n ###  rgn1 cords = %d  \n", region1->coords.size()); 
    }
  */

  // add restraint to simulation //

  res = new PmRestraint(name);
  res->setSimObjs(body1, body2);
  res->setRegions(region1, region2);
  res->setTime(tint);
  res->setSpringData(force_const, distance);

  if (dlist.getStringList("power", strl)) {
    convToFloat (strl[0], power_params[0]);
    convToFloat (strl[1], power_params[1]);
    res->setPowerParams(power_params);
    }

  if (dlist.getString("type", str)) {
    PmRestraint::convType(str, type);

    if (type == PM_RESTRAINT_UNKNOWN) {
      pm_ErrorReport (PM, "unknown restraint type \"%s\".", "*", str.c_str());
      return;
      }

    res->setType(type);
    }

  if (dlist.getFloat("ramp", ramp)) {
    res->setRampParams(ramp);
    }

  if (dlist.getBoolean("absolute_distance", abs_dist)) {
    res->setAbsoluteDistance(abs_dist);
    }

  if (dlist.getBoolean("compute_energy", comp_energy)) {
    res->setCompEnergy(comp_energy);
    }

  if (dlist.getVector("color", color)) { 
    res->setColor(color);
    }

  if (dlist.getFloat("width", width)) {
    res->setWidth(width);
    }

  if (dlist.getBoolean("show", show)) {
    res->setShow(show);
    }

  sim->addRestraint(res);
  pm_PrintMsg (CMDPR, "restraint \"%s\" added \n", name.c_str());
  }

//*============================================================*
//*==========              pm_CmdRigidSim            ==========*
//*============================================================*
// process rigid simualtion command.

void
pm_CmdRigidSim (PmSimulation *sim, PmCmdDataList& dlist)
  {
  string dv, name, bname, jname;
  vector<string> body_names; 
  vector<string> joint_names; 
  PmCmdData data, tdata;
  PmBody *body;
  PmJoint *joint;
  PmModel *model;
  PmTimeInterval time; 
  string prefix;

  sim->getName(name);
  PmRigidSimulation *rsim = PmSimulation::getRigid(name);

  // process parameters //

  while (dlist.getNext(data, dv)) {

    //====== add ======//

    if (data.name == "add") {
      dlist.getNext(data, dv);

      if (data.name == "body") {
        pmSystem.getBody(dv, &body);

        if (!body) {
          pm_ErrorReport (PM, "no body named \"%s\".", "*", dv.c_str());
          return;
          }

        if (rsim->addBody(body)) {
          pm_PrintMsg (CMDPR, "body \"%s\" added to simulation \"%s\" ", 
                       dv.c_str(), name.c_str());
          }
        else {
          pm_ErrorReport (PM, "body \"%s\" not added", "*", dv.c_str());
          }
        }

      // add bodies using all or a list //

      else if (data.name == "bodies") {

        if (dv == "all") {
          vector<PmBody*> bodies;
          pmSystem.getBodies(bodies);
          pm_PrintMsg (CMDPR, "add all %d bodies. \n", bodies.size());

          for (unsigned int i = 0; i < bodies.size(); i++) {
            rsim->addBody(bodies[i]);
            }
          }

        else if (dv == "prefix") {
          size_t pos;
          stringstream dss;
          vector<PmBody*> bodies;

          if (!dlist.getString("prefix", prefix)) {
            pm_ErrorReport (PM, "no prefix given.", "*");
            return;
            }

          pmSystem.getBodies(bodies);

          for (unsigned int i = 0; i < bodies.size(); i++) {
            bodies[i]->getName(bname);
            pos = bname.find(prefix);

            if (pos == 0) {
              rsim->addBody(bodies[i]);
              }
            }
          }

        else if (dlist.getStringList("bodies", body_names)) {
          for (unsigned int i = 0; i < body_names.size(); i++) {
            bname = body_names[i];
            pmSystem.getBody(bname, &body);

            if (!body) {
              pm_ErrorReport (PM, "no body named \"%s\".", "*", bname.c_str());
              return;
              }

            if (rsim->addBody(body)) {
              pm_PrintMsg (CMDPR, "body \"%s\" added to simulation \"%s\" ", 
                           bname.c_str(), name.c_str());
              }
            }
          }
        }


      // add a force //

      else if (data.name == "force") {
        string fname, bname;
        PmForce *force;

        if (dv == "") {
          pm_ErrorReport (PM, "force name not specified.", "*");
          return;
          }

        pmSystem.getForce(dv, &force);

        if (!force) {
          pm_ErrorReport (PM, "no force named \"%s\".", "*", dv.c_str());
          return;
          }

        if (!dlist.getString("body", bname)) {
          pm_ErrorReport (PM, "force body not specified.", "*");
          return;
          }

        if (!rsim->getBody(bname, &body)) {
          pm_ErrorReport (PM, "no body \"%s\" in simulation.", "*", bname.c_str());
          return;
          }

        if (dlist.getData("time", tdata)) {
          pm_CmdProcTimeInterval(tdata, time);
          }

        if (rsim->addForce(bname, force, time)) {
          pm_PrintMsg (CMDPR, "force \"%s\" added to simulation \"%s\" ", dv.c_str(), 
                       name.c_str());
          }
        }

      // add a geometry //

      else if (data.name == "geometry") {
        string gname, bname;
        PmGeometry *geom;

        if (dv == "") {
          pm_ErrorReport (PM, "geometry name not specified.", "*");
          return;
          }

        pmSystem.getGeometry(dv, &geom);

        if (!geom) {
          pm_ErrorReport (PM, "no geometry named \"%s\".", "*", dv.c_str());
          return;
          }

        if (!dlist.getString("body", bname)) {
          pm_ErrorReport (PM, "geometry body not specified.", "*");
          return;
          }

        if (!rsim->getBody(bname, &body)) {
          pm_ErrorReport (PM, "no body \"%s\" in simulation.", "*", bname.c_str());
          return;
          }

        geom->setBody(body);
        rsim->addGeometry(geom);
        }

      // add an interaction //

      else if (data.name == "interaction") {
        pm_CmdRigidSimAddInteraction(dv, sim, dlist);
        }

      // add a joint //

      else if (data.name == "joint") {
        pmSystem.getJoint(dv, &joint);

        if (!joint) {
          pm_ErrorReport (PM, "no joint named \"%s\".", "*", dv.c_str());
          return;
          }

        if (rsim->addJoint(joint)) {
          pm_PrintMsg (CMDPR, "joint \"%s\" added to simulation \"%s\" ", dv.c_str(), 
                       name.c_str());
          }
        else {
          pm_ErrorReport (PM, "joint \"%s\" not added", "*", dv.c_str());
          }
        }

      // add joints //

      else if (data.name == "joints") {
        if (dv == "all") {
          vector<PmJoint*> joints;
          pmSystem.getJoints(joints);
          pm_PrintMsg (CMDPR, "add all %d joints. ", joints.size());

          for (unsigned int i = 0; i < joints.size(); i++) {
            rsim->addJoint(joints[i]);
            }
          }

        else if (dv == "prefix") {
          size_t pos;
          stringstream dss;
          vector<PmJoint*> joints;

          if (!dlist.getString("prefix", prefix)) {
            pm_ErrorReport (PM, "no prefix given.", "*");
            return;
            }

          pmSystem.getJoints(joints);

          for (unsigned int i = 0; i < joints.size(); i++) {
            joints[i]->getName(jname);
            pos = jname.find(prefix);

            if (pos == 0) {
              rsim->addJoint(joints[i]);
              }
            }
          }

        else if (!dlist.getStringList("joints", joint_names)) {
          pm_ErrorReport (PM, "wrong format for joint names.", "*");
          return;
          }

        for (unsigned int i = 0; i < joint_names.size(); i++) {
          jname = joint_names[i];
          pmSystem.getJoint(jname, &joint);

          if (!joint) {
            pm_ErrorReport (PM, "no joint named \"%s\".", "*", jname.c_str());
            return;
            }

          if (rsim->addJoint(joint)) {
            pm_PrintMsg (CMDPR, "joint \"%s\" added to simulation \"%s\" ", 
                         jname.c_str(), name.c_str());
            }
          }
        }

      // add a measurement //

      else if (data.name == "measurement") {
        PmMeasurement *msr; 

        if (dv == "") {
          pm_ErrorReport (PM, "measurement name not specified.", "*");
          return;
          }

        pmSystem.getMeasurement(dv, &msr);

        if (!msr) {
          pm_ErrorReport (PM, "no measurement named \"%s\".", "*", dv.c_str());
          return;
          }

        if (rsim->addMeasurement(msr)) {
          pm_PrintMsg (CMDPR, "measurement \"%s\" added to simulation \"%s\" ",
                       dv.c_str(), name.c_str());
          }
        }

      // add a motor //

      else if (data.name == "motor") {
        PmMotor *motor;

        if (dv == "") {
          pm_ErrorReport (PM, "motor name not specified.", "*");
          return;
          }

        pmSystem.getMotor(dv, &motor);

        if (!motor) {
          pm_ErrorReport (PM, "no motor named \"%s\".", "*", dv.c_str());
          return;
          }

        if (dlist.getData("time", tdata)) {
          pm_CmdProcTimeInterval(tdata, time);
          }

        if (rsim->addMotor(motor, time)) {
          pm_PrintMsg (CMDPR, "motor \"%s\" added to simulation \"%s\" ",
                       dv.c_str(), name.c_str());
          }
        }

      // add a restraint // 

      else if (data.name == "restraint") {
        pm_CmdRigidSimAddRestraint (dv, sim, dlist);
        }

      else if (data.name == "model") {
        pmSystem.getModel(dv, &model);
        }

      else if (data.name == "trace") {
        pm_CmdRigidSimAddTrace(dv, rsim, dlist);
        }
      }

    //====== contact ======//

    else if (data.name == "contact") {
      vector<PmInteraction*> intr;
      string name, pstype;
      PmPotentialType ptype;
      vector<PmPotential*> potentials;
      bool contact;
      int num_contact = 0;

      if (dlist.getString("check", dv)) {
        sim->getInteractions(intr);
        fprintf (stderr, "\n---------- check contact ---------- \n");

        for (unsigned int i = 0; i < intr.size(); i++) {
          intr[i]->getName(name);
          intr[i]->getPotentialType(ptype);
          intr[i]->getPotentials(potentials);

          if (ptype == PM_POTENTIAL_CONTACT) {
            PmContactPotential *cpot = dynamic_cast<PmContactPotential*>(potentials[0]);
            cpot->checkContact(contact);

            if (contact) {
              fprintf (stderr, ">>> \"%s\" in contact. \n", name.c_str());
              num_contact += 1;
              }
            }
          }

        if (num_contact == 0) {
          fprintf (stderr, ">>> no contact. \n");
          }
        }
      }

    //====== set damping on/off ======//

    else if (data.name == "damping") {
      bool flag = (dv == "on"); 
      rsim->setDamping(flag);
      }

    //====== initialize ======//

    else if (data.name == "initialize") {
      rsim->initialize();
      }

    //====== interactions ======//

    else if (data.name == "interactions") {
      vector<PmInteraction*> intr;
      string name, pstype;
      PmPotentialType ptype;
      vector<PmPotential*> potentials;
      int num_int;

      if (dlist.getString("print", dv)) {
        sim->getInteractions(intr);
        fprintf (stderr, "\n---------- interactions ---------- \n");
        for (unsigned int i = 0; i < intr.size(); i++) { 
          intr[i]->getName(name);
          intr[i]->getPotentialType(ptype);
          intr[i]->getPotentials(potentials);
          intr[i]->getNumInteractions(num_int);
          PmPotential::convType(ptype, pstype);
          fprintf (stderr, "%2d:  name=\"%s\"  \n", i+1, name.c_str());
          fprintf (stderr, "     type=\"%s\"  \n", pstype.c_str());

          if (ptype == PM_POTENTIAL_SPRING) {
            float force_const;
            bool ljspring;
            PmSpringPotential *spot = dynamic_cast<PmSpringPotential*>(potentials[0]);
            spot->getForceConst(force_const);
            spot->getLJSpring(ljspring);
            fprintf (stderr, "     number of springs=%d \n", num_int); 
            fprintf (stderr, "     force constant=%f \n", force_const); 
            if (ljspring) {
              fprintf (stderr, "     ljspring=true \n"); 
              }
            }

          else if (ptype == PM_POTENTIAL_CONTACT) {
            }

          fprintf (stderr, "\n\n");
          }
        }
      }

    //====== set momentum on/off ======//

    else if (data.name == "momentum") {
      bool flag = (dv == "off"); 
      rsim->setMomentumOff(flag);
      }

    //====== print ======//

    else if (data.name == "print") {
      string itype;
      bool pbodies, pjoints;
      dlist.getBoolean("bodies", pbodies); 
      dlist.getBoolean("joints", pjoints); 
      rsim->printModelInformation(pbodies, pjoints);
      }

    //====== read results ======//

    else if (data.name == "read") {
      string file_name;
      bool bval;

      if (!dlist.getString("file_name", file_name)) {
        pm_ErrorReport (PM, "file name must be given. ", "*");
        return;
        }

      // read state data //

      if (dlist.getBoolean("state", bval)) {
        if (bval) {
          rsim->readState(file_name);
          }
        }
      }

    //====== restraints ======//

    else if (data.name == "restraints") {
      vector<PmRestraint*> rests;
      string name, pstype;
      float dist, pdist;

      if (dlist.getString("check", dv)) {
        sim->getRestraints(rests);
        fprintf (stderr, "\n---------- check restraints ---------- \n");

        for (unsigned int i = 0; i < rests.size(); i++) {
          rests[i]->getName(name);
          rests[i]->getPointDistance(dist, pdist);
          fprintf (stderr, ">>> %s  distance=%f \n", name.c_str(), dist);
          }
        }
      }

    //====== step ======//

    else if (data.name == "step") {
      int nsteps;
      string bexp_name;
      bool step_silent;

      if (dv != "") {
        nsteps = data.getInt();
        rsim->setNumInc(nsteps);
        }
      else if (dlist.getString("while", bexp_name)) {
        PmBooleanExpression *bexp;
        //fprintf (stderr, "\n>>>>>>> pm_CmdSimulation:  while=%s \n", bexp_name.c_str());
        if (!pm_CmdGetBooleanExpression(bexp_name, &bexp)) {
          pm_ErrorReport (PM, "unknown boolean expression \"%s\".", "*", 
                          bexp_name.c_str());
          return;
          }

        rsim->setWhileExp(bexp);
        }
      else {
        nsteps = 1;
        rsim->setNumInc(1);
        }

      if (dlist.getBoolean("silent", step_silent)) {
        rsim->setSilentStep(step_silent);
        }

      if (!pmSystem.useGraphics()) {
        rsim->step(nsteps);
        }
      else {
        rsim->step(1);
        }
      }

    //====== solver ======//

    else if (data.name == "solver") {
      PmRigidBodySolverType stype;

      if (dv != "") {
        PmRigidBodySolver::convSolverType(dv, stype);

        if (stype == PM_RIGID_BODY_SOLVER_UNKNOWN) {
          pm_ErrorReport (PM, "unknown solver \"%s\".", "*", dv.c_str());
          }

        PmRigidBodySolver *solver = PmRigidBodySolver::create("ode", stype);
  
        if (solver) {
          rsim->setSolver(solver);
          pm_PrintMsg (CMDPR, "%s solver set for simulation [%s]. ", dv.c_str(),
                        name.c_str());
          }
        }

      // process solver parameters //

      while (dlist.getNext(data, dv)) {
        if (!rsim->setSolverParam(data.name, dv)) {
          pm_ErrorReport (PM, "unknown solver parameter \"%s\".", "*", data.name.c_str());
          }
        }
      }

    //====== time ======//

    else if (data.name == "time") {
      dlist.getNext(data, dv);

      if (data.name == "step") {
        float dt = data.getFloat();
        rsim->setTimeStep(dt);
        }
      }

    //====== set results write parameters ======//

    else if (data.name == "write") {
      string file_name, wname;
      bool binary = true;
      bool print;

      if (!dlist.getString("file_name", file_name)) {
        pm_ErrorReport (PM, "file name must be given. ", "*");
        return;
        }

      if (dlist.getString("measurement", wname)) {
        PmMeasurement *msr;
        rsim->getMeasurement(wname, &msr);

        if (!msr) {
          pm_ErrorReport (PM, "no measurement named \"%s\". ", "*", wname.c_str());
          return;
          }

        if (dlist.getBoolean("print", print)) {
          msr->setPrintFlag(print);
          }

        msr->setWriteResultsParams(file_name, binary);
        }

      // write joint data //

      else if (dlist.getString("joint", wname)) {
        PmJoint *joint;
        rsim->getJoint(wname, &joint);

        if (!joint) {
          pm_ErrorReport (PM, "no joint named \"%s\". ", "*", wname.c_str());
          return;
          }

        joint->setWriteResultsParams(file_name, binary);
        }

      // write motor data //

      else if (dlist.getString("motor", wname)) {
        PmMotor *motor;
        rsim->getMotor(wname, &motor);

        if (!motor) {
          pm_ErrorReport (PM, "no motor named \"%s\". ", "*", wname.c_str());
          return;
          }

        motor->setWriteResultsParams(file_name, binary);
        }

      // write trace data //

      else if (dlist.getString("trace", wname)) {
        PmTrace *trace;
        rsim->getTrace(wname, &trace);

        if (!trace) {
          pm_ErrorReport (PM, "no trace named \"%s\". ", "*", wname.c_str());
          return;
          }

        trace->setWriteResultsParams(file_name, binary);
        }

      // write simulation results //

      else {
        bool bval;

        if (dlist.getBoolean("contact", bval)) {
          rsim->setWriteContact(bval);
          }

        if (dlist.getBoolean("state", bval)) {
          rsim->setWriteState(bval);
          }

        if (dlist.getBoolean("energy", bval)) {
          rsim->setWriteEnergy(bval);
          }

        if (dlist.getBoolean("kinetic_energy", bval)) {
          rsim->setWriteKineticEnergy(bval);
          }

        if (dlist.getBoolean("strain", bval)) {
          rsim->setWriteStrain(bval);
          }

        if (dlist.getBoolean("domains", bval)) {
          rsim->setWriteDomains(bval);
          }

        if (dlist.getBoolean("joints", bval)) {
          rsim->setWriteJoints(bval);
          }

        rsim->setWriteResultsParams(file_name, binary);

        if (dlist.getBoolean("current", bval)) {
          rsim->setWriteCurrent(bval);
          rsim->writeResults();
          }
        }
      }
    }
  }

//*============================================================*
//*==========              pm_CmdSimulation          ==========*
//*============================================================*
// process simualtion command.

void
pm_CmdSimulation (PmCmdDataList& dlist)
  {

  string dv, name, dname, res, type_str, sim_name;
  PmSimulation *sim;
  PmSimulationType stype;
  PmCmdData data;

  dlist.getNext (data, dv);
  //fprintf (stderr, "\n>>>>>>> pm_CmdSimulation:  [%s] \n", data.name.c_str());


  // create a new simulation and add it to the system //

  if (data.name == "create") {
    dlist.getString("name", name);
    dlist.getString("type", type_str);

    if (type_str != "") {
      PmSimulation::convSimulationType(type_str, stype);

      if (stype == PM_SIMULATION_UNKNOWN) {
        pm_ErrorReport (PM, "unknown simulation type \"%s\".", "*", type_str.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no simulation type specified.", "*");
      return;
      }

    // create a simulation of a particular type and //
    // add it to the system.                        //

    sim = PmSimulation::create(name, stype);
    pmSystem.addSimulation (sim);
    sim_name = name;
    }

  // replay the results for a simulation //

  else if (data.name == "replay") {
    int inc;

    if (dlist.getNext(data, dv)) {
      inc = data.getInt();
      }
    else {
      inc = 1;
      }

    pmSystem.setReplaySimulations(true, 0, 1, inc);
    return;
    }

  // state information //

  else if (data.name == "state") {
    if (!dlist.getNext(data, dv)) {
      return;
      }

    if (data.name == "save_frequency") {
      int save_freq = data.getInt();
      pmSystem.setSimulationsStateSaveFreq(save_freq);
      }

    return;
    }

  // step all simulations //

  else if (data.name == "step") {
    dlist.getNext (data, dv);

    if (dv != "") {
      int num_inc = data.getInt();
      pmSystem.setSimulationsNumInc(num_inc);
      }
    else {
      pmSystem.setSimulationsNumInc(1);
      //pmSystem.stepSimulations();
      }

    return;
    }
  else {
    pmSystem.getSimulation (data.name, &sim);
    }

  if (!sim) {
    pm_ErrorReport (PM, "no simulation named \"%s\".", "*", data.name.c_str());
    return;
    }

  stype = sim->getType();

  if (stype == PM_SIMULATION_RIGID) {
    pm_CmdRigidSim (sim, dlist);
    }
  }
