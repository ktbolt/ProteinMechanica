
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
//              p o t e n t i a l  c o m m a n d s           //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========         pm_CmdPotentialParams          ==========*
//*============================================================*
// process potential parameters command.

void
pm_CmdPotentialParams (vector<PmPotential*> pots, PmPotentialType ptype, 
                       PmCmdDataList& dlist)
  {
  fprintf (stderr, "\n>>>>>> pm_CmdPotentialParams \n");
  string dv;
  PmCmdData data;
  PmVector3 color;
  vector<string> map_scale;
  float map_min, map_max, width;
  bool show, show_set, map_strain, map_energy, has_map_scale, color_set; 

  //===== process axis potential parameters =====//

  if (ptype == PM_POTENTIAL_AXIS) { 
    float rest_strength, rest_cutoff, val;
    vector<float> fconsts;
    string str;
    vector<string> fclist;

    if (!dlist.getFloat("rest_strength", rest_strength)) {
      rest_strength = 1.0;            
      }

    if (dlist.getFloat("rest_cutoff", rest_cutoff)) {
       rest_cutoff = 0.0;
       }

    if (dlist.getStringList("axes_strength", fclist)) {
      for (unsigned int i = 0; i < fclist.size(); i++) {
        str = fclist[i];

        if (!convToFloat(str, val)) {
          pm_ErrorReport (PM, "improper force constant \"%s\".", "*", str.c_str());
          pmSystem.setCmdError(true);
          return;
          }

        fconsts.push_back(val);
        }
      }

    for (unsigned int i = 0; i < pots.size(); i++) {
      PmAxisPotential *axpot = dynamic_cast<PmAxisPotential*>(pots[i]);
      axpot->setRestStrength(rest_strength);
      axpot->setRestCutoff(rest_cutoff);

      if (fconsts.size()) {
        axpot->setAxesStrength(fconsts);
        }
      }
    }

  //===== process bend potential parameters =====//

  else if (ptype == PM_POTENTIAL_BEND) { 
    float force_const;

    if (!dlist.getFloat("force_const", force_const)) {
      force_const = 1.0;
      }

    for (unsigned int i = 0; i < pots.size(); i++) {
      PmBendPotential *bpot = dynamic_cast<PmBendPotential*>(pots[i]);
      bpot->setForceConst(force_const);
      }
    }

  //===== process contact potential parameters =====//

  else if (ptype == PM_POTENTIAL_CONTACT) {
    float strength;

    if (!dlist.getFloat("strength", strength)) {
      strength = 1.0;
      }

    for (unsigned int i = 0; i < pots.size(); i++) {
      PmContactPotential *cpot = dynamic_cast<PmContactPotential*>(pots[i]);
      cpot->setContactStrength(strength);
      }
    }

  //===== process entropic spring potential parameters =====//

  else if (ptype == PM_POTENTIAL_ENTROPIC_SPRING) { 
    string str;
    PmEntropicSpringFuncType ftype;
    float val;
    vector<float> graph_data;
    float persistence_length=0.0, contour_length=0.0;

    // force-extension function //

    if (dlist.getString("function", str)) {
      PmEntropicSpringPotential::convFunctionType(str, ftype);

      if (ftype == PM_ENTROPIC_SPRING_FUNC_UNKNOWN) {
        pm_ErrorReport (PM, "unknown entropic spring function type \"%s\".", "*",
                        str.c_str());
        pmSystem.setCmdError(true);
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no entropic spring function type given", "*");
      pmSystem.setCmdError(true);
      return;
      }

    // graph //

    if (ftype == PM_ENTROPIC_SPRING_FUNC_GRAPH) {
      vector<string> list;

      if (dlist.getStringList("graph_data", list)) {
        if ((list.size() % 2) != 0) {
          pm_ErrorReport (PM, "graph data must contain even number of elements.", "*");
          pmSystem.setCmdError(true);
          return;
          }

        for (unsigned int i = 0; i < list.size(); i++) {
          str = list[i];

          if (!convToFloat(str, val)) {
            pm_ErrorReport (PM, "invalid graph element %s", "*", str.c_str());
            pmSystem.setCmdError(true);
            return;
            }

          graph_data.push_back(val);
          }
        }
      }

    // worm-like chain //

    else if (ftype == PM_ENTROPIC_SPRING_FUNC_WLC) {
      int num_res;

      if (!dlist.getFloat("persistence_length", persistence_length)) {
        persistence_length = 0.0;
        }

      if (!dlist.getFloat("contour_length", contour_length)) {
        contour_length = 0.0;
        }

      if (dlist.getInt("num_res", num_res)) {
        contour_length = num_res * 0.34;
        }
      }

    // set parameters for each potential //

    for (unsigned int i = 0; i < pots.size(); i++) {
      PmEntropicSpringPotential*epot = dynamic_cast<PmEntropicSpringPotential*>(pots[i]);
      epot->setForceFunctionType(ftype);

      if (graph_data.size()) {
        epot->setForceGraphData(graph_data);
        }

      if ((persistence_length != 0.0) && (contour_length != 0.0)) {
        epot->setWlcParams(persistence_length, contour_length);
        }
      }
    }

  //===== process molecular mechanics potential parameters =====//

  else if (ptype == PM_POTENTIAL_MOLECULAR_MECHANICS) { 
    vector<string> slist;
    PmPotentialMolMechTerms terms;
    PmPotentialMolMechParameters params;
    float maximum_force, cutoff, strength;
    bool report_energy;

    if (!dlist.getFloat("maximum_force", maximum_force)) {
      maximum_force = 100.0;
      }

    if (!dlist.getFloat("cutoff", cutoff)) {
      cutoff = 1.0;
      }

    if (!dlist.getFloat("strength", strength)) {
      strength = 1.0;
      }

    if (!dlist.getBoolean("report_energy", report_energy)) {
      report_energy = true;
      }

    if (dlist.getStringList("terms", slist)) {
      for (unsigned int i = 0; i < slist.size(); i++) {
        if (!terms.setTerm(slist[i])) {
          pm_ErrorReport (PM, "unknown term name = \"%s\".", "*", slist[i].c_str());
          pmSystem.setCmdError(true);
          return;
          }
        }
      }

    if (dlist.getStringList("parameters", slist)) {
      for (unsigned int i = 0; i < slist.size(); i++) {
        if (!params.setParameter(slist[i])) {
          pm_ErrorReport (PM, "unknown parameter name = \"%s\".", "*", slist[i].c_str());
          pmSystem.setCmdError(true);
          return;
          }
        }
      }

    for (unsigned int i = 0; i < pots.size(); i++) {
      PmMolecularMechanicsPotential *mpot = 
         dynamic_cast<PmMolecularMechanicsPotential*>(pots[i]);
      mpot->setTerms(terms);
      mpot->setParameters(params);
      mpot->setMaxForce(maximum_force);
      mpot->setCutoff(cutoff);
      mpot->setStrength(strength);
      mpot->setReportEnergy(report_energy);
      }
    }

  //===== process spring potential parameters =====//

  else if (ptype == PM_POTENTIAL_SPRING) { 
    fprintf (stderr, "\n>>>>>> add spring potential \n");
    float cutoff, force_const, max_strain, max_dist;
    bool ljspring;

    if (!dlist.getFloat("cutoff", cutoff)) {
      cutoff = 1.0;
      }

    if (!dlist.getFloat("max_strain", max_strain)) {
      max_strain = 0.0;
      }

    if (!dlist.getFloat("max_distance", max_dist)) {
      max_dist = 0.0;
      }

    if (!dlist.getFloat("force_const", force_const)) {
      force_const = 1.0;
      }

    if (!dlist.getBoolean("ljspring", ljspring)) {
      ljspring = false;
      }

    for (unsigned int i = 0; i < pots.size(); i++) {
      PmSpringPotential *spot = dynamic_cast<PmSpringPotential*>(pots[i]);
      spot->setCutoff(cutoff);
      spot->setMaxStrain(max_strain);
      spot->setMaxDistance(max_dist);
      spot->setForceConst(force_const);
      spot->setLJSpring(ljspring);
      }
    }

  //===== process torsion potential =====//

  else if (ptype == PM_POTENTIAL_TORSION) {
    float force_const;

    if (!dlist.getFloat("force_const", force_const)) {
      force_const = 1.0;
      }

    for (unsigned int i = 0; i < pots.size(); i++) {
      PmTorsionPotential *tpot = dynamic_cast<PmTorsionPotential*>(pots[i]);
      tpot->setForceConst(force_const);
      }
    }

  //===== process graphics attributes =====//

  has_map_scale = false;
  map_min = map_max = 0.0;
  show_set = false; 
  map_strain = false; 
  map_energy = false; 
  has_map_scale = false; 
  color_set = false; 
  width = 1.0;

  while (dlist.getNext(data, dv)) {
    if (data.name == "color") {
      data.getVector(color);
      color_set = true; 
      }
    else if (data.name == "map_energy") {
      map_energy = data.getBoolean();
      }
    else if (data.name == "map_scale") {
      data.getStringList(map_scale);
      has_map_scale = true;

      if (!convToFloat(map_scale[0], map_min)) {
        pm_ErrorReport (PM, "bad map scale value \"%s\" ", "*", map_scale[0].c_str());
        return;
        }

      if (!convToFloat(map_scale[1], map_max)) {
        pm_ErrorReport (PM, "bad map scale value \"%s\" ", "*", map_scale[1].c_str());
        return;
        }
      }
    else if (data.name == "map_strain") {
      map_strain = data.getBoolean();
      }
    else if (data.name == "show") {
      show = data.getBoolean();
      show_set = true;
      }
    else if (data.name == "width") {
      width = data.getFloat();
      }
    }

  for (unsigned int i = 0; i < pots.size(); i++) {
    if (color_set) {
      pots[i]->setColor(color);
      }

    pots[i]->setLineWidth(width);

    if (ptype == PM_POTENTIAL_SPRING) {
      PmSpringPotential *spot = dynamic_cast<PmSpringPotential*>(pots[i]);
      spot->setMapStrain(map_strain);
      spot->setMapEnergy(map_energy);

      if (has_map_scale) {
        spot->setMapScale(map_min, map_max);
        }
      }

    if (show_set) {
      pots[i]->displayGeometry(show);
      }
    }
  }

//*============================================================*
//*==========             pm_CmdPotential            ==========*
//*============================================================*
// process potential command.

void
pm_CmdPotential (PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, dname, res, type_str, rgn, gname;
  PmPotentialType ptype;
  float scale, ramp;
  PmCmdData pdata;
  PmPotential *pot;
  PmPotentialGeom *pgeom1, *pgeom2;
  PmVector3 bpos, color;
  PmPotentialGeomType gtype;
  PmPotentialParameters params;
  vector<PmPotential*> pots;

  dlist.getNext (data);

  // create a new potential and add it to the system //

  if (data.name == "create") {
    dlist.getString("name", name);
    dlist.getString("type", type_str);

    if (type_str != "") {
      PmPotential::convType(type_str, ptype);

      if (ptype == PM_POTENTIAL_UNKNOWN) {
        pm_ErrorReport (PM, "unknown potential type \"%s\".", "*", type_str.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no potential type specified.", "*");
      return;
      }

    // process potential geometry //

    if (dlist.getString("geometry1", gname)) {
      pmSystem.getPotentialGeometry(gname, &pgeom1);

      if (!pgeom1) {
        pm_ErrorReport (PM, "unknown potential geometry type \"%s\".", "*", gname.c_str());
        pmSystem.setCmdError(true);
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no potential geometry1 given.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    if (dlist.getString("geometry2", gname)) {
      pmSystem.getPotentialGeometry(gname, &pgeom2);

      if (!pgeom2) {
        pm_ErrorReport (PM, "unknown potential geometry type \"%s\".", "*", gname.c_str());
        pmSystem.setCmdError(true);
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no potential geometry2 given.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    // create a potential of a particular type and add it to the system //

    pot = PmPotential::create(name, ptype, pgeom1, pgeom2);
    pmSystem.addPotential(pot);
    }

  // process potential parameters //

  pots.push_back(pot);

  pm_CmdPotentialParams (pots, ptype, dlist);


  // process graphics attributes? //

/*
  bool show, show_set = false;

  while (dlist.getNext(data, dv)) {
    if (data.name == "color") {
      if (data.getVector(color)) {
        pot->setColor(color);
        }
      }

    else if (data.name == "show") {
      show = data.getBoolean();
      show_set = true;
      }
    }

  if (show_set) {
    pot->displayGeometry(show);
    }
*/
  }

