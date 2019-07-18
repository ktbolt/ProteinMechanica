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
//              j o i n t   c o m m a n d s                  //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========         pm_CmdJointBondAxis            ==========*
//*============================================================*
// process joint bond axis command.

void
pm_CmdJointBondAxis (PmJoint *joint, PmCmdData& data, PmCmdDataList& dlist)
  {
  string desc, dname, rname, aname, str;
  PmAtomFilter filter;
  PmVector3 axis, axis1, axis2, axis3, jloc;
  vector<string> bond_names;
  vector<string> names(3);

  // get atoms names if they are defined  //

  dlist.getStringList("atom_names", filter.names);

  // get joint type //

  //PmJointType jtype = joint->getType();

  if (!dlist.getStringList("bond", bond_names)) {
    pm_ErrorReport (PM, "bond not given.", "*");
    return;
    }

  // get the atoms for the bonds //

  int num_bonds = bond_names.size();
  PmVector3 pos[3];

  for (int i = 0; i < num_bonds && i < 3; i++) {
    str = bond_names[i];
    int dnum = 0;

    for (unsigned int j = 0; j < str.size() && dnum < 3; j++) {
      char c = str[j];

      if (c == ':') {
        dnum += 1;
        }
      else {
        names[dnum].push_back(c);
        }
      }

    if (dnum != 2) {
      pm_ErrorReport (PM, "bad format for bond.", "*");
      return;
      }

    PmMolecule *domain;
    pmSystem.getDomain (names[0], &domain);

    if (!domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", names[0].c_str());
      return;
      }

    if (!domain->getResidueCoords(names[1], names[2], pos[i])) { 
      pm_ErrorReport (PM, "can't get atoms for domain \"%s\".", "*", names[0].c_str());
      return;
      }

    //pm_PrintMsg (CMDPR, "pos = [ %g, %g, %g ] \n", pos[i][0], pos[i][1], pos[i][2]);

    for (int j = 0; j <= dnum; j++) {
      names[j].clear();
      }
    }

  // set hinge joint axis //

  if (PmHingeJoint *pjoint = dynamic_cast<PmHingeJoint*>(joint)) {
    axis = pos[0] - pos[1];

    if (axis.length() == 0.0) {
      pm_ErrorReport (PM, "axis length is zero.", "*");
      return;
      }

    pjoint->setAxis(axis);
    //pm_PrintMsg (CMDPR, "axis = [ %g, %g, %g ] \n", axis[0], axis[1], axis[2]);
    pjoint->setPosition2(pos[1]);
    }
  }

//*============================================================*
//*==========         pm_CmdJointRange               ==========*
//*============================================================*
// process joint range command. 

void
pm_CmdJointRange(PmJoint *joint, PmCmdData& data, PmCmdDataList& dlist)
  {
  string desc, dname, str;
  float r = data.getFloat();
  PmJointType jtype = joint->getType();

  if (jtype == PM_JOINT_UNIVERSAL) {
    //PmUniversalJoint *ujoint = dynamic_cast<PmUniversalJoint*>(joint);

    if (data.name == "range1") {
      //ujoint->setRange(1,k);
      }
    else if (data.name == "range2") {
      //ujoint->setRange(2,r);
      }
    }

  else if (jtype == PM_JOINT_HINGE) {
    PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);
    hjoint->setRange(r);
    }

  else if (jtype == PM_JOINT_BALL) {
    PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint);

    if (data.name == "range1") {
      bjoint->setRange(1, r);
      }
    else if (data.name == "range2") {
      bjoint->setRange(2, r);
      }
    else if (data.name == "range3") {
      bjoint->setRange(3, r);
      }
    }
  }

//*============================================================*
//*==========         pm_CmdJointForceConst          ==========*
//*============================================================*
// process joint force constant command. this command is only 
// valid for certain types of joints. 

void
pm_CmdJointForceConst (PmJoint *joint, PmCmdData& data, PmCmdDataList& dlist)
  {
  string desc, dname, str;
  float k = data.getFloat();
  PmJointType jtype = joint->getType();

  if (jtype == PM_JOINT_UNIVERSAL) {
    PmUniversalJoint *ujoint = dynamic_cast<PmUniversalJoint*>(joint);

    if (data.name == "force_const1") {
      ujoint->setForceConst(1,k);
      }
    else if (data.name == "force_const2") {
      ujoint->setForceConst(2,k);
      }
    }

  else if (jtype == PM_JOINT_HINGE) {
    PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);
    hjoint->setForceConst(k);
    }

  else if (jtype == PM_JOINT_BALL) {
    PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint);

    if (data.name == "force_const1") {
      bjoint->setForceConst(1, k);
      }
    else if (data.name == "force_const2") {
      bjoint->setForceConst(2, k);
      }
    else if (data.name == "force_const3") {
      bjoint->setForceConst(3, k);
      }
    }
  }

//*============================================================*
//*==========         pm_CmdJointPcaAxis             ==========*
//*============================================================*
// process joint pca axis command. 

void
pm_CmdJointPcaAxis (PmJoint *joint, PmCmdData& data, PmCmdDataList& dlist)
  {
  string desc, dname, str;
  PmAtomFilter filter;
  PmMolecule *domain;
  PmPcaResults pca;
  vector<PmVector3> coords;
  PmVector3 axis, axis1, axis2, axis3;
  float pca_rot, axes_sign;
  bool reverse_axes = false;

  // get atoms names if they are defined  //

  dlist.getStringList("atom_names", filter.names);

  // get joint type //

  //PmJointType jtype = joint->getType();

  // get domain //

  if ((dlist.getString("pca_domain", dname) || dlist.getString("domain", dname))) {
    pmSystem.getDomain (dname, &domain);

    if (!domain) {
      pm_ErrorReport (PM, "no pca domain named \"%s\".", "*", dname.c_str());
      return;
      }
    }

  // process pca residues //

  if (!dlist.getString("pca_res", desc)) {
    pm_ErrorReport (PM, "no pca residues given.", "*");
    return;
    }

  vector<string> names(3);
  char c;
  int n = 0;

  for (unsigned int i = 0; i < desc.size(); i++) {
    c = desc[i];

    if (c == ':') {
      n = n + 1;
      }
    else {
      names[n].push_back(c);
      }
    }

  if (n != 0) {
    dname = names[0];
    pmSystem.getDomain (dname, &domain);

    if (!domain) {
      pm_ErrorReport (PM, "no pca domain named \"%s\".", "*", dname.c_str());
      return;
      }

    if (n >= 1) {
      desc = names[1];
      }

    if (n == 2) {
      filter.names.push_back(names[2]);
      }
    }

  if (domain) {
    domain->getAtomCoords(desc, filter, coords);
    }

  if (coords.size() == 0) {
    pm_ErrorReport (PM, "domain \"%s\" contains no atoms specified.", "*", dname.c_str());
    return;
    }

  pm_MathPrincipalComponents(coords, pca);

  // get rotation about pca axis 1 //

  if (!dlist.getFloat("pca_rot", pca_rot)) { 
    pca_rot = 0.0;
    }

  dlist.getBoolean("reverse_axes", reverse_axes);

  if (reverse_axes) {
    axes_sign = -1.0;
    }
  else {
    axes_sign = 1.0;
    }

  // set appropriate axes //

  if (dlist.getString("axis", str)) {
    if (str == "pca1") {
      axis = axes_sign*pca.axis1;
      }
    else if (str == "pca2") {
      axis = axes_sign*pca.axis2;
      }
    else if (str == "pca3") {
      axis = axes_sign*pca.axis3;
      }
    else {
      pm_ErrorReport (PM, "unknown axis \"%s\" specified for pca.", "*", str.c_str());
      return;
      }
    }

  if (dlist.getString("axis1", str)) {
    if (str == "pca1") {
      axis1 = axes_sign*pca.axis1;
      }
    else if (str == "pca2") {
      axis1 = axes_sign*pca.axis2;
      }
    else if (str == "pca3") {
      axis1 = axes_sign*pca.axis3;
      }
    else {
      pm_ErrorReport (PM, "unknown axis1 \"%s\" specified for pca.", "*", str.c_str());
      return;
      }
    }

  if (dlist.getString("axis2", str)) {
    if (str == "pca1") {
      axis2 = axes_sign*pca.axis1;
      }
    else if (str == "pca2") {
      axis2 = axes_sign*pca.axis2;
      }
    else if (str == "pca3") {
      axis2 = axes_sign*pca.axis3;
      }
    }

  if (pca_rot != 0.0) {
    PmMatrix3x3 rmat;
    PmVector3 v1, v2; 
    pm_MathRotationAroundAxis (pca_rot, pca.axis1, rmat);

    if (axis1.length() != 0.0) {
      v1 = rmat*axis1;
      v2 = rmat*axis2;
      axis1 = v1;
      axis2 = v2;
      }

    else if (axis.length() != 0.0) {
      v1 = rmat*axis;
      axis = v1;
      }
    }

  // set hinge joint axis //

  if (PmHingeJoint *pjoint = dynamic_cast<PmHingeJoint*>(joint)) {
    if (axis.length() == 0.0) {
      pm_ErrorReport (PM, "no axis given for pca.", "*");
      return;
      }

    pjoint->setAxis(axis);
    pm_PrintMsg (CMDPR, "axis = (%g, %g, %g) ", axis[0], axis[1], axis[2]);
    }

  // set universal joint axes //

  else if (PmUniversalJoint *ujoint = dynamic_cast<PmUniversalJoint*>(joint)) {
    if (axis1.length() == 0.0) {
      pm_ErrorReport (PM, "no axis1 given for pca.", "*");
      return;
      }

    if (axis2.length() == 0.0) {
      pm_ErrorReport (PM, "no axis2 given for pca.", "*");
      return;
      }

    ujoint->setAxis(1,axis1);
    ujoint->setAxis(2,axis2);
    pm_PrintMsg (CMDPR, "axis 1 = (%g %g %g) ", axis1[0], axis1[1], axis1[2]);
    pm_PrintMsg (CMDPR, "axis 2 = (%g %g %g) ", axis2[0], axis2[1], axis2[2]);
    }

  // set ball joint axes //

  else if (PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint)) {
    bjoint->setAxes(axis1, axis2);
    }
  }

//*============================================================*
//*==========         pm_CmdJointAxis                ==========*
//*============================================================*
// process joint axis command. the axis command is only valid 
// for certain types of joints. 

void
pm_CmdJointAxis (PmJoint *joint, PmCmdData& data, PmCmdDataList& dlist)
  {
  string desc, dname, str;
  PmVector3 axis;

  data.getVector(axis);
  PmJointType jtype = joint->getType();

  if (jtype == PM_JOINT_UNIVERSAL) {
    PmUniversalJoint *ujoint = dynamic_cast<PmUniversalJoint*>(joint);

    if (!ujoint) { 
      pm_ErrorReport (PM, "can't set axis data for joint: wrong type", "*");
      return;
      }
  
    if (data.name == "axis1") {
      ujoint->setAxis(1,axis);
      }
    else if (data.name == "axis2") {
      ujoint->setAxis(2,axis);
      }
    }

  else if (jtype == PM_JOINT_HINGE) {
    PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);
    hjoint->setAxis(axis);
    }

  else if (jtype == PM_JOINT_BALL) {
    PmVector3 axis1, axis2;
    PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint);

    if (data.name == "axis1") {
      axis1 = axis;

      if (!dlist.getVector("axis2", axis2)) {
        pm_ErrorReport (PM, "no axis2 given.", "*");
        return;
        }
      }

    else if (data.name == "axis2") {
      axis2 = axis;

      if (!dlist.getVector("axis1", axis1)) {
        pm_ErrorReport (PM, "no axis2 given.", "*");
        return;
        }
      }

    else { 
      pm_ErrorReport (PM, "unknown axis \"%s\".", "*", data.name.c_str());
      return;
      }

    /*
    fprintf (stderr, "\n>>> data name = %s \n", data.name.c_str()); 
    fprintf (stderr, ">>> axis1 = %f %f %f \n", axis1[0], axis1[1], axis1[2]);
    fprintf (stderr, ">>> axis2 = %f %f %f \n", axis2[0], axis2[1], axis2[2]);
    */
    bjoint->setAxes(axis1, axis2);
    }
  }

//*============================================================*
//*==========         pm_CmdJoint                    ==========*
//*============================================================*
// process joint command.

void
pm_CmdJoint (PmCmdDataList& dlist)
  {
  string dv, name, dname, res, type_str, use_type, pstr;
  bool use_pca, use_bond;
  PmJoint *joint;
  PmJointType jtype;
  PmVector3 jpos;
  PmMolecule *domain;
  PmCmdData data, pdata;

  dlist.getNext (data, dv);

  // create a new joint and add it to the system //

  if (data.name == "create") {
    dlist.getString("name", name);
    dlist.getString("type", type_str);

    if (type_str != "") {
      PmJoint::convJointType(type_str, jtype);

      if (jtype == PM_JOINT_UNKNOWN) {
        pm_ErrorReport (PM, "unknown joint type \"%s\".", "*", type_str.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no joint type specified.", "*");
      return;
      }

    // create a joint of a particular type and //
    // add it to the system.                   //

    joint = PmJoint::create(name, jtype);
    pmSystem.addJoint (joint);

    //===== define joint position =====//

    if (dlist.getString("domain", dname)) {
      pmSystem.getDomain (dname, &domain);

      if (!domain) {
        pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
        return;
        }

      if (!dlist.getString("res", res) && !dlist.getString("residues", res)) {
        pm_ErrorReport (PM, "no residues given", "*");
        return;
        }

      dlist.getString("use", use_type);

      if (!res.empty() && !domain) {
        pm_ErrorReport (PM, "no domain specified for res \"%s\".", "*", res.c_str());
        return;
        }
      else if ((res != "") && domain) {
        if (!domain->getResidueCoords(res, use_type, jpos)) {
          pm_ErrorReport (PM, "invalid res \"%s\".", "*", res.c_str());
          return;
          }
        }

      joint->setPosition (jpos);
      pm_PrintMsg (CMDPR, "joint pos = (%g %g %g) ", jpos[0], jpos[1], jpos[2]);
      }

    else if (dlist.getString("cylinder", dname)) {
      PmCylinder *cyl;
      PmSolid::getCylinder(dname, &cyl);
      PmVector3 center, axis, pos;
      float radius, length, dist;

      if (!cyl) {
        pm_ErrorReport (PM, "no cylinder named \"%s\".", "*", dname.c_str());
        return;
        }

      if (!dlist.getFloat("distance", dist)) {
        pm_ErrorReport (PM, "no axis distance specified.", "*");
        return;
        }

      if (dist < -1.0) dist = -1.0;
      if (dist >  1.0) dist =  1.0;

      cyl->getParameters(center, radius, length, axis);
      jpos = center + dist*length*axis / 2.0;
      joint->setPosition (jpos);
      }

    else if (dlist.getString("position", pstr)) {
      dlist.getData("position", pdata);

      if (!pm_CmdProcPosition(pdata, jpos)) {
        pm_ErrorReport (PM, "couldn't process joint postion data \"%s\".", "*",
                        pstr.c_str());
        return;
        }

      //fprintf (stderr, " joint  pos = %f %f %f \n", jpos[0], jpos[1], jpos[2]);
      joint->setPosition (jpos);
      }

    else if (dlist.getVector("position", jpos)) {
      joint->setPosition (jpos);
      }

    // define joint position from shared residue in bodies //

    else {
      vector<string> blist;
      PmMolecule *dom1, *dom2;
      PmBody *body1, *body2;
      PmAtomFilter filter;
      vector<PmAtom> atoms1, atoms2;
      bool overlap;

      if (!dlist.getStringList("bodies", blist)) {
        pm_ErrorReport (PM, "no joint position given.", "*");
        }

      if (blist.size() != 2) {
        pm_ErrorReport (PM, "must specify two bodies.", "*");
        return;
        }

      pmSystem.getBody (blist[0], &body1);

      if (!body1) {
        pm_ErrorReport (PM, "no body named \"%s\".", "*", blist[0].c_str());
        return;
        }

      pmSystem.getBody (blist[1], &body2);

      if (!body2) {
        pm_ErrorReport (PM, "no body named \"%s\".", "*", blist[1].c_str());
        return;
        }

      body1->getMolecule(&dom1);
      body2->getMolecule(&dom2);

      if (!dom1 || !dom2) {
        pm_ErrorReport (PM, "bodies are not derived from domains.", "*");
        return;
        }

      filter.names.push_back("CA");
      dom1->getAtoms(filter, atoms1);
      dom2->getAtoms(filter, atoms2);
      overlap = false;

      for (unsigned int i = 0; i < atoms1.size(); i++) {
        for (unsigned int j = 0; j < atoms2.size(); j++) {
          if (atoms1[i].seq == atoms2[j].seq) {
            overlap = true;
            jpos = atoms1[i].pos;
            joint->setPosition (jpos);
            //fprintf (stderr, ">>> position=(%f, %f, %f) \n", jpos[0], jpos[1], jpos[2]); 
            break;
            }
          }
        }

      if (!overlap) {
        pm_ErrorReport (PM, "body domains don't overlap.", "*");
        return;
        }

      joint->setBodies (body1, body2);
      }

    // define joint axes using a bond //

    if (dlist.getBoolean("use_bond", use_bond)) {
      pm_CmdJointBondAxis (joint, data, dlist);
      }

    // define joint axes using pca //

    else if (dlist.getBoolean("use_pca", use_pca)) {
      pm_CmdJointPcaAxis (joint, data, dlist);
      }
    }
  else {
    pmSystem.getJoint (data.name, &joint);
    }

  if (!joint) {
    pm_ErrorReport (PM, "no joint named \"%s\".", "*", data.name.c_str());
    return;
    }

  // process joint parameters //

  bool show = false , show_set = false, show_axes = false;

  while (dlist.getNext(data)) {
    if (!use_pca && data.name.find("axis") != string::npos) {
      pm_CmdJointAxis (joint, data, dlist);
      }
    else if (data.name == "bodies") {
      vector<string> blist;
      PmBody *body1, *body2;
      dlist.getStringList ("bodies", blist);
   
      if (blist.size() != 2) {
        pm_ErrorReport (PM, "must specify two bodies.", "*");
        return;
        }

      pmSystem.getBody (blist[0], &body1);

      if (!body1) {
        pm_ErrorReport (PM, "no body named \"%s\".", "*", blist[0].c_str());
        return;
        }

      pmSystem.getBody (blist[1], &body2);

      if (!body2) {
        pm_ErrorReport (PM, "no body named \"%s\".", "*", blist[1].c_str());
        return;
        }

      joint->setBodies (body1, body2);
      }

    else if (data.name == "color") {
      PmVector3 color;
      data.getVector(color);
      joint->setColor(color);
      }

    else if (data.name == "display") {
      data.getString(dv);
      PmGeometryDisplayType dtype;
      PmGraphicsGeometry::convDisplayType(dv, dtype);
      joint->setDisplayType(dtype);
      }

    // look for force_const, force_const1, ... //

    else if (data.name.find("force_const") != string::npos) {
      pm_CmdJointForceConst(joint, data, dlist);
      }

    else if (data.name == "msize") {
      float msize = data.getFloat(); 
      joint->setMsize (msize);
      }

    else if (data.name.find("range") != string::npos) {
      pm_CmdJointRange(joint, data, dlist);
      }

    else if (data.name == "print") {
      string name, str1, str2;
      PmBody *body1, *body2;
      PmVector3 pos;
      PmJointType jtype;

      joint->getName(name);
      joint->getType(str1);
      jtype = joint->getType();
      joint->getBodies(&body1, &body2);
      joint->getPosition(pos);

      fprintf (stderr, "\n------ joint \"%s\" ------ \n", name.c_str()); 
      fprintf (stderr, "type = \"%s\" \n", str1.c_str()); 
      fprintf (stderr, "position = (%f, %f, %f) \n", pos[0], pos[1], pos[2]); 

      if (!body1 || !body2) {
        fprintf (stderr, ">>> bodies: *** not defined *** \n");
        }
      else {
        body1->getName(str1);
        body2->getName(str2);
        fprintf (stderr, "bodies = \"%s\"  \"%s\" \n", str1.c_str(), str2.c_str());
        }

      if (jtype == PM_JOINT_UNIVERSAL) {
        PmVector3 axis1, axis2;
        PmUniversalJoint *ujoint = dynamic_cast<PmUniversalJoint*>(joint);
        ujoint->getAxis(1, axis1);
        ujoint->getAxis(2, axis2);
        fprintf (stderr, "axis1 = (%f, %f, %f) \n", axis1[0], axis1[1], axis1[2]); 
        fprintf (stderr, "axis2 = (%f, %f, %f) \n", axis2[0], axis2[1], axis2[2]); 
        }

      else if (jtype == PM_JOINT_HINGE) {
        PmVector3 axis;
        PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);
        hjoint->getAxis(axis);
        fprintf (stderr, "axis = (%f, %f, %f) \n", axis[0], axis[1], axis[2]); 
        }
      }

    else if (data.name == "shading") {
      data.getString(dv);
      PmGeometryShadingType stype;
      PmGraphicsGeometry::convShadingType(dv, stype);
      joint->setShadingType(stype);
      }

    else if (data.name == "show") {
      show = data.getBoolean();
      show_set = true; 
      }

    else if (data.name == "show_axes") {
      show_axes = data.getBoolean();
      joint->setShowAxes(show_axes);
      }
    }

  if (show_set) {
    joint->setShow(show);
    joint->display(show);
    }
  }

