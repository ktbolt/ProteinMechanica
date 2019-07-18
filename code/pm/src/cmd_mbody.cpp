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
//          m u l t i b o d y     c o m m a n d s            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========          pm_CmdMultiBodyCalpha         ==========*
//*============================================================*
// process multibody command for a c-alpha chain.

void
pm_CmdMultiBodyCalpha (PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, dname, desc;
  PmMolecule *domain, *mol;

  int num_res;
  vector<PmResidue*> rlist;
  PmResidue *res, *nres;
  char chain_id;
  stringstream dss;
  PmMolecule *ca_domain; 
  PmAtomFilter sidechain_filter;
  PmAtom *atom;

  PmExtent extent;
  PmVector3 color;

  string bname, jname, pstr, nstr; 
  PmBody *body;
  PmJoint *joint;
  PmVector3 jpos, pos, axis;
  vector <PmBody*> bodies;
  vector <PmJoint*> joints;
  PmJointType joint_type;
  bool verbose = pmSystem.getCmdVerbose();
  bool first_peptide=true;

  if (verbose) {
    fprintf (stderr, "\n------ create c-alpha chain ------ \n");
    }

  if (!dlist.getString("name", name)) {  
    pm_ErrorReport (PM, "no name given.", "*");
    return;
    }

  //===== get current molecule =====//

  pmSystem.getMolecule (&mol);

  if (!mol) {
    pm_ErrorReport (PM, "no current molecule.", "*");
    return;
    }

  //===== get domain =====//

  if (dlist.getString("domain", dname)) {  
    pmSystem.getDomain (dname, &domain);

    if (!domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
      return;
      }
    }
  else {
    pm_ErrorReport (PM, "no domain name given.", "*");
    return;
    }

  joint_type = PM_JOINT_BALL;


  //===== get domain residues =====//

  domain->getResidues(desc, rlist);
  num_res = rlist.size();

  if (verbose) {
    fprintf (stderr, ">>> number of residues=%d \n", num_res);
    }

  sidechain_filter.sidechain_group = true;

  //===== create domains and rigid bodies =====//

  vector<PmMolecule*> ca_domains; 

  if (verbose) {
    fprintf (stderr, "\n------ create domains and bodies ------ \n");
    }

  for (int i = 0; i < num_res-1; i++) {
    res = rlist[i];
    nres = rlist[i+1];
    chain_id = res->getChainId();

    dss << name << "Ca" << res->id;
    dname = dss.str();
    dss.str(std::string());
    dss << chain_id << '[' << res->id << ']';
    desc = dss.str();
    dss.str(std::string());

    if (verbose) {
      fprintf (stderr, ">>> ca domain=%s residues=%s",dname.c_str(),desc.c_str());
      }

    ca_domain = mol->createDomain (dname, desc, sidechain_filter);

    if (!ca_domain) {
      string mol_name;
      mol->getName(mol_name);
      pm_ErrorReport (PM, "can't create domain with desc=%s from mol=%s ", "*", 
                      desc.c_str(), mol_name.c_str());
      return;
      }

    // add peptide and next residue N and CA atoms //

    res->getAtom ("O", &atom);
    ca_domain->addAtom (atom);
    res->getAtom ("C", &atom);
    ca_domain->addAtom (atom);
    nres->getAtom ("N", &atom);
    ca_domain->addAtom (atom);
    nres->getAtom ("CA", &atom);
    ca_domain->addAtom (atom);

    ca_domain->getExtent(extent);
    pmSystem.updateExtent (extent);

    if (i % 2) {
      color.set(0,1,0);
      }
    else {
      color.set(1,0,0);
      }

    ca_domain->setAtomColor(color);
    ca_domain->setBondAtomRenderType(PM_MOLECULE_RENDER_LINE);
    ca_domain->setDisplayBondAtoms (true);
    ca_domain->displayBonds("", true);

    vector<PmAtom> satoms;
    ca_domain->getAtoms(satoms);

    if (verbose) {
      fprintf (stderr, " atoms={");
      for (int ia = 0; ia < satoms.size(); ia++) {
        fprintf (stderr, "%d:%s ", satoms[ia].seq, satoms[ia].name); 
        }
      fprintf (stderr, "} ");
      }

    ca_domains.push_back(ca_domain);
    dss << name << "Ca" << res->id << "Body";
    bname = dss.str();
    dss.str(std::string());

    if (verbose) {
      fprintf (stderr, " body=%s \n", bname.c_str());
      }

    body = new PmBody(bname, PM_BODY_RIGID);
    pmSystem.addBody (body);
    body->setPhysicalObj(ca_domain);
    bodies.push_back(body);
    }

  //===== create joints =====//

  PmMolecule *next_ca_domain;
  PmVector3 jc;
  int res_count=0;

  if (verbose) {
    fprintf (stderr, "\n------ create joints ------ \n");
    }

  for (int i = 0; i < num_res-2; i++) {
    res = rlist[i];
    nres = rlist[i+1];
    chain_id = res->getChainId();

    ca_domain = ca_domains[i];

    dss << name << "CAJnt" << i+1;
    jname = dss.str();
    dss.str(std::string());
    joint = PmJoint::create(jname, PM_JOINT_BALL);
    pmSystem.addJoint (joint);

    dss << chain_id << '[' << nres->id << ']';
    pstr = dss.str();
    dss.str(std::string());

    if (!ca_domain->getResidueCoords(pstr, "CA", pos)) { 
      fprintf (stderr, "*** no CA for residue %s \n", pstr.c_str()); 
      }

    if (verbose) {
      fprintf (stderr, ">>> create joint=%s position=%s \n", jname.c_str(),
               pstr.c_str()); 
      }

    jc.set(1,1,0);
    joint->setPosition(pos);
    //joint->setMsize(0.1);
    joint->setMsize(0.01);
    joint->setColor(jc);
    joint->display(true);
    joints.push_back(joint);
    }

  for (unsigned int i = 0; i < joints.size(); i++) {
    joint = joints[i];
    joint->setBodies (bodies[i], bodies[i+1]);
    }
  }

//*============================================================
//*==========          pm_CmdMultiBodyDna            ==========
//*============================================================
// process multibody command for dna.

void
pm_CmdMultiBodyDna (PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, dname, desc;
  vector<string> domains_desc; 
  PmMolecule *full_domain, *domain, *mol;
  vector<PmMolecule*> domains; 
  vector<string> begin_res, end_res;
  char chain1, chain2;

  int num_bodies;
  int nres1, nres2, step1, step2;
  int bid1, eid1, bid2, eid2; 
  int imin1, imax1, imin2, imax2;
  int jid1, jid2;
  vector<PmResidue*> bres1, bres2, eres1, eres2;

  PmResidue *res, *nres;
  char chain_id;
  stringstream dss;
  PmAtomFilter sidechain_filter, filter, disp_filter;
  vector<PmVector3> coords; 
  vector<PmVector3> domains_com; 
  PmMassProperties props;
  vector<PmAtom> atoms;
  vector<string> atom_names; 

  PmExtent extent;
  PmVector3 dom_color, joint_color, atom_color, bond_color;
  float width, msize;
  bool show_axes, show_atoms, show_bonds, no_bond_color; 

  string bname, jname, pstr;
  PmBody *body;
  PmJoint *joint;
  PmVector3 jpos, pos; 
  vector <PmBody*> bodies;
  vector <PmJoint*> joints;
  PmJointType joint_type;
  PmPcaResults pca;
  PmVector3 axis, axis1, axis2, axis3, u, v, w;
  vector<string> force_consts;
  float fconsts[3]; 
  bool verbose = pmSystem.getCmdVerbose();

  if (verbose) {
    fprintf (stderr, "\n------ create dna chain ------\n");
    }

  if (!dlist.getString("name", name)) {
    pm_ErrorReport (PM, "no name given.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  // get current molecule //

  pmSystem.getMolecule (&mol);

  if (!mol) {
    pm_ErrorReport (PM, "no current molecule.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  // get domain //

  if (dlist.getString("domain", dname)) {
    pmSystem.getDomain (dname, &full_domain);

    if (!full_domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
      pmSystem.setCmdError(true);
      return;
      }
    }
  else {
    pm_ErrorReport (PM, "no domain name given.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  // get start residues for the dna chain //

  if (dlist.getStringList("begin", begin_res)) {
    if (begin_res.size() != 2) {
      pm_ErrorReport (PM, "number of begin residues should be two. ", "*");
      pmSystem.setCmdError(true);
      return;
      }

    full_domain->getResidues(begin_res[0], bres1);
    full_domain->getResidues(begin_res[1], bres2);

    if (bres1.size() != 1) {
      pm_ErrorReport (PM, "can't process begin residue \"%s\" ", "*", 
                      begin_res[0].c_str());
      pmSystem.setCmdError(true);
      return;
      }

    if (bres2.size() != 1) {
      pm_ErrorReport (PM, "can't process begin residue \"%s\" ", "*", 
                      begin_res[1].c_str());
      pmSystem.setCmdError(true);
      return;
      }
    }

  // get end residues for the dna chain //

  if (dlist.getStringList("end", end_res)) {
    if (end_res.size() != 2) {
      pm_ErrorReport (PM, "number of end residues should be two. ", "*");
      pmSystem.setCmdError(true);
      return;
      }

    full_domain->getResidues(end_res[0], eres1);
    full_domain->getResidues(end_res[1], eres2);

    if (eres1.size() != 1) {
      pm_ErrorReport (PM, "can't process end residue \"%s\" ", "*", end_res[0].c_str());
      pmSystem.setCmdError(true);
      return;
      }

    if (eres2.size() != 1) {
      pm_ErrorReport (PM, "can't process end residue \"%s\" ", "*", end_res[1].c_str());
      pmSystem.setCmdError(true);
      return;
      }
    }

  // get force constants //

  if (dlist.getStringList("force_constants", force_consts)) {
    if (force_consts.size() != 3) {
      pm_ErrorReport (PM, "need three force constants ", "*"); 
      pmSystem.setCmdError(true);
      return;
      }

    for (unsigned int i = 0; i < force_consts.size(); i++) {
      if (!convToFloat(force_consts[i], fconsts[i])) {
        pm_ErrorReport (PM, "can't process force constant \"%s\" ", "*", 
                          force_consts[i].c_str());
        pmSystem.setCmdError(true);
        return;
        }
      }
    }
  else {
    pm_ErrorReport (PM, "no force constants set ", "*");
    pmSystem.setCmdError(true);
    return;
    }

  // process display arguments //

  if (!dlist.getVector("domain_color", dom_color)) {
    dom_color.set(1,1,1);
    }

  if (!dlist.getVector("joint_color", joint_color)) {
    joint_color.set(0.6, 0.6, 0.6);
    }

  if (!dlist.getVector("bond_color", bond_color)) {
    no_bond_color = true; 
    }
  else { 
    no_bond_color = false; 
    }

  if (!dlist.getFloat("width", width)) {
    width = 1;
    }

  if (!dlist.getFloat("msize", msize)) {
    msize = 0.05;
    }

  dlist.getBoolean("show_axes", show_axes); 
  dlist.getBoolean("show_atoms", show_atoms);
  dlist.getBoolean("show_bonds", show_bonds);

  if (dlist.getStringList("atom_names", atom_names)) {
    disp_filter.names = atom_names; 
    }

  // compute stepping increments //

  nres1 = eres1[0]->id - bres1[0]->id + 1;
  step1 = 9;

  if (nres1 < 0) {
    nres1 = -nres1;
    step1 = -9;
    } 

  nres2 = eres2[0]->id - bres2[0]->id + 1;
  step2 = 9;

  if (nres2 < 0) {
    nres2 = -nres2;
    step2 = -9;
    } 

  num_bodies = nres1 / 9 + 1;
  bid1 = bres1[0]->id;
  eid1 = bid1 + 9;
  chain1 = begin_res[0][0];

  bid2 = bres2[0]->id;
  eid2 = bid2 + 9;
  chain2 = begin_res[1][0];

  if (verbose) {
  fprintf (stderr, ">>> num segments = %d \n", num_bodies); 
  fprintf (stderr, ">>> step1 = %d \n", step1); 
  fprintf (stderr, ">>> step2 = %d \n", step2); 
  }

  // create domains and bodies //

  for (int i = 0; i < num_bodies; i++) {
    if (bid1 > eid1) {
      imin1 = eid1;
      imax1 = bid1;
      }
    else {
      imin1 = bid1;
      imax1 = eid1;
      }

    if (bid2 > eid2) {
      imin2 = eid2;
      imax2 = bid2;
      }
    else {
      imin2 = bid2;
      imax2 = eid2;
      }

    // create domain //

    dss << name << i+1;
    dname = dss.str();
    dss.str(std::string());

    dss << chain1 << '[' << imin1 << '-' << imax1 << ']';
    dss << chain2 << '[' << imin2 << '-' << imax2 << ']';
    desc = dss.str();
    dss.str(std::string());
    domains_desc.push_back(desc);

    if (verbose) {
    fprintf (stderr, ">>> domain %s %s --->", dname.c_str(), desc.c_str());
    }
    domain = mol->createDomain (dname, desc, sidechain_filter);

    if (!domain) { 
      pm_ErrorReport (PM, "can't create domain %s", "*", desc.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    domains.push_back(domain);
    domain->getMassProps(props);
    domains_com.push_back(props.com);

    domain->getExtent(extent);
    pmSystem.updateExtent (extent);
    domain->setColor(dom_color);
    domain->setLineWidth(width);
    domain->displayBackbone("", true);

    if (show_bonds) { 
      domain->setAtomColorType (PM_ATOM_COLOR_ELEMENT);

      if (no_bond_color) { 
        domain->setBondColorType (PM_ATOM_COLOR_ELEMENT);
        }
      else {
        domain->setBondColor(bond_color);
        }

      if (show_atoms) { 
        domain->setMarkerSize(0.05);
        domain->setDisplayBondAtoms (show_atoms);
        }

      domain->displayBonds("", show_atoms);
      }

    else if (show_atoms) { 
      domain->setAtomColorType (PM_ATOM_COLOR_ELEMENT);
      domain->displayAtoms("", disp_filter, show_atoms);
      }

    // create rigid body //

    dss << dname << "Body";
    bname = dss.str();
    dss.str(std::string());
    if (verbose) {
    fprintf (stderr, " body %s \n", bname.c_str());
    }

    body = new PmBody(bname, PM_BODY_RIGID);
    pmSystem.addBody (body);
    body->setPhysicalObj(domain);
    bodies.push_back(body);

    // create joint //

    if (i < num_bodies-1) {
      dss << name << "Jnt" << i+1;
      jname = dss.str();
      dss.str(std::string());
      joint = PmJoint::create(jname, PM_JOINT_BALL);
      pmSystem.addJoint (joint);

      if (step1 > 0) {
        jid1 = imax1; 
        }
      else {
        jid1 = imin1; 
        }

      if (step2 > 0) {
        jid2 = imax2; 
        }
      else {
        jid2 = imin2; 
        }

      dss << chain1 << '[' << jid1 << ']';
      dss << chain2 << '[' << jid2 << ']';
      pstr = dss.str();
      dss.str(std::string());
      if (verbose) {
      fprintf (stderr, "    joint at %s \n", pstr.c_str());
      }

      coords.clear();
      domain->getAtomCoords(pstr, filter, coords);
      jpos.set(0,0,0);

      for (unsigned int j = 0; j < coords.size(); j++) {
        jpos = jpos + coords[j];
        }

      jpos = (1.0/coords.size())*jpos;
      joint->setPosition(jpos);
      joints.push_back(joint);
      }
  
    bid1 += step1;
    eid1 += step1;
    bid2 += step2;
    eid2 += step2;
    }

  // add joint bodies and define spring axes. //

  filter.names.push_back("P");

  for (unsigned int i = 0; i < joints.size(); i++) {
    joint = joints[i];
    joint->setBodies (bodies[i], bodies[i+1]);
    joint->getPosition(jpos);
    PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint);

    // compute axis //

    u = domains_com[i+1] - jpos;
    u .normalize();

    domains[i]->getAtoms(atoms);
    domains[i]->getAtomCoords("", filter, coords);
    pos = coords.back();
    v = pos - jpos;
    v.normalize();

    w = u.cross(v);  
    w.normalize();
    v = u.cross(w);

    axis1 = u;
    axis2 = v;
    bjoint->setAxes(axis1, axis2);

    for (int j = 0; j < 3; j++) {
      bjoint->setForceConst(j+1, fconsts[j]);
      }

    // display joint //

    joint->setMsize(msize);
    joint->setColor(joint_color);
    joint->setShowAxes(show_axes);
    joint->display(true);
    }
  }

//*============================================================
//*==========          pm_CmdMultiBodyResidue        ==========
//*============================================================
// process multibody command for residue.

void
pm_CmdMultiBodyResidue(PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, dname, desc;
  PmMolecule *domain, *mol;

  int num_res;
  vector<PmResidue*> rlist;
  PmResidue *res, *nres;
  char chain_id;
  stringstream dss;
  PmAtomFilter sidechain_filter;

  PmExtent extent;
  PmVector3 dom_color, joint_color;
  float width, msize;
  bool show_bonds;

  string bname, jname, pstr;
  PmBody *body;
  PmJoint *joint;
  PmVector3 jpos, pos, axis;
  vector <PmBody*> bodies;
  vector <PmJoint*> joints;
  PmJointType joint_type;

  fprintf (stderr, "\n------ create residue chain ------ \n");

  if (!dlist.getString("name", name)) {
    pm_ErrorReport (PM, "no name given.", "*");
    return;
    }

  // get current molecule //

  pmSystem.getMolecule (&mol);

  if (!mol) {
    pm_ErrorReport (PM, "no current molecule.", "*");
    return;
    }

  // get domain //

  if (dlist.getString("domain", dname)) {
    pmSystem.getDomain (dname, &domain);

    if (!domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
      return;
      }
    }
  else {
    pm_ErrorReport (PM, "no domain name given.", "*");
    return;
    }

  if (!dlist.getVector("domain_color", dom_color)) {
    dom_color.set(1,1,1);
    }

  if (!dlist.getVector("joint_color", joint_color)) {
    joint_color.set(0.6, 0.6, 0.6);
    }

  if (!dlist.getFloat("width", width)) {
    width = 1;
    }

  if (!dlist.getFloat("msize", msize)) {
    msize = 0.05;
    }

  dlist.getBoolean("show_bonds",show_bonds);
    

  // create multibody model //

  domain->getResidues(desc, rlist);
  num_res = rlist.size();
  fprintf (stderr, "    >>> number of residues = %d \n", num_res);
  fprintf (stderr, "    >>> create objects: ");
  sidechain_filter.sidechain_ca = true;

  for (int i = 0; i < num_res; i++) {
    res = rlist[i];

    if (i < num_res-1) {
      nres = rlist[i+1];
      }
    else {
      nres = NULL; 
      }

    chain_id = res->getChainId();
    fprintf (stderr, "------ residue %d ----- \n", res->id);

    // create domain //

    dss << name << i+1;
    dname = dss.str();
    dss.str(std::string());

    if (nres) {
      dss << chain_id << '[' << res->id << '-' << nres->id << ']';
      }
    else {
      dss << chain_id << '[' << res->id << ']';
      }

    desc = dss.str();
    dss.str(std::string());

    fprintf (stderr, ">>> domain %s %s ---> ", dname.c_str(), desc.c_str());
    domain = mol->createDomain (dname, desc, sidechain_filter);
    vector<PmAtom> atoms;
    domain->getAtoms(atoms);

    /*
    for (int i = 0; i < atoms.size(); i++) {
      fprintf (stderr, ">>> atoms %d %s %d \n", i+1, atoms[i].name, atoms[i].seq); 
      }
    */

    if (domain) { 
      domain->getExtent(extent);
      pmSystem.updateExtent (extent);
      domain->setColor(dom_color);
      domain->setLineWidth(width);
      domain->displayBackbone("", true);

      if (show_bonds) {
        domain->setAtomColor(dom_color);
        domain->setDisplayBondAtoms(true);
        domain->displayBonds("", true);
        }
      }
    else {
      pm_ErrorReport (PM, "can't create domain %s", "*", desc.c_str());
      return;
      }

    // create rigid body //

    dss << dname << "Body";
    bname = dss.str();
    dss.str(std::string());
    fprintf (stderr, " body %s \n", bname.c_str());

    body = new PmBody(bname, PM_BODY_RIGID);
    pmSystem.addBody (body);
    body->setPhysicalObj(domain);
    bodies.push_back(body);

    // create joint //

    if (nres) {
      dss << name << "Jnt" << i+1;
      jname = dss.str();
      dss.str(std::string());
      joint = PmJoint::create(jname, PM_JOINT_BALL);
      pmSystem.addJoint (joint);

      dss << chain_id << '[' << nres->id << ']';
      pstr = dss.str();
      dss.str(std::string());

      domain->getResidueCoords(pstr, "CA", jpos);
      joint->setPosition(jpos);

      joint->setMsize(msize);
      joint->setColor(joint_color);
      joint->display(true);
      joints.push_back(joint);
      }
    }


  // add joint bodies //

  for (unsigned int i = 0; i < joints.size(); i++) {
    joint = joints[i];
    joint->setBodies (bodies[i], bodies[i+1]);
    }
  }

//*============================================================
//*==========          pm_CmdMultiBodyDomains        ==========
//*============================================================
// process multibody command for domains.

void
pm_CmdMultiBodyDomains(PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, dname, desc, jstr;
  stringstream dss;

  vector<string> domain_list;
  vector<string> joint_list;

  PmMolecule *mol, *domain;
  vector<PmMolecule*> domains;
  PmAtomFilter filter;

  string bname, jname, pstr;
  PmBody *body;
  PmJoint *joint;
  PmVector3 jpos, pos, axis, center, jaxis;
  vector <PmBody*> bodies;
  vector <PmJoint*> joints;
  PmExtent extent;
  PmVector3 dom_color, joint_color;
  PmVector3 axis1, axis2; 
  PmJointType joint_type;
  int num_joints;
  float width, msize;

  fprintf (stderr, "\n    >>> create multibody from domains. \n");

  pmSystem.getMolecule (&mol);

  if (!dlist.getString("name", name)) {
    pm_ErrorReport (PM, "no name given.", "*");
    return;
    }

  if (!dlist.getStringList("bodies", domain_list)) {
    pm_ErrorReport (PM, "no bodies given.", "*");
    return;
    }

  if (!dlist.getStringList("joints", joint_list)) {
    pm_ErrorReport (PM, "no joints given.", "*");
    return;
    }

  if (dlist.getString("joint_type", jstr)) {
    PmJoint::convJointType(jstr, joint_type);

    if (joint_type == PM_JOINT_UNKNOWN) {
      pm_ErrorReport (PM, "unknown joint type \"%s\".", "*", jstr.c_str());
      return;
      }

    if (joint_type == PM_JOINT_HINGE) {
      if (!dlist.getVector("axis", axis1)) {
        axis1.set(0,0,1);
        }
      }
    else if (joint_type == PM_JOINT_UNIVERSAL) {
      if (!dlist.getVector("axis1", axis1)) {
        axis1.set(0,0,1);
        }

      if (!dlist.getVector("axis2", axis2)) {
        axis2.set(1,0,0);
        }
      }
    }
  else {
    joint_type = PM_JOINT_BALL; 
    jstr = "ball";
    }

  if (!dlist.getVector("domain_color", dom_color)) {
    dom_color.set(1,1,1);
    }

  if (!dlist.getVector("joint_color", joint_color)) {
    joint_color.set(0.6, 0.6, 0.6);
    }

  if (!dlist.getFloat("width", width)) {
    width = 1;
    }

  if (!dlist.getFloat("msize", msize)) {
    msize = 0.05;
    }

  //===== create domains =====//

  fprintf (stderr, "\n------ create domains ------ \n");

  for (int i = 0; i < domain_list.size(); i++) {
    dss << name << i+1;
    dname = dss.str();
    dss.str(std::string());
    desc = domain_list[i];
    domain = mol->createDomain (dname, desc, filter);

    if (domain) {
      domain->getExtent(extent);
      pmSystem.updateExtent (extent);
      fprintf (stderr, ">>> create domain name = %s \n", dname.c_str());
      }
    else {
      pm_ErrorReport (PM, "domain %s not created.", "*", dname.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    domains.push_back(domain);
    domain->setColor(dom_color);
    domain->displayBackbone("", true);
    }

  //===== create rigid bodies =====//

  fprintf (stderr, "\n------ create rigid bodies ------ \n");

  for (int i = 0; i < domains.size(); i++) {
    domain = domains[i];
    domain->getName(dname);
    dss << dname << "Body";
    bname = dss.str();
    dss.str(std::string());
    fprintf (stderr, ">>> create rigid body name=%s \n", bname.c_str());

    body = new PmBody(bname, PM_BODY_RIGID);
    pmSystem.addBody (body);
    body->setPhysicalObj(domain);
    bodies.push_back(body);
    }

  //===== create joints =====//

  fprintf (stderr, "\n------ create joints ------ \n");
  num_joints = joint_list.size();

  for (int i = 0; i < domains.size()-1; i++) {
    domain = domains[i];
    domain->getName(dname);

    if (i == num_joints) {
      pm_ErrorWarnReport (PM, "not enough joints for bodies. ", "*");
      break;
      }

    pstr = joint_list[i];

    if (!domain->getResidueCoords(pstr, "CA", jpos)) {
      pm_ErrorReport (PM, "can't process joint residues %s for domain \"%s\".", "*", 
         pstr.c_str(), dname.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    dss << dname << "Jnt" << i+1;
    jname = dss.str();
    dss.str(std::string());

    joint = PmJoint::create(jname, joint_type);
    pmSystem.addJoint (joint);
    fprintf (stderr, ">>> create %s joint name = %s \n", jstr.c_str(), jname.c_str());

    joint->setPosition(jpos);
    joint->setColor(joint_color);
    joint->setMsize(msize);
    joint->display(true);
    joint->setBodies (bodies[i], bodies[i+1]);

    if (joint_type == PM_JOINT_HINGE) {
      PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);
      hjoint->setAxis(axis1);
      }
    else if (joint_type == PM_JOINT_UNIVERSAL) {
      PmUniversalJoint *ujoint = dynamic_cast<PmUniversalJoint*>(joint);
      ujoint->setAxis(1,axis1);
      ujoint->setAxis(2,axis2);
      }
    }
  }

//*============================================================
//*==========          pm_CmdMultiBodySolid          ==========
//*============================================================
// process multibody command for solids.

void
pm_CmdMultiBodySolid(PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, sname; 
  stringstream dss;

  PmSolid *solid;
  PmSolidType type;

  int num_objs;
  string bname, jname, pstr;
  PmBody *body;
  PmJoint *joint;
  PmVector3 jpos, pos, axis, center, jaxis;
  vector <PmBody*> bodies;
  vector <PmJoint*> joints;
  float radius, length;
  PmExtent extent;
  PmVector3 color; 
  PmJointType joint_type;

  fprintf (stderr, "\n    >>> create multibody from solids. \n");

  if (!dlist.getString("name", name)) {
    pm_ErrorReport (PM, "no name given.", "*");
    return;
    }

  if (!dlist.getInt("number", num_objs)) {
    pm_ErrorReport (PM, "number of solids not given.", "*");
    return;
    }

  dlist.getString ("solid", dv);
  PmSolid::convSolidType(dv, type);

  if (type == PM_SOLID_UNKNOWN) {
    pm_ErrorReport ("pm> ", "unknown solid type \"%s\".", "*", dv.c_str());
    return;
    }

  if (dlist.getString("joint_type", dv)) {
    PmJoint::convJointType(dv, joint_type);

    if (joint_type == PM_JOINT_UNKNOWN) {
      pm_ErrorReport (PM, "unknown joint type \"%s\".", "*", dv.c_str());
      return;
      }
    }
  else {
    joint_type = PM_JOINT_HINGE;
    }

  dlist.getFloat("radius", radius);
  dlist.getFloat("length", length);
  dlist.getVector("center", center);
  axis.set(0,1,0);
  fprintf (stderr, "    >>> radius = %f \n", radius);
  fprintf (stderr, "    >>> length = %f \n", length);
  fprintf (stderr, "    >>> center = %f %f %f \n", center[0], center[1], center[2]);

  for (int i = 0; i < num_objs; i++) {
    dss << name << i+1;
    sname = dss.str();
    dss.str(std::string());
    solid = PmSolid::create(sname, type);
    pmSystem.addSolid(solid);
    fprintf (stderr, "    >>> solid %s \n", sname.c_str());

    if (type == PM_SOLID_CYLINDER) {
      PmCylinder *cyl;
      PmSolid::getCylinder(sname, &cyl);
      cyl->setRadius(radius);
      cyl->setLength(length);
      cyl->setAxis(axis);
      cyl->setCenter(center);
      }

    /*
    solid->getExtent(extent);
    fprintf (stderr, "extent  min (%g %g %g)  max (%g %g %g) \n",
             extent.min[0], extent.min[1], extent.min[2],
             extent.max[0], extent.max[1], extent.max[2]);

    PmMassProperties props;
    solid->getMassProps(props);
    fprintf (stderr, "        com = %g %g %g \n", props.com[0], 
             props.com[1], props.com[2]);
    */

    pmSystem.updateExtent (extent);

    if (i == 0) {
      color.set(0,1,0);
      }
    else if (i == 1) {
      color.set(0,1,1);
      }
    else {
      color.set(1,1,1);
      }

    solid->setColor(color);
    solid->display(true);

    // create bodies //

    dss << name << "Body" << i+1;
    bname = dss.str();
    dss.str(std::string());
    fprintf (stderr, "    >>> solid body %s \n", bname.c_str());

    body = new PmBody(bname, PM_BODY_RIGID);
    pmSystem.addBody (body);
    body->setPhysicalObj(solid);
    bodies.push_back(body);

    if (i == num_objs -1) {
      break;
      }

    // create joints //
    
    dss << name << "Jnt" << i+1;
    jname = dss.str();
    dss.str(std::string());
    joint = PmJoint::create(jname, joint_type);
    jpos = center + (0.5)*length*axis;
    joint->setPosition (jpos);

    if (joint_type == PM_JOINT_HINGE) {
      PmHingeJoint *pjoint = dynamic_cast<PmHingeJoint*>(joint);
      jaxis.set(0,0,1);
      pjoint->setAxis(jaxis);
      }

    color.set(1,0,0);
    joint->setColor(color);
    joint->setMsize(radius/4.0);
    joint->display(true);
    joints.push_back(joint);
    pmSystem.addJoint (joint);

    if (type == PM_SOLID_CYLINDER) {
      center[1] = center[1] + length;
      }
    }

  // add bodies to joints //

  for (unsigned int i = 0; i < joints.size(); i++) {
    joint = joints[i];
    joint->setBodies (bodies[i], bodies[i+1]);
    }
  }

//*============================================================*
//*==========      pm_CmdMultiBodySecondaryStruct    ==========*
//*============================================================*
// process multibody command for secondary structure. 

void
pm_CmdMultiBodySecondaryStruct(PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, dname, desc, sname, sdesc;
  PmMolecule *domain, *mol;
  int num_chains;
  vector<PmChain*> clist;
  PmTopology *topo;
  PmStructureList *slist, *prev_slist;

  int num_res;
  vector<PmResidue*> rlist;
  PmResidue *first_res, *last_res;
  char chain_id;
  stringstream dss;
  PmMolecule *sdomain; 
  PmAtomFilter filter; 
  int init_snum, term_snum;

  PmExtent extent;
  PmVector3 color, joint_color;
  float width, msize;

  int num_objs;
  string bname, jname, pstr;
  PmBody *body;
  PmJoint *joint;
  PmVector3 jpos, pos, axis;
  vector <PmMolecule*> domains;
  vector <PmBody*> bodies;
  vector <PmJoint*> joints;

  fprintf (stderr, "\n    >>> create multibody from secondary structure. \n");

  if (!dlist.getString("name", name)) {
    pm_ErrorReport (PM, "no name given.", "*");
    return;
    }

  // get current molecule //

  pmSystem.getMolecule (&mol);

  if (!mol) {
    pm_ErrorReport (PM, "no current molecule.", "*");
    return;
    }

  // get topology //

  mol->getTopology(&topo);

  if (!topo) {
    pm_ErrorReport (PM, "no topology for molecule.", "*");
    return;
    }

  // get domain //

  if (dlist.getString("domain", dname)) {
    pmSystem.getDomain (dname, &domain);

    if (!domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
      return;
      }
    }
  else {
    pm_ErrorReport (PM, "no domain name given.", "*");
    return;
    }

  domain->getResidues(desc, rlist);
  num_res = rlist.size();
  first_res = rlist[0];
  last_res = rlist[num_res-1];
  fprintf (stderr, "    >>> number of residues = %d \n", num_res);

  // check for missing residues. this could cause a problem //
  // when connecting rigid bodies.                          //

  if ((last_res->id - first_res->id + 1) != num_res) {
    pm_ErrorReport (PM, "domain has missing residues.", "*");
    return;
    }

  chain_id = first_res->getChainId();
  fprintf (stderr, "    >>> use chain %c \n", chain_id);
  fprintf (stderr, "    ---- create objects ---   \n");

  mol->getChains (clist);
  num_chains = clist.size();
  slist = NULL;

  // find topology for the domain chain  //

  for (int i = 0; i < num_chains; i++) {
    if (topo[i].chain_id == chain_id) {
      slist = topo[i].list; 
      break;
      }
    }

  // get display parameters //

  if (!dlist.getFloat("width", width)) {
    width = 1;
    }

  if (!dlist.getFloat("msize", msize)) {
    msize = 0.02;
    }

  if (!dlist.getVector("joint_color", joint_color)) { 
    joint_color.set(1,1,1);
    }

  // create a domain and rigid body for each secondary   //
  // structure in the topology that is within the domain. // 
  // add joints between pairs of bodies.

  num_objs = 0;
  prev_slist = NULL;

  while (slist) {
    init_snum = slist->init_seq_num;
    term_snum = slist->term_seq_num;

    if ((init_snum < first_res->id) || (term_snum > last_res->id)) { 
      slist = slist->next;
      continue;
      }

    if (slist->type == PM_SECONDARY_STRUCTURE_HELIX) {
      dss << "helix" << slist->id;
      color.set(1,0,0);
      }
    else if (slist->type == PM_SECONDARY_STRUCTURE_SHEET) {
      dss << "sheet" << slist->id;
      color.set(0,1,0);
      }

    // for a loop add one residue to the end to  //
    // ensure the domains overlap by one residue //

    else if (slist->type == PM_SECONDARY_STRUCTURE_LOOP) {
      dss << "loop" << slist->id;
      color.set(0,1,1);
      if (init_snum != first_res->id) init_snum -= 1;
      if (term_snum != last_res->id)  term_snum += 1;
      }

    sname = dss.str();
    dss.str(std::string());
    
    dss << name << '_' << sname;
    dname = dss.str();
    dss.str(std::string());

    dss << chain_id << '[' << init_snum << '-' << term_snum << ']';
    sdesc = dss.str();
    dss.str(std::string());

    //sdomain = mol->createDomain (dname, sname, filter);
    sdomain = mol->createDomain (dname, sdesc, filter);
    domains.push_back(sdomain);
    fprintf (stderr, "    >>> domain %s = %c[%d-%d] \n", dname.c_str(),
             chain_id, init_snum, term_snum);

    if (!sdomain) {
      return;
      }

    // display the structure //

    sdomain->getExtent(extent);
    pmSystem.updateExtent (extent);
    sdomain->setColor(color);
    sdomain->setLineWidth(width);
    sdomain->displayBackbone("", true);

    //sdomain->setAtomColor(color);
    //sdomain->setDisplayBondAtoms(true);
    //sdomain->displayBonds("", true);
    sdomain->displayBackbone("", true);

    // create a rigid body //

    dss << dname << "Body";
    bname = dss.str();
    dss.str(std::string());
    fprintf (stderr, "    >>> rigid body %s \n", bname.c_str());

    body = new PmBody(bname, PM_BODY_RIGID);
    pmSystem.addBody (body);
    body->setPhysicalObj(sdomain);
    bodies.push_back(body);

    // create a joint //

    if (prev_slist) {
      //dss << dname << "Jnt";
      dss << name << "Jnt" << num_objs+1;
      jname = dss.str();
      dss.str(std::string());

      joint = PmJoint::create(jname, PM_JOINT_BALL);
      //PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint);
      //joint = PmJoint::create(jname, PM_JOINT_HINGE);
      //PmHingeJoint *pjoint = dynamic_cast<PmHingeJoint*>(joint);
      pmSystem.addJoint (joint);

      // for a loop add one residue to the end to  //
      // ensure the domains overlap by one residue //

      term_snum = prev_slist->term_seq_num;

      if (prev_slist->type == PM_SECONDARY_STRUCTURE_LOOP) {
        if (term_snum != last_res->id)  term_snum += 1;
        }

      dss << chain_id << '[' << term_snum << ']';
      pstr = dss.str();
      dss.str(std::string());
      domains[num_objs-1]->getResidueCoords(pstr, "CA", jpos);

      joint->setPosition(jpos);
      //pjoint->setAxis(axis);
      joint->setMsize(msize);
      joint->setColor(joint_color);
      joint->display(true);
      joint->setBodies (bodies[num_objs-1], bodies[num_objs]);
      joints.push_back(joint);

      fprintf (stderr, "    >>> joint %s at %s  pos = [ %f, %f, %f ] \n", jname.c_str(),
               pstr.c_str(), jpos[0], jpos[1], jpos[2]);
      }

    num_objs += 1;
    fprintf (stderr, "    ------------------------- \n");
    prev_slist = slist;
    slist = slist->next;
    }
  }

//*============================================================*
//*==========          pm_CmdMultiBodyKchain         ==========*
//*============================================================*
// process multibody command for a kinematic chain.

void
pm_CmdMultiBodyKchain(PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, dname, desc;
  PmMolecule *domain, *mol;

  int num_res;
  vector<PmResidue*> rlist;
  PmResidue *res, *nres;
  char chain_id;
  stringstream dss;
  PmMolecule *peptide_domain, *sidechain_domain;
  PmAtomFilter peptide_filter, sidechain_filter;
  bool last_peptide = false;

  PmExtent extent;
  PmVector3 color;

  string bname, jname, pstr, nstr; 
  PmBody *body;
  PmJoint *joint;
  PmVector3 jpos, pos, axis;
  vector <PmBody*> bodies;
  vector <PmJoint*> joints;
  PmJointType joint_type;
  bool verbose = pmSystem.getCmdVerbose();
  bool first_peptide=true;

  if (verbose) {
    fprintf (stderr, "\n------ create kinematic chain ------ \n");
    }

  if (!dlist.getString("name", name)) {  
    pm_ErrorReport (PM, "no name given.", "*");
    return;
    }

  //===== get current molecule =====//

  pmSystem.getMolecule (&mol);

  if (!mol) {
    pm_ErrorReport (PM, "no current molecule.", "*");
    return;
    }

  //===== get domain =====//

  if (dlist.getString("domain", dname)) {  
    pmSystem.getDomain (dname, &domain);

    if (!domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
      return;
      }
    }
  else {
    pm_ErrorReport (PM, "no domain name given.", "*");
    return;
    }

  if (dlist.getString("joint_type", dv)) {  
    PmJoint::convJointType(dv, joint_type);

    if (joint_type == PM_JOINT_UNKNOWN) {
      pm_ErrorReport (PM, "unknown joint type \"%s\".", "*", dv.c_str());
      return;
      }
    }
  else {
    joint_type = PM_JOINT_HINGE;
    }

  //===== flag to make the last domain a peptide =====//

  if (dlist.getString("last", dv)) {  
    last_peptide = (dv == "peptide");
    }

  domain->getResidues(desc, rlist);
  num_res = rlist.size();

  if (verbose) {
    fprintf (stderr, ">>> number of residues=%d \n", num_res);
    }
  /*
  */
  peptide_filter.peptide = true;
  sidechain_filter.sidechain_group = true;

  //===== create domains and rigid bodies =====//

  vector<PmMolecule*> peptide_domains, sidechain_domains; 
  if (verbose) {
    fprintf (stderr, "\n------ create domains and bodies ------ \n");
    }

  for (int i = 0; i < num_res-1; i++) {
    res = rlist[i];
    nres = rlist[i+1];
    chain_id = res->getChainId();
    //fprintf (stderr, "--- residue %d ---\n", res->id);

    //----- create sidechain domain and body -----//

    if (first_peptide || (i != 0)) {
      dss << name << "Sc" << res->id;
      //dss << name << "Sc" << i+1;
      dname = dss.str();
      dss.str(std::string());
      dss << chain_id << '[' << res->id << ']';
      desc = dss.str();
      dss.str(std::string());

      if (verbose) {
        fprintf (stderr, ">>> sidechain domain=%s residues=%s",dname.c_str(),desc.c_str());
        }

      sidechain_domain = mol->createDomain (dname, desc, sidechain_filter);

      if (sidechain_domain) {
        sidechain_domain->getExtent(extent);
        pmSystem.updateExtent (extent);
        //sidechain_domain->setDisplayBondAtoms(true);
        color.set(0,1,0);
        sidechain_domain->setAtomColor(color);
        sidechain_domain->displayBonds("", true);
        vector<PmAtom> satoms;
        sidechain_domain->getAtoms(satoms);

        if (verbose) {
          fprintf (stderr, " atoms={");
          for (int ia = 0; ia < satoms.size(); ia++) {
            fprintf (stderr, "%s ", satoms[ia].name); 
            }
          fprintf (stderr, "} ");
          }
        }
      else {
        string mol_name;
        mol->getName(mol_name);
        pm_ErrorReport (PM, "can't create domain with desc=%s from mol=%s ", "*", 
                        desc.c_str(), mol_name.c_str());
        return;
        }

      sidechain_domains.push_back(sidechain_domain);
      dss << name << "Sc" << res->id << "Body";
      //dss << name << "ScBody" << i+1;
      bname = dss.str();
      dss.str(std::string());

      if (verbose) {
        fprintf (stderr, " body=%s \n", bname.c_str());
        }

      body = new PmBody(bname, PM_BODY_RIGID);
      pmSystem.addBody (body);
      body->setPhysicalObj(sidechain_domain);
      bodies.push_back(body);
      }

    //----- create peptide domain and rigid body -----//

    dss << name << "Pt" << res->id;
    //dss << name << "Pt" << i+1;
    dname = dss.str();
    dss.str(std::string());
    dss << chain_id << '[' << res->id << '-' << nres->id << ']';
    desc = dss.str();
    dss.str(std::string());

    if (verbose) {
      fprintf (stderr, ">>> peptide domain=%s residues=%s ", dname.c_str(), desc.c_str());
      }

    peptide_domain = mol->createDomain (dname, desc, peptide_filter);
    vector<PmAtom> pep_atoms;
    peptide_domain->getAtoms(pep_atoms);

    if (verbose) {
      fprintf (stderr, "atoms={");
      for (int ia = 0; ia < pep_atoms.size(); ia++) {
        fprintf (stderr, "%s ", pep_atoms[ia].name); 
        }
      fprintf (stderr, "} ");
      }

    if (peptide_domain) {
      peptide_domain->getExtent(extent);
      pmSystem.updateExtent (extent);
      //peptide_domain->setColor(color);
      //peptide_domain->setLineWidth(width);
      //peptide_domain->displayBackbone("", true);

      color.set(1,0,0);
      peptide_domain->setAtomColor(color);
      peptide_domain->setDisplayBondAtoms(true);
      peptide_domain->displayBonds("", true);
      }
    else {
      return;
      }

    peptide_domains.push_back(peptide_domain);
    dss << name << "Pt" << res->id << "Body";
    //dss << name << "PtBody" << i+1;
    bname = dss.str();
    dss.str(std::string());
    if (verbose) {
    fprintf (stderr, " body=%s\n", bname.c_str());
    }

    body = new PmBody(bname, PM_BODY_RIGID);
    pmSystem.addBody (body);
    body->setPhysicalObj(peptide_domain);
    bodies.push_back(body);
    }

  //===== create joints =====//

  PmMolecule *next_sidechain_domain;
  PmVector3 jc;
  int res_count=0;

  #ifdef pm_CmdMultiBodyKchain_pr_jnt
  fprintf (stderr, "\n------ create joints ------ \n");
  #endif

  for (int i = 0; i < num_res-1; i++) {
    res = rlist[i];
    nres = rlist[i+1];
    chain_id = res->getChainId();
    #ifdef pm_CmdMultiBodyKchain_pr_jnt
    fprintf (stderr, "--- residue %d ---\n", res->id);
    #endif
    sidechain_domain = sidechain_domains[i];
    peptide_domain = peptide_domains[i];

    //----- CA-C joint -----//

    dss << name << "CACJnt" << i+1;
    jname = dss.str();
    dss.str(std::string());
    joint = PmJoint::create(jname, joint_type);
    pmSystem.addJoint (joint);

    dss << chain_id << '[' << res->id << ']';
    pstr = dss.str();
    dss.str(std::string());
    dss << chain_id << '[' << nres->id << ']';
    nstr = dss.str();
    dss.str(std::string());

    if (!sidechain_domain->getResidueCoords(pstr, "CA", pos)) { 
      fprintf (stderr, "*** no CA for residue %s \n", pstr.c_str()); 
      }

    if (!peptide_domain->getResidueCoords(pstr, "C", jpos)) {
      fprintf (stderr, "*** no C for residue %s \n", pstr.c_str()); 
      }

    if (joint_type == PM_JOINT_HINGE) {
      PmHingeJoint *pjoint = dynamic_cast<PmHingeJoint*>(joint);
      axis = jpos - pos;
      axis.normalize(); 
      pjoint->setAxis(axis);
      pjoint->setPosition2(pos);
      }
    else if (joint_type == PM_JOINT_BALL) {
      }

    #ifdef pm_CmdMultiBodyKchain_pr_jnt
    fprintf (stderr, ">>> joint=%s pos=[%f, %f, %f] \n", jname.c_str(), 
             jpos[0], jpos[1], jpos[2]); 
    #endif

    jc.set(1,1,0);
    joint->setPosition(jpos);
    //joint->setMsize(0.1);
    joint->setMsize(0.05);
    joint->setColor(jc);
    joint->display(true);
    joints.push_back(joint);

    //----- N-CA joint -----//

    if (i == num_res-2) {
      continue;
      }

    //res = nres;
    //nres = rlist[i+2];
    dss << name << "NCAJnt" << i+1;
    jname = dss.str();
    dss.str(std::string());
    joint = PmJoint::create(jname, joint_type);
    pmSystem.addJoint (joint);

    dss << chain_id << '[' << res->id << ']';
    pstr = dss.str();
    dss.str(std::string());
    dss << chain_id << '[' << nres->id << ']';
    nstr = dss.str();
    dss.str(std::string());

    next_sidechain_domain = sidechain_domains[i+1];

    if (!peptide_domain->getResidueCoords(nstr, "N", pos)) {
      fprintf (stderr, "*** no N for residue %s \n", nstr.c_str()); 
      }

    if (!next_sidechain_domain->getResidueCoords(nstr, "CA", jpos)) { 
      fprintf (stderr, "*** no CA for N-CA joint for residue %s \n", nstr.c_str()); 
      }

    joint->setPosition(jpos);

    if (joint_type == PM_JOINT_HINGE) {
      PmHingeJoint *pjoint = dynamic_cast<PmHingeJoint*>(joint);
      axis = jpos - pos;
      axis.normalize(); 
      pjoint->setAxis(axis);
      pjoint->setPosition2(pos);
      }

    jc.set(1,0,1);
    joint->setColor(jc);
    //joint->setMsize(0.05);
    joint->setMsize(0.1);
    joint->display(true);
    joints.push_back(joint);

    #ifdef pm_CmdMultiBodyKchain_pr_jnt
    fprintf (stderr, ">>> joint=%s residue=%s pos=[%f, %f, %f] \n", jname.c_str(), 
             pstr.c_str(), jpos[0], jpos[1], jpos[2]); 
    #endif
    }

  for (unsigned int i = 0; i < joints.size(); i++) {
    joint = joints[i];
    joint->setBodies (bodies[i], bodies[i+1]);
    }
  }

//*============================================================*
//*==========             pm_CmdMultiBody            ==========*
//*============================================================*
// process multibody command.

void
pm_CmdMultiBody (PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, bname, dname, res, type_str, use_type;
  vector<string> body_names;
  vector<string> domain_names;
  PmVector3 bpos;

  // get the next data item //

  dlist.getNext (data);

  // create a new body and add it to the system //

  if (data.name == "create") {
    dlist.getString("name", name);

    if (!dlist.getString("type", type_str)) {
      pm_ErrorReport (PM, "no multibody type specified.", "*");
      return;
      }

    if (type_str == "residue") {
      pm_CmdMultiBodyResidue(dlist);
      }

    else if (type_str == "dna") {
      pm_CmdMultiBodyDna(dlist);
      }

    else if (type_str == "domains") {
      pm_CmdMultiBodyDomains(dlist);
      }

    else if (type_str == "kinematic_chain") {
      if (dlist.getString("solid", dv)) {
        pm_CmdMultiBodySolid(dlist);
        }
      else {
        pm_CmdMultiBodyKchain(dlist);
        }
      }

    else if (type_str == "secondary_structure") {
      pm_CmdMultiBodySecondaryStruct(dlist);
      }

    else if (type_str == "c_alpha_chain") {
      pm_CmdMultiBodyCalpha(dlist);
      }

    else {
      pm_ErrorReport (PM, "unknown multibody type \"%s\" ", "*", type_str.c_str());
      pmSystem.setCmdError(true);
      return;
      }
    }
  }

//*============================================================*
//*==========          pm_CmdMultiBodyDnaPChain      ==========*
//*============================================================*
// process multibody command for a dna p-chain.

void
pm_CmdMultiBodyDnaPChain(PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, dname, desc;
  PmMolecule *domain, *mol;

  int num_res;
  vector<PmResidue*> rlist;
  PmResidue *res, *nres;
  char chain_id;
  stringstream dss;
  PmMolecule *p_domain; 
  PmAtomFilter sidechain_filter;
  PmAtom *atom;

  PmExtent extent;
  PmVector3 color;

  string bname, jname, pstr, nstr; 
  PmBody *body;
  PmJoint *joint;
  PmVector3 jpos, pos, axis;
  vector <PmBody*> bodies;
  vector <PmJoint*> joints;
  PmJointType joint_type;
  bool verbose = pmSystem.getCmdVerbose();
  bool first_peptide=true;

  if (verbose) {
    fprintf (stderr, "\n------ create dna p-chain ------ \n");
    }

  if (!dlist.getString("name", name)) {  
    pm_ErrorReport (PM, "no name given.", "*");
    return;
    }

  //===== get current molecule =====//

  pmSystem.getMolecule (&mol);

  if (!mol) {
    pm_ErrorReport (PM, "no current molecule.", "*");
    return;
    }

  //===== get domain =====//

  if (dlist.getString("domain", dname)) {  
    pmSystem.getDomain (dname, &domain);

    if (!domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
      return;
      }
    }
  else {
    pm_ErrorReport (PM, "no domain name given.", "*");
    return;
    }

  joint_type = PM_JOINT_BALL;


  //===== get domain residues =====//

  domain->getResidues(desc, rlist);
  num_res = rlist.size();

  if (verbose) {
    fprintf (stderr, ">>> number of residues=%d \n", num_res);
    }

  sidechain_filter.sidechain_group = true;

  //===== create domains and rigid bodies =====//

  vector<PmMolecule*> p_domains; 

  if (verbose) {
    fprintf (stderr, "\n------ create domains and bodies ------ \n");
    }

  for (int i = 0; i < num_res-1; i++) {
    res = rlist[i];
    nres = rlist[i+1];
    chain_id = res->getChainId();

    dss << name << "P" << res->id;
    dname = dss.str();
    dss.str(std::string());
    dss << chain_id << '[' << res->id << ']';
    desc = dss.str();
    dss.str(std::string());

    if (verbose) {
      fprintf (stderr, ">>> P domain=%s residues=%s \n",dname.c_str(),desc.c_str());
      }

    p_domain = mol->createDomain(dname, desc, sidechain_filter);

    if (!p_domain) {
      string mol_name;
      mol->getName(mol_name);
      pm_ErrorReport (PM, "can't create domain with desc=%s from mol=%s ", "*", 
                      desc.c_str(), mol_name.c_str());
      return;
      }
    else {
      fprintf (stderr, ">>> P domain found \n");
      }

    // add peptide and next residue N and CA atoms //

    res->getAtom("P", &atom);
    if (!atom) { 
      fprintf(stderr, ">>> no P atom in res\n");
    }

    p_domain->addAtom(atom);

    /*
    res->getAtom ("C", &atom);
    ca_domain->addAtom (atom);
    nres->getAtom ("N", &atom);
    ca_domain->addAtom (atom);
    nres->getAtom ("CA", &atom);
    ca_domain->addAtom (atom);
    */

    p_domain->getExtent(extent);
    pmSystem.updateExtent(extent);

    if (i % 2) {
      color.set(0,1,0);
      }
    else {
      color.set(1,0,0);
      }

    p_domain->setAtomColor(color);
    p_domain->setBondAtomRenderType(PM_MOLECULE_RENDER_LINE);
    p_domain->setDisplayBondAtoms (true);
    p_domain->displayBonds("", true);

    vector<PmAtom> satoms;
    p_domain->getAtoms(satoms);

    if (verbose) {
      fprintf (stderr, " atoms={");
      for (int ia = 0; ia < satoms.size(); ia++) {
        fprintf (stderr, "%d:%s ", satoms[ia].seq, satoms[ia].name); 
        }
      fprintf (stderr, "} ");
      }

    p_domains.push_back(p_domain);
    dss << name << "P" << res->id << "Body";
    bname = dss.str();
    dss.str(std::string());

    if (verbose) {
      fprintf (stderr, " body=%s \n", bname.c_str());
      }

    body = new PmBody(bname, PM_BODY_RIGID);
    pmSystem.addBody (body);
    body->setPhysicalObj(p_domain);
    bodies.push_back(body);
    }

  //===== create joints =====//

  PmMolecule *next_p_domain;
  PmVector3 jc;
  int res_count=0;

  if (verbose) {
    fprintf (stderr, "\n------ create joints ------ \n");
    }

  for (int i = 0; i < num_res-2; i++) {
    res = rlist[i];
    nres = rlist[i+1];
    chain_id = res->getChainId();

    p_domain = p_domains[i];

    dss << name << "PJnt" << i+1;
    jname = dss.str();
    dss.str(std::string());
    joint = PmJoint::create(jname, PM_JOINT_BALL);
    pmSystem.addJoint (joint);

    dss << chain_id << '[' << nres->id << ']';
    pstr = dss.str();
    dss.str(std::string());

    if (!p_domain->getResidueCoords(pstr, "P", pos)) { 
      fprintf (stderr, "*** no P for residue %s \n", pstr.c_str()); 
      }

    if (verbose) {
      fprintf (stderr, ">>> create joint=%s position=%s \n", jname.c_str(),
               pstr.c_str()); 
      }

    jc.set(1,1,0);
    joint->setPosition(pos);
    //joint->setMsize(0.1);
    joint->setMsize(0.01);
    joint->setColor(jc);
    joint->display(true);
    joints.push_back(joint);
    }

  for (unsigned int i = 0; i < joints.size(); i++) {
    joint = joints[i];
    joint->setBodies (bodies[i], bodies[i+1]);
    }
  }
