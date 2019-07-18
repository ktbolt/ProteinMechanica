
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
//              d o m a i n   c o m m a n d s                //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========               pm_CmdDomainHelix        ==========*
//*============================================================*
// process domain helix command.

void
pm_CmdDomainHelix(PmMolecule *domain, PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name;
  bool show; 
  vector<PmHelixProps> helix_props;
  PmGraphicsAttributes atts;
  PmVector3 color;
  float scale;

  atts.setScale(0.01); 
  //atts.setMarker(false); 

  domain->getHelixProps(helix_props);

  while (dlist.getNext(data,dv)) {
    if (data.name == "print") {
      }

    else if (data.name == "color") {
      data.getVector(color);
      atts.setColor(color);
      }

    else if (data.name == "msize") {
      scale = data.getFloat();
      atts.setScale(scale); 
      }

    else if (data.name == "show") {
      atts.visible = data.getBoolean();
      domain->displayHelixProps(helix_props, atts);
      }
    }
  }

//*============================================================*
//*==========           pm_CmdDomainCheck            ==========*
//*============================================================*
// process domain check command.

void
pm_CmdDomainCheck(PmMolecule *domain, PmCmdDataList& dlist)
  {
  domain->check();
  }

//*============================================================*
//*==========           pm_CmdDomainHbonds           ==========*
//*============================================================*
// process domain hbond command.

void
pm_CmdDomainHbonds (PmMolecule *domain, PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name;
  PmVector3 color;
  bool flag, show, show_set;
  float width;

  show_set = true;
  show = true;

  if (dlist.getBoolean("show", flag)) {
    show_set = true;
    show = flag; 
    }

  if (!dlist.getVector("color", color)) {
    color.set(1,0,1);
    }

  if (!dlist.getFloat("width", width)) {
    width = 1.0;     
    }

  domain->displayHydrogenBonds(color, width, show);
  }

//*============================================================*
//*==========         pm_CmdDomainPrint              ==========*
//*============================================================*
// process domain print command.

void
pm_CmdDomainPrint (PmMolecule *domain, PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, str, name;

  domain->getName(name);

  if (dlist.getString("properties", str)) {
    if (str == "mass") {
      PmMassProperties mass_props;
      domain->getMassProps(mass_props);
      fprintf (stderr, "\n------ domain \"%s\" mass properties ------ \n", name.c_str());
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

    else if (str == "geometric") {
      PmPcaResults pca;
      PmAtomFilter filter;

      domain->compPrincipalComponents(filter, pca);
      fprintf (stderr, "\n------ domain \"%s\" geometric properties ------ \n", 
               name.c_str());
      fprintf (stderr, "center = %g %g %g \n", pca.com[0],pca.com[1],pca.com[2]);
      fprintf (stderr, "axis1 = %g %g %g \n", pca.axis1[0],pca.axis1[1],pca.axis1[2]);
      fprintf (stderr, "axis2 = %g %g %g \n", pca.axis2[0],pca.axis2[1],pca.axis2[2]);
      fprintf (stderr, "axis3 = %g %g %g \n", pca.axis3[0],pca.axis3[1],pca.axis3[2]);
      fprintf (stderr, "scales = %g %g %g \n", pca.s1, pca.s2, pca.s3);
      }
    }

  else { 
    domain->print();
    }
  }

//*============================================================*
//*==========         pm_CmdDomainContact            ==========*
//*============================================================*
// process domain contact command.

void
pm_CmdDomainContact(PmMolecule *domain, PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv;
  vector<string> list;
  PmMolecule *mol;
  bool use_radii = false; 
  bool show = true; 
  PmVector3 color(1,1,1);
  PmAtomFilter filter; 
  float cutoff = 1.0;
  float tol = 0.01;

  // get current molecule //

  pmSystem.getMolecule (&mol);

  if (!mol) {
    pm_ErrorReport (PM, "no current molecule.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  while (dlist.getNext(data,dv)) {
   if (data.name == "atom_names") {
     dlist.getStringList("atom_names", filter.names);
     }
   else if (data.name == "color") {
     data.getVector(color);
     }
   else if (data.name == "tolerance") {
     tol = data.getFloat();
     }
   else if (data.name == "show") {
     show = data.getBoolean();
     }
   else if (data.name == "use_radii") { 
     use_radii = data.getBoolean();
     }
   else if (data.name == "use_sidechains") { 
     domain->getBackboneAtomNames(filter.names);
     filter.exclude = true;
     }
   }

  if (dlist.getStringList("list", list)) {
    vector <PmMolecule*> dom_list;
    PmMolecule *dom;

    for (unsigned int i = 0; i < list.size(); i++) {
      pmSystem.getDomain(list[i], &dom);

      if (!dom) {
        pm_ErrorReport (PM, "no domain named \"%s\".", "*", list[i].c_str());
        pmSystem.setCmdError(true);
        return;
        }

      dom_list.push_back(dom);
      }

    domain->checkDomainContact(dom_list, tol, color, show);
    }
  else {
    domain->displayResContact(cutoff, use_radii, filter, color, show);
    }
  }

//*============================================================*
//*==========         pm_CmdDomainTopology           ==========*
//*============================================================*
// process domain topology command.

void
pm_CmdDomainTopology(PmMolecule *domain, PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv; 
  PmMolecule *mol;
  bool chain = false, contact = false;
  bool show = true; 
  PmVector3 color(1,1,1);

  // get current molecule //

  pmSystem.getMolecule (&mol);

  if (!mol) {
    pm_ErrorReport (PM, "no current molecule.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  while (dlist.getNext(data,dv)) {
    if (data.name == "chain") {
      chain = data.getBoolean();
      }
    else if (data.name == "contact") {
      contact = data.getBoolean();
      }
    else if (data.name == "color") {
      data.getVector(color);
      }
    else if (data.name == "show") {
      show = data.getBoolean();
      }
    }

  domain->displayTopology(chain, contact, color, show);
  }

//*============================================================*
//*==========         pm_CmdDomainJoin               ==========*
//*============================================================*
// process domain joint command.

void
pm_CmdDomainJoin (PmMolecule *domain, PmCmdDataList& dlist)
  {
  string dv, dname, name, desc, dst_desc;
  vector<string> chains;
  PmAtomFilter filter, dst_filter;
  bool nterm;
  PmMolecule *dst_domain, *join_domain;

  if (!dlist.getString("destination", dname)) {
    pm_ErrorReport (PM, "no destination domain given.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  if (!dlist.getString("name", name)) {
    pm_ErrorReport (PM, "no join domain given.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  pmSystem.getMolecule (name, &dst_domain);

  if (dst_domain) {
    pm_ErrorReport (PM, "domain named already exists \"%s\".", "*", name.c_str());
    pmSystem.setCmdError(true);
    return;
    }

  if (!dlist.getBoolean("nterm", nterm)) {
    nterm = false;
    }

  pmSystem.getDomain(dname, &dst_domain);

  if (!dst_domain) {
    pmSystem.getMolecule (dname, &dst_domain);

    if (!dst_domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
      pmSystem.setCmdError(true);
      return;
      }
    }

  // check for valid joining residues //

  if (dlist.getString("residue", desc)) {
    if (!domain->checkDescriptor(desc)) {
      pm_ErrorReport (PM, "residue descriptor \"%s\" is not valid for this domain.",
                      "*", desc.c_str());
      pmSystem.setCmdError(true);
      return;
      }
    }

  if (dlist.getString("destination_residue", dst_desc)) {
    if (!dst_domain->checkDescriptor(dst_desc)) {
      pm_ErrorReport (PM, "residue descriptor \"%s\" is not valid for domain \"%s\".",
                      "*", dst_desc.c_str(), dname.c_str());
      pmSystem.setCmdError(true);
      return;
      }
    }

  dlist.getStringList("chain_ids", chains);

  join_domain = dst_domain->joinDomains(name, dst_desc, domain, desc, nterm, chains);
  pmSystem.addDomain(join_domain);
  }

//*============================================================*
//*==========             pm_CmdDomainMode           ==========*
//*============================================================*
// process domain mode command.

void
pm_CmdDomainMode (PmMolecule *domain, PmCmdDataList& dlist)
  {
  string dv, name, atom_name;
  PmCmdData data;

  PmVector3 color;
  PmAtomFilter filter;
  float scale;
  int number;
  bool show = true, show_set = true, vectors = false;

  while (dlist.getNext(data,dv)) {
    if (data.name == "atom_names") {
      dlist.getStringList("atom_names", filter.names);
      }
    else if (data.name == "color") {
      data.getVector(color);
      domain->setModeColor(color);
      }
    else if (data.name == "number") {
      number = data.getInt();
      domain->setModeNumber(number);
      }
    else if (data.name == "scale") {
      scale = data.getFloat();
      domain->setModeScale(scale);
      }
    else if (data.name == "show") {
      show = data.getBoolean();
      }
    else if (data.name == "vectors") {
      vectors = true;
      }
    }

  if (show_set) {
    if (vectors) { 
      domain->displayModeVector(name, filter, show);
      }
    }
  }

//*============================================================*
//*==========         pm_CmdDomainDefRegion          ==========*
//*============================================================*
// process domain define region command.

void
pm_CmdDomainDefRegion (PmMolecule *domain, string name, PmCmdDataList& dlist)
  {
  //pm_PrintMsg (CMDPR, "region \"%s\" defined. ", name.c_str());
  string desc, rstr;
  PmAtomFilter filter;
  bool bval, use_spheres, show;
  PmMolRegionParameters params;
  PmRegion *rgn;
  PmMoleculeRenderType rtype;
  vector<string> atom_names;

  //===== use centers of residue side chains =====//

  if (dlist.getBoolean("use_sidechains", params.use_sidechains)) {
    if (params.use_sidechains) { 
      domain->getBackboneAtomNames(filter.names);
      filter.exclude = true;
      }
    else {
      domain->getBackboneAtomNames(filter.names);
      filter.sidechain = false;
      filter.exclude = false;
      }
    }

  //===== use pca =====//

  else if (dlist.getBoolean("use_pca", bval)) {
    params.use_pca = bval;
    }

  else if (dlist.getStringList("atom_names", atom_names)) {
    if (atom_names.size() == 1) {
      if (atom_names[0] == "mc") {
        domain->getBackboneAtomNames(filter.names);
        }
      else if (atom_names[0] == "sc") {
        domain->getBackboneAtomNames(filter.names);
        filter.exclude = true;
        }
      else {
        filter.names = atom_names;
        }
      }
    else if (atom_names.size() != 0) {
      filter.names = atom_names;
      }
    }

  dlist.getBoolean("use_surface", params.use_surface); 
  dlist.getFloat("tolerance", params.tolerance); 

  //===== define region based on protein residues =====//

  if (dlist.getString("residues", desc)) {
    //fprintf (stderr, "######### desc: %s \n", desc.c_str() );
    if (!domain->checkDescriptor(desc)) {
      pm_ErrorReport (PM, "residue descriptor \"%s\" is not valid for this domain.", 
                      "*", desc.c_str());
      pmSystem.setCmdError(true);
      return;
      }

    domain->defineRegion(name, filter, desc, params);

    if (dlist.getString("render", rstr)) {
      PmMolecule::convMoleculeRenderType(rstr, rtype);
      }
    else {
      rtype = PM_MOLECULE_RENDER_SOLID;
      }

    if (dlist.getBoolean("show", show) && show) {
      PmVector3 color;

      if (!dlist.getVector("color", color)) {
        color.set(1,1,1);
        }

      dlist.getBoolean("use_spheres", use_spheres);
      domain->displayRegion(name, color, rtype, use_spheres); 
      }
    }
  }

//*============================================================*
//*==========         pm_CmdDomainFit                ==========*
//*============================================================*
// process domain fit command.

void
pm_CmdDomainFit (PmMolecule *domain, PmCmdDataList& dlist)
  {
  string dv, seq, name, dname, dst_seq, dst_names;
  PmAtomFilter filter, dst_filter;
  bool rot;
  PmMolecule *dst_domain, *dom_copy;
  vector<string> domain_names;
  vector<float> wt;
  float weight;
  //fprintf (stderr, ">>>>>>> pm_CmdDomainFit: \n");

  if (!dlist.getString("destination", dname)) {
    pm_ErrorReport (PM, "no destination domain given.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  pmSystem.getDomain(dname, &dst_domain); 

  if (!dst_domain) {
    pmSystem.getMolecule (dname, &dst_domain);

    if (!dst_domain) {
      pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
      pmSystem.setCmdError(true);
      return;
      }
    }

  if (!dlist.getString("sequence", seq)) {
    pm_ErrorReport (PM, "no sequence given.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  if (!dlist.getString("destination_sequence", dst_seq)) {
    pm_ErrorReport (PM, "no destination sequence given.", "*");
    pmSystem.setCmdError(true);
    return;
    }

  // define which atoms to use in fitting  //

  if (!dlist.getStringList("atom_names", filter.names)) {
    filter.names.push_back("CA");
    }

  if (!dlist.getBoolean("mainchain", filter.mainchain)) {
    filter.mainchain = false;
    }

  if (!dlist.getBoolean("sidechain", filter.sidechain)) {
    filter.sidechain = false;
    }

  if (!dlist.getBoolean("rotation", rot)) {
    rot = true;
    }

  // get atoms coordinates for the given sequences //

  vector<PmVector3> coords, dst_coords; 
  domain->getAtomCoords(seq, filter, coords);

  if (!coords.size()) {
    pm_ErrorReport (PM, "could not process sequence [%s]. no coordinates.", "*", seq.c_str());
    pmSystem.setCmdError(true);
    return;
    }

  if (dlist.getStringList("destination_atom_names", dst_filter.names)) {
    dst_domain->getAtomCoords(dst_seq, dst_filter, dst_coords);
    }
  else {
    dst_domain->getAtomCoords(dst_seq, filter, dst_coords);
    }

  if (!dst_coords.size()) {
    pm_ErrorReport (PM, "could not process destination sequence [%s].", "*", 
                    dst_seq.c_str());
    pmSystem.setCmdError(true);
    return;
    }

  if (dst_coords.size() != coords.size()) {
    pm_ErrorReport (PM, "sequences are not equal: sequence[%s] has %d atoms; destination sequnce[%s] has %d atoms.", "*", seq.c_str(), coords.size(), dst_seq.c_str(), 
      dst_coords.size());
    pmSystem.setCmdError(true);
    return;
    }

  if (dlist.getFloat("weight", weight)) {
    wt.push_back(weight);

    for (int i = 1; i < coords.size(); i++) {
      wt.push_back(1.0);
      }
    }

  // perform fit //

  float err, angs[3], rmat[3][3];
  PmXform xform;

  pm_MathFitRms (dst_coords, coords, wt, xform, err);
  domain->getName(name);

  fprintf (stderr, "\n------- fit \"%s\" to \"%s\" ------\n",name.c_str(),dname.c_str()); 
  fprintf (stderr, ">>> fit using %d coordinates \n", coords.size());
  fprintf (stderr, ">>> fit error = %g \n", err);
  fprintf (stderr, ">>> fit translation = %g %g %g \n", xform.translation[0], 
           xform.translation[1], xform.translation[2]); 
  fprintf (stderr, ">>> fit center = %g %g %g \n", xform.center[0], xform.center[1], xform.center[2]); 

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rmat[i][j] = xform.matrix(i,j);
      }
    }

  pm_MathExtractRotations(rmat, angs);
  fprintf (stderr, ">>> fit rotation angles = %f %f %f  \n", angs[0], angs[1], angs[2]);

  // don't rotate  //

  if (!rot) {
    pm_PrintMsg (CMDPR, "don't rotate ");
    xform.resetRotation();
    }

  //===== create a copy of the fit =====//

  if (dlist.getString("copy", dname)) {
    pm_PrintMsg (CMDPR, "create copy");
    domain->copy(dname, &dom_copy);
    pmSystem.addDomain(dom_copy);
    dom_copy->xformAtoms(xform);
    }
  else {
    domain->xformAtoms(xform);
    }

  //===== additional domains to transform =====//

  if (dlist.getStringList("list", domain_names)) {
    string dname;
    PmMolecule *dom;
    PmMassProperties props;

    for (unsigned int i = 0; i < domain_names.size(); i++) {
      dname = domain_names[i];
      pmSystem.getDomain (dname, &dom);

      if (!dom) {
        pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
        pmSystem.setCmdError(true);
        return;
        }

      //dom->getMassProps(props);
      //xform.center = props.com;
      dom->xformAtoms(xform);
      }
    }
  }

//*============================================================*
//*==========         pm_CmdDomainAtoms              ==========*
//*============================================================*
// process domain atoms command.

void
pm_CmdDomainAtoms (PmMolecule *domain, PmCmdDataList& dlist)
  {
  string dv, name, atom_name;
  vector<string> atom_names; 
  PmCmdData data;
  PmVector3 color;
  bool show, show_set;
  PmAtomFilter filter;

  show_set = true;
  show = true;

  dlist.getStringList("atom_names", atom_names);

  if (atom_names.size() == 1) {
    if (atom_names[0] == "mc") {
      domain->getBackboneAtomNames(filter.names);
      }
    else if (atom_names[0] == "sc") {
      domain->getBackboneAtomNames(filter.names);
      filter.exclude = true;
      }
    else {
      filter.names = atom_names;
      }
    }
  else if (atom_names.size() != 0) {
    filter.names = atom_names;
    }

  while (dlist.getNext(data,dv)) {
    if (data.name == "color") {
      if (dv == "element") {
        domain->setAtomColorType (PM_ATOM_COLOR_ELEMENT);
        }
      else if (dv == "chain") {
        domain->setAtomColorType (PM_ATOM_COLOR_CHAIN);
        }
      else if (dv == "force") {
        domain->setAtomColorType (PM_ATOM_COLOR_FORCE);
        }
      else if (dv == "mode_energy") {
        domain->setAtomColorType (PM_ATOM_COLOR_MODE_ENERGY);
        }
      else {
        data.getVector(color);
        domain->setColor(color);
        }
      }

    // display the geometry as:

    else if (data.name == "display") {

      if (dv == "point") {
        domain->setDisplayType(PM_MOLECULE_DISPLAY_POINT);
        }
      else if (dv == "cross") {
        domain->setDisplayType(PM_MOLECULE_DISPLAY_CROSS);
        }
      else if (dv == "sphere") {
        domain->setDisplayType(PM_MOLECULE_DISPLAY_SPHERE);
        }
      }

    // added name to geometry //

    else if (data.name == "name") {
      data.getString(name);
      }

    // print atoms  //

    else if (data.name == "print") {
      vector<PmAtom> atoms;
      string name;
      domain->getName(name);
      domain->getAtoms(atoms);
      fprintf (stderr, "\n---------- domain \"%s\" ---------- \n", name.c_str());
      fprintf (stderr, "number of atoms = %d \n\n", atoms.size());
      fprintf (stderr, "#     resid     name  \n");
      fprintf (stderr, "--------------------- \n");
      for (unsigned int i = 0; i < atoms.size(); i++) {
        fprintf (stderr, "%d      %d     %s \n", i+1, atoms[i].seq, atoms[i].name);
        }
      return;
      }

    // render sphere as: point | line | solid //

    else if (data.name == "render") {
      if (dv == "point") {
        domain->setRenderType(PM_MOLECULE_RENDER_POINT);
        }
      else if (dv == "line") {
        domain->setRenderType(PM_MOLECULE_RENDER_LINE);
        }
      else if (dv == "solid") {
        domain->setRenderType(PM_MOLECULE_RENDER_SOLID);
        }
      }

    // shading = none | flat | color | normal //

    /*
    else if (data.name== "shading") {
      PmGeometryShadingType stype;
      PmGraphicsGeometry::convShadingType (dv, stype);

      if (stype == PM_GEOMETRY_SHADING_UNKNOWN) {
        pm_ErrorReport (PM, "unknown shading type \"%s\".", "*", dv.c_str());
        }
      else {
        domain->setShadingType(stype);
        }
      }
     */

    else if (data.name == "show") {
      show = data.getBoolean();
      }
    }

  if (show_set) {
    domain->displayAtoms (name, filter, show);
    }
  }

//*============================================================*
//*==========         pm_CmdDomainBonds              ==========*
//*============================================================*
// process domain bonds command.

void
pm_CmdDomainBonds (PmMolecule *domain, PmCmdDataList& dlist)
  {
  string dv, name; 
  PmVector3 color, bond_color, atom_color;
  float width; 
  bool show, show_set;
  PmCmdData data;

  show_set = true;
  show = true;

  while (dlist.getNext (data,dv)) {
    if (data.name == "atom_color") {
      if (dv == "element") {
        domain->setAtomColorType (PM_ATOM_COLOR_ELEMENT);
        }
      else if (dv == "chain") {
        domain->setAtomColorType (PM_ATOM_COLOR_CHAIN);
        }
      else if (dv == "force") {
        domain->setAtomColorType (PM_ATOM_COLOR_FORCE);
        }
      else if (dv == "mode_energy") {
        domain->setAtomColorType (PM_ATOM_COLOR_MODE_ENERGY);
        }
      else {
        data.getVector(color);
        domain->setAtomColor(color);
        }
      }

    else if (data.name == "atoms") {
      bool bond_atoms = data.getBoolean();
      domain->setDisplayBondAtoms (bond_atoms);
      }

    else if (data.name == "cbonly") {
      //bool cbonly = data.getBoolean();  
      }

    else if (data.name == "base_planes") {
      //bool base_planes = data.getBoolean();  
      }

    else if (data.name == "bond_color") {
      if (dv == "element") {
        domain->setBondColorType (PM_ATOM_COLOR_ELEMENT);
        }
      else if (dv == "charge") {
        domain->setBondColorType (PM_ATOM_COLOR_CHARGE);
        }
      else if (dv == "charge-polar") {
        domain->setBondColorType (PM_ATOM_COLOR_CHARGE_POLAR);
       }
      else if (dv == "force") {
        domain->setBondColorType (PM_ATOM_COLOR_FORCE);
        }
      else {
        data.getVector(bond_color);
        domain->setBondColor(bond_color);
        }
      }

    // display = point | cross | sphere //

    else if (data.name== "display") {
      PmMoleculeDisplayType dtype;
      PmMolecule::convMoleculeDisplayType(dv, dtype);

      if (dtype == PM_MOLECULE_DISPLAY_UNKNOWN) {
        pm_ErrorReport (PM, "unknown display type \"%s\".", "*", dv.c_str());
        pmSystem.setCmdError(true);
        }
      else {
        domain->setDisplayType(dtype);
        }
      }

    else if (data.name == "line_width") {
      width = data.getFloat();
      domain->setLineWidth(width);
      }

    else if (data.name == "msize") {
      float size = data.getFloat();
      domain->setMarkerSize(size);
      }

    else if (data.name == "name") {
      data.getString(name);
      }

    // render = point | line | solid //

    else if (data.name == "render") {
      PmMoleculeRenderType rtype;
      PmMolecule::convMoleculeRenderType(dv, rtype);

      if (rtype == PM_MOLECULE_RENDER_UNKNOWN) {
        pm_ErrorReport (PM, "unknown render type \"%s\".", "*", dv.c_str());
        pmSystem.setCmdError(true);
        }
      else {
        domain->setBondAtomRenderType(rtype);
        }
      }

    else if (data.name == "show") {
      show = data.getBoolean();  
      }
    }

  if (show_set) {
    domain->displayBonds(name, show);
    }
  }

//*============================================================*
//*==========           pm_CmdDomainBackbone         ==========*
//*============================================================*
// process domain backbone command.

void
pm_CmdDomainBackbone (PmMolecule *domain, PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name;
  PmVector3 color;
  float width, msize;
  bool planes, angles, show, show_set, tube;

  planes = false;
  angles = false;
  show_set = true;
  show = true;
  fprintf (stderr, ">>>>>>> pm_CmdDomainBackbone: \n");
  color.set(1,1,1);

  while (dlist.getNext(data,dv)) {
    if (data.name == "color") {
      if (dv == "element") {
        domain->setAtomColorType (PM_ATOM_COLOR_ELEMENT);
        }
      else if (dv == "chain") {
        domain->setAtomColorType (PM_ATOM_COLOR_CHAIN);
        }
      else if (dv == "force") {
        domain->setAtomColorType (PM_ATOM_COLOR_FORCE);
        }
      else if (dv == "mode_energy") {
        domain->setAtomColorType (PM_ATOM_COLOR_MODE_ENERGY);
        }
      else {
        data.getColor(color);
        //data.getVector(color);
        domain->setColor(color);
        }
      }

    else if (data.name == "planes") {
      planes = data.getBoolean();
      }

    else if (data.name == "tube") {
      tube = data.getBoolean();
      domain->setTubeDisplay(tube);
      fprintf (stderr, ">>>>>>> pm_CmdDomainBackbone: set tube display \n");
      }

    else if (data.name == "plane_color") {
      data.getColor(color);
      }

    else if (data.name == "msize") {
      msize = data.getFloat();
      domain->setMarkerSize (msize);
      }

    else if (data.name == "name") {
      data.getString(name);
      }

    else if (data.name == "show") {
      show_set = true;
      show = data.getBoolean();
      }

    else if (data.name == "width") {
      width = data.getFloat();
      domain->setLineWidth(width);
      }
    }

  if (show_set) {
    if (planes ) {
      string res_desc;
      vector<PmResidue*> rlist;
      PmResidue* res = NULL;
      PmMolecule *mol;

      if (dlist.getString("residue", res_desc)) {
        pmSystem.getCurrentMolecule(&mol);
        mol->getResidues(res_desc, rlist);

        if (rlist.size() != 0) {
          res = rlist[0];
          }
        }

      domain->displayBackbonePlanes (name, res, color, show);
      }
    else {
      domain->displayBackbone (name, show);
      }
    }
  }

//*============================================================*
//*==========           pm_CmdDomainPca              ==========*
//*============================================================*
// process domain pca command.

void
pm_CmdDomainPca (PmMolecule *domain, PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name;

  bool show = true;
  PmPcaResults pca;
  PmAtomFilter filter;

  // get atoms names for computation.  //

  dlist.getStringList("atom_names", filter.names);

  domain->compPrincipalComponents(filter, pca);
  pm_PrintMsg (CMDPR, "center = %g %g %g", pca.com[0],pca.com[1],pca.com[2]);
  pm_PrintMsg (CMDPR, "axis1 = %g %g %g ", pca.axis1[0],pca.axis1[1],pca.axis1[2]);
  pm_PrintMsg (CMDPR, "axis2 = %g %g %g ", pca.axis2[0],pca.axis2[1],pca.axis2[2]);
  pm_PrintMsg (CMDPR, "axis3 = %g %g %g ", pca.axis3[0],pca.axis3[1],pca.axis3[2]);
  pm_PrintMsg (CMDPR, "scales = %g %g %g ", pca.s1, pca.s2, pca.s3);

  if (dlist.getBoolean("show", show)) {

    if (show) {
      domain->displayPca(pca, show);
      }
    }
  }

//*============================================================*
//*==========         pm_CmdDomainSurface            ==========*
//*============================================================*
// process domain atoms command.

void
pm_CmdDomainSurface (PmMolecule *domain, PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name;
  PmSurface *surf;
  PmDbInterfaceSelect db_sel;
  PmVector3 color;
  float msize;
  int show_set = true;
  int show = true;

  // get domain surface //

  domain->getSurface(&surf);

  while (dlist.getNext (data,dv)) {
    if (surf && data.name == "color") {

      if (dv == "charge") {
        surf->setMapCharge(true);
        }
      else {
        data.getVector(color);
        surf->setColor(color);
        }
      }

    // display = point | line | solid //

    else if (surf && data.name== "display") {
      PmGeometryDisplayType dtype;
      PmGraphicsGeometry::convDisplayType (dv, dtype);

      if (dtype == PM_GEOMETRY_DISPLAY_UNKNOWN) {
        pm_ErrorReport (PM, "unknown display type \"%s\".", "*", dv.c_str());
        pmSystem.setCmdError(true);
        }
      else {
        surf->setDisplayType(dtype);
        }
      }

    else if (surf && data.name == "lighting") {
      if (dv == "two-sided") {
        surf->setLighting(true);
        }
      }

    else if (surf && data.name == "msize") {
      msize = data.getFloat();
      domain->setMarkerSize (msize);
      surf->setMarkerSize (msize);
      }

    // shading = none | flat | color | normal //

    else if (surf && data.name== "shading") {
      PmGeometryShadingType stype;
      PmGraphicsGeometry::convShadingType (dv, stype);

      if (stype == PM_GEOMETRY_SHADING_UNKNOWN) {
        pm_ErrorReport (PM, "unknown shading type \"%s\".", "*", dv.c_str());
        }
      else {
        surf->setShadingType(stype);
        }
      }

    else if (data.name == "name") {
      data.getString(name);
      }

    else if (data.name == "read") {
      string fname, name, dname, format;
      domain->getName(dname);
      dlist.getString("file", fname);
      PmDbInterface *db = db_sel.create(dname, PM_DB_PM);
      db->open(fname, PM_DB_MODE_READ, "");
      db->read(name);
      db->close();
      surf = db->getSurface();

      if (!surf) {
        pm_ErrorReport (PM, "surface not processed.", "*");
        return;
        }

      if (dlist.getString("format", format)) {
        }

      domain->setSurface (surf);
      }
    else if (data.name == "show") {
      show = data.getBoolean();
      }
    }

  if (show_set) {
    domain->displaySurface(name, show);
    }
  }

//*============================================================*
//*==========         pm_CmdCheckOutputDesc          ==========*
//*============================================================*
// check an output descriptor.

bool
pm_CmdCheckOutputDesc(string desc)
  {
  string src, dst;
  int n = desc.length();
  size_t found = desc.find('=');

  if (found != string::npos) {
    src.assign(desc, 0, found);
    dst.assign(desc, found+1, n-found+1);
    }
  else {
    src = desc;
    }

  for (unsigned int i = 0; i < src.length(); i++) {
    char c = src[i];

    if ( !isdigit(c) && !isalpha(c) && (c != ']') && (c != '[') && (c != '-') 
         && (c != ',') ) {
      pm_ErrorReport (PM, "invalid chain descriptor char %c in %s.", "*", c, src.c_str());
      return false;
      }
    }

  for (unsigned int i = 0; i < dst.length(); i++) {
    char c = dst[i];

    if ( !isdigit(c) && !isalpha(c) && (c != ']') && (c != '[')) {
      pm_ErrorReport (PM, "invalid chain descriptor char %c in %s.", "*", c, dst.c_str());
      return false;
      }
    }

  return true;
  }

//*============================================================*
//*==========              pm_CmdDomains             ==========*
//*============================================================*
// process molecular domains commands.                        

void
pm_CmdDomains (PmCmdDataList& dlist)
  {
  string dv, str, one_chain, start_res, fname, atypes; 
  PmMolecule *domain; 
  vector<string> domain_names;
  vector<PmMolecule*> domains; 
  string chain_id, sidechains, res_names, model_name;
  bool mainchain, no_hydrogen;
  vector<string> chains;
  PmAtomFilter filter;
  PmDbModeType mtype; 
  PmCmdData data, pdata;
  int model=-1;

  dlist.getNext (data, dv);

  //===== write or append a domain =====//

  if ((data.name == "write") || (data.name == "append")) {
    if (!dlist.getString("file", fname)) {
      return;
      }

    // get write mode //

    mtype = PmDbInterface::convModeType(data.name);

    if (!dlist.getStringList("list", domain_names)) {
      pm_ErrorReport (PM, "domain names not given.", "*");
      pmSystem.setCmdError(true);
      return;
      }

    //===== get domains to write =====//

    domains.clear();

    for (unsigned int i = 0; i < domain_names.size(); i++) {
      string dname = domain_names[i];
      pmSystem.getDomain (dname, &domain);

      if (!domain) {
        pm_ErrorReport (PM, "no domain named \"%s\".", "*", dname.c_str());
        pmSystem.setCmdError(true);
        return;
        }
 
      domains.push_back(domain);
      }

    //===== write domain atoms =====//

    if (dlist.getString("atoms", str)) { 
      dlist.getString("types", atypes);
      dlist.getStringList("chains", chains);

      if (chains.size()) {
        for (unsigned int i = 0; i < chains.size(); i++) {
          string desc = chains[i];

          if (!pm_CmdCheckOutputDesc(desc)) {
            pm_ErrorReport (PM, "improper chain descriptor \"%s\".", "*", desc.c_str());
            pmSystem.setCmdError(true);
            return;
            }
          }
        }

      dlist.getString("one_chain", one_chain);
      dlist.getBoolean("mainchain", mainchain);
      dlist.getBoolean("no_hydrogen", no_hydrogen);
      dlist.getString("start_res", start_res);
      dlist.getString("res_names", res_names);
      dlist.getInt("model", model);
      dlist.getString("model_name", model_name);

      PmDbInterfaceSelect db_sel;
      PmDbInterface *db = db_sel.create("tmp", PM_DB_PDB);
      db->open(fname, mtype, "");
      db->writeDomainAtoms(domains, chains, one_chain, no_hydrogen, start_res, res_names,
                           model, model_name);
      db->close();
      }

    //===== write domain sequences =====//

    else if (dlist.getString("sequence", str)) { 
      PmDbInterfaceSelect db_sel;
      PmDbInterface *db = db_sel.create("tmp", PM_DB_PDB);
      db->open(fname, mtype, "");
      db->writeDomainSequence(domains);
      db->close();
      }
    }

  //===== create =====//

  else if (data.name == "create") { 
    bool flag;

    if (dlist.getBoolean("helices", flag)) {
      string prefix, dname, hname;
      vector<PmHelix> helices;
      PmMolecule *mol, *domain;
      stringstream dss;
      PmExtent extent;
      PmVector3 color;
      float width;

      if (!dlist.getString("prefix", prefix)) {
        pm_ErrorReport (PM, "no prefix given.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      pmSystem.getCurrentMolecule(&mol);
      mol->getHelices(helices);

      if (!helices.size()) {
        pm_ErrorReport (PM, "no helices for molecule.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      if (!dlist.getVector("color", color)) {
        color.set(1,1,1);
        }

      if (!dlist.getFloat("width", width)) {
        width = 1.0;
        }

      for (unsigned int i = 0; i < helices.size(); i++) {
        dss << prefix << i+1;
        dname = dss.str();
        //dss.str().erase();
        dss.str(std::string());

        dss << "helix" << helices[i].id;
        hname = dss.str();
        //dss.str().erase();
        dss.str(std::string());
        /*
        fprintf (stderr, "    >>> domain \"%s\"  \"%s\" \n", dname.c_str(), 
                 hname.c_str());
        */

        domain = mol->createDomain (dname, hname, filter);

        if (domain) {
          domain->getExtent(extent);
          pmSystem.updateExtent (extent);
          domain->setColor(color);
          domain->setLineWidth(width);
          domain->displayBackbone("", true);
          }
        }
      }

    // create domains from all loops //

    else if (dlist.getBoolean("loops", flag)) {
      string prefix, dname, hname;
      vector<PmLoop> loops;
      PmMolecule *mol, *domain;
      stringstream dss;
      PmExtent extent;
      PmVector3 color;
      float width;

      if (!dlist.getString("prefix", prefix)) {
        pm_ErrorReport (PM, "no prefix given.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      pmSystem.getCurrentMolecule(&mol);
      mol->getLoops(loops);

      if (!dlist.getVector("color", color)) {
        color.set(1,1,1);
        }

      if (!dlist.getFloat("width", width)) {
        width = 1.0;
        }

      for (unsigned int i = 0; i < loops.size(); i++) {
        dss << prefix << i+1;
        dname = dss.str();
        //dss.str().erase();
        dss.str(std::string());

        dss << "loop" << loops[i].id;
        hname = dss.str();
        //dss.str().erase();
        dss.str(std::string());
        /*
        fprintf (stderr, "    >>> domain \"%s\"  \"%s\" \n", dname.c_str(), 
                 hname.c_str());
        */

        domain = mol->createDomain (dname, hname, filter);

        if (domain) {
          domain->getExtent(extent);
          pmSystem.updateExtent (extent);
          domain->setColor(color);
          domain->setLineWidth(width);
          domain->displayBackbone("", true);
          }
        }
      }

    // create domains from all sheets //

    else if (dlist.getBoolean("sheets", flag)) {
      string prefix, dname, hname;
      vector<PmSheet> sheets;
      PmMolecule *mol, *domain;
      stringstream dss;
      PmExtent extent;
      PmVector3 color;
      float width;

      if (!dlist.getString("prefix", prefix)) {
        pm_ErrorReport (PM, "no prefix given.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      pmSystem.getCurrentMolecule(&mol);
      mol->getSheets(sheets);

      if (!sheets.size()) {
        pm_ErrorReport (PM, "no sheets for molecule.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      if (!dlist.getVector("color", color)) {
        color.set(1,1,1);
        }

      if (!dlist.getFloat("width", width)) {
        width = 1.0;
        }

      for (unsigned int i = 0; i < sheets.size(); i++) {
        dss << prefix << i+1;
        dname = dss.str();
        //dss.str().erase();
        dss.str(std::string());

        dss << "sheet" << sheets[i].id << ':' << sheets[i].num; 
        hname = dss.str();
        //dss.str().erase();
        dss.str(std::string());
        /*
        fprintf (stderr, "    >>> domain \"%s\"  \"%s\" \n", dname.c_str(), 
                 hname.c_str());
        */

        domain = mol->createDomain (dname, hname, filter);

        if (domain) {
          domain->getExtent(extent);
          pmSystem.updateExtent (extent);
          domain->setColor(color);
          domain->setLineWidth(width);
          domain->displayBackbone("", true);
          }
        }
      }
    }

  //===== pick =====//

  else if (data.name == "pick") {
    //pm_GrPickActivate (PM_ENT_DOMAIN, true);
    return;
    }
  }

//*============================================================*
//*==========         pm_CmdDomainXform              ==========*
//*============================================================*
// process domain xform command.

void
pm_CmdDomainXform (PmMolecule *domain, PmCmdDataList& dlist)
  {
  string dv;
  PmCmdData data, pdata;
  PmVector3 translation, rotation, axis, point;
  PmXform xform;
  PmMassProperties props;
  float angle = 0.0;
  bool rot_about_axis = false, use_pca;
  PmMolecule *dom_copy = NULL;

  // set center of xform  //

  domain->getMassProps(props);
  xform.center = props.com;

  while (dlist.getNext (data, dv)) {
    if (data.name == "angle") {
      angle = data.getFloat();
      }

    else if (data.name == "axis") {
      data.getVector(axis);
      rot_about_axis = true;
      }

    else if (data.name == "atoms") {
      vector<string> atoms_desc;
      PmAtom atom1, atom2;

      if (!dlist.getStringList("atoms", atoms_desc)) {
        pm_ErrorReport (PM, "atoms is not a string list.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      if (atoms_desc.size() != 4) {
        pm_ErrorReport (PM, 
          "wrong format for atoms:  <resname> <atom> <resname> <atom> ", "*");
        pmSystem.setCmdError(true);
        return;
        }

      if (!domain->getAtoms(atoms_desc, atom1, atom2)) {
        pm_ErrorReport (PM, "atoms are not valid.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      axis = atom2.pos - atom1.pos;
      axis.normalize();
      rot_about_axis = true;
      xform.center = 0.5*(atom1.pos + atom2.pos); 
      xform.center = atom2.pos; 
      xform.center = atom1.pos; 
      pm_PrintMsg (CMDPR, "atom1 = (%g %g %g) ", atom1.pos[0], atom1.pos[1], 
                   atom1.pos[2]); 
      pm_PrintMsg (CMDPR, "atom2 = (%g %g %g) ", atom2.pos[0], atom2.pos[1], 
                  atom2.pos[2]); 
      pm_PrintMsg (CMDPR, "axis = (%g %g %g) ", axis[0], axis[1], axis[2]); 
      }

    else if (data.name == "copy") {
      domain->copy(dv, &dom_copy);
      pmSystem.addDomain(dom_copy);
      }


    // xform point to origin //

    else if ((data.name == "origin")) { 
      dlist.getData("origin", pdata);

      if (!pm_CmdProcPosition(pdata, point)) {
        pm_ErrorReport (PM, "couldn't process point data.", "*");
        pmSystem.setCmdError(true);
        return;
        }

      xform.translation = -point;
      }

    else if ((data.name == "translation") || 
             (data.name == "translate")) {
      data.getVector(xform.translation);

      if (dlist.getData("point", pdata)) {
        if (!pm_CmdProcPosition(pdata, point)) {
          pm_ErrorReport (PM, "couldn't process point data.", "*");
          pmSystem.setCmdError(true);
          return;
          }

        xform.translation = xform.translation - point;
        }
      }

    else if ((data.name == "rotation") || 
             (data.name == "rotate")) {
      PmVector3 angles;
      data.getVector(angles);
      xform.setRotation(angles[0], angles[1], angles[2]);
      }

    // use pca //

    else if (data.name == "use_pca") { 
      string dname, desc;
      PmMolecule *domain;
      PmAtomFilter filter;
      PmPcaResults pca;
      vector<PmVector3> coords;
      float wd[3], mag;
      PmExtent proj;
      PmVector3 vec;
      int pca_axis;

      use_pca = data.getBoolean();
      if (!use_pca) continue; 

      if (!dlist.getString("pca_domain", dname)) {
        pm_ErrorReport (PM, "no pca domain given.", "*");
        pmSystem.setCmdError(true);
        return;
        }
     
      pmSystem.getDomain (dname, &domain);

      if (!domain) {
        pm_ErrorReport (PM, "no pca domain named \"%s\".", "*", dname.c_str());
        pmSystem.setCmdError(true);
        return;
        }

      dlist.getString("pca_res", desc); 
      dlist.getStringList("atom_names", filter.names);
      domain->getAtomCoords(desc, filter, coords);

      if (coords.size() == 0) {
        pm_ErrorReport (PM, "no coordinates for pca.", "*");
        return;
        }

      pm_MathPrincipalComponents(coords, pca);
      pm_MathPrincipalComponentsProj(coords, pca, proj, wd);

      // align the domains principle axes with the vector //

      if (dlist.getVector("align_pca", vec)) { 
        PmVector3 v1, v2, v3, a1, a2, a3;
        PmMatrix3x3 mat;
        float ang1, ang2, ang3;
        xform.center = props.com;

        v1 = vec;
        pm_MathBasisCompute (v1, v2, v3);

        a1 = pca.axis1;
        a2 = pca.axis2;
        a3 = pca.axis3;

        a1.normalize();
        a2.normalize();
        a3.normalize();

        pm_MathFrameRotation(a1, a2, a3, v1, v2, v3, mat);
        xform.setRotation(mat);
        }

      //===== translate =====//

      else {
        if (dlist.getVector("translate", vec)) { 
          xform.translation = vec[0]*pca.axis1 + vec[1]*pca.axis2 + vec[2]*pca.axis3;
          }

        if (dlist.getInt("pca_axis", pca_axis)) { 
          if (pca_axis == 1) {
            axis = pca.axis1;
            }
          else if (pca_axis == 2) {
            axis = pca.axis2;
            }
          else if (pca_axis == 3) {
            axis = pca.axis3;
            }
          else {
            pm_ErrorReport (PM, "pca axis must be between 1 and 3.", "*");
            pmSystem.setCmdError(true);
            return;
            }

          rot_about_axis = true;
          xform.center = props.com;
          }
        }
      }
    }

  //===== apply the transformation =====//

  if (rot_about_axis) {
    xform.setRotation(angle, axis);
    }

  if (dom_copy) {
    dom_copy->xformAtoms (xform);
    }
  else {
    domain->xformAtoms (xform);
    }
  }

//*============================================================*
//*==========             pm_CmdDomain               ==========*
//*============================================================*
// process molecular domain commands.

void
pm_CmdDomain (PmCmdDataList& dlist)
  {

  string dv, name, domain_name, residues; 
  vector<string> domain_names; 
  PmAtomFilter filter;
  bool sidechain, peptide, backbone, sidechain_group, last_nca;
  PmMolecule *mol, *domain; 
  PmExtent extent;
  PmCmdData data;

  dlist.getNext (data, dv);

  //fprintf (stderr, ">>>>>> pm_CmdDomain: [%s] \n", data.name.c_str());

  // ===== create ===== //

  if (data.name == "create") {
    if (dlist.getString("molecule", name)) {
      pmSystem.getMolecule (name, &mol); 

      if (!mol) {
        pm_ErrorReport (PM, "no molecule named \"%s\".", "*", name.c_str());
        pmSystem.setCmdError(true);
        return;
        }
      }
    else { 
      pmSystem.getMolecule (&mol); 
      }

    if (!mol) {
      pm_ErrorReport (PM, "no current molecule", "*");
      pmSystem.setCmdError(true);
      return;
      }

    // create domain from name = <value> format //

    else if (dlist.getString("name", domain_name)) {

      // create a domain from two or more domains //

      if (dlist.getStringList("list", domain_names)) {
        domain = mol->createDomain (domain_name, domain_names);
        }

      else if (!dlist.getString("residues", residues)) {
        pm_ErrorReport (PM, "no residues given", "*");
        pmSystem.setCmdError(true);
        return;
        }
      }

    // command of the form <name> = <res> //

    else {
      dlist.getNext (data, residues);
      domain_name = data.name;
      }

    if (mol) {
      dlist.getStringList("atom_names", filter.names);

      if (dlist.getBoolean("sidechain", sidechain)) {
        mol->getBackboneAtomNames(filter.names);
        filter.exclude = true;
        }

      else if (dlist.getBoolean("backbone", backbone)) {
        mol->getBackboneAtomNames(filter.names);
        }
      }

    if (dlist.getBoolean("peptide", peptide)) {
      filter.peptide = peptide; 
      }

    else if (dlist.getBoolean("sidechain_group", sidechain_group)) {
      filter.sidechain_group = sidechain_group; 
      }

    else if (dlist.getBoolean("last_nca_atoms", last_nca)) {
      filter.last_nca_atoms = last_nca; 
      }

    if (domain_name.size() && residues.size()) {
      domain = mol->createDomain (domain_name, residues, filter);
      //fprintf (stderr, "\n>>> create domain from residues  \n");
      }

    // set extent //

    if (domain) {
      domain->getExtent(extent);
      pmSystem.updateExtent (extent);
      }
    else {
      pm_ErrorReport (PM, "domain not created.", "*");
      pmSystem.setCmdError(true);
      }

    //dlist.getNext(data,dv); 
    }

  else {
    pmSystem.getDomain (data.name, &domain);
    }

  if (!domain) {
    pm_ErrorReport (PM, "no domain named \"%s\".", "*", data.name.c_str());
    pmSystem.setCmdError(true);
    return;
    }

  pmSystem.setCurrentDomain (domain);

  while (dlist.getNext(data,dv)) {
     //fprintf (stderr, "data.name = %s \n", data.name.c_str());

    // ===== atoms ===== //

    if (data.name == "atoms") {
      pm_CmdDomainAtoms (domain, dlist);
      return;
      }

    // ===== backbone ===== //

    else if (data.name == "backbone") {
      pm_CmdDomainBackbone (domain, dlist);
      return;
      }

    // ===== bonds ===== //

    else if (data.name == "bonds") {
      pm_CmdDomainBonds (domain, dlist);
      return;
      }

    // ===== check ===== //

    else if (data.name == "check") {
      pm_CmdDomainCheck(domain, dlist);
      return;
      }

    // ===== contact ===== //

    else if (data.name == "contact") {
      pm_CmdDomainContact(domain, dlist);
      return;
      }

    // ===== copy ===== //

    else if (data.name == "copy") {
      PmMolecule *dom_copy; 
      domain->copy(dv, &dom_copy);
      pmSystem.addDomain(dom_copy);
      return;
      }

     // ===== define ===== //

    else if (data.name == "define") {
      dlist.getNext(data, dv);

      if (data.name == "region") {
        pm_CmdDomainDefRegion(domain, dv, dlist);
        }

      return;
      }

    // ===== fit ===== //

    else if (data.name == "fit") {
      pm_CmdDomainFit (domain, dlist);
      return;
      }

    // ===== hbonds ===== //

    else if (data.name == "hbonds") {
      pm_CmdDomainHbonds(domain, dlist);
      return;
      }

    // ===== helix ===== //

    else if (data.name == "helix") {
      pm_CmdDomainHelix(domain, dlist);
      return;
      }

    // ===== join ===== //

    else if (data.name == "join") {
      pm_CmdDomainJoin (domain, dlist);
      return;
      }

    // ===== mode ===== //

    else if (data.name == "mode") {
      pm_CmdDomainMode(domain, dlist);
      return;
      }

    // ===== pca ===== //

    else if (data.name == "pca") {
      pm_CmdDomainPca (domain, dlist);
      return;
      }

    // ===== print ===== //

    else if (data.name == "print") {
      pm_CmdDomainPrint (domain, dlist);
      return;
      }

    // ===== surface ===== //

    else if (data.name == "surface") {
      pm_CmdDomainSurface (domain, dlist);
      return;
      }

    // ===== topology ===== //

    else if (data.name == "topology") {
      pm_CmdDomainTopology(domain, dlist);
      return;
      }


    // ===== xform ===== //

    else if (data.name == "xform") {
      pm_CmdDomainXform (domain, dlist);
      return;
      }
    }
  }

