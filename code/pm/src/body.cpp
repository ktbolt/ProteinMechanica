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

//*============================================================*
//* body:                    b o d y                           *
//*============================================================*
// rigid body object.                       

#include "body.h"

namespace ProteinMechanica {

typedef struct PmBodyContact {
  vector<PmVector3> coords;
  vector<float> rads;
  vector<PmAtom> atoms;
  PmVector3 center;
  float radius;
  } PmBodyContact;

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*
// create a body object from a name and type.

PmBody::PmBody(const string name, PmBodyType type)
  {
  this->name = name;
  this->type = type;
  pobj = NULL;
  has_damping = false;
  damping_factor = 1.0; 
  damping_length = 1.0;
  }

//*============================================================*
//*==========            convBodyType                ==========*
//*============================================================*
// convert string type to enum.

void
PmBody::convBodyType(const string str, PmBodyType& type)
  {
  if (str == "ground") {
    type = PM_BODY_GROUND;
    }
  else if (str == "static") {
    type = PM_BODY_STATIC;
    }
  else if (str == "rigid") {
    type = PM_BODY_RIGID;
    }
  else if (str == "deformed") {
    type = PM_BODY_DEFORM;
    }
  else { 
    type = PM_BODY_UNKNOWN;
    }
  }

//*============================================================*
//*==========            getType                     ==========*
//*============================================================*
// get body type string.

void 
PmBody::getType(string& str)
  {
  switch (type) {
    case PM_BODY_GROUND: 
      str = "ground"; 
    break;

    case PM_BODY_STATIC:
      str = "static" ;
    break;

    case PM_BODY_RIGID:
      str = "rigid" ;
    break;

    case PM_BODY_DEFORM:
      str = "deformed" ;
    break;

    default:
      str = "unknown" ;
    break;
    }
  }

//*============================================================*
//*==========               isGround                 ==========*
//*============================================================*
// check for a ground body.

bool 
PmBody::isGround() { 
  return type == PM_BODY_GROUND; 
  }

//*============================================================*
//*==========               isStatic                 ==========*
//*============================================================*
// check for a static body.

bool 
PmBody::isStatic() { 
  return type == PM_BODY_STATIC; 
  }

//*============================================================*
//*==========               getXform                 ==========*
//*============================================================*

void 
PmBody::getXform(PmXform& xform) {
  xform = this->xform;
  }

//*============================================================*
//*==========               checkContact             ==========*
//*============================================================*
// check if all bodies are in contact. use bounding sphere.  

void
PmBody::checkContact(PmVector3& color, bool show)
  {
  string bname1, bname2;
  int num_bodies, n;
  vector<PmBody*> bodies;
  PmBody* body, *body1, *body2;
  PmPhysicalObj *pobj;
  PmVector3 center, center1, center2, v, pt1, pt2;
  PmExtent extent;
  float d, dx, dy, dz, rad, srad1, srad2, r1, r2;
  PmBodyType type;
  int num_contact;
  vector<int> ids1, ids2;

  fprintf (stderr, "\n---------- check contact ----------\n");
  pmSystem.getBodies(bodies);
  num_bodies = bodies.size();

  if (num_bodies == 0) {
    return;
    }

  PmBodyContact *bcon = new PmBodyContact[num_bodies];

  for (int i = 0; i < num_bodies; i++) {
    body = bodies[i];
    type = body->getType();

    if (type != PM_BODY_RIGID) {
      continue;
      }

    body->getPhysicalObj(&pobj);

    if (!pobj) {
      continue;
      }

    PmMolecule *mol = dynamic_cast<PmMolecule*>(pobj);

    if (mol) {
      mol->getAtoms(bcon[i].atoms);
      }

    // compute bounding sphere of current (transformed) coordinates for each body  //

    body->getCurrentCoordinats (bcon[i].coords);
    pobj->getRadii(bcon[i].rads);
    center.set(0,0,0);
    extent.min[0] =  1e6; extent.min[1] =  1e6; extent.min[2] = 1e6;
    extent.max[0] = -1e6; extent.max[1] = -1e6; extent.max[2] = -1e6;
    n = bcon[i].coords.size();

    for (int j = 0; j < n; j++) {
      center = center + bcon[i].coords[j];
      extent.update(bcon[i].coords[j]);
      }

    center = center / n;
    dx = extent.max[0] - extent.min[0];
    dy = extent.max[1] - extent.min[1];
    dz = extent.max[2] - extent.min[2];

    if (dx > dy) {
      rad = dx;
      }
    else {
      rad = dy;
      }

    if (dz > rad) {
      rad = dz;
      }

    rad = rad / 2.0 + 0.2;      // add 2 angstroms for atom radii //
    bcon[i].center = center;
    bcon[i].radius = rad;
    }

  //===== now check for contact =====//

  for (int i = 0; i < num_bodies-1; i++) {
    body1 = bodies[i];
    type = body1->getType();

    if (type != PM_BODY_RIGID) {
      continue;
      }

    body1->getName(bname1);
    center1 = bcon[i].center;
    srad1 = bcon[i].radius;
    fprintf (stderr, ">>> check contact of %s with: \n", bname1.c_str());

    for (int j = i+1; j < num_bodies; j++) {
      body2 = bodies[j];
      type = body2->getType();

      if (type != PM_BODY_RIGID) {
        continue;
        }

      ids1.clear();
      ids2.clear();

      body2->getName(bname2);
      center2 = bcon[j].center;
      srad2 = bcon[j].radius;
      fprintf (stderr, "    %s  ", bname2.c_str());

      num_contact = 0;
      v = center1 - center2;
      d = v.length();

      if ((d - srad1 - srad2) < 0.0) {
        fprintf (stderr, " overlaps ");

        for (unsigned int ic = 0; ic < bcon[i].coords.size(); ic++) {
          pt1 = bcon[i].coords[ic];
          r1 = bcon[i].rads[ic];

          for (unsigned int jc = 0; jc < bcon[j].coords.size(); jc++) {
            pt2 = bcon[j].coords[jc];
            r2 = bcon[j].rads[jc];
            v = pt1 - pt2;
            d = v.length();

            if ((d - r1 - r2) < 0.0) {
              ids1.push_back(ic);
              ids2.push_back(jc);
              num_contact += 1;
              }
            }
          }

        fprintf (stderr, "    number of contacts = %d \n", num_contact); 

        if (num_contact && show) {
          PmGraphicsAttributes atts;
          PmVector3 ccolor;
          int maxn, n1, n2;
          string geom_name = "contact";
          geom_name = geom_name + '[' + bname1 + ']' + "points";

          n1 = ids1.size();
          n2 = ids2.size();
          maxn = (n1 > n2 ? n1 : n2);

          PmVector3 *verts = new PmVector3[maxn];
          float *srads = new float[maxn];

          for (unsigned int ic = 0; ic < n1; ic++) {
            verts[ic] = bcon[i].coords[ids1[ic]];
            srads[ic] = bcon[i].rads[ids1[ic]];
            }

          PmGraphicsSphere *geom = new PmGraphicsSphere(geom_name, n1, verts, srads);
          atts.setColor(color);
          atts.setDisplayType(PM_GEOMETRY_DISPLAY_SOLID);
          geom->setAttributes(atts);
          geom->display();

          // in contact with //

          string cgeom_name = "contact";
          cgeom_name = cgeom_name + '[' + bname2 + ']' + "points";

          for (unsigned int ic = 0; ic < n2; ic++) {
            verts[ic] = bcon[j].coords[ids2[ic]];
            srads[ic] = bcon[j].rads[ids2[ic]];
            }

          for (int ic = 0; ic < 3; ic++) {
            ccolor[ic] = color[ic] - 0.5;
            if (ccolor[ic] < 0.0) ccolor[ic] = 0.0;
            }

          geom = new PmGraphicsSphere(cgeom_name, n2, verts, srads);
          atts.setColor(ccolor);
          atts.setDisplayType(PM_GEOMETRY_DISPLAY_SOLID);
          geom->setAttributes(atts);
          geom->display();

          delete[] verts;
          delete[] srads;
          }
        }
      else {
        fprintf (stderr, " no overlap \n");
        }
      }

    fprintf (stderr, "\n");
    }

  delete[] bcon;
  }

//*============================================================*
//*==========               checkContact             ==========*
//*============================================================*
// check if two bodies are in contact. use the coordinates and
// radii for each body.

void
PmBody::checkContact(PmBody *cbody, float tol, PmVector3& color1, PmVector3& color2,
                     bool show)
  {
  PmPhysicalObj *cpobj;
  vector<PmVector3> coords, ccoords;
  vector<float> rads, crads;
  PmXform xform, cxform;
  PmRigidSimResults simr, csimr;
  PmVector3 com, ccom;
  int num_contact;
  PmVector3 p1, p2, v;
  float r1, r2, d, overlap, max_overlap;
  vector<int> ids, cids;
  string cname;

  cbody->getPhysicalObj(&cpobj);

  if (!cpobj) {
    return;
    }

  cbody->getName(cname);
  fprintf (stderr, "\n---------- checking contact ---------- \n");
  fprintf (stderr, ">>> body=\"%s\" with body=\"%s\" \n", this->name.c_str(),
           cname.c_str());
  fprintf (stderr, ">>> tolerance=%f \n", tol);


  // get radii and current (transformed) coordinates for each body's  //

  this->getCurrentCoordinats (coords);
  pobj->getRadii(rads);
  cbody->getCurrentCoordinats (ccoords);
  cpobj->getRadii(crads);

  PmMolecule *mol = dynamic_cast<PmMolecule*>(pobj);
  PmMolecule *cmol = dynamic_cast<PmMolecule*>(cpobj);
  vector<PmAtom> atoms, catoms;
  PmAtomFilter filter;

  if (mol && cmol) {
    mol->getAtoms(atoms);
    cmol->getAtoms(catoms);
    mol->getBackboneAtomNames(filter.names);
    }

  // compute contact // 

  num_contact = 0;
  max_overlap = 0.0;

  for (unsigned int i = 0; i < coords.size(); i++) {
    p1 = coords[i];
    r1 = rads[i];

    for (unsigned int j = 0; j < ccoords.size(); j++) {
      p2 = ccoords[j];
      r2 = crads[j];
      v = p1 - p2;
      d = v.length();
      overlap = d - r1 - r2; 

      if (overlap < tol) {
        if (!atoms.empty() && !catoms.empty()) {
          if (atoms[i].seq != catoms[j].seq) {

            if (!filter.hasName(atoms[i].name) && !filter.hasName(catoms[j].name)) {
              fprintf (stderr, " atoms:seq=(%d %d) names=(%s %s) overlap=%g\n", 
                       atoms[i].seq, catoms[j].seq, atoms[i].name, catoms[j].name,
                       fabs(overlap));
              num_contact += 1;
              ids.push_back(i);
              cids.push_back(j);
              if (fabs(overlap) > max_overlap) max_overlap = fabs(overlap);
              }
            }
          }

        else {
          num_contact += 1;
          ids.push_back(i);
          cids.push_back(j);
          if (fabs(overlap) > max_overlap) max_overlap = fabs(overlap);
          }
        }
      }
    }

  fprintf (stderr, ">>> number of contacts=%d  max overlap=%g \n", num_contact, 
           max_overlap);

  if (!show || !num_contact) {
    return;
    }

  // show this body's contact //
 
  PmGraphicsAttributes atts;
  string geom_name = "contact";
  geom_name = geom_name + '[' + name + ']' + "points";

  int num_verts = ids.size();
  PmVector3 *verts = new PmVector3[num_verts];
  float *srads = new float[num_verts];

  for (unsigned int i = 0; i < ids.size(); i++) {
    verts[i] = coords[ids[i]];
    srads[i] = rads[ids[i]];
    }

  PmGraphicsSphere *geom = new PmGraphicsSphere(geom_name, num_verts, verts, srads);
  atts.setColor(color1);
  atts.setDisplayType(PM_GEOMETRY_DISPLAY_SOLID);
  geom->setAttributes(atts);
  geom->display();

  // show target body's contact //

  string cgeom_name = "contact";
  cgeom_name = cgeom_name + '[' + cname + ']' + "points";

  int num_cverts = cids.size();
  PmVector3 *cverts = new PmVector3[num_cverts];
  float *csrads = new float[num_cverts];

  for (unsigned int i = 0; i < cids.size(); i++) {
    cverts[i] = ccoords[cids[i]];
    csrads[i] = crads[cids[i]];
    }

  PmGraphicsSphere *cgeom = new PmGraphicsSphere(cgeom_name, num_cverts, cverts, csrads);
  atts.setColor(color2);
  //atts.setDisplayType(PM_GEOMETRY_DISPLAY_SOLID);
  atts.setDisplayType(PM_GEOMETRY_DISPLAY_LINE);
  cgeom->setAttributes(atts);
  cgeom->display();
  /*
  PmGraphicsPoint *geom = new PmGraphicsPoint(geom_name, num_verts, verts);
  atts.setMarker(false);
  atts.setScale(2.0);
  atts.setColor(color);
  geom->setAttributes(atts);
  geom->display();
  */

  delete[] verts;
  delete[] srads;
  delete[] cverts;
  delete[] csrads;
  }

//*============================================================*
//*==========       getDomainCoords                  ==========*
//*============================================================*

void 
PmBody::getDomainCoords(const string seq, vector<PmVector3>& coords)
  {
  fprintf (stderr, "\n>>>>>> PmBody::getDomainCoords %s \n", this->name.c_str());
  PmMolecule *mol = dynamic_cast<PmMolecule*>(pobj);

  if (!mol) {
    return;
    }

  PmAtomFilter filter;
  PmVector3 com, disp;
  PmMatrix3x3 mat;
  PmRigidSimResults simr;

  //===== get coords =====//

  coords.clear();
  filter.names.push_back("CA");
  mol->getAtomCoords(seq, filter, coords);
  fprintf (stderr, ">>> num coords=%d \n", coords.size()); 

  //===== transform coords =====//

  this->getSimResult(simr);
  this->getCenterOfMass(com);
  mat = simr.rot_mat;
  disp = simr.disp;

  for (unsigned int i = 0; i < coords.size(); i++) {
    coords[i] = (mat*(coords[i]-com) + disp + com);
    }
  }

//*============================================================*
//*==========       getCurrentCoordinats             ==========*
//*============================================================*
// get the current coordinates for a body.

void
PmBody::getCurrentCoordinats(vector<PmVector3>& cur_coords)
  {
  vector<PmVector3> coords; 
  PmRigidSimResults simr;
  PmVector3 com, disp;
  PmMatrix3x3 mat;

  // get simulation result //

  this->getSimResult(simr);

  // get coordinates from body's physical object //

  pobj->getCoordinates(coords);

  // transform coords //

  cur_coords.clear();
  this->getCenterOfMass(com);
  mat = simr.rot_mat;
  disp = simr.disp;

  for (unsigned int i = 0; i < coords.size(); i++) {
    cur_coords.push_back(mat*(coords[i]-com) + disp + com);
    }
  }

//*============================================================*
//*==========              getDampingFactor          ==========*
//*============================================================*
// get damping factors. 

bool 
PmBody::getDampingFactor(float& factor, float& length)
  {
  if (!has_damping) {
    factor = 0.0;
    length = 0.0;
    return false;
    }

  factor = damping_factor;
  length = damping_length;
  return true;
  }

//*============================================================*
//*==========              setDampingFactor          ==========*
//*============================================================*
// set damping factors. 

void 
PmBody::setDampingFactor(const float val) 
  { 
  has_damping = true;
  damping_factor = val;
  vector<float> dims;

  // get damping length from char length of object //

  if (pobj) {
    pobj->getDimensions(dims);

    if (dims.size()) {
      damping_length = dims[0];
      }
    }
  else {
    damping_length = 1.0;
    }
  }

//*============================================================*
//*==========              getCenterOfMass           ==========*
//*============================================================*
// get center of mass.

void
PmBody::getCenterOfMass(PmVector3& com)
  {
  if (!mass_properties.set) {
    PmMassProperties props;
    getMassProps(props);
    }

  com = mass_properties.com;
  }

//*============================================================*
//*==========              getMassProps              ==========*
//*============================================================*
// get mass properties.

void
PmBody::getMassProps(PmMassProperties& props)
  {
  #define ndbg_PmBody_getMassProps
  #ifdef dbg_PmBody_getMassProps 
  fprintf (stderr, ">>>>>> PmBody::getMassProps:  [%s] \n", name.c_str());
  #endif

  if (pobj) {
    pobj->getMassProps (mass_properties);
    }
  else {
    fprintf (stderr, "   **** WARNING: body \"%s\" does not have a pobj \n", 
             name.c_str());
    return;
    }

  props = mass_properties;

  #ifdef dbg_PmBody_getMassProps 
  fprintf (stderr, "   >>> mass [%g] \n", mass_properties.mass);
  fprintf (stderr, "   >>> center of mass (%g %g %g) \n",
           mass_properties.com[0],
           mass_properties.com[1],
           mass_properties.com[2]);

  fprintf (stderr, "   >>> inertia (%g %g %g) \n",
           mass_properties.inertia(0,0),
           mass_properties.inertia(0,1),
           mass_properties.inertia(0,2));

  fprintf (stderr, "               (%g %g %g) \n",
           mass_properties.inertia(1,0),
           mass_properties.inertia(1,1),
           mass_properties.inertia(1,2));

  fprintf (stderr, "               (%g %g %g) \n",
           mass_properties.inertia(2,0),
           mass_properties.inertia(2,1),
           mass_properties.inertia(2,2));
  #endif

  mass_properties.set = true;
  }

//*============================================================*
//*==========              addKineticEnergy          ==========*
//*============================================================*
// add the kinetic energy for a body. 

void 
PmBody::addKineticEnergy(bool upd, float ke)
  {
  if (upd) {
    kinetic_energy.push_back(ke);
    }
  }

bool 
PmBody::getKineticEnergy(const int step, float& ke)
  {
  ke = 0.0;

  if ((this->type == PM_BODY_GROUND) || (this->type == PM_BODY_STATIC)) {
    return false;
    }

  if (step >= (int)kinetic_energy.size()) {
    return false;
    }

  ke = kinetic_energy[step];
  return true;
  }

//*============================================================*
//*==========              getPosition               ==========*
//*============================================================*
// get the current position for a body.

void 
PmBody::getPosition(PmVector3& pos)
  {
  PmVector3 com;
  this->getCenterOfMass(com);
  pos = com + current_state.displacement;
  }

//*============================================================*
//*==========              addSimResults             ==========*
//*============================================================*
// add the simulation results for a body. note that we still 
// need to update the transformation for a body to update the 
// geometry of any potentials defined for the body.
// -------------------------------------------------------------
// upd: if true then store results for a body
// disp: body displacement

void 
PmBody::addSimResults(bool upd, float time, PmVector3& disp, PmVector3& rot, 
                      PmVector3& vel, PmVector3& avel, PmMatrix3x3& mat, 
                      PmQuaternion& q) 
  {
  /*
  fprintf (stderr, ">>>>>> PmBody::addSimResults [%s] \n", name.c_str());
  fprintf (stderr, "   >>> step = %d  \n", simulation_results.size()); 
  fprintf (stderr, "   >>> disp %g %g %g \n", disp[0], disp[1], disp[2]); 
  fprintf (stderr, "   >>> mat  = %f %f %f \n", mat(0,0), mat(0,1), mat(0,2));
  fprintf (stderr, "              %f %f %f \n", mat(1,0), mat(1,1), mat(1,2));
  fprintf (stderr, "              %f %f %f \n", mat(2,0), mat(2,1), mat(2,2));
  */

  // store results into object //

  if (upd) {
    PmRigidSimResults res;
    res.disp = disp;
    res.rot = rot;
    res.vel = vel;
    res.avel = avel;
    res.rot_mat = mat;
    simulation_results.push_back(res);
    /*
    fprintf (stderr, "\n>>>>>> PmBody::addSimResults [%s] \n", name.c_str());
    fprintf (stderr, "   >>> time = %g \n", time); 
    fprintf (stderr, "   >>> sim res size[%d] \n", simulation_results.size());
    */
    }

  // set xform that will be used for the associated //
  // physical objects graphics xform.               //

  PmVector3 com;
  xform.set = true;
  getCenterOfMass (com);
  xform.center = com;
  xform.translation = disp;

  // convert quaternion into rotation matrix //
  //pm_MathQuaternionToMatrix(q, xform.matrix);
  xform.setRotation(mat);

  if (!pobj) {
    return;
    }

  // set the current xform for the associated physical object //

  pobj->setXform (xform);

  // set current xform for any potentials defined //

  for (unsigned int i = 0; i < potential_geometries.size(); i++) {
    potential_geometries[i]->setXform(xform);
    }

  // set the current xform for xform physical objects //

  for (unsigned int i = 0; i < xform_pobjs.size(); i++) {
    xform_pobjs[i]->setXform (xform);
    }
  }

//*============================================================*
//*==========              getSimResult              ==========*
//*============================================================*
// get the simulation results for a body for the given time step
// <step>. This is used for simulation results output.

void
PmBody::getSimResult(const int step, PmRigidSimResults& result)
  {
  /*
  fprintf (stderr, " body %s  results size=%d step=%d \n", this->name.c_str(),
           simulation_results.size(), step);
  */

  if (type == PM_BODY_STATIC) {
    result.disp.set(0,0,0);
    result.rot.set(0,0,0); 
    result.rot_mat.diag(1.0);
    result.vel.set(0,0,0);
    result.avel.set(0,0,0);
    return;
    }

  if (step >= (int)simulation_results.size()) {
    pm_ErrorWarnReport (PM, "step %d too large for getting results for body \"%s\" .", 
                        "*", step, name.c_str());
    return;
    }

  if (!pobj) {
    result.disp.set(0,0,0);
    result.rot.set(0,0,0); 
    result.rot_mat.diag(1.0);
    result.vel.set(0,0,0);
    result.avel.set(0,0,0);
    return;
    }

  result.disp = simulation_results[step].disp;
  result.rot = simulation_results[step].rot;
  result.rot_mat = simulation_results[step].rot_mat;
  result.vel = simulation_results[step].vel;
  result.avel = simulation_results[step].avel;
  }

//*============================================================*
//*==========              getSimResult              ==========*
//*============================================================*
// get the current simulation result for a body.

void
PmBody::getSimResult(PmRigidSimResults& result)
  {
  if (simulation_results.size() == 0) {
    return;
    }

  if (!pobj) {
    return;
    }

  result = simulation_results.back();
  }

//*============================================================*
//*==========              setSimResults             ==========*
//*============================================================*
// set the simulation results for a body for the given time step
// <step>. This is used for replaying a simulation.

void
PmBody::setSimResults(const int step)
  {
  /*
  fprintf (stderr, ">>>>>> PmBody::setSimResults [%s] \n", name.c_str());
  fprintf (stderr, "   >>> step[%d] \n", step);
  fprintf (stderr, "   >>> sim res size[%d] \n", simulation_results.size());
  */

  if (step >= (int)simulation_results.size()) {
    return;
    }

  if (!pobj) {
    return;
    }

  PmVector3 disp = simulation_results[step].disp;
  PmMatrix3x3 mat = simulation_results[step].rot_mat;

  /*
  fprintf (stderr, "   >>> mat  = %f %f %f \n", mat(0,0), mat(0,1), mat(0,2));
  fprintf (stderr, "              %f %f %f \n", mat(1,0), mat(1,1), mat(1,2));
  fprintf (stderr, "              %f %f %f \n", mat(2,0), mat(2,1), mat(2,2));
  */

  PmXform pxform;
  PmVector3 com;
  pxform.set = true;
  getCenterOfMass (com);
  pxform.center = com;
  pxform.translation = disp;
  pxform.setRotation(mat);

  // set the current xform for the associated physical object //

  pobj->setXform(pxform);
  }

//*============================================================*
//*==========              get/setState              ==========*
//*============================================================*

void 
PmBody::getState(PmRigidBodyState& state)
  {
  state = current_state;
  }

void 
PmBody::setState(const PmRigidBodyState& state) 
  {
  last_state = current_state;
  current_state = state;
  }

//*============================================================*
//*==========              getVelocity               ==========*
//*============================================================*
// get the velocity at the center of mass.

void 
PmBody::getVelocity(PmVector3& vel)
  {
  vel = current_state.velocity;
  }

//*============================================================*
//*==========              getVelocity               ==========*
//*============================================================*
// get the velocity at the given point.   

void 
PmBody::getVelocity(const PmVector3& point, PmVector3& vel) 
  {
  PmVector3 r, com;
  this->getCenterOfMass(com);
  r = point - com + current_state.displacement;
  vel = current_state.velocity + current_state.angular_velocity.cross(r);
  }

}



