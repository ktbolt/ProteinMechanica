
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
//              f o r c e   c o m m a n d s                  //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========             pm_CmdForce                ==========*
//*============================================================*
// process force command.

void
pm_CmdForce (PmCmdDataList& dlist)
  {

  PmCmdData data;
  string dv, name, dname, res, type_str, jnt_name;
  PmForceType ftype;
  PmVector3 point, dir;
  float scale, ramp, mag;
  PmCmdData pdata;
  PmForce *force;
  PmVector3 bpos, color;
  bool use_pca, interactive;

  dlist.getNext (data);

  //===== create a new force and add it to the system =====//

  if (data.name == "create") {
    dlist.getString("name", name);
    dlist.getString("type", type_str);

    if (type_str != "") {
      PmForce::convForceType(type_str, ftype);

      if (ftype == PM_FORCE_UNKNOWN) {
        pm_ErrorReport (PM, "unknown force type \"%s\".", "*", type_str.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no force type specified.", "*");
      return;
      }

    // create a force of a particular type and add it to the system //

    force = PmForce::create(name, ftype);
    pmSystem.addForce(force);
    dlist.getBoolean("use_pca", use_pca);

    if (ftype == PM_FORCE_EXPLICIT) {
      bool global_frame, torque;
      PmExplicitForce *eforce = dynamic_cast<PmExplicitForce*>(force);

      // Note [6Jan2010] change to use global or world terms //

      if (dlist.getBoolean("global_frame", global_frame)) {
        eforce->setGlobalFrame(global_frame);
        }

      if (dlist.getBoolean("world_frame", global_frame)) {
        eforce->setGlobalFrame(global_frame);
        }

      if (dlist.getBoolean("torque", torque)) {
        eforce->setTorque(torque);
        }

      if (dlist.getBoolean("interactive", interactive)) {
        eforce->setInteractive(interactive);
        }

      if (!dlist.getFloat("magnitude", mag)) {
        mag = 1.0;
        }

      if (!torque && !use_pca) {
        if (!dlist.getData("point", pdata)) {
          pm_ErrorReport (PM, "force point not specified.", "*");
          return;
          }

        if (!pm_CmdProcPosition(pdata, point)) {
          pm_ErrorReport (PM, "couldn't process force point.", "*");
          return;
          }
        }

      // set the torque on a joint axis //

      if (dlist.getString("joint", jnt_name)) {
        PmJoint *joint;
        int axis_id;
        PmJointType jtype;
        pmSystem.getJoint(jnt_name, &joint);

        if (!joint) {
          pm_ErrorReport (PM, "no joint named \"%s\".", "*", jnt_name.c_str());
          return;
          }

        jtype = joint->getType();

        if ((jtype == PM_JOINT_FREE) || (jtype == PM_JOINT_WELD)) {
          pm_ErrorReport (PM, "joint \"%s\" does not have an axis.", "*", 
                          jnt_name.c_str());
          return;
          }

        if ((jtype != PM_JOINT_HINGE) && !dlist.getInt("axis", axis_id)) {
          pm_ErrorReport (PM, "no axis id given.", "*");
          return;
          }

        if (jtype == PM_JOINT_HINGE) {
          PmHingeJoint *hjoint = dynamic_cast<PmHingeJoint*>(joint);
          hjoint->getAxis(dir);
          }
        else if (jtype == PM_JOINT_BALL) {
          if ((axis_id < 0) || (axis_id > 3)) {
            pm_ErrorReport (PM, "axis id must be between 1-3.", "*");
            return;
            }

          PmBallJoint *bjoint = dynamic_cast<PmBallJoint*>(joint);
          bjoint->getAxis(axis_id, dir);
          }
        else if (jtype == PM_JOINT_UNIVERSAL) {
          if ((axis_id < 0) || (axis_id > 2)) {
            pm_ErrorReport (PM, "axis id must be 1 or 2.", "*");
            return;
            }

          PmUniversalJoint *ujoint = dynamic_cast<PmUniversalJoint*>(joint);
          ujoint->getAxis(axis_id, dir);
          }

        dir = dir*mag;
        }

      //===== set axis using pca =====//

      else if (dlist.getBoolean("use_pca", use_pca) && use_pca) {
        string dname, desc;
        PmMolecule *domain;
        PmAtomFilter filter;
        PmPcaResults pca;
        vector<PmVector3> coords;
        float wd[3];
        PmExtent proj;
        int pca_axis;

        if (dlist.getString("pca_domain", dname)) {
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
        //fprintf (stderr, ">>> pca number of coords=%d \n", coords.size());
        //fprintf (stderr, ">>> pca wd=%f %f %f\n", wd[0], wd[1], wd[2]); 

        if (!dlist.getInt("direction", pca_axis)) {
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

        eforce->getGlobalFrame(global_frame);

        if (global_frame) {
          point = pca.com;
          }
        else {
          point.set(0,0,0);
          }

        dir = mag*dir;

        if (interactive) {
          PmVector3 com, axis1, axis2, axis3;
          float val, s=1.2;
          bool use_sphere;
          //axis1 = pca.s1*pca.axis1;
          //axis2 = pca.s2*pca.axis2;
          //axis3 = pca.s3*pca.axis3;
          axis1 = s*wd[0]*pca.axis1;
          axis2 = s*wd[0]*pca.axis2;
          axis3 = s*wd[0]*pca.axis3;
          eforce->setAxes(axis1, axis2, axis3);

          if (dlist.getFloat("strength", val)) {
            eforce->setStrength(val);
            }

          if (dlist.getBoolean("use_sphere", use_sphere)) {
            eforce->setUsePickingSphere(use_sphere);
            eforce->setPickingSphereRadius(1.1*wd[0]);
            }
          }
        }

      //===== use vector =====//

      else if (!dlist.getVector("direction", dir)) {
        string vname;
        PmGeometry *geom;
        PmGeometryVector *vgeom;
        PmVector3 pts[2];

        if (!dlist.getString("direction", vname)) {
          pm_ErrorReport (PM, "force direction not specified.", "*");
          return;
          }

        pmSystem.getGeometry(vname, &geom);
        vgeom = dynamic_cast<PmGeometryVector*>(geom);

        if (!vgeom) {
          pm_ErrorReport (PM, "no vector named \"%s\".", "*", vname.c_str());
          return;
          }

        vgeom->getCurrVector(dir);
        dir.normalize();
        }

      dir = mag*dir;
      eforce->setPoint(point);
      eforce->setDirection(dir);
      }

    else if (ftype == PM_FORCE_RANDOM) {
      float val;
      int ival;
      PmRandomForce *rforce = dynamic_cast<PmRandomForce*>(force);

      if (dlist.getFloat("mean", val)) {
        rforce->setMean(val);
        }

      if (dlist.getFloat("sd", val)) {
        rforce->setStandardDev(val);
        }

      if (dlist.getInt("seed", ival)) {
        rforce->setSeed(ival);
        }
      }

    if (dlist.getFloat("scale", scale)) {
      force->setScale(scale);
      }
    }

  else {
    pmSystem.getForce(data.name, &force);
    }

  //fprintf (stderr, "\n###### force=%x \n", force);

  if (!force) {
    pm_ErrorReport (PM, "no force named \"%s\".", "*", data.name.c_str());
    return;
    }

  // process force parameters //

  bool show, show_set = false;

  while (dlist.getNext(data, dv)) {
    if (data.name == "color") {
      if (data.getVector(color)) {
        force->setColor(color);
        }
      }

   else if (data.name == "ramp") {
     ramp = data.getFloat();
     force->setRampParams(ramp);
     }

    else if (data.name == "show") {
      show = data.getBoolean();
      show_set = true;
      }
    }

  if (show_set) {
    force->setShow(show);
    force->display(point);
    }
  }

