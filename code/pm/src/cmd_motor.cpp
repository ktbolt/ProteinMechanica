
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
//              m o t o r   c o m m a n d s                  //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========             pm_CmdMotor                ==========*
//*============================================================*
// process motor command.

void
pm_CmdMotor (PmCmdDataList& dlist)
  {

  PmCmdData data;
  string dv, name, jname, type_str;
  PmMotorType mtype;
  PmVector3 point, dir;
  PmCmdData pdata;
  PmMotor *motor;
  PmVector3 bpos;
  PmJoint *joint;

  dlist.getNext (data);

  // create a new motor and add it to the system //

  if (data.name == "create") {
    dlist.getString("name", name);
    dlist.getString("type", type_str);

    if (type_str != "") {
      PmMotor::convMotorType(type_str, mtype);

      if (mtype == PM_MOTOR_UNKNOWN) {
        pm_ErrorReport (PM, "unknown motor type \"%s\".", "*", type_str.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no motor type specified.", "*");
      return;
      }

    if (dlist.getString("joint", jname)) {
      pmSystem.getJoint(jname, &joint);

      if (!joint) {
        pm_ErrorReport (PM, "not joint named \"%s\".", "*", jname.c_str());
        return;
        }
      }
    else {
      pm_ErrorReport (PM, "no motor type specified.", "*");
      return;
      }

    // create a motor of a particular type and add it to the system //

    motor = PmMotor::create(name, mtype);
    pmSystem.addMotor(motor);

    if (mtype == PM_MOTOR_ANGULAR) {
      vector<string> list;
      int num_axes, axes[3];
      float val;

      if (!dlist.getStringList("axes", list)) {
        pm_ErrorReport (PM, "no joint axes given.", "*");
        return;
        }

      num_axes = list.size();

      if (num_axes > 3) {
        pm_ErrorReport (PM, "too many axes ids given.", "*");
        return;
        }

      for (int i = 0; i < num_axes; i++) {
        axes[i] = atoi(list[i].c_str());

        if ((axes[i] < 1) || (axes[i] > 3)) { 
          pm_ErrorReport (PM, "axis id not in range 1-3.", "*");
          return;
          }
        }

      PmAngularMotor *amotor = dynamic_cast<PmAngularMotor*>(motor);
      amotor->setJoint(joint);
      amotor->setAxes(num_axes, axes);
      joint->addMotor(motor);
      fprintf (stderr, "    >>> angular motor \"%s\" created. \n", name.c_str());

      // set maximum velocities //

      if (dlist.getStringList("max_vel", list)) {
        if ((int)list.size() != num_axes) {
          pm_ErrorReport (PM, "number of velocities not equal to number of axes.", "*");
          return;
          }

        for (int i = 0; i < num_axes; i++) {
          val = atof(list[i].c_str());
          amotor->setMaxVelocity(i+1, val);
          }
        }

      // set maximum forces //

      if (dlist.getStringList("max_force", list)) {
        if ((int)list.size() != num_axes) {
          pm_ErrorReport (PM, "number of forces not equal to number of axes.", "*");
          return;
          }

        for (int i = 0; i < num_axes; i++) {
          val = atof(list[i].c_str());
          amotor->setMaxForce(i+1, val);
          }
        }

      // set maximum angles //

      if (dlist.getStringList("max_angle", list)) {
        if ((int)list.size() != num_axes) {
          pm_ErrorReport (PM, "number of angles not equal to number of axes.", "*");
          return;
          }

        for (int i = 0; i < num_axes; i++) {
          val = atof(list[i].c_str());
          amotor->setMaxAngle(i+1, val);
          //fprintf (stderr, "    >>> set max angle = %f \n", val);
          }
        }

      // get angle increments //

      if (dlist.getStringList("angle_inc", list)) {
        if ((int)list.size() != num_axes) {
          pm_ErrorReport (PM, "number of angle inc not equal to number of axes.", "*");
          return;
          }

        for (int i = 0; i < num_axes; i++) {
          val = atof(list[i].c_str());
          amotor->setAngleInc(i+1, val);
          }
        }
      }
    }

  else {
    pmSystem.getMotor(data.name, &motor);
    }

  if (!motor) {
    pm_ErrorReport (PM, "no motor named \"%s\".", "*", data.name.c_str());
    return;
    }

  // process motor parameters //

  bool show_set = true;

  while (dlist.getNext(data)) {
    }

  if (show_set) {
    //body->display (show);
    }
  }

