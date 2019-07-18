
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
//              m o d e l   c o m m a n d s                  //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========              pm_CmdModel               ==========*
//*============================================================*
// process model command.                                 

void
pm_CmdModel (PmCmdDataList& dlist)
  {

  string dv, fname, name, format_str;
  PmModel *model;
  PmJoint *joint;
  PmDbInterfaceSelect db_sel;
  PmExtent extent;
  PmDbType format;
  PmCmdData data;

  dlist.getNext (data);

  // create a new model and add it to the system //

  if (data.name == "create") {
    dlist.getString ("name", name);
    model = new PmModel(name);
    pmSystem.addModel(model);
    return;
    }

  else if (data.name == "read") {
    dlist.getString ("name", name);
    dlist.getString ("format", format_str);
    dlist.getString ("file", fname);

    format = pm_DbFormatConv (format_str);
    PmDbInterface *db = db_sel.create ("tmp", format);
    db->open(fname, PM_DB_MODE_READ, "");
    db->read(name);
    db->close();
    //model = db->getModel();

    // add model to system
    pmSystem.addModel(model);

    //mol->getExtent(extent);
   
    fprintf (stderr, "    >>> extent  min (%g %g %g)  max (%g %g %g) \n", 
             extent.min[0], extent.min[1], extent.min[2],
             extent.max[0], extent.max[1], extent.max[2]);

    // set extent
    pmSystem.setExtent (extent);
    return;
    }

  else if (data.name == "set") {
    dlist.getString ("name", name);
    pmSystem.getModel(name, &model);
    pmSystem.setCurrentModel (model);
    return;
    }

  pmSystem.getModel(data.name, &model);

  if (!model) {
    pm_ErrorReport (PM, "no model named \"%s\".", "*", data.name.c_str());
    return;
    }

  // process model commands 

  while (dlist.getNext(data)) { 
    if (data.name == "add") {
      dlist.getNext(data); 

      if (data.name == "joint") {
        data.getString (dv);
        pmSystem.getJoint (dv, &joint);

        if (!joint) {
          pm_ErrorReport (PM, "no joint named \"%s\".", "*", dv.c_str());
          return;
          }

        model->addJoint (joint);
        fprintf (stderr, "    >>> joint \"%s\" added to model \n", dv.c_str());
        }

      else if (data.name == "ground") {
        dlist.getNext(data); 

        if (data.name == "body") {
          PmBody *body;
          data.getString (dv);
          pmSystem.getBody (dv, &body);

          if (!body) {
            pm_ErrorReport (PM, "no body named \"%s\".", "*", dv.c_str());
            return;
            }

          if (model->addBody(body)) {
            fprintf (stderr, "    >>> body \"%s\" added to model \n", dv.c_str());
            }
          }
        }

      else if (data.name == "body") {
        PmBody *body;
        data.getString (dv);
        pmSystem.getBody (dv, &body);

        if (!body) {
          pm_ErrorReport (PM, "no body named \"%s\".", "*", dv.c_str());
          return;
          }

        if (model->addBody(body)) { 
          fprintf (stderr, "    >>> body \"%s\" added to model \n", dv.c_str());
          }
        }
      }

    else if (data.name == "bodies") {
      bool show_set = true;
      bool show = true;

      while (dlist.getNext(data)) {
        if (data.name == "msize") {
          float msize = data.getFloat();
          model->setBodyMsize (msize);
          }
        }

      if (show_set) {
        model->displayBodies(show);
        }
      }

    else if (data.name == "joints") {
      bool show_set = true;
      bool show = true;

      while (dlist.getNext(data)) {
        if (data.name == "msize") {
          float msize = data.getFloat();
          model->setJointMsize (msize);
          }
        }

      if (show_set) {
        model->displayJoints(show);
        }
      }
    }
  }

