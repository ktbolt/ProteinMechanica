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
//*==========              pm_CmdParticle            ==========*
//*============================================================*
// process particle command.                                 

void
pm_CmdParticle (PmCmdDataList& dlist)
  {
  string dv, dbname, name, type; 
  PmParticle *particle;
  PmExtent extent;
  PmCmdData data;

  dlist.getNext (data);

  if (data.name == "read") {
    string fname, format_str;
    PmDbInterfaceSelect db_sel;
    PmDbType format;

    dlist.getString ("name", name);
    dlist.getString ("format", format_str);
    dlist.getString ("file", fname);
    format = pm_DbFormatConv (format_str);
    PmDbInterface *db = db_sel.create ("tmp", format);

    if (!db) {
      pm_ErrorReport ("pm> ", "unknown format \"%s\".", "*", format_str.c_str());
      return;
      }

    db->open(fname, PM_DB_MODE_READ, "");
    db->read(name);
    db->close();

    PmDbPmInterface* pmdb = dynamic_cast<PmDbPmInterface*>(db);
    particle = pmdb->getParticle();

    if (!particle) {
      pm_ErrorReport (PM, "could not read particle.", "*");
      return;
      }

    // add particle to system //

    pmSystem.addParticle(particle);
    particle->getExtent(extent);
   
    fprintf (stderr, "    >>> extent  min (%g %g %g)  max (%g %g %g) \n", 
             extent.min[0], extent.min[1], extent.min[2],
             extent.max[0], extent.max[1], extent.max[2]);

    // set extent //

    pmSystem.setExtent (extent);
    return;
    }

  else {
    pmSystem.getParticle(data.name, &particle);
    }

  if (!particle) {
    pm_ErrorReport (PM, "no particle named \"%s\".", "*", data.name.c_str());
    return;
    }

  while (dlist.getNext(data)) {

    if (data.name == "show") {
      bool show = data.getBoolean();
      particle->display(show);
      }

    else if (data.name == "color") {
      PmVector3 color;
      data.getVector(color);
      particle->setColor(color);
      }

    else if (data.name == "spheres") {
      bool spheres = data.getBoolean();
      particle->setDisplaySpheres(spheres);
      }

    else if (data.name == "map_mass") {
      bool val = data.getBoolean();
      particle->setMapMass(val);
      }

    else if (data.name == "xform") {
      }
    }
  }

