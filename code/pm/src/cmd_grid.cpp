
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
//*==========              pm_CmdGridSlice           ==========*
//*============================================================*
// process grid slice command.

void
pm_CmdGridSlice(PmGrid *grid, PmCmdDataList& dlist)
  {
  PmCmdData data;
  float vmax, vmin;
  int id = -1, dim = -1, index = -1;

  grid->getDataRange(vmin, vmax);

  while (dlist.getNext(data)) {
    if (data.name == "id") {
      id = data.getInt();
      }

    else if (data.name == "i") {
      index = data.getInt();
      dim = 1;
      }

    else if (data.name == "j") {
      index = data.getInt();
      dim = 2;
      }

    else if (data.name == "k") {
      index = data.getInt();
      dim = 3;
      }

    else if (data.name == "vmax") {
      vmax = data.getFloat();
      }

    else if (data.name == "vmin") {
      vmin = data.getFloat();
      }
    }

  grid->displaySlice(id, dim, index, vmin, vmax);
  }

//*============================================================*
//*==========        pm_CmdGridIsosurface            ==========*
//*============================================================*
// process grid isosurface command.

void
pm_CmdGridIsosurface(PmGrid *grid, PmCmdDataList& dlist)
  {
  float level = 0.0;
  PmCmdData data;
  bool seeded = false;
  bool write = false;
  string dv, fname;
  PmVector3 point, color;

  PmDbInterface *db = NULL; 
  PmDbPmInterface *pmdb = NULL; 
  PmDbInterfaceSelect db_sel;
  PmDbType format;

  while (dlist.getNext(data, dv)) {
    if (data.name == "color") {
      data.getVector(color);
      grid->setColor(color);
      }
    else if (data.name == "level") {
      level = data.getFloat();
      }

    else if (data.name == "point") {
      data.getVector(point);
      seeded = true;
      }

    else if (data.name == "file") {
      data.getString(fname);
      db = db_sel.create ("tmp", PM_DB_PM);
      db->open(fname, PM_DB_MODE_WRITE, "");
      pmdb = dynamic_cast<PmDbPmInterface*>(db); 
      }

    else if (data.name== "shading") {
      PmGeometryShadingType stype;
      PmGraphicsGeometry::convShadingType (dv, stype);

      if (stype == PM_GEOMETRY_SHADING_UNKNOWN) {
        pm_ErrorReport (PM, "unknown shading type \"%s\".", "*", dv.c_str());
        }
      else {
        grid->setShadingType(stype);
        }
      }

    else if (data.name== "display") {
      PmGeometryDisplayType dtype;
      PmGraphicsGeometry::convDisplayType (dv, dtype);

      if (dtype == PM_GEOMETRY_DISPLAY_UNKNOWN) {
        pm_ErrorReport (PM, "unknown display type \"%s\".", "*", dv.c_str());
        pmSystem.setCmdError(true);
        }
      else {
        grid->setDisplayType(dtype);
        }
      }
    }

  grid->displayIsosurface(level, seeded, point, pmdb);

  if (db) {
    db->close();
    }
  }

//*============================================================*
//*==========        pm_CmdGridVolume                ==========*
//*============================================================*
// process grid volume command.

void
pm_CmdGridVolume(PmGrid *grid, PmCmdDataList& dlist)
  {
  float vmax = 0.0, vmin = 0.0;
  int num_samples =  1;
  PmCmdData data;

  grid->getDataRange(vmin, vmax);

  while (dlist.getNext(data)) {
    if (data.name == "samples") {
      num_samples = data.getInt();
      }
    else if (data.name == "vmax") {
      vmax = data.getFloat();
      }
    else if (data.name == "vmin") {
      vmin = data.getFloat();
      }
    }

  grid->displayVolume(vmin, vmax, num_samples);
  }

//*============================================================*
//*==========              pm_CmdGridPoints          ==========*
//*============================================================*
// process grid points command.

void
pm_CmdGridPoints(PmGrid *grid, PmCmdDataList& dlist)
  {
  float vmax = 0.0, vmin = 0.0;
  PmCmdData data;
  string fname;
  PmDbInterface *db = NULL;
  PmDbPmInterface *pmdb = NULL;
  PmDbInterfaceSelect db_sel;
  PmDbType format;

  while (dlist.getNext(data)) {
    if (data.name == "vmax") {
      vmax = data.getFloat();
      }

    else if (data.name == "vmin") {
      vmin = data.getFloat();
      }

    else if (data.name == "file") {
      data.getString(fname);
      db = db_sel.create ("tmp", PM_DB_PM);
      db->open(fname, PM_DB_MODE_WRITE, "");
      pmdb = dynamic_cast<PmDbPmInterface*>(db);
      }
    }

  grid->displayPoints(vmin, vmax, pmdb);

  if (db) {
    db->close();
    }

  }

//*============================================================*
//*==========              pm_CmdGrid                ==========*
//*============================================================*
// process grid command.                                 

void
pm_CmdGrid (PmCmdDataList& dlist)
  {
  string dv, dbname, name, type; 
  PmGrid *grid;
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
    grid = db->getGrid();

    if (!grid) {
      pm_ErrorReport (PM, "could not read grid.", "*");
      return;
      }

    //grid->setName(name);

    // add grid to system //

    pmSystem.addGrid(grid);
    grid->getExtent(extent);
   
    fprintf (stderr, "    >>> extent  min (%g %g %g)  max (%g %g %g) \n", 
             extent.min[0], extent.min[1], extent.min[2],
             extent.max[0], extent.max[1], extent.max[2]);

    // set extent //

    pmSystem.setExtent (extent);
    return;
    }

  else if (data.name == "set") {
    dlist.getString("name", name);
    //pmSystem.getMolecule(name, &mol);
    //pmSystem.setCurrentMolecule(mol);
    return;
    }

  pmSystem.getGrid(data.name, &grid);

  if (!grid) {
    pm_ErrorReport (PM, "no grid named \"%s\".", "*", data.name.c_str());
    return;
    }

  while (dlist.getNext(data)) {

    // ===== points ===== //

    if (data.name == "points") {
      pm_CmdGridPoints(grid, dlist);
      }

    // ===== xform ===== //

    else if (data.name == "xform") {
      }

    // ===== isosurface ===== //

    else if (data.name == "isosurface") {
      pm_CmdGridIsosurface(grid, dlist);
      }

    // ===== slice ===== //

    else if (data.name == "slice") {
      pm_CmdGridSlice(grid, dlist);
      }

    // ===== volume ===== //

    else if (data.name == "volume") {
      pm_CmdGridVolume(grid, dlist);
      }
    }
  }

