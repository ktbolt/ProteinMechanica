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
//              g r a p h i c s   c o m m a n d s            //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========            pm_CmdGraphics              ==========*
//*============================================================*
// process graphics commands.                                 

void
pm_CmdGraphics (PmCmdDataList& dlist)
  {
  PmCmdData data;
  string dv, name, type, obj, format;

  // get graphics interface //

  PmGraphicsInterface *graphics = pmSystem.getGraphicsInterace();

  if (!graphics) {
    return;
    }

  // get next item //

  while (dlist.getNext(data, dv)) { 

    /*===== background =====*/

    if (data.name == "background") {
      PmVector3 color;

      if (dlist.getVector("color", color)) {
        graphics->setWindowColor(color);
        }
      }

    /*===== center =====*/

    else if (data.name == "center") {
      PmVector3 center;
      PmMolecule *domain;
      dlist.getNext(data, dv);

      if (data.name == "pick") {
        graphics->getCurrentPick(center);
        graphics->setCenter(center);
        }
      else if (data.name == "geometry") {
        //gr_GeomGet (dv, &geom);
        //gr_GeomCenterGet (geom, center);
        //gr_SceneXformCenterSet (man.scene, center);
        }

      else if (data.name == "point") {
        if (!data.getVector(center)) { 
          string pstr;
          PmCmdData pdata;
          dlist.getString("point", pstr);
          dlist.getData("point", pdata);

          if (!pm_CmdProcPosition(pdata, center)) {
            pm_ErrorReport (PM, "couldn't process point data \"%s\".", "*",
                            pstr.c_str());
           return;
           }
         }

        graphics->setCenter(center);
        }

      else if (data.name == "domain") {
        string desc;
        PmAtomFilter filter;
        pmSystem.getDomain (dv, &domain);

        if (!domain) {
          pm_ErrorReport (PM, "no domain named \"%s\".", "*", dv.c_str());
          return;
          }

        domain->getCenter(desc, filter, center);
        fprintf (stderr, "    >>> center set to (%f %f %f) \n", center[0], center[1], 
                 center[2]); 
        graphics->setCenter(center);
        }
      }

    /*===== grid =====*/

    else if (data.name == "grid") {
      PmCmdData pdata;
      string pstr, name;
      PmVector3 center, spacing, number;
      PmVector3 color;
      bool show;
      PmGraphicsAttributes atts;
      float msize;

      dlist.getData("center", pdata);
      dlist.getString("point", pstr);

      if (!pm_CmdProcPosition(pdata, center)) {
        pm_ErrorReport (PM, "couldn't process point data \"%s\".", "*", pstr.c_str());
        return;
        }

      if (!dlist.getString("name", name)) {
        pm_ErrorReport (PM, "no grid name given.", "*");
        return;
        }

      if (dlist.getVector("color", color)) {
        atts.setColor(color);
        }

      if (dlist.getFloat("msize", msize)) {
        atts.setScale(msize);
        }

      if (!dlist.getVector("number", number)) {
        pm_ErrorReport (PM, "no grid number given.", "*");
        return;
        }

      if (!dlist.getVector("spacing", spacing)) {
        pm_ErrorReport (PM, "no grid spacing given.", "*");
        return;
        }

      dlist.getBoolean("show", show);

      PmGraphicsGrid *grid = new PmGraphicsGrid(name, center, number, spacing);

      grid->setAttributes(atts);
      grid->display(show);
      }


    /*===== plane =====*/

    else if (data.name == "plane") {
      PmCmdData pdata;
      string pstr, name;
      PmVector3 center, normal; 
      float spacing;
      int number;
      PmVector3 color;
      bool show;
      PmGraphicsAttributes atts;
      float msize;

      dlist.getData("center", pdata);
      dlist.getString("point", pstr);

      if (!pm_CmdProcPosition(pdata, center)) {
        pm_ErrorReport (PM, "couldn't process point data \"%s\".", "*", pstr.c_str());
        return;
        }

      if (!dlist.getString("name", name)) {
        pm_ErrorReport (PM, "no grid name given.", "*");
        return;
        }

      if (!dlist.getVector("normal", normal)) {
        pm_ErrorReport (PM, "no plane normal given.", "*");
        return;
        }

      if (dlist.getVector("color", color)) {
        atts.setColor(color);
        }

      if (dlist.getFloat("msize", msize)) {
        atts.setScale(msize);
        }

      if (!dlist.getInt("number", number)) {
        pm_ErrorReport (PM, "no plane number given.", "*");
        return;
        }

      if (!dlist.getFloat("spacing", spacing)) {
        pm_ErrorReport (PM, "no plane spacing given.", "*");
        return;
        }

      dlist.getBoolean("show", show);

      PmGraphicsPlane *plane = new PmGraphicsPlane(name, center, normal, number, spacing);

      plane->setAttributes(atts);
      plane->display(show);
      }
  
    /*===== pick =====*/

    else if (data.name == "pick") {
      /*
      pm_CmdDataValGet (&dlist, "obj", &obj);
      pm_CmdDataValGet (&dlist, "type", &type);
      pm_CmdDataIntGet (&dlist, "active", &active);

      while (pm_CmdDataNext(&dlist, &dn, &dv)) { 
        if (!strcmp(dn, "size")) {
          man.graphics.pick_size = atof (dv); 
          }
        else if (!strcmp(dn, "width")) {
          man.graphics.pick_width = atof (dv); 
          }
        else if (!strcmp(dn, "color")) {
          pm_CmdDataFvecGet (&dlist, "color", man.graphics.pick_color);
          }
        else if (!strcmp(dn, "radius")) {
          man.graphics.pick_radius = atof (dv); 
          }
        }
      */
      }

    /*===== point =====*/

    else if (data.name == "point") {
      /*
      pm_CmdDataValGet (&dlist, "name", &name);
      pm_CmdDataValGet (&dlist, "type", &type);
      pm_CmdDataFvecGet (&dlist, "pos", pos);
      pm_CmdDataFloatGet (&dlist, "size", &size);
      pm_CmdDataFvecGet (&dlist, "color", color);
      pm_GrPointShow (man.win, man.scene, name, type, pos, size, color); 
      */
      }

    /*===== reset =====*/

    else if (data.name == "reset") {
      //gr_SceneReset (man.win, man.scene, GR_RENDER_MODE_SURFACE);
      }

    /*===== record =====*/

    else if (data.name == "record") {
      //pm_CmdDataNext (&dlist, &dn, &dv);
      //pm_GrWinRecord (dn);

      if (!dlist.getString("name", name)) {
        return;
        }

      if (!dlist.getString("format", format)) {
        format = "jpeg";
        }

      graphics->recordWindow(name, format);
      fprintf (stderr, "    >>> graphics record on. \n");
      }

    /*===== rotate =====*/

    else if (data.name == "rotate") {
      PmVector3 rot;
      data.getVector(rot);
      graphics->rotateScene(rot);
      }

    /*===== scale =====*/

    else if (data.name == "scale") {
      float scale = data.getFloat();
      graphics->scaleScene(scale);
      }

    /*===== translate =====*/

    else if (data.name == "translate") {
      PmVector3 trans;
      data.getVector(trans);
      graphics->translateScene(trans);
      }


    /*===== xform =====*/

    else if (data.name == "xform") {
      PmXform xform;
      graphics->getSceneXform(xform);
      fprintf (stderr, "\n------ graphics scene transformation ------\n");
      fprintf (stderr, ">>> translation = %f %f %f \n", xform.translation[0], 
                                                        xform.translation[1], 
                                                        xform.translation[2]);
      fprintf (stderr, ">>> rotation    = %f %f %f \n", xform.angles[0], 
                                                        xform.angles[1], 
                                                        xform.angles[2]);
      }

    /*===== write =====*/

    else if (data.name == "write") {
      if (!dlist.getString("name", name)) {
        return;
        }

      if (!dlist.getString("format", format)) {
        format = "jpeg";
        }

      graphics->writeWindow(name, format);
      }
    }
  }
