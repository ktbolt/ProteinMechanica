
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
// * graphics:           g r a p h i c s     A P I             *
//*============================================================*

#include "graphics.h"
#include "graphics_prv.h"

namespace ProteinMechanica {

void pm_CmdProc(char*);

#ifdef USE_GRAPHICS
static GrGeometryDisplayType gr_ren_map[PM_MOLECULE_RENDER_SIZE+1];
#endif

//*============================================================*
//*==========              init                      ==========*
//*============================================================*
// initialize graphics.

void
PmGraphicsInterface::api_init ()
  {

  fprintf (stderr, ">>>>>> PmGraphicsInterface::api_init  this [%x] \n", this);
  fprintf (stderr, "   >>> prv_data [%x] \n", prv_data);

  // set cmd processing function.
  grSystem.setCommandCallback (pm_CmdProc);

  // set mapping from pm to gr symbols
  gr_ren_map[PM_MOLECULE_RENDER_UNKNOWN] = GR_GEOMETRY_DISPLAY_SOLID;
  gr_ren_map[PM_MOLECULE_RENDER_POINT] = GR_GEOMETRY_DISPLAY_POINT;
  gr_ren_map[PM_MOLECULE_RENDER_LINE] = GR_GEOMETRY_DISPLAY_LINE;
  gr_ren_map[PM_MOLECULE_RENDER_SOLID] = GR_GEOMETRY_DISPLAY_SOLID;
  }

//*============================================================*
//*==========            api_create_window           ==========*
//*============================================================*
// create window.             

void
PmGraphicsInterface::api_create_window (char *name, int x, int y, int win_width, 
                                        int win_height)
  { 
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;

#ifdef USE_GRAPHICS

  GrExtent extent(-10.0, 10.0, -10.0, 10.0, -200.0, 200.0);

  // create a window

  int wtype = GR_WINDOW_TYPE_DOUBLE | GR_WINDOW_TYPE_RGB | GR_WINDOW_TYPE_DEPTH |
              GR_WINDOW_TYPE_DIRECT;

  GrWindow *win = new GrWindow (name, x, y, win_width, win_height, wtype);

  win->setExtent (extent);
  win->open ();


  // create a scene

  GrScene *scene = new GrScene (name);
  scene->setWindow (win);
  scene->setViewport (win_width, win_height, extent);
  scene->setExtent (extent);

  // create a light

  GrLight *light = new GrLight ("light 1", GR_LIGHT_DIRECTIONAL);
  scene->addLight (light);


  // attach the scene to the window
  win->setScene (scene);

  //win->processEvents (NULL, false, "test>");

  prvd->window = win;
  prvd->scene = scene;

  fprintf (stderr, " done. \n");

#else

  prvd->window = NULL;
  prvd->scene = NULL;

#endif

  }

//*============================================================*
//*==========            api_proc_events             ==========*
//*============================================================*
// process event loop.

void
PmGraphicsInterface::api_proc_events (const char *prompt, const char *script) 
  {
  char cmd[100];
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;

#ifdef USE_GRAPHICS

  GrWindow *win = prvd->window;

  fprintf (stderr, "\n>>>>>> PmGraphicsInterface::api_proc_events: \n");
  fprintf (stderr, "   >>> win [%x] \n", win);

  if (script) {
    sprintf (cmd, "read %s", script);
    }

  win->processEvents (NULL, true, prompt, cmd);

#endif
  }

//*============================================================*
//*==========            api_set_extent              ==========*
//*============================================================*
// process event loop.

void 
PmGraphicsInterface::api_set_extent (PmExtent& extent) 
  {
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  float dx, dy, dz, max_dim, cx, cy, cz, xmin, xmax, ymin, ymax,
        zmin, zmax;

  xmin = extent.min[0]; ymin = extent.min[1]; zmin = extent.min[2];
  xmax = extent.max[0]; ymax = extent.max[1]; zmax = extent.max[2];

  dx = xmax - xmin;
  dy = ymax - ymin;
  dz = zmax - zmin;
  max_dim = (dx > dy ? dx : dy);
  max_dim = (dz > max_dim ? dz : max_dim);

  cx = (xmax + xmin) / 2.0;
  cy = (ymax + ymin) / 2.0;
  cz = (zmax + zmin) / 2.0;

  xmin = cx - max_dim;
  xmax = cx + max_dim;
  ymin = cy - max_dim;
  ymax = cy + max_dim;
  zmin = cz - 100.0*max_dim;
  zmax = cz + 100.0*max_dim;

  GrExtent gextent(xmin, xmax, ymin, ymax, zmin, zmax);
  prvd->window->setExtent (gextent);
  prvd->scene->setExtent (gextent);
#endif
  }

//*============================================================*
//*==========            api_get_context             ==========*
//*============================================================*
// create and initialize a graphics context.

void
PmGraphicsInterface::api_get_context (void **ptr)
  {
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  GraphicsContext *context = new GraphicsContext;
  context->window = prvd->window;
  context->scene = prvd->scene;
  context->geom = NULL;
  *ptr = context;
#else
  *ptr = NULL;
#endif
  }

//*============================================================*
//*==========            api_line_display            ==========*
//*============================================================*

void
PmGraphicsBackbone::api_line_display (void *gc) 
  {
#ifdef USE_GRAPHICS
  GraphicsContext *context = (GraphicsContext*)gc;
  GrScene *scene = context->scene;
  GrLine *line = (GrLine *)context->geom;
  line->setColor(color);
  line->setWidth(line_width);
  scene->render();
#endif
  }

//*============================================================*
//*==========            api_line_create             ==========*
//*============================================================*

void
PmGraphicsBackbone::api_line_create (void *gc, const string name, int num_verts, 
                                     PmVector3 *verts) 
  {

#ifdef USE_GRAPHICS

  fprintf (stderr, "\n>>>>>> PmGraphicsBackbone::api_line_display \n");
  fprintf (stderr, "   >>> name      [%s] \n", name.c_str());
  fprintf (stderr, "   >>> num verts [%d] \n", num_verts);

  GraphicsContext *context = (GraphicsContext*)gc;

  GrWindow *win = context->window;
  GrScene *scene = context->scene;
  fprintf (stderr, "   >>> win [%x] \n", win);
  fprintf (stderr, "   >>> scene [%x] \n", scene);
  
  GrLine *line = new GrLine (name, num_verts, verts);
  line->setColor(color);
  line->setWidth(line_width);
  scene->addGeometry (line);
  context->geom = dynamic_cast<GrVertexGeometry*>(line);
#endif
  }

//*============================================================*
//*==========            api_spheres_display         ==========*
//*============================================================*

void
PmGraphicsAtoms::api_spheres_display (void *gc) 
  {
#ifdef USE_GRAPHICS
  fprintf (stderr, "\n>>>>>> PmGraphicsAtoms::api_spheres_display \n");
  GraphicsContext *context = (GraphicsContext*)gc;
  GrScene *scene = context->scene;
  GrSphere *sphere;
  sphere = dynamic_cast<GrSphere*>(context->geom);
  sphere->setColor(color);
  sphere->setDisplayType (gr_ren_map[ren_type]);

  if (radii) {
    sphere->setRadius (radii);
    }

  scene->render();
#endif
  }

//*============================================================*
//*==========            api_spheres_create          ==========*
//*============================================================*

void
PmGraphicsAtoms::api_spheres_create (void *gc, const string name, int num_verts,
                                     PmVector3 *verts, PmVector3 *colors)
  {
#ifdef USE_GRAPHICS
  fprintf (stderr, "\n>>>>>> PmGraphicsAtoms::api_spheres_create \n");
  fprintf (stderr, "   >>> name      [%s] \n", name.c_str());
  fprintf (stderr, "   >>> num verts [%d] \n", num_verts);

  GraphicsContext *context = (GraphicsContext*)gc;
  GrWindow *win = context->window;
  GrScene *scene = context->scene;
  fprintf (stderr, "   >>> win [%x] \n", win);
  fprintf (stderr, "   >>> scene [%x] \n", scene);

  GrSphere *sphere = new GrSphere (name, num_verts, verts);
  sphere->setColor(color);
  sphere->setRadius (0.2);
  fprintf (stderr, "   >>> radii [%x] \n", radii);
  fprintf (stderr, "   >>> colors [%x] \n", colors);

  if (colors) {
    sphere->setColors (num_verts, colors);
    }

  sphere->setShadingType (GR_GEOMETRY_SHADING_COLOR);
  scene->addGeometry (sphere);
  context->geom = dynamic_cast<GrVertexGeometry*>(sphere);
#endif
  }


};

