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
//* graphics:           g r a p h i c s                        *
//*============================================================*

#include "graphics.h"
#include "graphics_prv.h"
#include "bez.h"

namespace ProteinMechanica {

const char* pmGraphicsGeomAtoms          = "atoms";
const char* pmGraphicsGeomBackbone       = "backbone";
const char* pmGraphicsGeomBackbonePlanes = "backbonePlanes";
const char* pmGraphicsGeomBonds          = "bonds";
const char* pmGraphicsGeomHbonds         = "hbonds";

void pm_CmdProcess(FILE *fp, const string cmd);

#ifdef USE_GRAPHICS

static GrGeometryDisplayType gr_molren_map[PM_MOLECULE_RENDER_SIZE+1];
static GrGeometryShadingType gr_shading_map[PM_GEOMETRY_SHADING_SIZE+1];
static GrGeometryDisplayType gr_display_map[PM_GEOMETRY_DISPLAY_SIZE+1];
void PM_EXPORT pm_SystemIdle();

void pickProc(GrScene *scene, GrPickResult& pick_res);
void pickDisplay(PmVector3& point);

static void
graphics_ColorMapHsvToRgb (float h, float s, float v, float *r, float *g, float *b);

#endif

// debugging symbols
#define ndbg_PmGraphicsBackbone 


///////////////////////////////////////////////////////////////
//       G r a p h i c s    I n t e r f a c e               //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsInterface::PmGraphicsInterface(bool enable_graphics) 
  { 
  //fprintf (stderr, ">>>>>> PmGraphicsInterface: ctor \n");
  graphics_enabled = enable_graphics;
  win_width = 800; win_height = 800;
  win_width = 256; win_height = 256;
  win_width = 512; win_height = 512;
  //win_width = 1100; win_height = 1100;
  pick_geom = NULL;
  }

//*============================================================*
//*==========              init                      ==========*
//*============================================================*
// initialize graphics.

void
PmGraphicsInterface::init() 
  { 
  //fprintf (stderr, ">>>>>> PmGraphicsInterface: init \n");

#ifdef USE_GRAPHICS

  prv_data = new GraphicsPrvData;
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;


  // set cmd processing function //

  grSystem.setCommandCallback (pm_CmdProcess);

  // set mapping from pm to gr symbols
  gr_molren_map[PM_MOLECULE_RENDER_UNKNOWN] = GR_GEOMETRY_DISPLAY_SOLID;
  gr_molren_map[PM_MOLECULE_RENDER_POINT] = GR_GEOMETRY_DISPLAY_POINT;
  gr_molren_map[PM_MOLECULE_RENDER_LINE] = GR_GEOMETRY_DISPLAY_LINE;
  gr_molren_map[PM_MOLECULE_RENDER_SOLID] = GR_GEOMETRY_DISPLAY_SOLID;

  gr_shading_map[PM_GEOMETRY_SHADING_NONE] = GR_GEOMETRY_SHADING_NONE;
  gr_shading_map[PM_GEOMETRY_SHADING_FLAT] = GR_GEOMETRY_SHADING_FLAT;
  gr_shading_map[PM_GEOMETRY_SHADING_COLOR] = GR_GEOMETRY_SHADING_COLOR;
  gr_shading_map[PM_GEOMETRY_SHADING_NORMAL] = GR_GEOMETRY_SHADING_NORMAL;

  gr_display_map[PM_GEOMETRY_DISPLAY_UNKNOWN] = GR_GEOMETRY_DISPLAY_POINT;
  gr_display_map[PM_GEOMETRY_DISPLAY_POINT] = GR_GEOMETRY_DISPLAY_POINT;
  gr_display_map[PM_GEOMETRY_DISPLAY_LINE] = GR_GEOMETRY_DISPLAY_LINE;
  gr_display_map[PM_GEOMETRY_DISPLAY_SOLID] = GR_GEOMETRY_DISPLAY_SOLID;

  // set initial extent // 

  GrExtent extent(-10.0, 10.0, -10.0, 10.0, -200.0, 200.0);

  // create a window //

  int wtype = GR_WINDOW_TYPE_DOUBLE | GR_WINDOW_TYPE_RGB | GR_WINDOW_TYPE_DEPTH |
              GR_WINDOW_TYPE_DIRECT;

  string name = "Protein Mechanica";
  int x = 0, y = 0;
  GrWindow *win = new GrWindow (name, x, y, win_width, win_height, wtype);

  win->setIdleFunction(pm_SystemIdle);
  win->setExtent (extent);
  win->open ();

  // create a scene //

  GrScene *scene = new GrScene (name);
  scene->setWindow (win);
  scene->setViewport (win_width, win_height, extent);
  scene->setExtent (extent);

  // create a light //

  GrLight *light = new GrLight ("light 1", GR_LIGHT_DIRECTIONAL);
  scene->addLight (light);

  // attach the scene to the window //

  win->setScene (scene);

  prvd->window = win;
  prvd->scene = scene;

  // initialize picking // 

  initPicking();

#else

  //prvd->window = NULL;
  //prvd->scene = NULL;

#endif
  }

//*============================================================*
//*==========              mapRenderType             ==========*
//*============================================================*

PmGeometryDisplayType 
PmGraphicsInterface::mapRenderType(PmMoleculeRenderType rtype) 
  {
  switch (rtype) {
    case PM_MOLECULE_RENDER_UNKNOWN:
      return (PM_GEOMETRY_DISPLAY_UNKNOWN);
    break;

    case PM_MOLECULE_RENDER_POINT:
      return (PM_GEOMETRY_DISPLAY_POINT);
    break;

    case PM_MOLECULE_RENDER_LINE:
      return (PM_GEOMETRY_DISPLAY_LINE);
    break;

    case PM_MOLECULE_RENDER_SOLID:
      return (PM_GEOMETRY_DISPLAY_SOLID);
    break;
    } 
  }

//*============================================================*
//*==========              update                    ==========*
//*============================================================*
// update graphics.

void
PmGraphicsInterface::update() {
  //fprintf (stderr, "\n>>>>>> PmGraphicsInterface::update \n");
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  prvd->scene->render();
#endif
  }

//*============================================================*
//*==========              setCenter                 ==========*
//*============================================================*
// set graphics rotational center 

void
PmGraphicsInterface::setCenter(PmVector3& center)
  {
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  prvd->scene->setXformCenter(center);
  prvd->scene->render();
#endif
  }

//*============================================================*
//*==========              rotateScene               ==========*
//*============================================================*
// rotate scene.

void 
PmGraphicsInterface::rotateScene(PmVector3& rot)
  {
#ifdef USE_GRAPHICS
  GrVector3 grot;
  grot[0] = rot[0];
  grot[1] = rot[1];
  grot[2] = rot[2];
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  prvd->scene->rotate(false, rot);
  prvd->scene->render();
#endif
  }

//*============================================================*
//*==========              scaleScene                ==========*
//*============================================================*
// scale scene.

void
PmGraphicsInterface::scaleScene(const float scale)
  {
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  prvd->scene->scale(false, scale);
  prvd->scene->render();
#endif
  }

//*============================================================*
//*==========              translateScene            ==========*
//*============================================================*
// translate scene.

void
PmGraphicsInterface::translateScene(PmVector3& vec)
  {
#ifdef USE_GRAPHICS
  GrVector3 gvec;
  gvec[0] = vec[0];
  gvec[1] = vec[1];
  gvec[2] = vec[2];
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  prvd->scene->translate(false, vec);
  prvd->scene->render();
#endif
  }

//*============================================================*
//*==========              getSceneXform             ==========*
//*============================================================*
// get the current scene xform.

void
PmGraphicsInterface::getSceneXform(PmXform& xform)
  {
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  prvd->scene->getXform(xform);
#endif
  }

//*============================================================*
//*==========              setExtent                 ==========*
//*============================================================*
// set graphics extent.               

void 
PmGraphicsInterface::setExtent(PmExtent& extent) 
  {
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsInterface::setExtent \n");
  fprintf (stderr, ">>> extent x = %f %f \n", extent.min[0], extent.max[0]);
  fprintf (stderr, ">>>        y = %f %f \n", extent.min[1], extent.max[1]);
  fprintf (stderr, ">>>        z = %f %f \n", extent.min[2], extent.max[2]);
  */

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

  //fprintf (stderr, ">>> center %f %f %f \n", cx, cy, cz); 
  //fprintf (stderr, ">>> max dim = %f \n", max_dim); 

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
//*==========              procEvents                ==========*
//*============================================================*
// start process events loop.

void
PmGraphicsInterface::procEvents(const char *script)
  {
  char cmd[100];

#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  GrWindow *win = prvd->window;
  char *prompt = "pm> ";

  //fprintf (stderr, "\n>>>>>> PmGraphicsInterface::procEvents: \n");
  //fprintf (stderr, "   >>> win [%x] \n", win);

  if (script) {
    sprintf (cmd, "read %s", script);
    }
  else {
    cmd[0] = '\0';
    }

  win->processEvents (NULL, true, prompt, cmd);

#else

  if (script) {
    sprintf (cmd, "read %s", script);
    pm_CmdProcess (NULL, cmd);
    }

#endif
  }

//*============================================================*
//*==========              getContext                ==========*
//*============================================================*

void
PmGraphicsInterface::getContext (void **gc) {
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  GraphicsContext *context = new GraphicsContext;
  context->window = prvd->window;
  context->scene = prvd->scene;
  context->geom = NULL;
  context->pmgr_geom = NULL;
  *gc = context;
#else
  *gc = NULL;
#endif
  }

//*============================================================*
//*==========              setWindowColor            ==========*
//*============================================================*
// set the window background color.

void
PmGraphicsInterface::setWindowColor(const PmVector3& c) {
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  GrWindow *win = prvd->window;
  GrColor color;
  color[0] = c[0];
  color[1] = c[1];
  color[2] = c[2];
  win->setColor(color);
#endif
  }

//*============================================================*
//*==========              recordWindow              ==========*
//*============================================================*

void
PmGraphicsInterface::recordWindow(const string fname, const string fstr)
  {
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  //fprintf (stderr, ">>>>> PmGraphicsInterface::recordWindow  prvd[%x] \n", prv_data);

  if (!prvd) {
    return;
    }

  GrWindow *win = prvd->window;
  GrImageFormat format;
  //fprintf (stderr, "   >>> win[%x] \n", win);

  if (!win) {
    return;
    }

  if (fstr == "jpeg") {
    format = GR_IMAGE_FORMAT_JPEG;
    }

  GrWindowRecordParams params;
  params.format = format;
  params.prefix = fname;

  win->setRecordParams(params);
  win->setEnableRecord(true);

#endif
  }

//*============================================================*
//*==========              writeWindow               ==========*
//*============================================================*

void
PmGraphicsInterface::writeWindow(const string fname, const string fstr)
  {
  //fprintf (stderr, "\n>>>>>> PmGraphicsInterface::writeWindow \n"); 
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;

  if (!prvd) {
    return;
    }

  GrWindow *win = prvd->window;
  GrImageFormat format; 
  int status;
  //fprintf (stderr, "   >>> win[%x] \n", win); 

  if (!win) {
    return;
    }

  if (fstr == "jpeg") {
    format = GR_IMAGE_FORMAT_JPEG;
    }
  else if (fstr == "ppm") {
    format = GR_IMAGE_FORMAT_PPM;
    }
  else {
    format = GR_IMAGE_FORMAT_JPEG;
    }

  win->writeImage(fname, format, &status);

#endif
  }

///////////////////////////////////////////////////////////////
//                 c o l o r   m a p p i n g                //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========            mapDataToColors             ==========*
//*============================================================*

void
PmGraphicsInterface::mapDataToColors(float vmin, float vmax, int num_data, float *data, 
                                     PmVector3 *colors)
  {
#ifdef USE_GRAPHICS
  vector<PmVector3> cmap;
  float val, f, dv;
  int ci, num_colors;

  // generate color map //

  PmGraphicsInterface::getSpectrumColorMap(16, cmap);
  num_colors = cmap.size();
  dv = vmax - vmin;

  if (dv != 0.0) {
    f = (num_colors-1) / dv;

    for (int i = 0; i < num_data; i++ ) {
      val = data[i];
      ci = (int)(f*(val - vmin));
      if (ci < 0) ci = 0;
      if (ci > num_colors-1) ci = num_colors-1;
      colors[i] = cmap[ci];
      }
    }
  else{
    int cmax, cmin;
    cmin = 0;
    cmax = num_colors - 1;

    for (int i = 0; i < num_data; i++ ) {
      val = data[i];

      if (val < vmax) {
        ci = cmin;
        }
      else {
        ci = cmax;
        }

      colors[i] = cmap[ci];
      }
    }
#endif
  }

//*============================================================*
//*==========      grid_ColorMapSpectrumCreate       ==========*
//*============================================================*

void
PmGraphicsInterface::getSpectrumColorMap(int num_colors, vector<PmVector3>& colors)
  {
#ifdef USE_GRAPHICS
  float inc, r, g, b, hue;
  PmVector3 color;
  int i;

  inc = 240.0 / (float)(num_colors - 1);

  for (i = 0, hue = 240.0; i < num_colors; i++, hue -= inc) {
    graphics_ColorMapHsvToRgb (hue, 1.0, 1.0, &r, &g, &b);
    if (r > 1.0) r = 1.0;
    if (r < 0.0) r = 0.0;
    if (g > 1.0) g = 1.0;
    if (g < 0.0) g = 0.0;
    if (b > 1.0) b = 1.0;
    if (b < 0.0) b = 0.0;
    color[0] = r;
    color[1] = g;
    color[2] = b;
    colors.push_back(color);
    }
#endif
  }

//*============================================================*
//*==========          graphics_ColorMapHsvToRgb     ==========*
//*============================================================*

static void
graphics_ColorMapHsvToRgb (float h, float s, float v, float *r, float *g, float *b)
  {

  float f, p, q, t;
  int i;

  if (s == 0.0) {
    *r = v;
    *g = v;
    *b = v;
    }
  else {
    if (h == 360.0) {
      h = 0.0;
      }

    h = h / 60.0;
    i = (int)h;
    f = h - i;
    p = v * (1.0 - s);
    q = v * (1.0 - (s*f));
    t = v * (1.0 - (s * (1.0 - f)));

    switch (i) {
      case 0:
        *r = v, *g = t, *b = p;
      break;

      case 1:
        *r = q, *g = v, *b = p;
      break;

      case 2:
        *r = p, *g = v, *b = t;
      break;

      case 3:
        *r = p, *g = q, *b = v;
      break;

      case 4:
        *r = t, *g = p, *b = v;
      break;

      case 5:
        *r = v, *g = p, *b = q;
      break;
      }
    }
  }

///////////////////////////////////////////////////////////////
//                 p i c k i n g                            //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========              initPicking               ==========*
//*============================================================*

void
PmGraphicsInterface::initPicking() 
  {
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;
  GrPickAtts pick_atts;

  pick_atts.type = GR_PICK_POINT;
  pick_atts.type = GR_PICK_VERTEX;
  pick_atts.tol = 1.0;

  prvd->scene->setPickAtts (pick_atts);
  prvd->scene->setPickActive (true);
  prvd->scene->requestPick (1, GR_PICK_POINT, pickProc);
  current_picked_point.set(0,0,0);
#endif
  }

#ifdef USE_GRAPHICS

//*============================================================*
//*==========              pickProc                  ==========*
//*============================================================*

void
pickProc(GrScene *scene, GrPickResult& pick_res)
  {
  #define ndbg_pickProc
  #ifdef dbg_pickProc
  fprintf (stderr, "\n>>>>>> pick_proc:  num sel [%d] \n", pick_res.geom_list.size());
  #endif
  string name;
  int entity, index;
  float min_dist;

  #ifdef dbg_pickProc
  for (int i = 0; i < pick_res.geom_list.size(); i++) {
    GrVertexGeometry *geom = pick_res.geom_list[i].geom;
    geom->getName (name);
    geom->getPickedEntity(entity);
    fprintf (stderr, "   >>> geom \"%s\"  entity [%d] \n", name.c_str(), entity);
    fprintf (stderr, "       dist [%g] \n", pick_res.geom_list[i].dist);
    fprintf (stderr, "       intr [%d] \n", pick_res.geom_list[i].intersected);
    fprintf (stderr, "       ipt (%g %g %g) \n", pick_res.geom_list[i].int_pt[0],
                                                 pick_res.geom_list[i].int_pt[1],
                                                 pick_res.geom_list[i].int_pt[2]);
    }
  #endif

  // first find the closest geometry that was intersected, if any //

  min_dist = 1e6;
  index = -1;
  
  for (unsigned int i = 0; i < pick_res.geom_list.size(); i++) {
    if (pick_res.geom_list[i].intersected && (pick_res.geom_list[i].dist < min_dist)) {
      index = i;
      min_dist = pick_res.geom_list[i].dist;
      }
    }

  // next just find closest geometry //

  if (index == -1) {
    min_dist = 1e6;
  
    for (unsigned int i = 0; i < pick_res.geom_list.size(); i++) {
      if (pick_res.geom_list[i].dist < min_dist) {
        index = i;
        min_dist = pick_res.geom_list[i].dist;
        }
      }
    }

  if (index == -1) {
    return;
    }

  GrPickGeom *pgeom = &pick_res.geom_list[index];
  GrVertexGeometry *geom = pgeom->geom;

  int n;
  string obj_type;

  geom->getName (name);
  geom->getPickedEntity(entity);
  fprintf (stderr, "\n");
  fprintf (stderr, ">>> picked geometry \"%s\" \n", name.c_str());
  fprintf (stderr, "    entity [%d] \n", entity);
  fprintf (stderr, "    point (%g %g %g) \n", pgeom->int_pt[0],
                                              pgeom->int_pt[1],
                                              pgeom->int_pt[2]);
  fprintf (stderr, "    ---------------------------------------\n");

  // determine if object was created by something that can be queried //

  n = name.find("[");

  if (n == (int)string::npos) {
    return;
    }

  // send query to appropriate object //

  PmQuery query;
  query.name = name;
  query.point = pgeom->int_pt;
  query.entity = entity;

  pmSystem.procQuery(query);

  // show the picked point //

  PmGraphicsInterface *grint = pmSystem.getGraphicsInterace();

  grint->displayPick(pgeom->int_pt);
  }

#endif

//*============================================================*
//*==========              displayPick               ==========*
//*============================================================*

void
PmGraphicsInterface::displayPick(PmVector3& point)
  {
#ifdef USE_GRAPHICS
  GraphicsPrvData *prvd = (GraphicsPrvData*)prv_data;

  if (!pick_geom) { 
    PmGraphicsAttributes atts;
    atts.setMarker(true);
    atts.setScale(0.05);
    string geom_name = "pm:pick_point";
    PmVector3 *verts = new PmVector3[1];
    verts[0] = point;
    pick_geom = new PmGraphicsPoint(geom_name, 1, verts);
    pick_geom->setAttributes(atts);
    pick_geom->display();
    }
  else {
    PmVector3 verts[1]; 
    verts[0] = point;
    pick_geom->update(1, &point);
    }

  prvd->scene->render();
  current_picked_point = point;
#endif
  }

//*============================================================*
//*==========              getCurrentPick            ==========*
//*============================================================*

void
PmGraphicsInterface::getCurrentPick(PmVector3& point)
  {
  point = current_picked_point;
  }

///////////////////////////////////////////////////////////////
//          G r a p h i c s     G e o m e t r y             //
/////////////////////////////////////////////////////////////

PmGraphicsGeometry::PmGraphicsGeometry() 
  {
  num_vertices = 0;
  vertices = NULL;
  colors = NULL;
  copy_vertices = false;
  }

PmGraphicsGeometry::PmGraphicsGeometry (string name) 
  {
  this->name = name;
  }

void 
PmGraphicsGeometry::setAttributes(PmGraphicsAttributes& atts) 
  { 
  this->attributes = atts; 
  }

//*============================================================*
//*==========              setColors                 ==========*
//*============================================================*

void 
PmGraphicsGeometry::setColors(int num_colors, PmVector3 *colors) 
  { 
  this->num_colors = num_colors; 
  this->colors = colors; 
  }

void 
PmGraphicsGeometry::setColors(vector<PmVector3>& cvals)
  {
  if (cvals.size() == 0) {
    return;
    }

  if (cvals.size() != num_vertices) {
    pm_ErrorWarnReport (PM, "colors size (%d) not equal to vertex size (%d)", "*",
                        cvals.size(), num_vertices); 
    return;
    }

  this->num_colors = cvals.size(); 
  this->colors = new PmVector3[cvals.size()]; 

  for (unsigned int i = 0; i < cvals.size(); i++) {
    this->colors[i] = cvals[i]; 
    }
  }

//*============================================================*
//*==========              convDisplayType           ==========*
//*============================================================*
// convert a string into a display type.

void 
PmGraphicsGeometry::convDisplayType (const string s, PmGeometryDisplayType& dtype)
  {
  if (s == "point") {
    dtype = PM_GEOMETRY_DISPLAY_POINT;
    }
  else if (s == "line") {
    dtype = PM_GEOMETRY_DISPLAY_LINE;
    }
  else if (s == "solid") {
    dtype = PM_GEOMETRY_DISPLAY_SOLID;
    }
  else {
    dtype = PM_GEOMETRY_DISPLAY_UNKNOWN;
    }
  }

//*============================================================*
//*==========              convShadingType           ==========*
//*============================================================*
// convert a string into a shading type.

void
PmGraphicsGeometry::convShadingType (const string s, PmGeometryShadingType& stype)
  {
  if (s == "none") {
    stype = PM_GEOMETRY_SHADING_NONE;
    }
  else if (s == "flat") {
    stype = PM_GEOMETRY_SHADING_FLAT;
    }
  else if (s == "color") {
    stype = PM_GEOMETRY_SHADING_COLOR;
    }
  else {
    stype = PM_GEOMETRY_SHADING_UNKNOWN;
    }
  }

//*============================================================*
//*==========               update                   ==========*
//*============================================================*
// update geometry with new vertices.

void
PmGraphicsGeometry::update(int num_verts, PmVector3 *verts)
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GrGeometry *geom = (GrGeometry*)gc->geom;
  geom->setVertices(num_verts, verts);
#endif
  }

//*============================================================*
//*==========               updateColors             ==========*
//*============================================================*
// update geometry with new colors.

void 
PmGraphicsGeometry::updateColors(int num_colors, PmVector3 *colors)
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GrGeometry *geom = (GrGeometry*)gc->geom;
  geom->setColors(num_colors, colors);
#endif
  }

/////////////////////////////////////////////////////////////////
//            G r a p h i c s    A t o m s                    //
///////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsAtoms::PmGraphicsAtoms(string name, int num_verts, PmVector3 *verts,
                                 PmVector3 *colors)
  {
  this->name = name;
  this->num_vertices = num_verts;
  this->vertices = verts;
  this->colors = colors;
  radii = NULL;
  disp_type = PM_MOLECULE_DISPLAY_SPHERE;
  ren_type = PM_MOLECULE_RENDER_SOLID;
  attributes.color.set(1,1,1);
  pmSystem.getGraphicsContext (&this->context);
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display atoms.

void
PmGraphicsAtoms::display()
  {
#ifdef USE_GRAPHICS
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsAtoms::display \n");
  fprintf (stderr, "   >>> name [%s] \n", this->name.c_str());
  fprintf (stderr, "   >>> has radii \n");
  fprintf (stderr, "   >>> atts scale = %f \n", attributes.scale);
  */
  GraphicsContext *gc = (GraphicsContext*)context;
  GrSphere *sphere;

  if (!gc->geom) {
    GrScene *scene = gc->scene;
    sphere = new GrSphere (name, num_vertices, vertices);
    //sphere->setOptimized(false);

    if (radii) {
      sphere->setRadius (radii);
      }

    if (colors) {
      sphere->setColors (num_vertices, colors);
      }

    gc->geom = dynamic_cast<GrVertexGeometry*>(sphere);
    scene->addGeometry(sphere);
    }
  else {
    sphere = dynamic_cast<GrSphere*>(gc->geom);
    }

  sphere->setRadius(attributes.scale);
  sphere->setColor(attributes.color);
  sphere->setDisplayType (gr_molren_map[ren_type]);
  sphere->setResolution(16);

  if (xform.set) {
    sphere->setXform(xform);
    }

  sphere->setActive(attributes.visible);
#endif
  }

/////////////////////////////////////////////////////////////////
//            G r a p h i c s    A x e s                      //
///////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsAxes::PmGraphicsAxes(const string name, PmVector3 center, PmVector3 axis1,
                               PmVector3 axis2, PmVector3 axis3)
  {
  this->name = name;
  num_vertices = 6;
  vertices = new PmVector3[6];
  vertices[0] = center - 0.5*axis1;
  vertices[1] = center + 0.5*axis1;
  vertices[2] = center - 0.5*axis2;
  vertices[3] = center + 0.5*axis2;
  vertices[4] = center - 0.5*axis3;
  vertices[5] = center + 0.5*axis3;
  pmSystem.getGraphicsContext (&this->context);
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display axes geometry.

void
PmGraphicsAxes::display()
  {
#ifdef USE_GRAPHICS
  //fprintf (stderr, "\n>>>>>> PmGraphicsAxes::display\n");
  GraphicsContext *gc = (GraphicsContext*)context;
  GrScene *scene = gc->scene;
  GrLine *line;

  if (!gc->geom) { 
    line = new GrLine (name, num_vertices, vertices);
    line->setDisjoint(true);
    scene->addGeometry(line);
    gc->geom = dynamic_cast<GrVertexGeometry*>(line);
    }
  else {
    line = (GrLine*)gc->geom;
    }

  line->setColor(attributes.color);
  line->setWidth(attributes.line_width);
#endif
  }

/////////////////////////////////////////////////////////////////
//            G r a p h i c s    B a c k b o n e              //
///////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsBackbone::PmGraphicsBackbone() {
  num_vertices = 0;
  vertices = NULL;
  connectivity = NULL;
  }

PmGraphicsBackbone::PmGraphicsBackbone(const string name, int num_verts, 
                                       PmVector3 *verts) 
  {
  this->name = name;
  this->num_vertices = num_verts; 
  this->vertices = verts; 
  this->number = 0;
  connectivity = NULL;
  pmSystem.getGraphicsContext (&this->context);
  }

PmGraphicsBackbone::PmGraphicsBackbone(const string name, int num, PmConn *conn,
                                       int num_verts, PmVector3 *verts) 
  {
  this->name = name;
  this->num_vertices = num_verts; 
  this->vertices = verts; 
  this->number = num;
  this->connectivity = conn;
  pmSystem.getGraphicsContext (&this->context);
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create backbone geometry.

void
PmGraphicsBackbone::create()
  {
#ifdef USE_GRAPHICS
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsBackboneTube::create \n");
  fprintf (stderr, "   >>> name [%s] \n", this->name.c_str());
  */
  GraphicsContext *gc = (GraphicsContext*)context;
  if (gc->geom) return; 
  GrScene *scene = gc->scene;
  GrLine *line;

  if (connectivity) {
    line = new GrLine (name, number, *connectivity, num_vertices, vertices);
    }
  else {
    line = new GrLine (name, num_vertices, vertices);
    }

  line->setColor(attributes.color);
  line->setWidth(attributes.line_width);
  scene->addGeometry (line);
  gc->geom = dynamic_cast<GrVertexGeometry*>(line);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display backbone.

void
PmGraphicsBackbone::display()
  {
#ifdef USE_GRAPHICS
  #ifdef dbg_PmGraphicsBackbone
  fprintf (stderr, "\n>>>>>> PmGraphicsBackbone::display \n");
  fprintf (stderr, "   >>> name [%s] \n", this->name.c_str());
  fprintf (stderr, "   >>> xform[%d] \n", this->xform.set);
  #endif
  GraphicsContext *gc = (GraphicsContext*)context;
  //GrScene *scene = gc->scene;

  create();

  GrLine *line = (GrLine*)gc->geom;
  line->setColor(attributes.color);
  line->setWidth(attributes.line_width);
  line->setActive(attributes.visible);
  //fprintf (stderr, "\n>>>>>> PmGraphicsBackbone::disp visible=%d\n", attributes.visible);

  if (xform.set) {
    #ifdef dbg_PmGraphicsBackbone
    fprintf (stderr, "   >>> xform: translate (%g %g %g) \n",
             this->xform.translation[0], this->xform.translation[1],
             this->xform.translation[2]);
    #endif
    line->setXform(xform);
    }

#endif
  }

/////////////////////////////////////////////////////////////////
//        G r a p h i c s    T u b e   B a c k b o n e        //
///////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsBackboneTube::PmGraphicsBackboneTube() {
  num_vertices = 0;
  vertices = NULL;
  connectivity = NULL;
  resolution = 4;
  width = 0.01;
  }

PmGraphicsBackboneTube::PmGraphicsBackboneTube(const string name, int num_verts, 
                                               PmVector3 *verts, int res, float width) 
  {
  this->name = name;
  this->num_vertices = num_verts; 
  this->vertices = verts; 
  this->number = 0;
  this->resolution = res;
  this->width = width;
  connectivity = NULL;
  pmSystem.getGraphicsContext (&this->context);
  }

PmGraphicsBackboneTube::PmGraphicsBackboneTube(const string name, int num, PmConn *conn,
                                               int num_verts, PmVector3 *verts, int res, 
                                               float width) 
  {
  this->name = name;
  this->num_vertices = num_verts; 
  this->vertices = verts; 
  this->number = num;
  this->connectivity = conn;
  this->resolution = res;
  this->width = width;
  pmSystem.getGraphicsContext (&this->context);
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create backbone geometry.

void
PmGraphicsBackboneTube::create()
  {
#ifdef USE_GRAPHICS
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsBackboneTube::create \n");
  fprintf (stderr, ">>> name=%s \n", this->name.c_str());
  fprintf (stderr, ">>> num_vertices=%d \n", num_vertices);
  */
  GraphicsContext *gc = (GraphicsContext*)context;
  if (gc->geom) return; 
  GrScene *scene = gc->scene;
  GrLine *line;

  /*
  if (connectivity) {
    line = new GrLine (name, number, *connectivity, num_vertices, vertices);
    }
  else {
    line = new GrLine (name, num_vertices, vertices);
    }

  line->setColor(attributes.color);
  line->setWidth(attributes.line_width);
  scene->addGeometry (line);
  gc->geom = dynamic_cast<GrVertexGeometry*>(line);
  */

  PmVector3 u, v, w;
  PmVector3 p1, p2, p3, p;
  float ang, dang, r;
  int n;
  vector<PmVector3> pts;
  float error, dp;
  stringstream dss;
  string pname;
  int num_curves;
  PmVector3 *curves;
  int num_pts;
  int num_sects;
  int num_div;
  PmVector3 *cpts, *ctang, *cnorm;
  PmVector3 *tpts, *npts, tang, norm, norm1;

  error = 0.001;
  error = 0.01;
  r = 0.02;

  pm_BezierFit (num_vertices, vertices, error, &num_curves, &curves);
  //fprintf (stderr, "\n>>> num curves=%d \n", num_curves);

  num_div = 50;
  num_div = 10;
  cpts = (PmVector3*)malloc(sizeof(PmVector3) * num_div * num_curves);
  ctang = (PmVector3*)malloc(sizeof(PmVector3) * num_div * num_curves);
  cnorm = (PmVector3*)malloc(sizeof(PmVector3) * num_div * num_curves);
  tpts = (PmVector3*)malloc(sizeof(PmVector3) * 2*num_div * num_curves);
  npts = (PmVector3*)malloc(sizeof(PmVector3) * 2*num_div * num_curves);

  for (int i = 0; i < num_curves; i++) {
    pm_BezierEval (4, curves+4*i, num_div, cpts+i*num_div, ctang+i*num_div,
                   cnorm+i*num_div);
    }

  num_sects = num_curves * num_div;
  float s = 0.01;

  // make sure adjacent normals are in the same direction //

  for (int i = 1; i < num_sects; i++) {
    norm = cnorm[i];
    norm.normalize();
    norm1 = cnorm[i-1];
    norm1.normalize();
    dp = norm*norm1;

    if (dp < 0.0) {
      cnorm[i] = -norm;
      }

    if (fabs(dp) < 0.7) {
      cnorm[i] = norm1;
      }

    if (i % (num_div) == 0) {
      norm = 0.5*(cnorm[i] + cnorm[i-1]);
      norm.normalize();
      cnorm[i-1] = norm;
      }

    //fprintf (stderr, ">>> %d: dp=%f \n", i, dp);
    }

  for (int i = 0; i < num_sects; i++) {
    tang = ctang[i];
    tang.normalize();
    tpts[2*i] = cpts[i];
    tpts[2*i+1] = cpts[i] + r*tang;

    norm = cnorm[i];
    norm.normalize();
    npts[2*i] = cpts[i];
    npts[2*i+1] = cpts[i] + r*norm;

    //fprintf (stderr, ">>> %d norm=%f %f %f\n", i, norm[0], norm[1], norm[2]);
    //fprintf (stderr, "    norm*tang=%f \n", norm*tang);
    }

  /*
  dss << name << '_' << "Bez";
  pname = dss.str();
  dss.str(std::string());
  PmVector3 color;

  GrPoint *geom_pts = new GrPoint(pname, num_sects, cpts);
  scene->addGeometry (geom_pts);

  GrLine *geom_line = new GrLine(pname, num_sects, cpts);
  scene->addGeometry (geom_line);

  GrLine *geom_tline = new GrLine(pname, 2*num_sects, tpts);
  geom_tline->setDisjoint(true);
  color.set(0,1,1);
  geom_tline->setColor(color);
  scene->addGeometry (geom_tline);

  GrLine *geom_nline = new GrLine(pname, 2*num_sects, npts);
  color.set(1,0,0);
  geom_nline->setDisjoint(true);
  geom_nline->setColor(color);
  scene->addGeometry (geom_nline);
  */

  n = 32;
  dang = 2.0*M_PI / n;

  for (int i = 0; i < num_sects; i++) {
    if ((i != 0) && (i % (num_div) == 0)) {
      //fprintf (stderr, ">>> skip %d \n", i);
      continue;
      }

    p1 = cpts[i];
    tang = ctang[i];
    tang.normalize();
    norm = cnorm[i];
    norm.normalize();

    u = norm;
    v = norm.cross(tang);
    ang = 0.0;
 
    //fprintf (stderr, "%d u=%f %f %f \n", i, u[0], u[1], u[2]); 
    //fprintf (stderr, "  v=%f %f %f \n", v[0], v[1], v[2]); 
    //fprintf (stderr, "  u*v=%f \n", u*v); 
    //fprintf (stderr, "  p1=%f %f %f \n", p1[0], p1[1], p1[2]); 

    for (int j = 0; j < n; j++) {
      p = p1 + r*cos(ang)*u + r*sin(ang)*v;
      pts.push_back(p);
      //fprintf (stderr, "  p=%f %f %f \n", p[0], p[1], p[2]); 
      ang += dang;
      }
    }

  int num_pverts = pts.size();
  PmVector3 *pverts = new PmVector3[num_pverts];

  for (int j = 0; j < num_pverts; j++) {
    pverts[j] = pts[j];
    //fprintf (stderr, "%d %f %f %f \n", j, pverts[j][0],  pverts[j][1],  pverts[j][2]); 
    }

  PmConn *conn = new PmConn(5*n*num_curves*(num_div-1));
  int *pconn = conn->vals;
  int m=0;
  int off = 0;
  int k = 0;
  int num_poly=0;

  for (int i = 0; i < num_sects-1; i++) {
    if ((i != 0) && (i % (num_div) == 0)) {
      continue;
      }

    for (int j = 0; j < n; j++) {
      k = j + 1;

      if (k == n) {
        k = 0;
        }

      pconn[m++] = 4;
      pconn[m++] = off + j;
      pconn[m++] = off + j + n;
      pconn[m++] = off + k + n;
      pconn[m++] = off + k;
      num_poly++;
      }

    off += n; 
    }

  dss << name << '_' << "tube";
  pname = dss.str();
  dss.str(std::string());

  //GrPoint *pgeom_pts = new GrPoint(pname, num_pverts, pverts);
  //scene->addGeometry (pgeom_pts);

  /*
  fprintf (stderr, ">>> num poly=%d \n", num_poly); 
  fprintf (stderr, ">>> num pverts=%d \n", num_pverts); 
  fprintf (stderr, ">>> num_sects=%d \n", num_sects); 
  fprintf (stderr, ">>> m=%d \n", m); 
  fprintf (stderr, ">>> 5*n*num_curves*(num_div-1)-1=%d \n", 
           5*n*num_curves*(num_div-1)); 
  */

  GrPolygon *poly;
  poly = new GrPolygon(pname, num_poly, *conn, num_pverts, pverts, 0);
  scene->addGeometry (poly);
  //gc->geom = dynamic_cast<GrVertexGeometry*>(poly);

#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display backbone.                             

void
PmGraphicsBackboneTube::display() 
  {
#ifdef USE_GRAPHICS
  #ifdef dbg_PmGraphicsBackbone 
  fprintf (stderr, "\n>>>>>> PmGraphicsBackbone::display \n");
  fprintf (stderr, "   >>> name [%s] \n", this->name.c_str());
  fprintf (stderr, "   >>> xform[%d] \n", this->xform.set);
  #endif 
  GraphicsContext *gc = (GraphicsContext*)context;
  //GrScene *scene = gc->scene;

  create();

  GrLine *line = (GrLine*)gc->geom;
  line->setColor(attributes.color);
  line->setWidth(attributes.line_width);
  line->setActive(attributes.visible);
  //fprintf (stderr, "\n>>>>>> PmGraphicsBackbone::disp visible=%d\n", attributes.visible);

  if (xform.set) {
    #ifdef dbg_PmGraphicsBackbone 
    fprintf (stderr, "   >>> xform: translate (%g %g %g) \n", 
             this->xform.translation[0], this->xform.translation[1],
             this->xform.translation[2]);
    #endif
    line->setXform(xform);
    }

#endif
  }

///////////////////////////////////////////////////////////////
//          G r a p h i c s    E l l ip s o i d             //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsEllipsoid::PmGraphicsEllipsoid(const string name, int num_verts, 
                                         PmVector3 *verts, PmEllipsoidData *data)
  {
  this->name = name;
  this->num_vertices = num_verts;
  this->vertices = verts;
  this->axes_data = data;
  attributes.display_type = PM_GEOMETRY_DISPLAY_SOLID;
  attributes.display_type = PM_GEOMETRY_DISPLAY_LINE;
  pmSystem.getGraphicsContext (&this->context);
  attributes.color.set(1,0,1);
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display an ellipsoid.

void
PmGraphicsEllipsoid::display()
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GrScene *scene = gc->scene;
  GrEllipsoid *ellipsoid;

  if (!gc->geom) {
    GrEllipsoidData *data = new GrEllipsoidData[num_vertices];

    for (int i = 0; i < num_vertices; i++) {
      data[i].rad1 = axes_data[i].radius1;
      data[i].axis1 = axes_data[i].axis1;
      data[i].rad2 = axes_data[i].radius2;
      data[i].axis2 = axes_data[i].axis2;
      data[i].rad3 = axes_data[i].radius3;
      data[i].axis3 = axes_data[i].axis3;
      }

    ellipsoid = new GrEllipsoid(name, num_vertices, vertices, data);
    ellipsoid->setColor(attributes.color);
    ellipsoid->setDisplayType(gr_display_map[attributes.display_type]);
    ellipsoid->setShadingType(gr_shading_map[attributes.shading_type]);
    ellipsoid->setLineWidth(attributes.line_width);

    scene->addGeometry (ellipsoid);
    gc->geom = dynamic_cast<GrVertexGeometry*>(ellipsoid);
    }
  else {
    ellipsoid= (GrEllipsoid*)gc->geom;
    }

  if (xform.set) {
    ellipsoid->setXform(xform);
    }

#endif
  }

///////////////////////////////////////////////////////////////
//          G r a p h i c s    J o i n t                    //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsJoint::PmGraphicsJoint(const string name, int n, PmVector3 *verts) 
  {
  this->name = name;
  num_vertices = n;
  vertices = verts;
  attributes.scale = 0.5;
  attributes.color.set(1,0,1);
  pmSystem.getGraphicsContext (&this->context);
  pmSystem.getGraphicsContext (&this->aux_context);
  }

void
PmGraphicsJoint::update(const PmVector3& pos)
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  PmVector3 verts[1];
  verts[0] = pos;
  gc->geom->setVertices(1, verts);
#endif
  }

//*============================================================*
//*=================  b a l l  j o i n t  =====================*
//*============================================================*

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create ball joint.

void
PmGraphicsBallJoint::create() 
  {
#ifdef USE_GRAPHICS
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsBallJoint::create\n");
  fprintf (stderr, "   >>> name   [%s] \n", this->name.c_str());
  */
  GraphicsContext *gc = (GraphicsContext*)context;
  if (gc->geom) return;

  float length, radius;
  GrScene *scene = gc->scene;
  GrSphere *sphere = new GrSphere (name, 1, vertices);

  length = 4.0*attributes.scale;
  radius = attributes.scale / 3.0;

  sphere->setRadius(attributes.scale);
  sphere->setActive(attributes.visible); 
  scene->addGeometry (sphere);
  gc->geom = dynamic_cast<GrVertexGeometry*>(sphere);

  GraphicsContext *cgc = (GraphicsContext*)aux_context;

  for (int i = 0; i < 3; i++) {
    vertices[i] = vertices[i] + attributes.scale*axes[i];
    }

  GrCylinder *cyl = new GrCylinder (name, 3, vertices);
  cyl->setAxis (axes);
  cyl->setLength(length);
  cyl->setRadius(radius);
  cyl->setShadingType(gr_shading_map[attributes.shading_type]);
  cyl->setColor(attributes.color);
  cyl->setActive(attributes.visible_aux); 
  scene->addGeometry (cyl);
  cgc->geom = dynamic_cast<GrVertexGeometry*>(cyl);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display ball joint.

void
PmGraphicsBallJoint::display() 
  {
#ifdef USE_GRAPHICS
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsBallJoint::display\n");
  fprintf (stderr, "   >>> name   [%s] \n", this->name.c_str());
  fprintf (stderr, "   >>> scale [%f] \n", attributes.scale);
  fprintf (stderr, "   >>> vis    [%d] \n", attributes.visible);
  fprintf (stderr, "   >>> aux    [%d] \n", attributes.visible_aux);
  */

  GraphicsContext *gc = (GraphicsContext*)context;
  GraphicsContext *cgc = (GraphicsContext*)aux_context;

  create();

  GrSphere *sphere = dynamic_cast<GrSphere*>(gc->geom);
  sphere->setColor(attributes.color);
  sphere->setDisplayType(gr_display_map[attributes.display_type]);
  sphere->setShadingType(gr_shading_map[attributes.shading_type]);
  sphere->setRadius (attributes.scale);
  sphere->setActive(attributes.visible); 

  GrCylinder *cyl = dynamic_cast<GrCylinder*>(cgc->geom);
  cyl->setRadius(attributes.scale/3.0);
  cyl->setLength(4.0*attributes.scale);
  cyl->setActive(attributes.visible_aux); 
#endif
  }

//*============================================================*
//*==========              xform                     ==========*
//*============================================================*
// xform a ball joint.

void
PmGraphicsBallJoint::xform(PmXform& xform)
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GrSphere *sphere = dynamic_cast<GrSphere*>(gc->geom);
  GraphicsContext *cgc = (GraphicsContext*)aux_context;
  GrCylinder *cyl = dynamic_cast<GrCylinder*>(cgc->geom);

  if (xform.set) {
    sphere->setXform(xform);
    cyl->setXform(xform);
    }
#endif
  }

//*============================================================*
//*==========              update                    ==========*
//*============================================================*
// update the position of a ball joint.

/*
void
PmGraphicsBallJoint::update(const PmVector3& pos)
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GrSphere *sphere = dynamic_cast<GrSphere*>(gc->geom);
  PmVector3 verts[1];
  verts[0] = pos;
  sphere->setVertices(1, verts);
#endif
  }
*/

//*============================================================*
//*=================  free joint  =============================*
//*============================================================*

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create free joint.

void
PmGraphicsFreeJoint::create() 
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  if (gc->geom) return;
  GrScene *scene = gc->scene;
  GrSphere *sphere = new GrSphere (name, 1, vertices);
  scene->addGeometry (sphere);
  gc->geom = dynamic_cast<GrVertexGeometry*>(sphere);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display free joint.

void
PmGraphicsFreeJoint::display() 
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  create();
  GrSphere *sphere = dynamic_cast<GrSphere*>(gc->geom);
  sphere->setColor(attributes.color);
  sphere->setDisplayType(gr_display_map[attributes.display_type]);
  sphere->setShadingType(gr_shading_map[attributes.shading_type]);
  sphere->setRadius (attributes.scale);
  sphere->setActive(attributes.visible);
#endif
  }

//*============================================================*
//*=============  w e l d   j o i n t  ========================*
//*============================================================*

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create weld joint.

void
PmGraphicsWeldJoint::create() 
  {
#ifdef USE_GRAPHICS
  //fprintf (stderr, "\n>>>>>> PmGraphicsWeldJoint::create \n");
  GraphicsContext *gc = (GraphicsContext*)context;
  if (gc->geom) return;

  PmVector3 pverts[8], vx, vy, vz;
  GrIndex conn(30);
  int n;

  GrScene *scene = gc->scene;

  // create box geometry //

  vx.set(1,0,0);
  vy.set(0,1,0);
  vz.set(0,0,-1);

  pverts[0][0] = vertices[0][0] - 0.5;
  pverts[0][1] = vertices[0][1] - 0.5;
  pverts[0][2] = vertices[0][2] + 0.5;
  pverts[1] = pverts[0] + vx;
  pverts[2] = pverts[0] + vx + vy;
  pverts[3] = pverts[0] + vy;

  pverts[4] = pverts[0] + vz;
  pverts[5] = pverts[0] + vx + vz;
  pverts[6] = pverts[0] + vx + vy + vz;
  pverts[7] = pverts[0] + vy + vz;

  // create box topology //

  n = 0;
  conn[n++] = 4; conn[n++] = 0; conn[n++] = 1; conn[n++] = 2;  conn[n++] = 3;
  conn[n++] = 4; conn[n++] = 4; conn[n++] = 7; conn[n++] = 6;  conn[n++] = 5;
  conn[n++] = 4; conn[n++] = 1; conn[n++] = 5; conn[n++] = 6;  conn[n++] = 2;
  conn[n++] = 4; conn[n++] = 0; conn[n++] = 3; conn[n++] = 7;  conn[n++] = 4;
  conn[n++] = 4; conn[n++] = 0; conn[n++] = 4; conn[n++] = 5;  conn[n++] = 1;
  conn[n++] = 4; conn[n++] = 3; conn[n++] = 2; conn[n++] = 6;  conn[n++] = 7;

  GrPolygon *ipoly = new GrPolygon (name, 6, conn, 8, pverts);
  scene->addGeometry(ipoly);
  gc->geom = dynamic_cast<GrVertexGeometry*>(ipoly);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display weld joint.

void
PmGraphicsWeldJoint::display() 
  {
#ifdef USE_GRAPHICS
  //fprintf (stderr, "\n>>>>>> PmGraphicsWeldJoint::display \n");
  PmVector3 pverts[8], vx, vy, vz;
  float s;
  GraphicsContext *gc = (GraphicsContext*)context;

  create();

  GrScene *scene = gc->scene;

  // create box geometry //

  s = attributes.scale;
  vx.set(s,0,0);
  vy.set(0,s,0);
  vz.set(0,0,-s);

  pverts[0][0] = vertices[0][0] - 0.5*s;
  pverts[0][1] = vertices[0][1] - 0.5*s;
  pverts[0][2] = vertices[0][2] + 0.5*s;
  pverts[1] = pverts[0] + vx;
  pverts[2] = pverts[0] + vx + vy;
  pverts[3] = pverts[0] + vy;

  pverts[4] = pverts[0] + vz;
  pverts[5] = pverts[0] + vx + vz;
  pverts[6] = pverts[0] + vx + vy + vz;
  pverts[7] = pverts[0] + vy + vz;

  GrPolygon *ipoly = dynamic_cast<GrPolygon*>(gc->geom);
  ipoly->setVertices(8, pverts);

  ipoly->setColor(attributes.color);
  ipoly->setDisplayType(gr_display_map[attributes.display_type]);
  ipoly->setShadingType(gr_shading_map[attributes.shading_type]);
  ipoly->setActive(attributes.visible);
#endif
  }

//*============================================================*
//*==========              xform                     ==========*
//*============================================================*
// xform a weld joint.

void
PmGraphicsWeldJoint::xform(PmXform& xform)
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;

  if (!xform.set || !gc->geom) {
    return;
    }

  GrPolygon *ipoly = dynamic_cast<GrPolygon*>(gc->geom);
  ipoly->setXform(xform);
#endif
  }

//*============================================================*
//*=============  h i n g e   j o i n t  ======================*
//*============================================================*

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create hinge joint.

void
PmGraphicsHingeJoint::create() 
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  if (gc->geom) return;
  //fprintf (stderr, "\n>>>>>> PmGraphicsHingeJoint::create \n");
  //fprintf (stderr, ">>> num vertices=%d \n", num_vertices); 

  float length, radius;
  GrScene *scene = gc->scene;
  GrCylinder *cyl;

  if (num_vertices == 1) {
    cyl = new GrCylinder (name, 1, vertices);
    length = attributes.scale;
    radius = attributes.scale / 10.0;
    }
  else {
    PmVector3 v;
    v = vertices[1] - vertices[0];
    length = v.length();
    vertices[1] = (vertices[1] + vertices[0]) / 2.0;
    cyl = new GrCylinder (name, 1, &vertices[1]);
    radius = 0.005;
    }

  //fprintf (stderr, ">>> length=%f \n", length); 
  //fprintf (stderr, ">>> radius=%f \n", radius); 

  cyl->setAxis (axis);
  cyl->setLength(length);
  cyl->setRadius(radius);
  scene->addGeometry (cyl);
  gc->geom = dynamic_cast<GrVertexGeometry*>(cyl);

  GraphicsContext *cgc = (GraphicsContext*)aux_context;
  GrSphere *sphere = new GrSphere (name, 1, vertices);
  sphere->setRadius(attributes.scale / 8.0);
  sphere->setActive(attributes.visible);
  scene->addGeometry (sphere);
  cgc->geom = dynamic_cast<GrVertexGeometry*>(sphere);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display hinge joint.

void
PmGraphicsHingeJoint::display() 
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  //fprintf (stderr, "\n>>>>>> PmGraphicsHingeJoint::display\n");
  //fprintf (stderr, ">>> visible=%d \n", attributes.visible); 

  create();

  GrCylinder *cyl = dynamic_cast<GrCylinder*>(gc->geom);
  cyl->setColor(attributes.color);
  cyl->setDisplayType(gr_display_map[attributes.display_type]);
  cyl->setShadingType(gr_shading_map[attributes.shading_type]);
  //cyl->setRadius(attributes.scale/10.0);
  //cyl->setLength(attributes.scale);
  cyl->setActive(attributes.visible);

  GraphicsContext *cgc = (GraphicsContext*)aux_context;
  GrSphere *sphere = dynamic_cast<GrSphere*>(cgc->geom);
  sphere->setColor(attributes.color);
  sphere->setDisplayType(gr_display_map[attributes.display_type]);
  sphere->setShadingType(gr_shading_map[attributes.shading_type]);
  sphere->setRadius (attributes.scale / 8.0);
  sphere->setActive(attributes.visible);
#endif
  }

//*============================================================*
//*==========              xform                     ==========*
//*============================================================*
// xform a hinge joint.

void
PmGraphicsHingeJoint::xform(PmXform& xform)
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GraphicsContext *cgc = (GraphicsContext*)aux_context;
  GrCylinder *cyl = dynamic_cast<GrCylinder*>(gc->geom);
  GrSphere *sphere = dynamic_cast<GrSphere*>(cgc->geom);

  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsHingeJoint::xform\n");
  fprintf (stderr, ">>> center      = %f %f %f \n", xform.center[0],
                                                    xform.center[1],
                                                    xform.center[2]);

  fprintf (stderr, ">>> translation = %f %f %f \n", xform.translation[0],
                                                    xform.translation[1],
                                                    xform.translation[2]);
  fprintf (stderr, ">>> rotation    = %f %f %f \n", xform.angles[0],
                                                    xform.angles[1],
                                                    xform.angles[2]);
  */

  if (xform.set) {
    cyl->setXform(xform);
    sphere->setXform(xform);
    }
#endif
  }

//*============================================================*
//*===========  u n i v e r s a l  j o i n t  =================*
//*============================================================*

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create universal joint.

void
PmGraphicsUniversalJoint::create() 
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  if (gc->geom) return;
  GrScene *scene = gc->scene;
  float length, radius, radii[2];

  GrCylinder *cyl = new GrCylinder (name, 2, vertices);
  cyl->setAxis (axes);

  length = attributes.scale;
  cyl->setLength(length);
  cyl->setRadius(radius);

  radii[0] = attributes.scale / 10.0;
  radii[1] = attributes.scale / 10.0;
  cyl->setRadius(radii);
  scene->addGeometry (cyl);
  gc->geom = dynamic_cast<GrVertexGeometry*>(cyl);

  GraphicsContext *cgc = (GraphicsContext*)aux_context;
  GrSphere *sphere = new GrSphere (name, 1, vertices);
  sphere->setRadius(0.01);
  sphere->setActive(attributes.visible);
  scene->addGeometry (sphere);
  cgc->geom = dynamic_cast<GrVertexGeometry*>(sphere);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display universal joint.

void
PmGraphicsUniversalJoint::display()  
  {
#ifdef USE_GRAPHICS
  GrColor colors[2];

  create();

  GraphicsContext *gc = (GraphicsContext*)context;
  GrCylinder *cyl = dynamic_cast<GrCylinder*>(gc->geom);
  colors[0] = attributes.color;
  colors[1] = attributes.color;
  //colors[1][0] = attributes.color[0] - 0.5;
  //colors[1][1] = attributes.color[1] - 0.5;
  //colors[1][2] = attributes.color[2] - 0.5;

  for (int i = 0; i < 3; i++) {
    if (colors[1][i] < 0.0) colors[1][i] = 0.0; 
    }

  cyl->setColor(attributes.color);
  cyl->setColors(2, colors);
  cyl->setDisplayType(gr_display_map[attributes.display_type]);
  cyl->setShadingType(gr_shading_map[attributes.shading_type]);
  cyl->setActive(attributes.visible);

  GraphicsContext *cgc = (GraphicsContext*)aux_context;
  GrSphere *sphere = dynamic_cast<GrSphere*>(cgc->geom);
  sphere->setColor(attributes.color);
  sphere->setDisplayType(gr_display_map[attributes.display_type]);
  sphere->setShadingType(gr_shading_map[attributes.shading_type]);
  sphere->setRadius (attributes.scale / 8.0);
  sphere->setActive(attributes.visible);
#endif
  }

//*============================================================*
//*==========              xform                     ==========*
//*============================================================*
// xform a universal joint.

void 
PmGraphicsUniversalJoint::xform(PmXform& xform) 
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GraphicsContext *cgc = (GraphicsContext*)aux_context;
  GrCylinder *cyl = dynamic_cast<GrCylinder*>(gc->geom);
  GrSphere *sphere = dynamic_cast<GrSphere*>(cgc->geom);

  if (xform.set) {
    cyl->setXform(xform);
    sphere->setXform(xform);
    }
#endif
  }

///////////////////////////////////////////////////////////////
//            G r a p h i c s    B o n d s                  //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsBonds::PmGraphicsBonds(const string name, int num_verts, PmVector3 *verts, 
                                 int num_bonds, PmConn *bonds, PmVector3 *colors)
  {
#ifdef USE_GRAPHICS
  this->name = name;
  this->num_vertices = num_verts;
  this->vertices = verts;
  this->num_bonds = num_bonds;
  this->bonds = bonds;
  this->atom_geom = NULL;
  this->colors = colors;

  atom_color_type = PM_ATOM_COLOR_RGB;
  bond_color_type = PM_ATOM_COLOR_RGB;
  atom_render = PM_MOLECULE_RENDER_SOLID;
  atom_color.set(1,1,1);
  bond_color.set(0.6, 0.6, 0.6);
  //bond_color.set(1,1,1);
  bond_atoms = false;
  base_planes = false;
  cbonly = false;
  //iatom_radius = 0.04;
  attributes.scale = 0.04;

  pmSystem.getGraphicsContext (&this->context);
  pmSystem.getGraphicsContext (&this->sphere_context);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display bonds.    

void
PmGraphicsBonds::display()
  {
#ifdef USE_GRAPHICS
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsBonds::display \n");
  fprintf (stderr, "   >>> name   [%s] \n", this->name.c_str());
  fprintf (stderr, "   >>> colors [%x] \n", colors);
  fprintf (stderr, "   >>> bond color (%f %f %f)  \n", bond_color[0], bond_color[1],
           bond_color[2]);
  fprintf (stderr, "   >>> bond atoms [%d] \n", bond_atoms);
  fprintf (stderr, "   >>> bond color [%d] \n", bond_color_type);
  */

  GraphicsContext *gc = (GraphicsContext*)context;
  GraphicsContext *sgc = (GraphicsContext*)sphere_context;
  GrScene *scene = gc->scene;
  GrLine *line;

  // create geometry //

  if (!gc->geom && !gc->pmgr_geom) {
    if (num_bonds) {
      line = new GrLine (name, num_bonds, *bonds, num_vertices, vertices);
      line->setColor(bond_color);
      line->setWidth(attributes.line_width);
      scene->addGeometry (line);
      gc->geom = dynamic_cast<GrVertexGeometry*>(line);
      }
    else {
      line = NULL;
      }

    PmGraphicsAttributes atts;
    string aname = name + "_atoms";
    atom_geom = new PmGraphicsAtoms(aname, num_vertices, vertices, colors);
    atts.setColor(atom_color);
    atts.setVisible(false);
    atom_geom->setAttributes(attributes);
    gc->pmgr_geom = atom_geom;

    if (line && colors && (bond_color_type == PM_ATOM_COLOR_ELEMENT)) {
      line->setVertexColors(num_vertices, colors);
      }
    }
  else {
    line = (GrLine*)gc->geom;
    atom_geom = (PmGraphicsAtoms*)gc->pmgr_geom;
    }

  //line->setWidth(attributes.line_width);

  //===== set atom attributes =====//

  if (bond_atoms) { 
    PmGraphicsAttributes atts;
    //atts.setScale(atom_radius);
    atts.setVisible(true);
    atts.setColor(atom_color);
    atts.setScale(attributes.scale);
    atom_geom->setAttributes(atts);
    atom_geom->setRender(atom_render);
    atom_geom->display();
    }
  else {
    PmGraphicsAttributes atts;
    atts.setVisible(false);
    atom_geom->setRender(atom_render);
    atom_geom->setAttributes(atts);
    atom_geom->display();
    }

  //===== xform lines and atoms =====//

  if (xform.set) {
    if (line) { 
      line->setXform(xform);
      }

    atom_geom->setXform(xform);
    atom_geom->display();
    }
#endif
  }


///////////////////////////////////////////////////////////////
//            G r a p h i c s    P l a n e s                //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsPlanes::PmGraphicsPlanes(const string name, int num_verts, PmVector3 *verts, 
                                   int num_planes, PmConn *conn, PmVector3 *colors)
  {
#ifdef USE_GRAPHICS
  this->name = name;
  this->num_vertices = num_verts;
  this->vertices = verts;
  this->num_planes = num_planes;
  this->connectivity = conn;
  this->atom_geom = NULL;
  this->colors = colors;

  atom_color_type = PM_ATOM_COLOR_RGB;
  bond_color_type = PM_ATOM_COLOR_RGB;
  atom_render = PM_MOLECULE_RENDER_SOLID;
  atom_color.set(1,1,1);
  bond_color.set(1,1,1);
  bond_atoms = false;
  base_planes = false;
  cbonly = false;
  atom_radius = 0.02;
  atom_radius = 0.04;

  pmSystem.getGraphicsContext (&this->context);
  pmSystem.getGraphicsContext (&this->sphere_context);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display planes.   

void
PmGraphicsPlanes::display()
  {
#ifdef USE_GRAPHICS
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsPlanes::display \n");
  fprintf (stderr, "   >>> name   [%s] \n", this->name.c_str());
  */
  GraphicsContext *gc = (GraphicsContext*)context;
  GrScene *scene = gc->scene;
  GrPolygon *poly;

  // create geometry //

  if (!gc->geom) {
    poly = new GrPolygon (name, num_planes, *connectivity, num_vertices, vertices, 0);
    scene->addGeometry (poly);
    gc->geom = dynamic_cast<GrVertexGeometry*>(poly);
    }
  else {
    poly = (GrPolygon*)gc->geom;
    }

  poly->setColor(attributes.color);

  // xform lines and atoms //

  if (xform.set) {
    poly->setXform(xform);
    }
#endif
  }

///////////////////////////////////////////////////////////////
//                    c y l i n d e r                       //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsCylinder::PmGraphicsCylinder(const string name, PmVector3 center, 
                                       float radius, float length, PmVector3 axis)
  {
#ifdef USE_GRAPHICS
  this->name = name;
  this->center = center;
  this->radius = radius;
  this->length = length;
  this->axis = axis;
  pmSystem.getGraphicsContext (&this->context);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display cylinder.  

void
PmGraphicsCylinder::display() {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GrCylinder *cyl;
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsCylinder::display \n");
  fprintf (stderr, "   >>> radius[%f] \n", radius);
  fprintf (stderr, "   >>> length[%f] \n", length);
  */

  if (!gc->geom) {
    GrScene *scene = gc->scene;
    cyl = new GrCylinder (name, 1, &center);
    cyl->setAxis(axis);
    cyl->setLength(length);
    cyl->setRadius(radius);
    scene->addGeometry(cyl);
    gc->geom = dynamic_cast<GrVertexGeometry*>(cyl);
    }
  else {
    cyl = (GrCylinder*)gc->geom; 
    }

  cyl->setColor(attributes.color);
  cyl->setShadingType(gr_shading_map[attributes.shading_type]);
  cyl->setDisplayType(gr_display_map[attributes.display_type]);

  if (xform.set) {
    cyl->setXform(xform);
    }
#endif
  }

/////////////////////////////////////////////////////////////////
//                       m o d e l                            //
///////////////////////////////////////////////////////////////

//*============================================================*
//*=================    model body     ========================*
//*============================================================*

PmGraphicsModelBody::PmGraphicsModelBody(const string name, int num_coords, 
                                         PmVector3 *coords)
  {
  this->name = name;
  this->num_vertices = num_coords;
  this->vertices = coords;
  pmSystem.getGraphicsContext (&this->context);
  pmSystem.getGraphicsContext (&this->sphere_context);
#ifdef USE_GRAPHICS
#endif
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create model body geometry.

void
PmGraphicsModelBody::create()
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  if (gc->geom) return;
  GrScene *scene = gc->scene;
  GrPoint *pts = new GrPoint(name, num_vertices, vertices);
  pts->setMarkerType (GR_GEOMETRY_MARKER_CROSS);
  pts->setColor(attributes.color);
  scene->addGeometry (pts);
  gc->geom = dynamic_cast<GrVertexGeometry*>(pts);

  gc = (GraphicsContext*)sphere_context;
  GrSphere *sphere = new GrSphere (name, num_vertices, vertices);
  PmVector3 sphere_color(0.7, 0.7, 0.7);
  sphere->setColor(sphere_color);
  sphere->setDisplayType(GR_GEOMETRY_DISPLAY_LINE);
  scene->addGeometry(sphere);
  gc->geom = dynamic_cast<GrVertexGeometry*>(sphere);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display model body as spheres.

void
PmGraphicsModelBody::display()
  {
#ifdef USE_GRAPHICS
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsModelBody::display \n");
  fprintf (stderr, "   >>> name [%s] \n", this->name.c_str());
  fprintf (stderr, "   >>> msize [%g] \n", this->msize);
  */
  GraphicsContext *gc = (GraphicsContext*)context;
  create();

  GrPoint *pts = dynamic_cast<GrPoint*>(gc->geom);
  pts->setSize (attributes.scale);
  pts->setMarkerSize(attributes.scale);

  gc = (GraphicsContext*)sphere_context;
  GrSphere *sphere = dynamic_cast<GrSphere*>(gc->geom);
  sphere->setRadius (attributes.scale);
#endif
  }

//*============================================================*
//*=================    model joint    ========================*
//*============================================================*

PmGraphicsModelJoint::PmGraphicsModelJoint(const string name, int num_coords, 
                                           PmVector3 *coords)
  {
  this->name = name;
  this->num_vertices = num_coords;
  this->vertices = coords;
  pmSystem.getGraphicsContext (&this->context);
  attributes.color.set(1,0,1);
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create model joints geometry.

void
PmGraphicsModelJoint::create()
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  if (gc->geom) return;
  GrScene *scene = gc->scene;
  GrPoint *pts = new GrPoint(name, num_vertices, vertices);
  pts->setMarkerType (GR_GEOMETRY_MARKER_CROSS);
  pts->setColor(attributes.color);
  scene->addGeometry (pts);
  gc->geom = dynamic_cast<GrVertexGeometry*>(pts);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display model joints as a cross marker for points.

void
PmGraphicsModelJoint::display()
  {
#ifdef USE_GRAPHICS
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsModelJoint::display \n");
  fprintf (stderr, "   >>> name [%s] \n", this->name.c_str());
  fprintf (stderr, "   >>> msize [%g] \n", this->msize);
  */
  GraphicsContext *gc = (GraphicsContext*)context;
  create();
  GrPoint *pts = dynamic_cast<GrPoint*>(gc->geom);
  pts->setSize (attributes.scale);
  pts->setMarkerSize(attributes.scale);
#endif
  }

///////////////////////////////////////////////////////////////
//                       L i n e                            //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsLine::PmGraphicsLine(const string name, int num_verts, PmVector3 *verts, 
                               bool copy)
  {
  this->name = name;
  this->num_vertices = num_verts;
  this->vertices = verts;
  this->copy_vertices = copy;
  pmSystem.getGraphicsContext (&this->context);
  attributes.color.set(1,0,1);
  }

//*============================================================*
//*==========              setDisjoint               ==========*
//*============================================================*

void
PmGraphicsLine::setDisjoint(const bool flag)
  {
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display lines.  

void
PmGraphicsLine::display()
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GrScene *scene = gc->scene;
  GrLine *line;

  if (!gc->geom) {
    line = new GrLine(name, num_vertices, vertices, this->copy_vertices);

    if (attributes.disjoint) { 
      line->setDisjoint(true);
      }

    scene->addGeometry(line);
    gc->geom = dynamic_cast<GrVertexGeometry*>(line);
    }
  else {
    line = (GrLine*)gc->geom;
    }

  line->setDisplayType(gr_display_map[attributes.display_type]);
  line->setColor(attributes.color);
  line->setWidth(attributes.line_width);
  line->setActive(attributes.visible);

  if (colors) {
    line->setColors(num_colors, colors);
    }

  if (xform.set) {
    line->setXform(xform);
    }

#endif
  }

//*============================================================*
//*==========              update                    ==========*
//*============================================================*
// update lines.

/*
void
PmGraphicsLine::update(int num_verts, PmVector3 *vertices)
  {
  }
*/


///////////////////////////////////////////////////////////////
//                       P o i n t s                        //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsPoint::PmGraphicsPoint(const string name, int num_verts, PmVector3 *verts)
  {
  this->name = name;
  this->num_vertices = num_verts;
  this->vertices = verts;
  pmSystem.getGraphicsContext (&this->context);
  attributes.color.set(1,0,1);
  attributes.scale = 4.0;
  attributes.scale = 8.0;
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display points. 

void
PmGraphicsPoint::display()
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GrScene *scene = gc->scene;
  GrPoint *point;

  if (!gc->geom) { 
    point = new GrPoint(name, num_vertices, vertices);
    point->setColor(attributes.color);
    point->setSize(attributes.scale);
    point->setSize(2.0);
    scene->addGeometry(point);
    gc->geom = dynamic_cast<GrVertexGeometry*>(point);

    if (attributes.marker) {
      point->setMarkerType (GR_GEOMETRY_MARKER_CROSS);
      point->setSize(attributes.scale);
      }
    }
  else {
    point = (GrPoint*)gc->geom; 
    }

  if (colors) {
    point->setColors(num_vertices, colors);
    }

  if (xform.set) {
    point->setXform(xform);
    }

  point->setActive(attributes.visible);
#endif
  }

///////////////////////////////////////////////////////////////
//                       p o l y g o n                      //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsPolygon::PmGraphicsPolygon(const string name, int num_poly, PmConn *conn, 
                                     int hgeom, int num_verts, PmVector3 *verts)
  {
  this->name = name;
  this->num_polygons = num_poly;
  this->nodes_per_poly = hgeom;
  this->connectivity = conn;
  this->num_vertices = num_verts;
  this->vertices = verts;
  pmSystem.getGraphicsContext (&this->context);
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display polygons.

void
PmGraphicsPolygon::display()
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GrScene *scene = gc->scene;
  GrPolygon *poly;

  if (!gc->geom) { 
    poly = new GrPolygon (name, num_polygons, *connectivity, num_vertices, vertices, 
                          nodes_per_poly);
    scene->addGeometry(poly);
    gc->geom = dynamic_cast<GrVertexGeometry*>(poly);
    }
  else {
    poly = dynamic_cast<GrPolygon*>(gc->geom);
    }

  if (colors) {
    poly->setVertexColors(num_vertices, colors);
    }

  poly->setColor(attributes.color);
  poly->setDisplayType(gr_display_map[attributes.display_type]);
  poly->setShadingType(gr_shading_map[attributes.shading_type]);

  if (xform.set) {
    poly->setXform(xform);
    }

  poly->setActive(attributes.visible);
#endif
  }

///////////////////////////////////////////////////////////////
//                       r e s t r a i n t                  //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsRestraint::PmGraphicsRestraint(const string name, PmVector3& pt1, 
                                         PmVector3& pt2) 
  {
#ifdef USE_GRAPHICS
  //fprintf (stderr, "\n>>>>>> PmGraphicsRestraint::ctor name[%s] \n", name.c_str());
  this->name = name;
  vertices = new PmVector3[2];
  vertices[0] = pt1; 
  vertices[1] = pt2; 
  num_vertices = 2;
  pmSystem.getGraphicsContext (&this->context);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display a restraint as a line connectined two points.

void
PmGraphicsRestraint::display()
  {
#ifdef USE_GRAPHICS
  //fprintf (stderr, "\n>>>>>> PmGraphicsRestraint::display\n");
  GraphicsContext *gc = (GraphicsContext*)context;
  GrScene *scene = gc->scene;
  GrLine *line;

  if (!gc->geom) {
    line = new GrLine (name, num_vertices, vertices);
    scene->addGeometry(line);
    gc->geom = dynamic_cast<GrVertexGeometry*>(line);
    }
  else {
    line = (GrLine*)gc->geom;
    }

  line->setColor(attributes.color);
  line->setWidth(attributes.line_width);
#endif
  }

//*============================================================*
//*==========               update                   ==========*
//*============================================================*
// update a restraint with new end points.

void
PmGraphicsRestraint::update(PmVector3& pt1, PmVector3& pt2)
  {
#ifdef USE_GRAPHICS
  vertices[0] = pt1;
  vertices[1] = pt2;
  GraphicsContext *gc = (GraphicsContext*)context;
  GrLine *line = (GrLine*)gc->geom;
  line->setVertices(num_vertices, vertices);
#endif
  }

///////////////////////////////////////////////////////////////
//                       S p h e r e                        //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsSphere::PmGraphicsSphere(const string name, int num_verts, PmVector3 *verts,
                                   float *radius)
  {
  this->name = name;
  this->num_vertices = num_verts;
  this->vertices = verts;
  this->radius = radius;
  pmSystem.getGraphicsContext (&this->context);
  attributes.color.set(1,0,1);
  attributes.display_type = PM_GEOMETRY_DISPLAY_LINE;
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display spheres.

void
PmGraphicsSphere::display()
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  GrScene *scene = gc->scene;
  GrSphere *sphere; 

  if (!gc->geom) { 
    sphere = new GrSphere (name, num_vertices, vertices);
    sphere->setColor(attributes.color);
    //sphere->setRadius (0.02);
    sphere->setRadius (radius);
    sphere->setDisplayType(gr_display_map[attributes.display_type]);
    sphere->setShadingType(gr_shading_map[attributes.shading_type]);
    sphere->setLineWidth(attributes.line_width);
    sphere->setResolution(32);

    if (colors) {
      sphere->setColors (num_vertices, colors);
      }

    scene->addGeometry (sphere);
    gc->geom = dynamic_cast<GrVertexGeometry*>(sphere);

    /*
    sphere1 = new GrSphere (name, num_vertices, vertices);
    color.set(1,0,0);
    sphere1->setColor(color);
    sphere1->setRadius (radius);
    sphere1->setOptimized(false);
    sphere1->setShadingType(gr_shading_map[PM_GEOMETRY_SHADING_COLOR]);
    scene->addGeometry (sphere1);
    gc->geom1 = dynamic_cast<GrVertexGeometry*>(sphere1);
    */
    }
  else {
    sphere = (GrSphere*)gc->geom; 
    //sphere1 = (GrSphere*)gc->geom1; 
    }

  if (xform.set) {
    sphere->setXform(xform);
    //sphere1->setXform(xform);
    }

  //sphere->setOptimized(false);

#endif
  }

///////////////////////////////////////////////////////////////
//                       S u r f a c e                      //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========      constructors / destructor         ==========*
//*============================================================*

PmGraphicsSurface::PmGraphicsSurface(const string name, int num_verts, PmVector3 *verts,
                                     int num_polys, PmConn *conn, int hgeom)
  {
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsSurface:: ctor \n");
  fprintf (stderr, "   >>> name [%s] \n", this->name.c_str());
  */

  this->name = name;
  this->num_vertices = num_verts;
  this->vertices = verts;
  this->num_polygons = num_polys;
  this->connectivity = conn;
  this->hgeom = hgeom;
  pmSystem.getGraphicsContext (&this->context);
  }

//*============================================================*
//*==========              create                    ==========*
//*============================================================*
// create surface geometry.

void
PmGraphicsSurface::create()
  {
#ifdef USE_GRAPHICS
  GraphicsContext *gc = (GraphicsContext*)context;
  if (gc->geom) return;
  GrScene *scene = gc->scene;
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsSurface::create \n");
  fprintf (stderr, "   >>> name [%s] \n", this->name.c_str());
  fprintf (stderr, "   >>> hgeom[%d] \n", this->hgeom);
  */

  GrPolygon *poly = new GrPolygon (name, num_polygons, *connectivity, num_vertices,
                                   vertices, hgeom);
  poly->setColor(attributes.color);
  poly->setDisplayType(gr_display_map[attributes.display_type]);
  poly->setShadingType(gr_shading_map[attributes.shading_type]);
  poly->setLighting(attributes.twosided_lighting);

  if (colors) {
    poly->setVertexColors(num_vertices, colors);
    }

  scene->addGeometry(poly);
  gc->geom = dynamic_cast<GrVertexGeometry*>(poly);
#endif
  }

//*============================================================*
//*==========              display                   ==========*
//*============================================================*
// display surface. 

void
PmGraphicsSurface::display()
  {
#ifdef USE_GRAPHICS
  /*
  fprintf (stderr, "\n>>>>>> PmGraphicsSurface::display \n");
  fprintf (stderr, "   >>> name   [%s] \n", this->name.c_str());
  */
  GraphicsContext *gc = (GraphicsContext*)context;
  create();

  GrPolygon *poly = (GrPolygon*)gc->geom;
  poly->setColor(attributes.color);
  poly->setDisplayType(gr_display_map[attributes.display_type]);
  poly->setShadingType(gr_shading_map[attributes.shading_type]);
  poly->setLighting(attributes.twosided_lighting);

  if (xform.set) {
    poly->setXform(xform);
    }
#endif
  }


}
