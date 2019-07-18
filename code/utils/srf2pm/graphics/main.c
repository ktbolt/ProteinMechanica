/*----------------------------------------------------------------------*
 *                                                                      *
 *                        ****  gr test  ****                           *
 *                                                                      *
 *----------------------------------------------------------------------*/

#include "dm.h"
#include "gr.h"
/*
#include "/usr/include/GL/gl.h"
#include "glut.h"
*/


DmObj *gwin, *gscene;

void 
GrExtentSet (DmObj *geom);

void
show_poly (int p, DmObj *geom);

void
GrPickShow (DmPoint3 pt);

void
GrPickInit (void (*proc_func)());

void
GrPickProc (DmObj *scene, GrPickResult *pick_res);


/*----------------------------------------------------------------------*
 *                                                                      *
 *                          ****  main  ****                            *
 *                                                                      *
 *----------------------------------------------------------------------*/

main ()
  {

  int win_x, win_y;

  int win_width, win_height;

  char *win_name;

  int type;

  DmExtent extent;

  char *scene_name;

  DmObj *scene;

  int rep;

  int r;

  float t;

  float dt;

  int *conn;

  DmObj *win;

  GrGeometryType geom_type; 

  GrGeometryOptions geom_opts;

  DmObj *geom;

  int num_verts;

  DmPoint3 *verts;

  int num_polys;

  int i;

  DmPoint3 origin;

  float dims[3];

  float radius, height;

  DmObj *light;

  GrSystemType system;

  GrRenderMode render;

  DmObj *sys;

  int status;

  GrColor color;

 /**************
  ***  body  ***
  **************/

  gr_Init ();

  system = GR_SYSTEM_OPENGL;
  render = GR_RENDER_MODE_SURFACE;

  gr_SystemGet (system, render, &sys);


  /*  create and open a window.  */

  win_name = "test";
  win_x = 0;
  win_y = 0;
  win_width = 500;
  win_height = 500;

  type = GR_WINDOW_TYPE_DOUBLE | GR_WINDOW_TYPE_RGB | GR_WINDOW_TYPE_DEPTH |
         GR_WINDOW_TYPE_DIRECT;

  gr_WinCreate (win_name, win_x, win_y, win_width, win_height, type, sys, &win);
  gr_WinOpen (win);


  /*  set the coordinates and viewport.  */

  dm_ExtentSet (extent, 0.0, 10.0, 0.0, 10.0, -20.0, 20.0);
  gr_WinExtentSet (win, extent);


  /*  create a scene and set it for the window.  */

  scene_name = "simple test";
  gr_SceneCreate (scene_name, &scene);
  gr_WinSceneSet (win, scene);
  gr_SceneExtentSet (scene, extent);
  gr_ViewportSet (win, sys, scene, win_width, win_height, extent);


  /*  create a light.  */

  gr_LightCreate ("light 1", GR_LIGHT_DIRECTIONAL, &light);
  gr_SceneLightAdd (sys, scene, light);


  /*  read surf  */ 

  read ("t.srf", &geom);

  gr_SceneGeomAdd (scene, geom);

  gwin = win;
  gscene = scene;

  GrExtentSet (geom);

  PolyReorient (geom);

  show_poly (3657, geom);
  /*
  show_poly (4070, geom);
  show_poly (4199, geom);
  show_poly (1138, geom);
  show_poly (1139, geom);
  show_poly (1160, geom);
  show_poly (1161, geom);
  show_poly (1910, geom);
  show_poly (1911, geom);
  show_poly (1909, geom);
  */

  GrPickInit (GrPickProc); 


  /*  render the scene.  */

  gr_SceneRender (win, scene, GR_RENDER_MODE_SURFACE);
  gr_WinEventProcLoop (win, 0, DM_TRUE, "cmd>");
  }


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  GrExtentSet  ****                 *
 *                                                            *
 *------------------------------------------------------------*/

void 
GrExtentSet (DmObj *geom)
  {

  DmExtent extent;

  float xmin, xmax, ymin, ymax, zmin, zmax;

  float dx, dy, dz, max_dim;

  float cx, cy, cz;

  float x, y, z;

  DmExtent gr_extent;

  int num_verts; 

  DmPoint3 *verts;

  int num_polys, *conn;

  int i;

 /**************
  ***  body  ***
  **************/

  gr_GeomIndexPolygonsGet (geom, &num_verts, &verts, &num_polys, &conn);

  for (i = 0; i < num_verts; i++) {
    x = verts[i][0];
    y = verts[i][1];
    z = verts[i][2];

    if (i == 0) {
      xmin = xmax = x;
      ymin = ymax = y;
      zmin = zmax = z;
      }
    else {
      xmin = (x < xmin ? x : xmin);
      ymin = (y < ymin ? y : ymin);
      zmin = (z < zmin ? z : zmin);
      xmax = (x > xmax ? x : xmax);
      ymax = (y > ymax ? y : ymax);
      zmax = (z > zmax ? z : zmax);
      }
    }

  /*
  fprintf (stderr, " ------  ExtentSet  ------- \n");
  fprintf (stderr, " >>>>>>  xmin [%f]  xmax [%f] \n", xmin, xmax);
  fprintf (stderr, " >>>>>>  ymin [%f]  ymax [%f] \n", ymin, ymax);
  fprintf (stderr, " >>>>>>  zmin [%f]  zmax [%f] \n", zmin, zmax);
  */

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
  zmin = cz - 4.0*max_dim;
  zmax = cz + 4.0*max_dim;
  zmin = cz - 10.0*max_dim;
  zmax = cz + 10.0*max_dim;

  dm_ExtentSet (gr_extent, xmin, xmax, ymin, ymax, zmin, zmax);
  gr_WinExtentSet (gwin, gr_extent);
  gr_SceneExtentSet (gscene, gr_extent);
  gr_WinSet (gwin);
  }



/*------------------------------------------------------------*
 *                                                            *
 *                      ****  read  ****                      *
 *                                                            *
 *------------------------------------------------------------*/

void
show_poly (int p, DmObj *geom)
  {

  GrGeometryType geom_type;

  GrGeometryOptions geom_opts;

  DmObj *ngeom;

  float color[3];

  int num_verts;

  DmPoint3 *verts;

  DmPoint3 *nverts;

  int num_polys;

  int *conn, *pconn;

  char name[80];

  int i, j, k, n;

 /**************
  ***  body  ***
  **************/

  gr_GeomIndexPolygonsGet (geom, &num_verts, &verts, &num_polys, &conn);

  geom_type = GR_GEOMETRY_LINES;
  geom_opts = GR_GEOMETRY_OPTION_NONE;

  sprintf (name, "poly:line:%d", p);
  gr_GeomCreate (name, geom_type, geom_opts, &ngeom);
  n = 8;
  nverts = dm_MemAlloc (0, DmPoint3, 2*n);
  n = 0;
  fprintf (stderr, "\n\n");

  for (i = 0; i < 3; i++) {
    j = conn[3*p+i];
    dm_PointCopy3 (nverts[n], verts[j]);
    fprintf (stderr, " %f %f %f \n", dm_Plist3(nverts[n]));
    n += 1;

    if (i < 2) {
      j = i + 1;
      }
    else {
      j = 0;
      }

    k = conn[3*p+j];
    dm_PointCopy3 (nverts[n], verts[k]);
    fprintf (stderr, " %f %f %f \n", dm_Plist3(nverts[n]));
    n += 1;
    }

  dm_PointCopy3 (nverts[n], verts[k]);
  n += 1;
  dm_PointSet3 (nverts[n], 0, 0, 0);

  gr_GeomLinesSet (ngeom, n, nverts);

  gr_GeomLineWidthSet (ngeom, 4.0);

  color[0] = 1.0, color[1] = 0.0, color[2] = 0.0;

  gr_GeomColorSet (ngeom, color);

  gr_SceneGeomAdd (gscene, ngeom);
  }


void
GrPickShow (DmPoint3 pt)
  {

  char geom_name[80];

  GrGeometryType geom_type;

  GrGeometryOptions geom_opts;

  DmObj *gr_geom;

  int n;

  DmPoint3 *pts;

  GrColor color;

  GrGeometryMarkerType mtype;

 /**************
  ***  body  ***
  **************/

  /*  show picked point.  */

  sprintf (geom_name, "pick");
  geom_type = GR_GEOMETRY_MARKERS;
  geom_opts = GR_GEOMETRY_OPTION_NONE;
  gr_GeomGet (geom_name, &gr_geom);

  if (!gr_geom) {
    pts = (DmPoint3*)malloc(sizeof(DmPoint3) * 1);
    gr_GeomCreate (geom_name, geom_type, geom_opts, &gr_geom);
    color[0] = 1.0,  color[1] = 1.0, color[2] = 1.0;
    gr_GeomColorSet (gr_geom, color);
    gr_GeomMarkersScaleSet (gr_geom, 1.5);
    gr_SceneGeomAdd (gscene, gr_geom);
    }
  else {
    gr_GeomMarkersGet (gr_geom, &mtype, &n, &pts);
    }

  pts[0][0] = pt[0];
  pts[0][1] = pt[1];
  pts[0][2] = pt[2];

  gr_GeomMarkersSet (gr_geom, GR_GEOMETRY_MARKER_CROSS, 1, pts);
  gr_SceneRender (gwin, gscene, GR_RENDER_MODE_SURFACE);
  }


/*------------------------------------------------------------*
 *                                                            *
 *                  ****  pm_GrPickInit  ****                 *
 *                                                            *
 *------------------------------------------------------------*/

void
GrPickInit (void (*proc_func)())
  {

  GrPickAtts pick_atts;

 /**************
  ***  body  ***
  **************/

  pick_atts.type = GR_PICK_POINT;
  pick_atts.tol = 1.0e-01;
  pick_atts.tol = 1.0;
  gr_ScenePickAttsSet (gscene, &pick_atts);
  gr_ScenePickActiveSet (gscene, DM_TRUE);
  /*
  gr_ScenePickRequest (man.scene, 1, GR_PICK_VERTEX, proc_func);
  */
  gr_ScenePickRequest (gscene, 1, GR_PICK_POINT, proc_func);
  }



void
GrPickProc (DmObj *scene, GrPickResult *pick_res)
  {

  GrPickGeom *pick_list;

  DmObj *geom;

  int i;

  DmObj *body;

  int num;

  char *s, str[80], body_name[80];

  int entity;

 /**************
  ***  body  ***
  **************/

  /*
  fprintf (stderr, "\n---------- dyn_GrBodyPickProc ---------- \n");
  */
  pick_list = pick_res->geom_list;

  while (pick_list) {
    geom = pick_list->geom;
    GrPickShow (pick_list->int_pt);
    break;
    /*
    fprintf (stderr, ">>>>>> geom name [%s] \n", geom->name);
    */
    pick_list = pick_list->next;
    }

  gr_GeomPickEntityGet (geom, &entity);
  fprintf (stderr, ">>>>>> ent [%d]. \n", entity);
  }
