
#include "gr.h"
#include "win.h"
#include "scene.h"
#include "light.h"
#include "geom.h"
#include "sphere.h"

using namespace PmGraphics;

int
main (int nargs, char **args) 
  {

  int wtype;

  // create a window

  char *win_name = "test";
  int win_x = 0;
  int win_y = 0;
  int win_width = 500;
  int win_height = 500;

  wtype = GR_WINDOW_TYPE_DOUBLE | GR_WINDOW_TYPE_RGB | GR_WINDOW_TYPE_DEPTH |
          GR_WINDOW_TYPE_DIRECT;

  GrWindow *win = new GrWindow (win_name, win_x, win_y, win_width, win_height, wtype);
  GrExtent extent(-10.0, 10.0, -10.0, 10.0, -20.0, 20.0);
  win->setExtent (extent);
  win->open ();


  // create a scene

  GrScene *scene = new GrScene ("scene");
  scene->setWindow (win);
  scene->setViewport (win_width, win_height, extent);
  scene->setExtent (extent);


  // create a light 

  GrLight *light = new GrLight ("light 1", GR_LIGHT_DIRECTIONAL);
  scene->addLight (light);


  // create some spheres 

  GrVector3 verts[7];
  int n = 0;
  verts[n++].set(-6, 0, 0);
  verts[n++].set(-4, 0, 0);
  verts[n++].set(-2, 0, 0);
  verts[n++].set(0, 0, 0);
  verts[n++].set(2, 0, 0);
  verts[n++].set(4, 0, 0);
  verts[n++].set(6, 0, 0);

  GrColor colors[7];
  n = 0;
  colors[n++].set(1,0,0);
  colors[n++].set(0,1,0);
  colors[n++].set(0,0,1);
  colors[n++].set(1,1,0);
  colors[n++].set(1,0,1);
  colors[n++].set(0,0,1);
  colors[n++].set(0,1,1);

  float radii[7];
  n = 0;
  radii[n++] = 1.0;
  radii[n++] = 0.4;
  radii[n++] = 0.2;
  radii[n++] = 0.1;
  radii[n++] = 0.2;
  radii[n++] = 0.4;
  radii[n++] = 1.0;

  GrSphere *sp1 = new GrSphere ("sphere1", 7, verts);
  sp1->setShadingType (GR_GEOMETRY_SHADING_NONE);
  sp1->setShadingType (GR_GEOMETRY_SHADING_COLOR);
  sp1->setColor(1,0,0); 
  sp1->setColors(7, colors); 

  sp1->setRadius (0.2); 
  sp1->setRadius (radii); 

  scene->addGeometry (sp1);
  scene->render();

  grSystem.setCommandCallback1 ();


  // process events

  win->setScene (scene);
  win->processEvents (NULL, false, "test>", NULL);
  }


