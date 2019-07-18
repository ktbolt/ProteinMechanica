
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

  int n = 2;
  int m = 2;
  int size = n*m;
  int k = 0;

  GrVector3 verts[size];
  float dx = 4.0;
  float dy = 4.0;
  float x, y, z;

  z = 0.0;
  y = -5.0;

  for (int i = 0; i < n; i++) {
    x = -5.0;

    for (int j = 0; j < m; j++) {
      x += dx;
      verts[k++].set(x, y, z);
      }

    y += dy;
    }

  GrSphere *sp1 = new GrSphere ("sphere1", size, verts);
  sp1->setShadingType (GR_GEOMETRY_SHADING_NONE);
  sp1->setShadingType (GR_GEOMETRY_SHADING_COLOR);
  sp1->setDisplayType (GR_GEOMETRY_DISPLAY_LINE);
  sp1->setColor(1,0,0); 

  sp1->setRadius (1.15); 
  scene->addGeometry (sp1);


  GrXform xform;
  xform.set = true;
  //xform.angles[0] = 45.0;
  xform.angles[1] = 5.0;
  xform.angles[2] = 45.0;

  xform.translation[0] = 1.0;

  sp1->setXform(xform); 
  //sp1->setOptimized(false); 


  // process events

  win->setScene (scene);
  win->processEvents (NULL, false, "test>", NULL);
  }


