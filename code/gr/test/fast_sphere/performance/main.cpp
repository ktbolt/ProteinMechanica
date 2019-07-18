
#include "gr.h"
#include "win.h"
#include "scene.h"
#include "light.h"
#include "geom.h"
#include "sphere.h"

using namespace PmGraphics;

float RandNumber(long *seed);

void
pm_CmdProcess (FILE *fp, const string line)
  {
  //printf (">>> proc cmd = %s \n", cmd.c_str());
  }

int
main (int nargs, char **args) 
  {

  int wtype;

  // create a window

  char *win_name = "test";
  int win_x = 0;
  int win_y = 0;
  int win_width = 400;
  int win_height = 400;

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

  //const int num = 25000;
  const int num = 500000;
  GrVector3 *verts = new GrVector3[num];
  GrColor colors[num];

  long seed;
  float dx, dy, dz;
  float x, y, z;
  float xmin, xmax, ymin, ymax, zmin, zmax;

  xmin = -6.5; xmax =  6.5;
  ymin = -6.5; ymax =  6.5;
  zmin = -6.5; zmax =  6.5;

  dx = (xmax - xmin) / 2.0;
  dy = (ymax - ymin) / 2.0;
  dz = (zmax - zmin) / 2.0;
  seed = 1357931;

  for (int i = 0; i < num; i++) {
    x = dx*(2*RandNumber(&seed) - 1.0);
    y = dy*(2*RandNumber(&seed) - 1.0);
    z = dz*(2*RandNumber(&seed) - 1.0);

    if (num == 1) {
      x = y = z = 0.0;
      }

    verts[i].set(x,y,z);
    }

  GrSphere *sp1 = new GrSphere ("sphere1", num, verts);
  sp1->setShadingType (GR_GEOMETRY_SHADING_NONE);
  sp1->setShadingType (GR_GEOMETRY_SHADING_COLOR);
  //sp1->setDisplayType (GR_GEOMETRY_DISPLAY_LINE);
  sp1->setColor(1,0,0); 

  sp1->setRadius (0.5); 
  scene->addGeometry (sp1);

   grSystem.setCommandCallback (pm_CmdProcess);


  // process events

  win->setScene (scene);
  win->processEvents (NULL, false, "test>", NULL);


  }


float 
RandNumber(long *seed)
  {

  float rnum;

  int k;

  float am;

  static long a = 16807;

  static long m = 2147483647;

  static long q = 127773;

  static long mask = 123459876;

  static long r = 2836;

 /**************
  ***  body  ***
  **************/

  am = 1.0 / m;
  *seed ^= mask;

  k = (*seed) / q;
  *seed = a * (*seed - k*q) - r*k;

  if (*seed < 0) {
    *seed += m;
    }

  rnum = am * (float)(*seed);
  *seed ^= mask;
  return (rnum);
  }
