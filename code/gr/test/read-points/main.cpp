
#include "gr.h"
#include "win.h"
#include "scene.h"
#include "light.h"
#include "geom.h"
#include "sphere.h"
#include "point.h"
#include "line.h"
#include "poly.h"

using namespace PmGraphics;

int
main (int nargs, char **args) 
  {

  //-------------- read data --------------//
  GrVector3 verts[200000];
  int n = 0;
  char *s, fline[100];
  float x, y, z;
  int j, index;

  //char *fn = "t.dat";
  char *fn = "t1.dat";
  FILE *fp = fopen(fn, "r");
  fgets(fline, 2000, fp);
  int num_verts = atoi(fline);
  printf( "num verts: %d\n", num_verts);
  GrExtent vextent(1e6,-1e6,1e6,-1e6,1e6,-1e6);

  n = 0;
  for (int i = 0; i < num_verts; i++) {
     fgets(fline, 2000, fp);
     sscanf(fline, "%d: %f %f %f \n", &n, &x, &y, &z);
     //printf ( "%d: %f %f %f \n", n, x, y, z);
     verts[n].set(x, y, z);
     vextent.update( verts[n] );
     n += 1;
  }

  printf ( "extent: %f %f %f \n", vextent.min[0], vextent.min[1], vextent.min[2] );
  printf ( "        %f %f %f \n", vextent.max[0], vextent.max[1], vextent.max[2] );

/*
  fgets(fline, 2000, fp);
  int num_index = atoi(fline);
  int *indexes = new int[100000];
  printf( "num index: %d\n", num_index);
  n = 0;
  for (int i = 0; i < num_index; i++) {
     fgets(fline, 2000, fp);
     sscanf(fline, "%d: %d  \n", &j, &index );
     //printf ( "%d: %d  \n", j, index);
     indexes[n] = index;
     n += 1;
  }

*/

  int wtype;

  //---------- create a window ----------//
  char *win_name = "test";
  int win_x = 0;
  int win_y = 0;
  int win_width = 1000;
  int win_height = 1000;

  wtype = GR_WINDOW_TYPE_DOUBLE | GR_WINDOW_TYPE_RGB | GR_WINDOW_TYPE_DEPTH |
          GR_WINDOW_TYPE_DIRECT;

  GrWindow *win = new GrWindow (win_name, win_x, win_y, win_width, win_height, wtype);
  //GrExtent extent(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0);
  GrExtent extent(-5.0, 5.0, -8.0, -2.0, -100.0, 100.0);
  win->setExtent(extent);
  win->open();

  // create a scene //

  GrScene *scene = new GrScene ("scene");
  scene->setWindow(win);
  scene->setViewport(win_width, win_height, extent);
  scene->setExtent(extent);

  // create a light 
  GrLight *light = new GrLight ("light 1", GR_LIGHT_DIRECTIONAL);
  scene->addLight (light);

  // add geom //
  GrPoint *geom_pts = new GrPoint("verts", num_verts, verts);
  geom_pts->setMarkerSize( 2.0 );
  scene->addGeometry(geom_pts);

  GrSphere *sp1 = new GrSphere ("sphere1", num_verts, verts);
  sp1->setShadingType (GR_GEOMETRY_SHADING_NONE);
  sp1->setShadingType (GR_GEOMETRY_SHADING_COLOR);
  sp1->setColor(1,0,0);
  sp1->setRadius (0.05);
  sp1->setResolution(32);
  //scene->addGeometry (sp1);

/*
  GrLine *line = new GrLine("line", num_verts, verts);
  line->setColor(0,1,0);
  scene->addGeometry(line);

  int num_poly = num_index / 3;
  num_poly = 36; 
  GrIndex conn(4*num_poly);
  n = 0;

  for (int i = 0; i < num_poly*3; i += 3 ) {
     conn[4*n] = 3;
     conn[4*n+1] = indexes[i];
     conn[4*n+2] = indexes[i+1];
     conn[4*n+3] = indexes[i+2];
     printf ( "%d: %d %d %d \n", n, conn[4*n+1], conn[4*n+2], conn[4*n+3] ); 
     n = n + 1;
  }

  GrPolygon *poly = new GrPolygon("poly", num_poly, conn, num_verts, verts);
  poly->setColor(1,1,0);
  scene->addGeometry(poly);

*/
 
  scene->render();

  // process events //
  win->setScene (scene);
  win->processEvents (NULL, false, "test>", NULL);
  }


