
#include "gr.h"
#include "win.h"
#include "scene.h"
#include "light.h"
#include "geom.h"
#include "line.h"
#include "mth.h"

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

  // create a scene //

  GrScene *scene = new GrScene ("scene");
  scene->setWindow (win);
  scene->setViewport (win_width, win_height, extent);
  scene->setExtent (extent);

  // create a light 

  GrLight *light = new GrLight ("light 1", GR_LIGHT_DIRECTIONAL);
  scene->addLight (light);

  // axes //
  GrVector3 verts[6];
  int n = 0;
  verts[n++].set(0, 0, 0); verts[n++].set(1, 0, 0);
  verts[n++].set(0, 0, 0); verts[n++].set(0, 1, 0);
  verts[n++].set(0, 0, 0); verts[n++].set(0, 0, 1);

  GrColor colors[3];
  colors[0].set (1,0,0);
  colors[1].set (0,1,0);
  colors[2].set (0,1,1);

  GrLine *line = new GrLine ("line1", 6, verts);
  line->setDisjoint(true);
  line->setColors(3, colors);
  scene->addGeometry (line);

  // axes to transform //

  PmVector3 v1, v2, v3, a1, a2, a3;
  PmMatrix3x3 mat, zmat, mat2, mat1;
  float ang1, ang2, ang3;
  GrVector3 averts[6];

  a1.set(1,1,0); 
  a1.set(0,1,1); 
  pm_MathBasisCompute (a1, a2, a3);

  n = 0;
  averts[n++].set(0, 0, 0); averts[n++] = a1; 
  averts[n++].set(0, 0, 0); averts[n++] = a2; 
  averts[n++].set(0, 0, 0); averts[n++] = a3; 

  GrLine *aline = new GrLine ("aline", 6, averts);
  aline->setDisjoint(true);

  GrColor acolors[3];
  acolors[0].set (0.5,0,0);
  acolors[1].set (0,0.5,0);
  acolors[2].set (0,0.5,0.5);
  aline->setColors(3, acolors);
  scene->addGeometry (aline);

  // compute transformation //

  a1.normalize(); a2.normalize(); a3.normalize();

  fprintf (stderr, "\n");
  fprintf (stderr, ">>> a1=(%f %f %f) \n", a1[0], a1[1], a1[2]);
  fprintf (stderr, ">>> a2=(%f %f %f) \n", a2[0], a2[1], a2[2]);
  fprintf (stderr, ">>> a3=(%f %f %f) \n", a3[0], a3[1], a3[2]);

  v1.set(1,0,0); v2.set(0,1,0); v3.set(0,0,1);
  v1.set(3,5,1); 
  v1.normalize(); 
  pm_MathBasisCompute (v1, v2, v3);
  v2.normalize(); v3.normalize();

  fprintf (stderr, "\n");
  fprintf (stderr, ">>> v1=(%f %f %f) \n", v1[0], v1[1], v1[2]);
  fprintf (stderr, ">>> v2=(%f %f %f) \n", v2[0], v2[1], v2[2]);
  fprintf (stderr, ">>> v3=(%f %f %f) \n", v3[0], v3[1], v3[2]);

  // mat rot //

  pm_MathFrameRotation(a1, a2, a3, v1, v2, v3, mat);

  PmVector3 r1, r2, r3;
  r1 = mat*a1;
  r2 = mat*a2;
  r3 = mat*a3;

  fprintf (stderr, "\n rot\n");
  fprintf (stderr, ">>> r1=(%f %f %f) \n", r1[0], r1[1], r1[2]);
  fprintf (stderr, ">>> r2=(%f %f %f) \n", r2[0], r2[1], r2[2]);
  fprintf (stderr, ">>> r3=(%f %f %f) \n", r3[0], r3[1], r3[2]);

  // dir cosines //

  mat(0,0) = v1*a1; mat(0,1) = v1*a2; mat(0,2) = v1*a3;
  mat(1,0) = v2*a1; mat(1,1) = v2*a2; mat(1,2) = v2*a3;
  mat(2,0) = v3*a1; mat(2,1) = v3*a2; mat(2,2) = v3*a3;

  mat1(0,0) = a1[0]; mat1(0,1) = a2[0]; mat1(0,2) = a3[0];
  mat1(1,0) = a1[1]; mat1(1,1) = a2[1]; mat1(1,2) = a3[1];
  mat1(2,0) = a1[2]; mat1(2,1) = a2[2]; mat1(2,2) = a3[2];

  mat2(0,0) = v1[0]; mat2(0,1) = v2[0]; mat2(0,2) = v3[0];
  mat2(1,0) = v1[1]; mat2(1,1) = v2[1]; mat2(1,2) = v3[1];
  mat2(2,0) = v1[2]; mat2(2,1) = v2[2]; mat2(2,2) = v3[2];

  mat1.transpose();

  r1 = mat2*mat1*a1;
  r2 = mat2*mat1*a2;
  r3 = mat2*mat1*a3;

  fprintf (stderr, "\n dir cos \n");
  fprintf (stderr, ">>> r1=(%f %f %f) \n", r1[0], r1[1], r1[2]);
  fprintf (stderr, ">>> r2=(%f %f %f) \n", r2[0], r2[1], r2[2]);
  fprintf (stderr, ">>> r3=(%f %f %f) \n", r3[0], r3[1], r3[2]);

  // render scene //

  scene->render();

  // process events

  win->setScene (scene);
  win->processEvents (NULL, false, "test>", NULL);
  }


