
/*------------------------------------------------------------*
 *                                                            *
 *                      ****  srf2pm  ****                    *
 *                                                            *
 * convert SwissViewer srf file to pm file.                   *
 *------------------------------------------------------------*/

#include "dm.h" 
#include "gr.h" 
#include "gc.h" 
#include <stdio.h>
#include <math.h>

typedef struct Tri {
  int n1, n2, n3;
  short c;
  } Tri;


/*------------------------------------------------------------*
 *                                                            *
 *                      ****  read  ****                      *
 *                                                            *
 *------------------------------------------------------------*/

void
read (char *name, DmObj **p_geom)
  {

  FILE *fp;

  char c, *str, s[10000], fname[80];

  int i, j, n;

  char blanks[80];

  int num_verts;

  float *verts;

  int vert_start;
 
  int num_tri, *conn, tri_start, *nconn;

  Tri *tri;

  int nodes[3];

  short d;

  char date[1000];

  GrGeometryType geom_type;

  GrGeometryOptions geom_opts;

  DmObj *geom;

  float color[3];

  float x, y, z;

  float dist, dx, dy, dz;

  float tol, x1, y1, z1;

  float area; 

  DmPoint3 pts[3];

  int n1, n2, n3;

 /**************
  ***  body  ***
  **************/

  sprintf (fname, "%s", name); 
  fp = fopen (fname, "r");

  fgets (s, 10000, fp);
  fprintf (stderr, "\n s [%s] \n", s);

  fgets (s, 10000, fp);
  fprintf (stderr, "\n s [%s] \n", s);
  sscanf (s, "NbOfVertices %d \n", &num_verts);
  fprintf (stderr, " num verts [%d] \n", num_verts);
  verts = (float*)malloc(sizeof(float) * 3 * num_verts);

  fgets (s, 10000, fp);
  fprintf (stderr, "\n s [%s] \n", s);
  sscanf (s, "NbOfTriangles %d \n", &num_tri);
  fprintf (stderr, " num tri [%d] \n", num_tri);
  conn = (int*)malloc(sizeof(int) * 3 * num_tri);
  tri = (Tri*)malloc(sizeof(Tri) * num_tri);

  fgets (s, 10000, fp);
  fprintf (stderr, " str [%s] \n", s);
  sscanf (s, "StartVertices %d \n", &vert_start);
  fprintf (stderr, " vert start [%d] \n", vert_start);

  fgets (s, 10000, fp);
  fprintf (stderr, " str [%s] \n", s);

  fgets (s, 10000, fp);
  fprintf (stderr, " str [%s] \n", s);

  fgets (s, 10000, fp);
  fprintf (stderr, " str [%s] \n", s);
  sscanf (s, "StartTriangle %d \n", &tri_start);
  fprintf (stderr, " tri start [%d] \n", tri_start);


  fseek (fp, vert_start, SEEK_SET);
  fread (verts, sizeof(float), num_verts*3, fp);
  /*
  for (i = 0; i < num_verts; i++) {
    fprintf (stderr, "%d  (%f %f %f) \n", i, verts[3*i], verts[3*i+1], verts[3*i+2]);
    }
  */

  fseek (fp, tri_start, SEEK_SET);
  /*
  fread (conn, sizeof(int), num_tri*3, fp);
  fprintf (stderr, " sizeof(Tri)  [%d] \n", sizeof(Tri)); 
  fprintf (stderr, " sizeof(short)  [%d] \n", sizeof(short)); 
  fprintf (stderr, " sizeof(int*3)  [%d] \n", sizeof(int)*3); 
  fread (tri, sizeof(Tri), num_tri, fp);
  */

  for (i = 0; i < num_tri; i++) {
    fread (conn+3*i, sizeof(int), 3, fp);
    fread (&d, sizeof(short), 1, fp);
    }


  /*  check verts  */

/*
  tol = 1.0e-06;

  for (i = 0; i < num_verts-1; i++) {
    x1 = verts[3*i]; 
    y1 = verts[3*i+1]; 
    z1 = verts[3*i+2];

    for (j = i+1; j < num_verts; j++) {
      x = verts[3*j]; 
      y = verts[3*j+1]; 
      z = verts[3*j+2];
      dx = x - x1;
      dy = y - y1;
      dz = z - z1;
      dist = sqrt (dx*dx + dy*dy + dz*dz);

      if (dist < tol) {
        fprintf (stderr, " dupe verts  %d  %d \n", i, j); 
        exit (0);
        }

      }
    }
*/


  /*  check area */

  nconn = (int*)malloc(sizeof(int) * 3 * num_tri);
  tol = 1.0e-06;
  fprintf (stderr, "\n--------- area check ---------- \n"); 
  n = 0;

  for (i = 0, j = 0; i < num_tri; i++) {
    n1 = conn[3*i]; 
    n2 = conn[3*i+1]; 
    n3 = conn[3*i+2];
    pts[0][0] = verts[3*n1]; pts[0][1] = verts[3*n1+1]; pts[0][2] = verts[3*n1+2];
    pts[1][0] = verts[3*n2]; pts[1][1] = verts[3*n2+1]; pts[1][2] = verts[3*n2+2];
    pts[2][0] = verts[3*n3]; pts[2][1] = verts[3*n3+1]; pts[2][2] = verts[3*n3+2];

    gc_PtsAreaComp (3, pts, &area);

    if (area < tol) {
      fprintf (stderr, " poly [%d]   area [%f] \n", i, area); 
      }
    else {
      nconn[3*n]   = n1; 
      nconn[3*n+1] = n2; 
      nconn[3*n+2] = n3;
      n += 1;
      }
    }

  /*
  num_tri = n + 1;
  fprintf (stderr, "\n >>>>>> new num tri [%d] \n", num_tri);
  */

   nconn = conn;


  /*  create geom  */

  geom_type = GR_GEOMETRY_INDEX_POLYGONS;
  geom_opts = GR_GEOMETRY_OPTION_TRIANGLES;

  gr_GeomCreate ("surf", geom_type, geom_opts, &geom);

  gr_GeomIndexPolygonsSet (geom, num_verts, verts, num_tri, nconn);

  color[0] = 0.0, color[1] = 1.0, color[2] = 0.0;

  gr_GeomColorSet (geom, color);

  gr_GeomShadingSet (geom, GR_GEOMETRY_SHADING_FLAT);

  /*
  gr_GeomDisplaySet (geom, GR_GEOMETRY_DISPLAY_LINE);
  */

  gr_GeomPickEnableSet (geom, DM_TRUE);

  *p_geom = geom;



  }


/*------------------------------------------------------------*
 *                                                            *
 *                      ****  write  ****                     *
 *                                                            *
 *------------------------------------------------------------*/

void
write (char *fname, DmObj *geom)
  {

 /**************
  ***  body  ***
  **************/

#ifdef use 

  fprintf (stdout, "# pm surface file \n"); 
  fprintf (stdout, "# source file \"%s\"\n", fname); 
  dm_TimeStrGet (date);
  fprintf (stdout, "# created %s \n", date); 
  fprintf (stdout, "number of vertices %d \n", num_verts);
  fprintf (stdout, "vertices \n");

  for (i = 0; i < num_verts; i++) {
    fprintf (stdout, "%d  %f %f %f \n", i, verts[3*i], verts[3*i+1], verts[3*i+2]);
    }

  fprintf (stdout, "end vertices \n\n");
  fprintf (stdout, "number of triangles %d \n", num_tri);
  fprintf (stdout, "connectivity \n");

  for (i = 0; i < num_tri/*num_tri*/; i++) {
    fprintf (stdout, "%d  %d %d %d \n", i, conn[3*i], conn[3*i+1], conn[3*i+2]);
    }

  fprintf (stdout, "end connectivity \n");

#endif


  }

