/*------------------------------------------------------------*
 *                                                            *
 *                      ****  msms2pm  ****                   *
 *                                                            *
 * convert msms surface files into pm file.                   *
 *------------------------------------------------------------*/

#include <stdio.h>

typedef struct Tri {
  int n1, n2, n3;
  short c;
  } Tri;

typedef  float  (DmPoint3)[3];
typedef  int    (DmIntPoint3)[3];

int
main (int num_args, char **args)
  {
  FILE *fp;
  char c, *str, *name, s[10000], fname[80];
  int i, j, n;
  char blanks[80];

  int num_verts;
  float *verts;
  int vert_start;
 
  int num_tri, *conn, tri_start;
  Tri *tri;

  int nodes[3];
  short d;
  char date[1000];
  int num_sphere;
  float x, y, z; 
  float f1, f2, f3;
  int i1, i2, i3;
  int binary;

  binary = 1;
  name = args[1];
  sprintf (fname, "%s.vert", name); 
  fp = fopen (fname, "r");

  fgets (s, 10000, fp);
  fprintf (stderr, "\n >>>> %s ", s);
  fgets (s, 10000, fp);
  fprintf (stderr, "\n >>>> %s ", s);


  /*  read in verts  */

  fscanf (fp, "%d %d %f %f \n", &num_verts, &num_sphere, &f1, &f2);
  fprintf (stderr, "\n >>>>>> num verts [%d] \n", num_verts);
  verts = (float*)malloc(sizeof(float) * 3 * num_verts);

  for (i = 0; i < num_verts; i++) {
    fscanf (fp, "%f %f %f %f %f %f %d %d %d \n", &x, &y, &z, 
            &f1, &f2, &f3, &i1, &i2, &i3);
    verts[3*i] = x; 
    verts[3*i+1] = y; 
    verts[3*i+2] = z; 
    }

  fclose (fp);


  /*  read in faces  */

  sprintf (fname, "%s.face", name); 
  fp = fopen (fname, "r");
  fgets (s, 10000, fp);
  fprintf (stderr, "\n >>>> %s ", s);
  fgets (s, 10000, fp);
  fscanf (fp, "%d %d %f %f \n", &num_tri, &num_sphere, &f1, &f2);
  fprintf (stderr, "\n >>>>>> num faces [%d] \n", num_tri);
  conn = (int*)malloc(sizeof(int) * 3 * num_tri);

  for (i = 0; i < num_tri; i++) {
    fscanf (fp, "%d %d %d %d %d \n", &conn[3*i], &conn[3*i+1], &conn[3*i+2], &i1, &i2);
    conn[3*i] -= 1;  
    conn[3*i+1] -= 1;
    conn[3*i+2] -= 1; 
    }


  /* write pm surface file  */

  fprintf (stdout, "# pm surface file \n"); 
  fprintf (stdout, "# source file \"%s\"\n", fname); 
  fprintf (stdout, "surface \"%s\" \n", name);

  if (binary) {
    fprintf (stdout, "binary\n");
    fwrite (&num_verts, sizeof(int), 1, stdout);
    fwrite (&num_tri, sizeof(int), 1, stdout);
    fwrite (verts, sizeof(DmPoint3), num_verts, stdout);
    fwrite (conn, sizeof(DmIntPoint3), num_tri, stdout);
    return 0;
    }

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
  fprintf (stdout, "end surface \n");

  return 0;
  }

