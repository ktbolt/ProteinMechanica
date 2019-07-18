/*------------------------------------------------------------*
 *                                                            *
 *                      ****  sfc2pm  ****                    *
 *                                                            *
 * convert Swiss PDBViewer sfc file to pm file.               *
 *------------------------------------------------------------*/

#include <stdio.h>
#include "sfc2pm.h" 

/*------------------------------------------------------------*
 *                                                            *
 *                      ****  main  ****                      *
 *                                                            *
 *------------------------------------------------------------*/

int
main (int num_args, char **args)
  {
  FILE *fp;
  char c, *str, *name, s[10000], fname[80];
  int i, j, n;
  char blanks[80];
  int num_verts, start_verts;
  float *verts;
  int vert_start;
  int num_tri, *conn, start_tri;
  int start_norm, start_charge;
  int nodes[3];
  int fsize;
  float *charge;
  short d;
  char date[1000];
  int num_sphere;
  float x, y, z; 
  float f1, f2, f3;
  int i1, i2, i3;
  int binary;
  short s1, s2, s3;

  binary = 1;
  name = args[1];
  sprintf (fname, "%s.sfc", name); 
  //sprintf (fname, "%s.", name); 
  fp = fopen (fname, "r");

  fgets (s, 10000, fp);
  fprintf (stderr, "\n >>>> %s ", s);


  /*  read in verts  */

  fscanf (fp, "NbOfVertices %d\n", &num_verts);
  fprintf (stderr, "\n >>>>>> num verts [%d] \n", num_verts);

  fscanf (fp, "NbOfTriangles %d\n", &num_tri);
  fprintf (stderr, "\n >>>>>> num tri [%d] \n", num_tri);

  fscanf (fp, "StartVertices %d\n", &start_verts);
  fprintf (stderr, "\n >>>>>> start verts [%d] \n", start_verts);

  fscanf (fp, "StartNormals%d\n", &start_norm);
  fprintf (stderr, "\n >>>>>> start norm [%d] \n", start_norm);

  fscanf (fp, "StartCharges%d\n", &start_charge);
  fprintf (stderr, "\n >>>>>> start charge [%d] \n", start_charge);

  fscanf (fp, "StartTriangle%d\n", &start_tri);
  fprintf (stderr, "\n >>>>>> start tri [%d] \n", start_tri);

  fsize = 3*num_verts * sizeof(float) +  
          3*num_tri   * sizeof(int)   +
            num_tri   * sizeof(short) +
          3*num_verts * sizeof(float) +  
          num_verts   * sizeof(float);  

  verts = (float*)malloc(sizeof(float) * 3 * num_verts);
  charge = (float*)malloc(sizeof(float) * num_verts);
  conn = (int*)malloc(sizeof(int) * 3 * num_tri);

  rewind (fp);
  fseek(fp, start_verts, SEEK_SET);
  fread (verts, sizeof(float), 3*num_verts, fp);

  rewind (fp);
  fseek(fp, start_tri, SEEK_SET);

  for (i = 0; i < num_tri; i++) {
    fread (&i1, sizeof(int), 1, fp); 
    fread (&i2, sizeof(int), 1, fp); 
    fread (&i3, sizeof(int), 1, fp); 
    fread (&s1, sizeof(short), 1, fp); 

    conn[3*i]   = i1;
    conn[3*i+1] = i2;
    conn[3*i+2] = i3;

    if (s1 != 0) {
      //fprintf (stderr, "s1 = %d \n", s1);
      }
    }

  rewind (fp);
  fseek(fp, start_charge, SEEK_SET);
  fread (charge, sizeof(float), num_verts, fp);
  fclose (fp);

  /*
  fprintf (stderr, "\n---------- charge ---------- \n");
  for (i = 0; i < 30; i++) {
    fprintf (stderr, "%f \n", charge[i]);
    }
  */

  /*
  fprintf (stderr, "\n>>> reorient mesh ... ");
  PolyReorient (num_tri, conn, num_verts, verts);
  fprintf (stderr, "done \n");
  */

  /* write pm surface file  */

  fprintf (stderr, "\n>>> write pm file ... \n");
  fprintf (stdout, "# pm surface file \n");
  fprintf (stdout, "# source file \"%s\"\n", fname);
  fprintf (stdout, "surface \"%s\" \n", name);
  fprintf (stdout, "charge\n");

  if (binary) {
    fprintf (stdout, "binary\n");
    fwrite (&num_verts, sizeof(int), 1, stdout);
    fwrite (&num_tri, sizeof(int), 1, stdout);
    fwrite (verts, sizeof(float), 3*num_verts, stdout);
    fwrite (conn, sizeof(int), 3*num_tri, stdout);
    fwrite (charge, sizeof(float), num_verts, stdout);
    return;
    }
  }

