/*------------------------------------------------------------*
 *                                                            *
 *                      ****  mrc2pm  ****                    *
 *                                                            *
 * convert mrcm density map into a pm file.                   *
 *------------------------------------------------------------*/

#include <stdio.h>
#include "mrc.h"

/*------------------------------------------------------------*
 *                                                            *
 *                      ****  main  ****                      *
 *                                                            *
 *------------------------------------------------------------*/

int 
main (int num_args, char **args)
  {

  FILE *fp;

  int err;

  char *fn;

  MRCHeader hdr;

  MRC mrc;

  int size;

  float dx, dy, dz;

  char *data;

  float *fdata;

  static char *modes[] = { "byte", "short", "float", "short complex", "float complex",
                           "?", "unsigned short"};

  fn = args[1];

  loadMRC (fn, &mrc);

  fprintf (stderr, "\n>>> nc = %d  nr = %d  ns = %d \n", mrc.header.mapc, mrc.header.mapr,
           mrc.header.maps);

  fprintf (stderr, "\n>>> nx = %d  ny = %d  nz = %d \n", mrc.header.nx, mrc.header.ny,
           mrc.header.nz);

  fprintf (stderr, "\n>>> mx = %d  my = %d  mz = %d \n", mrc.header.mx, mrc.header.my,
           mrc.header.mz);

  fprintf (stderr, "\n>>> lx = %f  ly = %f  lz = %f \n", mrc.header.x_length,
           mrc.header.y_length, mrc.header.z_length);

  fprintf (stderr, "\n>>> amin = %f  amax = %f  amean = %f \n", mrc.header.amin, 
           mrc.header.amax, mrc.header.amean);

  fprintf (stderr, "\n>>> mode = %s  \n", modes[mrc.header.mode]);


  /* write pm grid file  */

  fprintf (stdout, "# pm grid file \n");
  fprintf (stdout, "# source file \"%s\"\n", fn);
  fprintf (stdout, "grid \"%s\" \n", fn);

  dx = mrc.header.x_length /  mrc.header.nx;
  dy = mrc.header.y_length /  mrc.header.ny;
  dz = mrc.header.z_length /  mrc.header.nz;
  dx /= 10.0;
  dy /= 10.0;
  dz /= 10.0;

  fprintf (stdout, "binary\n");
  fwrite (&mrc.header.nx, sizeof(int), 1, stdout);
  fwrite (&mrc.header.ny, sizeof(int), 1, stdout);
  fwrite (&mrc.header.nz, sizeof(int), 1, stdout);

  fwrite (&dx, sizeof(float), 1, stdout);
  fwrite (&dy, sizeof(float), 1, stdout);
  fwrite (&dz, sizeof(float), 1, stdout);

  data = mrc.pbyData;
  size = mrc.header.nx * mrc.header.ny * mrc.header.nz; 
  fwrite (data, sizeof(float), size, stdout);
  }


