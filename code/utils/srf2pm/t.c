
#include <stdio.h>

main ()
  {


  int n;

  float a[4];

  FILE *fp;

  char line[100];

  a[0] = 0.0;
  a[1] = 1.0;
  a[2] = 2.0;
  a[3] = 3.0;
  n = 4;

  fprintf (stdout, " bin file \n");
  fwrite (&n, sizeof(int), 1, stdout);
  fwrite (a, sizeof(float), 4, stdout);

  }


