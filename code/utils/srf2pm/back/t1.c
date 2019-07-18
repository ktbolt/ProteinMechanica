
#include <stdio.h>

main ()
  {


  int n;

  float a[4];

  FILE *fp;

  char line[100];

  fp = fopen ("tb", "r");
  fgets (line, 100, fp);
  fread (&n, sizeof(int), 1, fp);
  fread(a, sizeof(float), 4, fp);

  printf ( "\n\n %f %f %f %f \n", a[0], a[1], a[2], a[3]);
  }


