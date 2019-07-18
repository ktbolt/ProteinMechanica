
#include "pm/pm.h"

using namespace ProteinModeler;


int
main (int nargs, char **args) 
  {

  // matrix mult //

  PmMatrix3x3 A, B, C;

  A(0,0) = 1; A(0,1) = 2; A(0,2) = 3;
  A(1,0) = 4; A(1,1) = 5; A(1,2) = 6;
  A(2,0) = 7; A(2,1) = 8; A(2,2) = 9;

  B(0,0) = 3;  B(0,1) = 5;  B(0,2) = 8;
  B(1,0) = 11; B(1,1) = 13; B(1,2) = 9;
  B(2,0) = 1;  B(2,1) = 2;  B(2,2) = 4;

  C = A * B;
  C = B * A * B;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf (stderr, " %g ", C(i,j)); 
      }

    fprintf (stderr, "\n");
    }

  fprintf (stderr, "---------------\n");

  C = 2.0*A;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf (stderr, " %g ", C(i,j));
      }

    fprintf (stderr, "\n");
    }

  }


