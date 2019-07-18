
#include <stdio.h>
#include <math.h>

int
main(int nargs, char **args)
  {
  double eps, sig;
  double dx, x, energy, xb, xe, v, err, perr;
  double c12, c6; 
  double r0, k, r, xmax, xmin;
  int n, i;

  sig = 3.39967e-01; 
  eps = 3.59824e-01;
  eps = eps / 1.654;

  r0 = pow(2.0, 1.0/6.0) * sig;
  k = 36.0*eps /(pow(2.0, 2.0/3.0) * sig*sig);
  fprintf (stderr, "\n");
  fprintf (stderr, ">>> k=%g\n", k);

  n = 500;
  xb = 0.3;
  xe = 1.0;
  dx = (xe - xb) / n;
  x = xb;
  xmax = -1e6;
  xmin =  1e6;

  for (i = 0; i < n; i++) {
    c6 = pow(sig / x, 6.0);
    c12 = pow(sig / x, 12.0);
    energy = 4.0*eps*(c12 - c6);
    r = x - r0;
    v = -eps + 0.5*k*r*r;
    err = fabs(energy-v);
    perr = fabs(err) / fabs(energy);
    fprintf (stdout, "%g %g %g %g\n", x, energy, v, 100.0*perr);

    if (perr < 0.2) {
      fprintf (stderr, ">>> x=%g perr=%g energy=%g v=%g \n", x, 100.0*perr, energy, v);
      if (x < xmin) xmin = x;
      if (x > xmax) xmax = x;
      }

    x += dx;
    }

  fprintf (stderr, "\n");
  fprintf (stderr, ">>> r0=%g\n", r0);
  fprintf (stderr, ">>> xmin=%g  xmax=%g  dx=%g \n", xmin, xmax, xmax-xmin);

  /*
  c6 = 4.0*eps*pow(sig,6.0);
  c12 = 4.0*eps*pow(sig,12.0);
  fprintf (stderr, ">>> c6=%g  c12=%g \n", c6, c12); 
  */
  }

