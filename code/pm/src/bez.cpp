
/*------------------------------------------------------------*
 *                                                            *
 *                    ****  bez  ****                         *
 *                                                            *
 * bezier curve module.                                       *
 *------------------------------------------------------------*/

#include "bez.h"
#include "bez_prv.h"

namespace ProteinMechanica {

/*------------------------------------------------------------*
 *                                                            *
 *                    ****  pm_BezierFit  ****                *
 *                                                            *
 * fit a bezier curve to a set of points.                     *
 *------------------------------------------------------------*/

void 
pm_BezierFit (int num_pts, PmVector3 *pts, float error, int *p_num_curves,
              PmVector3 **curves)
  {
  PmVector3 tang_1, tang_2;

  /*
  fprintf (stderr, "\n>>>>>>  bez:pm_BezierFit \n");
  fprintf (stderr, ">>> fit error=%f \n", error);
  fprintf (stderr, ">>> num pts=%d \n", num_pts);
  fprintf (stderr, ">>> num pts x 4=%d \n", num_pts*4);
  */

  *p_num_curves = 0;
  *curves = (PmVector3*)malloc(sizeof(PmVector3) * num_pts * 4);

  fit_left_tangent_comp (0, pts, tang_1);
  //fprintf (stderr, ">>> tang_1=%f %f %f\n", tang_1[0], tang_1[1], tang_1[2]);

  fit_right_tangent_comp (num_pts - 1, pts, tang_2);
  //fprintf (stderr, ">>> tang_2=%f %f %f\n", tang_2[0], tang_2[1], tang_2[2]);

  fit_cubic (pts, 0, num_pts - 1, tang_1, tang_2, error, p_num_curves, *curves);
  }

/*------------------------------------------------------------*
 *                                                            *
 *                  ****  pm_BezierEval  ****                 *
 *                                                            *
 * evaluate a bezier curve.                                   *
 *------------------------------------------------------------*/

void
pm_BezierEval (int num_control_pts, PmVector3 *control_pts, int num_pts, 
               PmVector3 *pts, PmVector3 *tang, PmVector3 *norm)
  {
  //fprintf (stderr, "\n>>>>>> pm_BezierEval: \n");
  //fprintf (stderr, ">>> num_control_pts=%d \n", num_control_pts);

  int i, j, degree;
  float coef[10], tcoef[10], ncoef[10];
  float u, du, x, y, z, u2, tx, ty, tz, nx, ny, nz; 

  degree = 3; 
  du = 1.0 / (float)(num_pts-1);
  u = 0.0;

  for (i = 0; i < num_pts; i++) {
    if (i == (num_pts - 1)) {
      u = 1.0;
      }

    pm_BezierBernsteinGet (degree, u, coef);

    u2 = u*u;
    x = y = z = 0.0; 
    tx = ty = tz = 0.0; 
    nx = ny = nz = 0.0; 
    tcoef[0] = -3.0 + 6.0*u - 3.0*u2;
    tcoef[1] = 3.0*(1.0 - 4.0*u + 3.0*u2); 
    tcoef[2] = 3.0*(2.0*u - 3.0*u2); 
    tcoef[3] = 3.0*u*u;

    ncoef[0] = 6.0 - 6.0*u; 
    ncoef[1] = 3.0*(-4.0 + 6.0*u);
    ncoef[2] = 3.0*(2.0 - 6.0*u);
    ncoef[3] = 6.0*u;

    /*
    fprintf (stderr, ">>> %d: u=%f coef=%f %f %f %f\n", i, u, coef[0], coef[1],
             coef[2], coef[3]);
    fprintf (stderr, "             tcoef=%f %f %f %f\n", tcoef[0], tcoef[1],
             tcoef[2], tcoef[3]);
    fprintf (stderr, "             ncoef=%f %f %f %f\n", ncoef[0], ncoef[1],
             ncoef[2],ncoef[3]);
    */

    for (j = 0; j <= degree; j++) {
      x += coef[j] * control_pts[j][0];
      y += coef[j] * control_pts[j][1];
      z += coef[j] * control_pts[j][2];

      tx += tcoef[j] * control_pts[j][0];
      ty += tcoef[j] * control_pts[j][1];
      tz += tcoef[j] * control_pts[j][2];

      nx += ncoef[j] * control_pts[j][0];
      ny += ncoef[j] * control_pts[j][1];
      nz += ncoef[j] * control_pts[j][2];
      }

     pts[i][0] = x;   pts[i][1] = y;   pts[i][2] = z;
    tang[i][0] = tx; tang[i][1] = ty; tang[i][2] = tz;
    norm[i][0] = nx; norm[i][1] = ny; norm[i][2] = nz;

    /*
    fprintf (stderr, ">>> %d: u=%f tang=%f %f %f\n", i, u, tang[i][0], tang[i][1], 
             tang[i][2]);
    */

    u += du;
    }
  }

/*------------------------------------------------------------*
 *                                                            *
 *           ****  pm_BezierBernsteinGet  ****                *
 *                                                            *
 * get the bernstein coeficients for a given parameter.       *
 *------------------------------------------------------------*/

void
pm_BezierBernsteinGet (int degree, float u, float *coef) 
  {
  float u1, last, tmp;
  int i, j;

  coef[0] = 1.0;
  u1 = 1.0 - u;

  for (i = 1; i <= degree; i++) {
    last = 0.0;

    for (j = 0; j < i; j++) {
      tmp = coef[j];
      coef[j] = last + u1 * tmp;
      last = u * tmp;
      }

    coef[i] = last;
    }
  }

/*------------------------------------------------------------*
 *                    internal functions                      *
 *------------------------------------------------------------*/

/*------------------------------------------------------------*
 *                                                            *
 *               ****  fit_left_tangent_comp  ****            *
 *                                                            *
 * approximate unit tangent at left endpoint.                 * 
 *------------------------------------------------------------*/

static void
fit_left_tangent_comp (int end, PmVector3 *pts, PmVector3& tang)	
  {
  tang = pts[end+1] - pts[end];
  tang.normalize();
  }

/*------------------------------------------------------------*
 *                                                            *
 *               ****  fit_right_tangent_comp  ****           *
 *                                                            *
 * approximate unit tangent at right endpoint.                * 
 *------------------------------------------------------------*/

static void
fit_right_tangent_comp (int end, PmVector3 *pts, PmVector3& tang)	
  {
  tang = pts[end-1] - pts[end];
  tang.normalize();
  }

/*------------------------------------------------------------*
 *                                                            *
 *                    ****  fit_cubic  ****                   *
 *                                                            *
 * fit a bezier curve to a (sub)set of points.                *
 *------------------------------------------------------------*/

static void 
fit_cubic (PmVector3 *pts, int first, int last, PmVector3 tang_1, 
           PmVector3 tang_2, float error, int *p_num_curves, PmVector3 *curves)
  {

  PmVector3 *bezcurve; 
  float	*u, *u_prime, max_error;
  int split_pt, num_pts;	
  float	iter_error; 
  int max_iter; 
  PmVector3 tang_center;   	
  PmVector3 stang_1, stang_2, dv;
  int i;		
  float dist; 

  //fprintf (stderr, "\n>>>>>>  bez:fit_cubic \n");
  //fprintf (stderr, ">>> fit error=%f \n", error);
  //fprintf (stderr, ">>> first=%d \n", first);
  //fprintf (stderr, ">>> last=%d \n", last);

  max_iter = 10; 
  iter_error = error * error;
  num_pts = last - first + 1;
  //fprintf (stderr, ">>> num_pts=%d \n", num_pts);

  if (num_pts == 2) {
    //fprintf (stderr, ">>> first=%d  last=%d\n", first, last);
    dv = pts[last] - pts[first]; 
    dist = dv.length();
    dist /= 3.0;
    //fprintf (stderr, "   >>> dv=%f %f %f\n", dv[0], dv[1], dv[2]);

    bezcurve = (PmVector3 *)malloc(4 * sizeof(PmVector3));
    bezcurve[0] = pts[first];
    bezcurve[3] = pts[last];

    stang_1 = dist * tang_1;
    bezcurve[1] = stang_1 + bezcurve[0];

    stang_2 = dist * tang_2;
    bezcurve[2] = stang_2 + bezcurve[3];

    fit_output_curve (3, bezcurve, p_num_curves, curves);

    return;
    }

  //===== parameterize points, and attempt to fit curve =====//

  fit_chord_length_parameterize (pts, first, last, &u);

  generate_bezier (pts, first, last, u, tang_1, tang_2, &bezcurve);


  //===== find max deviation of points to fitted curve =====//

  fit_max_error_comp (pts, first, last, bezcurve, u, &split_pt, &max_error);

  //fprintf (stderr, ">>> max error [%f] \n", max_error);

  if (max_error < error) {
    /*
    DrawBezierCurve(3, bezcurve);
    */
    fit_output_curve (3, bezcurve,  p_num_curves, curves);
    free((void *)u);
    free((void *)bezcurve);
    return;
    }

  /*  if error not too large, try some  
      reparameterization and iteration.  */ 

  if (max_error < iter_error) {
    //fprintf (stderr, ">>> re-parameterize \n");

    for (i = 0; i < max_iter; i++) {
      //fprintf (stderr, "\n  fit_cubic: iter [%d] \n", i+1);

      fit_reparameterize (pts, first, last, u, bezcurve, &u_prime);

      generate_bezier (pts, first, last, u_prime, tang_1, tang_2, &bezcurve);

      fit_max_error_comp (pts, first, last, bezcurve, u_prime, &split_pt, 
                          &max_error);

      if (max_error < error) {
        fit_output_curve (3, bezcurve,  p_num_curves, curves);
        free((void *)u);
        free((void *)bezcurve);
        return;
        }

      free((void *)u);
      u = u_prime;
      }
    }


  /*  fitting failed so split at max 
      error point and fit recursively.   */

  //fprintf (stderr, ">>> split at max error point \n");
  free((void *)u);
  free((void *)bezcurve);

  fit_center_tangent_comp (pts, split_pt, tang_center);

  fit_cubic (pts, first, split_pt, tang_1, tang_center, error, p_num_curves, 
             curves);

  tang_center = -tang_center;

  fit_cubic (pts, split_pt, last, tang_center, tang_2, error, p_num_curves,
             curves);
  }


/*------------------------------------------------------------*
 *                                                            *
 *            ****  fit_chord_length_parameterize  ****       *
 *                                                            *
 * assign parameter values to the input points using relative *
 * distances between points.                                  *
 *------------------------------------------------------------*/

static void
fit_chord_length_parameterize (PmVector3 *pts, int first, int last, float **p_u) 
  {

  int i;	

  float	*u;

  float	dist;

  PmVector3 dv;

  u = (float *)malloc((unsigned)(last-first+1) * sizeof(float));
  u[0] = 0.0;

  for (i = first+1; i <= last; i++) {
    dv = pts[i] - pts[i-1]; 
    dist = dv.length();
    u[i-first] = u[i-first-1] + dist; 
    }

  for (i = first + 1; i <= last; i++) {
    u[i-first] = u[i-first] / u[last-first];
    }

  *p_u = u;
  }

/*------------------------------------------------------------*
 *                                                            *
 *                    ****  generate_bezier  ****             *
 *                                                            *
 * use least-squares method to find bezier control points for *
 * region.                                                    *
 *------------------------------------------------------------*/

static void
generate_bezier (PmVector3 *pts, int first, int last, float *u_prime, 
                 PmVector3 tang_1, PmVector3 tang_2, PmVector3 **p_curve)
  {

  int i;
  PmVector3 A[MAXPOINTS][2];
  int num_pts;	
  float C[2][2], X[2], det_C0_C1, det_C0_X, det_X_C1, alpha_l, alpha_r;
  PmVector3 tmp, *bezcurve, v1, v2;
  float s, dist; 
  PmVector3 sv1, sv2, res, av1, av2;

  bezcurve = (PmVector3 *)malloc(4 * sizeof(PmVector3));
  num_pts = last - first + 1;
 
  /*  compute the A's.  */

  for (i = 0; i < num_pts; i++) {
    v1 = tang_1;
    v2 = tang_2;

    s = bezmult (1, u_prime[i]);
    v1 = s * v1;

    s = bezmult (2, u_prime[i]);
    v2 = s * v2;

    A[i][0] = v1;
    A[i][1] = v2;

    if (i >= MAXPOINTS) {
      fprintf (stderr, "\n  **** error [bez fit] memory overflow. \n");
      }
    }


  /*  create the C and X matrices.  */

  C[0][0] = 0.0;
  C[0][1] = 0.0;
  C[1][0] = 0.0;
  C[1][1] = 0.0;
  X[0]    = 0.0;
  X[1]    = 0.0;

  for (i = 0; i < num_pts; i++) {
    C[0][0] += A[i][0]*A[i][0];
    C[0][1] += A[i][0]*A[i][1];
    C[1][0] = C[0][1];
    C[1][1] += A[i][1]*A[i][1];

    s = bezmult (2, u_prime[i]);
    sv1 = s * pts[last];

    s = bezmult (3, u_prime[i]);
    sv2 = s* pts[last];
 
    av1 = sv1 + sv2;

    s = bezmult (1, u_prime[i]);
    res = s * pts[first];
    av2 = av1 + res;

    s = bezmult (0, u_prime[i]);
    sv1 = s * pts[first];
    res = sv1 + av2;
    tmp = pts[first + i] - res;

    X[0] += A[i][0] * tmp;
    X[1] += A[i][1] * tmp;
    }

  /*  compute the determinants of C and X.  */

  det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
  det_C0_X  = C[0][0] * X[1]    - C[0][1] * X[0];
  det_X_C1  = X[0]    * C[1][1] - X[1]    * C[0][1];

  /*  finally, derive alpha values.  */

  if (det_C0_C1 == 0.0) {
    det_C0_C1 = (C[0][0] * C[1][1]) * 10e-12;
    }

  alpha_l = det_X_C1 / det_C0_C1;
  alpha_r = det_C0_X / det_C0_C1;

  /*  if alpha negative, use the Wu/Barsky heuristic.  */

  if (alpha_l < 0.0 || alpha_r < 0.0) {
    PmVector3 dv;
    dv = pts[last] - pts[first];
    dist = dv.length();
    dist /= 3.0;

    bezcurve[0] = pts[first];
    bezcurve[3] = pts[last];

    sv1 = dist * tang_1; 
    bezcurve[1] = bezcurve[0]  + sv1;

    sv2 = dist * tang_2; 
    bezcurve[2] = bezcurve[3] + sv2;

    *p_curve = bezcurve;
    return;
    }


  /*  first and last control points of the Bezier curve are 
      positioned exactly at the first and last data points.
      control points 1 and 2 are positioned an alpha distance 
      out on the tangent vectors, left and right, respectively.   */

  bezcurve[0] = pts[first];
  bezcurve[3] = pts[last];

  sv1 = alpha_l * tang_1; 
  bezcurve[1] = bezcurve[0] + sv1;

  sv2 = alpha_r * tang_2; 
  bezcurve[2] = bezcurve[3] + sv2;

  *p_curve = bezcurve;
  }

/*------------------------------------------------------------*
 *                                                            *
 *                    ****  bezmult  ****                     *
 *                                                            *
 *------------------------------------------------------------*/

static float 
bezmult (int n, float u)
  {
  float tmp, coef; 

  switch (n) {
    case 0:
      tmp = 1.0 - u;
      coef = tmp * tmp * tmp;;
    break;

    case 1:
      tmp = 1.0 - u;
      coef = 3.0 * u * (tmp * tmp);
    break;

    case 2:
      tmp = 1.0 - u;
      coef = 3.0 * u * u * tmp;
    break;

    case 3:
      tmp = 1.0 - u;
      coef = u * u * u;
    break;
    }

  return (coef);
  }

/*------------------------------------------------------------*
 *                                                            *
 *                    ****  fit_max_error_comp  ****          *
 *                                                            *
 * find the maximum squared distance of digitized points to   *
 * fitted curve.                                              *
 *------------------------------------------------------------*/

static void
fit_max_error_comp (PmVector3 *pts, int first, int last, PmVector3 *bezcurve,
                    float *u, int *split_pt, float *error)
  {

  int i;
  float	max_dist;
  float	dist;
  PmVector3 P;
  PmVector3 v;

  //fprintf (stderr, "\n>>>>>> fit_max_error_comp: \n");

  *split_pt = (last - first + 1) / 2;
  max_dist = 0.0;

  for (i = first + 1; i < last; i++) {
    fit_bezier_eval (3, bezcurve, u[i-first], P);
    v = P - pts[i];
    dist = v * v; 
    //fprintf (stderr, ">>> dist=%f \n", dist);

    if (dist >= max_dist) {
      max_dist = dist;
      *split_pt = i;
      }
    }

  *error = sqrt(max_dist);
  }

/*------------------------------------------------------------*
 *                                                            *
 *                    ****  fit_bezier_eval  ****             *
 *                                                            *
 * evaluate a bezier curve at a particular parameter value.   *
 *------------------------------------------------------------*/

static void
fit_bezier_eval (int degree, PmVector3 *V, float t, PmVector3& pt)
  {

  int i, j;		
  PmVector3 *Vtemp;

  /*  copy array.  */

  Vtemp = (PmVector3*)malloc((unsigned)((degree+1) * sizeof (PmVector3)));

  for (i = 0; i <= degree; i++) {
    Vtemp[i] = V[i];
    }

  /*  triangle computation.  */

  for (i = 1; i <= degree; i++) {	
    for (j = 0; j <= degree-i; j++) {
      Vtemp[j][0] = (1.0 - t) * Vtemp[j][0] + t * Vtemp[j+1][0];
      Vtemp[j][1] = (1.0 - t) * Vtemp[j][1] + t * Vtemp[j+1][1];
      }
    }

  pt[0] = Vtemp[0][0];
  pt[1] = Vtemp[0][1];
  pt[2] = Vtemp[0][2];
  free((void *)Vtemp);
  }

/*------------------------------------------------------------*
 *                                                            *
 *               ****  fit_reparameterize  ****               *
 *                                                            *
 * given set of points and their parameterization, try to     *
 * find a better parameterization.                            *
 *------------------------------------------------------------*/

static void
fit_reparameterize (PmVector3 *pts, int first, int last, float *u, 
                    PmVector3 *bezcurve, float **new_u)
  {

  int num_pts;

  int i;

  float	*u_prime;

  num_pts = last - first + 1;	
  u_prime = (float *)malloc(num_pts * sizeof(float));

  for (i = first; i <= last; i++) {
    u_prime[i-first] = fit_root_find (bezcurve, pts[i], u[i-first]);
    }

  *new_u = u_prime;
  }


/*------------------------------------------------------------*
 *                                                            *
 *                    ****  fit_root_find  ****               *
 *                                                            *
 * use newton-raphson iteration to find better root.          *
 *------------------------------------------------------------*/

static float 
fit_root_find (PmVector3 *Q, PmVector3 P, float u)
  {

  float numerator, denominator;
  PmVector3 Q1[4], Q2[3];	
  PmVector3 Q_u, Q1_u, Q2_u; 
  float u_prime;
  int i;

  fit_bezier_eval (3, Q, u, Q_u);
    

  /*  generate control vertices for Q'.  */

  for (i = 0; i <= 2; i++) {
    Q1[i][0] = (Q[i+1][0] - Q[i][0]) * 3.0;
    Q1[i][1] = (Q[i+1][1] - Q[i][1]) * 3.0;
    }
    

  /*  generate control vertices for Q''.  */

  for (i = 0; i <= 1; i++) {
    Q2[i][0] = (Q1[i+1][0] - Q1[i][0]) * 2.0;
    Q2[i][1] = (Q1[i+1][1] - Q1[i][1]) * 2.0;
    }
    

  /*  compute Q'(u) and Q''(u).  */

  fit_bezier_eval (2, Q1, u, Q1_u);
  fit_bezier_eval (1, Q2, u, Q2_u);
    

  /*  compute f(u)/f'(u).  */

  numerator = (Q_u[0] - P[0]) * (Q1_u[0]) + (Q_u[1] - P[1]) * (Q1_u[1]);
  denominator = (Q1_u[0])*(Q1_u[0]) + (Q1_u[1])*(Q1_u[1]) + 
                (Q_u[0] - P[0])*(Q2_u[0]) + (Q_u[1] - P[1])*(Q2_u[1]);
    

  /* u = u - f(u)/f'(u) */

  u_prime = u - (numerator / denominator);

  return (u_prime);
  }


/*------------------------------------------------------------*
 *                                                            *
 *               ****  fit_center_tangent_comp  ****          *
 *                                                            *
 *------------------------------------------------------------*/

static void
fit_center_tangent_comp (PmVector3 *pts, int center, PmVector3& tang) 
  {

  PmVector3 V1, V2;

  V1 = pts[center-1] - pts[center];
  V2 = pts[center] - pts[center+1];

  tang[0] = (V1[0] + V2[0]) / 2.0;
  tang[1] = (V1[1] + V2[1]) / 2.0;
  tang[2] = (V1[2] + V2[2]) / 2.0;
  // Note: this is orig statement //
  tang[2] = 0.0; 
 
  tang.normalize();
  }

/*------------------------------------------------------------*
 *                                                            *
 *                    ****  fit_output_curve  ****            *
 *                                                            *
 *------------------------------------------------------------*/

static void
fit_output_curve (int n, PmVector3 *bezcurve, int *p_num_curves, PmVector3 *curves)
  {

  int i, m;

  m = *p_num_curves * 4; 

  for (i = 0; i < n+1; i++) {
    curves[m][0] = bezcurve[i][0]; 
    curves[m][1] = bezcurve[i][1]; 
    curves[m][2] = bezcurve[i][2]; 
    m += 1;
    }

  *p_num_curves += 1; 
  }

};
