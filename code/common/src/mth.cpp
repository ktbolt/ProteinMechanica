/* -------------------------------------------------------------------------- *
 *                              Protein Mechanica                             *
 * -------------------------------------------------------------------------- *
 * This is part of the Protein Mechanica coarse-grained molecular motor       *
 * modeling application originating from Simbios, the NIH National Center for *
 * Physics-Based Simulation of Biological Structures at Stanford, funded      *
 * under the NIH Roadmap for Medical Research, grant U54 GM072970.            *
 * See https://simtk.org.                                                     *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: David Parker                                                      *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

//*============================================================*
//* mth:                      m a t h                          *
//*============================================================*
// numerical and math functions.
//
// The Ecf functions implement ellipsoid contact using the elliptic contact function
// (ECF) approach (Perram J. Comp. Phys. 58, 409 1985).
//
// The calculation of a closed polygon mesh mass properties is performed by computing
// surface integrals over the mesh and equating that with a volume inegral using
// Gauss's Theorem. (A. M. Messner and G. Q. Taylor, ACM TOMS, Vol. 6, No. 1, 1980).


#include "pm/mth.h"
#include <limits.h>
#include <algorithm>

namespace ProteinMechanica {

#define PM "PM math" 

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

// rotation macro for eigenvalue computation // 

#define MTH_ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
        a[k][l]=h+s*(g-h*tau);

#define MTH_MAX_ROTATIONS  50


// edge list for computing closed polygonal mesh //

typedef struct PmMathEdgeList {
  int count;
  int id;
  struct PmMathEdgeList *next;
  } PmMathEdgeList;


// constants for random number generation //

static
float   snorm_a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
    };

static
float   snorm_d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
    };

static
float   snorm_t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
    };

static
float   snorm_h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
    };

//*============================================================*
//*==========             pm_MathSign                ==========*
//*============================================================*
// return the sign of a number.

double
pm_MathSign(double a) 
 {
 if (a >= 0.0) {
   return 1.0 ;
   }
 else {
   return -1.0;
   }   
 }

//*============================================================*
//*==========            pm_MathBasisComp            ==========*
//*============================================================*
// compute an orthonormal basis from the given normal vector.

void
pm_MathBasisCompute (PmVector3& normal, PmVector3& u, PmVector3& v)
  {
  int i, j, n;
  PmVector3 q1, q2, q3;

  u.set(0,0,0);
  v.set(0,0,0);
  n = 0;

  for (i = 0; i < 3; i++) {
    if (normal[i] != 0.0) {
      j = (i + 1) % 3; 
      u[j] = 1.0;
      n += 1;
      if (n == 2) break;
      }
    }

  v = normal.cross(u);
  pm_MathGramSchmidt (normal, u, v, q1, q2, q3);
  u = q2;
  v = q3;
  }

///////////////////////////////////////////////////////////////
//             e l l i p s o i d   c o n t a c t            //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========            pm_MathEcfCompContact       ==========*
//*============================================================*
// compute contact for two ellipsoids.  

void
pm_MathEcfCompContact (PmMathContactEllipsoids& ellipsoids, bool& contact)
  {
  float *r1, *r2, dist, s1, s2; 
  float max_r1, max_r2, d;
  PmVector3 R; 

  contact = false;
  R = ellipsoids.pos2 - ellipsoids.pos1;

  // compute quadratic forms (A and B) for ellipsoids //

  pm_MathEcfQuadForm (ellipsoids.r1, ellipsoids.u1, ellipsoids.v1, ellipsoids.w1, 
                      ellipsoids.A);
  pm_MathEcfQuadForm (ellipsoids.r2, ellipsoids.u2, ellipsoids.v2, ellipsoids.w2, 
                      ellipsoids.B);

  // determine if ellipsoids are in contact //

  pm_MathEcfContactFuncSolve (ellipsoids.A, ellipsoids.B, R, &ellipsoids.lambda, 
                              &ellipsoids.cfunc);

  if (ellipsoids.cfunc > 1.0) {
    return;
    }

  // compute contact point //

  pm_MathEcfContactPtComp (ellipsoids.A, ellipsoids.B, ellipsoids.pos1, ellipsoids.pos2, 
                           ellipsoids.lambda, ellipsoids.cfunc, ellipsoids.npt1, 
                           ellipsoids.npt2, ellipsoids.cpt);
  contact = true;
  }

//*============================================================*
//*==========            pm_MathEcfDerivEval         ==========*
//*============================================================*
// evaluate ECF derivatives.            

float
pm_MathEcfDerivEval (PmMatrix3x3& A, PmMatrix3x3& B, PmVector3& R, float val)
  {

  PmMatrix3x3 G, C, D;
  PmVector3 X, V;
  bool singular;
  float v1, v2, result;

  pm_MathEcfInterp (A, B, 1.0-val, val, G);

  pm_MathMatrix3x3Solve (G, R, X, &singular);

  if (singular) {
    fprintf (stderr, "**** PM Error: singular matrix in ECF deriv computation.");
    return (0.0);
    }

  v1 = (1.0-val)*(1.0-val);
  v2 = -val*val;
  pm_MathEcfInterp (A, B, v1, v2, C);

  V = C * X;
  result = V * X;
  return (result);
  }

//*============================================================*
//*==========            pm_MathEcfEval              ==========*
//*============================================================*
// evaluate the ECF function.

float
pm_MathEcfEval (PmMatrix3x3& A, PmMatrix3x3& B, PmVector3& R, float val)
  {

  PmMatrix3x3 G, C;
  PmVector3 V;
  bool singular;
  float result;

  pm_MathEcfInterp (A, B, 1.0-val, val, G);

  pm_MathMatrix3x3Inverse (G, C, &singular);

  if (singular) {
    fprintf (stderr, "**** PM Error: singular matrix in ECF eval computation.");
    return (0.0);
    }

  V = C * R;
  result = val*(1.0-val) * (V*R);
  return (result);
  }

//*============================================================*
//*==========            pm_MathEcfInterp            ==========*
//*============================================================*
// interpolate matrices for ECF derivatives.

void 
pm_MathEcfInterp (PmMatrix3x3& A, PmMatrix3x3& B, float v1, float v2, PmMatrix3x3& G)
  {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      G(i,j) = v1*A(i,j) + v2*B(i,j);
      }
    }
  }

//*============================================================*
//*==========        pm_MathEcfContactFuncSolve      ==========*
//*============================================================*
// solve for the contact function.

void
pm_MathEcfContactFuncSolve (PmMatrix3x3& A, PmMatrix3x3& B, PmVector3& R, 
                            float *p_lambda, float *p_cfunc)
  {
  PmMatrix3x3 G;
  float a, b, fa, fb, tmp, c, fc, d;
  int mflag;
  static double tol = 1e-5;
  static int num_it = 100;

  a = 0.1;
  b = 0.9;

  fa = pm_MathEcfDerivEval (A, B, R, a);
  fb = pm_MathEcfDerivEval (A, B, R, b);

  for (int i = 0; i < num_it; i++) {
    if ((fb == 0.0) || (fabs(b-a) < tol)) {
      break;
      }

    c = (a + b) / 2.0;
    fc = pm_MathEcfDerivEval (A, B, R, c);

    if (fa*fc > 0) {
      a = c;
      fa = pm_MathEcfDerivEval (A, B, R, a);
      }
    else {
      b = c;
      }
    }

  *p_lambda = b;
  *p_cfunc = pm_MathEcfEval (A, B, R, b);
  }

//*============================================================*
//*==========        pm_MathEcfContactPtComp         ==========*
//*============================================================*
// compute the contact points for two ellipsoids.

void
pm_MathEcfContactPtComp (PmMatrix3x3& A, PmMatrix3x3& B, PmVector3& r, PmVector3& s, 
                         float lambda, float F, PmVector3& xa, PmVector3& xb, 
                         PmVector3& cpt)
  {
  int i;
  PmMatrix3x3 G, Ginv;
  PmVector3 r1, s1, r2, s2;
  PmVector3 R, ar, bs;
  float f;
  bool singular;

  // compute contact point //

  pm_MathEcfInterp (A, B, 1.0-lambda, lambda, G);
  pm_MathMatrix3x3Inverse (G, Ginv, &singular);

  if (singular) {
    fprintf (stderr, "**** PM Error: singular matrix in ECF pt computation.");
    xa.set(0,0,0);
    xb.set(0,0,0);
    cpt.set(0,0,0);
    return;
    }

  if (1) {
    R = s - r;
    r1 = Ginv * R;
    r2 = A * r1;
    cpt = (1.0-lambda)*r2 + r;
    }
  else {
    ar = A*r;
    r1 = Ginv * ar;
    bs = B*s;
    s1 = Ginv * bs;
    cpt = (1.0-lambda)*r1 + lambda*s1;
    }

  // compute points on ellipsoids //

  f = 1.0 / sqrt(F);
  xa  = r + f*(cpt - r);
  xb  = s + f*(cpt - s);
  }

//*============================================================*
//*==========            pm_MathEcfQuadForm          ==========*
//*============================================================*
// comute quadratic form A for an ellipsoid.

void
pm_MathEcfQuadForm (float *r, PmVector3& u, PmVector3& v, PmVector3& w, PmMatrix3x3& A)
  {
  float rs[3]; 

  for (int i = 0; i < 3; i++) {
    rs[i] = r[i]*r[i];

    for (int j = 0; j < 3; j++) {
      A(i,j) = 0.0;
      }
    }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      A(i,j) += rs[0]*u[i]*u[j] + rs[1]*v[i]*v[j] + rs[2]*w[i]*w[j];
      }
    }
  }

///////////////////////////////////////////////////////////////
//                         m a t r i x                      //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========        pm_MathEigMatrix3x3             ==========*
//*============================================================*
// compute eigenvalues and vectors for a 3x3 matrix.
// eigen vectors returned in columns of <v>.

void
pm_MathEigMatrix3x3 (float a[3][3], float w[3], float v[3][3], bool& singular)
  {
  int i, j, k, iq, ip, numPos;
  float tresh, theta, tau, t, sm, s, h, g, c;
  float b[3], z[3], tmp;
  static int n = 3;
  singular = false;

  for (ip = 0; ip < n; ip++) {
    for (iq = 0; iq < n; iq++) {
      v[ip][iq] = 0.0;
      }

    v[ip][ip] = 1.0;
    }

  for (ip = 0; ip < n; ip++) {
    b[ip] = w[ip] = a[ip][ip];
    z[ip] = 0.0;
    }

  // begin rotation sequence // 

  for (i = 0; i< MTH_MAX_ROTATIONS; i++) {
    sm = 0.0;

    for (ip = 0; ip < n-1; ip++) {
      for (iq = ip+1; iq < n; iq++) {
	sm += fabs(a[ip][iq]);
	}
      }

    if (sm == 0.0) {
      break;
      }

    if (i < 3) {
      tresh = 0.2*sm / (n*n);
      }
    else {
      tresh = 0.0;
      }

    for (ip = 0; ip < n-1; ip++) {
      for (iq = ip+1; iq < n; iq++) {
        g = 100.0*fabs(a[ip][iq]);

        if (i > 3 && (fabs(w[ip])+g) == fabs(w[ip]) && 
           (fabs(w[iq])+g) == fabs(w[iq])) {
          a[ip][iq] = 0.0;
          }
        else if (fabs(a[ip][iq]) > tresh) {
          h = w[iq] - w[ip];

          if ( (fabs(h)+g) == fabs(h)) {
	    t = (a[ip][iq]) / h;
	    }
          else {
            theta = 0.5*h / (a[ip][iq]);
            t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0)
	      {
	      t = -t;
	      }
            }

          c = 1.0 / sqrt(1+t*t);
          s = t*c;
          tau = s/(1.0+c);
          h = t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          w[ip] -= h;
          w[iq] += h;
          a[ip][iq]=0.0;

          // ip already shifted left by 1 unit

          for (j = 0;j <= ip-1;j++) {
            MTH_ROTATE(a,j,ip,j,iq)
            }

          // ip and iq already shifted left by 1 unit

          for (j = ip+1;j <= iq-1;j++) {
            MTH_ROTATE(a,ip,j,j,iq)
            }

          // iq already shifted left by 1 unit

          for (j = iq+1; j < n; j++) {
            MTH_ROTATE(a,ip,j,iq,j)
            }

          for (j = 0; j < n; j++) {
            MTH_ROTATE(v,j,ip,j,iq)
            }
          }
        }
      }

    for (ip = 0; ip < n; ip++) {
      b[ip] += z[ip];
      w[ip] = b[ip];
      z[ip] = 0.0;
      }
    }

  if ( i >= MTH_MAX_ROTATIONS ) {
    return;
    }

  //  sort eigenfunctions // 

  for (j = 0; j < n-1; j++) {
    k = j;
    tmp = w[k];

    for (i = j+1; i < n; i++) {
      if (w[i] >= tmp)    {
        k = i;
        tmp = w[k];
        }
      }

    if (k != j) {
      w[k] = w[j];
      w[j] = tmp;

      for (i = 0; i < n; i++) {
        tmp = v[i][j];
        v[i][j] = v[i][k];
        v[i][k] = tmp;
        }
      }
    }

  // Jacobi can compute vectors that are negative of       //
  // one another so select the most positive eigenvector.  //

  for (j = 0; j < n; j++) {
    for (numPos = 0, i = 0; i < n; i++) {
      if (v[i][j] >= 0.0) {
	numPos++;
	}
      }

    if (numPos < ceil(n/2.0)) {
      for(i=0; i<n; i++) {
	v[i][j] *= -1.0;
	}
      }
    }
  }

//*============================================================*
//*==========        pm_MathGramSchmidt              ==========*
//*============================================================*
// orthonormalize the vectors <v1>, <v2> and <v3>. 

void
pm_MathGramSchmidt (PmVector3& v1, PmVector3& v2, PmVector3& v3, PmVector3& q1, 
                    PmVector3& q2, PmVector3& q3)
  {
  int i;
  float len1, len2, dp1, dp2;

  dp1 = v1*v2;
  len1 = v1*v1;

  for (i = 0; i < 3; i++) {
    q1[i] = v1[i];
    q2[i] = v2[i] - dp1 / len1 * v1[i];
    }

  dp1 = q1*v3;
  dp2 = q2*v3;
  len2 = q2*q2;

  for (i = 0; i < 3; i++) {
    q3[i] = v3[i] - (dp1/len1)*q1[i] - (dp2/len2)*q2[i];
    }

  q1.normalize();
  q2.normalize();
  q3.normalize();
  }

//*============================================================*
//*==========        pm_MathMatrix3x3Inverse         ==========*
//*============================================================*
// compute inverse of a 3x3 matrix. 

void
pm_MathMatrix3x3Inverse(PmMatrix3x3& src, PmMatrix3x3& res, bool *singular)
  {

  float denom;

  float a, b, c, d, e, f, g, h, i;

  *singular = false;
  a = src(0,0), b = src(1,0), c = src(2,0);
  d = src(0,1), e = src(1,1), f = src(2,1);
  g = src(0,2), h = src(1,2), i = src(2,2);
  denom = c * (d*h - e*g) + b * (f*g - d*i) + a * (e*i - f*h);

  if (denom != 0.0) {
    denom = 1.0 / denom;
    }
  else {
    *singular = true;
    return;
    }

  res(0,0) = (e*i - f*h) * denom;
  res(1,0) = (c*h - b*i) * denom;
  res(2,0) = (b*f - c*e) * denom;

  res(0,1) = (f*g - d*i) * denom;
  res(1,1) = (a*i - c*g) * denom;
  res(2,1) = (c*d - a*f) * denom;

  res(0,2) = (d*h - e*g) * denom;
  res(1,2) = (b*g - a*h) * denom;
  res(2,2) = (a*e - b*d) * denom;
  }

//*============================================================*
//*==========        pm_MathMatrix3x3Solve           ==========*
//*============================================================*
// solve a 3x3 system of equations.                      

void
pm_MathMatrix3x3Solve (PmMatrix3x3& src, PmVector3& rhs, PmVector3& res, bool *singular)
  {
  PmMatrix3x3 inv;

  pm_MathMatrix3x3Inverse(src, inv, singular);

  if (*singular) {
    return;
    }

  res = inv * rhs;
  }

//*============================================================*
//*==========        pm_MathMatrix3x3DirCos          ==========*
//*============================================================*
// compute direction cosines for the coord frames given by
// (u1, u2, u3) to that given by (v1, v2, v3).

void
pm_MathMatrix3x3DirCos (PmVector3& u1, PmVector3& u2, PmVector3& u3, 
                        PmVector3& v1, PmVector3& v2, PmVector3& v3, 
                        PmMatrix3x3& mat)
  {
  mat(0,0) = u1*v1; mat(0,1) = u1*v2; mat(0,2) = u1*v3;
  mat(1,0) = u2*v1; mat(1,1) = u2*v2; mat(1,2) = u2*v3;
  mat(2,0) = u3*v1; mat(2,1) = u3*v2; mat(2,2) = u3*v3;
  }

//*============================================================*
//*==========        pm_MathFrameRotation            ==========*
//*============================================================*
// compute rotation matrix transforming axes (u1, u2, u3) into 
// axes (v1, v2, v3).    

void
pm_MathFrameRotation(PmVector3& u1, PmVector3& u2, PmVector3& u3,
                     PmVector3& v1, PmVector3& v2, PmVector3& v3, PmMatrix3x3& rmat)
  {
  PmMatrix3x3 mat1, mat2;

  mat1(0,0) = u1[0]; mat1(0,1) = u2[0]; mat1(0,2) = u3[0];
  mat1(1,0) = u1[1]; mat1(1,1) = u2[1]; mat1(1,2) = u3[1];
  mat1(2,0) = u1[2]; mat1(2,1) = u2[2]; mat1(2,2) = u3[2];

  mat2(0,0) = v1[0]; mat2(0,1) = v2[0]; mat2(0,2) = v3[0];
  mat2(1,0) = v1[1]; mat2(1,1) = v2[1]; mat2(1,2) = v3[1];
  mat2(2,0) = v1[2]; mat2(2,1) = v2[2]; mat2(2,2) = v3[2];

  mat1.transpose();
  rmat = mat2*mat1;
  }

//*============================================================*
//*==========        pm_MathMatrixToQuaternion       ==========*
//*============================================================*
// extract a quaternion from a 3x3 rotation matrix.

void
pm_MathMatrixToQuaternion(PmMatrix3x3& m, PmQuaternion q)
  {

  float tr, s; 

  int i;

  tr = m(0,0) + m(1,1) + m(2,2);

  if (tr >= 0) {
    s = sqrt(tr + 1);
    q[0] = 0.5 * s;
    s = 0.5 / s;
    q[1] = (m(2,1) - m(1,2)) * s;
    q[2] = (m(0,2) - m(2,0)) * s;
    q[3] = (m(1,0) - m(0,1)) * s;
    }
  else {
    i = 0;

    if (m(1,1) > m(0,0)) {
      i = 1;
      }

    if (m(2,2) > m(i,i)) {
      i = 2;
      }

    switch (i) {
      case 0:
        s = sqrt((m(0,0) - (m(1,1) + m(2,2))) + 1);
        q[1] = 0.5 * s;
        s = 0.5 / s;
        q[2] = (m(0,1) + m(1,0)) * s;
        q[3] = (m(2,0) + m(0,2)) * s;
        q[0] = (m(2,1) - m(1,2)) * s;
      break;

      case 1:
        s = sqrt((m(1,1) - (m(2,2) + m(0,0))) + 1);
        q[2] = 0.5 * s;
        s = 0.5 / s;

        q[3] = (m(1,2) + m(2,1)) * s;
        q[1] = (m(0,1) + m(1,0)) * s;
        q[0] = (m(0,2) - m(2,0)) * s;
      break;

      case 2:
        s = sqrt((m(2,2) - (m(0,0) + m(1,1))) + 1);
        q[3] = 0.5 * s;
        s = 0.5 / s;
        q[1] = (m(2,0) + m(0,2)) * s;
        q[2] = (m(1,2) + m(2,1)) * s;
        q[0] = (m(1,0) - m(0,1)) * s;
      break;
      }
    }
  }

//*============================================================*
//*==========        pm_MathQuaternionToMatrix       ==========*
//*============================================================*
// convert a quaternion into a 3x3 rotation matrix.

void
pm_MathQuaternionToMatrix(PmQuaternion q, PmMatrix3x3& m)
  {

  float s, x, y, z;

  s  = q[0];
  x = q[1];
  y = q[2];
  z = q[3];

  m(0,0) = 1.0 - 2*(y*y + z*z);
  m(0,1) = 2.0*(x*y - s*z);
  m(0,2) = 2.0*(x*z + s*y);

  m(1,0) = 2.0*(x*y + s*z);
  m(1,1) = 1.0 - 2.0*(x*x + z*z);
  m(1,2) = 2.0*(y*z - s*x);

  m(2,0) = 2.0*(x*z - s*y);
  m(2,1) = 2.0*(y*z + s*x);
  m(2,2) = 1.0 - 2.0*(x*x + y*y);
  }

///////////////////////////////////////////////////////////////
//             r a n d o m   n u m b e r s                  //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========            pm_MathRandNumber           ==========*
//*============================================================*
// compute a uniform random number in range (0, 1).          

float 
pm_MathRandNumber (int long *seed)
  {

  float rnum;
  int k;
  float am;
  static int long a = 16807;
  static int long m = 2147483647;
  static int long q = 127773;
  static int long mask = 123459876;
  static int long r = 2836;

  am = 1.0 / m;
  *seed ^= mask;

  k = (*seed) / q;
  *seed = a * (*seed - k*q) - r*k;

  if (*seed < 0) {
    *seed += m;
    }

  rnum = am * (float)(*seed);
  *seed ^= mask;
  return rnum;
  }

//*============================================================*
//*==========            pm_MathRandSnormNumber      ==========*
//*============================================================*
// generate a random number from a normal distribution. 

float
pm_MathRandSnormNumber(int long *seed)
  {
  int long i;
  float snorm, u, s, ustar, aa, w, y, tt;

  u = pm_MathRandNumber(seed);
  s = 0.0;

  if (u > 0.5) s = 1.0;

  u += (u-s);
  u = 32.0*u;
  i = (int long)(u);

  if (i == 32) i = 31;
  if (i == 0) goto S100;

  // start center  //

  ustar = u - (float)i;
  aa = *(snorm_a+i-1);

S40:
  if (ustar <= *(snorm_t+i-1)) goto S60;
  w = (ustar-*(snorm_t+i-1))**(snorm_h+i-1);

S50:

  // exit // 

  y = aa+w;
  snorm = y;

  if (s == 1.0) snorm = -y;
    return (snorm);

S60:

  // center //

  u =  pm_MathRandNumber(seed); 
  w = u*(*(snorm_a+i)-aa);
  tt = (0.5*w+aa)*w;
  goto S80;


S70:
  tt = u;
  ustar =  pm_MathRandNumber(seed); 

S80:
  if(ustar > tt) goto S50;
  u =  pm_MathRandNumber(seed);     
  if(ustar >= u) goto S70;
  ustar =  pm_MathRandNumber(seed); 
  goto S40;

S100:

  /* START TAIL */

  i = 6;
  aa = *(snorm_a+31);
  goto S120;

S110:
  aa += *(snorm_d+i-1);
  i += 1;

S120:
  u += u;
  if(u < 1.0) goto S110;
  u -= 1.0;

S140:
  w = u**(snorm_d+i-1);
  tt = (0.5*w+aa)*w;
  goto S160;

S150:
  tt = u;

S160:
  ustar = pm_MathRandNumber(seed); 
  if(ustar > tt) goto S50;
  u =  pm_MathRandNumber(seed); 
  if(ustar >= u) goto S150;
  u =  pm_MathRandNumber(seed); 
  goto S140;
  }

//*============================================================*
//*==========            pm_MathRandGaussNumber      ==========*
//*============================================================*
// generate a random number with a gaussian distribution with
// the given mean and standard deviation.

float 
pm_MathRandGaussNumber(const float mean, const float sd, int long *seed)
  {
  return sd * pm_MathRandSnormNumber(seed) + mean;
  }

//*============================================================*
//*==========            pm_MathRandGaussVector      ==========*
//*============================================================*
// generate a random vector.                          

void
pm_MathRandGaussVector(const float mean, const float sd, PmVector3& vec, int long *seed)
  {
  vec[0] = pm_MathRandGaussNumber(mean, sd, seed);
  vec[1] = pm_MathRandGaussNumber(mean, sd, seed);
  vec[2] = pm_MathRandGaussNumber(mean, sd, seed);
  }

//*============================================================*
//*==========        pm_MathExtractRotations         ==========*
//*============================================================*
// extract the equivalent rotations from a 3x3 rotation matrix.

#define use_new_pm_MathExtractRotations
#ifdef use_new_pm_MathExtractRotations

void 
pm_MathExtractRotations (float mat[3][3], float rot[3])
  {

  float quot, angle, th1, th2, th3, costh2, temp[10];
  float r11, r13, r22, r23, r32;

  r11 = mat[0][0];
  r22 = mat[1][1];
  r13 = mat[0][2];
  r32 = mat[2][1];
  r23 = mat[1][2];

  if (fabs(fabs(r13)-1.0) <= 1e-10) {
    if (r13 > 0.0) {
      temp[0] = M_PI / 2.0; 
      } 
    else {
      temp[0] = -M_PI / 2.0; 
      }

    th2 = temp[0];

    if (r32 > 1.0) {
      temp[0] = 1.0;
      } 
    else {
      if (r32 < -1.0) {
        temp[1] = -1.0;
        } 
      else {
        temp[1] = r32; 
        }

      temp[0] = temp[1];
      }

    angle = asin(temp[0]);

    if (r22 >= 0.0) {
      temp[0] = angle;
      } 
    else {
      temp[0] = (M_PI - angle);
      }

    th1 = temp[0];
    th3 = 0.0;
    } 
  else {
    if (r13 > 1.0) {
      temp[0] = 1.;
      } 
    else {
      if (r13 < -1.0) {
        temp[1] = -1.0;
        } 
      else {
        temp[1] = r13; 
        }

      temp[0] = temp[1];
      }

    th2 = asin(temp[0]);
    costh2 = cos(th2);
    quot = -r23 / costh2;

    if ((quot > 1.)  ) {
      temp[0] = 1.;
      } 
    else {
      if ((quot < -1.)  ) {
        temp[1] = -1.;
        } 
      else {
        temp[1] = quot;
        }

      temp[0] = temp[1];
      }

    angle = asin(temp[0]);

    if ((mat[2][2] >= 0.)  ) {
      temp[0] = angle;
      } 
    else {
      temp[0] = (3.14159265358979-angle);
      }

    th1 = temp[0];
    quot = ((-mat[0][1])/costh2);

    if ((quot > 1.)  ) {
      temp[0] = 1.;
      } 
    else {
      if ((quot < -1.)  ) {
        temp[1] = -1.;
        } 
      else {
        temp[1] = quot;
        }

      temp[0] = temp[1];
      }

    angle = asin(temp[0]);

    if (r11 >= 0.0) {
      temp[0] = angle;
      } 
    else {
      temp[0] = (M_PI - angle);
      }

    th3 = temp[0];
    }

  if (th1 > M_PI) {
    temp[0] = (th1 - 2.0*M_PI);
    } 
  else {
    temp[0] = th1;
    }

  rot[0] = temp[0];
  rot[1] = th2;

  if (th3 > M_PI) {
    temp[0] = (th3 - 2.0*M_PI);
    } 
  else {
    temp[0] = th3;
    }

  rot[2] = temp[0];
  }

#else

//*============================================================*
//*==========        pm_MathExtractRotations         ==========*
//*============================================================*
// extract the equivalent rotations from a 3x3 rotation matrix.

void 
pm_MathExtractRotations (float mat[3][3], float rot[3])
  {
  float xr, yr, zr, arg, tmp, cy;
  float r11, r12, r13, r22, r23, r33;
  static float tol = 1.0e-10;

  r11 = mat[0][0];
  r12 = mat[0][1];
  r13 = mat[0][2];
  r22 = mat[1][1];
  r23 = mat[1][2];
  r33 = mat[2][2];

  xr = yr = zr = 0.0;
  yr = asin (-r13);
  cy = cos (yr);

  if (fabs(cy) < tol) {
    zr = 0.0;
    xr = acos (r22); 
    tmp = sin(xr) * sin(yr);

    if (pm_MathSign(tmp) != pm_MathSign(r33)) {
      xr = -xr;
      }
    }
  else {
    arg = (r23 / cy);
    if (arg > 1.0) arg = 1.0;
    xr = asin (arg);
    tmp = cos (xr) * cy;

    if (pm_MathSign(tmp) != pm_MathSign(r33)) {
      xr += M_PI;
      tmp = sin(xr) * cy;

      if (pm_MathSign(tmp) != pm_MathSign(r23)) {
        xr = -xr;
        }
      }

    arg = (r12 / cy);
    if (arg > 1.0) arg = 1.0;
    zr = asin (arg);
    tmp = cos(zr) * cy;

    if (pm_MathSign(tmp) != pm_MathSign(r11)) {
      zr += M_PI;
      tmp = sin(zr) * cy;

      if (pm_MathSign(tmp) != pm_MathSign(r12)) {
        zr = -zr;
        }
      }
    }

  rot[0] = xr;
  rot[1] = yr;
  rot[2] = zr;
  }
#endif

///////////////////////////////////////////////////////////////
//                      f i t t i n g                       //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========             pm_MathFitRms              ==========*
//*============================================================*
// fit two sets of coordinates. returns the transformation 
// that will translate and rotate v2[] to v1[].                                
// Implements method from Kabsch (Acta Cryst. (1978) A34, 827-828). 

void
pm_MathFitRms (vector<PmVector3>& v1, vector<PmVector3>& v2, vector<float> wt, 
               PmXform& xform, float& p_err)
  {
  double m[3][3], aa[3][3], x[3], xx[3];
  double sumwt, tol, sig, gam;
  double sg, bb, cc, err, etmp, tmp;
  int maxiter, iters, ix, iy, iz;
  double t1[3], t2[3];

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      m[i][j] = 0.0;
      aa[i][j] = 0.0;
      }

    m[i][i] = 1.0;
    t1[i] = 0.0;
    t2[i] = 0.0;
    }

  sumwt = 0.0;
  tol = 0.00001;
  maxiter = 1000;

  // calculate center // 

  int n = v1.size();
  fprintf (stderr, "---------------------------------------------------------------------------- \n");

  if (!wt.empty()) {
    fprintf (stderr, " use weights \n");
    for (int k = 0; k < n; k++) {
      for (int i = 0; i< 3; i++) {
        t1[i] += wt[k]*v1[k][i];
        t2[i] += wt[k]*v2[k][i];
        }

      if (wt[k] != 0.0) {
        sumwt = sumwt + wt[k];
        } 
      else {
        sumwt = sumwt + 1.0;
        }
      }
    } 

  else {
    for (int k = 0; k < n; k++) {
      for (int i = 0; i < 3; i++) {
        t1[i] += v1[k][i];
        t2[i] += v2[k][i];
        }

      sumwt += 1.0;
      }
    }

  if (sumwt == 0.0) {
    sumwt = 1.0;
    }

  for (int i = 0; i < 3; i++) {
    t1[i] /= sumwt;
    t2[i] /= sumwt;
    }

  fprintf (stderr, "sumwt %g \n", sumwt); 
  fprintf (stderr, "center t1 %g %g %g \n", t1[0], t1[1], t1[2]);
  fprintf (stderr, "center t2 %g %g %g \n", t2[0], t2[1], t2[2]);

  // create correlation matrix // 

  for (int k = 0; k < n; k++) {
    if (!wt.empty()) {
      for (int i = 0; i < 3; i++) {
        x[i]  = wt[k]*v1[k][i] - t1[i];
        xx[i] = wt[k]*v2[k][i] - t2[i];
        }
      } 
    else {
      for (int i = 0; i < 3; i++) {
        x[i]  = v1[k][i] - t1[i];
        xx[i] = v2[k][i] - t2[i];
        }
      }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        aa[i][j] = aa[i][j] + xx[i]*x[j];
        }
      }
    }

  fprintf (stderr, "\n------ aa ------ \n");
  for (int i = 0; i < 3; i++) {
    fprintf (stderr, "%g %g %g \n", aa[i][0], aa[i][1], aa[i][2]);
    }

  // perform svd // 

  if (n > 1) {
    iters = 0;

    while (1) {
      iz = (iters+1) % 3;
      iy = (iz+1) % 3;
      ix = (iy+1) % 3;
      sig = aa[iz][iy] - aa[iy][iz];
      gam = aa[iy][iy] + aa[iz][iz];

      if (iters >= maxiter) {
        fprintf (stderr, "\n no convergence after %d iterations. \n", iters);
        break;
        }

      tmp = sig*sig + gam*gam;
      sg = sqrt(tmp);

      if ((sg!=0.0F) &&(fabs(sig)>(tol*fabs(gam)))) {
        sg = 1.0 / sg;

        for (int i = 0; i < 3; i++) {
          bb = gam*aa[iy][i] + sig*aa[iz][i];
          cc = gam*aa[iz][i] - sig*aa[iy][i];
          aa[iy][i] = bb*sg;
          aa[iz][i] = cc*sg;
            
          bb = gam*m[iy][i] + sig*m[iz][i];
          cc = gam*m[iz][i] - sig*m[iy][i];
          m[iy][i] = bb*sg;
          m[iz][i] = cc*sg;
          }
        } 
      else {
        break;
        }

      iters++;
      }
    }

  fprintf (stderr, "---------------------- M ---------------------- \n");
  for (int i = 0; i < 3; i++) {
    fprintf (stderr, "%g %g %g \n", m[i][0], m[i][1], m[i][2]);
    }


  // calculate the weighted rms error //

  err = 0.0;

  for (int i = 0; i < 3; i++) {
    double mag = sqrt(m[i][0]*m[i][0] + m[i][1]*m[i][1] + m[i][2]*m[i][2]);
    m[i][0] /= mag;
    m[i][1] /= mag;
    m[i][2] /= mag;
    }

  fprintf (stderr, "---------------------- M weighted ---------------------- \n");
  for (int i = 0; i < 3; i++) {
    fprintf (stderr, "%g %g %g \n", m[i][0], m[i][1], m[i][2]);
    }

  /*
  fprintf (stderr, "\n------ m ------ \n");
  for (int i = 0; i < 3; i++) {
    fprintf (stderr, "%g %g %g \n", m[i][0], m[i][1], m[i][2]);
    }
  */

  for (int j = 0; j < n; j++) {
    etmp = 0.0;

    for (int i = 0; i < 3; i++) {
      tmp = m[i][0]*(v2[j][0]-t2[0]) + m[i][1]*(v2[j][1]-t2[1]) + 
            m[i][2]*(v2[j][2]-t2[2]);
      tmp = (v1[j][i]-t1[i])-tmp;
      etmp += tmp*tmp;
      }

    if (!wt.empty()) {
      err += wt[j] * etmp;
      }
    else {
      err += etmp;
      }
    }

  err = err / sumwt;
  err = sqrt(err);

  for (int i = 0; i < 3; i++) {
    xform.translation[i] = t1[i] - t2[i];
    xform.center[i] = t2[i];

    for (int j = 0; j < 3; j++) {
      xform.matrix(i,j) = m[i][j];
      }
    }

  fprintf (stderr, "\n-------------------- xform coords ------------------ \n");
   
   float xc[3];

   for (int j = 0; j < n; j++) {
      for (int i = 0; i< 3; i++) {
          tmp = m[i][0]*(v2[j][0]-t2[0]) + m[i][1]*(v2[j][1]-t2[1]) + 
                m[i][2]*(v2[j][2]-t2[2]);
          xc[i] = tmp +  xform.translation[i];
      } 

      fprintf (stderr, ">>> coord %g %g %g \n", v2[j][0], v2[j][1], v2[j][2]);
      fprintf (stderr, "    xc    %g %g %g \n", xc[0], xc[1], xc[2]);
   }


  p_err = err; 
  }

///////////////////////////////////////////////////////////////
//         p r i n c i p a l   c o m p o n e n t s          //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========        compPrincipalComponents         ==========*
//*============================================================*
// compute principal components.

void
pm_MathPrincipalComponents(vector<PmVector3>& coords, PmPcaResults& pca)
  {
  int n = coords.size();

  if (!n) {
    return;
    }

  PmVector3 *ca = new PmVector3[n];

  for (int i = 0; i < coords.size(); i++) {
    ca[i] = coords[i];
    }

  pm_MathPrincipalComponents(n, ca, pca);
  delete[] ca;
  }

//*============================================================*
//*==========        compPrincipalComponents         ==========*
//*============================================================*
// compute principal components.

void
pm_MathPrincipalComponents(int num_coords, PmVector3 *coords, PmPcaResults& pca)
  {
  float x, y, z, cx, cy, cz;
  float xx, xy, xz, yy, yz, zz;
  PmVector3 com, csum;
  float var[3][3], eignvecs[3][3], eignvals[3]; 
  bool status;

  cx = cy = cz = 0.0;

  for (int i = 0; i < num_coords; i++) {
    cx += coords[i][0];
    cy += coords[i][1];
    cz += coords[i][2];
    }

  com[0] = cx / num_coords;
  com[1] = cy / num_coords;
  com[2] = cz / num_coords;
  csum[0] = csum[1] = csum[2] = 0.0;

  for (int i = 0; i < num_coords; i++) {
    csum[0] += coords[i][0] - com[0];
    csum[1] += coords[i][1] - com[1];
    csum[2] += coords[i][2] - com[2];
    }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      var[i][j] = 0.0;
      }
    }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < num_coords; k++) {
        var[i][j] += (coords[k][i] - com[i]) * (coords[k][j] - com[j]);
        }

      var[i][j] = (var[i][j] - csum[i]*csum[j] / num_coords) / (num_coords - 1);
      }
    }
   
  pm_MathEigMatrix3x3 (var, eignvals, eignvecs, status);

  if (status) {
    fprintf (stderr, "\n**** PM Error: PCA eigenvalue computation failed. \n");
    return;
    }

  pca.com = com;
  pca.s1 = eignvals[0];
  pca.s2 = eignvals[1];
  pca.s3 = eignvals[2];
  
  for (int i = 0; i < 3; i++) {
    pca.axis1[i] = eignvecs[i][0];
    pca.axis2[i] = eignvecs[i][1];
    pca.axis3[i] = eignvecs[i][2];
    }
  }

//*============================================================*
//*==========     pm_MathPrincipalComponentsProj     ==========*
//*============================================================*
// project coordinats onto principal components.

void
pm_MathPrincipalComponentsProj(vector<PmVector3>& coords, const PmPcaResults& pca, 
                               PmExtent& proj, float w[3])
  {
  int n = coords.size();

  if (!n) {
    return;
    }

  PmVector3 *ca = new PmVector3[n];

  for (int i = 0; i < coords.size(); i++) {
    ca[i] = coords[i];
    }

  pm_MathPrincipalComponentsProj(n, ca, pca, proj, w);
  delete[] ca;
  }

//*============================================================*
//*==========     pm_MathPrincipalComponentsProj     ==========*
//*============================================================*
// project coordinats onto principal components.

void
pm_MathPrincipalComponentsProj(const int num_coords, const PmVector3 *coords, 
                               const PmPcaResults& pca, PmExtent& proj, float widths[3])
  {
  PmVector3 u, v, w, center, a;
  float p1, p2, p3;
  float p1_min, p1_max, p2_min, p2_max, p3_min, p3_max;
  vector<float>wd(3);

  u = pca.axis1;
  v = pca.axis2;
  w = pca.axis3;
  center = pca.com;
  
  for (int i = 0; i < num_coords; i++) {
    a = coords[i] - center;
    p1 = a * u;
    p2 = a * v;
    p3 = a * w;

    if (i == 0) {
      p1_max = p1_min = p1;
      p2_max = p2_min = p2;
      p3_max = p3_min = p3;
      }
    else {
      p1_max = (p1 > p1_max ? p1 : p1_max); 
      p2_max = (p2 > p2_max ? p2 : p2_max); 
      p3_max = (p3 > p3_max ? p3 : p3_max); 
      p1_min = (p1 < p1_min ? p1 : p1_min); 
      p2_min = (p2 < p2_min ? p2 : p2_min); 
      p3_min = (p3 < p3_min ? p3 : p3_min); 
      }
    }

  proj.min[0] = p1_min;
  proj.min[1] = p2_min;
  proj.min[2] = p3_min;

  proj.max[0] = p1_max;
  proj.max[1] = p2_max;
  proj.max[2] = p3_max;

  // compute maximum widths for each axis //

  wd[0] = (p1_max - p1_min) / 2.0;
  wd[1] = (p2_max - p2_min) / 2.0;
  wd[2] = (p3_max - p3_min) / 2.0;

  sort (wd.begin(), wd.end());

  for (int i = 0; i < 3; i++) {
    widths[2-i] = wd[i]; 
    }
  }

///////////////////////////////////////////////////////////////
//         m a t r i x    r o t a t i o n s                 //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========       pm_MathRotationAroundAxis        ==========*
//*============================================================*
// construct a 3x3 rotation matrix for a rotation about an axis.

void
pm_MathRotationAroundAxis (float ang, PmVector3 axis, PmMatrix3x3& mat)
  {

  float cosa, sina, mag, d; 
  PmVector3 u;
  PmMatrix3x3 s;

  ang = ang * M_PI / 180.0;
  cosa = cos (ang);
  sina = sin (ang);
  mag = axis.length(); 

  if (mag != 0.0) {
    u = (1.0 / mag)*axis;
    }
  else {
    u[0] = u[1] = u[2] = 0.0;
    }

  s(0,0) =     0;  s(0,1) = -u[2];  s(0,2) =  u[1];
  s(1,0) =  u[2];  s(1,1) =     0;  s(1,2) = -u[0];
  s(2,0) = -u[1];  s(2,1) =  u[0];  s(2,2) =     0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      d = (i == j);
      mat(i,j) = u[i]*u[j]  +  cosa*(d - u[i]*u[j])  -  sina*s(i,j);
      }
    }
  }

//*============================================================*
//*==========       pm_MathRotationFromAngles        ==========*
//*============================================================*
// construct a 3x3 rotation matrix from three angles.

void
pm_MathRotationFromAngles (float xr, float yr, float zr, PmMatrix3x3& mat)
  {

  float sx, sy, sz, cx, cy, cz, f;

  f = (4.0 * atan(1.0)) / 180.0;
  sx = sin(f*xr);
  sy = sin(f*yr);
  sz = sin(f*zr);

  cx = cos(f*xr);
  cy = cos(f*yr);
  cz = cos(f*zr);

  mat(0,0) = cy * cz;
  mat(0,1) = -cy * sz;
  mat(0,2) = sy;

  mat(1,0) = sx * sy * cz + cx * sz;
  mat(1,1) = -sx * sy * sz + cx * cz;
  mat(1,2) = -sx * cy;

  mat(2,0) = -cx * sy * cz + sx * sz;
  mat(2,1) = cx * sy * sz + sx * cz;
  mat(2,2) = cx * cy;
  }


////////////////////////////////////////////////////////////////
//               p o l y g o n   m e s h                     //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========             pm_MathPolyClosed          ==========*
//*============================================================*
// determine if a polygon mesh is closed.

void
pm_MathPolyClosed(int num_polys, int *conn, int hgeom, int num_verts, bool& closed)
  {

  if (num_polys == 1) {
    closed = true;
    return;
    }

  PmVector3 pts[10];
  PmMathEdgeList *ptr, *lptr;
  int ids[10], mnid, mxid, id1, id2;
  int i, j, n;

  PmMathEdgeList **etab = new PmMathEdgeList*[num_verts];

  for (i = 0; i < num_verts; i++) {
    etab[i] = NULL;
    }

  // if homogeneous polygons //

  if (hgeom) {
    n = hgeom; 
    }

  for (i = 0; i < num_polys; i++) {
    if (!hgeom) {
      n = *conn;
      conn++;
      }

    for (j = 0; j < n; j++) {
      ids[j] = *conn;
      conn++;
      }

    for (j = 0; j < n; j++) {
      id1 = ids[j];

      if (j == n-1) {
        id2 = ids[0];
        }
      else {
        id2 = ids[j+1];
        }

      mnid = (id1 < id2 ? id1 : id2);
      mxid = (id1 > id2 ? id1 : id2);
      ptr = etab[mnid]; 

      while (ptr) { 
        if (ptr->id == mxid) {
          ptr->count += 1; 
          /*
          if (mnid == 0) {
            fprintf (stderr, ">>>>>> dup edge (%d %d) \n", mnid, ptr->id);
            }
          */
          break;
          }

        ptr = ptr->next; 
        }

      if (!ptr) {
        ptr = new PmMathEdgeList;
        ptr->id = mxid; 
        ptr->count = 1; 
        ptr->next = etab[mnid]; 
        etab[mnid] = ptr; 
        }
      }
    }

  closed = true;

  for (i = 0; i < num_verts; i++) {
    ptr = etab[i];

    while (ptr) { 
      if (ptr->count != 2) {
        closed = false;
        //fprintf (stderr, ">>>>>> unmatched edge (%d %d) \n", i, ptr->id);
        }

      lptr = ptr;
      ptr = ptr->next; 
      delete lptr;
      }
    }

  delete[] etab;
  }

//*============================================================*
//*==========             gc_PolyMassProps           ==========*
//*============================================================*
// get mass properties of a closed polygon mesh.

void
pm_MathPolyMassProps (int num_poly, int *conn, int hgeom, int num_verts, 
                      PmVector3 *verts, PmVector3 *face_norms, float density, 
                      float& mass, PmVector3& com, float J[3][3])
  {

  PmMathVolIntData vd;

  // compute volume integral //

  pm_MathPolyVolIntComp(num_poly, conn, hgeom, num_verts, verts, face_norms, vd);

  mass = density * vd.T0;

  // center of mass //

  com[0] = vd.T1[0] / vd.T0;
  com[1] = vd.T1[1] / vd.T0;
  com[2] = vd.T1[2] / vd.T0;

  // inertia tensor //

  J[0][0] = density * (vd.T2[1] + vd.T2[2]);
  J[1][1] = density * (vd.T2[2] + vd.T2[0]);
  J[2][2] = density * (vd.T2[0] + vd.T2[1]);
  J[0][1] = J[1][0] = - density * vd.TP[0];
  J[1][2] = J[2][1] = - density * vd.TP[1];
  J[2][0] = J[0][2] = - density * vd.TP[2];

  // translate inertia tensor to center of mass //

  J[0][0] -= mass * (com[1]*com[1] + com[2]*com[2]);
  J[1][1] -= mass * (com[2]*com[2] + com[0]*com[0]);
  J[2][2] -= mass * (com[0]*com[0] + com[1]*com[1]);
  J[0][1] = J[1][0] += mass * com[0] * com[1];
  J[1][2] = J[2][1] += mass * com[1] * com[2];
  J[2][0] = J[0][2] += mass * com[2] * com[0];
  }

//*============================================================*
//*==========         pm_MathPolyVolIntComp          ==========*
//*============================================================*
// compute volume integral for a closed polygon mesh.

void
pm_MathPolyVolIntComp (int num_poly, int *conn, int hgeom, int num_verts, 
                       PmVector3 *verts, PmVector3 *fnorms,  PmMathVolIntData& vd)
  {

  double nx, ny, nz;

  int i, j, A, B, C, n; 

  vd.T0 = vd.T1[0] = vd.T1[1] = vd.T1[2]
        = vd.T2[0] = vd.T2[1] = vd.T2[2]
        = vd.TP[0] = vd.TP[1] = vd.TP[2] = 0;

  for (i = 0; i < num_poly; i++) {
    if (hgeom) {
      n = hgeom-1;
      j = *conn;
      }
    else {
      n = *conn;
      j = *(conn+1);
      }

    nx = fabs(fnorms[i][0]);
    ny = fabs(fnorms[i][1]);
    nz = fabs(fnorms[i][2]);

    if (nx + ny + nz == 0.0) {
      continue;
      }

    if ((nx > ny) && (nx > nz)) {
      C = 0;
      }
    else {
      C = (ny > nz) ? 1 : 2;
      }

    A = (C + 1) % 3;
    B = (A + 1) % 3;

    vd.A = A;
    vd.B = B;
    vd.C = C;
    vd.conn = conn;
    vd.w = -fnorms[i][0]*verts[j][0] -fnorms[i][1]*verts[j][1] -fnorms[i][2]*verts[j][2];

    pm_MathPolyFaceIntComp (verts, fnorms, i, hgeom, vd);

    vd.T0 += fnorms[i][0] * ((A == 0) ? vd.Fa : ((B == 0) ? vd.Fb : vd.Fc));

    vd.T1[A] += fnorms[i][A] * vd.Faa;
    vd.T1[B] += fnorms[i][B] * vd.Fbb;
    vd.T1[C] += fnorms[i][C] * vd.Fcc;
    vd.T2[A] += fnorms[i][A] * vd.Faaa;
    vd.T2[B] += fnorms[i][B] * vd.Fbbb;
    vd.T2[C] += fnorms[i][C] * vd.Fccc;
    vd.TP[A] += fnorms[i][A] * vd.Faab;
    vd.TP[B] += fnorms[i][B] * vd.Fbbc;
    vd.TP[C] += fnorms[i][C] * vd.Fcca;

    conn += n+1;
    }

  vd.T1[0] /= 2; vd.T1[1] /= 2; vd.T1[2] /= 2;
  vd.T2[0] /= 3; vd.T2[1] /= 3; vd.T2[2] /= 3;
  vd.TP[0] /= 2; vd.TP[1] /= 2; vd.TP[2] /= 2;
  }

//*============================================================*
//*==========         pm_MathPolyFaceIntComp         ==========*
//*============================================================*
// compute integral over polygon face.

void
pm_MathPolyFaceIntComp (PmVector3 *verts, PmVector3 *fnorms, int poly, 
                        int hgeom, PmMathVolIntData& vd)
  {
  double w;
  double k1, k2, k3, k4;
  PmVector3 n;
  int A, B, C;

  n = fnorms[poly];
  w = vd.w;
  A = vd.A;
  B = vd.B;
  C = vd.C;

  pm_MathPolyProjIntComp (verts, hgeom, vd);

  if (n[C] == 0.0) {
    return;
    }

  k1 = 1 / n[C]; 
  k2 = k1 * k1; 
  k3 = k2 * k1; 
  k4 = k3 * k1;

  vd.Fa = k1 * vd.Pa;
  vd.Fb = k1 * vd.Pb;
  vd.Fc = -k2 * (n[A]*vd.Pa + n[B]*vd.Pb + w*vd.P1);

  vd.Faa = k1 * vd.Paa;
  vd.Fbb = k1 * vd.Pbb;
  vd.Fcc = k3 * (SQR(n[A])*vd.Paa + 2*n[A]*n[B]*vd.Pab + SQR(n[B])*vd.Pbb
         + w*(2*(n[A]*vd.Pa + n[B]*vd.Pb) + w*vd.P1));

  vd.Faaa = k1 * vd.Paaa;
  vd.Fbbb = k1 * vd.Pbbb;
  vd.Fccc = -k4 * (CUBE(n[A])*vd.Paaa + 3*SQR(n[A])*n[B]*vd.Paab
           + 3*n[A]*SQR(n[B])*vd.Pabb + CUBE(n[B])*vd.Pbbb
           + 3*w*(SQR(n[A])*vd.Paa + 2*n[A]*n[B]*vd.Pab + SQR(n[B])*vd.Pbb)
           + w*w*(3*(n[A]*vd.Pa + n[B]*vd.Pb) + w*vd.P1));

  vd.Faab = k1 * vd.Paab;
  vd.Fbbc = -k2 * (n[A]*vd.Pabb + n[B]*vd.Pbbb + w*vd.Pbb);
  vd.Fcca = k3 * (SQR(n[A])*vd.Paaa + 2*n[A]*n[B]*vd.Paab + SQR(n[B])*vd.Pabb
         + w*(2*(n[A]*vd.Paa + n[B]*vd.Pab) + w*vd.Pa));
  }

//*============================================================*
//*==========       pm_MathPolyProjIntComp           ==========*
//*============================================================*

void
pm_MathPolyProjIntComp (PmVector3 *verts, int hgeom, PmMathVolIntData& vd)
  {
  double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;
  double a0, a1, da;
  double b0, b1, db;
  double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
  double a1_2, a1_3, b1_2, b1_3;
  double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
  double Cab, Kab, Caab, Kaab, Cabb, Kabb;

  int i, j;
  int A, B, C;
  int n, *conn;

  P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;
  A = vd.A;
  B = vd.B;
  C = vd.C;
  conn = vd.conn;

  if (hgeom) {
    n = hgeom;
    }
  else {
    n = conn[0];
    conn++;
    }

  for (i = 0; i < n; i++) {
    j = conn[i];
    a0 = verts[j][A];
    b0 = verts[j][B];

    if (i+1 > n-1) {
      j = conn[0];
      }  
    else {
      j = conn[i+1];
      }  

    a1 = verts[j][A];
    b1 = verts[j][B];

    da = a1 - a0;
    db = b1 - b0;
    a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
    b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
    a1_2 = a1 * a1; a1_3 = a1_2 * a1;
    b1_2 = b1 * b1; b1_3 = b1_2 * b1;

    C1 = a1 + a0;
    Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
    Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
    Cab = 3*a1_2 + 2*a1*a0 + a0_2; Kab = a1_2 + 2*a1*a0 + 3*a0_2;
    Caab = a0*Cab + 4*a1_3; Kaab = a1*Kab + 4*a0_3;
    Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
    Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

    P1 += db*C1;
    Pa += db*Ca;
    Paa += db*Caa;
    Paaa += db*Caaa;
    Pb += da*Cb;
    Pbb += da*Cbb;
    Pbbb += da*Cbbb;
    Pab += db*(b1*Cab + b0*Kab);
    Paab += db*(b1*Caab + b0*Kaab);
    Pabb += da*(a1*Cabb + a0*Kabb);
    }

  vd.P1 = P1 / 2.0;
  vd.Pa = Pa / 6.0;
  vd.Paa = Paa / 12.0;
  vd.Paaa = Paaa / 20.0;
  vd.Pb = Pb / -6.0;
  vd.Pbb = Pbb / -12.0;
  vd.Pbbb = Pbbb / -20.0;
  vd.Pab = Pab / 24.0;
  vd.Paab = Paab / 60.0;
  vd.Pabb = Pabb / -60.0;
  }

///////////////////////////////////////////////////////////////
//                         s t a t i s t i c s              //
/////////////////////////////////////////////////////////////

//*============================================================*
//*==========       pm_MathStandardDeviationComp     ==========*
//*============================================================*

void
pm_MathStandardDeviationComp (const vector<float>& vals, float& mean, float& sd)
  {
  float sum, s, m, d, x;
  int n = vals.size();
  sum = 0.0;
  x = 0.0;

  for (unsigned i = 0; i < n; i++) {
    sum += vals[i];
    }

  mean = sum / n;
  sum = 0.0;

  for (unsigned i = 0; i < n; i++) {
    x = vals[i] - mean;
    sum += x*x;
    }

  s = sum / (n-1);
  sd = sqrt(s);
  }


};


