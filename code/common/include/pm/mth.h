
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

#ifndef _MATH_PM_H_
#define _MATH_PM_H_

#include "common.h"
#include <math.h>

namespace ProteinMechanica {

#ifndef M_PI
#define M_PI  3.14159265358979323846 
#endif

class PmPcaResults {
  public:
    PmPcaResults() { s1 = 0.0; s2 = 0.0; s3 = 0.0; }
    float s1, s2, s3;
    PmVector3 com, axis1, axis2, axis3;
  };

typedef struct PmMathVolIntData {
  // projection integrals // 
  double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;
  // face integrals //
  double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;
  // volume integrals //
  double T0, T1[3], T2[3], TP[3];
  int A, B, C;
  int *conn;
  float w;
  } PmMathVolIntData;

typedef struct PmMathContactEllipsoids {
  PmVector3 pos1, u1, v1, w1;
  PmVector3 pos2, u2, v2, w2;
  float r1[3], max_r1, r2[3], max_r2;     
  float lambda, cfunc;
  PmVector3 npt1, npt2, cpt; 
  PmMatrix3x3 A, B;
  } PmMathContactEllipsoids;

void 
pm_MathBasisCompute (PmVector3& normal, PmVector3& u, PmVector3& v);

void 
pm_MathEcfCompContact (PmMathContactEllipsoids& ellipsoids, bool& contact);

float
pm_MathEcfDerivEval (PmMatrix3x3& A, PmMatrix3x3& B, PmVector3& R, float val);

float
pm_MathEcfEval (PmMatrix3x3& A, PmMatrix3x3& B, PmVector3& R, float val);

void
pm_MathEcfInterp (PmMatrix3x3& A, PmMatrix3x3& B, float v1, float v2, PmMatrix3x3& G);

void
pm_MathEcfContactFuncSolve (PmMatrix3x3& A, PmMatrix3x3& B, PmVector3& R,
                            float *p_lambda, float *p_cfunc);

void
pm_MathEcfContactPtComp (PmMatrix3x3& A, PmMatrix3x3& B, PmVector3& r, PmVector3& s,
                         float lambda, float F, PmVector3& xa, PmVector3& xb, 
                         PmVector3& cpt);

void
pm_MathEcfQuadForm (float *r, PmVector3& u, PmVector3& v, PmVector3& w, PmMatrix3x3& A);

void
pm_MathEigMatrix3x3 (float a[3][3], float w[3], float v[3][3], bool& singular);

void
pm_MathGramSchmidt (PmVector3& v1, PmVector3& v2, PmVector3& v3, PmVector3& q1,
                    PmVector3& q2, PmVector3& q3);

void
pm_MathMatrix3x3Inverse(PmMatrix3x3& src, PmMatrix3x3& res, bool *singular);

void
pm_MathMatrix3x3Solve (PmMatrix3x3& src, PmVector3& rhs, PmVector3& res, bool *singular);

void 
pm_MathExtractRotations(float mat[3][3], float rot[3]);

void
pm_MathFitRms (vector<PmVector3>& v1, vector<PmVector3>& v2, vector<float> wt,
               PmXform& xform, float& p_err);

void
pm_MathMatrixToQuaternion(PmMatrix3x3& m, PmQuaternion q);

void
pm_MathMatrix3x3DirCos (PmVector3& u1, PmVector3& u2, PmVector3& u3,
                        PmVector3& v1, PmVector3& v2, PmVector3& v3,
                        PmMatrix3x3& mat);

void
pm_MathFrameRotation(PmVector3& u1, PmVector3& u2, PmVector3& u3,
                     PmVector3& v1, PmVector3& v2, PmVector3& v3, PmMatrix3x3& mat);

void
pm_MathQuaternionToMatrix(PmQuaternion q, PmMatrix3x3& m);

void
pm_MathPrincipalComponents(vector<PmVector3>& coords, PmPcaResults& pca);

void
pm_MathPrincipalComponents(int num_coords, PmVector3 *coords, PmPcaResults& pca);

void
pm_MathPrincipalComponentsProj(vector<PmVector3>& coords, const PmPcaResults& pca,
                               PmExtent& proj, float w[3]);

void
pm_MathPrincipalComponentsProj(const int num_coords, const PmVector3 *coords,
                               const PmPcaResults& pca, PmExtent& proj, float w[3]);

float
pm_MathRandNumber (int long *seed);

void
pm_MathRandGaussVector(const float mean, const float sd, PmVector3& vec, int long *seed);

void 
pm_MathRotationAroundAxis (float ang, PmVector3 axis, PmMatrix3x3& mat);

void
pm_MathRotationFromAngles (float xr, float yr, float zr, PmMatrix3x3& mat);

double
pm_MathSign(double a);

void
pm_MathPolyClosed(int num_polys, int *conn, int hgeom, int num_verts, bool& closed);

void 
pm_MathPolyMassProps (int num_poly, int *conn, int hgeom, int num_verts,
                      PmVector3 *verts, PmVector3 *face_norms, float density,
                      float& mass, PmVector3& com, float J[3][3]);

void 
pm_MathPolyVolIntComp (int num_poly, int *conn, int hgeom, int num_verts,
                       PmVector3 *verts, PmVector3 *fnorms, PmMathVolIntData& vd);

void 
pm_MathPolyFaceIntComp (PmVector3 *verts, PmVector3 *fnorms, int poly,
                        int hgeom, PmMathVolIntData& vd);

void 
pm_MathPolyProjIntComp (PmVector3 *verts, int hgeom, PmMathVolIntData& vd);

void
pm_MathStandardDeviationComp (const vector<float>& vals, float& mean, float& sd);


}


#endif



