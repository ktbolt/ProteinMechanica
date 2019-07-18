
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
//* gc:         g e o m e t r i c    c o m p u t a t i o n     *
//*============================================================*

#ifndef _GC_PM_H_
#define _GC_PM_H_

#include "pm/pm.h"

namespace ProteinMechanica {

typedef struct PmGcTri {
  PmVector3 p1, p2, p3;
  } PmGcTri;

typedef struct PmGcVolIntData {
  // projection integrals //
  double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;
  // face integrals //
  double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;
  // volume integrals //
  double T0, T1[3], T2[3], TP[3];
  int A, B, C;
  int *conn;
  float w;
  } PmGcVolIntData;

void PM_EXPORT
pm_GcBasisComp (PmVector3& normal, PmVector3& u, PmVector3& v);

void PM_EXPORT
pm_GcCylLineInt (float radius, float length, PmVector3& origin, PmVector3& axis,
               PmVector3 line[2], float tol, bool& intersect, PmVector3& ipt);

void PM_EXPORT
pm_GcLineLineDist(PmVector3& p1, PmVector3& v1, PmVector3& p2, PmVector3& v2,
                PmVector3& ipt, float& p_dist);

void PM_EXPORT
pm_GcLinePointProj (PmVector3 line[2], PmVector3& pt, float& dist, PmVector3& proj_pt);

void PM_EXPORT
pm_GcPolyFaceNormsComp (int num_poly, int *conn, int hgeom, int num_verts,
                        PmVector3 *verts, PmVector3 **p_fnorms);

void PM_EXPORT
pm_GcPolyVertexNormsComp (int num_poly, int *conn, int hgeom, int num_verts,
                          PmVector3 *verts, PmVector3 *fnorms, PmVector3 **p_vnorms);

void PM_EXPORT
pm_GcPolyLineInt (int n, PmVector3 *pts, PmVector3& norm, PmVector3 line[2], float tol,
                bool& p_intersect, PmVector3& ipt);

void PM_EXPORT
pm_GcPolyNormComp (int n, PmVector3 *pts, PmVector3& norm);

void PM_EXPORT
pm_GcPolyPointClassify (int n, PmVector3 *pts, PmVector3& norm, float tol, 
                        PmVector3& pt, bool& res);

void PM_EXPORT
pm_GcPolyPtIn3 (PmVector3 *tri, PmVector3& pt, int i0, int i1, int i2, bool& p_in);

void PM_EXPORT
pm_GcPolyWindingComp (int num_pts, PmVector3 *pts, PmVector3& normal, PmVector3& pt,
                      float& wnum);

void PM_EXPORT
pm_GcSphereGen (PmGcTri *f, int iterations, int *num);

void PM_EXPORT
pm_GcTriGetIndexes (int num_tri, PmVector3 *verts, int **p_conn, int *p_num_iverts,
                    PmVector3 **p_iverts);

void PM_EXPORT
pm_GcVecMaxProj (PmVector3& normal, int& i0, int& i1, int& i2);

void
pm_GcTriInscribedCircle(PmVector3 pts[3], PmVector3& center, float& radius, float& area);


}

#endif



