
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

#ifndef _GC_GR_H_
#define _GC_GR_H_

#include "pm/gr/gr.h"

namespace PmGraphics {

typedef struct GrGcTri {
  GrVector3 p1, p2, p3;
  } GrGcTri;

typedef struct GrGcVolIntData {
  // projection integrals //
  double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;
  // face integrals //
  double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;
  // volume integrals //
  double T0, T1[3], T2[3], TP[3];
  int A, B, C;
  int *conn;
  float w;
  } GrGcVolIntData;

void PM_EXPORT
gc_BasisComp (GrVector3& normal, GrVector3& u, GrVector3& v);

void PM_EXPORT
gc_CylLineInt (float radius, float length, GrVector3& origin, GrVector3& axis,
               GrVector3 line[2], float tol, bool& intersect, GrVector3& ipt);

void PM_EXPORT
gc_LineLineDist(GrVector3& p1, GrVector3& v1, GrVector3& p2, GrVector3& v2,
                GrVector3& ipt, float& p_dist);

void PM_EXPORT
gc_LinePointProj (GrVector3 line[2], GrVector3& pt, float& dist, GrVector3& proj_pt);

void PM_EXPORT
gc_PolyFaceNormsComp (int num_poly, int *conn, int hgeom, int num_verts,
                       GrVector3 *verts, GrVector3 **p_fnorms);

void PM_EXPORT
gc_PolyVertexNormsComp (int num_poly, int *conn, int hgeom, int num_verts,
                        GrVector3 *verts, GrVector3 *fnorms, GrVector3 **p_vnorms);

void PM_EXPORT
gc_PolyLineInt (int n, GrVector3 *pts, GrVector3& norm, GrVector3 line[2], float tol,
                bool& p_intersect, GrVector3& ipt);

void PM_EXPORT
gc_PolyNormComp (int n, GrVector3 *pts, GrVector3& norm);

void PM_EXPORT
gc_PolyPointClassify (int n, GrVector3 *pts, GrVector3& norm, float tol, GrVector3& pt, 
                      bool& res);

void PM_EXPORT
gc_PolyPtIn3 (GrVector3 *tri, GrVector3& pt, int i0, int i1, int i2, bool& p_in);

void PM_EXPORT
gc_PolyWindingComp (int num_pts, GrVector3 *pts, GrVector3& normal, GrVector3& pt,
                    float& wnum);

void PM_EXPORT
gc_SphereGen (GrGcTri *f, int iterations, int *num);

void PM_EXPORT
gc_TriGetIndexes (int num_tri, GrVector3 *verts, int **p_conn, int *p_num_iverts,
                  GrVector3 **p_iverts);

void PM_EXPORT
gc_VecMaxProj (GrVector3& normal, int& i0, int& i1, int& i2);


}

#endif



