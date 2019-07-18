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

/*============================================================*
 * common:          c o m m o n   t y p e s                   *
 *============================================================*/

#ifndef _PM_COMMON_H_
#define _PM_COMMON_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace ProteinMechanica {

// set api declarations needed for windows

//#define _PM_CYGWIN_

#if defined(_PM_WIN32_) || defined (_PM_CYGWIN_)
  #define  PM_EXPORT __declspec( dllexport )
  #define  PM_IMPORT __declspec(dllimport)
#else
  #define  PM_EXPORT
#endif


//*==============================================================*
//*==========                3d vector                ===========*
//*==============================================================*

class PmVector3 {
  private:
    float val[3];
  public:
    PmVector3 () { val[0] = 0.0; val[1] = 0.0; val[2] = 0.0; }
    PmVector3 (float x, float y, float z) { val[0] = x; val[1] = y; val[2] = z; }
    PmVector3 (const PmVector3& v) { val[0] = v[0]; val[1] = v[1]; val[2] = v[2]; }
    void set  (float x, float y, float z) { val[0] = x; val[1] = y; val[2] = z; }
    void set  (float v[3]) { val[0] = v[0]; val[1] = v[1]; val[2] = v[2]; }

    float &operator [] (const int i) { return val[i]; }
    const float &operator [] (const int i) const { return val[i]; }

    PmVector3 operator - () { return PmVector3(-val[0], -val[1], -val[2]); }

    PmVector3 operator+(PmVector3 o) {
      return PmVector3(val[0]+o[0], val[1]+o[1], val[2]+o[2]);
      }

    /*
    PmVector3 &operator=(const PmVector3& v) {
       val[0] = v[0], val[1] = v[1], val[2] = v[2];
       return *this;
       }
    */

    PmVector3 operator - (PmVector3 o) const {
      return PmVector3(val[0]-o[0], val[1]-o[1], val[2]-o[2]);
      }

    PmVector3 operator * (float s) {
      return PmVector3(s*val[0], s*val[1], s*val[2]);
      }

    friend inline PmVector3 operator * (const float s, PmVector3 v) { return ( v*s ); }
    float length() { return sqrt(val[0]*val[0] + val[1]*val[1] + val[2]*val[2]); }

    PmVector3 operator / (float s) {
      return PmVector3(val[0]/s, val[1]/s, val[2]/s);
      }

    friend inline PmVector3 operator / (const float s, PmVector3 v) { return ( v/s ); }

    // dot product //
    float operator*(PmVector3 o) const {
      return val[0]*o[0]+val[1]*o[1]+val[2]*o[2];}

    void normalize() {
      float mag = length();
      if (mag == 0.0) return;
      val[0] /= mag; val[1] /= mag; val[2] /= mag;
      }

    PmVector3 cross(const PmVector3& v) const {
      return PmVector3(val[1]*v[2] - val[2]*v[1], val[2]*v[0] - val[0]*v[2], 
                       val[0]*v[1] - val[1]*v[0]);
      }
   };

//*==============================================================*
//*==========                quaternion               ===========*
//*==============================================================*

class PmQuaternion {
  private:
    float val[4];
  public:
    PmQuaternion() { val[0] = 1.0; val[1] = 0.0; val[2] = 0.0; val[3] = 0.0; }
    float &operator [] (const int i) { return val[i]; }
    const float &operator [] (const int i) const { return val[i]; }
    /*
    PmQuaternion &operator=(const PmQuaternion& q) {
       val[0] = q[0], val[1] = q[1], val[2] = q[2]; val[3] = q[3];
       }
    */
  };

//*==============================================================*
//*==========               3x3 matrix                ===========*
//*==============================================================*

class PmMatrix3x3 {

   public:

     float val[3][3];

     PmMatrix3x3() {
       val[0][0] = 0.0; val[0][1] = 0.0; val[0][2] = 0.0;
       val[1][0] = 0.0; val[1][1] = 0.0; val[1][2] = 0.0;
       val[2][0] = 0.0; val[2][1] = 0.0; val[2][2] = 0.0;
       }

     PmMatrix3x3(const float v) {
       val[0][0] = v; val[0][1] = v; val[0][2] = v;
       val[1][0] = v; val[1][1] = v; val[1][2] = v;
       val[2][0] = v; val[2][1] = v; val[2][2] = v;
       }

     PmMatrix3x3(const float d1, const float d2, const float d3) {
       val[0][0] = d1;  val[0][1] = 0.0; val[0][2] = 0.0;
       val[1][0] = 0.0; val[1][1] = d2;  val[1][2] = 0.0;
       val[2][0] = 0.0; val[2][1] = 0.0; val[2][2] = d3;  
       }

     float &operator () (const int i, const int j) { return val[i][j]; }
     const float &operator () (const int i, const int j) const { return val[i][j]; }

     /*
     PmMatrix3x3 &operator=(const PmMatrix3x3& rhs) {
       val[0][0] = rhs(0,0); val[0][1] = rhs(0,1); val[0][2] = rhs(0,2);
       val[1][0] = rhs(1,0); val[1][1] = rhs(1,1); val[1][2] = rhs(1,2);
       val[2][0] = rhs(2,0); val[2][1] = rhs(2,1); val[2][2] = rhs(2,2);
       }
     */

     void diag(const float c) {
       val[0][0] =   c; val[0][1] = 0.0; val[0][2] = 0.0;
       val[1][0] = 0.0; val[1][1] =   c; val[1][2] = 0.0;
       val[2][0] = 0.0; val[2][1] = 0.0; val[2][2] =   c;
       }

     void setRow(int c, float c1, float c2, float c3 ) {
       val[c][0] = c1; val[c][1] = c2; val[c][2] = c2;
       }

     void transpose() {
       float tmp;
       tmp = val[0][1]; val[0][1] = val[1][0]; val[1][0] = tmp;
       tmp = val[0][2]; val[0][2] = val[2][0]; val[2][0] = tmp;
       tmp = val[1][2]; val[1][2] = val[2][1]; val[2][1] = tmp;
       }

     // matrix / vector multiply //

     PmVector3 operator * (const PmVector3& v) const {
       float v1 = v[0]*val[0][0] + v[1]*val[0][1] + v[2]*val[0][2];
       float v2 = v[0]*val[1][0] + v[1]*val[1][1] + v[2]*val[1][2];
       float v3 = v[0]*val[2][0] + v[1]*val[2][1] + v[2]*val[2][2];
       return PmVector3(v1, v2, v3);
       }

     // scalar multiply //

     PmMatrix3x3 operator * (const float s) const {
       PmMatrix3x3 C;
       for (int i = 0; i < 3; i++) 
         for (int j = 0; j < 3; j++) 
           C(i,j) = s*val[i][j];
       return C; 
       }

     friend PmMatrix3x3 operator * (const float s, const PmMatrix3x3& rhs) {
       return (rhs*s);
       }

     // matrix multiply //

     PmMatrix3x3 operator * (const PmMatrix3x3& rhs) const {
        PmMatrix3x3 C;
        C(0,0) = val[0][0]*rhs(0,0) + val[0][1]*rhs(1,0) + val[0][2]*rhs(2,0); 
        C(0,1) = val[0][0]*rhs(0,1) + val[0][1]*rhs(1,1) + val[0][2]*rhs(2,1); 
        C(0,2) = val[0][0]*rhs(0,2) + val[0][1]*rhs(1,2) + val[0][2]*rhs(2,2); 
        C(1,0) = val[1][0]*rhs(0,0) + val[1][1]*rhs(1,0) + val[1][2]*rhs(2,0); 
        C(1,1) = val[1][0]*rhs(0,1) + val[1][1]*rhs(1,1) + val[1][2]*rhs(2,1); 
        C(1,2) = val[1][0]*rhs(0,2) + val[1][1]*rhs(1,2) + val[1][2]*rhs(2,2); 
        C(2,0) = val[2][0]*rhs(0,0) + val[2][1]*rhs(1,0) + val[2][2]*rhs(2,0); 
        C(2,1) = val[2][0]*rhs(0,1) + val[2][1]*rhs(1,1) + val[2][2]*rhs(2,1); 
        C(2,2) = val[2][0]*rhs(0,2) + val[2][1]*rhs(1,2) + val[2][2]*rhs(2,2); 
        return C;
        }

   };


//*==============================================================*
//*==========               3d extent                 ===========*
//*==============================================================*

class PmExtent {

  public:

    PmVector3 min, max;

    PmExtent() : min(0,0,0), max(1,1,1) { }

    PmExtent(const float& xmin, const float& xmax, const float& ymin, const float& ymax,
             const float& zmin, const float& zmax) : min(xmin,ymin,zmin), 
                                                     max(xmax, ymax, zmax) { }

    PmExtent(PmVector3& vmin, PmVector3& vmax) : min(vmin), max(vmax) { }

    void set(PmVector3& min, PmVector3& max) {
       this->min = min;
       this->max = max;
       }

    void get(PmVector3& min, PmVector3& max) {
       min = this->min;
       max = this->max;
       }

    void update(PmVector3& v) {
      if (v[0] < min[0]) min[0] = v[0];
      if (v[1] < min[1]) min[1] = v[1];
      if (v[2] < min[2]) min[2] = v[2];
      if (v[0] > max[0]) max[0] = v[0];
      if (v[1] > max[1]) max[1] = v[1];
      if (v[2] > max[2]) max[2] = v[2];
      }

    void update(PmExtent& src) {
      if (src.min[0] < min[0]) min[0] = src.min[0];
      if (src.min[1] < min[1]) min[1] = src.min[1];
      if (src.min[2] < min[2]) min[2] = src.min[2];

      if (src.max[0] > max[0]) max[0] = src.max[0];
      if (src.max[1] > max[1]) max[1] = src.max[1];
      if (src.max[2] > max[2]) max[2] = src.max[2];
      }
  };

//*==============================================================*
//*==========                 xform                   ===========*
//*==============================================================*

class PmXform {
  public:
    PmXform () : set(false), center(0,0,0), matrix(1,1,1), scale(1,1,1), 
                 translation(0,0,0), angles(0, 0, 0) {} 

    bool set;
    //bool rot_order_zyx;
    PmVector3 center;
    PmMatrix3x3 matrix;
    PmVector3 scale;
    PmVector3 translation;
    PmVector3 angles;

    void init() {
      set = false;
      //rot_order_zyx = false;
      scale.set (1, 1, 1);
      translation.set (0, 0, 0);
      angles.set (0, 0, 0);
      matrix.diag(1.0);
      }

    void setRotation(float x, float y, float z);
    void setRotation(float angle, PmVector3& axis);
    void setRotation(PmMatrix3x3& mat) { set = true; matrix = mat; }

    void resetRotation() {
      matrix.diag(1.0);
      }
  };

//*==============================================================*
//*==========                 conn                    ===========*
//*==============================================================*
// connectivity.

class PmConn {
   public:
     PmConn(int n) {
       size = n;
       vals = new int [n];
       }

      int &operator [] (const int i) { return vals[i]; }
      const int &operator [] (const int i) const { return vals[i]; }

     int size;
     int *vals;
  };

}

#endif


