
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
//* gr:                g r a p h i c s                         *
//*============================================================*

#ifndef _GR_H_
#define _GR_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <vector>

#include "pm/common.h"

using namespace std; 
using namespace ProteinMechanica; 

namespace PmGraphics {

typedef PmMatrix3x3 GrMatrix3x3;
typedef PmVector3 GrVector3;
typedef PmExtent GrExtent;
typedef PmXform GrXform;
typedef PmConn GrIndex;


//////////////////////////////////////////////////////////////////
//                   b a s i c    t y p e s                    //
////////////////////////////////////////////////////////////////

// window types  
#define GR_WINDOW_TYPE_TEMPLATE    -1
#define GR_WINDOW_TYPE_RGB         0
#define GR_WINDOW_TYPE_RGBA        GR_WINDOW_TYPE_RGB
#define GR_WINDOW_TYPE_INDEX       1
#define GR_WINDOW_TYPE_SINGLE      0
#define GR_WINDOW_TYPE_DOUBLE      2
#define GR_WINDOW_TYPE_DIRECT      0
#define GR_WINDOW_TYPE_INDIRECT    4

#define GR_WINDOW_TYPE_ACCUM       8
#define GR_WINDOW_TYPE_ALPHA       16
#define GR_WINDOW_TYPE_DEPTH       32
#define GR_WINDOW_TYPE_STENCIL     64
#define GR_WINDOW_TYPE_VGR_WINDOW_TYPE         128

typedef enum {
  GR_GEOMETRY_UNKNOWN,
  GR_GEOMETRY_LINES,
  GR_GEOMETRY_QUADS,
  GR_GEOMETRY_TRIANGLES,
  GR_GEOMETRY_TRISTRIP,
  GR_GEOMETRY_GRID
  } GrGeometryType;

typedef enum {
  GR_RENDER_MODE_UNKNOWN,
  GR_RENDER_MODE_NONE,
  GR_RENDER_MODE_SURFACE,
  GR_RENDER_MODE_VOLUME,
  GR_RENDER_MODE_RAYTRACE
  } GrRenderMode;

typedef enum {
  GR_GEOMETRY_DISPLAY_UNKNOWN,
  GR_GEOMETRY_DISPLAY_POINT,
  GR_GEOMETRY_DISPLAY_LINE,
  GR_GEOMETRY_DISPLAY_SOLID,
  GR_GEOMETRY_DISPLAY_OUTLINE
  } GrGeometryDisplayType;

typedef enum {
  GR_GEOMETRY_SHADING_NONE,
  GR_GEOMETRY_SHADING_FLAT,
  GR_GEOMETRY_SHADING_COLOR,
  GR_GEOMETRY_SHADING_NORMAL
  } GrGeometryShadingType;

// line styles

typedef enum {
  GR_LINE_STYLE_UNKNOWN,
  GR_LINE_STYLE_SOLID,
  GR_LINE_STYLE_DOT,
  GR_LINE_STYLE_DASH,
  GR_LINE_STYLE_DOT_DASH,
  GR_LINE_STYLE_MAX
  } GrLineStyleType;

typedef enum {
  GR_IMAGE_FORMAT_UNKNOWN,
  GR_IMAGE_FORMAT_RGB,
  GR_IMAGE_FORMAT_PPM,
  GR_IMAGE_FORMAT_JPEG
  } GrImageFormat;

typedef enum {
  GR_GEOMETRY_MARKER_NONE,
  GR_GEOMETRY_MARKER_POINT,
  GR_GEOMETRY_MARKER_CROSS,
  GR_GEOMETRY_MARKER_SQUARE,
  GR_GEOMETRY_MARKER_BOX,
  GR_GEOMETRY_MARKER_CIRCLE,
  GR_GEOMETRY_MARKER_CIRCLE_FILLED
  } GrGeometryMarkerType;

class GrWindowRecordParams {
  public:
    GrWindowRecordParams() { format = GR_IMAGE_FORMAT_UNKNOWN; }
    GrImageFormat format;
    string prefix;
    GrWindowRecordParams &operator=(const GrWindowRecordParams& rhs) {
      format = rhs.format;
      prefix = rhs.prefix;
      }
  };

typedef float (GrReal3)[3];

typedef char (GrCmdPrompt)[20];

//*==============================================================*
//*==========              color                      ===========*
//*==============================================================*

class GrColor {
  public:
    float val[3];

    GrColor() { val[0] = 1.0; val[1] = 1.0; val[2] = 1.0; }
    GrColor (const float v[3]) { val[0] = v[0]; val[1] = v[1]; val[2] = v[2]; }
    GrColor(float r, float g, float b) { val[0] = r; val[1] = g; val[2] = b; }

    void set (float r, float g, float b) { val[0] = r; val[1] = g; val[2] = b; }
    void set (const GrVector3 v) { val[0] = v[0]; val[1] = v[1]; val[2] = v[2]; }
    void get (float v[3]) { v[0] = val[0]; v[1] = val[1]; v[2] = val[2]; }

    GrColor &operator=(const GrVector3& v) { 
       val[0] = v[0], val[1] = v[1], val[2] = v[2]; 
       }

    float &operator [] (const int i) { return val[i]; }
    const float &operator [] (const int i) const { return val[i]; }
  };

class GrAlphaColor {
  public:
    float val[4];

    GrAlphaColor() { val[0] = 1.0; val[1] = 1.0; val[2] = 1.0; val[3] = 1.0; }
    GrAlphaColor (float r, float g, float b, float a) { 
      val[0] = r; val[1] = g; val[2] = b; val[3] = a; }

    GrAlphaColor (const float v[4]) { 
       val[0]=v[0]; val[1]=v[1]; val[2]= v[2]; val[3] = v[3]; }

    void get (float v[4]) {
      v[0] = val[0]; v[1] = val[1]; v[2] = val[2]; v[3] = val[3]; }
  };

typedef vector<GrColor> GrColorVector;

const GrColor grWhite(1,1,1);
const GrColor grRed(1,0,0);
const GrColor grGreen(0,1,0);
const GrColor grBlue(0,0,1);

//////////////////////////////////////////////////////////////////
//                         p i c k i n g                       //
////////////////////////////////////////////////////////////////

typedef enum {
  GR_PICK_UNKNOWN,
  GR_PICK_SELECT,
  GR_PICK_POINT,
  GR_PICK_VERTEX
  } GrPickType;

typedef struct GrPickAtts {
  float tol;
  GrPickType type;
  } GrPickAtts;

typedef struct GrPickGeom {
  class GrVertexGeometry *geom;
  bool intersected;
  GrVector3 int_pt;
  float dist;
  } GrPickGeom;

typedef struct GrPickResult {
  int id;
  GrVector3 line[2];
  GrPickAtts atts;
  //GrPickGeom *geom_list;
  vector<class GrPickGeom> geom_list;
  } GrPickResult;

typedef void  (*GrPickReturnFunc)(class GrScene *scene, GrPickResult& pick_res);

//////////////////////////////////////////////////////////////////
//                         l i g h t                           //
////////////////////////////////////////////////////////////////

// light types  

typedef enum {
  GR_LIGHT_DIRECTIONAL = 1,
  GR_LIGHT_POSITIONAL  = 2
  } GrLightType;


//////////////////////////////////////////////////////////////////
//                         g r                                 //
////////////////////////////////////////////////////////////////


// no copy class 

class GrNoCopy {
  protected:
    GrNoCopy() {};
    ~GrNoCopy() {};

  private:
    GrNoCopy(const GrNoCopy&);
    GrNoCopy& operator=(const GrNoCopy&);
  };


// graphics system 

class PM_EXPORT GrSystem {
  public:
     GrSystem();
     void addWindow (class GrWindow *win);
     void setCommandCallback (void (*func)(FILE *fp, const string cmd));
     void setCommandCallback1 (void (*func)(const string cmd));
     void processCommand (const string line);

  private:
     vector<class GrWindow*> windows;
     void (*cmd_callback_func)(FILE *fp, const string cmd);
     void (*cmd_callback_func1)(const string cmd);
  };

extern PM_EXPORT GrSystem grSystem;

void gr_RotationFromAngles (float xr, float yr, float zr, PmMatrix3x3& mat);

}


#endif



