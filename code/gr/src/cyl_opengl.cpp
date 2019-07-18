
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
//* cyl_opengl:        o p e n g l   c y l i n d e r           *
//*============================================================*

#include "pm/gr/cyl.h"
#include "gr_opengl.h"

namespace PmGraphics {

typedef struct CylinderGeom {
  int num_polys, *conn, num_verts; 
  GrVector3 *verts, *vnorms, *fnorms;
  } CylinderGeom;

typedef struct CylinderPrvData {
  CylinderGeom *geoms;
  GrVector3 *verts;
  } CylinderPrvData;

//*============================================================*
//*==========              api_setup                 ==========*
//*============================================================*

void
GrCylinder::api_setup()
  {
  //fprintf (stderr, "\n>>>>>>> GrCylinder::api_setup \n");

  // create private data  //

  if (!prv_data) {
    CylinderPrvData *prvd = new CylinderPrvData;

    float r, l; 
    GrVector3 center, dir, *verts, *vnorms, *fnorms, pt;
    int num_polys, *conn, num_verts;
    CylinderGeom *geoms = new CylinderGeom[num_vertices];

    for (int i = 0; i < num_vertices; i++) {
      pt = vertices[i];

      if (radii) {
        r = radii[i];
        }
      else {
        r = radius;
        }

      if (lengths) {
        l = lengths[i];
        }
      else {
        l = length;
        }

      if (axes) {
        dir = axes[i];
        }
      else {
        dir = axis;
        }

      /*
      fprintf (stderr, "   >>>> %d: pt[%g %g %g] r[%f] l[%f] dir[%g %g %g] \n", i, 
               pt[0], pt[1], pt[2], r, l, dir[0], dir[1], dir[2]); 
      */
    
      genCylinderGeom (r, l, pt, dir, &num_polys, &conn, &num_verts, &verts, &vnorms,
                       &fnorms);

      geoms[i].num_polys = num_polys; 
      geoms[i].conn = conn; 
      geoms[i].num_verts = num_verts;
      geoms[i].verts = verts, 
      geoms[i].fnorms = fnorms;
      geoms[i].vnorms = vnorms;
      }

    prvd->geoms = geoms; 
    prv_data = prvd;
    }

  glDisable (GL_LIGHTING);
  //glLineWidth (width);
  glColor3f (color[0], color[1], color[2]);
  }

//*============================================================*
//*==========              api_reset_state           ==========*
//*============================================================*

void
GrCylinder::api_reset_state()
  {
  glEnable (GL_LIGHTING);
  glLineWidth (1.0);
  glDisable (GL_LINE_STIPPLE);
  }

//*============================================================*
//*==========              api_render                ==========*
//*============================================================*
// render spheres.   

void
GrCylinder::api_render()
  {

  //fprintf (stderr, "\n>>>>>>> GrCylinder::api_render \n");
  GrVector3 v; 
  CylinderPrvData *prvd = (CylinderPrvData*)prv_data;

  bool vcolors;
  GLenum mode;
  int hgeom = 3;

  // get the display mode: point, line, fill.
  opengl_getPolygonMode (display, mode);

  //fprintf (stderr, "    >>> colors [%x] \n", colors); 
  //fprintf (stderr, "    >>> norms  [%x] \n", norms); 

  for (int i = 0; i < num_vertices; i++) {
    int num_polys = prvd->geoms[i].num_polys; 
    int *conn = prvd->geoms[i].conn;
    int num_verts = prvd->geoms[i].num_verts; 
    GrVector3 *verts = prvd->geoms[i].verts; 
    GrVector3 *norms; 

    // get normals for shading
    if (shading == GR_GEOMETRY_SHADING_COLOR) {
      glEnable (GL_LIGHTING);
      glShadeModel (GL_SMOOTH);
      norms = prvd->geoms[i].vnorms;
      }
    else if (shading == GR_GEOMETRY_SHADING_FLAT) {
      glEnable (GL_LIGHTING);
      norms = prvd->geoms[i].fnorms;
      }
    else if (shading == GR_GEOMETRY_SHADING_NONE) {
      norms = NULL;
      }

    if (colors) {
      glColor3fv ((GLfloat*)(colors+i));
      }

    if (norms) {
      opengl_renderPolygonsShading (mode, shading, num_polys, conn, verts, norms, 
                                    hgeom);
      }

    else {
      opengl_renderPolygons (mode, num_polys, conn, verts, hgeom);
      }
    }
  }


}  // PmGraphics

