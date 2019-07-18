
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
//* sphere_opengl:        o p e n g l   s p h e r e            *
//*============================================================*

#include "pm/gr/ellipsoid.h"
#include "gr_opengl.h"
#include "gc.h"

namespace PmGraphics {

typedef struct EllipsoidGeom {
  int num_polys, *conn, num_verts; 
  GrVector3 *verts, *fnorms, *vnorms, *tmp_verts;
  } EllipsoidGeom;

typedef struct EllipsoidPrvData {
  EllipsoidGeom *geom;
  } EllipsoidPrvData;

//*============================================================*
//*==========              api_setup                 ==========*
//*============================================================*

void
GrEllipsoid::api_setup()
  {
  //fprintf (stderr, "\n>>>>>>> GrEllipsoid::api_setup \n");

  // create private data which will be the geometry for each ellipsoid //

  if (!prv_data) {
    EllipsoidPrvData *prvd = new EllipsoidPrvData;
    prvd->geom = new EllipsoidGeom[num_vertices];
    //fprintf (stderr, "   >>>> gen ellipsoids ... ");

    for (int i = 0; i < num_vertices; i++) {
      genEllipsoidGeom(axes_data[i], vertices[i], &prvd->geom[i].num_polys, 
                       &prvd->geom[i].conn, &prvd->geom[i].num_verts, 
                       &prvd->geom[i].verts, &prvd->geom[i].fnorms);

      gc_PolyVertexNormsComp(prvd->geom[i].num_polys, prvd->geom[i].conn, 0,
                             prvd->geom[i].num_verts, prvd->geom[i].verts, 
                             prvd->geom[i].fnorms, &prvd->geom[i].vnorms);
      }

    //fprintf (stderr, " done \n");
    prv_data = prvd;
    }

  glDisable (GL_LIGHTING);
  glLineWidth (line_width);
  glColor3f (color[0], color[1], color[2]);
  }

//*============================================================*
//*==========              api_reset_state           ==========*
//*============================================================*

void
GrEllipsoid::api_reset_state()
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
GrEllipsoid::api_render()
  {

  //fprintf (stderr, "\n>>>>>>> GrEllipsoid::api_render \n");

  GrVector3 v; 
  EllipsoidPrvData *prvd = (EllipsoidPrvData*)prv_data;
  bool vcolors;
  GLenum mode;
  GrVector3 *norms;

  // get the display mode: point, line, fill.
  opengl_getPolygonMode (display, mode);
  int hgeom = 0;

  for (int i = 0; i < num_vertices; i++) {
    int num_polys = prvd->geom[i].num_polys; 
    int *conn = prvd->geom[i].conn;
    int num_verts = prvd->geom[i].num_verts; 
    GrVector3 *verts = prvd->geom[i].verts; 
    GrVector3 *fnorms = prvd->geom[i].fnorms;
    GrVector3 *vnorms = prvd->geom[i].vnorms;

    // set lighting for shading
    if (shading == GR_GEOMETRY_SHADING_COLOR) { 
      glEnable (GL_LIGHTING);
      glShadeModel (GL_SMOOTH);
      norms = vnorms;
      }
    else if (shading == GR_GEOMETRY_SHADING_FLAT) { 
      glEnable (GL_LIGHTING);
      glShadeModel (GL_FLAT);
      norms = fnorms;
      }
    else if (shading == GR_GEOMETRY_SHADING_NONE) { 
      norms = NULL;
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

