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
 * part:             p a r t i c l e   s y s t e m            *
 *============================================================*/

#ifndef _PARTICLE_SYSTEM_PM_H_
#define _PARTICLE_SYSTEM_PM_H_

#include "pm/pm.h"
#include "pobj.h"
#include "graphics.h"

namespace ProteinMechanica {

// PmParticle
// ----------

 class PM_EXPORT PmParticle : public PmPhysicalObj {

    public:
       PmParticle(const string name, int num_verts, PmVector3 *coords);
       ~PmParticle();
       void getCoordinates(vector<PmVector3>& coords){};
       void getDimensions(vector<float>& dims) {};
       void getExtent (PmExtent& extent ) { extent = this->extent; }
       void setExtent (PmExtent& extent ) { this->extent = extent; }
       void setMasses(float *masses);
       void getMassProps (PmMassProperties& props);
       void getName (string& name);
       void setName (const string name);
       void getRadii(vector<float>& rads) {};
       void defineRegion(const string name, const vector<int>& ids){};

       void setColor(const PmVector3& color);
       void display(const bool show);
       void setDisplaySpheres(const bool spheres);
       void setMapMass(const bool val);
       void setXform (PmXform& xform);
       void getXform (PmXform& xform);

    protected:

    private:
      string name;
      int num_coordinates;
      PmVector3 *coordinates;
      float *masses, *data, *radii;
      PmMassProperties mass_properties;
      PmExtent extent;
      PmXform xform;

      bool display_spheres, map_mass;
      PmGraphicsPoint  *point_geom;
      PmGraphicsSphere *sphere_geom;
      PmVector3 color, *map_colors;
      PmGeometryDisplayType display_type;
      PmGeometryShadingType shading_type;

   };

}

#endif

