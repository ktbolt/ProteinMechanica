
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
//* rest:                     r e s t r a i n t                *
//*============================================================*

#ifndef _RESTRAINT_PM_H_
#define _RESTRAINT_PM_H_

#include "pm/pm.h"

namespace ProteinMechanica {

class PmGraphicsGeometry;
class PmSimulationObj;

typedef enum PmRestraintType {
  PM_RESTRAINT_UNKNOWN,
  PM_RESTRAINT_ALL,
  PM_RESTRAINT_CENTER
  } PmRestraintType;


// PmRestraint
// -----------

class PM_EXPORT PmRestraint {

    public:
       PmRestraint(const string name);
       bool isActive(float time);
       void getName(string& name);
       bool hasEnergy();
       void setCompEnergy(const bool flag);
       void getEnergy(float& energy);
       void getRegionCenter(const int rgn_num, PmVector3& center);
       void getResidual(float& residual);
       void getPointDistance(float& dist, float& pdist); 
       void compForces(const bool upd, const float time);
       bool getPowerParams(float power_params[2]);
       void setPowerParams(const float power_params[2]);
       void setRampParams(const float param);
       bool getRampParams(const float time, float& dramp);
       void setRegions(const string rgn1, const string rgn2);
       void getSimObjs(PmSimulationObj **sobj1, PmSimulationObj **sobj2);
       void setSimObjs(PmSimulationObj *sobj1, PmSimulationObj *sobj2);
       void getSpringData(float *k, float *dist); 
       void setSpringData(const float k, const float dist); 
       static void convType(const string str, PmRestraintType& type);
       void setType(const PmRestraintType type);
       void getAbsoluteDistance(bool& flag);
       void setAbsoluteDistance(const bool flag); 
       void setTime(PmTimeInterval& time);
       void setColor(const PmVector3 color); 
       void display();
       void setWidth(const float val); 
       void setShow(const bool val); 

    private:
       string name;
       PmRestraintType type;
       bool initialized;
       bool active;
       PmSimulationObj *sobjs[2];
       string regions[2]; 
       float distance, force_const;
       bool compute_energy; 
       float current_energy;
       float current_residual;
       bool use_absolute_distance;
       bool has_power_params;
       float power_params[2];
       bool has_ramp, ramp_init;
       float ramp_start, dramp, ramp_time; 
       vector<PmVector3> points1, points2, current_points1, current_points2;

       PmTimeInterval time_interval;
       PmGraphicsGeometry *line_geometry, *point_geometry; 
       PmVector3 *geometry_vertices; 
 
       PmVector3 color;
       float width;
       bool show;

       void xformPoints();
       void init();
   };

}

#endif

