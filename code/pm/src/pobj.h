
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
 * pobj:             p h y s i c a l   o b j e c t            *
 *============================================================*/

#ifndef _PHYSICAL_OBJECT_PM_H_
#define _PHYSICAL_OBJECT_PM_H_

#include "pm/pm.h"
#include "rgn.h"

namespace ProteinMechanica {

class PmSimulationObj;


// PmPhysicalObj
// -------------

 class PM_EXPORT PmPhysicalObj {

    public:

       PmPhysicalObj();
       ~PmPhysicalObj();
       static vector<PmPhysicalObj*> physical_objects;

       void getName(string& name);
       PmObjectType getObjectType() { return obj_type; };
       void getObjectType(string& str);
       virtual void getDimensions(vector<float>& dims) = 0;
       virtual void getCoordinates(vector<PmVector3>& coords) = 0;
       virtual void getRadii(vector<float>& rads) = 0;
       virtual void getMassProps(PmMassProperties& props) = 0;
       virtual void display(const bool show) = 0;
       virtual void getXform (PmXform& xform) = 0;
       virtual void setXform (PmXform& xform) = 0;

       void createHbondRegions(const string rgn_name, PmPhysicalObj *pobj, 
                               const float cutoff, const string chain, 
                               vector<PmRegion*>& rgns); 

       void createRegions(const string rgn_name, PmPhysicalObj *pobj, 
                          const float cutoff, bool use_sidechains, PmAtomFilter& filter,
                          vector<PmRegion*>& rgns); 

       virtual void defineRegion(const string name, const vector<int>& ids)=0;

       void getRegion(const string rgn_name, PmRegion **rgn); 
       void getRegionCenter(const string rgn_name, PmVector3& center); 
       void getRegionCoords(const string rgn_name, vector<PmVector3>& coords);
       void getRegionData(const string rgn_name, vector<float>& data);
       void getRegionIndices(const string rgn_name, vector<int>& ids);
       void getRegionRadii(const string rgn_name, vector<float>& rads);
       bool hasRegion(const string rgn_name); 

       void setSimulationObj(PmSimulationObj *sobj);
       void getSimulationObj(PmSimulationObj **sobj);

       static bool get(string& name, PmPhysicalObj **pobj);

    protected:
       PmObjectType obj_type;
       vector <PmRegion*> regions;
       PmSimulationObj *sobj;
   };

}

#endif

