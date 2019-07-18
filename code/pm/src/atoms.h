
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
 * atoms:                 a t o m s                           *
 *============================================================*/

#include "pm/pm.h"
#include "pm/mth.h"

#ifndef _ATOMS_PM_H_
#define _ATOMS_PM_H_

namespace ProteinMechanica {

#include <vector>  

using namespace std; 


// PmAtoms
// -----------------

 class PM_EXPORT PmAtoms {

   public:
     PmAtoms();
     PmAtoms(int n);
     ~PmAtoms();

     virtual void addAtom (const int id, PmAtomName name, PmAtomElementType type, 
                           PmVector3 pos, char chain_id, PmCompoundType ctype, 
                           PmMoleculeType mtype, int res_seq);
     virtual void addAtom (PmAtom *atom); 

     int getNumAtoms();
     int getNumAtoms (PmAtomName name);

     int getNumChains();

     void getAtom (const int id, PmAtom **atom);
     void getAtoms (PmAtomName name, vector<PmAtom*>& sel_atoms);
     void getAtoms (vector<PmAtom*>& sel_atoms);

     void getAtomCoords (int *num, PmVector3 **coords);
     void getAtomCoords (PmAtomName name, int *num_coords, PmVector3 **coords);
     void getAtomCoords (vector<PmVector3> &coords);
     void getAtomCoords (PmAtomName name, vector<PmVector3> &coords);

     void getMassProps (PmMassProperties& props);

     void compPrincipalComponents(PmAtomFilter& filter, PmPcaResults& pca);

    protected:
      int number_of_atoms; 
      vector<PmAtom> atoms;
      PmExtent extent;

    private:
      int num_chains;
      void resize ();

   };

}

#endif

