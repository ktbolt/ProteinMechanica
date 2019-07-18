
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
 * atom:                 a t o m                              *
 *============================================================*/

#ifndef _ATOM_PRV_PM_H_
#define _ATOM_PRV_PM_H_

namespace ProteinMechanica {

//===== atom colors =====//

static float
atom_colors[][3] = {
  {0.0, 0.0, 0.0},
  {0.7, 0.7, 0.7},    /*  carbon     */
  {1.0, 1.0, 1.0},    /*  hydrogen   */
  {0.6, 0.6, 1.0},    /*  nitrogen   */
  {0.9, 0.0, 0.0},    /*  oxygen     */
  {1.0, 0.8, 0.2},    /*  sulfur     */
  {1.0, 0.6, 0.0},    /*  phosphor   */
  {1.0, 0.6, 0.0},    /*  iron       */
  {1.0, 0.8, 0.2}};   /*  selinium   */


//===== atom properties =====//
// radii in nm,   mass in amu 
// these are very approximate VanderWaals radii and
// are used for drawing atoms as balls or for calculating 
// atomic overlap.

static const float Crad = 0.15; 
static const float Hrad = 0.04; 
static const float Nrad = 0.11; 
static const float Orad = 0.105; 
static const float Prad = 0.16; 
static const float Srad = 0.16; 
static const float Irad = 0.12; 

static PmAtomicProps 
atomic_props[] = {{0.00,  0.0,   0.000},            /*  undefined  */
                  {Crad,  0.0,  12.011},            /*  carbon     */
                  {Hrad,  0.0,   1.00794},          /*  hydrogen   */
                  {Nrad,  0.0,  14.0067},           /*  nitrogen   */
                  {Orad,  0.0,  15.9994},           /*  oxygen     */
                  {Srad,  0.0,  32.066},            /*  sulfur     */
                  {Prad,  0.0,  30.9738},           /*  phosphorus */
                  {Irad,  0.0,  55.847},            /*  iron       */
                  {0.115,  0.0,  78.96 }};           /*  selenium   */
}

#endif


