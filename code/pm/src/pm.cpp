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
//* pm:             p r o t e i n   m o d e l e r              *
//*============================================================*

#include "pm/pm.h"

namespace ProteinMechanica {

////////////////////////////////////////////////////////////////
//              t i m e    i n t e r v a l                   //
//////////////////////////////////////////////////////////////

//*============================================================*
//*==========               hasTime                  ==========*
//*============================================================*
// test if given time is within an interval.

bool
PmTimeInterval::hasTime(float t)
  {
  if (!interval_set) {
    return true;
    }
  
  for (unsigned int i = 0; i < begin_intervals.size(); i++) {
    if ((t >= begin_intervals[i]) && (t <=  end_intervals[i])) {
      return true;
      }
    }

  return false;
  }

//*============================================================*
//*==========               addInterval              ==========*
//*============================================================*
// add a time interval.

void 
PmTimeInterval::addInterval(float begin, float end)
  {
  float tmp;

  if (begin > end) {
    tmp = end;
    end = begin;
    begin = tmp;
    }

  begin_intervals.push_back(begin);
  end_intervals.push_back(end);
  interval_set = true;
  }


}

