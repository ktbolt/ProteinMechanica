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
 * bexp:        b o o l e a n   e x p r e s s i o n           *
 *============================================================*/

#include "pm/pm.h" 
#include "bexp.h" 

namespace ProteinMechanica {

  bool convToFloat (string s, float& val);

//*============================================================*
//*==========      PmBooleanExpression::evaluate     ==========*
//*============================================================*

bool 
PmBooleanExpression::evaluate()
  {
  //fprintf (stderr, ">>>>>> PmBooleanExpression::evaluate \n"); 
  PmBooleanCondition cond;
  string var1, var2;
  float val1, val2;
  vector<bool> cond_vals;
  bool bval;

  for (unsigned int i = 0; i < conditions.size(); i++) {
    var1 = conditions[i].var1;
    var2 = conditions[i].var2;
    //fprintf (stderr, ">>> cond %d var1=%s var2=%s \n", i, var1.c_str(), var2.c_str()); 

    if (isalpha(var1[0])) {
      val1 = pmSystem.getVariableValue(var1);
      }
    else {
      convToFloat(var1, val1);
      }

    if (isalpha(var2[0])) {
      val2 = pmSystem.getVariableValue(var2);
      }
    else {
      convToFloat(var2, val2);
      }

    switch (conditions[i].op) {
      case PM_BOOLEAN_RELOP_LT:
        bval = val1 < val2;
      break;

      case PM_BOOLEAN_RELOP_LE:
        bval = val1 <= val2;
      break;

      case PM_BOOLEAN_RELOP_GT:
        bval = val1 > val2;
      break;

      case PM_BOOLEAN_RELOP_GE:
        bval = val1 >= val2;
      break;

      case PM_BOOLEAN_RELOP_EQ:
        bval = val1 == val2;
      break;

      case PM_BOOLEAN_RELOP_NE:
        bval = val1 != val2;
      break;
      }

    cond_vals.push_back(bval);
    }

  current_value = cond_vals[0];

  for (unsigned int i = 1; i < conditions.size(); i++) {
    if (ops[i-1] == PM_BOOLEAN_OP_AND) {
      current_value = current_value && cond_vals[i];
      }
    else if (ops[i-1] == PM_BOOLEAN_OP_OR) {
      current_value = current_value || cond_vals[i];
      }
    }

  return current_value;
  }

}

