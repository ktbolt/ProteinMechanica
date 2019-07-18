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

#ifndef _BEXP_PM_H_
#define _BEXP_PM_H_

namespace ProteinMechanica {

//-------------------------------------------------------------
//                    PmBooleanRelOpType                      |
//-------------------------------------------------------------

typedef enum PmBooleanRelOpType {
  PM_BOOLEAN_RELOP_UNKNOWN,
  PM_BOOLEAN_RELOP_LT,
  PM_BOOLEAN_RELOP_LE,
  PM_BOOLEAN_RELOP_GT,
  PM_BOOLEAN_RELOP_GE,
  PM_BOOLEAN_RELOP_EQ,
  PM_BOOLEAN_RELOP_NE
  } PmBooleanRelOpType;

//-------------------------------------------------------------
//                    PmBooleanOpType                         |
//-------------------------------------------------------------

typedef enum PmBooleanOpType {
  PM_BOOLEAN_OP_UNKNOWN,
  PM_BOOLEAN_OP_AND,
  PM_BOOLEAN_OP_OR,
  PM_BOOLEAN_OP_NOT
  } PmBooleanOpType;

//-------------------------------------------------------------
//                    PmBooleanCondition                      |
//-------------------------------------------------------------

typedef struct PmBooleanCondition {
  string var1, var2;
  PmBooleanRelOpType op;
  } PmBooleanCondition;

//-------------------------------------------------------------
//                    PmBooleanExpression                     |
//-------------------------------------------------------------

class PmBooleanExpression {
  public:
     string name;
     vector<PmBooleanOpType> ops;
     vector<PmBooleanCondition> conditions;
     bool current_value;
     bool evaluate();
  };
}

#endif



