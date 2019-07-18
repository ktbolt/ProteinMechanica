
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
//* grid:                    g r i d                           *
//*============================================================*

#ifndef _GRID_PRV_PM_H_
#define _GRID_PRV_PM_H_

#include "pm/pm.h"
#include "graphics.h"
#include "pm/mth.h"
#include "grid.h"

namespace ProteinMechanica {

static bool iso_init = false;

#define Subset unsigned long

#define SubsetSize (sizeof(Subset)*8)

#define NilSubset NULL 

#define SubsetTest(block, mask) \
            (((block) == NilSubset) || ((mask) & 0x1) == 1)

#define SubsetSet(block) (((block) == NilSubset) || (*(block) & 0x1) == 1)


typedef struct LocalIndex {
  int i, j, k;
  } LocalIndex;

typedef struct IsoCellTable {
   LocalIndex list[3];
   } IsoCellTable;

 typedef struct EdgePt {
   int node;
   float index;
   struct EdgePt *next_pt;
   } EdgePt;

 typedef struct EdgePtBlockList {
   int size;
   EdgePt *ptr_block;
   struct EdgePtBlockList *next_block;
   } EdgePtBlockList;

static IsoCellTable cell_table[8][8];


static int
edge_list[12][2] = {{0, 1},
                    {0, 3},
                    {0, 4},
                    {1, 2},
                    {1, 5},
                    {2, 3},
                    {2, 6},
                    {3, 7},
                    {4, 5},
                    {4, 7},
                    {5, 6},
                    {6, 7}};


//*------------------------------------------------------------*
//* tables and such for extracting isosurfaces from cells.     *
//*------------------------------------------------------------*

#define  NUM_ISO_CASES  256

int *iso_edge_table[NUM_ISO_CASES];
int *iso_poly_table[NUM_ISO_CASES];

int case_0_edges[] = {0};
int case_0_polys[] = {0};

int case_1_edges[] = {3, 0, 1, 2};
int case_1_polys[] = {1, 3, 0, 1, 2};

int case_2_edges[] = {3, 0, 3, 4};
int case_2_polys[] = {1, 3, 0, 2, 1}; 

int case_3_edges[] = {4, 1, 2, 3, 4}; 
int case_3_polys[] = {2, 3, 1, 3, 2, 3, 0, 1, 2};

int case_4_edges[] = {3, 3, 5, 6}; 
int case_4_polys[] = {1, 3, 0, 2, 1}; 

int case_5_edges[] = {6, 0, 1, 2, 3, 5, 6};
int case_5_polys[] = {2, 3, 0, 1, 2, 3, 3, 5, 4}; 

int case_6_edges[] = {4, 0, 4, 5, 6}; 
int case_6_polys[] = {2, 3, 0, 1, 2, 3, 1, 3, 2}; 

int case_7_edges[] = {5, 1, 2, 4, 5, 6}; 
int case_7_polys[] = {3, 3, 0, 4, 3, 3, 0, 2, 4, 3, 0, 1, 2}; 

int case_8_edges[] = {3, 1, 5, 7}; 
int case_8_polys[] = {1, 3, 0, 1, 2}; 

int case_9_edges[] = {4, 0, 2, 5, 7}; 
int case_9_polys[] = {2, 3, 0, 2, 1, 3, 1, 2, 3}; 

int case_10_edges[] = {6, 0, 1, 3, 4, 5, 7}; 
int case_10_polys[] = {2, 3, 0, 3, 2, 3, 1, 4, 5}; 

int case_11_edges[] = {5, 2, 3, 4, 5, 7}; 
int case_11_polys[] = {3, 3, 0, 2, 4, 3, 1, 4, 2, 3, 1, 3, 4}; 

int case_12_edges[] = {4, 1, 3, 6, 7}; 
int case_12_polys[] = {2, 3, 0, 1, 2, 3, 0, 2, 3}; 

int case_13_edges[] = {5, 0, 2, 3, 6, 7}; 
int case_13_polys[] = {3, 3, 0, 2, 3, 3, 0, 3, 1, 3, 1, 3, 4}; 

int case_14_edges[] = {5, 7,  1,  6,  0,  4};
int case_14_polys[] = {3, 3, 0, 4, 2, 3, 0, 1, 4, 3, 1, 3, 4};

int case_15_edges[] = {4, 7,  6,  2,  4}; 
int case_15_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3};

int case_16_edges[] = {3, 2, 8, 9};
int case_16_polys[] = {1, 3, 0, 2, 1};

int case_17_edges[] = {4, 0, 1, 8, 9}; 
int case_17_polys[] = {2, 3, 0, 1, 3, 3, 0, 3, 2}; 

int case_18_edges[] = {6, 0, 3, 4, 2, 8, 9}; 
int case_18_polys[] = {2, 3, 0, 2, 1, 3, 3, 5, 4}; 

int case_19_edges[] = {5, 9, 8, 1, 4, 3};
int case_19_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_20_edges[] = {6, 2, 3, 5, 6, 8, 9}; 
int case_20_polys[] = {2, 3, 1, 3, 2, 3, 0, 5, 4}; 

int case_21_edges[] = {7, 0, 1, 3, 5, 6, 8, 9}; 
int case_21_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_22_edges[] = {7,  0,  4,  2,  8,  9,  5,  6}; 
int case_22_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_23_edges[] = {6, 1, 4, 5, 6, 8, 9}; 
int case_23_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_24_edges[] = {6, 1, 5, 7, 2, 8, 9}; 
int case_24_polys[] = {2, 3, 0, 1, 2, 3, 3, 5, 4}; 

int case_25_edges[] = {5, 8, 9, 7, 0, 5}; 
int case_25_polys[] = {3, 3, 0, 3, 4, 3, 0, 4, 2, 3, 0, 2, 1}; 

int case_26_edges[] = {9, 7,  9,  5,  1,  3,  8,  4,  2,  0}; 
int case_26_polys[] = {3, 3, 0, 3, 2, 3, 4, 8, 6, 3, 1, 5, 7}; 

int case_27_edges[] = {6, 3, 4, 5, 7, 8, 9}; 
int case_27_polys[] = {4, 3, 2, 3, 5, 3, 2, 5, 4, 3, 0, 4, 1, 3, 0, 2, 4}; 

int case_28_edges[] = {7, 3,  2,  1,  8,  9,  6,  7};
int case_28_polys[] = {3, 3, 1, 4, 3, 3, 0, 6, 2, 3, 0, 5, 6}; 

int case_29_edges[] = {6, 7,  0,  3,  6,  9,  8};
int case_29_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_30_edges[] = {8, 7,  1,  6,  9,  0,  2,  4,  8}; 
int case_30_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_31_edges[] = {5, 4,  8,  6,  9,  7};
int case_31_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_32_edges[] = {3, 4, 8, 10}; 
int case_32_polys[] = {1, 3, 0, 1, 2}; 

int case_33_edges[] = {6, 0, 1, 2, 4, 8, 10}; 
int case_33_polys[] = {2, 3, 0, 1, 2, 3, 3, 4, 5}; 

int case_34_edges[] = {4, 0, 3, 8, 10}; 
int case_34_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3}; 

int case_35_edges[] = {5, 1, 2, 3, 8, 10};
int case_35_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_36_edges[] = {6, 3, 5, 6, 4, 8, 10};
int case_36_polys[] = {2, 3, 0, 2, 1, 3, 3, 4, 5}; 

int case_37_edges[] = {9, 10,  6,  8,  4,  2,  5,  1,  3,  0}; 
int case_37_polys[] = {3, 3, 0, 3, 2, 3, 4, 8, 6, 3, 1, 5, 7}; 

int case_38_edges[] = {5, 8,  10,  0,  6,  5}; 
int case_38_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_39_edges[] = {6, 10,  8,  6,  5,  2,  1}; 
int case_39_polys[] = {4, 3, 2, 3, 5, 3, 2, 5, 4, 3, 0, 4, 1, 3, 0, 2, 4}; 

int case_40_edges[] = {6, 1, 5, 7, 4, 8, 10}; 
int case_40_polys[] = {2, 3, 0, 1, 2, 3, 3, 4, 5}; 

int case_41_edges[] = {7, 2,  0,  8,  4,  10,  7,  5}; 
int case_41_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_42_edges[] = {7, 8,  1,  0,  7,  5,  10,  3}; 
int case_42_polys[] = {3, 3, 1, 4, 3, 3, 0, 6, 2, 3, 0, 5, 6}; 

int case_43_edges[] = {6, 2, 3, 5, 7, 8, 10}; 
int case_43_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_44_edges[] = {7, 7,  10,  6,  8,  4,  1,  3}; 
int case_44_polys[] = {3, 3, 1, 4, 3, 3, 0, 6, 2, 3, 0, 5, 6}; 

int case_45_edges[] = {8, 2, 0, 7, 8, 3, 4, 6, 10};
int case_45_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_46_edges[] = {6, 0,  6,  1,  7,  10,  8}; 
int case_46_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_47_edges[] = {5, 6,  10,  7,  8,  2};
int case_47_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_48_edges[] = {4, 2, 9, 4, 10}; 
int case_48_polys[] = {2, 3, 0, 1, 3, 3, 0, 3, 2}; 

int case_49_edges[] = {5,  10,  4,  9,  0,  1}; 
int case_49_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_50_edges[] = {5, 3, 0, 10, 2, 9}; 
int case_50_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_51_edges[] = {4, 1, 3, 9, 10}; 
int case_51_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3}; 

int case_52_edges[] = {7, 2, 3, 4, 5, 6, 9, 10}; 
int case_52_polys[] = {3, 3, 1, 4, 3, 3, 0, 6, 2, 3, 0, 5, 6}; 

int case_53_edges[] = {8, 10,  4,  9,  6,  0,  3,  1,  5}; 
int case_53_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_54_edges[] = {6, 0,  10,  6,  5,  2,  9}; 
int case_54_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_55_edges[] = {5, 10, 6, 9, 5, 1};
int case_55_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_56_edges[] = {7, 2,  9,  1,  7,  5,  4,  10}; 
int case_56_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_57_edges[] = {6, 9, 0, 7, 5, 4, 10}; 
int case_57_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_58_edges[] = {8, 3,  0,  10,  5,  2,  1,  9,  7}; 
int case_58_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_59_edges[] = {5, 3, 5, 7, 10, 9};
int case_59_polys[] = {3, 3, 0, 4, 3, 3, 0, 2, 4, 3, 0, 1, 2};

int case_60_edges[] = {8,  3,  4,  6,  1,  7,  2,  10,  9}; 
int case_60_polys[] = {4, 3, 0, 2, 4, 3, 0, 4, 3, 3, 1, 5, 7, 3, 1, 7, 6}; 

int case_61_edges[] = {7, 9,  4,  10,  0,  3,  7,  6};
int case_61_polys[] = {3, 3, 1, 3, 4, 3, 0, 2, 6, 3, 0, 6, 5};

int case_62_edges[] = {7, 9,  7,  2,  1,  0,  10,  6};
int case_62_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_63_edges[] = {4, 6, 10, 7, 9};
int case_63_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3};

int case_64_edges[] = {3, 6, 10, 11}; 
int case_64_polys[] = {1, 3, 0, 1, 2}; 

int case_65_edges[] = {6, 0, 1, 2, 6, 10, 11}; 
int case_65_polys[] = {2, 3, 0, 1, 2, 3, 3, 4, 5};

int case_66_edges[] = {6, 0, 3, 4, 6, 10, 11}; 
int case_66_polys[] = {2, 3, 0, 2, 1, 3, 3, 4, 5}; 

int case_67_edges[] = {7, 1,  6,  3,  11,  10,  2,  4}; 
int case_67_polys[] = {3, 3, 1, 4, 3, 3, 0, 6, 2, 3, 0, 5, 6}; 

int case_68_edges[] = {4, 3, 5, 10, 11}; 
int case_68_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3}; 

int case_69_edges[] = {7, 10,  0,  3,  2,  1,  11,  5}; 
int case_69_polys[] = {3, 3, 1, 4, 3, 3, 0, 6, 2, 3, 0, 5, 6}; 

int case_70_edges[] = {5, 0,  4,  5,  10,  11};
int case_70_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_71_edges[] = {6, 4,  5,  1,  2,  10,  11}; 
int case_71_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_72_edges[] = {6, 1, 5, 7, 6, 10, 11}; 
int case_72_polys[] = {2, 3, 0, 1, 2, 3, 3, 4, 5}; 

int case_73_edges[] = {7, 5,  7,  6,  11,  10,  0,  2}; 
int case_73_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_74_edges[] = {9, 1,  0,  7,  5,  11,  4,  10,  3,  6};
int case_74_polys[] = {3, 3, 0, 3, 2, 3, 4, 8, 6, 3, 1, 5, 7}; 

int case_75_edges[] = {8, 4,  3,  2,  10,  5,  6,  7,  11}; 
int case_75_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_76_edges[] = {5, 1, 7, 3, 11, 10}; 
int case_76_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1}; 

int case_77_edges[] = {6, 7,  3,  11,  10,  0,  2}; 
int case_77_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_78_edges[] = {6, 4,  0,  10,  11,  1,  7}; 
int case_78_polys[] = {4, 3, 0, 2, 3, 3, 0, 3, 1, 3, 1, 3, 5, 3, 1, 5, 4};

int case_79_edges[] = {5, 4, 10, 11, 2, 7};
int case_79_polys[] = {3, 3, 0, 4, 3, 3, 0, 2, 4, 3, 0, 1, 2};

int case_80_edges[] = {6, 2, 8, 9, 6, 10, 11}; 
int case_80_polys[] = {2, 3, 0, 2, 1, 3, 3, 4, 5}; 

int case_81_edges[] = {7, 0,  10,  8,  6,  11,  1,  9};
int case_81_polys[] = {3, 3, 1, 4, 3, 3, 0, 6, 2, 3, 0, 5, 6}; 

int case_82_edges[] = {9, 6,  3,  11,  10,  9,  0,  2,  4,  8}; 
int case_82_polys[] = {3, 3, 0, 3, 2, 3, 4, 8, 6, 3, 1, 5, 7}; 

int case_83_edges[] = {8, 9,  8,  1,  11,  4,  10,  3,  6}; 
int case_83_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_84_edges[] = {7, 10,  11,  8,  9,  2,  3,  5}; 
int case_84_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_85_edges[] = {8, 5,  1,  3,  11,  10,  9,  0,  8}; 
int case_85_polys[] = {4, 3, 0, 2, 4, 3, 0, 4, 3, 3, 1, 5, 7, 3, 1, 7, 6}; 

int case_86_edges[] = {8, 0,  4,  5,  2,  10,  8,  11,  9}; 
int case_86_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_87_edges[] = {7, 1,  8,  9,  4,  10,  5,  11};
int case_87_polys[] = {3, 3, 1, 3, 4, 3, 0, 2, 6, 3, 0, 6, 5};

int case_88_edges[] = {9, 10,  8,  6,  11,  5,  2,  1,  9,  7}; 
int case_88_polys[] = {3, 3, 0, 3, 2, 3, 4, 8, 6, 3, 1, 5, 7}; 

int case_89_edges[] = {8, 5,  7,  0,  6,  9,  11,  8,  10}; 
int case_89_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_90_edges[] = {12, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; 
int case_90_polys[] = {4, 3, 0, 4, 3, 3, 6, 10, 11, 3, 2, 9, 8, 3, 1, 5, 7}; 

int case_91_edges[] = {9, 3, 4, 5, 6, 7, 8,9, 10, 11};
int case_91_polys[] = {3, 3, 0, 2, 3, 3, 4, 6, 8, 3, 1, 7, 5};

int case_92_edges[] = {8, 10,  11,  3,  8,  7,  9,  1,  2}; 
int case_92_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_93_edges[] = {7, 8,  10,  9,  11,  7,  0,  3};
int case_93_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_94_edges[] = {9, 4, 0, 10, 8, 11, 1, 7, 2, 9};
int case_94_polys[] = {3, 3, 0, 2, 3, 3, 4, 6, 8, 3, 1, 7, 5};

int case_95_edges[] = {6, 4, 8, 10, 7, 9, 11};
int case_95_polys[] = {2, 3, 0, 2, 1, 3, 3, 4, 5};

int case_96_edges[] = {4, 4,  6,  8,  11}; 
int case_96_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3};

int case_97_edges[] = {7, 4,  8,  0,  2,  1,  6,  11}; 
int case_97_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_98_edges[] = {5, 0, 3, 6, 8, 11}; 
int case_98_polys[] = {3, 3, 0, 3, 4, 3, 0, 4, 2, 3, 0, 2, 1}; 

int case_99_edges[] = {6, 3,  8,  6,  11,  2,  1}; 
int case_99_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_100_edges[] = {5, 5,  3,  11,  4,  8}; 
int case_100_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_101_edges[] = {8, 5,  3,  11,  1,  4,  0,  8,  2}; 
int case_101_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_102_edges[] = {4, 0,  5,  8,  11}; 
int case_102_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3}; 

int case_103_edges[] = {5, 8,  2,  11,  1,  5};
int case_103_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_104_edges[] = {7, 11,  6,  7,  5,  1,  8,  4}; 
int case_104_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_105_edges[] = {8, 2,  8,  0,  7,  5,  11,  4,  6}; 
int case_105_polys[] = {4, 3, 0, 2, 4, 3, 0, 4, 3, 3, 1, 5, 7, 3, 1, 7, 6}; 

int case_106_edges[] = {8, 11,  6,  8,  7,  3,  5,  0,  1}; 
int case_106_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_107_edges[] = {7, 7,  11,  5,  6,  3,  2,  8};
int case_107_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_108_edges[] = {6, 1, 3, 4, 7, 8, 11}; 
int case_108_polys[] = {4, 3, 1, 2, 5, 3, 2, 4, 5,  3, 0, 1, 5, 3, 0, 5, 3};

int case_109_edges[] = {7, 8,  2,  4,  0,  3,  11,  7};
int case_109_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_110_edges[] = {5, 0,  1,  8,  7,  11};
int case_110_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_111_edges[] = {4, 7,  2,  11,  8};
int case_111_polys[] = {2, 3, 0, 1, 3, 3, 0, 3, 2};

int case_112_edges[] = {5, 2,  9,  4,  11,  6}; 
int case_112_polys[] = {3, 3, 0, 4, 2, 3, 0, 1, 4, 3, 1, 3, 4};

int case_113_edges[] = {6, 4,  9,  11,  6,  0,  1}; 
int case_113_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_114_edges[] = {6, 9,  2,  11,  6,  0,  3}; 
int case_114_polys[] = {4, 3, 2, 3, 5, 3, 2, 5, 4, 3, 0, 4, 1, 3, 0, 2, 4}; 

int case_115_edges[] = {5, 9, 11, 1, 6, 3};
int case_115_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_116_edges[] = {6, 4,  11,  3,  5,  9,  2}; 
int case_116_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_117_edges[] = {7, 1,  5,  0,  3,  4,  9,  11};
int case_117_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_118_edges[] = {5, 11,  9,  5,  2,  0};
int case_118_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_119_edges[] = {4, 1, 5, 9, 11};
int case_119_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3};

int case_120_edges[] = {8, 2,  9,  4,  1,  11,  7,  6,  5};
int case_120_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_121_edges[] = {7, 6,  5,  11,  7,  9,  4,  0};
int case_121_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_122_edges[] = {9, 9,  2,  11,  7,  6,  0,  3,  1,  5};
int case_122_polys[] = {3, 3, 0, 2, 3, 3, 4, 6, 8, 3, 1, 7, 5};

int case_123_edges[] = {6, 3, 5, 6, 7, 9, 11};
int case_123_polys[] = {2, 3, 0, 1, 2, 3, 3, 4, 5};

int case_124_edges[] = {7, 4,  9,  2,  11,  7,  3,  1};
int case_124_polys[] = {3, 3, 1, 3, 4, 3, 0, 2, 6, 3, 0, 6, 5};

int case_125_edges[] = {6, 0, 3, 4, 7, 9, 11};
int case_125_polys[] = {2, 3, 0, 1, 2, 3, 3, 4, 5};

int case_126_edges[] = {6, 0, 1, 2, 7, 9, 11};
int case_126_polys[] = {2, 3, 0, 2, 1, 3, 3, 4, 5};

int case_127_edges[] = {3, 7, 9, 11};
int case_127_polys[] = {1, 3, 0, 1, 2};

int case_128_edges[] = {3, 7, 9, 11}; 
int case_128_polys[] = {1, 3, 0, 2, 1}; 

int case_129_edges[] = {6, 0, 1, 2, 7, 9, 11}; 
int case_129_polys[] = {2, 3, 0, 1, 2, 3, 3, 5, 4}; 

int case_130_edges[] = {6, 0, 3, 4, 7, 9, 11}; 
int case_130_polys[] = {2, 3, 0, 2, 1, 3, 3, 5, 4};

int case_131_edges[] = {7, 4,  9,  2,  11,  7,  3,  1}; 
int case_131_polys[] = {3, 3, 1, 4, 3, 3, 0, 6, 2, 3, 0, 5, 6}; 

int case_132_edges[] = {6, 3, 5, 6, 7, 9, 11}; 
int case_132_polys[] = {2, 3, 0, 2, 1, 3, 3, 5, 4}; 

int case_133_edges[] = {9, 9,  2,  11,  7,  6,  0,  3,  1,  5}; 
int case_133_polys[] = {3, 3, 0, 3, 2, 3, 4, 8, 6, 3, 1, 5, 7}; 

int case_134_edges[] = {7, 6,  5,  11,  7,  9,  4,  0}; 
int case_134_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_135_edges[] = {8, 6,  5,  4,  11,  1,  7,  2,  9}; 
int case_135_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_136_edges[] = {4, 1, 5, 9, 11}; 
int case_136_polys[] = {2, 3, 0, 1, 3, 3, 0, 3, 2}; 

int case_137_edges[] = {5, 11,  9,  5,  2,  0}; 
int case_137_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_138_edges[] = {7, 1,  5,  0,  3,  4,  9,  11};
int case_138_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_139_edges[] = {6, 5,  2,  3,  4,  9,  11};
int case_139_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_140_edges[] = {5, 9, 11, 1, 6, 3}; 
int case_140_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1}; 

int case_141_edges[] = {6, 6,  3,  11,  9,  0,  2};
int case_141_polys[] = {4, 3, 2, 3, 5, 3, 2, 5, 4, 3, 0, 4, 1, 3, 0, 2, 4}; 

int case_142_edges[] = {6, 1,  6,  11,  9,  0,  4}; 
int case_142_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_143_edges[] = {5, 2,  9,  4,  11,  6};
int case_143_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 1, 3, 1, 4, 3};

int case_144_edges[] = {4, 7,  2,  11,  8}; 
int case_144_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3};

int case_145_edges[] = {5, 0,  1,  8,  7,  11}; 
int case_145_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_146_edges[] = {7, 8,  2,  4,  0,  3,  11,  7}; 
int case_146_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_147_edges[] = {6, 8,  1,  7,  11,  4,  3}; 
int case_147_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_148_edges[] = {7, 7,  11,  5,  6,  3,  2,  8}; 
int case_148_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_149_edges[] = {8, 0,  1, 3,  5,  6, 7,  8,  11};
int case_149_polys[] = {4, 3, 2, 4, 3, 3, 1, 5, 7, 3, 0, 1, 7, 3, 0, 7, 6};

int case_150_edges[] = {8, 0, 2, 4, 5, 6, 7, 8, 11}; 
int case_150_polys[] = {4, 3, 0, 2, 4, 3, 0, 4, 3, 3, 1, 5, 7, 3, 1, 7, 6}; 

int case_151_edges[] = {7, 11,  6,  7,  5,  1,  8,  4};
int case_151_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_152_edges[] = {5, 8,  2,  11,  1,  5}; 
int case_152_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_153_edges[] = {4, 5,  0,  11,  8}; 
int case_153_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3}; 

int case_154_edges[] = {8, 8,  2,  11,  4,  1,  0,  5,  3}; 
int case_154_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_155_edges[] = {5, 5,  3,  11,  4,  8};
int case_155_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_156_edges[] = {6, 1,  11,  2,  8,  6,  3}; 
int case_156_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_157_edges[] = {5, 0, 3, 6, 8, 11};
int case_157_polys[] = {3, 3, 0, 4, 3, 3, 0, 2, 4, 3, 0, 1, 2};

int case_158_edges[] = {7, 4,  8,  0,  2,  1,  6,  11};
int case_158_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_159_edges[] = {4, 4,  6,  8,  11};
int case_159_polys[] = {2, 3, 0, 1, 3, 3, 0, 3, 2};

int case_160_edges[] = {6, 4, 8, 10, 7, 9, 11}; 
int case_160_polys[] = {2, 3, 0, 1, 2, 3, 3, 5, 4}; 

int case_161_edges[] = {9, 4, 0, 10, 8, 11, 1, 7, 2, 9}; 
int case_161_polys[] = {3, 3, 0, 3, 2, 3, 4, 8, 6, 3, 1, 5, 7}; 

int case_162_edges[] = {7, 8,  10,  9,  11,  7,  0,  3}; 
int case_162_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_163_edges[] = {8, 1, 2, 3, 7, 8, 9, 10, 11}; 
int case_163_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_164_edges[] = {9, 3, 4, 5, 6, 7, 8,9, 10, 11}; 
int case_164_polys[] = {3, 3, 0, 3, 2, 3, 4, 8, 6, 3, 1, 5, 7}; 

int case_165_edges[] = {12, 2, 9, 8, 1, 0, 7, 5, 11, 4, 10, 3, 6}; 
int case_165_polys[] = {4, 3, 0, 4, 3, 3, 6, 10, 11, 3, 2, 9, 8, 3, 1, 5, 7}; 

int case_166_edges[] = {8, 8,  10,  0,  9,  6,  11,  5,  7}; 
int case_166_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 5, 7};

int case_167_edges[] = {9, 10,  8,  6,  11,  5,  2,  1,  9,  7};
int case_167_polys[] = {3, 3, 0, 2, 3, 3, 4, 6, 8, 3, 1, 7, 5};

int case_168_edges[] = {7, 1,  8,  9,  4,  10,  5,  11}; 
int case_168_polys[] = {3, 3, 1, 4, 3, 3, 0, 6, 2, 3, 0, 5, 6}; 

int case_169_edges[] = {8, 11,  9,  5,  10,  2,  8,  0,  4};
int case_169_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_170_edges[] = {8,  8,  9,  10,  0,  3,  1,  11,  5}; 
int case_170_polys[] = {4, 3, 0, 2, 4, 3, 0, 4, 3, 3, 1, 5, 7, 3, 1, 7, 6}; 

int case_171_edges[] = {7, 10,  11,  8,  9,  2,  3,  5};
int case_171_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_172_edges[] = {8, 3,  6,  1,  4,  11,  10,  9,  8}; 
int case_172_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_173_edges[] = {9, 6,  3,  11,  10,  9,  0,  2,  4,  8};
int case_173_polys[] = {3, 3, 0, 2, 3, 3, 4, 6, 8, 3, 1, 7, 5};

int case_174_edges[] = {7, 0,  10,  8,  6,  11,  1,  9};
int case_174_polys[] = {3, 3, 1, 3, 4, 3, 0, 2, 6, 3, 0, 6, 5};

int case_175_edges[] = {6, 2, 8, 9, 6, 10, 11};
int case_175_polys[] = {2, 3, 0, 1, 2, 3, 3, 5, 4};

int case_176_edges[] = {5, 4, 10, 11, 2, 7}; 
int case_176_polys[] = {3, 3, 0, 3, 4, 3, 0, 4, 2, 3, 0, 2, 1}; 

int case_177_edges[] = {6, 1,  0,  7,  11,  4,  10}; 
int case_177_polys[] = {4, 3, 2, 3, 5, 3, 2, 5, 4, 3, 0, 4, 1, 3, 0, 2, 4}; 

int case_178_edges[] = {6, 10,  2,  11,  7,  0,  3};
int case_178_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_179_edges[] = {5, 1, 7, 3, 11, 10};
int case_179_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_180_edges[] = {8, 7,  11,  2,  5,  10,  6,  4,  3}; 
int case_180_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_181_edges[] = {9, 1,  0,  7,  5,  11,  4,  10,  3,  6};
int case_181_polys[] = {3, 3, 0, 2, 3, 3, 4, 6, 8, 3, 1, 7, 5};

int case_182_edges[] = {7, 5,  7,  6,  11,  10,  0,  2};
int case_182_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_183_edges[] = {6, 1, 5, 7, 6, 10, 11};
int case_183_polys[] = {2, 3, 0, 2, 1, 3, 3, 5, 4};

int case_184_edges[] = {6, 11,  2,  1,  5,  10,  4};
int case_184_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_185_edges[] = {5, 0,  4,  5,  10,  11};
int case_185_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_186_edges[] = {7, 10,  0,  3,  2,  1,  11,  5};
int case_186_polys[] = {3, 3, 1, 3, 4, 3, 0, 2, 6, 3, 0, 6, 5};

int case_187_edges[] = {4, 3, 5, 10, 11};
int case_187_polys[] = {2, 3, 0, 1, 3, 3, 0, 3, 2};

int case_188_edges[] = {7, 1,  6,  3,  11,  10,  2,  4};
int case_188_polys[] = {3, 3, 1, 3, 4, 3, 0, 2, 6, 3, 0, 6, 5};

int case_189_edges[] = {6, 0, 3, 4, 6, 10, 11};
int case_189_polys[] = {2, 3, 0, 1, 2, 3, 3, 5, 4};

int case_190_edges[] = {6, 0, 1, 2, 6, 10, 11};
int case_190_polys[] = {2, 3, 0, 2, 1, 3, 3, 5, 4};

int case_191_edges[] = {3, 6, 10, 11};
int case_191_polys[] = {1, 3, 0, 2, 1};

int case_192_edges[] = {4, 6, 10, 7, 9}; 
int case_192_polys[] = {2, 3, 0, 1, 3, 3, 0, 3, 2}; 

int case_193_edges[] = {7, 9,  7,  2,  1,  0,  10,  6}; 
int case_193_polys[] = {3, 3, 2, 4, 3, 3, 0, 1, 6, 3, 0, 6, 5}; 

int case_194_edges[] = {7, 9,  4,  10,  0,  3,  7,  6}; 
int case_194_polys[] = {3, 3, 1, 4, 3, 3, 0, 6, 2, 3, 0, 5, 6}; 

int case_195_edges[] = {8,  1,  7,  2,  3,  4,  6,  9,  10};
int case_195_polys[] = {4, 3, 0, 2, 4, 3, 0, 4, 3, 3, 1, 5, 7, 3, 1, 7, 6}; 

int case_196_edges[] = {5, 3, 5, 7, 10, 9};
int case_196_polys[] = {3, 3, 0, 3, 4, 3, 0, 4, 2, 3, 0, 2, 1}; 

int case_197_edges[] = {8, 9,  7,  10,  2,  5,  1,  3,  0}; 
int case_197_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_198_edges[] = {6, 10,  5,  4,  0,  7,  9}; 
int case_198_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_199_edges[] = {7, 2,  9,  1,  7,  5,  4,  10};
int case_199_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_200_edges[] = {5, 10, 6, 9, 5, 1}; 
int case_200_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1}; 

int case_201_edges[] = {6, 9,  5,  6,  10,  2,  0};
int case_201_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_202_edges[] = {8, 1,  5,  9,  0,  6,  3,  10,  4};
int case_202_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_203_edges[] = {7, 2, 3, 4, 5, 6, 9, 10};
int case_203_polys[] = {3, 3, 1, 3, 4, 3, 0, 2, 6, 3, 0, 6, 5};

int case_204_edges[] = {4, 9,  10,  1,  3};
int case_204_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3}; 

int case_205_edges[] = {5, 3, 0, 10, 2, 9};
int case_205_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_206_edges[] = {5,  10,  4,  9,  0,  1};
int case_206_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_207_edges[] = {4, 2, 9, 4, 10};
int case_207_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3};

int case_208_edges[] = {5, 6,  10,  7,  8,  2}; 
int case_208_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_209_edges[] = {6, 7,  8,  1,  0,  10,  6}; 
int case_209_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_210_edges[] = {8, 6,  10,  7,  3,  8,  4,  2,  0}; 
int case_210_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_211_edges[] = {7, 7,  10,  6,  8,  4,  1,  3};
int case_211_polys[] = {3, 3, 1, 3, 4, 3, 0, 2, 6, 3, 0, 6, 5};

int case_212_edges[] = {6, 7,  10,  8,  2,  5,  3};
int case_212_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_213_edges[] = {7, 8,  1,  0,  7,  5,  10,  3};
int case_213_polys[] = {3, 3, 1, 3, 4, 3, 0, 2, 6, 3, 0, 6, 5};

int case_214_edges[] = {7, 2,  0,  8,  4,  10,  7,  5};
int case_214_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_215_edges[] = {6, 1, 5, 7, 4, 8, 10};
int case_215_polys[] = {2, 3, 0, 2, 1, 3, 3, 5, 4};

int case_216_edges[] = {6, 10,  6,  8,  2,  5,  1}; 
int case_216_polys[] = {4, 3, 2, 3, 5, 3, 2, 5, 4, 3, 0, 4, 1, 3, 0, 2, 4}; 

int case_217_edges[] = {5, 8,  10,  0,  6,  5};
int case_217_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_218_edges[] = {9, 10,  6,  8,  4,  2,  5,  1,  3,  0};
int case_218_polys[] = {3, 3, 0, 2, 3, 3, 4, 6, 8, 3, 1, 7, 5};

int case_219_edges[] = {6, 3, 5, 6, 4, 8, 10};
int case_219_polys[] = {2, 3, 0, 1, 2, 3, 3, 5, 4};

int case_220_edges[] = {5, 1, 2, 3, 8, 10};
int case_220_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_221_edges[] = {4, 0, 3, 8, 10};
int case_221_polys[] = {2, 3, 0, 1, 3, 3, 0, 3, 2};

int case_222_edges[] = {6, 0, 1, 2, 4, 8, 10};
int case_222_polys[] = {2, 3, 0, 2, 1, 3, 3, 5, 4};

int case_223_edges[] = {3, 4, 8, 10};
int case_223_polys[] = {1, 3, 0, 2, 1};

int case_224_edges[] = {5, 4,  8,  6,  9,  7}; 
int case_224_polys[] = {3, 3, 0, 4, 2, 3, 0, 3, 4, 3, 0, 1, 3};

int case_225_edges[] = {8, 4,  8,  6,  0,  9,  2,  7,  1}; 
int case_225_polys[] = {4, 3, 0, 6, 2, 3, 0, 1, 6, 3, 1, 4, 6, 3, 3, 7, 5};

int case_226_edges[] = {6, 6,  8,  9,  7,  3,  0}; 
int case_226_polys[] = {4, 3, 1, 2, 3, 3, 1, 3, 4, 3, 0, 4, 3, 3, 1, 4, 5}; 

int case_227_edges[] = {7, 3,  2,  1,  8,  9,  6,  7};
int case_227_polys[] = {3, 3, 1, 3, 4, 3, 0, 2, 6, 3, 0, 6, 5};

int case_228_edges[] = {6, 8,  4,  9,  7,  3,  5}; 
int case_228_polys[] = {4, 3, 0, 2, 3, 3, 0, 3, 1, 3, 1, 3, 5, 3, 1, 5, 4};

int case_229_edges[] = {9, 7,  9,  5,  1,  3,  8,  4,  2,  0};
int case_229_polys[] = {3, 3, 0, 2, 3, 3, 4, 6, 8, 3, 1, 7, 5};

int case_230_edges[] = {5, 8, 9, 7, 0, 5};
int case_230_polys[] = {3, 3, 0, 4, 3, 3, 0, 2, 4, 3, 0, 1, 2};

int case_231_edges[] = {6, 1, 5, 7, 2, 8, 9};
int case_231_polys[] = {2, 3, 0, 2, 1, 3, 3, 4, 5};

int case_232_edges[] = {6, 9,  6,  8,  4,  5,  1};
int case_232_polys[] = {4, 3, 1, 3, 4, 3, 0, 4, 3, 3, 0, 3, 2, 3, 0, 5, 4}; 

int case_233_edges[] = {7,  0,  4,  2,  8,  9,  5,  6};
int case_233_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_234_edges[] = {7, 0, 1, 3, 5, 6, 8, 9};
int case_234_polys[] = {3, 3, 2, 3, 4, 3, 0, 6, 1, 3, 0, 5, 6};

int case_235_edges[] = {6, 2, 3, 5, 6, 8, 9};
int case_235_polys[] = {2, 3, 1, 2, 3, 3, 0, 4, 5};

int case_236_edges[] = {5, 9, 8, 1, 4, 3};
int case_236_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 3, 3, 0, 3, 1};

int case_237_edges[] = {6, 0, 3, 4, 2, 8, 9};
int case_237_polys[] = {2, 3, 0, 1, 2, 3, 3, 4, 5};

int case_238_edges[] = {4, 0, 1, 8, 9};
int case_238_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3};

int case_239_edges[] = {3, 2, 8, 9};
int case_239_polys[] = {1, 3, 0, 1, 2};

int case_240_edges[] = {4, 2,  4,  7,  6};
int case_240_polys[] = {2, 3, 0, 3, 1, 3, 0, 2, 3}; 

int case_241_edges[] = {5, 7,  1,  6,  0,  4};
int case_241_polys[] = {3, 3, 0, 2, 4, 3, 0, 4, 1, 3, 1, 4, 3};

int case_242_edges[] = {5, 0, 2, 3, 6, 7};
int case_242_polys[] = {3, 3, 0, 3, 2, 3, 0, 1, 3, 3, 1, 4, 3};

int case_243_edges[] = {4, 1, 3, 6, 7};
int case_243_polys[] = {2, 3, 0, 2, 1, 3, 0, 3, 2};

int case_244_edges[] = {5, 2, 3, 4, 5, 7};
int case_244_polys[] = {3, 3, 0, 4, 2, 3, 1, 2, 4, 3, 1, 4, 3};

int case_245_edges[] = {6, 0, 1, 3, 4, 5, 7};
int case_245_polys[] = {2, 3, 0, 2, 3, 3, 1, 5, 4};

int case_246_edges[] = {4, 0, 2, 5, 7};
int case_246_polys[] = {2, 3, 0, 1, 2, 3, 1, 3, 2};

int case_247_edges[] = {3, 1, 5, 7};
int case_247_polys[] = {1, 3, 0, 2, 1};

int case_248_edges[] = {5, 1, 2, 4, 5, 6};
int case_248_polys[] = {3, 3, 0, 3, 4, 3, 0, 4, 2, 3, 0, 2, 1};

int case_249_edges[] = {4, 0, 4, 5, 6};
int case_249_polys[] = {2, 3, 0, 2, 1, 3, 1, 2, 3};

int case_250_edges[] = {6, 0, 1, 2, 3, 5, 6};
int case_250_polys[] = {2, 3, 0, 2, 1, 3, 3, 4, 5};

int case_251_edges[] = {3, 3, 5, 6};
int case_251_polys[] = {1, 3, 0, 1, 2};

int case_252_edges[] = {4, 1, 2, 3, 4};
int case_252_polys[] = {2, 3, 1, 2, 3, 3, 0, 2, 1};

int case_253_edges[] = {3, 0, 3, 4};
int case_253_polys[] = {1, 3, 0, 1, 2};

int case_254_edges[] = {3, 0, 1, 2};
int case_254_polys[] = {1, 3, 0, 2, 1};



static void
grid_IsoInit();

static void
grid_IsoExtract (int dims[3], float cell_sizes[3], float *data, float level,
                 vector<PmVector3>& verts, vector<int>& polys);


static void
grid_IsoSeededExtract (int dims[3], float cell_sizes[3], float *data, float level,
                       int iseed[3], vector<PmVector3>& verts, vector<int>& polys);

static void
grid_SliceExtract (int dims[3], float cell_sizes[3], float *data, int dim, int index,
                   vector<PmVector3>& verts, vector<int>& polys, vector<float>& vals);

static bool
grid_SubsetBitIsSet (int num_blocks, Subset *blocks, int n);

static void
grid_SubsetBitSet (int num_blocks, Subset *blocks, int n);

static void
grid_SubsetCreate (int size, int *p_num_blocks, Subset **p_blocks);

static void
grid_MapData(float vmin, float vmax, vector<float>& data, vector<PmVector3>& colors);

static void
grid_ColorMapSpectrumCreate (int num_colors, vector<PmVector3>& colors);

static void
grid_ColorMapHsvToRgb (float h, float s, float v, float *r, float *g, float *b);



  }

#endif



