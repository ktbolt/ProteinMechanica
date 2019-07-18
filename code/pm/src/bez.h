
/*------------------------------------------------------------*
 *                                                            *
 *                    ****  bez  ****                         *
 *                                                            *
 * include file for bezier code.                              *
 *------------------------------------------------------------*/

#ifndef _BEZ_PM_H_
#define _BEZ_PM_H_

#include "pm/pm.h"

namespace ProteinMechanica {

void
pm_BezierEval (int num_control_pts, PmVector3 *control_pts, int num_pts,
               PmVector3 *pts, PmVector3 *tang, PmVector3 *norm);

void
pm_BezierFit (int num_pts, PmVector3 *pts, float error, int *p_num_curves,
              PmVector3 **curves);

void
pm_BezierBernsteinGet (int degree, float u, float *coef);

}

#endif


