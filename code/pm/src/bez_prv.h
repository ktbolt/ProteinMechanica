
/*------------------------------------------------------------*
 *                                                            *
 *                    ****  bez  ****                         *
 *                                                            *
 * include file for bezier code.                              *
 *------------------------------------------------------------*/

#ifndef _BEZ_PRV_PM_H_
#define _BEZ_PRV_PM_H_

#define MAXPOINTS  1000
#include "pm/pm.h"

namespace ProteinMechanica {

static void
fit_left_tangent_comp (int end, PmVector3 *pts, PmVector3& tang);

static void
fit_right_tangent_comp (int end, PmVector3 *pts, PmVector3& tang);

static void
fit_cubic (PmVector3 *pts, int first, int last, PmVector3 tang_1,
           PmVector3 tang_2, float error, int *p_num_curves, PmVector3 *curves);

static void
fit_chord_length_parameterize (PmVector3 *pts, int first, int last, float **p_u);

static void
generate_bezier (PmVector3 *pts, int first, int last, float *uPrime, 
                 PmVector3 tang_1, PmVector3 tang_2, PmVector3 **curve);

static float
bezmult (int n, float u);

static void
fit_bezier_eval (int degree, PmVector3 *V, float t, PmVector3& pt);

static void
fit_max_error_comp (PmVector3 *pts, int first, int last,  PmVector3 *bezcurve,
                    float *u, int *split_pt, float *error);


static void
fit_reparameterize (PmVector3 *pts, int first, int last, float *u, 
                    PmVector3 *bezcurve, float **new_u);

static float
fit_root_find (PmVector3 *Q, PmVector3 P, float u);

static void
fit_center_tangent_comp (PmVector3 *pts, int center, PmVector3& tang);


static void
fit_output_curve (int n, PmVector3 *bezcurve, int *p_num_curves, PmVector3 *curves);

void
pm_BezierBernsteinGet (int degree, float u, float *coef);
}

#endif


