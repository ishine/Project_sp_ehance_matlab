/*
 * File: nlms.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Apr-2019 10:59:21
 */

#ifndef NLMS_H
#define NLMS_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "nlms_types.h"

/* Function Declarations */
extern void nlms(const float d[128], const float u[384], float Pest[256], float
                 A_st[256], double e[128]);

#endif

/*
 * File trailer for nlms.h
 *
 * [EOF]
 */
