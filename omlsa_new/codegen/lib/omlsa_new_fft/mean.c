/*
 * File: mean.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 09-Nov-2018 14:14:28
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "omlsa_new_fft.h"
#include "mean.h"

/* Function Definitions */

/*
 * Arguments    : const double x[128]
 * Return Type  : double
 */
double mean(const double x[128])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 127; k++) {
    y += x[k + 1];
  }

  y /= 128.0;
  return y;
}

/*
 * File trailer for mean.c
 *
 * [EOF]
 */
