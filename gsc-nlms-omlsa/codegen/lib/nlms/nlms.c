/*
 * File: nlms.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Apr-2019 10:59:21
 */

/* Include Files */
#include <string.h>
#include "rt_nonfinite.h"
#include "nlms.h"

/* Function Definitions */

/*
 * Arguments    : const float d[128]
 *                const float u[384]
 *                float Pest[256]
 *                float A_st[256]
 *                double e[128]
 * Return Type  : void
 */
void nlms(const float d[128], const float u[384], float Pest[256], float A_st
          [256], double e[128])
{
  int k;
  float Y_Frame_Block_data[383];
  float y;
  int loop_ub;
  float yout;
  float b_Pest;
  for (k = 0; k < 128; k++) {
    memcpy(&Y_Frame_Block_data[0], &u[k], (unsigned int)(256 * (int)sizeof(float)));

    /*     %% apply the coeff              */
    y = 0.0F;
    for (loop_ub = 0; loop_ub < 256; loop_ub++) {
      y += A_st[loop_ub] * Y_Frame_Block_data[loop_ub];
    }

    yout = d[k] - y;
    e[k] = yout;

    /*     %% update the coeff   */
    y = 0.0F;
    for (loop_ub = 0; loop_ub < 256; loop_ub++) {
      y += Y_Frame_Block_data[loop_ub] * Y_Frame_Block_data[loop_ub];
    }

    y *= 0.102F;
    for (loop_ub = 0; loop_ub < 256; loop_ub++) {
      b_Pest = 0.898F * Pest[loop_ub] + y;
      Pest[loop_ub] = b_Pest;
      A_st[loop_ub] += 0.072F / (b_Pest + 0.01F) * yout *
        Y_Frame_Block_data[loop_ub];
    }
  }
}

/*
 * File trailer for nlms.c
 *
 * [EOF]
 */
