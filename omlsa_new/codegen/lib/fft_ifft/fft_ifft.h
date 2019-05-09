/*
 * File: fft_ifft.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Nov-2018 11:09:30
 */

#ifndef FFT_IFFT_H
#define FFT_IFFT_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "fft_ifft_types.h"

/* Function Declarations */
extern void b_r2br_r2dit_trig(const creal32_T x[1024], const float costab[513],
  const float sintab[513], creal32_T y[1024]);
extern void fft_ifft(const float y[320], float x[320]);
extern void r2br_r2dit_trig(const creal32_T x[639], const float costab[513],
  const float sintab[513], creal32_T y[1024]);
extern void r2br_r2dit_trig_impl(const creal32_T x[320], int xoffInit, const
  float costab[513], const float sintab[513], creal32_T y[1024]);

#endif

/*
 * File trailer for fft_ifft.h
 *
 * [EOF]
 */
