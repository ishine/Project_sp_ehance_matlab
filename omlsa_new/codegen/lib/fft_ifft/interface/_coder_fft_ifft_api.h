/*
 * File: _coder_fft_ifft_api.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Nov-2018 11:09:30
 */

#ifndef _CODER_FFT_IFFT_API_H
#define _CODER_FFT_IFFT_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_fft_ifft_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void fft_ifft(real32_T y[320], real32_T x[320]);
extern void fft_ifft_api(const mxArray *prhs[1], const mxArray *plhs[1]);
extern void fft_ifft_atexit(void);
extern void fft_ifft_initialize(void);
extern void fft_ifft_terminate(void);
extern void fft_ifft_xil_terminate(void);

#endif

/*
 * File trailer for _coder_fft_ifft_api.h
 *
 * [EOF]
 */
