/*
 * File: _coder_omlsa_new_fft_api.h
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 09-Nov-2018 14:14:28
 */

#ifndef _CODER_OMLSA_NEW_FFT_API_H
#define _CODER_OMLSA_NEW_FFT_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_omlsa_new_fft_api.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void omlsa_new_fft(real_T Y_BUF[512000], emxArray_real_T *y, real_T out
  [512]);
extern void omlsa_new_fft_api(const mxArray *prhs[1], const mxArray *plhs[2]);
extern void omlsa_new_fft_atexit(void);
extern void omlsa_new_fft_initialize(void);
extern void omlsa_new_fft_terminate(void);
extern void omlsa_new_fft_xil_terminate(void);

#endif

/*
 * File trailer for _coder_omlsa_new_fft_api.h
 *
 * [EOF]
 */
