/*
 * File: _coder_nlms_api.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Apr-2019 10:59:21
 */

#ifndef _CODER_NLMS_API_H
#define _CODER_NLMS_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_nlms_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void nlms(real32_T d[128], real32_T u[384], real32_T Pest[256], real32_T
                 A_st[256], real_T e[128]);
extern void nlms_api(const mxArray * const prhs[5], int32_T nlhs, const mxArray *
                     plhs[1]);
extern void nlms_atexit(void);
extern void nlms_initialize(void);
extern void nlms_terminate(void);
extern void nlms_xil_terminate(void);

#endif

/*
 * File trailer for _coder_nlms_api.h
 *
 * [EOF]
 */
