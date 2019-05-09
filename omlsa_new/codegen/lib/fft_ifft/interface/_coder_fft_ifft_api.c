/*
 * File: _coder_fft_ifft_api.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Nov-2018 11:09:30
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_fft_ifft_api.h"
#include "_coder_fft_ifft_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true, false, 131434U, NULL, "fft_ifft", NULL,
  false, { 2045744189U, 2170104910U, 2743257031U, 4284093946U }, NULL };

/* Function Declarations */
static real32_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[320];
static real32_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[320];
static real32_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *y, const
  char_T *identifier))[320];
static const mxArray *emlrt_marshallOut(const real32_T u[320]);

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real32_T (*)[320]
 */
static real32_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[320]
{
  real32_T (*y)[320];
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real32_T (*)[320]
 */
  static real32_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[320]
{
  real32_T (*ret)[320];
  static const int32_T dims[1] = { 320 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "single", false, 1U, dims);
  ret = (real32_T (*)[320])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *y
 *                const char_T *identifier
 * Return Type  : real32_T (*)[320]
 */
static real32_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *y, const
  char_T *identifier))[320]
{
  real32_T (*b_y)[320];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_y = b_emlrt_marshallIn(sp, emlrtAlias(y), &thisId);
  emlrtDestroyArray(&y);
  return b_y;
}
/*
 * Arguments    : const real32_T u[320]
 * Return Type  : const mxArray *
 */
  static const mxArray *emlrt_marshallOut(const real32_T u[320])
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[1] = { 0 };

  static const int32_T iv1[1] = { 320 };

  y = NULL;
  m0 = emlrtCreateNumericArray(1, iv0, mxSINGLE_CLASS, mxREAL);
  mxSetData((mxArray *)m0, (void *)u);
  emlrtSetDimensions((mxArray *)m0, iv1, 1);
  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const mxArray *prhs[1]
 *                const mxArray *plhs[1]
 * Return Type  : void
 */
void fft_ifft_api(const mxArray *prhs[1], const mxArray *plhs[1])
{
  real32_T (*x)[320];
  real32_T (*y)[320];
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  x = (real32_T (*)[320])mxMalloc(sizeof(real32_T [320]));
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, false, -1);

  /* Marshall function inputs */
  y = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "y");

  /* Invoke the target function */
  fft_ifft(*y, *x);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*x);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void fft_ifft_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  fft_ifft_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void fft_ifft_initialize(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void fft_ifft_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_fft_ifft_api.c
 *
 * [EOF]
 */
