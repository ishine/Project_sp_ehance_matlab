/*
 * File: _coder_nlms_api.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Apr-2019 10:59:21
 */

/* Include Files */
#include "tmwtypes.h"
#include "_coder_nlms_api.h"
#include "_coder_nlms_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131466U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "nlms",                              /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

/* Function Declarations */
static real32_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[128];
static real32_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const char_T *identifier))[384];
static real32_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[384];
static real32_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Pest,
  const char_T *identifier))[256];
static real32_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *d, const
  char_T *identifier))[128];
static const mxArray *emlrt_marshallOut(const real_T u[128]);
static real32_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[256];
static real32_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[128];
static real32_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[384];
static real32_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[256];

/* Function Definitions */

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real32_T (*)[128]
 */
static real32_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[128]
{
  real32_T (*y)[128];
  y = g_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const char_T *identifier
 * Return Type  : real32_T (*)[384]
 */
  static real32_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const char_T *identifier))[384]
{
  real32_T (*y)[384];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(u), &thisId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real32_T (*)[384]
 */
static real32_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[384]
{
  real32_T (*y)[384];
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *Pest
 *                const char_T *identifier
 * Return Type  : real32_T (*)[256]
 */
  static real32_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Pest,
  const char_T *identifier))[256]
{
  real32_T (*y)[256];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(Pest), &thisId);
  emlrtDestroyArray(&Pest);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *d
 *                const char_T *identifier
 * Return Type  : real32_T (*)[128]
 */
static real32_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *d, const
  char_T *identifier))[128]
{
  real32_T (*y)[128];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(d), &thisId);
  emlrtDestroyArray(&d);
  return y;
}
/*
 * Arguments    : const real_T u[128]
 * Return Type  : const mxArray *
 */
  static const mxArray *emlrt_marshallOut(const real_T u[128])
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[1] = { 0 };

  static const int32_T iv1[1] = { 128 };

  y = NULL;
  m0 = emlrtCreateNumericArray(1, iv0, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m0, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m0, iv1, 1);
  emlrtAssign(&y, m0);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real32_T (*)[256]
 */
static real32_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[256]
{
  real32_T (*y)[256];
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real32_T (*)[128]
 */
  static real32_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[128]
{
  real32_T (*ret)[128];
  static const int32_T dims[1] = { 128 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "single", false, 1U, dims);
  ret = (real32_T (*)[128])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real32_T (*)[384]
 */
static real32_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[384]
{
  real32_T (*ret)[384];
  static const int32_T dims[1] = { 384 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "single", false, 1U, dims);
  ret = (real32_T (*)[384])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real32_T (*)[256]
 */
  static real32_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[256]
{
  real32_T (*ret)[256];
  static const int32_T dims[1] = { 256 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "single", false, 1U, dims);
  ret = (real32_T (*)[256])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const mxArray * const prhs[5]
 *                int32_T nlhs
 *                const mxArray *plhs[1]
 * Return Type  : void
 */
void nlms_api(const mxArray * const prhs[5], int32_T nlhs, const mxArray *plhs[1])
{
  real_T (*e)[128];
  const mxArray *prhs_copy_idx_2;
  const mxArray *prhs_copy_idx_3;
  static const uint32_T N[4] = { 3889729348U, 2275829390U, 1900332195U,
    1174496241U };

  real32_T (*d)[128];
  real32_T (*u)[384];
  real32_T (*Pest)[256];
  real32_T (*A_st)[256];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  (void)nlhs;
  st.tls = emlrtRootTLSGlobal;
  e = (real_T (*)[128])mxMalloc(sizeof(real_T [128]));
  prhs_copy_idx_2 = emlrtProtectR2012b(prhs[2], 2, false, -1);
  prhs_copy_idx_3 = emlrtProtectR2012b(prhs[3], 3, false, -1);

  /* Check constant function inputs */
  emlrtCheckArrayChecksumR2014a(&st, "N", N, prhs[4], false);

  /* Marshall function inputs */
  d = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "d");
  u = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "u");
  Pest = e_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_2), "Pest");
  A_st = e_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_3), "A_st");

  /* Invoke the target function */
  nlms(*d, *u, *Pest, *A_st, *e);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*e);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void nlms_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  nlms_xil_terminate();
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void nlms_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

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
void nlms_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * File trailer for _coder_nlms_api.c
 *
 * [EOF]
 */
