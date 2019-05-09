/*
 * File: main.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 24-Apr-2019 10:59:21
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include Files */
#include "rt_nonfinite.h"
#include "nlms.h"
#include "main.h"
#include "nlms_terminate.h"
#include "nlms_initialize.h"

/* Function Declarations */
static void argInit_128x1_real32_T(float result[128]);
static void argInit_256x1_real32_T(float result[256]);
static void argInit_384x1_real32_T(float result[384]);
static float argInit_real32_T(void);
static void main_nlms(void);

/* Function Definitions */

/*
 * Arguments    : float result[128]
 * Return Type  : void
 */
static void argInit_128x1_real32_T(float result[128])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 128; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real32_T();
  }
}

/*
 * Arguments    : float result[256]
 * Return Type  : void
 */
static void argInit_256x1_real32_T(float result[256])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 256; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real32_T();
  }
}

/*
 * Arguments    : float result[384]
 * Return Type  : void
 */
static void argInit_384x1_real32_T(float result[384])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 384; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real32_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : float
 */
static float argInit_real32_T(void)
{
  return 0.0F;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_nlms(void)
{
  float fv0[128];
  float fv1[384];
  float fv2[256];
  float fv3[256];
  double e[128];

  /* Initialize function 'nlms' input arguments. */
  /* Initialize function input argument 'd'. */
  /* Initialize function input argument 'u'. */
  /* Initialize function input argument 'Pest'. */
  /* Initialize function input argument 'A_st'. */
  /* Call the entry-point 'nlms'. */
  argInit_128x1_real32_T(fv0);
  argInit_384x1_real32_T(fv1);
  argInit_256x1_real32_T(fv2);
  argInit_256x1_real32_T(fv3);
  nlms(fv0, fv1, fv2, fv3, e);
}

/*
 * Arguments    : int argc
 *                const char * const argv[]
 * Return Type  : int
 */
int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  nlms_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_nlms();

  /* Terminate the application.
     You do not need to do this more than one time. */
  nlms_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
