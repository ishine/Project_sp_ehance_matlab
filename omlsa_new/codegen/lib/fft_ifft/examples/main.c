/*
 * File: main.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 20-Nov-2018 11:09:30
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
#include "fft_ifft.h"
#include "main.h"
#include "fft_ifft_terminate.h"
#include "fft_ifft_initialize.h"

/* Function Declarations */
static void argInit_320x1_real32_T(float result[320]);
static float argInit_real32_T(void);
static void main_fft_ifft(void);

/* Function Definitions */

/*
 * Arguments    : float result[320]
 * Return Type  : void
 */
static void argInit_320x1_real32_T(float result[320])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 320; idx0++) {
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
static void main_fft_ifft(void)
{
  float fv0[320];
  float x[320];

  /* Initialize function 'fft_ifft' input arguments. */
  /* Initialize function input argument 'y'. */
  /* Call the entry-point 'fft_ifft'. */
  argInit_320x1_real32_T(fv0);
  fft_ifft(fv0, x);
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
  fft_ifft_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_fft_ifft();

  /* Terminate the application.
     You do not need to do this more than one time. */
  fft_ifft_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
