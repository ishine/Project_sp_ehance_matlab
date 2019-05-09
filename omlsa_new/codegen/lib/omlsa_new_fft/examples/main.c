/*
 * File: main.c
 *
 * MATLAB Coder version            : 3.1
 * C/C++ source code generated on  : 09-Nov-2018 14:14:28
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
#include "omlsa_new_fft.h"
#include "main.h"
#include "omlsa_new_fft_terminate.h"
#include "omlsa_new_fft_emxAPI.h"
#include "omlsa_new_fft_initialize.h"

/* Function Declarations */
static void argInit_512000x1_real_T(double result[512000]);
static double argInit_real_T(void);
static void main_omlsa_new_fft(void);

/* Function Definitions */

/*
 * Arguments    : double result[512000]
 * Return Type  : void
 */
static void argInit_512000x1_real_T(double result[512000])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 512000; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
static void main_omlsa_new_fft(void)
{
  emxArray_real_T *y;
  static double dv5[512000];
  double out[512];
  emxInitArray_real_T(&y, 1);

  /* Initialize function 'omlsa_new_fft' input arguments. */
  /* Initialize function input argument 'Y_BUF'. */
  /* Call the entry-point 'omlsa_new_fft'. */
  argInit_512000x1_real_T(dv5);
  omlsa_new_fft(dv5, y, out);
  emxDestroyArray_real_T(y);
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
  omlsa_new_fft_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_omlsa_new_fft();

  /* Terminate the application.
     You do not need to do this more than one time. */
  omlsa_new_fft_terminate();
  return 0;
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
