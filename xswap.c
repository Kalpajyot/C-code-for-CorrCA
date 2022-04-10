/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xswap.c
 *
 * Code generation for function 'xswap'
 *
 */

/* Include files */
#include "xswap.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void xswap(double x[6400], int ix0, int iy0)
{
  double temp;
  int i;
  int k;
  int temp_tmp;
  for (k = 0; k < 80; k++) {
    temp_tmp = (ix0 + k) - 1;
    temp = x[temp_tmp];
    i = (iy0 + k) - 1;
    x[temp_tmp] = x[i];
    x[i] = temp;
  }
}

/* End of code generation (xswap.c) */
