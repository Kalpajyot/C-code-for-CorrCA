/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xaxpy.c
 *
 * Code generation for function 'xaxpy'
 *
 */

/* Include files */
#include "xaxpy.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void b_xaxpy(int n, double a, const double x[6400], int ix0, double y[80],
             int iy0)
{
  int i;
  int i1;
  int k;
  if (!(a == 0.0)) {
    i = n - 1;
    for (k = 0; k <= i; k++) {
      i1 = (iy0 + k) - 1;
      y[i1] += a * x[(ix0 + k) - 1];
    }
  }
}

void c_xaxpy(int n, double a, const double x[80], int ix0, double y[6400],
             int iy0)
{
  int i;
  int i1;
  int k;
  if (!(a == 0.0)) {
    i = n - 1;
    for (k = 0; k <= i; k++) {
      i1 = (iy0 + k) - 1;
      y[i1] += a * x[(ix0 + k) - 1];
    }
  }
}

void xaxpy(int n, double a, int ix0, double y[6400], int iy0)
{
  int i;
  int i1;
  int k;
  if (!(a == 0.0)) {
    i = n - 1;
    for (k = 0; k <= i; k++) {
      i1 = (iy0 + k) - 1;
      y[i1] += a * y[(ix0 + k) - 1];
    }
  }
}

/* End of code generation (xaxpy.c) */
