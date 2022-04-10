/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * permute.c
 *
 * Code generation for function 'permute'
 *
 */

/* Include files */
#include "permute.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void permute(const double a[640000], double b[640000])
{
  int b_k;
  int b_tmp;
  int c_k;
  int k;
  for (k = 0; k < 10; k++) {
    b_tmp = 64000 * k;
    for (b_k = 0; b_k < 80; b_k++) {
      for (c_k = 0; c_k < 10; c_k++) {
        memcpy(&b[(b_k * 80 + c_k * 6400) + b_tmp],
               &a[(b_k * 800 + c_k * 80) + b_tmp], 80U * sizeof(double));
      }
    }
  }
}

/* End of code generation (permute.c) */
