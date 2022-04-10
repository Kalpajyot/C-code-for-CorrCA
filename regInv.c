/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * regInv.c
 *
 * Code generation for function 'regInv'
 *
 */

/* Include files */
#include "regInv.h"
#include "rt_nonfinite.h"
#include "svd.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void regInv(const double R[6400], double K, double invR[6400])
{
  static double U[6400];
  static double V[6400];
  double b_data[6400];
  double y_data[6400];
  double s[80];
  double z_data[80];
  double bkj;
  int aoffset;
  int b_i;
  int b_loop_ub;
  int b_size_idx_0_tmp;
  int boffset;
  int coffset;
  int i;
  int inner;
  int j;
  int k;
  int loop_ub;
  boolean_T p;
  /* invR = regInv(R, K) */
  /*    PCA regularized inverse of square symmetric positive definite matrix R
   */
  p = true;
  for (k = 0; k < 6400; k++) {
    if ((!p) || (rtIsInf(R[k]) || rtIsNaN(R[k]))) {
      p = false;
    }
  }
  if (p) {
    b_svd(R, U, s, V);
  } else {
    for (loop_ub = 0; loop_ub < 6400; loop_ub++) {
      U[loop_ub] = rtNaN;
    }
    for (i = 0; i < 80; i++) {
      s[i] = rtNaN;
    }
    for (loop_ub = 0; loop_ub < 6400; loop_ub++) {
      V[loop_ub] = rtNaN;
    }
  }
  if (1.0 > K) {
    b_loop_ub = 0;
  } else {
    b_loop_ub = (int)K;
  }
  for (loop_ub = 0; loop_ub < b_loop_ub; loop_ub++) {
    z_data[loop_ub] = 1.0 / s[loop_ub];
  }
  b_size_idx_0_tmp = (signed char)b_loop_ub;
  loop_ub = (signed char)b_loop_ub * (signed char)b_loop_ub;
  if (0 <= loop_ub - 1) {
    memset(&b_data[0], 0, loop_ub * sizeof(double));
  }
  for (j = 0; j < b_loop_ub; j++) {
    b_data[j + (signed char)b_loop_ub * j] = z_data[j];
  }
  if (1.0 > K) {
    b_i = 0;
  } else {
    b_i = (int)K;
  }
  inner = b_i - 1;
  for (j = 0; j < b_size_idx_0_tmp; j++) {
    coffset = j * 80;
    boffset = j * (signed char)b_loop_ub;
    memset(&y_data[coffset], 0, 80U * sizeof(double));
    for (k = 0; k <= inner; k++) {
      aoffset = k * 80;
      bkj = b_data[boffset + k];
      for (i = 0; i < 80; i++) {
        loop_ub = aoffset + i;
        b_i = coffset + i;
        y_data[b_i] += U[loop_ub % 80 + 80 * (loop_ub / 80)] * bkj;
      }
    }
  }
  for (j = 0; j < 80; j++) {
    coffset = j * 80;
    memset(&invR[coffset], 0, 80U * sizeof(double));
    for (k = 0; k < b_size_idx_0_tmp; k++) {
      aoffset = k * 80;
      loop_ub = k * 80 + j;
      bkj = V[loop_ub % 80 + 80 * (loop_ub / 80)];
      for (i = 0; i < 80; i++) {
        loop_ub = coffset + i;
        invR[loop_ub] += y_data[aoffset + i] * bkj;
      }
    }
  }
}

/* End of code generation (regInv.c) */
