/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlarf.c
 *
 * Code generation for function 'xzlarf'
 *
 */

/* Include files */
#include "xzlarf.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Declarations */
static int div_nde_s32_floor(int numerator, int denominator);

/* Function Definitions */
static int div_nde_s32_floor(int numerator, int denominator)
{
  int b_numerator;
  if (((numerator < 0) != (denominator < 0)) &&
      (numerator % denominator != 0)) {
    b_numerator = -1;
  } else {
    b_numerator = 0;
  }
  return numerator / denominator + b_numerator;
}

void xzlarf(int m, int n, int iv0, double tau, double C[6400], int ic0,
            double work[80])
{
  double c;
  int b_i;
  int exitg1;
  int i;
  int ia;
  int iac;
  int lastc;
  int lastv;
  boolean_T exitg2;
  if (tau != 0.0) {
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C[i - 2] == 0.0)) {
      lastv--;
      i--;
    }
    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      i = ic0 + lastc * 80;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);
      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = -1;
  }
  if (lastv > 0) {
    if (lastc + 1 != 0) {
      if (0 <= lastc) {
        memset(&work[0], 0, (lastc + 1) * sizeof(double));
      }
      b_i = ic0 + 80 * lastc;
      for (iac = ic0; iac <= b_i; iac += 80) {
        c = 0.0;
        i = (iac + lastv) - 1;
        for (ia = iac; ia <= i; ia++) {
          c += C[ia - 1] * C[((iv0 + ia) - iac) - 1];
        }
        i = div_nde_s32_floor(iac - ic0, 80);
        work[i] += c;
      }
    }
    if (!(-tau == 0.0)) {
      i = ic0;
      for (iac = 0; iac <= lastc; iac++) {
        if (work[iac] != 0.0) {
          c = work[iac] * -tau;
          b_i = lastv + i;
          for (ia = i; ia < b_i; ia++) {
            C[ia - 1] += C[((iv0 + ia) - i) - 1] * c;
          }
        }
        i += 80;
      }
    }
  }
}

/* End of code generation (xzlarf.c) */
