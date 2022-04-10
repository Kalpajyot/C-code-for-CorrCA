/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rank.c
 *
 * Code generation for function 'rank'
 *
 */

/* Include files */
#include "rank.h"
#include "rt_nonfinite.h"
#include "svd.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
int local_rank(const double A[6400])
{
  double s[80];
  double absx;
  int exponent;
  int i;
  int irank;
  boolean_T p;
  irank = 0;
  p = true;
  for (i = 0; i < 6400; i++) {
    if ((!p) || (rtIsInf(A[i]) || rtIsNaN(A[i]))) {
      p = false;
    }
  }
  if (p) {
    svd(A, s);
  } else {
    for (i = 0; i < 80; i++) {
      s[i] = rtNaN;
    }
  }
  absx = fabs(s[0]);
  if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
    if (absx <= 2.2250738585072014E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = ldexp(1.0, exponent - 53);
    }
  } else {
    absx = rtNaN;
  }
  absx *= 80.0;
  i = 0;
  while ((i < 80) && (s[i] > absx)) {
    irank++;
    i++;
  }
  return irank;
}

/* End of code generation (rank.c) */
