/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sort.c
 *
 * Code generation for function 'sort'
 *
 */

/* Include files */
#include "sort.h"
#include "rt_nonfinite.h"
#include "sortIdx.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void sort(double x[80], int idx[80])
{
  double xwork[80];
  double x4[4];
  double d;
  double d1;
  int iwork[80];
  int i1;
  int i2;
  int i3;
  int i4;
  int ib;
  int k;
  int nNaNs;
  signed char idx4[4];
  signed char perm[4];
  x4[0] = 0.0;
  idx4[0] = 0;
  x4[1] = 0.0;
  idx4[1] = 0;
  x4[2] = 0.0;
  idx4[2] = 0;
  x4[3] = 0.0;
  idx4[3] = 0;
  memset(&idx[0], 0, 80U * sizeof(int));
  memset(&xwork[0], 0, 80U * sizeof(double));
  nNaNs = 0;
  ib = 0;
  for (k = 0; k < 80; k++) {
    if (rtIsNaN(x[k])) {
      idx[79 - nNaNs] = k + 1;
      xwork[79 - nNaNs] = x[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = (signed char)(k + 1);
      x4[ib - 1] = x[k];
      if (ib == 4) {
        ib = k - nNaNs;
        if (x4[0] >= x4[1]) {
          i1 = 1;
          i2 = 2;
        } else {
          i1 = 2;
          i2 = 1;
        }
        if (x4[2] >= x4[3]) {
          i3 = 3;
          i4 = 4;
        } else {
          i3 = 4;
          i4 = 3;
        }
        d = x4[i1 - 1];
        d1 = x4[i3 - 1];
        if (d >= d1) {
          d = x4[i2 - 1];
          if (d >= d1) {
            perm[0] = (signed char)i1;
            perm[1] = (signed char)i2;
            perm[2] = (signed char)i3;
            perm[3] = (signed char)i4;
          } else if (d >= x4[i4 - 1]) {
            perm[0] = (signed char)i1;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)i2;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)i1;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)i2;
          }
        } else {
          d1 = x4[i4 - 1];
          if (d >= d1) {
            if (x4[i2 - 1] >= d1) {
              perm[0] = (signed char)i3;
              perm[1] = (signed char)i1;
              perm[2] = (signed char)i2;
              perm[3] = (signed char)i4;
            } else {
              perm[0] = (signed char)i3;
              perm[1] = (signed char)i1;
              perm[2] = (signed char)i4;
              perm[3] = (signed char)i2;
            }
          } else {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)i4;
            perm[2] = (signed char)i1;
            perm[3] = (signed char)i2;
          }
        }
        idx[ib - 3] = idx4[perm[0] - 1];
        idx[ib - 2] = idx4[perm[1] - 1];
        idx[ib - 1] = idx4[perm[2] - 1];
        idx[ib] = idx4[perm[3] - 1];
        x[ib - 3] = x4[perm[0] - 1];
        x[ib - 2] = x4[perm[1] - 1];
        x[ib - 1] = x4[perm[2] - 1];
        x[ib] = x4[perm[3] - 1];
        ib = 0;
      }
    }
  }
  if (ib > 0) {
    perm[1] = 0;
    perm[2] = 0;
    perm[3] = 0;
    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] >= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] >= x4[1]) {
      if (x4[1] >= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] >= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] >= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] >= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }
    for (k = 0; k < ib; k++) {
      i2 = perm[k] - 1;
      i1 = ((k - nNaNs) - ib) + 80;
      idx[i1] = idx4[i2];
      x[i1] = x4[i2];
    }
  }
  ib = (nNaNs >> 1) + 80;
  for (k = 0; k <= ib - 81; k++) {
    i2 = (k - nNaNs) + 80;
    i1 = idx[i2];
    idx[i2] = idx[79 - k];
    idx[79 - k] = i1;
    x[i2] = xwork[79 - k];
    x[79 - k] = xwork[i2];
  }
  if ((nNaNs & 1) != 0) {
    ib -= nNaNs;
    x[ib] = xwork[ib];
  }
  memset(&iwork[0], 0, 80U * sizeof(int));
  i4 = 80 - nNaNs;
  if (80 - nNaNs > 1) {
    i3 = (80 - nNaNs) >> 2;
    i1 = 4;
    while (i3 > 1) {
      if ((i3 & 1) != 0) {
        i3--;
        ib = i1 * i3;
        i2 = 80 - (nNaNs + ib);
        if (i2 > i1) {
          merge(idx, x, ib, i1, i2 - i1, iwork, xwork);
        }
      }
      ib = i1 << 1;
      i3 >>= 1;
      for (k = 0; k < i3; k++) {
        merge(idx, x, k * ib, i1, i1, iwork, xwork);
      }
      i1 = ib;
    }
    if (80 - nNaNs > i1) {
      merge(idx, x, 0, i1, 80 - (nNaNs + i1), iwork, xwork);
    }
  }
  if ((nNaNs > 0) && (80 - nNaNs > 0)) {
    for (k = 0; k < nNaNs; k++) {
      ib = (k - nNaNs) + 80;
      xwork[k] = x[ib];
      iwork[k] = idx[ib];
    }
    for (k = i4; k >= 1; k--) {
      ib = (nNaNs + k) - 1;
      x[ib] = x[k - 1];
      idx[ib] = idx[k - 1];
    }
    if (0 <= nNaNs - 1) {
      memcpy(&x[0], &xwork[0], nNaNs * sizeof(double));
      memcpy(&idx[0], &iwork[0], nNaNs * sizeof(int));
    }
  }
}

/* End of code generation (sort.c) */
