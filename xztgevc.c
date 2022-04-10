/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xztgevc.c
 *
 * Code generation for function 'xztgevc'
 *
 */

/* Include files */
#include "xztgevc.h"
#include "C_corrca_data.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void xztgevc(const creal_T A[6400], creal_T V[6400])
{
  creal_T work1[80];
  creal_T work2[80];
  double rworka[80];
  double acoeff;
  double ai;
  double anorm;
  double ascale;
  double brm;
  double d_im;
  double d_re;
  double dmin;
  double salpha_im;
  double salpha_re;
  double scale;
  double temp;
  double xmx;
  double y;
  int b_i;
  int i;
  int j;
  int je;
  int jr;
  int re_tmp;
  int x_tmp_tmp;
  boolean_T lscalea;
  boolean_T lscaleb;
  memset(&rworka[0], 0, 80U * sizeof(double));
  anorm = fabs(A[0].re) + fabs(A[0].im);
  for (j = 0; j < 79; j++) {
    for (i = 0; i <= j; i++) {
      re_tmp = i + 80 * (j + 1);
      rworka[j + 1] += fabs(A[re_tmp].re) + fabs(A[re_tmp].im);
    }
    re_tmp = (j + 80 * (j + 1)) + 1;
    y = rworka[j + 1] + (fabs(A[re_tmp].re) + fabs(A[re_tmp].im));
    if (y > anorm) {
      anorm = y;
    }
  }
  y = anorm;
  if (2.2250738585072014E-308 > anorm) {
    y = 2.2250738585072014E-308;
  }
  ascale = 1.0 / y;
  for (je = 0; je < 80; je++) {
    x_tmp_tmp = 80 * (79 - je);
    re_tmp = (x_tmp_tmp - je) + 79;
    xmx = A[re_tmp].re;
    scale = A[re_tmp].im;
    y = (fabs(xmx) + fabs(scale)) * ascale;
    if (1.0 > y) {
      y = 1.0;
    }
    temp = 1.0 / y;
    salpha_re = ascale * (temp * xmx);
    salpha_im = ascale * (temp * scale);
    acoeff = temp * ascale;
    if ((temp >= 2.2250738585072014E-308) &&
        (acoeff < 8.0166734400358911E-291)) {
      lscalea = true;
    } else {
      lscalea = false;
    }
    xmx = fabs(salpha_re) + fabs(salpha_im);
    if ((xmx >= 2.2250738585072014E-308) && (xmx < 8.0166734400358911E-291)) {
      lscaleb = true;
    } else {
      lscaleb = false;
    }
    scale = 1.0;
    if (lscalea) {
      y = anorm;
      if (1.2474001934592E+290 < anorm) {
        y = 1.2474001934592E+290;
      }
      scale = 8.0166734400358911E-291 / temp * y;
    }
    if (lscaleb) {
      y = 8.0166734400358911E-291 / xmx;
      if (y > scale) {
        scale = y;
      }
    }
    if (lscalea || lscaleb) {
      y = acoeff;
      if (1.0 > acoeff) {
        y = 1.0;
      }
      if (xmx > y) {
        y = xmx;
      }
      y = 1.0 / (2.2250738585072014E-308 * y);
      if (y < scale) {
        scale = y;
      }
      if (lscalea) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }
      salpha_re *= scale;
      salpha_im *= scale;
    }
    memset(&work1[0], 0, 80U * sizeof(creal_T));
    work1[79 - je].re = 1.0;
    work1[79 - je].im = 0.0;
    dmin = 2.2204460492503131E-16 * acoeff * anorm;
    y = 2.2204460492503131E-16 * (fabs(salpha_re) + fabs(salpha_im));
    if (y > dmin) {
      dmin = y;
    }
    if (2.2250738585072014E-308 > dmin) {
      dmin = 2.2250738585072014E-308;
    }
    b_i = 78 - je;
    for (jr = 0; jr <= b_i; jr++) {
      re_tmp = jr + x_tmp_tmp;
      work1[jr].re = acoeff * A[re_tmp].re;
      work1[jr].im = acoeff * A[re_tmp].im;
    }
    work1[79 - je].re = 1.0;
    work1[79 - je].im = 0.0;
    b_i = (int)(((-1.0 - ((-(double)je + 80.0) - 1.0)) + 1.0) / -1.0);
    for (j = 0; j < b_i; j++) {
      i = 78 - (je + j);
      re_tmp = i + 80 * i;
      d_re = acoeff * A[re_tmp].re - salpha_re;
      d_im = acoeff * A[re_tmp].im - salpha_im;
      if (fabs(d_re) + fabs(d_im) <= dmin) {
        d_re = dmin;
        d_im = 0.0;
      }
      brm = fabs(d_re);
      y = fabs(d_im);
      xmx = brm + y;
      if (xmx < 1.0) {
        scale = fabs(work1[i].re) + fabs(work1[i].im);
        if (scale >= 5.6177910464447375E+305 * xmx) {
          temp = 1.0 / scale;
          re_tmp = 79 - je;
          for (jr = 0; jr <= re_tmp; jr++) {
            work1[jr].re *= temp;
            work1[jr].im *= temp;
          }
        }
      }
      temp = -work1[i].re;
      ai = -work1[i].im;
      if (d_im == 0.0) {
        if (ai == 0.0) {
          y = temp / d_re;
          xmx = 0.0;
        } else if (temp == 0.0) {
          y = 0.0;
          xmx = ai / d_re;
        } else {
          y = temp / d_re;
          xmx = ai / d_re;
        }
      } else if (d_re == 0.0) {
        if (temp == 0.0) {
          y = ai / d_im;
          xmx = 0.0;
        } else if (ai == 0.0) {
          y = 0.0;
          xmx = -(temp / d_im);
        } else {
          y = ai / d_im;
          xmx = -(temp / d_im);
        }
      } else if (brm > y) {
        scale = d_im / d_re;
        xmx = d_re + scale * d_im;
        y = (temp + scale * ai) / xmx;
        xmx = (ai - scale * temp) / xmx;
      } else if (y == brm) {
        if (d_re > 0.0) {
          scale = 0.5;
        } else {
          scale = -0.5;
        }
        if (d_im > 0.0) {
          xmx = 0.5;
        } else {
          xmx = -0.5;
        }
        y = (temp * scale + ai * xmx) / brm;
        xmx = (ai * scale - temp * xmx) / brm;
      } else {
        scale = d_re / d_im;
        xmx = d_im + scale * d_re;
        y = (scale * temp + ai) / xmx;
        xmx = (scale * ai - temp) / xmx;
      }
      work1[i].re = y;
      work1[i].im = xmx;
      if (i + 1 > 1) {
        if (fabs(work1[i].re) + fabs(work1[i].im) > 1.0) {
          temp = 1.0 / (fabs(work1[i].re) + fabs(work1[i].im));
          if (acoeff * rworka[i] >= 5.6177910464447375E+305 * temp) {
            re_tmp = 79 - je;
            for (jr = 0; jr <= re_tmp; jr++) {
              work1[jr].re *= temp;
              work1[jr].im *= temp;
            }
          }
        }
        d_re = acoeff * work1[i].re;
        d_im = acoeff * work1[i].im;
        for (jr = 0; jr < i; jr++) {
          re_tmp = jr + 80 * i;
          work1[jr].re += d_re * A[re_tmp].re - d_im * A[re_tmp].im;
          work1[jr].im += d_re * A[re_tmp].im + d_im * A[re_tmp].re;
        }
      }
    }
    memset(&work2[0], 0, 80U * sizeof(creal_T));
    b_i = 79 - je;
    for (i = 0; i <= b_i; i++) {
      xmx = work1[i].re;
      scale = work1[i].im;
      for (jr = 0; jr < 80; jr++) {
        re_tmp = jr + 80 * i;
        work2[jr].re += V[re_tmp].re * xmx - V[re_tmp].im * scale;
        work2[jr].im += V[re_tmp].re * scale + V[re_tmp].im * xmx;
      }
    }
    xmx = fabs(work2[0].re) + fabs(work2[0].im);
    for (jr = 0; jr < 79; jr++) {
      y = fabs(work2[jr + 1].re) + fabs(work2[jr + 1].im);
      if (y > xmx) {
        xmx = y;
      }
    }
    if (xmx > 2.2250738585072014E-308) {
      temp = 1.0 / xmx;
      for (jr = 0; jr < 80; jr++) {
        b_i = jr + x_tmp_tmp;
        V[b_i].re = temp * work2[jr].re;
        V[b_i].im = temp * work2[jr].im;
      }
    } else {
      memset(&V[x_tmp_tmp], 0, 80U * sizeof(creal_T));
    }
  }
}

/* End of code generation (xztgevc.c) */
