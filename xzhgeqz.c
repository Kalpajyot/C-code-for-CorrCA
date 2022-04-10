/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhgeqz.c
 *
 * Code generation for function 'xzhgeqz'
 *
 */

/* Include files */
#include "xzhgeqz.h"
#include "C_corrca_data.h"
#include "rt_nonfinite.h"
#include "sqrt.h"
#include "xzlartg.h"
#include <math.h>

/* Function Definitions */
void xzhgeqz(creal_T A[6400], int ilo, int ihi, creal_T Z[6400], int *info,
             creal_T alpha1[80], creal_T beta1[80])
{
  creal_T b_ascale;
  creal_T ctemp;
  creal_T shift;
  double ad22_im;
  double ad22_re;
  double anorm;
  double ascale;
  double b_atol;
  double colscale;
  double colssq;
  double eshift_im;
  double eshift_re;
  double scale;
  double ssq;
  double t;
  double t1_im;
  double t1_im_tmp;
  double t1_re;
  int b_i;
  int col;
  int exitg1;
  int i;
  int ifirst;
  int iiter;
  int ilast;
  int ilastm1;
  int istart;
  int j;
  int jiter;
  int nm1;
  int row;
  boolean_T b_guard1 = false;
  boolean_T exitg2;
  boolean_T failed;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  *info = 0;
  for (i = 0; i < 80; i++) {
    alpha1[i].re = 0.0;
    alpha1[i].im = 0.0;
    beta1[i].re = 1.0;
    beta1[i].im = 0.0;
  }
  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  anorm = 0.0;
  if (ilo <= ihi) {
    scale = 3.3121686421112381E-170;
    ssq = 0.0;
    nm1 = ihi - ilo;
    for (j = 0; j <= nm1; j++) {
      colscale = 3.3121686421112381E-170;
      colssq = 0.0;
      col = (ilo + j) - 1;
      i = j + 1;
      if (i >= nm1) {
        i = nm1;
      }
      b_i = ilo + i;
      for (row = ilo; row <= b_i; row++) {
        i = (row + 80 * col) - 1;
        anorm = fabs(A[i].re);
        if (anorm > colscale) {
          t = colscale / anorm;
          colssq = colssq * t * t + 1.0;
          colscale = anorm;
        } else {
          t = anorm / colscale;
          colssq += t * t;
        }
        anorm = fabs(A[i].im);
        if (anorm > colscale) {
          t = colscale / anorm;
          colssq = colssq * t * t + 1.0;
          colscale = anorm;
        } else {
          t = anorm / colscale;
          colssq += t * t;
        }
      }
      if (scale >= colscale) {
        t = colscale / scale;
        ssq += t * t * colssq;
      } else {
        t = scale / colscale;
        ssq = colssq + t * t * ssq;
        scale = colscale;
      }
    }
    anorm = scale * sqrt(ssq);
  }
  t = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (t > 2.2250738585072014E-308) {
    b_atol = t;
  }
  t = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    t = anorm;
  }
  ascale = 1.0 / t;
  failed = true;
  b_i = ihi + 1;
  for (j = b_i; j < 81; j++) {
    alpha1[j - 1] = A[(j + 80 * (j - 1)) - 1];
  }
  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 0;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1) - 1) {
        b_guard1 = false;
        if (ilast + 1 == ilo) {
          goto60 = true;
          b_guard1 = true;
        } else {
          b_i = ilast + 80 * ilastm1;
          if (fabs(A[b_i].re) + fabs(A[b_i].im) <= b_atol) {
            A[b_i].re = 0.0;
            A[b_i].im = 0.0;
            goto60 = true;
            b_guard1 = true;
          } else {
            j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (j + 1 >= ilo)) {
              if (j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else {
                b_i = j + 80 * (j - 1);
                if (fabs(A[b_i].re) + fabs(A[b_i].im) <= b_atol) {
                  A[b_i].re = 0.0;
                  A[b_i].im = 0.0;
                  guard3 = true;
                  exitg2 = true;
                } else {
                  j--;
                  guard3 = false;
                }
              }
            }
            if (guard3) {
              ifirst = j + 1;
              goto70 = true;
            }
            if (goto70) {
              b_guard1 = true;
            } else {
              for (i = 0; i < 80; i++) {
                alpha1[i].re = rtNaN;
                alpha1[i].im = 0.0;
                beta1[i].re = rtNaN;
                beta1[i].im = 0.0;
              }
              for (b_i = 0; b_i < 6400; b_i++) {
                Z[b_i].re = rtNaN;
                Z[b_i].im = 0.0;
              }
              *info = 1;
              exitg1 = 1;
            }
          }
        }
        if (b_guard1) {
          if (goto60) {
            goto60 = false;
            alpha1[ilast] = A[ilast + 80 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              jiter++;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              if (iiter - iiter / 10 * 10 != 0) {
                i = ilastm1 + 80 * ilastm1;
                anorm = ascale * A[i].re;
                t = ascale * A[i].im;
                if (t == 0.0) {
                  shift.re = anorm / 0.11180339887498948;
                  shift.im = 0.0;
                } else if (anorm == 0.0) {
                  shift.re = 0.0;
                  shift.im = t / 0.11180339887498948;
                } else {
                  shift.re = anorm / 0.11180339887498948;
                  shift.im = t / 0.11180339887498948;
                }
                i = ilast + 80 * ilast;
                anorm = ascale * A[i].re;
                t = ascale * A[i].im;
                if (t == 0.0) {
                  ad22_re = anorm / 0.11180339887498948;
                  ad22_im = 0.0;
                } else if (anorm == 0.0) {
                  ad22_re = 0.0;
                  ad22_im = t / 0.11180339887498948;
                } else {
                  ad22_re = anorm / 0.11180339887498948;
                  ad22_im = t / 0.11180339887498948;
                }
                t1_re = 0.5 * (shift.re + ad22_re);
                t1_im = 0.5 * (shift.im + ad22_im);
                t1_im_tmp = t1_re * t1_im;
                i = ilastm1 + 80 * ilast;
                anorm = ascale * A[i].re;
                t = ascale * A[i].im;
                if (t == 0.0) {
                  colscale = anorm / 0.11180339887498948;
                  colssq = 0.0;
                } else if (anorm == 0.0) {
                  colscale = 0.0;
                  colssq = t / 0.11180339887498948;
                } else {
                  colscale = anorm / 0.11180339887498948;
                  colssq = t / 0.11180339887498948;
                }
                i = ilast + 80 * ilastm1;
                anorm = ascale * A[i].re;
                t = ascale * A[i].im;
                if (t == 0.0) {
                  ssq = anorm / 0.11180339887498948;
                  anorm = 0.0;
                } else if (anorm == 0.0) {
                  ssq = 0.0;
                  anorm = t / 0.11180339887498948;
                } else {
                  ssq = anorm / 0.11180339887498948;
                  anorm = t / 0.11180339887498948;
                }
                t = shift.re * ad22_re - shift.im * ad22_im;
                scale = shift.re * ad22_im + shift.im * ad22_re;
                shift.re = ((t1_re * t1_re - t1_im * t1_im) +
                            (colscale * ssq - colssq * anorm)) -
                           t;
                shift.im = ((t1_im_tmp + t1_im_tmp) +
                            (colscale * anorm + colssq * ssq)) -
                           scale;
                b_sqrt(&shift);
                if ((t1_re - ad22_re) * shift.re +
                        (t1_im - ad22_im) * shift.im <=
                    0.0) {
                  shift.re += t1_re;
                  shift.im += t1_im;
                } else {
                  shift.re = t1_re - shift.re;
                  shift.im = t1_im - shift.im;
                }
              } else {
                i = ilast + 80 * ilastm1;
                anorm = ascale * A[i].re;
                t = ascale * A[i].im;
                if (t == 0.0) {
                  colscale = anorm / 0.11180339887498948;
                  colssq = 0.0;
                } else if (anorm == 0.0) {
                  colscale = 0.0;
                  colssq = t / 0.11180339887498948;
                } else {
                  colscale = anorm / 0.11180339887498948;
                  colssq = t / 0.11180339887498948;
                }
                eshift_re += colscale;
                eshift_im += colssq;
                shift.re = eshift_re;
                shift.im = eshift_im;
              }
              j = ilastm1;
              i = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                nm1 = j + 80 * j;
                ctemp.re = ascale * A[nm1].re - shift.re * 0.11180339887498948;
                ctemp.im = ascale * A[nm1].im - shift.im * 0.11180339887498948;
                anorm = fabs(ctemp.re) + fabs(ctemp.im);
                i += 80 * j;
                t = ascale * (fabs(A[i].re) + fabs(A[i].im));
                scale = anorm;
                if (t > anorm) {
                  scale = t;
                }
                if ((scale < 1.0) && (scale != 0.0)) {
                  anorm /= scale;
                  t /= scale;
                }
                b_i = j + 80 * (j - 1);
                if ((fabs(A[b_i].re) + fabs(A[b_i].im)) * t <= anorm * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  i = j;
                  j--;
                }
              }
              if (!goto90) {
                istart = ifirst;
                nm1 = (ifirst + 80 * (ifirst - 1)) - 1;
                ctemp.re = ascale * A[nm1].re - shift.re * 0.11180339887498948;
                ctemp.im = ascale * A[nm1].im - shift.im * 0.11180339887498948;
              }
              goto90 = false;
              i = istart + 80 * (istart - 1);
              b_ascale.re = ascale * A[i].re;
              b_ascale.im = ascale * A[i].im;
              b_xzlartg(ctemp, b_ascale, &anorm, &shift);
              j = istart;
              nm1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  col = j + 80 * nm1;
                  xzlartg(A[col - 1], A[col], &anorm, &shift,
                          &A[(j + 80 * nm1) - 1]);
                  A[col].re = 0.0;
                  A[col].im = 0.0;
                }
                for (nm1 = j; nm1 < 81; nm1++) {
                  row = j + 80 * (nm1 - 1);
                  ad22_re = anorm * A[row - 1].re +
                            (shift.re * A[row].re - shift.im * A[row].im);
                  ad22_im = anorm * A[row - 1].im +
                            (shift.re * A[row].im + shift.im * A[row].re);
                  t = A[row - 1].im;
                  scale = A[row - 1].re;
                  A[row].re = anorm * A[row].re - (shift.re * A[row - 1].re +
                                                   shift.im * A[row - 1].im);
                  A[row].im =
                      anorm * A[row].im - (shift.re * t - shift.im * scale);
                  A[row - 1].re = ad22_re;
                  A[row - 1].im = ad22_im;
                }
                shift.re = -shift.re;
                shift.im = -shift.im;
                nm1 = j;
                if (ilast + 1 < j + 2) {
                  nm1 = ilast - 1;
                }
                for (i = 1; i <= nm1 + 2; i++) {
                  row = (i + 80 * (j - 1)) - 1;
                  col = (i + 80 * j) - 1;
                  ad22_re = anorm * A[col].re +
                            (shift.re * A[row].re - shift.im * A[row].im);
                  ad22_im = anorm * A[col].im +
                            (shift.re * A[row].im + shift.im * A[row].re);
                  t = A[col].im;
                  scale = A[col].re;
                  A[row].re = anorm * A[row].re -
                              (shift.re * A[col].re + shift.im * A[col].im);
                  A[row].im =
                      anorm * A[row].im - (shift.re * t - shift.im * scale);
                  A[col].re = ad22_re;
                  A[col].im = ad22_im;
                }
                for (i = 0; i < 80; i++) {
                  row = i + 80 * (j - 1);
                  col = i + 80 * j;
                  ad22_re = anorm * Z[col].re +
                            (shift.re * Z[row].re - shift.im * Z[row].im);
                  ad22_im = anorm * Z[col].im +
                            (shift.re * Z[row].im + shift.im * Z[row].re);
                  t = Z[col].im;
                  scale = Z[col].re;
                  Z[row].re = anorm * Z[row].re -
                              (shift.re * Z[col].re + shift.im * Z[col].im);
                  Z[row].im =
                      anorm * Z[row].im - (shift.re * t - shift.im * scale);
                  Z[col].re = ad22_re;
                  Z[col].im = ad22_im;
                }
                nm1 = j - 1;
                j++;
              }
            }
            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }
  if (guard2) {
    if (failed) {
      *info = ilast + 1;
      for (i = 0; i <= ilast; i++) {
        alpha1[i].re = rtNaN;
        alpha1[i].im = 0.0;
        beta1[i].re = rtNaN;
        beta1[i].im = 0.0;
      }
      for (b_i = 0; b_i < 6400; b_i++) {
        Z[b_i].re = rtNaN;
        Z[b_i].im = 0.0;
      }
    } else {
      guard1 = true;
    }
  }
  if (guard1) {
    for (j = 0; j <= ilo - 2; j++) {
      alpha1[j] = A[j + 80 * j];
    }
  }
}

/* End of code generation (xzhgeqz.c) */
