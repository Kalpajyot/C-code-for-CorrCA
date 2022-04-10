/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eig.c
 *
 * Code generation for function 'eig'
 *
 */

/* Include files */
#include "eig.h"
#include "C_corrca_data.h"
#include "C_corrca_rtwutil.h"
#include "rt_nonfinite.h"
#include "xdhseqr.h"
#include "xnrm2.h"
#include "xzggev.h"
#include "xzlarf.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void eig(const double A[6400], creal_T V[6400], creal_T D[6400])
{
  static creal_T At[6400];
  creal_T alpha1[80];
  creal_T beta1[80];
  double b_D[6400];
  double b_V[6400];
  double work[80];
  double tau[79];
  double absxk;
  double ai;
  double ar;
  double brm;
  double colnorm;
  double scale;
  double t;
  int alpha1_tmp;
  int b_i;
  int c_i;
  int exitg1;
  int i;
  int i1;
  int ia;
  int im1n_tmp;
  int in;
  int iv0;
  int j;
  int k;
  int knt;
  int lastc;
  int lastv;
  boolean_T exitg2;
  boolean_T p;
  p = true;
  for (k = 0; k < 6400; k++) {
    if ((!p) || (rtIsInf(A[k]) || rtIsNaN(A[k]))) {
      p = false;
    }
  }
  if (!p) {
    for (i = 0; i < 6400; i++) {
      V[i].re = rtNaN;
      V[i].im = 0.0;
      D[i].re = 0.0;
      D[i].im = 0.0;
    }
    for (k = 0; k < 80; k++) {
      i = k + 80 * k;
      D[i].re = rtNaN;
      D[i].im = 0.0;
    }
  } else {
    p = true;
    j = 0;
    exitg2 = false;
    while ((!exitg2) && (j < 80)) {
      b_i = 0;
      do {
        exitg1 = 0;
        if (b_i <= j) {
          if (!(A[b_i + 80 * j] == A[j + 80 * b_i])) {
            p = false;
            exitg1 = 1;
          } else {
            b_i++;
          }
        } else {
          j++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);
      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
    if (p) {
      p = true;
      for (k = 0; k < 6400; k++) {
        if ((!p) || (rtIsInf(A[k]) || rtIsNaN(A[k]))) {
          p = false;
        }
      }
      if (!p) {
        for (i = 0; i < 6400; i++) {
          b_V[i] = rtNaN;
        }
        c_i = 2;
        for (j = 0; j < 79; j++) {
          if (c_i <= 80) {
            memset(&b_V[(j * 80 + c_i) + -1], 0, (81 - c_i) * sizeof(double));
          }
          c_i++;
        }
        for (i = 0; i < 6400; i++) {
          b_D[i] = rtNaN;
        }
      } else {
        memcpy(&b_D[0], &A[0], 6400U * sizeof(double));
        memset(&work[0], 0, 80U * sizeof(double));
        for (b_i = 0; b_i < 79; b_i++) {
          im1n_tmp = b_i * 80 + 2;
          in = (b_i + 1) * 80;
          alpha1_tmp = (b_i + 80 * b_i) + 1;
          t = b_D[alpha1_tmp];
          if (b_i + 3 < 80) {
            c_i = b_i + 1;
          } else {
            c_i = 78;
          }
          j = c_i + im1n_tmp;
          tau[b_i] = 0.0;
          scale = xnrm2(78 - b_i, b_D, j);
          if (scale != 0.0) {
            absxk = rt_hypotd_snf(t, scale);
            if (t >= 0.0) {
              absxk = -absxk;
            }
            if (fabs(absxk) < 1.0020841800044864E-292) {
              knt = -1;
              i = (j - b_i) + 77;
              do {
                knt++;
                for (k = j; k <= i; k++) {
                  b_D[k - 1] *= 9.9792015476736E+291;
                }
                absxk *= 9.9792015476736E+291;
                t *= 9.9792015476736E+291;
              } while (!(fabs(absxk) >= 1.0020841800044864E-292));
              absxk = rt_hypotd_snf(t, xnrm2(78 - b_i, b_D, j));
              if (t >= 0.0) {
                absxk = -absxk;
              }
              tau[b_i] = (absxk - t) / absxk;
              scale = 1.0 / (t - absxk);
              i = (j - b_i) + 77;
              for (k = j; k <= i; k++) {
                b_D[k - 1] *= scale;
              }
              for (k = 0; k <= knt; k++) {
                absxk *= 1.0020841800044864E-292;
              }
              t = absxk;
            } else {
              tau[b_i] = (absxk - t) / absxk;
              scale = 1.0 / (t - absxk);
              i = (j - b_i) + 77;
              for (k = j; k <= i; k++) {
                b_D[k - 1] *= scale;
              }
              t = absxk;
            }
          }
          b_D[alpha1_tmp] = 1.0;
          iv0 = (b_i + im1n_tmp) - 1;
          k = in + 1;
          if (tau[b_i] != 0.0) {
            lastv = 78 - b_i;
            c_i = (iv0 - b_i) + 78;
            while ((lastv + 1 > 0) && (b_D[c_i] == 0.0)) {
              lastv--;
              c_i--;
            }
            lastc = 80;
            exitg2 = false;
            while ((!exitg2) && (lastc > 0)) {
              knt = in + lastc;
              ia = knt;
              do {
                exitg1 = 0;
                if (ia <= knt + lastv * 80) {
                  if (b_D[ia - 1] != 0.0) {
                    exitg1 = 1;
                  } else {
                    ia += 80;
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
            lastv = -1;
            lastc = 0;
          }
          if (lastv + 1 > 0) {
            if (lastc != 0) {
              if (0 <= lastc - 1) {
                memset(&work[0], 0, lastc * sizeof(double));
              }
              knt = iv0;
              i = (in + 80 * lastv) + 1;
              for (c_i = k; c_i <= i; c_i += 80) {
                i1 = (c_i + lastc) - 1;
                for (ia = c_i; ia <= i1; ia++) {
                  j = ia - c_i;
                  work[j] += b_D[ia - 1] * b_D[knt];
                }
                knt++;
              }
            }
            if (!(-tau[b_i] == 0.0)) {
              knt = in;
              for (j = 0; j <= lastv; j++) {
                scale = b_D[iv0 + j];
                if (scale != 0.0) {
                  scale *= -tau[b_i];
                  i = knt + 1;
                  i1 = lastc + knt;
                  for (c_i = i; c_i <= i1; c_i++) {
                    b_D[c_i - 1] += work[(c_i - knt) - 1] * scale;
                  }
                }
                knt += 80;
              }
            }
          }
          xzlarf(79 - b_i, 79 - b_i, b_i + im1n_tmp, tau[b_i], b_D,
                 (b_i + in) + 2, work);
          b_D[alpha1_tmp] = t;
        }
        memcpy(&b_V[0], &b_D[0], 6400U * sizeof(double));
        for (j = 78; j >= 0; j--) {
          ia = (j + 1) * 80;
          for (b_i = 0; b_i <= j; b_i++) {
            b_V[ia + b_i] = 0.0;
          }
          i = j + 3;
          for (b_i = i; b_i < 81; b_i++) {
            c_i = ia + b_i;
            b_V[c_i - 1] = b_V[c_i - 81];
          }
        }
        memset(&b_V[0], 0, 80U * sizeof(double));
        b_V[0] = 1.0;
        memset(&work[0], 0, 80U * sizeof(double));
        for (b_i = 78; b_i >= 0; b_i--) {
          c_i = (b_i + b_i * 80) + 81;
          if (b_i + 1 < 79) {
            b_V[c_i] = 1.0;
            xzlarf(79 - b_i, 78 - b_i, c_i + 1, tau[b_i], b_V, c_i + 81, work);
            j = c_i + 2;
            i = (c_i - b_i) + 79;
            for (k = j; k <= i; k++) {
              b_V[k - 1] *= -tau[b_i];
            }
          }
          b_V[c_i] = 1.0 - tau[b_i];
          for (j = 0; j < b_i; j++) {
            b_V[(c_i - j) - 1] = 0.0;
          }
        }
        eml_dlahqr(b_D, b_V);
        c_i = 4;
        for (j = 0; j < 77; j++) {
          if (c_i <= 80) {
            memset(&b_D[(j * 80 + c_i) + -1], 0, (81 - c_i) * sizeof(double));
          }
          c_i++;
        }
      }
      for (j = 0; j < 79; j++) {
        b_D[(j + 80 * j) + 1] = 0.0;
        for (b_i = 0; b_i <= j; b_i++) {
          b_D[b_i + 80 * (j + 1)] = 0.0;
        }
      }
      for (i = 0; i < 6400; i++) {
        V[i].re = b_V[i];
        V[i].im = 0.0;
        D[i].re = b_D[i];
        D[i].im = 0.0;
      }
    } else {
      for (i = 0; i < 6400; i++) {
        At[i].re = A[i];
        At[i].im = 0.0;
      }
      xzggev(At, &c_i, alpha1, beta1, V);
      for (c_i = 0; c_i <= 6320; c_i += 80) {
        colnorm = 0.0;
        scale = 3.3121686421112381E-170;
        knt = c_i + 80;
        for (k = c_i + 1; k <= knt; k++) {
          absxk = fabs(V[k - 1].re);
          if (absxk > scale) {
            t = scale / absxk;
            colnorm = colnorm * t * t + 1.0;
            scale = absxk;
          } else {
            t = absxk / scale;
            colnorm += t * t;
          }
          absxk = fabs(V[k - 1].im);
          if (absxk > scale) {
            t = scale / absxk;
            colnorm = colnorm * t * t + 1.0;
            scale = absxk;
          } else {
            t = absxk / scale;
            colnorm += t * t;
          }
        }
        colnorm = scale * sqrt(colnorm);
        for (j = c_i + 1; j <= knt; j++) {
          ar = V[j - 1].re;
          ai = V[j - 1].im;
          if (ai == 0.0) {
            absxk = ar / colnorm;
            scale = 0.0;
          } else if (ar == 0.0) {
            absxk = 0.0;
            scale = ai / colnorm;
          } else {
            absxk = ar / colnorm;
            scale = ai / colnorm;
          }
          V[j - 1].re = absxk;
          V[j - 1].im = scale;
        }
      }
      memset(&D[0], 0, 6400U * sizeof(creal_T));
      for (k = 0; k < 80; k++) {
        ar = alpha1[k].re;
        ai = alpha1[k].im;
        t = beta1[k].re;
        colnorm = beta1[k].im;
        if (colnorm == 0.0) {
          if (ai == 0.0) {
            i = k + 80 * k;
            D[i].re = ar / t;
            D[i].im = 0.0;
          } else if (ar == 0.0) {
            i = k + 80 * k;
            D[i].re = 0.0;
            D[i].im = ai / t;
          } else {
            i = k + 80 * k;
            D[i].re = ar / t;
            D[i].im = ai / t;
          }
        } else if (t == 0.0) {
          if (ar == 0.0) {
            i = k + 80 * k;
            D[i].re = ai / colnorm;
            D[i].im = 0.0;
          } else if (ai == 0.0) {
            i = k + 80 * k;
            D[i].re = 0.0;
            D[i].im = -(ar / colnorm);
          } else {
            i = k + 80 * k;
            D[i].re = ai / colnorm;
            D[i].im = -(ar / colnorm);
          }
        } else {
          brm = fabs(t);
          scale = fabs(colnorm);
          if (brm > scale) {
            absxk = colnorm / t;
            scale = t + absxk * colnorm;
            i = k + 80 * k;
            D[i].re = (ar + absxk * ai) / scale;
            D[i].im = (ai - absxk * ar) / scale;
          } else if (scale == brm) {
            if (t > 0.0) {
              absxk = 0.5;
            } else {
              absxk = -0.5;
            }
            if (colnorm > 0.0) {
              scale = 0.5;
            } else {
              scale = -0.5;
            }
            i = k + 80 * k;
            D[i].re = (ar * absxk + ai * scale) / brm;
            D[i].im = (ai * absxk - ar * scale) / brm;
          } else {
            absxk = t / colnorm;
            scale = colnorm + absxk * t;
            i = k + 80 * k;
            D[i].re = (absxk * ar + ai) / scale;
            D[i].im = (absxk * ai - ar) / scale;
          }
        }
      }
    }
  }
}

/* End of code generation (eig.c) */
