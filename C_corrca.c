/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * C_corrca.c
 *
 * Code generation for function 'C_corrca'
 *
 */

/* Include files */
#include "C_corrca.h"
#include "C_corrca_data.h"
#include "C_corrca_initialize.h"
#include "eig.h"
#include "permute.h"
#include "rank.h"
#include "regInv.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "sqrt.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void C_corrca(const double X[240000], creal_T ISC[3])
{
  static creal_T W[6400];
  static creal_T b_ISC[6400];
  static creal_T b_Rt[6400];
  static creal_T b_W[6400];
  static creal_T c_W[6400];
  static double Rkl[640000];
  static double dv[640000];
  static double x[240000];
  static double Rt[6400];
  static double Rw[6400];
  static double a[6400];
  static double b_a[6400];
  creal_T b_x[80];
  creal_T y[80];
  double a__1[80];
  double b_br;
  double bi;
  double brm;
  double im;
  double muj;
  double re;
  int iidx[80];
  int ar;
  int br;
  int i;
  int i1;
  int ic;
  int irank;
  int w;
  if (!isInitialized_C_corrca) {
    C_corrca_initialize();
  }
  /*  do compute W */
  /*  don't truncate the eigs */
  /*  exemplars, dimensions, subjects */
  /*  compute within- and between-subject covariances */
  memcpy(&x[0], &X[0], 240000U * sizeof(double));
  for (irank = 0; irank < 800; irank++) {
    muj = 0.0;
    for (br = 0; br < 300; br++) {
      muj += x[br + 300 * irank];
    }
    muj /= 300.0;
    for (br = 0; br < 300; br++) {
      ar = br + 300 * irank;
      x[ar] -= muj;
    }
  }
  memset(&dv[0], 0, 640000U * sizeof(double));
  for (irank = 0; irank <= 639200; irank += 800) {
    i = irank + 1;
    i1 = irank + 800;
    if (i <= i1) {
      memset(&dv[i + -1], 0, ((i1 - i) + 1) * sizeof(double));
    }
  }
  br = -1;
  for (irank = 0; irank <= 639200; irank += 800) {
    ar = -1;
    i = irank + 1;
    i1 = irank + 800;
    for (ic = i; ic <= i1; ic++) {
      muj = 0.0;
      for (w = 0; w < 300; w++) {
        muj += x[(w + ar) + 1] * x[(w + br) + 1];
      }
      dv[ic - 1] += 0.0033444816053511705 * muj;
      ar += 300;
    }
    br += 300;
  }
  permute(dv, Rkl);
  for (irank = 0; irank < 6400; irank++) {
    Rw[irank] = Rkl[(irank % 80 + 80 * (irank / 80 % 80)) +
                    6400 * (11 * (irank / 6400))];
  }
  for (ic = 0; ic < 9; ic++) {
    br = (ic + 1) * 6400;
    for (irank = 0; irank < 6400; irank++) {
      i = br + irank;
      Rw[irank] +=
          Rkl[(i % 80 + 80 * (i / 80 % 80)) + 6400 * (11 * (i / 6400))];
    }
  }
  /*  pooled over all subjects */
  memcpy(&Rt[0], &Rkl[0], 6400U * sizeof(double));
  for (ic = 0; ic < 99; ic++) {
    br = (ic + 1) * 6400;
    for (irank = 0; irank < 6400; irank++) {
      Rt[irank] += Rkl[br + irank];
    }
  }
  /*  pooled over all pairs of subjects */
  for (i = 0; i < 6400; i++) {
    Rt[i] = (Rt[i] - Rw[i]) / 9.0;
  }
  /*  find projections W that maximize ISC */
  irank = local_rank(Rw);
  /*  handle rank deficient data.  */
  if (irank >= 80) {
    irank = 80;
  }
  regInv(Rw, irank, a);
  for (i = 0; i < 80; i++) {
    for (i1 = 0; i1 < 80; i1++) {
      bi = 0.0;
      for (irank = 0; irank < 80; irank++) {
        bi += a[i + 80 * irank] * Rt[irank + 80 * i1];
      }
      b_a[i + 80 * i1] = bi;
    }
  }
  eig(b_a, W, b_ISC);
  /*  make sure they are sorted by ISC and W normalized */
  for (i = 0; i < 6400; i++) {
    a[i] = b_ISC[i].re;
  }
  for (ic = 0; ic < 80; ic++) {
    a__1[ic] = a[ic + 80 * ic];
  }
  sort(a__1, iidx);
  for (i = 0; i < 80; i++) {
    for (i1 = 0; i1 < 80; i1++) {
      b_W[i1 + 80 * i] = W[i1 + 80 * (iidx[i] - 1)];
    }
  }
  for (ic = 0; ic < 6400; ic++) {
    bi = b_W[ic].re;
    brm = b_W[ic].im;
    b_ISC[ic].re = bi * bi - brm * brm;
    bi *= brm;
    b_ISC[ic].im = bi + bi;
  }
  for (irank = 0; irank < 80; irank++) {
    ar = irank * 80;
    re = b_ISC[ar].re;
    im = b_ISC[ar].im;
    for (ic = 0; ic < 79; ic++) {
      br = (ar + ic) + 1;
      re += b_ISC[br].re;
      im += b_ISC[br].im;
    }
    y[irank].re = re;
    y[irank].im = im;
  }
  for (ic = 0; ic < 80; ic++) {
    b_sqrt(&y[ic]);
  }
  for (i = 0; i < 80; i++) {
    b_br = y[i].re;
    bi = y[i].im;
    if (bi == 0.0) {
      re = 1.0 / b_br;
      im = 0.0;
    } else if (b_br == 0.0) {
      re = 0.0;
      im = -(1.0 / bi);
    } else {
      brm = fabs(b_br);
      muj = fabs(bi);
      if (brm > muj) {
        brm = bi / b_br;
        muj = b_br + brm * bi;
        re = (brm * 0.0 + 1.0) / muj;
        im = (0.0 - brm) / muj;
      } else if (muj == brm) {
        if (b_br > 0.0) {
          b_br = 0.5;
        } else {
          b_br = -0.5;
        }
        if (bi > 0.0) {
          muj = 0.5;
        } else {
          muj = -0.5;
        }
        re = (b_br + 0.0 * muj) / brm;
        im = (0.0 * b_br - muj) / brm;
      } else {
        brm = b_br / bi;
        muj = bi + brm * b_br;
        re = brm / muj;
        im = (brm * 0.0 - 1.0) / muj;
      }
    }
    y[i].re = re;
    y[i].im = im;
  }
  memset(&b_ISC[0], 0, 6400U * sizeof(creal_T));
  for (irank = 0; irank < 80; irank++) {
    b_ISC[irank + 80 * irank] = y[irank];
  }
  for (i = 0; i < 80; i++) {
    for (i1 = 0; i1 < 80; i1++) {
      re = 0.0;
      im = 0.0;
      for (irank = 0; irank < 80; irank++) {
        br = i + 80 * irank;
        ar = irank + 80 * i1;
        bi = b_W[br].re;
        brm = b_W[br].im;
        muj = b_ISC[ar].re;
        b_br = b_ISC[ar].im;
        re += bi * muj - brm * b_br;
        im += bi * b_br + brm * muj;
      }
      irank = i + 80 * i1;
      W[irank].re = re;
      W[irank].im = im;
    }
  }
  /*  compute ISC for fixed W, or recompute with unregularized Rw */
  for (i = 0; i < 6400; i++) {
    b_W[i] = W[i];
    b_Rt[i].re = Rt[i];
    b_Rt[i].im = 0.0;
  }
  for (i = 0; i < 80; i++) {
    for (i1 = 0; i1 < 80; i1++) {
      re = 0.0;
      im = 0.0;
      for (irank = 0; irank < 80; irank++) {
        br = irank + 80 * i;
        muj = b_W[br].re;
        b_br = -b_W[br].im;
        br = irank + 80 * i1;
        bi = b_Rt[br].re;
        brm = b_Rt[br].im;
        re += muj * bi - b_br * brm;
        im += muj * brm + b_br * bi;
      }
      irank = i + 80 * i1;
      W[irank].re = re;
      W[irank].im = im;
    }
    for (i1 = 0; i1 < 80; i1++) {
      re = 0.0;
      im = 0.0;
      for (irank = 0; irank < 80; irank++) {
        br = i + 80 * irank;
        ar = irank + 80 * i1;
        bi = W[br].re;
        brm = W[br].im;
        muj = b_W[ar].re;
        b_br = b_W[ar].im;
        re += bi * muj - brm * b_br;
        im += bi * b_br + brm * muj;
      }
      irank = i + 80 * i1;
      b_ISC[irank].re = re;
      b_ISC[irank].im = im;
    }
  }
  for (i = 0; i < 6400; i++) {
    b_Rt[i].re = Rw[i];
    b_Rt[i].im = 0.0;
  }
  for (ic = 0; ic < 80; ic++) {
    for (i = 0; i < 80; i++) {
      re = 0.0;
      im = 0.0;
      for (i1 = 0; i1 < 80; i1++) {
        br = i1 + 80 * ic;
        muj = b_W[br].re;
        b_br = -b_W[br].im;
        br = i1 + 80 * i;
        bi = b_Rt[br].re;
        brm = b_Rt[br].im;
        re += muj * bi - b_br * brm;
        im += muj * brm + b_br * bi;
      }
      i1 = ic + 80 * i;
      c_W[i1].re = re;
      c_W[i1].im = im;
    }
    for (i = 0; i < 80; i++) {
      re = 0.0;
      im = 0.0;
      for (i1 = 0; i1 < 80; i1++) {
        br = ic + 80 * i1;
        ar = i1 + 80 * i;
        bi = c_W[br].re;
        brm = c_W[br].im;
        muj = b_W[ar].re;
        b_br = b_W[ar].im;
        re += bi * muj - brm * b_br;
        im += bi * b_br + brm * muj;
      }
      i1 = ic + 80 * i;
      W[i1].re = re;
      W[i1].im = im;
    }
    b_x[ic] = b_ISC[ic + 80 * ic];
  }
  for (ic = 0; ic < 80; ic++) {
    y[ic] = W[ic + 80 * ic];
  }
  if (y[0].im == 0.0) {
    if (b_x[0].im == 0.0) {
      ISC[0].re = b_x[0].re / y[0].re;
      ISC[0].im = 0.0;
    } else if (b_x[0].re == 0.0) {
      ISC[0].re = 0.0;
      ISC[0].im = b_x[0].im / y[0].re;
    } else {
      ISC[0].re = b_x[0].re / y[0].re;
      ISC[0].im = b_x[0].im / y[0].re;
    }
  } else if (y[0].re == 0.0) {
    if (b_x[0].re == 0.0) {
      ISC[0].re = b_x[0].im / y[0].im;
      ISC[0].im = 0.0;
    } else if (b_x[0].im == 0.0) {
      ISC[0].re = 0.0;
      ISC[0].im = -(b_x[0].re / y[0].im);
    } else {
      ISC[0].re = b_x[0].im / y[0].im;
      ISC[0].im = -(b_x[0].re / y[0].im);
    }
  } else {
    brm = fabs(y[0].re);
    muj = fabs(y[0].im);
    if (brm > muj) {
      brm = y[0].im / y[0].re;
      muj = y[0].re + brm * y[0].im;
      ISC[0].re = (b_x[0].re + brm * b_x[0].im) / muj;
      ISC[0].im = (b_x[0].im - brm * b_x[0].re) / muj;
    } else if (muj == brm) {
      if (y[0].re > 0.0) {
        b_br = 0.5;
      } else {
        b_br = -0.5;
      }
      if (y[0].im > 0.0) {
        muj = 0.5;
      } else {
        muj = -0.5;
      }
      ISC[0].re = (b_x[0].re * b_br + b_x[0].im * muj) / brm;
      ISC[0].im = (b_x[0].im * b_br - b_x[0].re * muj) / brm;
    } else {
      brm = y[0].re / y[0].im;
      muj = y[0].im + brm * y[0].re;
      ISC[0].re = (brm * b_x[0].re + b_x[0].im) / muj;
      ISC[0].im = (brm * b_x[0].im - b_x[0].re) / muj;
    }
  }
  if (y[1].im == 0.0) {
    if (b_x[1].im == 0.0) {
      ISC[1].re = b_x[1].re / y[1].re;
      ISC[1].im = 0.0;
    } else if (b_x[1].re == 0.0) {
      ISC[1].re = 0.0;
      ISC[1].im = b_x[1].im / y[1].re;
    } else {
      ISC[1].re = b_x[1].re / y[1].re;
      ISC[1].im = b_x[1].im / y[1].re;
    }
  } else if (y[1].re == 0.0) {
    if (b_x[1].re == 0.0) {
      ISC[1].re = b_x[1].im / y[1].im;
      ISC[1].im = 0.0;
    } else if (b_x[1].im == 0.0) {
      ISC[1].re = 0.0;
      ISC[1].im = -(b_x[1].re / y[1].im);
    } else {
      ISC[1].re = b_x[1].im / y[1].im;
      ISC[1].im = -(b_x[1].re / y[1].im);
    }
  } else {
    brm = fabs(y[1].re);
    muj = fabs(y[1].im);
    if (brm > muj) {
      brm = y[1].im / y[1].re;
      muj = y[1].re + brm * y[1].im;
      ISC[1].re = (b_x[1].re + brm * b_x[1].im) / muj;
      ISC[1].im = (b_x[1].im - brm * b_x[1].re) / muj;
    } else if (muj == brm) {
      if (y[1].re > 0.0) {
        b_br = 0.5;
      } else {
        b_br = -0.5;
      }
      if (y[1].im > 0.0) {
        muj = 0.5;
      } else {
        muj = -0.5;
      }
      ISC[1].re = (b_x[1].re * b_br + b_x[1].im * muj) / brm;
      ISC[1].im = (b_x[1].im * b_br - b_x[1].re * muj) / brm;
    } else {
      brm = y[1].re / y[1].im;
      muj = y[1].im + brm * y[1].re;
      ISC[1].re = (brm * b_x[1].re + b_x[1].im) / muj;
      ISC[1].im = (brm * b_x[1].im - b_x[1].re) / muj;
    }
  }
  if (y[2].im == 0.0) {
    if (b_x[2].im == 0.0) {
      ISC[2].re = b_x[2].re / y[2].re;
      ISC[2].im = 0.0;
    } else if (b_x[2].re == 0.0) {
      ISC[2].re = 0.0;
      ISC[2].im = b_x[2].im / y[2].re;
    } else {
      ISC[2].re = b_x[2].re / y[2].re;
      ISC[2].im = b_x[2].im / y[2].re;
    }
  } else if (y[2].re == 0.0) {
    if (b_x[2].re == 0.0) {
      ISC[2].re = b_x[2].im / y[2].im;
      ISC[2].im = 0.0;
    } else if (b_x[2].im == 0.0) {
      ISC[2].re = 0.0;
      ISC[2].im = -(b_x[2].re / y[2].im);
    } else {
      ISC[2].re = b_x[2].im / y[2].im;
      ISC[2].im = -(b_x[2].re / y[2].im);
    }
  } else {
    brm = fabs(y[2].re);
    muj = fabs(y[2].im);
    if (brm > muj) {
      brm = y[2].im / y[2].re;
      muj = y[2].re + brm * y[2].im;
      ISC[2].re = (b_x[2].re + brm * b_x[2].im) / muj;
      ISC[2].im = (b_x[2].im - brm * b_x[2].re) / muj;
    } else if (muj == brm) {
      if (y[2].re > 0.0) {
        b_br = 0.5;
      } else {
        b_br = -0.5;
      }
      if (y[2].im > 0.0) {
        muj = 0.5;
      } else {
        muj = -0.5;
      }
      ISC[2].re = (b_x[2].re * b_br + b_x[2].im * muj) / brm;
      ISC[2].im = (b_x[2].im * b_br - b_x[2].re * muj) / brm;
    } else {
      brm = y[2].re / y[2].im;
      muj = y[2].im + brm * y[2].re;
      ISC[2].re = (brm * b_x[2].re + b_x[2].im) / muj;
      ISC[2].im = (brm * b_x[2].im - b_x[2].re) / muj;
    }
  }
  /*  Y = zeros(T,D,N); */
  /*  projections into corrca space */
  /*  for l=N:-1:1 */
  /*      Y(:,:,l)=real(X(:,:,l)*W);  */
  /*  end */
}

/* End of code generation (C_corrca.c) */
