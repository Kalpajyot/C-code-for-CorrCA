/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * svd.c
 *
 * Code generation for function 'svd'
 *
 */

/* Include files */
#include "svd.h"
#include "rt_nonfinite.h"
#include "xaxpy.h"
#include "xdotc.h"
#include "xnrm2.h"
#include "xrot.h"
#include "xrotg.h"
#include "xswap.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void b_svd(const double A[6400], double U[6400], double s[80], double V[6400])
{
  double b_A[6400];
  double e[80];
  double work[80];
  double nrm;
  double rt;
  double scale;
  double sm;
  double snorm;
  double sqds;
  double ztest;
  int exitg1;
  int k;
  int m;
  int q;
  int qjj;
  int qp1;
  int qp1jj;
  int qq;
  boolean_T apply_transform;
  boolean_T exitg2;
  memcpy(&b_A[0], &A[0], 6400U * sizeof(double));
  memset(&s[0], 0, 80U * sizeof(double));
  memset(&e[0], 0, 80U * sizeof(double));
  memset(&work[0], 0, 80U * sizeof(double));
  memset(&U[0], 0, 6400U * sizeof(double));
  memset(&V[0], 0, 6400U * sizeof(double));
  for (q = 0; q < 79; q++) {
    qp1 = q + 2;
    qq = (q + 80 * q) + 1;
    apply_transform = false;
    nrm = xnrm2(80 - q, b_A, qq);
    if (nrm > 0.0) {
      apply_transform = true;
      if (b_A[qq - 1] < 0.0) {
        ztest = -nrm;
        s[q] = -nrm;
      } else {
        ztest = nrm;
        s[q] = nrm;
      }
      if (fabs(ztest) >= 1.0020841800044864E-292) {
        nrm = 1.0 / ztest;
        qp1jj = (qq - q) + 79;
        for (k = qq; k <= qp1jj; k++) {
          b_A[k - 1] *= nrm;
        }
      } else {
        qp1jj = (qq - q) + 79;
        for (k = qq; k <= qp1jj; k++) {
          b_A[k - 1] /= s[q];
        }
      }
      b_A[qq - 1]++;
      s[q] = -s[q];
    } else {
      s[q] = 0.0;
    }
    for (k = qp1; k < 81; k++) {
      qjj = q + 80 * (k - 1);
      if (apply_transform) {
        xaxpy(80 - q, -(xdotc(80 - q, b_A, qq, b_A, qjj + 1) / b_A[q + 80 * q]),
              qq, b_A, qjj + 1);
      }
      e[k - 1] = b_A[qjj];
    }
    for (k = q + 1; k < 81; k++) {
      qjj = (k + 80 * q) - 1;
      U[qjj] = b_A[qjj];
    }
    if (q + 1 <= 78) {
      nrm = b_xnrm2(79 - q, e, q + 2);
      if (nrm == 0.0) {
        e[q] = 0.0;
      } else {
        if (e[q + 1] < 0.0) {
          e[q] = -nrm;
        } else {
          e[q] = nrm;
        }
        nrm = e[q];
        if (fabs(e[q]) >= 1.0020841800044864E-292) {
          nrm = 1.0 / e[q];
          for (k = qp1; k < 81; k++) {
            e[k - 1] *= nrm;
          }
        } else {
          for (k = qp1; k < 81; k++) {
            e[k - 1] /= nrm;
          }
        }
        e[q + 1]++;
        e[q] = -e[q];
        for (k = qp1; k < 81; k++) {
          work[k - 1] = 0.0;
        }
        for (k = qp1; k < 81; k++) {
          b_xaxpy(79 - q, e[k - 1], b_A, (q + 80 * (k - 1)) + 2, work, q + 2);
        }
        for (k = qp1; k < 81; k++) {
          c_xaxpy(79 - q, -e[k - 1] / e[q + 1], work, q + 2, b_A,
                  (q + 80 * (k - 1)) + 2);
        }
      }
      for (k = qp1; k < 81; k++) {
        V[(k + 80 * q) - 1] = e[k - 1];
      }
    }
  }
  m = 78;
  s[79] = b_A[6399];
  e[78] = b_A[6398];
  e[79] = 0.0;
  memset(&U[6320], 0, 80U * sizeof(double));
  U[6399] = 1.0;
  for (q = 78; q >= 0; q--) {
    qp1 = q + 2;
    qq = q + 80 * q;
    if (s[q] != 0.0) {
      for (k = qp1; k < 81; k++) {
        qjj = (q + 80 * (k - 1)) + 1;
        xaxpy(80 - q, -(xdotc(80 - q, U, qq + 1, U, qjj) / U[qq]), qq + 1, U,
              qjj);
      }
      for (k = q + 1; k < 81; k++) {
        qjj = (k + 80 * q) - 1;
        U[qjj] = -U[qjj];
      }
      U[qq]++;
      for (k = 0; k < q; k++) {
        U[k + 80 * q] = 0.0;
      }
    } else {
      memset(&U[q * 80], 0, 80U * sizeof(double));
      U[qq] = 1.0;
    }
  }
  for (q = 79; q >= 0; q--) {
    if ((q + 1 <= 78) && (e[q] != 0.0)) {
      qp1 = q + 2;
      qjj = (q + 80 * q) + 2;
      for (k = qp1; k < 81; k++) {
        qp1jj = (q + 80 * (k - 1)) + 2;
        xaxpy(79 - q, -(xdotc(79 - q, V, qjj, V, qp1jj) / V[qjj - 1]), qjj, V,
              qp1jj);
      }
    }
    memset(&V[q * 80], 0, 80U * sizeof(double));
    V[q + 80 * q] = 1.0;
  }
  qq = 0;
  snorm = 0.0;
  for (q = 0; q < 80; q++) {
    ztest = s[q];
    if (ztest != 0.0) {
      rt = fabs(ztest);
      nrm = ztest / rt;
      s[q] = rt;
      if (q + 1 < 80) {
        e[q] /= nrm;
      }
      qjj = 80 * q;
      qp1jj = qjj + 80;
      for (k = qjj + 1; k <= qp1jj; k++) {
        U[k - 1] *= nrm;
      }
    }
    if (q + 1 < 80) {
      ztest = e[q];
      if (ztest != 0.0) {
        rt = fabs(ztest);
        nrm = rt / ztest;
        e[q] = rt;
        s[q + 1] *= nrm;
        qjj = 80 * (q + 1);
        qp1jj = qjj + 80;
        for (k = qjj + 1; k <= qp1jj; k++) {
          V[k - 1] *= nrm;
        }
      }
    }
    nrm = fabs(s[q]);
    ztest = fabs(e[q]);
    if ((nrm > ztest) || rtIsNaN(ztest)) {
      ztest = nrm;
    }
    if ((!(snorm > ztest)) && (!rtIsNaN(ztest))) {
      snorm = ztest;
    }
  }
  while ((m + 2 > 0) && (qq < 75)) {
    k = m;
    do {
      exitg1 = 0;
      q = k + 1;
      if (k + 1 == 0) {
        exitg1 = 1;
      } else {
        nrm = fabs(e[k]);
        if ((nrm <= 2.2204460492503131E-16 * (fabs(s[k]) + fabs(s[k + 1]))) ||
            (nrm <= 1.0020841800044864E-292) ||
            ((qq > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
          e[k] = 0.0;
          exitg1 = 1;
        } else {
          k--;
        }
      }
    } while (exitg1 == 0);
    if (k + 1 == m + 1) {
      qjj = 4;
    } else {
      qp1jj = m + 2;
      qjj = m + 2;
      exitg2 = false;
      while ((!exitg2) && (qjj >= k + 1)) {
        qp1jj = qjj;
        if (qjj == k + 1) {
          exitg2 = true;
        } else {
          nrm = 0.0;
          if (qjj < m + 2) {
            nrm = fabs(e[qjj - 1]);
          }
          if (qjj > k + 2) {
            nrm += fabs(e[qjj - 2]);
          }
          ztest = fabs(s[qjj - 1]);
          if ((ztest <= 2.2204460492503131E-16 * nrm) ||
              (ztest <= 1.0020841800044864E-292)) {
            s[qjj - 1] = 0.0;
            exitg2 = true;
          } else {
            qjj--;
          }
        }
      }
      if (qp1jj == k + 1) {
        qjj = 3;
      } else if (qp1jj == m + 2) {
        qjj = 1;
      } else {
        qjj = 2;
        q = qp1jj;
      }
    }
    switch (qjj) {
    case 1:
      ztest = e[m];
      e[m] = 0.0;
      qp1jj = m + 1;
      for (k = qp1jj; k >= q + 1; k--) {
        xrotg(&s[k - 1], &ztest, &sm, &sqds);
        if (k > q + 1) {
          rt = e[k - 2];
          ztest = -sqds * rt;
          e[k - 2] = rt * sm;
        }
        xrot(V, 80 * (k - 1) + 1, 80 * (m + 1) + 1, sm, sqds);
      }
      break;
    case 2:
      ztest = e[q - 1];
      e[q - 1] = 0.0;
      for (k = q + 1; k <= m + 2; k++) {
        xrotg(&s[k - 1], &ztest, &sm, &sqds);
        rt = e[k - 1];
        ztest = -sqds * rt;
        e[k - 1] = rt * sm;
        xrot(U, 80 * (k - 1) + 1, 80 * (q - 1) + 1, sm, sqds);
      }
      break;
    case 3:
      qjj = m + 1;
      nrm = s[m + 1];
      scale = fabs(nrm);
      ztest = fabs(s[m]);
      if ((!(scale > ztest)) && (!rtIsNaN(ztest))) {
        scale = ztest;
      }
      ztest = fabs(e[m]);
      if ((!(scale > ztest)) && (!rtIsNaN(ztest))) {
        scale = ztest;
      }
      ztest = fabs(s[q]);
      if ((!(scale > ztest)) && (!rtIsNaN(ztest))) {
        scale = ztest;
      }
      ztest = fabs(e[q]);
      if ((!(scale > ztest)) && (!rtIsNaN(ztest))) {
        scale = ztest;
      }
      sm = nrm / scale;
      nrm = s[m] / scale;
      ztest = e[m] / scale;
      sqds = s[q] / scale;
      rt = ((nrm + sm) * (nrm - sm) + ztest * ztest) / 2.0;
      nrm = sm * ztest;
      nrm *= nrm;
      if ((rt != 0.0) || (nrm != 0.0)) {
        ztest = sqrt(rt * rt + nrm);
        if (rt < 0.0) {
          ztest = -ztest;
        }
        ztest = nrm / (rt + ztest);
      } else {
        ztest = 0.0;
      }
      ztest += (sqds + sm) * (sqds - sm);
      nrm = sqds * (e[q] / scale);
      for (k = q + 1; k <= qjj; k++) {
        xrotg(&ztest, &nrm, &sm, &sqds);
        if (k > q + 1) {
          e[k - 2] = ztest;
        }
        nrm = e[k - 1];
        rt = s[k - 1];
        e[k - 1] = sm * nrm - sqds * rt;
        ztest = sqds * s[k];
        s[k] *= sm;
        xrot(V, 80 * (k - 1) + 1, 80 * k + 1, sm, sqds);
        s[k - 1] = sm * rt + sqds * nrm;
        xrotg(&s[k - 1], &ztest, &sm, &sqds);
        ztest = sm * e[k - 1] + sqds * s[k];
        s[k] = -sqds * e[k - 1] + sm * s[k];
        nrm = sqds * e[k];
        e[k] *= sm;
        xrot(U, 80 * (k - 1) + 1, 80 * k + 1, sm, sqds);
      }
      e[m] = ztest;
      qq++;
      break;
    default:
      if (s[q] < 0.0) {
        s[q] = -s[q];
        qjj = 80 * q;
        qp1jj = qjj + 80;
        for (k = qjj + 1; k <= qp1jj; k++) {
          V[k - 1] = -V[k - 1];
        }
      }
      qp1 = q + 1;
      while ((q + 1 < 80) && (s[q] < s[qp1])) {
        rt = s[q];
        s[q] = s[qp1];
        s[qp1] = rt;
        xswap(V, 80 * q + 1, 80 * (q + 1) + 1);
        xswap(U, 80 * q + 1, 80 * (q + 1) + 1);
        q = qp1;
        qp1++;
      }
      qq = 0;
      m--;
      break;
    }
  }
}

void svd(const double A[6400], double U[80])
{
  double b_A[6400];
  double e[80];
  double work[80];
  double nrm;
  double rt;
  double scale;
  double sm;
  double snorm;
  double sqds;
  double ztest;
  int A_tmp;
  int exitg1;
  int i;
  int iter;
  int k;
  int m;
  int q;
  int qp1;
  int qq;
  int qs;
  boolean_T apply_transform;
  boolean_T exitg2;
  memcpy(&b_A[0], &A[0], 6400U * sizeof(double));
  memset(&U[0], 0, 80U * sizeof(double));
  memset(&e[0], 0, 80U * sizeof(double));
  memset(&work[0], 0, 80U * sizeof(double));
  for (q = 0; q < 79; q++) {
    qp1 = q + 2;
    qs = q + 80 * q;
    qq = qs + 1;
    iter = 79 - q;
    apply_transform = false;
    nrm = xnrm2(80 - q, b_A, qs + 1);
    if (nrm > 0.0) {
      apply_transform = true;
      if (b_A[qs] < 0.0) {
        ztest = -nrm;
        U[q] = -nrm;
      } else {
        ztest = nrm;
        U[q] = nrm;
      }
      if (fabs(ztest) >= 1.0020841800044864E-292) {
        nrm = 1.0 / ztest;
        i = (qs - q) + 80;
        for (k = qq; k <= i; k++) {
          b_A[k - 1] *= nrm;
        }
      } else {
        i = (qs - q) + 80;
        for (k = qq; k <= i; k++) {
          b_A[k - 1] /= U[q];
        }
      }
      b_A[qs]++;
      U[q] = -U[q];
    } else {
      U[q] = 0.0;
    }
    for (m = qp1; m < 81; m++) {
      qq = q + 80 * (m - 1);
      if (apply_transform) {
        nrm = 0.0;
        for (k = 0; k <= iter; k++) {
          nrm += b_A[qs + k] * b_A[qq + k];
        }
        nrm = -(nrm / b_A[qs]);
        if (!(nrm == 0.0)) {
          for (k = 0; k <= iter; k++) {
            A_tmp = qq + k;
            b_A[A_tmp] += nrm * b_A[qs + k];
          }
        }
      }
      e[m - 1] = b_A[qq];
    }
    if (q + 1 <= 78) {
      nrm = b_xnrm2(79 - q, e, q + 2);
      if (nrm == 0.0) {
        e[q] = 0.0;
      } else {
        if (e[q + 1] < 0.0) {
          e[q] = -nrm;
        } else {
          e[q] = nrm;
        }
        nrm = e[q];
        if (fabs(e[q]) >= 1.0020841800044864E-292) {
          nrm = 1.0 / e[q];
          for (k = qp1; k < 81; k++) {
            e[k - 1] *= nrm;
          }
        } else {
          for (k = qp1; k < 81; k++) {
            e[k - 1] /= nrm;
          }
        }
        e[q + 1]++;
        e[q] = -e[q];
        for (A_tmp = qp1; A_tmp < 81; A_tmp++) {
          work[A_tmp - 1] = 0.0;
        }
        for (m = qp1; m < 81; m++) {
          ztest = e[m - 1];
          if (!(ztest == 0.0)) {
            qq = q + 80 * (m - 1);
            i = 78 - q;
            for (k = 0; k <= i; k++) {
              A_tmp = (q + k) + 1;
              work[A_tmp] += ztest * b_A[(qq + k) + 1];
            }
          }
        }
        for (m = qp1; m < 81; m++) {
          c_xaxpy(79 - q, -e[m - 1] / e[q + 1], work, q + 2, b_A,
                  (q + 80 * (m - 1)) + 2);
        }
      }
    }
  }
  m = 78;
  U[79] = b_A[6399];
  e[78] = b_A[6398];
  e[79] = 0.0;
  iter = 0;
  snorm = 0.0;
  for (q = 0; q < 80; q++) {
    ztest = U[q];
    sm = ztest;
    if (ztest != 0.0) {
      rt = fabs(ztest);
      sm = rt;
      U[q] = rt;
      if (q + 1 < 80) {
        e[q] /= ztest / rt;
      }
    }
    if ((q + 1 < 80) && (e[q] != 0.0)) {
      rt = fabs(e[q]);
      nrm = rt / e[q];
      e[q] = rt;
      U[q + 1] *= nrm;
    }
    nrm = fabs(sm);
    ztest = fabs(e[q]);
    if ((nrm > ztest) || rtIsNaN(ztest)) {
      ztest = nrm;
    }
    if ((!(snorm > ztest)) && (!rtIsNaN(ztest))) {
      snorm = ztest;
    }
  }
  while ((m + 2 > 0) && (iter < 75)) {
    A_tmp = m;
    do {
      exitg1 = 0;
      q = A_tmp + 1;
      if (A_tmp + 1 == 0) {
        exitg1 = 1;
      } else {
        nrm = fabs(e[A_tmp]);
        if ((nrm <=
             2.2204460492503131E-16 * (fabs(U[A_tmp]) + fabs(U[A_tmp + 1]))) ||
            (nrm <= 1.0020841800044864E-292) ||
            ((iter > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
          e[A_tmp] = 0.0;
          exitg1 = 1;
        } else {
          A_tmp--;
        }
      }
    } while (exitg1 == 0);
    if (A_tmp + 1 == m + 1) {
      qq = 4;
    } else {
      qs = m + 2;
      qq = m + 2;
      exitg2 = false;
      while ((!exitg2) && (qq >= A_tmp + 1)) {
        qs = qq;
        if (qq == A_tmp + 1) {
          exitg2 = true;
        } else {
          nrm = 0.0;
          if (qq < m + 2) {
            nrm = fabs(e[qq - 1]);
          }
          if (qq > A_tmp + 2) {
            nrm += fabs(e[qq - 2]);
          }
          ztest = fabs(U[qq - 1]);
          if ((ztest <= 2.2204460492503131E-16 * nrm) ||
              (ztest <= 1.0020841800044864E-292)) {
            U[qq - 1] = 0.0;
            exitg2 = true;
          } else {
            qq--;
          }
        }
      }
      if (qs == A_tmp + 1) {
        qq = 3;
      } else if (qs == m + 2) {
        qq = 1;
      } else {
        qq = 2;
        q = qs;
      }
    }
    switch (qq) {
    case 1:
      ztest = e[m];
      e[m] = 0.0;
      i = m + 1;
      for (k = i; k >= q + 1; k--) {
        xrotg(&U[k - 1], &ztest, &rt, &sqds);
        if (k > q + 1) {
          sm = e[k - 2];
          ztest = -sqds * sm;
          e[k - 2] = sm * rt;
        }
      }
      break;
    case 2:
      ztest = e[q - 1];
      e[q - 1] = 0.0;
      for (k = q + 1; k <= m + 2; k++) {
        xrotg(&U[k - 1], &ztest, &rt, &sqds);
        sm = e[k - 1];
        ztest = -sqds * sm;
        e[k - 1] = sm * rt;
      }
      break;
    case 3:
      qq = m + 1;
      nrm = U[m + 1];
      scale = fabs(nrm);
      ztest = fabs(U[m]);
      if ((!(scale > ztest)) && (!rtIsNaN(ztest))) {
        scale = ztest;
      }
      ztest = fabs(e[m]);
      if ((!(scale > ztest)) && (!rtIsNaN(ztest))) {
        scale = ztest;
      }
      ztest = fabs(U[q]);
      if ((!(scale > ztest)) && (!rtIsNaN(ztest))) {
        scale = ztest;
      }
      ztest = fabs(e[q]);
      if ((!(scale > ztest)) && (!rtIsNaN(ztest))) {
        scale = ztest;
      }
      sm = nrm / scale;
      nrm = U[m] / scale;
      ztest = e[m] / scale;
      sqds = U[q] / scale;
      rt = ((nrm + sm) * (nrm - sm) + ztest * ztest) / 2.0;
      nrm = sm * ztest;
      nrm *= nrm;
      if ((rt != 0.0) || (nrm != 0.0)) {
        ztest = sqrt(rt * rt + nrm);
        if (rt < 0.0) {
          ztest = -ztest;
        }
        ztest = nrm / (rt + ztest);
      } else {
        ztest = 0.0;
      }
      ztest += (sqds + sm) * (sqds - sm);
      nrm = sqds * (e[q] / scale);
      for (k = q + 1; k <= qq; k++) {
        xrotg(&ztest, &nrm, &rt, &sqds);
        if (k > q + 1) {
          e[k - 2] = ztest;
        }
        nrm = e[k - 1];
        sm = U[k - 1];
        e[k - 1] = rt * nrm - sqds * sm;
        ztest = sqds * U[k];
        U[k] *= rt;
        U[k - 1] = rt * sm + sqds * nrm;
        xrotg(&U[k - 1], &ztest, &rt, &sqds);
        ztest = rt * e[k - 1] + sqds * U[k];
        U[k] = -sqds * e[k - 1] + rt * U[k];
        nrm = sqds * e[k];
        e[k] *= rt;
      }
      e[m] = ztest;
      iter++;
      break;
    default:
      if (U[q] < 0.0) {
        U[q] = -U[q];
      }
      qp1 = q + 1;
      while ((q + 1 < 80) && (U[q] < U[qp1])) {
        rt = U[q];
        U[q] = U[qp1];
        U[qp1] = rt;
        q = qp1;
        qp1++;
      }
      iter = 0;
      m--;
      break;
    }
  }
}

/* End of code generation (svd.c) */
