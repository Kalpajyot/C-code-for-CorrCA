/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xdhseqr.c
 *
 * Code generation for function 'xdhseqr'
 *
 */

/* Include files */
#include "xdhseqr.h"
#include "rt_nonfinite.h"
#include "xdlanv2.h"
#include "xzlarfg.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
int eml_dlahqr(double h[6400], double z[6400])
{
  double v[3];
  double aa;
  double ab;
  double ba;
  double bb;
  double rt1r;
  double rt2r;
  double s;
  double s_tmp;
  double tst;
  int L;
  int b_i;
  int b_k;
  int hoffset;
  int i;
  int info;
  int its;
  int j;
  int k;
  int m;
  int nr;
  int sum1_tmp;
  boolean_T exitg1;
  boolean_T exitg2;
  boolean_T exitg3;
  boolean_T goto150;
  info = 0;
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;
  for (j = 0; j < 77; j++) {
    i = j + 80 * j;
    h[i + 2] = 0.0;
    h[i + 3] = 0.0;
  }
  h[6239] = 0.0;
  b_i = 79;
  exitg1 = false;
  while ((!exitg1) && (b_i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 2401)) {
      k = b_i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > L)) {
        i = k + 80 * (k - 1);
        ba = fabs(h[i]);
        if (ba <= 8.0166734400358911E-291) {
          exitg3 = true;
        } else {
          hoffset = k + 80 * k;
          bb = fabs(h[hoffset]);
          tst = fabs(h[i - 1]) + bb;
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              tst = fabs(h[(k + 80 * (k - 2)) - 1]);
            }
            if (k + 2 <= 80) {
              tst += fabs(h[hoffset + 1]);
            }
          }
          if (ba <= 2.2204460492503131E-16 * tst) {
            tst = fabs(h[hoffset - 1]);
            if (ba > tst) {
              ab = ba;
              ba = tst;
            } else {
              ab = tst;
            }
            tst = fabs(h[i - 1] - h[hoffset]);
            if (bb > tst) {
              aa = bb;
              bb = tst;
            } else {
              aa = tst;
            }
            s = aa + ab;
            tst = 2.2204460492503131E-16 * (bb * (aa / s));
            if ((8.0166734400358911E-291 > tst) || rtIsNaN(tst)) {
              tst = 8.0166734400358911E-291;
            }
            if (ba * (ab / s) <= tst) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }
      L = k + 1;
      if (k + 1 > 1) {
        h[k + 80 * (k - 1)] = 0.0;
      }
      if (k + 1 >= b_i) {
        goto150 = true;
        exitg2 = true;
      } else {
        if (its == 10) {
          hoffset = k + 80 * k;
          s = fabs(h[hoffset + 1]) + fabs(h[(k + 80 * (k + 1)) + 2]);
          tst = 0.75 * s + h[hoffset];
          ab = -0.4375 * s;
          aa = s;
          ba = tst;
        } else if (its == 20) {
          s = fabs(h[b_i + 80 * (b_i - 1)]) +
              fabs(h[(b_i + 80 * (b_i - 2)) - 1]);
          tst = 0.75 * s + h[b_i + 80 * b_i];
          ab = -0.4375 * s;
          aa = s;
          ba = tst;
        } else {
          hoffset = b_i + 80 * (b_i - 1);
          tst = h[hoffset - 1];
          aa = h[hoffset];
          ab = h[(b_i + 80 * b_i) - 1];
          ba = h[b_i + 80 * b_i];
        }
        s = ((fabs(tst) + fabs(ab)) + fabs(aa)) + fabs(ba);
        if (s == 0.0) {
          rt1r = 0.0;
          bb = 0.0;
          rt2r = 0.0;
          ba = 0.0;
        } else {
          tst /= s;
          aa /= s;
          ab /= s;
          ba /= s;
          bb = (tst + ba) / 2.0;
          tst = (tst - bb) * (ba - bb) - ab * aa;
          aa = sqrt(fabs(tst));
          if (tst >= 0.0) {
            rt1r = bb * s;
            rt2r = rt1r;
            bb = aa * s;
            ba = -bb;
          } else {
            rt1r = bb + aa;
            rt2r = bb - aa;
            if (fabs(rt1r - ba) <= fabs(rt2r - ba)) {
              rt1r *= s;
              rt2r = rt1r;
            } else {
              rt2r *= s;
              rt1r = rt2r;
            }
            bb = 0.0;
            ba = 0.0;
          }
        }
        m = b_i - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          hoffset = m + 80 * (m - 1);
          tst = h[hoffset];
          s_tmp = h[hoffset - 1];
          aa = s_tmp - rt2r;
          s = (fabs(aa) + fabs(ba)) + fabs(tst);
          ab = tst / s;
          hoffset = m + 80 * m;
          v[0] =
              (ab * h[hoffset - 1] + (s_tmp - rt1r) * (aa / s)) - bb * (ba / s);
          tst = h[hoffset];
          v[1] = ab * (((s_tmp + tst) - rt1r) - rt2r);
          v[2] = ab * h[hoffset + 1];
          s = (fabs(v[0]) + fabs(v[1])) + fabs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
          if (m == k + 1) {
            exitg3 = true;
          } else {
            i = m + 80 * (m - 2);
            if (fabs(h[i - 1]) * (fabs(v[1]) + fabs(v[2])) <=
                2.2204460492503131E-16 * fabs(v[0]) *
                    ((fabs(h[i - 2]) + fabs(s_tmp)) + fabs(tst))) {
              exitg3 = true;
            } else {
              m--;
            }
          }
        }
        for (b_k = m; b_k <= b_i; b_k++) {
          nr = (b_i - b_k) + 2;
          if (3 < nr) {
            nr = 3;
          }
          if (b_k > m) {
            hoffset = (b_k + 80 * (b_k - 2)) - 1;
            for (j = 0; j < nr; j++) {
              v[j] = h[j + hoffset];
            }
          }
          tst = v[0];
          bb = xzlarfg(nr, &tst, v);
          v[0] = tst;
          if (b_k > m) {
            h[(b_k + 80 * (b_k - 2)) - 1] = tst;
            i = b_k + 80 * (b_k - 2);
            h[i] = 0.0;
            if (b_k < b_i) {
              h[i + 1] = 0.0;
            }
          } else if (m > k + 1) {
            h[(b_k + 80 * (b_k - 2)) - 1] *= 1.0 - bb;
          }
          s = v[1];
          aa = bb * v[1];
          if (nr == 3) {
            rt1r = v[2];
            tst = bb * v[2];
            for (j = b_k; j < 81; j++) {
              sum1_tmp = b_k + 80 * (j - 1);
              ab = (h[sum1_tmp - 1] + s * h[sum1_tmp]) + rt1r * h[sum1_tmp + 1];
              h[sum1_tmp - 1] -= ab * bb;
              h[sum1_tmp] -= ab * aa;
              h[sum1_tmp + 1] -= ab * tst;
            }
            if (b_k + 3 < b_i + 1) {
              i = b_k + 2;
            } else {
              i = b_i;
            }
            for (j = 0; j <= i; j++) {
              sum1_tmp = j + 80 * (b_k - 1);
              hoffset = j + 80 * b_k;
              nr = j + 80 * (b_k + 1);
              ab = (h[sum1_tmp] + s * h[hoffset]) + rt1r * h[nr];
              h[sum1_tmp] -= ab * bb;
              h[hoffset] -= ab * aa;
              h[nr] -= ab * tst;
            }
            for (j = 0; j < 80; j++) {
              sum1_tmp = j + 80 * (b_k - 1);
              hoffset = j + 80 * b_k;
              nr = j + 80 * (b_k + 1);
              ab = (z[sum1_tmp] + s * z[hoffset]) + rt1r * z[nr];
              z[sum1_tmp] -= ab * bb;
              z[hoffset] -= ab * aa;
              z[nr] -= ab * tst;
            }
          } else if (nr == 2) {
            for (j = b_k; j < 81; j++) {
              sum1_tmp = b_k + 80 * (j - 1);
              tst = h[sum1_tmp - 1];
              ab = tst + s * h[sum1_tmp];
              h[sum1_tmp - 1] = tst - ab * bb;
              h[sum1_tmp] -= ab * aa;
            }
            for (j = 0; j <= b_i; j++) {
              sum1_tmp = j + 80 * (b_k - 1);
              hoffset = j + 80 * b_k;
              ab = h[sum1_tmp] + s * h[hoffset];
              h[sum1_tmp] -= ab * bb;
              h[hoffset] -= ab * aa;
            }
            for (j = 0; j < 80; j++) {
              sum1_tmp = j + 80 * (b_k - 1);
              tst = z[sum1_tmp];
              hoffset = j + 80 * b_k;
              ab = tst + s * z[hoffset];
              z[sum1_tmp] = tst - ab * bb;
              z[hoffset] -= ab * aa;
            }
          }
        }
        its++;
      }
    }
    if (!goto150) {
      info = b_i + 1;
      exitg1 = true;
    } else {
      if ((L != b_i + 1) && (L == b_i)) {
        i = b_i + 80 * b_i;
        s = h[i - 1];
        j = 80 * (b_i - 1);
        hoffset = b_i + j;
        rt1r = h[hoffset];
        tst = h[i];
        xdlanv2(&h[(b_i + 80 * (b_i - 1)) - 1], &s, &rt1r, &tst, &aa, &ab, &bb,
                &ba, &s_tmp, &rt2r);
        h[i - 1] = s;
        h[hoffset] = rt1r;
        h[i] = tst;
        if (80 > b_i + 1) {
          hoffset = 78 - b_i;
          nr = b_i + (b_i + 1) * 80;
          for (k = 0; k <= hoffset; k++) {
            sum1_tmp = nr + k * 80;
            tst = h[sum1_tmp];
            aa = h[sum1_tmp - 1];
            h[sum1_tmp] = s_tmp * tst - rt2r * aa;
            h[sum1_tmp - 1] = s_tmp * aa + rt2r * tst;
          }
        }
        hoffset = b_i * 80;
        if (b_i - 1 >= 1) {
          for (k = 0; k <= b_i - 2; k++) {
            nr = hoffset + k;
            sum1_tmp = j + k;
            tst = s_tmp * h[sum1_tmp] + rt2r * h[nr];
            h[nr] = s_tmp * h[nr] - rt2r * h[sum1_tmp];
            h[sum1_tmp] = tst;
          }
        }
        for (k = 0; k < 80; k++) {
          nr = hoffset + k;
          tst = z[nr];
          sum1_tmp = j + k;
          aa = z[sum1_tmp];
          z[nr] = s_tmp * tst - rt2r * aa;
          z[sum1_tmp] = s_tmp * aa + rt2r * tst;
        }
      }
      b_i = L - 2;
    }
  }
  return info;
}

/* End of code generation (xdhseqr.c) */
