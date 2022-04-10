/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzggev.c
 *
 * Code generation for function 'xzggev'
 *
 */

/* Include files */
#include "xzggev.h"
#include "C_corrca_rtwutil.h"
#include "rt_nonfinite.h"
#include "xzhgeqz.h"
#include "xzlartg.h"
#include "xztgevc.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void xzggev(creal_T A[6400], int *info, creal_T alpha1[80], creal_T beta1[80],
            creal_T V[6400])
{
  creal_T atmp;
  double a;
  double absxk;
  double anrm;
  double anrmto;
  double cto1;
  double ctoc;
  double stemp_im;
  int rscale[80];
  int A_tmp;
  int exitg2;
  int exitg3;
  int i;
  int ihi;
  int ii;
  int ilo;
  int j;
  int jcol;
  int jcolp1;
  int jrow;
  int nzcount;
  signed char b_I[6400];
  boolean_T exitg1;
  boolean_T exitg4;
  boolean_T guard1 = false;
  boolean_T ilascl;
  boolean_T notdone;
  *info = 0;
  anrm = 0.0;
  jcol = 0;
  exitg1 = false;
  while ((!exitg1) && (jcol < 6400)) {
    absxk = rt_hypotd_snf(A[jcol].re, A[jcol].im);
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }
      jcol++;
    }
  }
  if (rtIsInf(anrm) || rtIsNaN(anrm)) {
    for (i = 0; i < 80; i++) {
      alpha1[i].re = rtNaN;
      alpha1[i].im = 0.0;
      beta1[i].re = rtNaN;
      beta1[i].im = 0.0;
    }
    for (jrow = 0; jrow < 6400; jrow++) {
      V[jrow].re = rtNaN;
      V[jrow].im = 0.0;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    guard1 = false;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = true;
      guard1 = true;
    } else if (anrm > 1.4885657073574029E+138) {
      anrmto = 1.4885657073574029E+138;
      ilascl = true;
      guard1 = true;
    }
    if (guard1) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = true;
      while (notdone) {
        stemp_im = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((stemp_im > ctoc) && (ctoc != 0.0)) {
          a = 2.0041683600089728E-292;
          absxk = stemp_im;
        } else if (cto1 > absxk) {
          a = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          a = ctoc / absxk;
          notdone = false;
        }
        for (jrow = 0; jrow < 6400; jrow++) {
          A[jrow].re *= a;
          A[jrow].im *= a;
        }
      }
    }
    for (i = 0; i < 80; i++) {
      rscale[i] = 1;
    }
    ilo = 1;
    ihi = 80;
    do {
      exitg3 = 0;
      i = 0;
      j = 0;
      notdone = false;
      ii = ihi;
      exitg1 = false;
      while ((!exitg1) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = ihi;
        jcol = 0;
        exitg4 = false;
        while ((!exitg4) && (jcol <= ihi - 1)) {
          A_tmp = (ii + 80 * jcol) - 1;
          if ((A[A_tmp].re != 0.0) || (A[A_tmp].im != 0.0) ||
              (ii == jcol + 1)) {
            if (nzcount == 0) {
              j = jcol + 1;
              nzcount = 1;
              jcol++;
            } else {
              nzcount = 2;
              exitg4 = true;
            }
          } else {
            jcol++;
          }
        }
        if (nzcount < 2) {
          notdone = true;
          exitg1 = true;
        } else {
          ii--;
        }
      }
      if (!notdone) {
        exitg3 = 2;
      } else {
        if (i != ihi) {
          for (jcol = 0; jcol < 80; jcol++) {
            jcolp1 = (i + 80 * jcol) - 1;
            atmp = A[jcolp1];
            jrow = (ihi + 80 * jcol) - 1;
            A[jcolp1] = A[jrow];
            A[jrow] = atmp;
          }
        }
        if (j != ihi) {
          for (jcol = 0; jcol < ihi; jcol++) {
            jcolp1 = jcol + 80 * (j - 1);
            atmp = A[jcolp1];
            jrow = jcol + 80 * (ihi - 1);
            A[jcolp1] = A[jrow];
            A[jrow] = atmp;
          }
        }
        rscale[ihi - 1] = j;
        ihi--;
        if (ihi == 1) {
          rscale[0] = 1;
          exitg3 = 1;
        }
      }
    } while (exitg3 == 0);
    if (exitg3 != 1) {
      do {
        exitg2 = 0;
        i = 0;
        j = 0;
        notdone = false;
        jcol = ilo;
        exitg1 = false;
        while ((!exitg1) && (jcol <= ihi)) {
          nzcount = 0;
          i = ihi;
          j = jcol;
          ii = ilo;
          exitg4 = false;
          while ((!exitg4) && (ii <= ihi)) {
            A_tmp = (ii + 80 * (jcol - 1)) - 1;
            if ((A[A_tmp].re != 0.0) || (A[A_tmp].im != 0.0) || (ii == jcol)) {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                ii++;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              ii++;
            }
          }
          if (nzcount < 2) {
            notdone = true;
            exitg1 = true;
          } else {
            jcol++;
          }
        }
        if (!notdone) {
          exitg2 = 1;
        } else {
          if (i != ilo) {
            for (jcol = ilo; jcol < 81; jcol++) {
              jcolp1 = 80 * (jcol - 1);
              ii = (i + jcolp1) - 1;
              atmp = A[ii];
              jrow = (ilo + jcolp1) - 1;
              A[ii] = A[jrow];
              A[jrow] = atmp;
            }
          }
          if (j != ilo) {
            for (jcol = 0; jcol < ihi; jcol++) {
              jcolp1 = jcol + 80 * (j - 1);
              atmp = A[jcolp1];
              jrow = jcol + 80 * (ilo - 1);
              A[jcolp1] = A[jrow];
              A[jrow] = atmp;
            }
          }
          rscale[ilo - 1] = j;
          ilo++;
          if (ilo == ihi) {
            rscale[ilo - 1] = ilo;
            exitg2 = 1;
          }
        }
      } while (exitg2 == 0);
    }
    memset(&b_I[0], 0, 6400U * sizeof(signed char));
    for (jcol = 0; jcol < 80; jcol++) {
      b_I[jcol + 80 * jcol] = 1;
    }
    for (jrow = 0; jrow < 6400; jrow++) {
      V[jrow].re = b_I[jrow];
      V[jrow].im = 0.0;
    }
    if (ihi >= ilo + 2) {
      for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
        jcolp1 = jcol + 2;
        for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
          A_tmp = jrow + 80 * jcol;
          xzlartg(A[A_tmp - 1], A[A_tmp], &absxk, &atmp,
                  &A[(jrow + 80 * jcol) - 1]);
          A[A_tmp].re = 0.0;
          A[A_tmp].im = 0.0;
          for (j = jcolp1; j < 81; j++) {
            nzcount = jrow + 80 * (j - 1);
            ctoc = absxk * A[nzcount - 1].re +
                   (atmp.re * A[nzcount].re - atmp.im * A[nzcount].im);
            stemp_im = absxk * A[nzcount - 1].im +
                       (atmp.re * A[nzcount].im + atmp.im * A[nzcount].re);
            cto1 = A[nzcount - 1].im;
            a = A[nzcount - 1].re;
            A[nzcount].re =
                absxk * A[nzcount].re -
                (atmp.re * A[nzcount - 1].re + atmp.im * A[nzcount - 1].im);
            A[nzcount].im =
                absxk * A[nzcount].im - (atmp.re * cto1 - atmp.im * a);
            A[nzcount - 1].re = ctoc;
            A[nzcount - 1].im = stemp_im;
          }
          atmp.re = -atmp.re;
          atmp.im = -atmp.im;
          for (i = 1; i <= ihi; i++) {
            nzcount = (i + 80 * (jrow - 1)) - 1;
            ii = (i + 80 * jrow) - 1;
            ctoc = absxk * A[ii].re +
                   (atmp.re * A[nzcount].re - atmp.im * A[nzcount].im);
            stemp_im = absxk * A[ii].im +
                       (atmp.re * A[nzcount].im + atmp.im * A[nzcount].re);
            cto1 = A[ii].im;
            a = A[ii].re;
            A[nzcount].re = absxk * A[nzcount].re -
                            (atmp.re * A[ii].re + atmp.im * A[ii].im);
            A[nzcount].im =
                absxk * A[nzcount].im - (atmp.re * cto1 - atmp.im * a);
            A[ii].re = ctoc;
            A[ii].im = stemp_im;
          }
          for (i = 0; i < 80; i++) {
            nzcount = i + 80 * (jrow - 1);
            ii = i + 80 * jrow;
            ctoc = absxk * V[ii].re +
                   (atmp.re * V[nzcount].re - atmp.im * V[nzcount].im);
            stemp_im = absxk * V[ii].im +
                       (atmp.re * V[nzcount].im + atmp.im * V[nzcount].re);
            cto1 = V[ii].re;
            V[nzcount].re = absxk * V[nzcount].re -
                            (atmp.re * V[ii].re + atmp.im * V[ii].im);
            V[nzcount].im =
                absxk * V[nzcount].im - (atmp.re * V[ii].im - atmp.im * cto1);
            V[ii].re = ctoc;
            V[ii].im = stemp_im;
          }
        }
      }
    }
    xzhgeqz(A, ilo, ihi, V, info, alpha1, beta1);
    if (*info == 0) {
      xztgevc(A, V);
      if (ilo > 1) {
        for (i = ilo - 2; i + 1 >= 1; i--) {
          jcol = rscale[i] - 1;
          if (rscale[i] != i + 1) {
            for (j = 0; j < 80; j++) {
              jcolp1 = i + 80 * j;
              atmp = V[jcolp1];
              nzcount = jcol + 80 * j;
              V[jcolp1] = V[nzcount];
              V[nzcount] = atmp;
            }
          }
        }
      }
      if (ihi < 80) {
        jrow = ihi + 1;
        for (i = jrow; i < 81; i++) {
          ii = rscale[i - 1];
          if (ii != i) {
            for (j = 0; j < 80; j++) {
              jcolp1 = (i + 80 * j) - 1;
              atmp = V[jcolp1];
              nzcount = (ii + 80 * j) - 1;
              V[jcolp1] = V[nzcount];
              V[nzcount] = atmp;
            }
          }
        }
      }
      for (nzcount = 0; nzcount < 80; nzcount++) {
        absxk = fabs(V[80 * nzcount].re) + fabs(V[80 * nzcount].im);
        for (jcol = 0; jcol < 79; jcol++) {
          ii = (jcol + 80 * nzcount) + 1;
          ctoc = fabs(V[ii].re) + fabs(V[ii].im);
          if (ctoc > absxk) {
            absxk = ctoc;
          }
        }
        if (absxk >= 6.7178761075670888E-139) {
          absxk = 1.0 / absxk;
          for (jcol = 0; jcol < 80; jcol++) {
            jrow = jcol + 80 * nzcount;
            V[jrow].re *= absxk;
            V[jrow].im *= absxk;
          }
        }
      }
      if (ilascl) {
        notdone = true;
        while (notdone) {
          stemp_im = anrmto * 2.0041683600089728E-292;
          cto1 = anrm / 4.9896007738368E+291;
          if ((stemp_im > anrm) && (anrm != 0.0)) {
            a = 2.0041683600089728E-292;
            anrmto = stemp_im;
          } else if (cto1 > anrmto) {
            a = 4.9896007738368E+291;
            anrm = cto1;
          } else {
            a = anrm / anrmto;
            notdone = false;
          }
          for (jrow = 0; jrow < 80; jrow++) {
            alpha1[jrow].re *= a;
            alpha1[jrow].im *= a;
          }
        }
      }
    }
  }
}

/* End of code generation (xzggev.c) */
