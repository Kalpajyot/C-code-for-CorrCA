/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhgeqz.h
 *
 * Code generation for function 'xzhgeqz'
 *
 */

#ifndef XZHGEQZ_H
#define XZHGEQZ_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void xzhgeqz(creal_T A[6400], int ilo, int ihi, creal_T Z[6400], int *info,
             creal_T alpha1[80], creal_T beta1[80]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (xzhgeqz.h) */
