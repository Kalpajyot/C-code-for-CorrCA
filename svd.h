/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * svd.h
 *
 * Code generation for function 'svd'
 *
 */

#ifndef SVD_H
#define SVD_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_svd(const double A[6400], double U[6400], double s[80], double V[6400]);

void svd(const double A[6400], double U[80]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (svd.h) */
