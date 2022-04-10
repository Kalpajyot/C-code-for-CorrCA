/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * C_corrca_initialize.c
 *
 * Code generation for function 'C_corrca_initialize'
 *
 */

/* Include files */
#include "C_corrca_initialize.h"
#include "C_corrca_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void C_corrca_initialize(void)
{
  rt_InitInfAndNaN();
  isInitialized_C_corrca = true;
}

/* End of code generation (C_corrca_initialize.c) */
