#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _BISNR_BISN_missing(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BISNR_BISN_obsv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _BISNR_QUICParameterLearning(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_BISNR_BISN_missing",          (DL_FUNC) &_BISNR_BISN_missing,          8},
  {"_BISNR_BISN_obsv",             (DL_FUNC) &_BISNR_BISN_obsv,             6},
  {"_BISNR_QUICParameterLearning", (DL_FUNC) &_BISNR_QUICParameterLearning, 6},
  {NULL, NULL, 0}
};

void R_init_BISNR(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
