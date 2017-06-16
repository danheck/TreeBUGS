#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP TreeBUGS_betampt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TreeBUGS_loglikMPT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP TreeBUGS_simplempt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"TreeBUGS_betampt",   (DL_FUNC) &TreeBUGS_betampt,   8},
  {"TreeBUGS_loglikMPT", (DL_FUNC) &TreeBUGS_loglikMPT, 6},
  {"TreeBUGS_simplempt", (DL_FUNC) &TreeBUGS_simplempt, 8},
  {NULL, NULL, 0}
};


void R_init_TreeBUGS(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

