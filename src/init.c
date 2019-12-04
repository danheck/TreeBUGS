#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _TreeBUGS_betampt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeBUGS_loglikMPT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TreeBUGS_simplempt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_TreeBUGS_betampt",   (DL_FUNC) &_TreeBUGS_betampt,   10},
  {"_TreeBUGS_loglikMPT", (DL_FUNC) &_TreeBUGS_loglikMPT,  6},
  {"_TreeBUGS_simplempt", (DL_FUNC) &_TreeBUGS_simplempt, 10},
  {NULL, NULL, 0}
};

void R_init_TreeBUGS(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
