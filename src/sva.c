#include <R.h> 
#include <Rinternals.h> 
#include <Rmath.h>
#include <R_ext/Rdynload.h>

SEXP monotone(SEXP lfdr){ 
  int i,n;
  double *vec, *out;
  vec = REAL(lfdr);
  n = length(lfdr);
 
  SEXP Rlfdr; 
  PROTECT(Rlfdr = allocVector(REALSXP,n));
  out = REAL(Rlfdr);
  out[0] = vec[0];
  for (i = 1; i < n; i++){
    if(vec[i] < out[(i-1)]){
      out[i] = out[(i-1)];
    }else{
      out[i] = vec[i];
    }
  }  
  UNPROTECT(1);
  return(Rlfdr);
}
static const
R_CallMethodDef callMethods[]  = {
  {"monotone", (DL_FUNC) &monotone, 1},
  NULL
};

void R_init_sva(DllInfo *info)
{
     R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
     
void R_unload_sva(DllInfo *info)
{
    (void) info;
}
