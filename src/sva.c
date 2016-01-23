#include <R.h> 
#include <Rinternals.h> 
#include <Rmath.h>

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
