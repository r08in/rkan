#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
#define false 0
#define true 1
typedef int bool; 

SEXP Standardize(SEXP X_,SEXP Y_) {  
  
  //declaration
  int n = nrows(X_);
  int p = ncols(X_);  
  SEXP XX_, s_,YY_;
  PROTECT(XX_ = allocMatrix(REALSXP, n, p));
  PROTECT(YY_ = allocVector(REALSXP, n));
  PROTECT(s_ = allocVector(REALSXP, p));  
  double *X = REAL(X_);
  double *XX = REAL(XX_);
  double *Y=REAL(Y_);
  double *YY=REAL(YY_);
  double *s = REAL(s_);
  // bool hasIntercept=HasIntercept(X,n);
  bool hasIntercept=false;
  //scale each collum of x
  for (int j=0; j<p; j++) 
  {
    s[j] = 0;
    for (int i=0; i<n; i++) 
    {
      s[j] += pow(X[j*n+i], 2);
    }
    s[j] = sqrt(s[j]/n);
    if(s[j]!=0)
    {      
      
      for (int i=0; i<n; i++)
      {
        XX[j*n+i] = X[j*n+i]/s[j];
      }
    }
    else
    {
      for (int i=0; i<n; i++)
      {
        XX[j*n+i] =0;
      }
    }
  }
  
  //centralize y
  if(hasIntercept)
  {
    double sum=0;
    for(int i=0;i<n;i++)
    {
      sum+=Y[i];
    }
    double mean=sum/n;
    for(int i=0;i<n;i++)
    {
      YY[i]=Y[i]-mean;
    }
  }
  else
  {
    for(int i=0;i<n;i++)
      YY[i]=Y[i];
  }
  // Return list
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, XX_);
  SET_VECTOR_ELT(res, 1, YY_);
  SET_VECTOR_ELT(res, 2, s_);
  UNPROTECT(4);
  return(res);
}