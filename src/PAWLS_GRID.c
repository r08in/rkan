#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define false 0
#define true 1
typedef int bool; 

// get the Kth largest pos
int GetKLargest(int * next, int head, int expt, int k)
{
  int cur=head;
  if(head==expt)
  {
    cur=next[cur];
  }
    
  for(int i=1;i<k;i++)
  {
    cur=next[cur];
    if(cur==expt) 
    {
      cur=next[cur];
    }
  }
  return cur;
}
//put pos before cur
void modifyIndex(int * pre, int * next, int cur, int pos)
{
  next[pre[pos]]=next[pos];
  next[pre[cur]]=pos;
  pre[pos]=pre[cur];
  next[pos]=cur;
  pre[cur]=pos;
}

// maintain index
void maintainIndex(double * beta, int * pre, int * next, int pos, int head)
{
  int cur=head;
  double b=beta[pos];
  int status=0;
  while(next[cur]!=0)
  {
    if(cur==pos)
    {
      cur=next[cur];
      status=1;
    } else
    {
      if (status==0) // not pass pos yet
      {
        if(beta[cur] <= b)
        {
          modifyIndex(pre, next, cur, pos);
          break;
        }else
        {
          cur=next[cur];
        }
      }
      if(status==1) // already pass
      {
        if(beta[cur] <= b)
        {
          break;
        } else
        {
          modifyIndex(pre, next, pos, cur);
          cur=pos;
        }
       
      }
      
    }
  }
}


//calculate the value of ||Xj'Y\n||
double VectorProduct(double *x, double *y)
{
  int n=sizeof(y);
  double val=0;
  for(int i=0;i<n;i++)
  {
    val+=x[i]*y[i];
  }  
  return val;
}

double UpdateBeta(double z,double lambda,double c)
{
  if(z>lambda)
  {
    return (z-lambda)/c;
  }
  else if(z+lambda<0)
  {
    return(z+lambda)/c;
  }
  else
  {
    return 0;
  }
}

SEXP CleanupG(double *r, double *betaPre, double *wPre, double * shift, 
              double *c, double *Lam2, double *Lam1,
              SEXP beta_, SEXP w_, SEXP loss_, SEXP wloss_, SEXP iter_) 
{
  Free(r);
  Free(betaPre);
  Free(wPre);
  Free(shift);
  Free(c);
  Free(Lam2);
  Free(Lam1);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(res, 0, beta_);
  SET_VECTOR_ELT(res, 1, w_);
  SET_VECTOR_ELT(res, 2, wloss_);
  SET_VECTOR_ELT(res, 3, loss_);
  SET_VECTOR_ELT(res, 4, iter_);
  UNPROTECT(6);
  return(res);
}


SEXP PAWLS_GRID( SEXP X_, SEXP Y_, SEXP Penalty1_, SEXP Penalty2_, SEXP Lambda1_, SEXP Lambda2_,
               SEXP Beta0_, SEXP W0_, SEXP Delta_, SEXP MaxIter_, 
               SEXP Intercept_, SEXP StarBeta_, SEXP StarW_ )
{
  //data convert
  double *x=REAL(X_);
  double *y=REAL(Y_);
  const char *penalty2 = CHAR(STRING_ELT(Penalty2_, 0));
  const char *penalty1 = CHAR(STRING_ELT(Penalty1_, 0));
  double *lambda2=REAL(Lambda2_);
  double *lambda1=REAL(Lambda1_);
  double *beta0=REAL(Beta0_);
  double *w0=REAL(W0_);
  double *starBeta=NULL;
  if(StarBeta_!=NULL)
  {
    starBeta=REAL(StarBeta_);
  }
  double *starW=NULL;
  if(StarW_!=NULL)
  {
    starW=REAL(StarW_);
  }
  double delta=REAL(Delta_)[0];
  double maxIter = REAL(MaxIter_)[0];
  int intercept=REAL(Intercept_)[0];
  
  //data declare
  int n=nrows(X_);
  int m=ncols(X_);
  int L2=length(Lambda2_);
  int L1=length(Lambda1_);
  int lstart2=0, lstart1=0;
  
  //data return
  SEXP res_, beta_, w_, loss_, wloss_, iter_;
  PROTECT(beta_ = allocVector(REALSXP, L2*L1*m));
  PROTECT(w_ = allocVector(REALSXP, L2*L1*n));
  PROTECT(loss_ = allocVector(REALSXP, L2*L1));
  PROTECT(wloss_ = allocVector(REALSXP, L2*L1));
  PROTECT(iter_ = allocVector(REALSXP, L2*L1)); 
  double * beta=REAL(beta_);
  double * w=REAL(w_);
  double *loss=REAL(loss_);
  double *wloss=REAL(wloss_);
  double * iter=REAL(iter_);
  
  double *betaPre = Calloc(m, double);
  double *wPre=Calloc(n, double);
  double *r=Calloc(n, double);
  
  //initial
  if(StarBeta_==NULL)
  {
     for(int i=0;i<m;i++)
   {
     betaPre[i]=0;
   }
  }
  else
  {
    for(int i=0;i<m;i++)
   {
     betaPre[i]=starBeta[i];
   }
  }
  if(StarW_==NULL)
  {
    for(int i=0;i<n;i++)
   {
     wPre[i]=1;
   }
  }
  else
  {
    for(int i=0;i<n;i++)
   {
     wPre[i]=starW[i];
   }
  }
  
  for(int i=0;i<L2*L1;i++)
  {
    loss[i]=wloss[i]=iter[i]=0;
  }
  
  if(StarBeta_==NULL)
  {
    for(int i=0;i<n;i++)
    {
      r[i]=y[i];
    }
  }
  else
  {
    for(int i=0;i<n;i++)
   {
     double temp=0;
     for(int j=0;j<m;j++)
     {
       temp+=x[j*n+i]*betaPre[j];
     }
     r[i]=y[i]-temp;
   }     
  }
  
  
  //temp
  double *shift=Calloc(n+m, double);
  double *c=Calloc(m, double);
  double *Lam2=Calloc(n, double);
  double *Lam1=Calloc(m, double);
  
  //interation for each lambda2
  for(int l2=lstart2;l2<L2;l2++)
  {
    //initial
    for(int i=0;i<n+m;i++)
    {
      shift[i]=0;
    }
    for(int i=0;i<m;i++)
    {
      c[i]=0;
    }
    loss[0*n+l2]=VectorProduct(y,y); //initial loss[l2,1]
    
    if(strcmp(penalty2,"log")==0)
    {
      
      for(int i=0;i<n;i++)
      {
        Lam2[i]=sqrt(lambda2[l2]/fabs(log(w0[i]))*n) ;//init sqrt(lambda2/fabs(log(w0))n)
        
        
      }
    }
    else if(strcmp(penalty2,"L1")==0)//L1
    {
      for(int i=0;i<n;i++)
      {
        Lam2[i]=lambda2[l2]/fabs(1-w0[i])*n;//init sqrt(lambda2/fabs(log(w0))n)
      }
      
    }
    
    //iteration for each lambda1
    for(int l1=lstart1;l1<L1;l1++)
    {
      for(int i=0;i<m;i++)
      {
        Lam1[i]=lambda1[l1]/fabs(beta0[i]);
      }
      
      if(intercept==true)
      {
        Lam1[0]=0;
        //Lam1[1]=0;//for test
      }
      
      //iteration for all covariates
      while(iter[l1*L2+l2]<maxIter)
      {
        iter[l1*L2+l2]+=1;
        //calculate coefficient c
        for(int j=0;j<m;j++)
        {
          c[j]=0;
          for(int i=0;i<n;i++)
          {
            c[j]+=(x[j*n+i]*wPre[i])*(x[j*n+i]*wPre[i]);
            
          }
          c[j]=c[j]/n;
        }
        
        //iteration for each beta
        for(int j=0;j<m;j++)
        {
          //(1)calculate zj 
          double zj=0;
          for(int i=0;i<n;i++)
          {
            zj+=x[j*n+i]*wPre[i]*wPre[i]*r[i];
          }
          zj=zj/n+c[j]*betaPre[j];

          //(2)update betaj
          if(strcmp(penalty1,"L1")==0)
          {
            beta[j*L2*L1+l1*L2+l2]=UpdateBeta(zj,Lam1[j],c[j]);
          }
          else if(strcmp(penalty1,"RIDGE")==0)
          {
            beta[j*L2*L1+l1*L2+l2]=zj/(c[j]+Lam1[j]);
          }

          //(3)update r
          shift[j]=beta[j*L2*L1+l1*L2+l2]-betaPre[j];
          for(int i=0;i<n;i++)
          {
            r[i]-=x[j*n+i]*shift[j];
          }
        }

        //update w
        if(strcmp(penalty2,"L1")==0)//L1
        {
          double sqr=0;
          for(int i=0;i<n;i++)
          {
            sqr=r[i]*r[i];
            if(sqr>Lam2[i])
            {
              w[i*L2*L1+l1*L2+l2]=Lam2[i]/sqr;
            }
            else 
            {
              w[i*L2*L1+l1*L2+l2]=1;
            }
          }
        }
        
        for(int i=0;i<n;i++)
        {
          shift[m+i]=w[i*L2*L1+l1*L2+l2]-wPre[i];
        }
        
        //update betaPre and wPre for next iteration
        for(int i=0;i<m;i++)
        {
          betaPre[i]=beta[i*L2*L1+l1*L2+l2];
        }
        for(int i=0;i<n;i++)
        {
          wPre[i]=w[i*L2*L1+l1*L2+l2];
        }

        //Check for convergence
        if(VectorProduct(shift,shift)<delta)
        {
          break;
        }
        
      }//end for the inner loop
      
      //compute square of loss
      loss[l1*L2+l2]=VectorProduct(r,r);
      double temp=0;
      for(int i=0;i<n;i++)
      {
        temp+=r[i]*wPre[i]*r[i]*wPre[i];
      }
      wloss[l1*L2+l2]=temp;
      
    }//end iteration for each lambda1 fixed lambda2
    
  }//end iteration for each lambda2
  
  //clean and return
  res_=CleanupG(r, betaPre, wPre, shift, c, Lam2, Lam1,
            beta_, w_, loss_, wloss_, iter_);          
  return res_;       
}

