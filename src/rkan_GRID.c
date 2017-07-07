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
double  GetKLargest(int * rank, double * beta, int pos, int k, int m)
{
  //printf("m=%d pos=%d\n", m, pos);
  if(m==k)
    return 0;
  if(fabs(beta[rank[k-1]]) > fabs(beta[pos]))
  {
    //printf("beta[rank[k-1]]=%f k=%d\n", beta[rank[k-1]], k);
    return beta[rank[k-1]];
  }

  else 
    return beta[rank[k]];
}

// maintain index
void maintainRank(int * rank, double * beta, int pos,int m)
{
  int status=0;
  int temp1=0, temp2=0;
  //printf("m=%d\n", m);
  for(int i=0;i<m;i++)
  {
    
    //printf("beta[%d]=%f, beta[%d]=%f \n", rank[i], fabs(beta[rank[i]]), pos, fabs(beta[pos]));
    if(rank[i]==pos)
    {
      status=1;
    }
    else
    {
      if(status==0)
      {
        if(fabs(beta[rank[i]])<=fabs(beta[pos]))
        {
          
          temp1=rank[i];
          rank[i]=pos;
          while(temp1!=pos)
          {
            i++;
            temp2=rank[i];
            rank[i]=temp1;
            temp1=temp2;
          }
          break;
        }
      }
      else
      {
        if(fabs(beta[rank[i]])> fabs(beta[pos]))
        {
          //printf("beta[%d]=%f, beta[%d]=%f \n", rank[i], fabs(beta[rank[i]]), pos, fabs(beta[pos]));
          rank[i-1]=rank[i];
          rank[i]=pos;
        }
        else
        {
          break;
        }
      }
      
    }
  }
}

void getRank(int * rank, double *beta,int m)
{
  double * beta2=Calloc(m, double);
  double temp=0;
  int index=0;
  for(int i=0;i<m;i++)
  {
    beta2[i]=beta[i];
    rank[i]=i;
  }
  for(int i=0;i<m;i++)
  {
    index=i;
    for(int j=i+1;j<m;j++)
    {
      if(fabs(beta2[index])<fabs(beta2[j]))
      {
        index=j;
      }
    }
    temp=beta2[index];
    beta2[index]=beta2[i];
    beta2[i]=temp;
    //switch for rank
    temp=rank[index];
    rank[index]=rank[i];
    rank[i]=temp;
  }
  Free(beta2);
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

double Kthreshold(double z, double b, double lambda, double c)
{
  if(z>=lambda+b*c)
  {
    return(z-lambda)/c;
  }
  else if(z>=b*c)
  {
    return(b);
  }
  else if(z>=-b*c)
  {
    return z/c;
  }
  else if(z>=-b*c-lambda)
  {
    return -b;
  }
  else
  {
    return (z+lambda)/c;
  }
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

SEXP CleanupG(double *r, double *betaPre, double *wPre, int* rank, double * shift, 
              double *c, double *Lam2, double *Lam1,
              SEXP beta_, SEXP w_, SEXP loss_, SEXP wloss_, SEXP iter_) 
{
  Free(r);
  Free(betaPre);
  Free(wPre);
  Free(rank);
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


SEXP rkan_GRID( SEXP X_, SEXP Y_, SEXP Lambda1_, SEXP Lambda2_, SEXP k_, SEXP g_, 
                SEXP Beta0_, SEXP W0_, SEXP Delta_, SEXP MaxIter_, 
                //SEXP Intercept_, 
                SEXP StarBeta_, SEXP StarW_ )
{
  //printf("begin\n");
  //data convert
  double *x=REAL(X_);
  double *y=REAL(Y_);
  double *lambda2=REAL(Lambda2_);
  double *lambda1=REAL(Lambda1_);
  double *k=REAL(k_);
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
  //int intercept=REAL(Intercept_)[0];
  
  //data declare
  int n=nrows(X_);
  int m=ncols(X_);
  int L1=length(Lambda1_);
  int L2=length(Lambda2_);
  int L3=length(k_);
  
  //data return
  SEXP res_, beta_, w_, loss_, wloss_, iter_;
  PROTECT(beta_ = allocVector(REALSXP, L3*L2*L1*m));
  PROTECT(w_ = allocVector(REALSXP, L3*L2*L1*n));
  PROTECT(loss_ = allocVector(REALSXP, L3*L2*L1));
  PROTECT(wloss_ = allocVector(REALSXP, L3*L2*L1));
  PROTECT(iter_ = allocVector(REALSXP, L3*L2*L1)); 
  double * beta=REAL(beta_);
  double * w=REAL(w_);
  double *loss=REAL(loss_);
  double *wloss=REAL(wloss_);
  double * iter=REAL(iter_);
  
  int * rank=Calloc(m,int);
  double *betaPre = Calloc(m, double);
  double *wPre=Calloc(n, double);
  double *r=Calloc(n, double);

  double kl=0;
  
  for(int i=0;i<m;i++) //need revised later
  {
    rank[i]=i;
  }
  
  
  for(int i=0;i<L3*L2*L1;i++)
  {
    loss[i]=wloss[i]=iter[i]=0;
  }
  
  //temp
  double *shift=Calloc(n+m, double);
  double *c=Calloc(m, double);
  double *Lam2=Calloc(n, double);
  double *Lam1=Calloc(m, double);
  //initial
  for(int i=0;i<n+m;i++)
  {
    shift[i]=0;
  }
  for(int i=0;i<m;i++)
  {
    c[i]=0;
  }
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
  for (int l3=0; l3<L3;l3++)
  {
    int kth=floor(k[l3]*m);
    //interation for each lambda2
    for(int l2=0;l2<L2;l2++)
    {
      for(int i=0;i<n;i++)
      {
        Lam2[i]=lambda2[l2]/fabs(1-w0[i])*n;//init sqrt(lambda2/fabs(log(w0))n)
      }
      
      //iteration for each lambda1
      for(int l1=0;l1<L1;l1++)
      {
        for(int i=0;i<m;i++)
        {
          Lam1[i]=lambda1[l1]/fabs(beta0[i]);
        }
        
        /*if(intercept==true)
        {
        Lam1[0]=0;
        //Lam1[1]=0;//for test
        }*/
        
        //iteration for all covariates
        while(iter[l3*L1*L2+l2*L1+l1]<maxIter)
        {
          iter[l3*L1*L2+l2*L1+l1]+=1;
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
            kl=GetKLargest(rank, betaPre, j, kth,m); 
            //printf("k=%f klargest=%f \n",k[l3],kl);
            beta[j*L3*L2*L1+l3*L1*L2+l2*L1+l1]=Kthreshold(zj, fabs(kl), Lam1[j], c[j]);
            
            //(3)update r
            shift[j]=beta[j*L3*L2*L1+l3*L1*L2+l2*L1+l1]-betaPre[j];
            betaPre[j]=beta[j*L3*L2*L1+l3*L1*L2+l2*L1+l1];
            maintainRank(rank, betaPre, j,m);
           // getRank(rank,betaPre,m);
           /*
            printf("print beta:");
            for(int i=0;i<m;i++) //test
            {
              printf("%f ",betaPre[i]);
            }
            printf("\n");
            printf("print rank:");
            for(int i=0;i<m;i++) //test
            {
              printf("%d ",rank[i]);
            }
            printf("\n\n");*/
            //break;
              
            for(int i=0;i<n;i++)
            {
              r[i]-=x[j*n+i]*shift[j];
            }
          }
          //break;
          
          //update w
          double sqr=0;
          for(int i=0;i<n;i++)
          {
            sqr=r[i]*r[i];
            if(sqr>Lam2[i])
            {
              w[i*L3*L2*L1+l3*L1*L2+l2*L1+l1]=Lam2[i]/sqr;
            }
            else 
            {
              w[i*L3*L2*L1+l3*L1*L2+l2*L1+l1]=1;
            }
          }
          
          for(int i=0;i<n;i++)
          {
            shift[m+i]=w[i*L3*L2*L1+l3*L1*L2+l2*L1+l1]-wPre[i];
            wPre[i]=w[i*L3*L2*L1+l3*L1*L2+l2*L1+l1];
          }
          
          //Check for convergence
          if(VectorProduct(shift,shift)/(m+n)<delta)
          {
            break;
          }
          
        }//end for the inner loop
        
        //compute square of loss
        loss[l3*L1*L2+l2*L1+l1]=VectorProduct(r,r);
        double temp=0;
        for(int i=0;i<n;i++)
        {
          temp+=r[i]*wPre[i]*r[i]*wPre[i];
        }
        wloss[l3*L1*L2+l2*L1+l1]=temp;
        
      }//end iteration for each lambda1 fixed lambda1
      for(int i=0;i<m;i++)
      {
        betaPre[i]=beta[i*L3*L2*L1+l3*L1*L2+l2*L1+0];
      }
      getRank(rank,betaPre,m);
      for(int i=0;i<n;i++)
      {
        wPre[i]=w[i*L3*L2*L1+l3*L1*L2+l2*L1+0];
      }
      for(int i=0;i<n;i++)
      {
        double temp=0;
        for(int j=0;j<m;j++)
        {
          temp+=x[j*n+i]*betaPre[j];
        }
        r[i]=y[i]-temp;
      } 
      
  }//end iteration for each lambda2
    for(int i=0;i<m;i++)
    {
      betaPre[i]=beta[i*L3*L2*L1+l3*L1*L2+0*L1+0];
    }
    getRank(rank,betaPre,m);
    for(int i=0;i<n;i++)
    {
      wPre[i]=w[i*L3*L2*L1+l3*L1*L2+0*L1+0];
    }
    for(int i=0;i<n;i++)
    {
      double temp=0;
      for(int j=0;j<m;j++)
      {
        temp+=x[j*n+i]*betaPre[j];
      }
      r[i]=y[i]-temp;
    } 
  }//end iteration for each k
 
  
  //clean and return
  res_=CleanupG(r, betaPre, wPre, rank, shift, c, Lam2, Lam1,
            beta_, w_, loss_, wloss_, iter_);          
  return res_;       
}

