#include <Rcpp.h>
using namespace Rcpp;

double h(int i, double t);
double psi(int i, double t);

// [[Rcpp::export]]
/*-----Modified K-means clustering algorithm with known cluster number based on MADD------*/
Rcpp::NumericVector gMADD(int s_psi, int s_h, int n_clust, int lb, Rcpp::NumericMatrix M) 
{
  int N = M.nrow(), d = M.ncol();
  int i,j,k,l;
  double d0,min,n0,d00,maxD,minD,WC,PC,DI;
  int nclust,m1,m2,m11,m22;
  int change,iter;
  int B = d/lb;
  
  Rcpp::NumericVector kclust(N), clust_a(N), nc(n_clust), cdist(n_clust);
  Rcpp::NumericMatrix mdist(N,N), dist(N,N), cdist_kD(n_clust, n_clust), cD(n_clust, n_clust);
  Rcpp::NumericMatrix c(N,N), cdist_a(N,N);
  
  for(i=0;i<N-1;i++)
  {
    for(j=i+1;j<N;j++)
    {
      d0=0.0;
      for(k=0;k<B;k++)
      {
        d00 = 0.0;
        for(l=k*lb;l<(k+1)*lb;l++)
        {
          d00=d00+(M(i, l)-M(j, l))*(M(i, l)-M(j, l));
        } 
        d0=d0+psi(s_psi,sqrt(d00)/((double)lb)); 
      }
      dist(i,j)=h(s_h,d0/((double)B));
      dist(j,i)=dist(i,j);
    }
  }
  
  for(i=0;i<N-1;i++)
  {
    for(j=i+1;j<N;j++)
    {
      for(k=0;k<N;k++)
      {
        mdist(i,j)=mdist(i,j)+fabs(dist(i,k)-dist(j,k));
      }
      mdist(i,j)=(mdist(i,j)-2*dist(i,j))/((double) (N-2));
      mdist(j,i)=mdist(i,j);
    }
  }
  
  /* Agglomerative clustering algorithm */
  
  for(i=0;i<N-1;i++)
  {
    for(j=i+1;j<N;j++)
    {
      cdist_a(i,j)=mdist(i,j);
    }
  }
  
  for(i=0;i<N;i++)
  {
    clust_a[i]=i;
    /*initial cluster where all objects are individual clusters*/
  }
  
  for(l=0;l<N-n_clust;l++)  // N-l is the number of clusters to be formed
  {
    nclust=N-l;
    
    /*average linkage method*/
    min=cdist_a(0,1);
    m1=0; m2=1;
    
    for(i=0;i<nclust-1;i++)
    {
      for(j=i+1;j<nclust;j++)
      {
        if(cdist_a(i,j)<min)
        {
          m1=i; m2=j;
          min=cdist_a(i,j);
        }
      }
    }
    
    for(i=0;i<N;i++)
    {
      if(clust_a[i]==m2)
      {
        clust_a[i]=m1;
      }
      else if(clust_a[i]>m2)
      {
        clust_a[i]=clust_a[i]-1;
      }
    }
    
    /*updating distances*/
    for(i=0;i<N-1;i++)
    {
      for(j=i+1;j<N;j++)
      {
        cdist_a(i,j)=0.0;
        c(i,j)=0.0;
      }
    }
    
    for(i=0;i<N-1;i++)
    {
      m1=clust_a[i];
      for(j=i+1;j<N;j++)
      {
        m2=clust_a[j];
        if(m1<m2)
        {
          cdist_a(m1,m2)=cdist_a(m1,m2)+mdist(i,j);
          c(m1,m2)=c(m1,m2)+1.0;
        }
        else if(m1>m2)
        {
          cdist_a(m2,m1)=cdist_a(m2,m1)+mdist(i,j);
          c(m2,m1)=c(m2,m1)+1.0;
        }
      }
    }
    for(i=0;i<nclust-1;i++)
    {
      for(j=i+1;j<nclust;j++)
      {
        cdist_a(i,j)=cdist_a(i,j)/c(i,j);
      }
    }
    
  }
  
  /*K-means clustering algorithm*/
  
  for(i=0;i<N;i++)
  {
    kclust[i]=clust_a[i];
    /*output of average linkage as the initial partition*/
  }
  for(i=0;i<n_clust;i++)
  {
    nc[i]=0.0;
  }
  
  for(i=0;i<N;i++)
  {
    m1=kclust[i];
    nc[m1]=nc[m1]+1.0;
  }
  iter=0;
  change=1;
  while(change>0 && iter<20*N)
  {
    change=0;
    for(i=0;i<N;i++)
    {
      m1=kclust[i];
      n0=nc[m1];
      if(n0>1.0)
      {
        for(k=0;k<n_clust;k++)
        {
          cdist[k]=0.0;
        }
        
        for(j=0;j<N;j++)
        {
          if(j!=i)
          {
            m2=kclust[j];
            cdist[m2]=cdist[m2]+mdist(i,j);
          }
        }
        l=m1; min=cdist[m1]/(nc[m1]-1.0);
        for(k=0;k<n_clust;k++)
        {
          if(k!=m1)
          {
            d0=cdist[k]/nc[k];
            if(d0<min)
            {
              min=d0;
              l=k;
            }
          }
        }
        if(l!=m1)
        {
          nc[m1]=nc[m1]-1.0;
          nc[l]=nc[l]+1.0;
          kclust[i]=l;
          change=change+1;
          iter=iter+1;
        }
      }
    }
  }

  return kclust;              
}

/*-----h functions for MADD based modified K-means clustering algorithm-----*/
double h(int i, double t)
{
  double y;
  if(i==1){ y=sqrt(t); }
  else if(i==2){ y=t; }
  return y;
}
/*-------------------------------------------------------------------------------------------*/


/*-----Family of continuous functions for MADD based modified K-means clustering algorithm-----*/
double psi(int i, double t)
{
  double y;
  if(i==1){ y=pow(t,2.0); }
  else if(i==2){ y=1.0-exp(-t); }
  else if(i==3){ y=1.0-exp(-pow(t,2.0)); }
  else if(i==4){ y=log(1+t); }
  else if(i==5){ y=t; }
  return y;
}
/*-------------------------------------------------------------------------------------------*/


// [[Rcpp::export]]
/*-----Modified K-means clustering algorithm with estimated clust number based on MADD------*/
Rcpp::NumericMatrix gMADD_DI(int s_psi, int s_h, int kmax, int lb, Rcpp::NumericMatrix M) 
{
  int N = M.nrow(), d = M.ncol();
  int i,j,k,k1,l;
  double d0,min,n0,d00,maxD,minD;
  int nclust, n_clust, m1, m2, m11, m22;
  int change,iter;
  int B = d/lb;
  
  Rcpp::NumericVector WC(kmax),PC(kmax),DI(kmax),nc(kmax),cdist(kmax);
  Rcpp::NumericMatrix kclust(kmax,N+1), mdist(N,N), dist(N,N), cdist_kD(kmax,kmax), cD(kmax,kmax); 
  Rcpp::NumericMatrix c(N,N), cdist_a(N,N), clust_a(N,N);
  
  for(i=0;i<N-1;i++)
  {
    for(j=i+1;j<N;j++)
    {
      d0=0.0;
      for(k=0;k<B;k++)
      {
        d00 = 0.0;
        for(l=k*lb;l<(k+1)*lb;l++)
        {
          d00=d00+(M(i, l)-M(j, l))*(M(i, l)-M(j, l));
        } 
        d0=d0+psi(s_psi,sqrt(d00)/((double)lb)); 
      }
      dist(i,j)=h(s_h,d0/((double)B));
      dist(j,i)=dist(i,j);
    }
  }
  
  for(i=0;i<N-1;i++)
  {
    for(j=i+1;j<N;j++)
    {
      for(k=0;k<N;k++)
      {
        mdist(i,j)=mdist(i,j)+fabs(dist(i,k)-dist(j,k));
      }
      mdist(i,j)=(mdist(i,j)-2*dist(i,j))/((double) (N-2));
      mdist(j,i)=mdist(i,j);
    }
  }
  
  /* Agglomerative clustering algorithm */
  
  for(i=0;i<N-1;i++)
  {
    for(j=i+1;j<N;j++)
    {
      cdist_a(i,j)=mdist(i,j);
    }
  }
  
  for(i=0;i<N;i++)
  {
    clust_a(0,i)=i;
    /*initial cluster where all objects are individual clusters*/
  }
  
  for(l=0;l<N-1;l++)  // N-l is the number of clusters to be formed
  {
    nclust=N-l;
    
    /*average linkage method*/
    min=cdist_a(0,1);
    m1=0; m2=1;
    
    for(i=0;i<nclust-1;i++)
    {
      for(j=i+1;j<nclust;j++)
      {
        if(cdist_a(i,j)<min)
        {
          m1=i; m2=j;
          min=cdist_a(i,j);
        }
      }
    }
    
    for(i=0;i<N;i++)
    {
      if(clust_a(l,i)<m2)
      {
        clust_a(l+1,i)=clust_a(l,i);
      }
      else if(clust_a(l,i)==m2)
      {
        clust_a(l+1,i)=m1;
      }
      else
      {
        clust_a(l+1,i)=clust_a(l,i)-1;
      }
    }
    
    /*updating distances*/
    for(i=0;i<N-1;i++)
    {
      for(j=i+1;j<N;j++)
      {
        cdist_a(i,j)=0.0;
        c(i,j)=0.0;
      }
    }
    
    for(i=0;i<N-1;i++)
    {
      m1=clust_a(l+1,i);
      for(j=i+1;j<N;j++)
      {
        m2=clust_a(l+1,j);
        if(m1<m2)
        {
          cdist_a(m1,m2)=cdist_a(m1,m2)+mdist(i,j);
          c(m1,m2)=c(m1,m2)+1.0;
        }
        else if(m1>m2)
        {
          cdist_a(m2,m1)=cdist_a(m2,m1)+mdist(i,j);
          c(m2,m1)=c(m2,m1)+1.0;
        }
      }
    }
    for(i=0;i<nclust-1;i++)
    {
      for(j=i+1;j<nclust;j++)
      {
        cdist_a(i,j)=cdist_a(i,j)/c(i,j);
      }
    }
  }
  
  /*K-means clustering algorithm*/
  
  for(k1=2;k1<=kmax;k1++)
  {
    for(i=0;i<N;i++)
    {
      kclust(k1-1,i)=clust_a(N-k1,i);
      /*output of average linkage as the initial partition*/
    }
    for(i=0;i<k1;i++)
    {
      nc[i]=0.0;
    }
    
    for(i=0;i<N;i++)
    {
      m1=kclust(k1-1,i);
      nc[m1]=nc[m1]+1.0;
    }
    
    iter=0;
    change=1;
    while(change>0 && iter<20*N)
    {
      change=0;
      for(i=0;i<N;i++)
      {
        m1=kclust(k1-1,i);
        n0=nc[m1];
        if(n0>1.0)
        {
          for(k=0;k<k1;k++)
          {
            cdist[k]=0.0;
          }
          
          for(j=0;j<N;j++)
          {
            if(j!=i)
            {
              m2=kclust(k1-1,j);
              cdist[m2]=cdist[m2]+mdist(i,j);
            }
          }
          
          l=m1; min=cdist[m1]/(nc[m1]-1.0);
          for(k=0;k<k1;k++)
          {
            if(k!=m1)
            {
              d0=cdist[k]/nc[k];
              if(d0<min)
              {
                min=d0;
                l=k;
              }
            }
          }
          
          if(l!=m1)
          {
            nc[m1]=nc[m1]-1.0;
            nc[l]=nc[l]+1.0;
            kclust(k1-1,i)=l;
            change=change+1;
            iter=iter+1;
          }
        }
      }
    }
    
  }
  
  
  /* computation of DUNN-index based on k-Means*/
  for(l=1;l<kmax;l++)
  {
    n_clust=l+1;
    
    for(i=0;i<n_clust;i++)
    {
      for(j=i;j<n_clust;j++)
      {
        cdist_kD(i,j)=0.0; 
        cD(i,j)=0.0;
      }
    }
    
    for(i=0;i<N-1;i++)
    {
      m11=kclust(l,i);
      for(j=i+1;j<N;j++)
      {
        m22=kclust(l,j);
        if(m11<=m22)
        {
          cdist_kD(m11,m22)=cdist_kD(m11,m22)+mdist(i,j);
          cD(m11,m22)=cD(m11,m22)+1.0;
        }
        else
        {
          cdist_kD(m22,m11)=cdist_kD(m22,m11)+mdist(i,j);
          cD(m22,m11)=cD(m22,m11)+1.0;
        }
      }
    }
    for(i=0;i<n_clust;i++)
    {
      if(cD(i,i)>0)
      {
        cdist_kD(i,i)=cdist_kD(i,i)/cD(i,i);
      }
      else
      {
        cdist_kD(i,i)=0.0;
      }
    }
    
    for(i=0;i<(n_clust-1);i++)
    {
      for(j=i+1;j<n_clust;j++)
      {
        cdist_kD(i,j)=cdist_kD(i,j)/cD(i,j);
      }
    }
    
    maxD=cdist_kD(0,0);
    for(i=0;i<n_clust;i++)
    {
      if(cdist_kD(i,i)>maxD)
      {
        maxD=cdist_kD(i,i);
      }
    }
    WC[n_clust-1]=maxD;
    
    minD=cdist_kD(0,1);
    for(i=0;i<(n_clust-1);i++)
    {
      for(j=i+1;j<n_clust;j++)
      {
        if(cdist_kD(i,j)<minD)
        {
          minD=cdist_kD(i,j);
        }
      }
    }
    PC[n_clust-1]=minD;
    
  }
  
  for(i=1;i<kmax;i++)
  {
    DI[i]=PC[i]/WC[i];
  }
  
  kclust(_,N) = DI;
  
  return kclust;
}

/*--------------------------------------------------------------------------------------------*/

    
  // [[Rcpp::export]]
/* ------------------------Generating random RxC contigency tables-----------------------------*/
  
Rcpp::NumericMatrix rctab(Rcpp::NumericMatrix M) 
{  
  int nrow = M.nrow(), ncol = M.ncol();
  Rcpp::NumericVector nrowt(nrow), ncolt(ncol);
  bool done1, done2, lsm, lsp;
  int i, ia, iap, ib, ic, id, idp, ie, igp, ihp, ii, iip, j, jc, l, m,  nll, nlm, nlmp, *jwork;
  double r, sumprb, x, y, *fact;
  
  /* calculate rowSums*/
    for(i = 0; i < nrow; i++)
    {
      double row_total = 0;
      for(j = 0; j < ncol; j++)
      {
        row_total += M(i, j);
      }
      nrowt[i] = row_total;
    }
  
  /* calculate colSums*/
    for(j = 0; j < ncol; j++)
    {
      double col_total = 0;
      for (i = 0; i < nrow; i++)
      {
        col_total += M(i, j);
      }
      ncolt[j] = col_total;
    }  
  
  int ntotal = sum(nrowt);
  
  fact = (double*)malloc((ntotal+1)*sizeof(double)); 
  jwork = (int*)malloc(ncol*sizeof(int)); 
  
  for(i=0;i<= ntotal;i++)
  { 
    fact[i] = 0.0;
  }
  
  for(i=0;i< ncol;i++)
  { 
    jwork[i] = 0;
  }
  
  /* Calculate log-factorials.  fact[i] = lgamma(i+1) */
    
    fact[0] = fact[1] = 0.0;
  for(i = 2; i <= ntotal; i++)
    fact[i] = fact[i - 1] + log((double)(i));    
  
  /*  Construct a random matrix.*/
    
    for(i = 0; i < ncol - 1; i++)
    {
      jwork[i] = ncolt[i];
    }
  
  jc = ntotal;
  
  for(l = 0; l < nrow - 1; l++)
  {
    ia = nrowt[l];
    ic = jc;
    jc -= ia;
    
    for(m = 0; m < ncol - 1; m++)
    {
      id = jwork[m];
      ie = ic;
      ic -= id;
      ib = ie - ia;
      ii = ib - id;
      
      /* Test for zero entries in matrix.*/
        
        if(ie == 0)
        {
          ia = 0;
          for(j = m; j < ncol; j++)
          {
            M(l,j) = 0;
          }
          break;
        }
      
      /* Generate a pseudo-random number.*/
        
        r = unif_rand();
        
        /* Compute the conditional expected value of MATRIX(L,M).*/
          
          done1 = false;
          
          for( ; ; )
          {
            nlm = (int)((double)(ia*id)/(double)(ie) + 0.5);
            iap = ia + 1;
            idp = id + 1;
            igp = idp - nlm;
            ihp = iap - nlm;
            nlmp = nlm + 1;
            iip = ii + nlmp;
            x = exp(fact[iap-1] + fact[ib] + fact[ic] + fact[idp-1] - 
                      fact[ie] - fact[nlmp-1] - fact[igp-1] - fact[ihp-1] - fact[iip-1]);
            
            if(r<=x)
            {
              break;
            }
            
            sumprb = x;
            y = x;
            nll = nlm;
            lsp = false;
            lsm = false;
            
            /*Increment entry in row L, column M.*/
              
              while(!lsp)       
              {
                j = (id - nlm)*(ia - nlm);
                
                if(j == 0)
                {
                  lsp = true;
                }
                else
                {
                  nlm = nlm + 1;
                  x = x*(double)(j)/(double)(nlm*(ii + nlm));
                  sumprb += x;
                  
                  if(r<=sumprb)
                  {
                    done1 = true;
                    break;
                  }
                }
                
                done2 = false;
                
                while(!lsm)
                {
                  /* Decrement the entry in row L, column M.*/
                    
                    j = nll*(ii + nll);
                    
                    if(j == 0)
                    {
                      lsm = true;
                      break;
                    }
                    
                    nll = nll - 1;
                    y = y*(double)(j)/(double)((id - nll)*(ia - nll));
                    sumprb += y;
                    
                    if(r<=sumprb)
                    {
                      nlm = nll;
                      done2 = true;
                      break;
                    }
                    
                    if(!lsp)
                    {
                      break;
                    }
                    
                }  
                if(done2)
                {
                  break;
                }
              }
            
            if(done1)
            {
              break;
            }
            
            if(done2)
            {
              break;
            }
            r = sumprb*unif_rand();
          }
          
          M(l,m) = nlm;
          ia -= nlm;
          jwork[m] -= nlm;
    }
    M(l,(ncol-1)) = ia;
  }
  /* Compute the last row.*/
    
    for(j = 0; j < ncol - 1; j++)
    {
      M((nrow-1),j) = jwork[j];
    }
  M((nrow-1),(ncol-1)) = ib - M((nrow-1),(ncol-2));
  return M;
  }
 /*--------------------------------------------------------------------------------------------*/ 
