###This R-file contains the procedures that were used for simulations in the paper "The Online Closure Principle"

#Input: alpha, gamma, tau, lambda,lags, p, n. 

#alpha: Overall significance level. Number between 0 and 1.
#gamma: Weights for Alpha-Spending. Non-negative n-dimensional vector with sum less or equal than 1.
#w:     Edge weights used for Online-Graph. Non-negative n x n upper triangular matrix with row sum less or equal than 1.
#tau:   Used for discarding procedures. n-dimensional vector with values between 0 and 1. Besides, it can be chosen as a fixed number between 0 and 1.
#lambda:Used for the adaptive procedures. n-dimensional vector with values between 0 and tau_i. Besides, it can be chosen as a fixed number fulfilling these conditions.
#lags:  Representing the given local dependence structure. n-dimensional vector of natural numbers (including 0) or a fix natural number (including 0).
#p:     n-dimensional vector of p-values. 
#n:     Number of hypotheses.

#ADDIS-Spending under local dependence (Tian and Ramdas, 2021)
ADDIS_Spending=function(alpha, gamma, tau, lambda,lags, p, n){
  if(length(lags)== 1){lags=rep(lags,n)}
  if(length(tau)== 1){tau=rep(tau,n)}
  if(length(lambda)== 1){lambda=rep(lambda,n)}
  if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n| length(lags)!=n){
    warning("mismatching length")
  }
  S_C=rep(0,n)
  t=rep(1,n)
  alpha_ind=rep(0,n)
  for(i in 1:n){
    if(lags[i]>=i-1){
      t[i]=i
    }else{
      t[i]=1+lags[i]+sum(S_C[1:(i-lags[i]-1)])
    }
    alpha_ind[i]=alpha*gamma[t[i]]*(tau[i]-lambda[i])
    if(p[i]<=tau[i] & p[i]>lambda[i]){
      S_C[i]=1
    }
  }
  return(alpha_ind)
}

#Closed ADDIS-Spending under local dependence 
closed_ADDIS_Spending=function(alpha, gamma, tau, lambda,lags, p, n){
  if(length(lags)== 1){lags=rep(lags,n)}
  if(length(tau)== 1){tau=rep(tau,n)}
  if(length(lambda)== 1){lambda=rep(lambda,n)}
  if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n| length(lags)!=n){
    warning("mismatching length")
  }
  S_C=rep(0,n)
  R=rep(0,n)
  t=rep(1,n)
  alpha_ind=rep(0,n)
  for(i in 1:n){
    if(lags[i]>=1){
    if(lags[i]>=i-1){
      t[i]=1+sum(1-(R[(i-lags[i]):(i-1)]))
    }else{
      t[i]=1+sum(S_C[1:(i-lags[i]-1)])+sum(1-(R[(i-lags[i]):(i-1)]))
    }
    }else{
      t[i]=1+sum(S_C[1:(i-1)])
    }
    alpha_ind[i]=alpha*gamma[t[i]]*(tau[i]-lambda[i])
    if(p[i]<=tau[i] & p[i]>lambda[i]){
      S_C[i]=1
    }
    if(p[i]<=alpha_ind[i]){
      R[i]=1
    }
  }
  return(alpha_ind)
}


#Online Graph
Online_Graph=function(alpha, gamma, w, p, n){
  if(length(gamma)!=n | length(p)!=n){
    warning("mismatching length")
  }
  alpha_ind=rep(0,n)
  alpha_ind[1]=alpha*gamma[1]
  for(i in 2:n){
    alpha_ind[i]=alpha*gamma[i]+sum((alpha_ind[1:(i-1)]>=p[1:(i-1)])*alpha_ind[1:(i-1)]*w[1:(i-1),i])
  }
  return(alpha_ind)
}


#Alpha-Spending
alpha_Spending=function(alpha, gamma, p, n){
  if(length(gamma)!=n | length(p)!=n){
    warning("mismatching length")
  }
  alpha_ind=rep(0,n)
  for(i in 1:n){
    alpha_ind[i]=alpha*gamma[i]
  }
  return(alpha_ind)
} 

#closed Alpha-Spending
closed_Alpha_Spending=function(alpha, gamma, p, n){
  if(length(gamma)!=n | length(p)!=n){
    warning("mismatching length")
  }
  t=1
  alpha_ind=rep(0,n)
  for(i in 1:n){
    alpha_ind[i]=alpha*gamma[t]
    if(p[i]>alpha_ind[i]){
      t=t+1
    }
  }
  return(alpha_ind)
}
