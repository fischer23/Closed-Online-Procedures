#Used to calculate Table 1 of the paper "The Online Closure Principle"

rm(list=ls())
library(ggplot2)
library(MASS)
library(patchwork)

###load the procedures
source("OCP_Procedures.R")

mydata=read.csv("IMPC_Data.csv")

n=5000

mydata$my_id=seq(1,n)
mydata$lag=rep(0,n)

for(i in 1:n){
  mydata_sameexp=mydata[mydata$Experimental.Id==mydata$Experimental.Id[i],]
  first_id=min(mydata_sameexp$my_id)
  mydata$lag[i]=i-first_id
}


###Initialise Hyperparameters
tau=0.8
lambda=0.16
alpha=0.3

#Setting of lags and p-values
lags=mydata$lag
p=mydata$Sex.P.Val


###Predefine vectors for rejections, number of rejections and individual significance levels for the different procedures

alpha_ind_Spending=matrix(0,ncol=3,nrow=n) 
rej_Spending=matrix(0,ncol=3,nrow=n) 
n_rej_Spending=rep(0,3)

alpha_ind_closed_ADDIS=matrix(0,ncol=3,nrow=n) 
rej_closed_ADDIS=matrix(0,ncol=3,nrow=n)
n_rej_closed_ADDIS=rep(0,3)

alpha_ind_Graph_online=matrix(0,ncol=3,nrow=n)
rej_Graph_online=matrix(0,ncol=3,nrow=n)
n_rej_Graph_online=rep(0,3)

alpha_ind_closed_Alpha=matrix(0,ncol=3,nrow=n)
rej_closed_Alpha=matrix(0,ncol=3,nrow=n)
n_rej_closed_Alpha=rep(0,3)

alpha_ind_Alpha=matrix(0,ncol=3,nrow=n)
rej_Alpha=matrix(0,ncol=3,nrow=n)
n_rej_Alpha=rep(0,3)



for(b in 1:3){
if(b==1){
gamma=1/(3.93*(1:n)^{1.3})        #3.93 was chosen such that sum of (gamma)_{i\in \mathbb{N}} approximately equals 1
}else if(b==2){
gamma=1/(9.33*(1:n)^{1.1})           #9.33 was chosen such that sum of (gamma)_{i\in \mathbb{N}} approximately equals 1
}else{
gamma=1/(13.48*(1:n)^{1.05})       #13.48 was chosen such that sum of (gamma)_{i\in \mathbb{N}} approximately equals 1
}  

  w=abs(matrix(1:n-1 , nrow = n, ncol = n, byrow = TRUE) - (1:n-1))
  w[w==0]=1
  w=matrix(gamma[w],n,n)
  w[upper.tri(w)==0]=0
  
  
  ##ADDIS-Spending
  alpha_ind_Spending[,b]=ADDIS_Spending(alpha, gamma, tau,lambda,lags, p, n)
  rej_Spending[,b]=alpha_ind_Spending[,b]>=p
  n_rej_Spending[b]=sum(rej_Spending[,b])

  ##Closed ADDIS-Spending
  alpha_ind_closed_ADDIS[,b]=closed_ADDIS_Spending(alpha, gamma, tau,lambda,lags, p, n)
  rej_closed_ADDIS[,b]=alpha_ind_closed_ADDIS[,b]>=p
  n_rej_closed_ADDIS[b]=sum(rej_closed_ADDIS[,b])

  ##Online-Graph
  alpha_ind_Graph_online[,b]=Online_Graph(alpha, gamma, w, p, n)
  rej_Graph_online[,b]=alpha_ind_Graph_online[,b]>=p
  n_rej_Graph_online[b]=sum(rej_Graph_online[,b])
  
  #Closed Alpha-Spending
  alpha_ind_closed_Alpha[,b]=closed_Alpha_Spending(alpha, gamma, p, n)
  rej_closed_Alpha[,b]=alpha_ind_closed_Alpha[,b]>=p
  n_rej_closed_Alpha[b]=sum(rej_closed_Alpha[,b])
  
  #Alpha-Spending
  alpha_ind_Alpha[,b]=alpha_Spending(alpha, gamma, p, n)
  rej_Alpha[,b]=alpha_ind_Alpha[,b]>=p
  n_rej_Alpha[b]=sum(rej_Alpha[,b])
  
  
}

