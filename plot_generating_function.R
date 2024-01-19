###This function can generate plots according to the simulation setting in the
###paper "The Online Closure Principle"

#Input : mu_A, mu_N, n

#mu_A  : Strength of the alternative. Real number greater than 0.
#mu_N  : Conservativeness of null p-values
#n     : Number of hypotheses. 

#Output: Plot of FWER and power as described in Section 5 of the paper "The Online Closure Principle"

plot_generator=function(n, mu_A, mu_N){

###Gaussian testing problem for exact p-values
m=2000   #Number of Trials
n=n    #Number of Hypotheses per Trial
mu_A=mu_A    #Strength of the alternative
mu_N=mu_N   #Conservativeness of null p-values (<0 for conservative null p-values)
#pi_A is defined in the loop below

###Initialise Hyperparameters
gamma=6/(pi^2*((1:n)^2))
alpha=0.2
tau=0.8
lambda=0.3

#Correlation within batch
corr=0.8

###Set seed to make the results reproducible
set.seed(12345)

#Set the different batch sizes
batch_sizes=c(1,10,25,100)

###Predefine vectors for FWER and power of the different procedures
FWER_Spending=matrix(0,9,length(batch_sizes))
power_Spending=matrix(0,9,length(batch_sizes))
FWER_cSpending=matrix(0,9,length(batch_sizes))
power_cSpending=matrix(0,9,length(batch_sizes))
b=1 #Counter 
for(batch_size in batch_sizes){
###Parameters for a local dependence structure given by batches
batch_number=ceiling(n/batch_size)
lags=rep(seq(0,(batch_size-1),1),batch_number)
lags=lags[1:n]
sigma=matrix(corr,batch_size,batch_size)+diag((1-corr),batch_size)
mu=rep(0,batch_size)

###Generate p-values and compute FWER and power for the desired procedures
for(l in 1:9){
  ##Fast way for p-values in batches.
  pi_A=l/10
  p=matrix(,nrow=n,ncol=m)
  hypo=matrix(,nrow=n,ncol=m)
  for(j in 1:m){
    hypo[,j]=rbinom(n,1,pi_A)
    X=rep(0,batch_number*batch_size)
    for(k in 1:batch_number){
    X[((k-1)*batch_size+1):(k*batch_size)]=mvrnorm(1,mu,sigma)
    }
    X=X[1:n]
    Z=mu_N*(hypo[,j]-1)*(-1)+mu_A*hypo[,j]+X
    p[,j]=pnorm(-Z)
  }

##ADDIS-Spending
V=rep(0,m)      #Indicates, whether there was at least one type 1 error in a trial
power=rep(0,m)  #Power within each trial

for(j in 1:m){
  alpha_ind=ADDIS_Spending(alpha, gamma, tau, lambda,lags, p[,j], n)
  hypo_est=alpha_ind>=p[,j]
  V[j]=max((hypo[,j]==0 & hypo_est==1))
  D=(hypo[,j]==1 & hypo[,j]==hypo_est)
  power[j]=sum(D)/sum(hypo[,j])
}
FWER_Spending[l,b]=mean(V,na.rm=TRUE)
power_Spending[l,b]=mean(power,na.rm=TRUE)

##Closed ADDIS-Spending
V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
power=rep(0,m)  #Power within each trial

for(j in 1:m){
  alpha_ind=closed_ADDIS_Spending(alpha, gamma, tau,lambda,lags, p[,j], n)
  hypo_est=alpha_ind>=p[,j]
  V[j]=max((hypo[,j]==0 & hypo_est==1))
  D=(hypo[,j]==1 & hypo[,j]==hypo_est)
  power[j]=sum(D)/sum(hypo[,j])
}
FWER_cSpending[l,b]=mean(V,na.rm=TRUE)
power_cSpending[l,b]=mean(power,na.rm=TRUE)

}
b=b+1
}
###Create Plot for ADDIS-Spending

results_df=data.frame(seq(0.1,0.9,0.1),power_Spending,FWER_Spending)

p_spending=ggplot(results_df, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "1",colour = "2")) + 
  geom_point(aes(y = X1,colour = "2")) +
  geom_line(aes(y = X2, linetype = "2",colour = "2")) + 
  geom_point(aes(y = X2,colour = "2")) +
  geom_line(aes(y = X3, linetype = "3",colour = "2")) + 
  geom_point(aes(y = X3,colour = "2")) +
  geom_line(aes(y = X4, linetype = "4",colour = "2")) + 
  geom_point(aes(y = X4,colour = "2")) +
  scale_colour_manual(guide="none", values=c("1"="limegreen", "2"="#f84f4f"))+
  scale_linetype_manual(name  ="Batch-size", values = c("1"="solid","2"="longdash","3"="dashed","4"="dotted"), labels=c("1","10","25", "100"))+
  xlab(expression(pi[A]))+
  ylab("Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle("ADDIS-Spending") 

###Create Plot for closed ADDIS-Spending
results_df_cSpending=data.frame(seq(0.1,0.9,0.1), power_cSpending, FWER_cSpending)

p_cSpending=ggplot(results_df_cSpending, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "1",colour = "4")) + 
  geom_point(aes(y = X1,colour = "4")) +
  geom_line(aes(y = X2, linetype = "2",colour = "4")) + 
  geom_point(aes(y = X2,colour = "4")) +
  geom_line(aes(y = X3, linetype = "3",colour = "4")) + 
  geom_point(aes(y = X3,colour = "4")) +
  geom_line(aes(y = X4, linetype = "4",colour = "4")) + 
  geom_point(aes(y = X4,colour = "4")) +
  scale_colour_manual(guide="none", values=c("3"="limegreen", "4"="cornflowerblue"))+
  scale_linetype_manual(name  ="Batch-size", values = c("1"="solid","2"="longdash","3"="dashed","4"="dotted"), labels=c("1","10","25", "100"))+
  xlab(expression(pi[A]))+
  ylab("Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle("Closed ADDIS-Spending")

combined <- p_spending + p_cSpending  & theme(legend.position = "bottom")
return(combined + plot_layout(guides = "collect"))
}


