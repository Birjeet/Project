library(MASS)
library(binhf)
#Function for mean vector generator 
Mean_cal <- function(p,a){
  mu=rep(a,times=p)  
  return(mu)
}
#Function to create correlation matrix
Variance_cal <- function(p,r){
  r_i=c()
  for (i in 1:p){
    r_i=c(r_i,r^(i-1))
  }
  s=r_i
  for(i in 1:(p-1)){
    s=rbind(s,shift(r_i, places = i))
  }
  return(s)
}
#Simulating with multivariate normal
n=100
p=2
Sample_1=mvrnorm(n, mu=Mean_cal(p,+1.4246), Sigma = Variance_cal(p,0.5))
Sample_0=mvrnorm(n, mu=Mean_cal(p,-1.4246), Sigma = Variance_cal(p,0.5))



# 1st Process for data mislabeling
sum_vec=function(Smpl,p){
  sum_=c()
  for (i in 1:n){
    sum_=c(sum_,sum(Smpl[i,]))
  }
  return(sum_/sqrt(p)) 
}
sqrt(p)
p_exp_o=function(d_i,a,b){
  func=1/(1+exp(-a-exp(b)*d_i))
  return(func)
}
generate_out=function(prob){
  success<-rbinom(1,size=1,prob=prob)
  c(success)
}


#Probabilities for mis_classifications
p_1=data.frame(p_exp_o(abs(sum_vec(Sample_1,p)),-5,0.45))
p_0=data.frame(p_exp_o(abs(sum_vec(Sample_0,p)),-6,0.45))

p_vec_1=c(t(apply(p_1,1,generate_out)))
p_vec_0=c(t(apply(p_0,1,generate_out)))

Data_1<- data.frame(
  Xp=Sample_1,
  p_out=p_vec_1
)
Data_0<- data.frame(
  Xp=Sample_0,
  p_out=p_vec_0
)

Data_n_1=Data_1[order(Data_1$p_out), ]
Data_n_0=Data_0[order(Data_0$p_out), ]

a=rep(1,times=100)
b=rep(0,times=100)

Data_n_1=cbind(Data_n_1[,-3],z=a,y=b)
Data_n_0=cbind(Data_n_0[,-3],z=b,y=a)
Data=rbind(Data_n_0,Data_n_1)
Data_clean=Data[,-4]


model1=glm(z~ (Xp.s+Xp.V2),family = 'binomial' ,data=Data)
summary(model1)

Data_validation=rbind(Data_n_0[86:100,],Data_n_1[91:100,])
Data_nonvalidation=rbind(Data_n_0[1:85,],Data_n_1[1:90,])

Data_0_=Data_validation
Data_1_=Data_nonvalidation

log_i<- function(theta){
  B_0=theta[1]
  B1=theta[2]
  B2=theta[3]
  a0=theta[4]
  b0=theta[5]
  a1=theta[6]
  b1=theta[7]
  sum_data_0=0
  for (i in 1:nrow(Data_0_)){
    c=2^(1/2)
    s_f_1=1/(1+exp(-B_0-(B1*Data_0_[i,]$Xp.s+B2*Data_0_[i,]$Xp.V2)))
    s_f_0=1/(1+exp(+B_0+(B1*Data_0_[i,]$Xp.s+B2*Data_0_[i,]$Xp.V2)))
    g_i_0=1/(1+exp(-a0-exp(b0)*(B_0+B1*Data_0_[i,]$Xp.s+B2*Data_0_[i,]$Xp.V2)/c))
    g_i_1=1/(1+exp(-a1-exp(b1)*(B_0+B1*Data_0_[i,]$Xp.s+B2*Data_0_[i,]$Xp.V2)/c))
    sum_data_0=sum_data_0+Data_0_[i,]$z*log(s_f_1)+(1-Data_0_[i,]$z)*log(s_f_0)+Data_0_[i,]$y*Data_0_[i,]$z*log(1-g_i_1)+(1-Data_0_[i,]$y)*Data_0_[i,]$z*log(g_i_1)+(1-Data_0_[i,]$z)*Data_0_[i,]$y*log(g_i_0)+(1-Data_0_[i,]$y)*(1-Data_0_[i,]$z)*log(g_i_0)
  }
  sum_data_1=0
  for (i in 1:nrow(Data_1_)){
    c=2^(1/2)
    s_f_1=1/(1+exp(-B_0-(B1*Data_1_[i,]$Xp.s+B2*Data_1_[i,]$Xp.V2)))
    s_f_0=1/(1+exp(+B_0+(B1*Data_1_[i,]$Xp.s+B2*Data_1_[i,]$Xp.V2)))
    g_i_0=1/(1+exp(-a0-exp(b0)*(B_0+B1*Data_1_[i,]$Xp.s+B2*Data_1_[i,]$Xp.V2)/c))
    g_i_1=1/(1+exp(-a1-exp(b1)*(B_0+B1*Data_1_[i,]$Xp.s+B2*Data_1_[i,]$Xp.V2)/c))
    sum_data_1=sum_data_1+log((g_i_0^Data_1_[i,]$y)*(1-g_i_0)^(1-Data_1_[i,]$y)*s_f_0+(g_i_1^(1-Data_1_[i,]$y))*(1-g_i_1)^(Data_1_[i,]$y)*s_f_1)
  }
  loglikelihood=sum_data_0+sum_data_1
  return(-loglikelihood)
  #sum(Data_0_$z*log(s_f_1)+(1-Data_0_$z)*log(s_f_0)+Data_0_$y*Data_0_$z*log(1-g_i_1)+(1-Data_0_$y)*Data_0_$z*log(g_i_1)+(1-Data_0_$z)*Data_0_$y*log(g_i_0)+(1-Data_0_$y)*(1-Data_0_$z)*log(g_i_0))+sum(log((g_i_0_^Data_1_$y)*(1-g_i_0_)^(1-Data_1_$y)*s_f_0_+(g_i_1_^(1-Data_1_$y))*(1-g_i_1_)^(Data_1_$y)*s_f_1_))
}
start_vals <- c(0,0,0,0,0,0,0)
result <- optim(par = start_vals, fn = log_i, method = "BFGS")

Data_fin_1=rbind(Data_n_1[1:90,],Data_n_0[86:100,])
Data_fin_0=rbind(Data_n_0[1:85,],Data_n_1[91:100,])
p=rep(1,times=105)
q=rep(0,times=95)
Data_fin_1=cbind(Data_fin_1,p=p)
Data_fin_0=cbind(Data_fin_0,p=q)
Data_noised=rbind(Data_fin_0,Data_fin_1)
model2=glm(p~ Xp.s+Xp.V2,family = 'binomial' ,data=Data_noised)
summary(model2)
func_ldf=function(vector){
  if (sum(vector)>0){
    return(1)
  }else{
    return(0)
  }
}

B_0=result$par[1]
B_1=result$par[2]
B_2=result$par[3]

func_our=function(vector){
  if (B_0+B_1*vector[1]+B_2*vector[2]>0){
    return(1)
  }else{
    return(0)
  }
}
z_ldf=apply(Data_noised[,c(1,2)],1,func_ldf)
z_our=apply(Data_noised[,c(1,2)],1,func_our)
result


