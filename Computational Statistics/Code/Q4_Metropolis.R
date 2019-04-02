
#Method 2 Two-dimensional Metropolis Algorithm

#install.packages("plotMCMC")
#install.packages("MASS")
#install.packages("plot3D")
library(plotMCMC)
library(MASS)
library(plot3D)

#function calculating running time
runtime=function(fn){
  start_time=Sys.time()
  fn
  end_time=Sys.time()
  return(end_time-start_time)
}

runtime=function(fn){
  start_time=Sys.time()
  fn
  end_time=Sys.time()
  return(end_time-start_time)
}
#Q4 MCMC sampling
f=function(x,y){
  #if(x<0|x>1|y<0|y>1) return("x,y out of range")
  result=(sin(pi*x)*sin(pi*x^2)^20+sin(pi*y)*sin(2*pi*y^2)^20)
  return(result)
}

#proposal distribution, bivariate normal
p=function(x,y){
  result=c(-1,-1)
  while(result[1]<0|result[1]>1|result[2]<0|result[2]>1){
    mean=c(x,y)
    sigma=matrix(c(1,0,0,1),2,2)
    result=mvrnorm(n=1,mu=mean,sigma)
  }
  result1=result
  return(result1)
}


Metropolis=function(f,p,x0,y0,n){
  x=rep(NA,n)
  y=rep(NA,n)
  accept=0
  i=1
  while(i<=n){
    xy=p(x0,y0)
    x_c=xy[1]
    y_c=xy[2]
    u=runif(1,0,1)
    a=f(x_c,y_c)/f(x0,y0)
    if(u<=a){
      x[i]=x_c
      y[i]=y_c
      accept=accept+1
    }
    else{
      x[i]=x0
      y[i]=y0
    }
    x0=x[i]
    y0=y[i]
    i=i+1
  }
  chain=cbind(x,y)
  result=list("chain"=chain,"acceptRate"=accept/n)
  return(result)
}


chainxy=Metropolis(f,p,0.25,0.35,1e5)
chainxy$acceptRate
x=chainxy$chain[(1+1e3):1e5,1]
y=chainxy$chain[(1+1e3):1e5,2]
z=rep(NA,length(x))
for(i in 1:length(x)){
  z[i]=f(x[i],y[i])
}
scatter3D(x,y,z,colvar=z,pch=19,cex=0.5)
plotAuto(chainxy$chain)
plotAuto(chainxy$chain,thin=5)
runtime_MH=runtime(Metropolis(f,p,0.5,0.5,1e5))
runtime_MH


#repeated trails for m=100 times, needs 1.5 h
m=100
n=1e5
burn=1e3
acceptanceRates_MH=rep(NA,m)
corrs_MH=rep(NA,m)
xmatrix_MH=matrix(rep(NA,m*(n-burn)),nrow=(n-burn),ncol=m)
ymatrix_MH=matrix(rep(NA,m*(n-burn)),nrow=(n-burn),ncol=m)

starttime=Sys.time()
for( j in 1:m){
  chainxy=Metropolis(f,p,0.5,0.5,n)
  x=chainxy$chain[(1+burn):n,1]
  y=chainxy$chain[(1+burn):n,2]
  xmatrix_MH[,i]=x
  ymatrix_MH[,i]=y
  acceptanceRates_MH[j]=chainxy$acceptRate
  corrs_MH[j]=cor(x,y)
}
endtime=Sys.time()
runtime_MH_100=endtime-starttime
runtime_MH_100
runtime_MH=runtime(Metropolis(f,p,0.5,0.5,n))
runtime_MH

plot(acceptanceRates_MH)
mean(acceptanceRates_MH)

plot(corrs_MH)
mean(corrs_MH)
sd(corrs_MH)
range(corrs_MH)
max(corrs_MH)-min(corrs_MH)
LCI_MH=mean(corrs_Mh)-qnorm(0.975,mean=0,sd=1)*sd(corrs_Mh)
UCI_MH=mean(corrs_MH)+qnorm(0.975,mean=0,sd=1)*sd(corrs_Mh)
LCI_MH
UCI_MH
LenCI_MH=UCI_MH-LCI_MH
LenCI_MH

#Plots for the result
par(mfrow=c(1,1))
plot(corrs_MH)
abline(h=mean(corrs_MH),col="red")
abline(h=LCI_MH,col="blue")
abline(h=UCI_MH,col="blue")

#save the result
write.csv(xmatrix_MH,file="MH_xMatrix.csv")
write.csv(ymatrix_MH,file="MH_yMatrix.csv")
write.csv(corrs_MH,file="MH_correlation.csv")
write.csv(c(runtime_MH,runtime_MH_100),file="MH_runtime.csv")





#thin the chain by 5
MH_xMatrix=read.csv(file="MH_xMatrix.csv") 
MH_yMatrix=read.csv(file="MH_yMatrix.csv")


x=MH_xMatrix[,2:101]
y=MH_yMatrix[,2:101]


chain=cbind(x[,1],y[,1])
colnames(chain)=c("x","y")
par(mfrow=c(2,1))
plotAuto(chain)


x_thin5=x[seq(from=1,to=99000,by=5),]
y_thin5=y[seq(from=1,to=99000,by=5),]
dim(x_thin5)


corrs_MH_thin5=rep(NA,100)
for(i in 1:100){
  corrs_MH_thin5[i]=cor(x_thin5[,i],y_thin5[,i])
}


mean(corrs_MH_thin5)
sd(corrs_MH_thin5)
LCI_MH_thin5=mean(corrs_MH_thin5)-qnorm(0.975,mean=0,sd=1)*sd(corrs_MH_thin5)
UCI_MH_thin5=mean(corrs_MH_thin5)+qnorm(0.975,mean=0,sd=1)*sd(corrs_MH_thin5)

par(mfrow=c(1,1))
plot(corrs_MH_thin5)
abline(h=mean(corrs_MH_thin5),col="red")
abline(h=LCI_MH_thin5,col="blue")
abline(h=UCI_MH_thin5,col="blue")



chain_thin5=cbind(x[,1],y[,1])
colnames(chain_thin5)=c("x","y")

plotAuto(chain)
plotAuto(chain, thin=5)
plotAuto(chain_thin5)

#save the result
write.csv(corrs_MH_thin5,file="correlation_MH_thin5.csv")
write.csv(x_thin5,file="MH_thin5_xMatrix.csv")
write.csv(y_thin5,file="MH_thin5_yMatrix.csv")
