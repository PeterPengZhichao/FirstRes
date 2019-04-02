#Method 3: Gibbs Sampling with Rejection Method 

install.packages("plotMCMC")
install.packages("MASS")
install.packages("plot3D")
library(plotMCMC)
library(MASS)
library(plot3D)

#runtime funti
runtime=function(fn){
  start_time=Sys.time()
  fn
  end_time=Sys.time()
  return(end_time-start_time)
}

#Define functions
f=function(x,y){
  #if(x<0|x>1|y<0|y>1) return("x,y out of range")
  result=(sin(pi*x)*sin(pi*x^2)^20+sin(pi*y)*sin(2*pi*y^2)^20)
  return(result)
}

fx=function(x){
  #if(x<0|x>1) return("x,y out of range")
  result=sin(pi*x)*sin(pi*x^2)^20
  return(result)
}

fxx=function(x){
  #if(x<0|x>1) return("x,y out of range")
  result=pi*cos(pi*x)*sin(pi*x^2)^20+20*sin(pi*x^2)^19*cos(pi*x^2)*2*pi*x
  return(result)
}

fxy=function(x,y0){
  #if(x<0|x>1) return("x,y out of range")
  result=(sin(pi*x)*sin(pi*x^2)^20+fy(y0))/(Itgfx+fy(y0))
  return(result)
}


fy=function(y){
  #if(y<0|y>1) return("x,y out of range")
  result=sin(pi*y)*sin(2*pi*y^2)^20
  return(result)
}

fyy=function(y){
  #if(y<0|y>1) return("x,y out of range")
  result=pi*cos(pi*y)*sin(2*pi*y^2)^20+20*sin(2*pi*y^2)^19*cos(2*pi*y^2)*4*pi*y
  return(result)
}

fyx=function(y,x0){
  #if(y<0|y>1) return("x,y out of range")
  result=(sin(pi*y)*sin(2*pi*y^2)^20+fx(x0))/(Itgfy+fx(x0))
  return(result)
}

#Inverse CDF Method, to find the integral of fx and fy
n=1e6
ux=runif(n,0,1)
uy=runif(n,0,1)
Itgfx=mean(fx(ux))
Itgfy=mean(fy(uy))

#Rejection Method
Rejection=function(f,g,a,b,mx,n){
  i=0
  accept=0
  x=rep(NA,n)
  while(i<=n){
    u=runif(1,0,1);
    x_rnd=runif(1,a,b);
    if(u<(f(x_rnd)/mx*g(x_rnd))){
      x[i]=x_rnd;
      i=i+1
    }
    
  }
  return(x)
}

px=function(x){
  a=0
  b=1
  result=1/(b-a)
  return(result)
}

#find maxfx and maxfy
#Biseciton Method
bisection=function(fx,a,b,error,step){
  i=0
  if(fx(a)*fx(b)>0) return("Invalid Interval [a,b]")
  
  while(b-a>error){
    c=(a+b)/2
    ifelse(fx(a)*fx(c)<0,b<-c,a<-c);
    i=i+1
    if(i>=step) return("Steps has been used up")
  }
  root=(a+b)/2
  return(root)
}

max_x=fx(bisection(fxx,0.2,0.8,1e-10,100))
max_y=fy(bisection(fyy,0,0.8,1e-10,100))


#Gibbs with rejection method
Gibbs_RJ=function(fx,fy,x0,n){
  x=rep(NA,n)
  y=rep(NA,n)
  for (i in 1:n){
    fyx=function(y){
      #if(y<0|y>1) return("x,y out of range")
      result=(sin(pi*y)*sin(2*pi*y^2)^20+fx(x0))/(Itgfy+fx(x0))
      return(result)
    }
    mx=(fx(x0)+max_y)/(fx(x0)+Itgfy)
    y[i]=Rejection(fyx,px,0,1,mx,1)
    y0=y[i]
    fxy=function(x){
      #if(x<0|x>1) return("x,y out of range")
      result=(sin(pi*x)*sin(pi*x^2)^20+fy(y0))/(Itgfx+fy(y0))
      return(result)
    }
    mx=(max_x+fy(y0))/(Itgfx+fy(y0))
    x[i]=Rejection(fxy,px,0,1,mx,1)
    x0=x[i]
  }
  result=cbind(x,y)
  return(result)
}

chainxy=Gibbs_RJ(fx,fy,0.5,1e5)
x=chainxy[(1+1e3):1e5,1]
y=chainxy[(1+1e3):1e5,2]
cor(x,y)
z=rep(NA,length(x))
for(i in 1:length(x)){
  z[i]=f(x[i],y[i])
}
scatter3D(x,y,z,colvar=z,pch=19,cex=0.5)
plotAuto(chainxy)

#repeated trails for m=100 times, needs 50 mins
m=100
n=1e5
burn=1e3
corrs_GRJ=rep(NA,m)
xmatrix_GRJ=matrix(rep(NA,m*(n-burn)),nrow=(n-burn),ncol=m)
ymatrix_GRJ=matrix(rep(NA,m*(n-burn)),nrow=(n-burn),ncol=m)
starttime=Sys.time()
for (i in 1:m){
  chainxy=Gibbs_RJ(fx,fy,0.5,n)
  x=chainxy[(1+burn):n,1]
  y=chainxy[(1+burn):n,2]
  xmatrix_GRJ[,i]=x
  ymatrix_GRJ[,i]=y
  corrs_GRJ[i]=cor(x,y)
}
endtime=Sys.time()
runtime_GRJ_100=endtime-starttime
runtime_GRJ_100
runtime_GRJ=runtime(Gibbs_RJ(fx,fy,0.5,n))
runtime_GRJ

mean(corrs_GRJ)
sd(corrs_GRJ)
range(corrs_GRJ)
max(corrs_GRJ)-min(corrs_GRJ)

LCI_GRJ=mean(corrs_GRJ)-qnorm(0.975,mean=0,sd=1)*sd(corrs_GRJ)
UCI_GRJ=mean(corrs_GRJ)+qnorm(0.975,mean=0,sd=1)*sd(corrs_GRJ)
LCI_GRJ
UCI_GRJ
LenCI_GRJ=UCI_GRJ-LCI_GRJ
LenCI_GRJ

#Plots for the result
par(mfrow=c(1,1))
plot(corrs_GRJ)
abline(h=mean(corrs_GRJ),col="red")
abline(h=LCI_GRJ,col="blue")
abline(h=UCI_GRj,col="blue")

#save the result
write.csv(xmatrix_GRJ,file="Gibbs_RJ_xMatrix.csv")
write.csv(ymatrix_GRJ,file="Gibbs_RJ_yMatrix.csv")
write.csv(corrs_GRJ,file="Gibbs_RJ_correlation.csv")
write.csv(c(runtime_GRJ,runtime_GRJ_100),file="Gibbs_RJ_runtime.csv")