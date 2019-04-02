

#function calculating running time
runtime=function(fn){
  start_time=Sys.time()
  fn
  end_time=Sys.time()
  return(end_time-start_time)
}


#Method 1: Using definition of correlation  

#Method 1.1: Directly using package cubature to calculate integrals
install.packages("cubature")
library(cubature)
f<-function(x){(sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20)}
fx<-function(x){x[1]*((sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20))}
fy<-function(x){x[2]*((sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20))}
fxx<-function(x){(x[1]^2)*((sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20))}
fyy<-function(x){(x[2]^2)*((sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20))}
fxy<-function(x){x[2]*x[1]*((sin(x[1]*pi)*sin(pi*x[1]^2)^20+sin(x[2]*pi)*sin(2*pi*x[2]^2)^20))}

Corr_def_true=function(){
  start_time=Sys.time()
  c<- as.numeric(hcubature(f,c(0,0),c(1,1))['integral'])
  Ex<- as.numeric(hcubature(fx,c(0,0),c(1,1))['integral'])/c
  Ey<- as.numeric(hcubature(fy,c(0,0),c(1,1))['integral'])/c
  Exx<- as.numeric(hcubature(fxx,c(0,0),c(1,1))['integral'])/c
  Eyy<- as.numeric(hcubature(fyy,c(0,0),c(1,1))['integral'])/c
  Exy<- as.numeric(hcubature(fxy,c(0,0),c(1,1))['integral'])/c
  corr<-(Exy-Ex*Ey)/sqrt((Exx-Ex^2)*(Eyy-Ey^2))
  end_time=Sys.time()
  runtime=end_time-start_time
  result=list("corr"=corr,"runtime"=runtime)
  return(result)
}
corr_true=Corr_def_true()
corr_true$corr
corr_true$runtime
#Result using cubature package: -0.06095707, assume it to be the true answer


#Method 1.2: Using Monte Carlo integration to calculate the integrals 
MC <-function(f,n){
  x<-rep(NA,n);
  i=1;
  r <- runif(n);
  while(i<=n){
    x[i] <- f(r[i]);
    i=i+1;
  }
  result=mean(x);
  return(result);
}

#let f(x,y)= f(x) + f(y)
fx<-function(x){sin(x*pi)*sin(pi*x^2)^20}
fy<-function(y){sin(y*pi)*sin(2*pi*y^2)^20}
xfx<-function(x){x*sin(x*pi)*sin(pi*x^2)^20}
yfy<-function(y){y*sin(y*pi)*sin(2*pi*y^2)^20}
xxfx<-function(x){(x^2)*sin(x*pi)*sin(pi*x^2)^20}
yyfy<-function(y){(y^2)*sin(y*pi)*sin(2*pi*y^2)^20}


Corr_def_Q4=function(N){  
  start_time=Sys.time()
  mfx=MC(fx,N)
  mxfx=MC(xfx,N)
  mxxfx=MC(xxfx,N)
  mfy=MC(fy,N)
  myfy=MC(yfy,N)
  myyfy=MC(yyfy,N)
  
  c=mfx + mfy
  
  #E(XY) = 1/2F(xf(x)) + 1/2F(yf(y))
  #E(X) = F(xf(x)) + 1/2F(f(y))
  #E(X^2) = F(x^2f(x)) + 1/3F(f(y))
  #E(Y) = F(yf(y)) + 1/2F(f(x))
  #E(Y^2) = F(y^2f(y)) + 1/3F(f(x))
  
  Ex=(mxfx + mfy/2)/c
  Ey=(myfy + mfx/2)/c
  Exx=(mxxfx + mfy/3)/c
  Eyy=(myyfy + mfx/3)/c
  Exy=(myfy/2 + mxfx/2)/c
  
  corr=(Exy-Ex*Ey)/sqrt((Exx-Ex^2)*(Eyy-Ey^2))
  end_time=Sys.time()
  runtime=end_time-start_time
  result=list("corr"=corr,"runtime"=runtime)
  return(result)
}

Corr_def_Q4(1e8)

#repeated trails for m=100 times, need 3-4mins
m=100
corrs_def=rep(NA,m)
corrs_def_runtime=rep(NA,m)
start_time=Sys.time()
for(i in 1:m){
  result=Corr_def_Q4(1e5)
  corrs_def[i]=result$corr
  corrs_def_runtime[i]=result$runtime
}
end_time=Sys.time()
def_runtime_100=end_time-start_time


mean(corrs_def)
sd(corrs_def)

LCI_def=mean(corrs_def)-qnorm(0.975,mean=0,sd=1)*sd(corrs_def)
UCI_def=mean(corrs_def)+qnorm(0.975,mean=0,sd=1)*sd(corrs_def)

#Plots for the result
plot(corrs_def,xlab)
abline(h=mean(corrs_def),col="red")
abline(h=LCI_def,col="blue")
abline(h=UCI_def,col="blue")

#save the result
write.csv(corrs_def,file="def_correlation.csv")
write.csv(c(mean(corrs_def_runtime),def_runtime_100),file="def_runtime.csv")






