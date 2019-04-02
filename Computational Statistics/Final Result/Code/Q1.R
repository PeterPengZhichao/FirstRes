#Q1-Bisection & Newton Method

#f(x)=log(x+x^2)/(1+x^3)

#Be aware that the domain of x is (-inf,-1)U(0,+inf)

f=function(x){
  result=log(x+(x^2))/(1+(x^3))
  return(result)
}

#f'(x)=((1/(x+x^2))*(1+2x)(1+x^3)-3x^2*log(x+x^2))/(1+x^3)^2

#find f'(x)=0 equal to find its nominator g(x)= 0#Q1-Bisection & Newton Method
#f(x)=log(x+x^2)/(1+x^3)
#Be aware that the domain of x is (-inf,-1)U(0,+inf)

f=function(x){
  result=log(x+(x^2))/(1+(x^3))
  return(result)
}

#f'(x)=((1/(x+x^2))*(1+2x)(1+x^3)-3x^2*log(x+x^2))/(1+x^3)^2
#find f'(x)=0 equal to find its nominator g(x)= 0
g=function(x){
  result=((1/(x+x^2))*(1+2*x)*(1+(x^3))-3*(x^2)*log(x+(x^2)))/(1+x^3)^2
  return(result)
}

#set f''x nominator h(x)
h=function(x){
  result=(x^2*(x*(6*(2*x^3-1)*log(x*(x+1))-14*x^3+8*x^2-9*x-6)-1)-1)/(x^2*(x+1)^3*(x^2-x+1)^3)
  return(result)
}

#Biseciton Method
bisection=function(fx,a,b,error,step){
  i<-0
  if(f(a)*f(b)>0) return("Invalid Interval [a,b]")
  while(b-a>error){
    c=(a+b)/2
    ifelse(fx(a)*fx(c)<0,b<-c,a<-c);
    i=i+1
    if(i>=step) return("Steps has been used up")
  }
  root=(a+b)/2
  return(root)
}

#Newton's Method
newton=function(fx,dfx,x_n,error,step){
  j<-0
  x_n_1=x_n+1
  while(abs(x_n_1-x_n)>error){
    x_c=x_n-fx(x_n)/dfx(x_n)
    x_n_1=x_n
    x_n=x_c
    j=j+1
    if(j>=step) return("Steps has been used up")
  }
  root=x_n_1
  return(root)
}

#obtain educated guess by ploting
xneg=seq(from=-5,to=-1, by=0.1)
xpos=seq(from=0,to=2, by=0.1)
x=c(xneg,xpos)
plot(x,f(x),type="l")

#Impletement bisection method and Newton's method
b<-bisection(g,0.5,1.5,10e-10,1000000)
b
f(b)
n<-newton(g,h,1,10e-10,1000000)
n
f(n)
#number of iterations
bisection_step=function(fx,a,b,error,step){
  i<-0
  if(f(a)*f(b)>0) return("Invalid Interval [a,b]")
  while(b-a>error){
    c=(a+b)/2
    ifelse(fx(a)*fx(c)<0,b<-c,a<-c);
    i=i+1
    if(i>=step) return("Steps has been used up")
  }
  return(i)
}
bisection_step(g,0.5,1.5,10e-10,1000000)
newton_step=function(fx,dfx,x_n,error,step){
  j<-0
  x_n_1=x_n+1
  while(abs(x_n_1-x_n)>error){
    x_c=x_n-fx(x_n)/dfx(x_n)
    x_n_1=x_n
    x_n=x_c
    j=j+1
    if(j>=step) return("Steps has been used up")}
  return(j)}
newton_step(g,h,1,10e-10,1000000)