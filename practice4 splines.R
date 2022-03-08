#1
m=4;k=3
library(splines)
x=(1:1000)/1001
knotlist=(1:3)/4
ord=m
knot_all=c(rep(0,ord),knotlist,rep(1,ord))
nb=length(knotlist)+ord
plot(x,x,type="n",ylab="")
X5=c();X6=c(0);X7=c();ans=c()
x1=rep(1,1000);x2=function(x){x}
x3=function(x){x^2};x4=function(x){x^3}
x5=function(x){
  if(x>=knotlist[1]){(x-knotlist[1])^3}
  else{return(0)}
}
for (i in 1:1000){
  X5[i]=x5(x[i])
}
x6=function(x){
  if(x>=knotlist[2]){(x-knotlist[2])^3}
  else{return(0)}
}
for (i in 1:1000){
  X6[i]=x6(x[i])
}
x7=function(x){
  if(x>=knotlist[3]){(x-knotlist[3])^3}
  else{return(0)}
}
for (i in 1:1000){
  X7[i]=x7(x[i])
}
X=cbind(x1,x2(x),x3(x),x4(x),X5,X6,X7)
for ( i in 1:nb){
  y = knot_all[i:(i+ord)]
  f1=function(x){
    bx=bs(x,deg=ord-1, knots=knotlist, Boundary.knots=c(0,1), intercept=T)
    return(bx[,i])
  }
  lines(x,f1(x), type="l")
  ans[i]=summary(lm(f1(x)~X-1))$r.squared
}

#2a
f0=function(x){x*sin(20*x)}
n=1000;x=seq(0,1,length=n);y=f0(x)
order=4;ISE=c();ans=c()
for (i in 1:7){
  k=i
  knotlist=(1:k)/(k+1)
  fhat=function(x0){
    bx=bs(x,deg=order-1, knots=knotlist, Boundary.knots=c(0,1), intercept=T)
    z1=lm(y~bx-1);bx0=bs(x0,deg=order-1, knots=knotlist, Boundary.knots=c(0,1), intercept=T)
    return((bx0 %*% z1$coefficients)[,1])
  }
  f2=function(x){(fhat(x)-f0(x))^2}
  ans=integrate(f2,0,1)$value
  ISE[i]=ans
}

#2b
m=1;ISE1=1
f0=function(x){x*sin(20*x)}
n=1000;x=seq(0,1,length=n);y=f0(x)
while(ISE1>=min(ISE)){
  m=m+1
  x1=matrix(0,1000,m+1)
  x1[,1]=rep(1,1000)
  for (i in 2:(m+1)){
    x1[,i]=x^(i-1)
  }
  fhat2=function(x0){
    lm1=lm(y~x1-1)
    x1=matrix(0,length(x0),m+1)
    x1[,1]=rep(1,length(x0))
    for (i in 2:(m+1)){
      x1[,i]=x0^(i-1)
    }
    return((x1 %*% lm1$coefficients)[,1])
  }
    f3=function(x){(fhat2(x)-f0(x))^2}
    ISE1=integrate(f3,0,1)$value
}

#3(a)
f=function(x){x*sin(20*x)}
n=1000;x=seq(0,1,length=n);sig=0.2
knotlist=(1:6)/7;order=4;ans1=c();ans2=c()
for (i in 1:1000){
  set.seed(i)
  e=rnorm(n,0,sig);y=f(x)+e
  fhat=function(x0){
    bx=bs(x,deg=order-1, knots=knotlist, Boundary.knots=c(0,1), intercept=T)
    z1=lm(y~bx-1);bx0=bs(x0,deg=order-1, knots=knotlist, Boundary.knots=c(0,1), intercept=T)
    return((bx0 %*% z1$coefficients)[,1])
  }
  f2=function(x){(fhat(x)-f(x))^2}
  ans1[i]=integrate(f2,0,1)$value
}
#3(b)
x1=matrix(0,1000,12)
x1[,1]=rep(1,1000)
for (i in 2:12){
  x1[,i]=x^(i-1)
}
for (i in 1:1000){
  set.seed(i)
  e=rnorm(n,0,sig);y=f(x)+e
  lm1=lm(y~x1-1);lm2=lm1$coefficients
  fhat2=function(x0){
    x2=matrix(0,length(x0),12)
    x2[,1]=rep(1,length(x0))
    for (i in 2:12){
      x2[,i]=x0^(i-1)
    }
    return((x2 %*% lm2)[,1])
  }
  f3=function(x){(fhat2(x)-f(x))^2}
  ans2[i]=integrate(f3,0,1)$value
}
IMSE1=sum(ans1)/1000;IMSE2=sum(ans2)/1000

#4(a)
f=function(x){x*sin(20*x)}
n=1000;x=seq(0,1,length=n);sig=0.2;order=4
ans=c()
for (i in 1:1000){
  set.seed(i)
  e=rnorm(n,0,sig)
  y=f(x)+e
  rsscv=c()
  for (k in 1:7){
    knotlist=(1:k)/(k+1)
    bx=bs(x,deg=order-1,knots=knotlist,Boundary.knots=c(0,1),intercept=T)
    Z=cbind(rep(1,n),bx)
    mod=lm(y~Z-1)
    hii.v=lm.influence(mod)$hat
    rsscv[k]=sum(mod$resid^2/(1-hii.v)^2)
  }
  k=which.min(rsscv)
  knotlist=(1:k)/(k+1)
  bx=bs(x,deg=order-1, knots=knotlist, Boundary.knots=c(0,1), intercept=T)
  z1=lm(y~bx-1)
  fhat=function(x0){
  bx0=bs(x0,deg=order-1, knots=knotlist, Boundary.knots=c(0,1), intercept=T)
    return((bx0 %*% z1$coefficients)[,1])
  }
  f2=function(x){(fhat(x)-f(x))^2}
  ans[i]=integrate(f2,0,1)$value
}
IMSE3=sum(ans)/1000

#4(b)
ans1=c()
for (i in 1:1000){
  set.seed(i)
  e=rnorm(n,0,sig)
  y=f(x)+e
  rsscv=c()
  for (m in 1:12){
    Z=matrix(0,1000,m+1)
    Z[,1]=rep(1,1000)
    for (j in 2:(m+1)){
      Z[,j]=x^(j-1)
    }
    mod=lm(y~Z-1)
    hii.v=lm.influence(mod)$hat
    rsscv[m]=sum(mod$resid^2/(1-hii.v)^2)
  }
  m=which.min(rsscv)
  Z=matrix(0,1000,m+1)
  Z[,1]=rep(1,1000)
  for (j in 2:(m+1)){
    Z[,j]=x^(j-1)
  }
  lm1=lm(y~Z-1);lm2=lm1$coefficients
  fhat2=function(x0){
    x2=matrix(0,length(x0),m+1)
    x2[,1]=rep(1,length(x0))
    for (i in 2:(m+1)){
      x2[,i]=x0^(i-1)
    }
    return((x2 %*% lm2)[,1])
  }
  f3=function(x){(fhat2(x)-f(x))^2}
  ans1[i]=integrate(f3,0,1)$value
}
IMSE4=sum(ans1)/1000
