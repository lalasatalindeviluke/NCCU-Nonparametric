#1(a)
library(splines)
n=150
set.seed(1)
x=runif(n)
y=sin(2*x)+runif(n,-0.1,0.1)*10
lm1=lm(y~x)
knotlist=0.5
bx=bs(x, knots = knotlist, Boundary.knots = c(0,1), intercept=T)
w.fun=function(x,y,resid=F){
  n=length(y)
  rss0 <- sum((y-lm1$coefficients[1]-lm1$coefficients[2]*x)^2)
  y.lm <- lm(y~bx-1)
  rss <- sum(y.lm$resid^2)
  w <- n*log(rss0/rss)
  if (resid) { ans <- list(w, y.lm$resid); return(ans) } else { return(w) }
}
pvalue1=1-pchisq(w.fun(x,y),3)

#1(b)
set.seed(1)
pv.fun <- function(x,y,m){
  ans <- w.fun(x,y,resid=T)
  w.obs <- ans[[1]]
  resid <- ans[[2]]
  n <- length(y)
  w <- rep(0,m)
  for (j in 1:m){
    e <- sample(resid,n)
    ynew <-lm1$coefficients[1]+lm1$coefficients[2]*x+e
    w[j] <- w.fun(x,ynew)
  }
  return(length(w[w>w.obs])/m)
}
pvalue2=pv.fun(x,y,200) 

#2
pvalue=c()
for (i in 1:500){
  set.seed(i)
  n = 150
  x = runif(n)
  y =1+x+runif(n, -0.1, 0.1)*5
  lm2=lm(y~x)
  knotlist=0.5
  bx=bs(x, knots = knotlist, Boundary.knots = c(0,1), intercept=T)
  w.fun=function(x,y,resid=F){
    n=length(y)
    rss0 <- sum((y-lm2$coefficients[1]-lm2$coefficients[2]*x)^2)
    y.lm <- lm(y~bx-1)
    rss <- sum(y.lm$resid^2)
    w <- n*log(rss0/rss)
    if (resid) { ans <- list(w, y.lm$resid); return(ans) } else { return(w) }
  }
  pvalue[i]=1-pchisq(w.fun(x,y),3)
}
ks.test(pvalue,punif)
#ans for 2a	One-sample Kolmogorov-Smirnov test
#data:  pvalue
#D = 0.060764, p-value = 0.04983
#alternative hypothesis: two-sided
#Reject H0
pvalue1=c()
for (i in 1:500){
  set.seed(i)
  n = 5000
  x = runif(n)
  y =1+x+runif(n, -0.1, 0.1)*5
  lm2=lm(y~x)
  knotlist=0.5
  bx=bs(x, knots = knotlist, Boundary.knots = c(0,1), intercept=T)
  w.fun=function(x,y,resid=F){
    n=length(y)
    rss0 <- sum((y-lm2$coefficients[1]-lm2$coefficients[2]*x)^2)
    y.lm <- lm(y~bx-1)
    rss <- sum(y.lm$resid^2)
    w <- n*log(rss0/rss)
    if (resid) { ans <- list(w, y.lm$resid); return(ans) } else { return(w) }
  }
  pvalue1[i]=1-pchisq(w.fun(x,y),3)
}
ks.test(pvalue1,punif)
#ans for 2b		One-sample Kolmogorov-Smirnov test
#data:  pvalue1
#D = 0.023597, p-value = 0.9435
#alternative hypothesis: two-sided
#Do not reject H0