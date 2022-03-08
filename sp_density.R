###generate data of size 1000 (stored in x) from density f
set.seed(1)
mu1 = .2
mu2 = .7
n <- 10000
m <- n*10
z <- rnorm(m, mean=mu1, sd=.1); x <- z[(z>0)&(z<1)]
x <- x[1:n]
z <- rnorm(m, mean=mu2, sd=.2); x2 <- z[(z>0)&(z<1)]
x2 <- x2[1:n]
z <- sample(0:1, size = n, replace=T)
x[z==1] <- x2[z==1]

#### compute the matrix whose (i,j)th element is the integral of B_i*B_j
require("splines")
knotlist <- (1:8)/9
nb <- length(knotlist)+4
M <- matrix(0, nb, nb)
for (i in 1:nb){
  for (j in i:nb){
    tem <- function(u){
      bx <- bs(u, knots=knotlist,
               Boundary.knots = c(0,1), intercept=T)
      return(bx[,i]*bx[,j])
    }
    M[i,j] <- integrate(tem, 0, 1)$value
    if (j > i) { M[j,i] <- M[i,j] }
  }
}

#### compute fhat, the estimator of f using method of moments
moments <- apply(bs(x, knots=knotlist, Boundary.knots=c(0,1), intercept=T),
                 2, mean)
ahat <- solve( M, moments)
fhat <- function(u){
  ans <- bs(u, knots = knotlist,
            Boundary.knots = c(0, 1), intercept=T) %*% ahat
  return( as.numeric(ans) )
}

##### compare fhat with the true density f
k0 = pnorm(1, mean=mu1, sd=.1) - pnorm(0, mean=mu1, sd=.1)
k1 = pnorm(1, mean=mu2, sd=.2) - pnorm(0, mean=mu2, sd=.2)
f <- function(x){
  ans <- 0.5*dnorm(x, mean=mu1, sd=.1)/k0 + 0.5*dnorm(x, mean=mu2, sd=.2)/k1
  ans[x>1] = 0
  ans[x<0] = 0
  return(ans)
}
curve(f, 0, 1)
curve(fhat, 0, 1, add=T, col=2)
## compute ISE
tem <- function(u){ (fhat(u)-f(u))^2 }
integrate(tem, 0, 1)

####### obtain a normalized version
k2 <- integrate(fhat, 0, 1)$value
fhat1 <- function(u){ fhat(u)/k2 }
curve(fhat1, 0, 1, add=T, col=3)

###### Check spline approximation accuracy using the given basis functions
x0 <- (1:1000)/1001
y1 <- f(x0)
bx <- bs(x0, knots=knotlist, Boundary.knots=c(0,1),
         intercept = T)
y1.lm <- lm(y1~bx-1)
lines(x0, y1.lm$fitted, col=4)
f.reg <- function(u){
  bx <- bs(u, knots=knotlist, Boundary.knots = c(0,1),
           intercept = T)
  ans <- bx %*% y1.lm$coefficients
  return(ans[,1])
}
tem <- function(u){ (f.reg(u)-f(u))^2 }
integrate(tem, 0, 1)

