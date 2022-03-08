knots = c(1/4, 2/4, 3/4)
order = 4
x = (1:1000)/1001

require(splines)
bs_basis = bs(x, deg=order-1, knots=knots, intercept=TRUE, Boundary.knots = c(0,1))


basis_matrix = matrix(1, length(x), 7)
for (i in 2:4){
  basis_matrix[, i] = x^(i-1)
}

f <- function(x, a){
  ans <- (x-a)^3
  ans[x<a] <- 0
  return(ans)
}

for (i in 1:3){
  basis_matrix[, i+4] = f(x, i/4)
}

for (i in 1:7){
  model <- lm(bs_basis[,i]~basis_matrix-1)
  print(summary(model)$r.squared)
}
