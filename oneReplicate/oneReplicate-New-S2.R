X =  as.numeric(gl(2, n / 2)) - 1
Z <- cbind(rep(1, n), rnorm(n))
b <- c(0, -0.5)

v1 = 2/(1+r1)
v2 = r1 * v1 
v3 = r1 * v2 


Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et

Z.r1.a <- Z[1 : (n / 2), ]
Et.r1.a <- Et[1 : (n / 2)]
Y.r1.a <- a3 + Z.r1.a %*% b + r2 * sqrt(v3) * Et.r1.a/sqrt(v1) +
  rnorm(n/2, sd = sqrt( (1 - r2^2) * v3 ) )


Z.r2.a <- Z[(n/2+1) : n, , drop = F]
Et.r2.a <- Et[(n/2+1) : n]
Y.r2.a <- a0 + a3 + Z.r2.a %*% b + r3 * sqrt(v3) * Et.r2.a/sqrt(v2) + 
  rnorm(n/2, sd = sqrt( (1 - r3^2)*v3 ) )


ind.r1 <- 1:n1 
Y.r1 = Y.r1.a[ind.r1]

ind.r2 <- 1:n2 
Y.r2 = Y.r2.a[ind.r2]
