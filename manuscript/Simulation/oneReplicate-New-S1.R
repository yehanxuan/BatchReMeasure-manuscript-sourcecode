#### oneReplicate-new #### 
X =  as.numeric(gl(2, n / 2)) - 1
Z <- cbind(rep(1, n), rnorm(n))
b <- c(0, -0.5)
v1 = r1^2
v2 = 1
Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
Z.r.a <- Z[1 : (n / 2), ]
Et.r.a <- Et[1 : (n / 2)]
Y.r.a <- a1 + Z.r.a %*% b + r2 * sqrt(v2) * Et.r.a/ sqrt(v1) +
  rnorm(n/2, sd = sqrt( (1 - r2^2) * v2 ) )

ind.r <- 1:n1
Y.r = Y.r.a[ind.r]

alpha = 0.05



