#### oneReplicate-new #### 

X =  as.numeric(gl(2, n / 2)) - 1
Z <- cbind(rep(1, n), rnorm(n))
b <- c(0, -0.5)

# Old version of data generating 
#vb = 1/(1 + r2)
#v2 = r2 * vb
#v1 = r1 * r2 * vb
v1 = 2/(1+r1)
v2 = r1 * v1

#Eb <- rnorm(n, sd = sqrt(vb))
Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
#Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Eb + Et
Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et


Z.r.a <- Z[1 : (n / 2), ]
#Eb.r.a <- Eb[1 : (n / 2)]
#Et.r.a <- rnorm(n / 2, sd =  sqrt(v2))
Et.r.a <- Et[1 : (n / 2)]
Y.r.a <- a1 + Z.r.a %*% b + r2 * sqrt(v2) * Et.r.a/ sqrt(v1) +
  rnorm(n/2, sd = sqrt( (1 - r2^2) * v2 ) )

ind.r <- 1:n1
Y.r = Y.r.a[ind.r]



