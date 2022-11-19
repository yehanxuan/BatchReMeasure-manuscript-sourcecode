### Replications when noise does not follow Gaussian distribution ##
X =  as.numeric(gl(2, n / 2)) - 1
Z <- cbind(rep(1, n), rnorm(n))
b <- c(0, -0.5)

v1 = r1^2
v2 = 1
Z.r.a <- Z[1 : (n / 2), ]

if (NoiseType == "t") {
  Et <- rt(n, df = 6) * ifelse (X == 0, sqrt(v1), sqrt(v2))/sqrt(3/2)
  Et.r.a <- Et[1 : (n / 2)]
  
  Y.r.a <- a1 + Z.r.a %*% b + r2 * sqrt(v2) * Et.r.a/ sqrt(v1) + 
    rt(n/2, df = 6) * sqrt( (1 - r2^2) * v2 )/sqrt(3/2)
} else if (NoiseType == "Gamma") {
  shape = 2
  scale = 1
  Et <- ( rgamma(n, shape=shape, scale=scale) - shape*scale ) * 
    ifelse (X == 0, sqrt(v1), sqrt(v2))/(sqrt(shape)*scale)
  Et.r.a <- Et[1 : (n / 2)]
  Y.r.a <- a1 + Z.r.a %*% b + r2 * sqrt(v2) * Et.r.a/ sqrt(v1) + 
    ( rgamma(n/2, shape=shape, scale=scale) - shape*scale )*sqrt( (1 - r2^2) * v2 )/
    (sqrt(shape)*scale)
}

Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
ind.r <- 1:n1
Y.r = Y.r.a[ind.r]
alpha = 0.05
