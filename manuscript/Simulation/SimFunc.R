################################################################################
#																			   #
#         Simulation functions for all scenarios in the manuscript             #
#                                                                              #
################################################################################

# The following simulation function is for Figure 2, 3, S1, S3, S4, S11, S12, S13 
SimulateData <- function(
    total.size = 100,
    remeasure.size = 10, 
    scale = 2,
    correlation = c('Weak', 'Moderate', 'Strong'),
    true.effect = c('Weak', "Moderate", "Strong"),
    batch.effect = c('Negative', 'Positive')
  ) {
  
  a0s = c(Weak = 0, Moderate=0.5, High=0.8)
  a1s = c(Negative = -0.5, Positive = 0.5)
  r2s = c(Weak = 0.3, Moderate = 0.6, Strong = 0.9)
  
  X = as.numeric(gl(2, total.size/2)) - 1
  Z = cbind(rep(1, total.size), rnorm(total.size))
  b <- c(0, -0.5)
  v1 = scale^2 
  v2 = 1
  Et <- rnorm(total.size, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
  
  a0 = a0s[true.effect]
  a1 = a1s[batch.effect]
  Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
  Z.r.a <- Z[1 : (total.size / 2), ]
  Et.r.a <- Et[1 : (total.size / 2)]
  r2 = r2s[correlation]
  
  Y.r.a <- a1 + Z.r.a %*% b + r2 * sqrt(v2) * Et.r.a/ sqrt(v1) +
    rnorm(total.size/2, sd = sqrt( (1 - r2^2) * v2 ) )
  if (remeasure.size > total.size/2) {
    cat('The remeasured sample size is larger than the total control size! Input a smaller number.')
  }
  ind.r <- 1:remeasure.size
  Y.r = Y.r.a[ind.r]
  return(list(case = X, design = Z, res = Y, remeasure.ind = ind.r, remeasure.res = Y.r))
}
