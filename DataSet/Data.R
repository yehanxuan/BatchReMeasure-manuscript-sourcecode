### Generate data ###
# nc1 samplesize of control in batch one 
# nt2 samplesize of case in batch two 
# nc2 samplesize of control in batch two 
# rho correlation 
# p dimension of covariate 
# a0 treatment effect 
# a1 Batch effect
library(MASS)
Generate_S1 = function(nc1, nt2, nc2, p, beta, a0, a1, sigma1, sigma2, rho) {
  if (nc1 < nc2){
    cat("The remeasure size should be smaller than the original control size")
    break()
  }
  
  Zc1 = matrix( rnorm(nc1*p), nc1, p)
  Zt2 = matrix( rnorm(nt2*p), nt2, p)
  # Uniform sampling scheme
  Index  = sample(nc1, nc2, replace = FALSE)
  Zc2 = Zc1[Index, ,drop = F]
  # generate noise 
  #Sigma = matrix(c(sigma1^2, rho*sigma1*sigma2, rho*sigma1*sigma2, sigma2^2), 2, 2)
  #etemp = mvrnorm(n = nc2, mu = c(0,0), Sigma = Sigma)
  #ec1 = c(etemp[, 1])
  ec1 = rnorm(nc1, mean = 0, sd = sigma1)
  Yc1 = ec1 + Zc1%*%beta
  et2 = rnorm(nt2, mean = 0, sd = sigma2)
  Yt2 = a0 + a1 + et2 + Zt2%*%beta
  ecInd = rnorm(nc2, mean = 0, sd = sigma1)
  ec2 = (rho * ec1[Index] + sqrt(1 - rho^2) * ecInd) * sigma2/sigma1
  Yc2 = a1 + ec2 + Zc2%*%beta
  return(list("Yc1" = Yc1, "Yt2" = Yt2, "Yc2" = Yc2, "Zc1" = Zc1, "Zt2" = Zt2, "Zc2" = Zc2, "Index" = Index))
}


Generate_S2 = function(nc1, nt2, nc3, nt3, p, beta, a0, a1, a3, sigma1, sigma2, sigma3, rho1, rho2) {
  if (nc1 < nc3){
    cat("The remeasure size should be smaller than the original control size")
    break()
  } 
  
  if (nt2 < nt3){
    cat("The remeasure size should be smaller than the original treatment size")
    break()
  }
  
  Zc1 = matrix( rnorm(nc1*p), nc1, p)
  Zt2 = matrix( rnorm(nt2*p), nt2, p)
  
  Index_C = sample(nc1, nc3, replace = FALSE)
  Index_T = sample(nt2, nt3, replace = FALSE)
  
  Zc3 = Zc1[Index_C, , drop = F]
  Zt3 = Zt2[Index_T, ,drop = F]
  ec1 = rnorm(nc1, mean = 0, sd = sigma1)
  Yc1 = ec1 + Zc1%*%beta
  et2 = rnorm(nt2, mean = 0, sd = sigma2)
  Yt2 = a0 + a1 + et2 + Zt2%*%beta
  
  ecInd = rnorm(nc3, mean = 0, sd = sigma1)
  ec3 = (rho1 * ec1[Index_C] + sqrt(1 - rho1^2) * ecInd) * sigma3/sigma1
  Zc3 = Zc1[Index_C, , drop = F]
  Yc3 = a3 + Zc3 %*% beta + ec3 
  
  etInd = rnorm(nt3, mean = 0, sd = sigma2)
  et3 = (rho2 * et2[Index_T] + sqrt(1 - rho2^2) * etInd) * sigma3/sigma2
  Zt3 = Zt2[Index_T, ,drop = F]
  Yt3 = a0 + a3 + Zt3 %*% beta + et3
  
  return(list("Yc1" = Yc1, "Yt2" = Yt2, "Yc3" = Yc3, "Yt3" = Yt3, "Zc1" = Zc1, 
              "Zt2" = Zt2, "Zc3" = Zc3, "Zt3" = Zt3, "Index_C" = Index_C, "Index_T" = Index_T))
}




