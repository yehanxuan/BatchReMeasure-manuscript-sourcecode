### Rewrite BC1 in our computer include Bootstrap ###

Wrap_func = function(part, paras) {
  set.seed(part)
  
  n <- paras$n
  b <- paras$b 
  
  niter <- paras$niter
  
  a0s <- paras$a0s
  a1s <- paras$a1s
  n1s <- paras$n1s
  r1s <- paras$r1s
  r2s <- paras$r2s
  
  resdir <- paras$resdir
  prefix <- paras$prefix
  prefix2 <- paras$prefix2
  
  resdir <- gsub('/$', '', resdir)
  source(file.path(resdir, "ReMeasure_S1.R"))
  sink(file.path(resdir, paste(prefix2, "_",  part, ".log", sep="")))
  cat(date(), '\n')
  X <- as.numeric(gl(2, n / 2)) - 1
  Z <- cbind(rep(1, n), rnorm(n))
  methods <- c('B2 only', "Ye", "Ye_Boot_res", "Ye_Boot_sim")
  pvs <- array(NA, c(length(a0s), length(a1s), length(r1s), length(r2s), length(n1s), length(methods), niter), 
               dimnames = list(a0 = paste(a0s), a1 = paste(a1s), r1 = paste(r1s), r2 = paste(r2s), n1 = paste(n1s),
                               method = methods, iter = paste(1 : niter)))
  
  for (i in 1 : niter) {
    for (a0 in a0s) {
      cat('.')
      for (a1 in a1s) {
        cat('*')
        for (r1 in r1s) {
          cat('!')
          for (r2 in r2s){
            vb <-  1 / (1 + r2)
            v2 <- r2 * vb
            v1 <- r1 * r2 * vb
            
            Eb <- rnorm(n, sd = sqrt(vb))
            Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
            Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Eb + Et
            Z.r.a <- Z[1 : (n / 2), ]
            Eb.r.a <- Eb[1 : (n / 2)]
            Et.r.a <- rnorm(n / 2, sd =  sqrt(v2))
            Y.r.a <- a1 + Z.r.a %*% b + Eb.r.a + Et.r.a  # 1/(1+r2) correlation 
            for (n1 in n1s) {
              ind.r <- 1 : n1
              Y.r <- Y.r.a[ind.r]
              
              pvs[paste(a0), paste(a1), paste(r1), paste(r2), paste(n1), 'Ye', i] <- 
                batch.ReMeasure.S1(Y, X, Z, ind.r, Y.r)$p.value
            }
          }
        }
      }
    }
  }
  res <- pvs
  cat('\n', date(), '\n')
  sink()
  save(res, file=file.path(resdir, paste(prefix2, "_res",  part, ".Rdata", sep="")))
  return(res)
}


prefix <- 'BC1'
resdir <- paste0("./", prefix)
source(paste0("./", prefix, "/Cluster_mayo_plus.R"))

paras.list <- list()
tempdir.list <- list()
dat.list <- list()
func.list <- list()


paras <- list()
paras$prefix <- prefix
paras$prefix2 <- 'Sim1'
paras$resdir <- resdir
paras$b <- c(0, -0.5)
paras$n <- 100
paras$niter <- 2
paras$a0s <- c(0, 0.25, 0.5)
paras$a1s <- c(0, 0.5, 1)
paras$r1s <- c(0.5, 1, 2)
paras$r2s <- c(0.5, 1, 2)
paras$n1s <- seq(5, 50, len = 10)


paras.list[[paras$prefix2]] <- paras
tempdir.list[[paras$prefix2]] <- file.path(resdir, paras$prefix2)
dat.list[[paras$prefix2]] <- 1:50
func.list[[paras$prefix2]] <- Wrap_func
clsapply.plus(dat.list, func.list, paras.list, tempdir.list, resdir = resdir, queque="1-hour")

