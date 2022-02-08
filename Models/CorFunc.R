# Need to consider which the level of X corresponds control/case
require(lmvar)
require(MASS)

loglik.r1 <- function (paras, Y, X, Z, ind.r, Y.r) {
	n <- length(Y)
	n1 <- length(ind.r)
	p <- ncol(Z)
	X <- as.numeric(factor(X)) - 1
	XX <- cbind(X, X)
	a0 <- paras[1]
	a1 <- paras[2]
	b <- paras[3 : (p + 2)]
	v1 <- exp(paras[p + 3])
	v2 <- exp(paras[p + 4])	
	rho <- (1 - exp(paras[p + 5])) / (1 + exp(paras[p + 5]))
	
	mu <- Z %*% b + XX %*% c(a0, a1)
	s2 <- c(v1, v2)[X + 1]
	
	Z.r <- Z[ind.r, ]
	Y.o <- Y[ind.r]
	mu.r <- a1 + Z.r %*% b + rho * sqrt(v2 / v1) * (Y.o - Z.r %*% b)
	s2.r <- rep((1 - rho^2) * v2, n1)
	
	-sum(stats::dnorm(c(Y, Y.r), mean = c(mu, mu.r), sd = sqrt(c(s2, s2.r)), log = TRUE))
	
}


### Generic MLE
batch.correct.r1 <- function (Y, X, Z, ind.r, Y.r, start = NULL, method = 'BFGS') {
	# Z include the intercept
  # X treatment and case 
  start_S1 = proc.time()[1]
	call <- match.call()
	if (is.null(start)) {
		
		Y.o <- Y[ind.r]
		c2 <-  mean(Y.r - Y.o) ## a1
		lm.obj <- lm(Y ~ X + Z - 1)
		lm.coef <- coef(lm.obj)
		c1 <- lm.coef[1] - c2
		c3 <- lm.coef[2 : length(lm.coef)]
		c4 <- c5 <-  log(sigma(lm.obj)^2)   # sigma residual standard error 
		c6 <- cor(Y.o - resid(lm.obj)[ind.r], Y.r - resid(lm.obj)[ind.r])
		c6 <- log((1 - c6) / (1 + c6))
		
		start <- c(c1, c2, c3, c4, c5, c6)  # a0, a1, beta, sigma, rho, here is slightly
		#c5 not c6
	}
	
	oout <- stats::optim(start, loglik.r1, method = method, hessian = TRUE, 
			Y = Y, X = X, Z = Z, ind.r = ind.r, Y.r = Y.r)
	coef <- oout$par
	try.obj <- try({	vcov <- solve(oout$hessian)
				stat <- coef[1] / sqrt(vcov[1, 1])
				pv <- 2 * stats::pnorm(-abs(stat))
	#			pv <- 2 * stats::pt(-abs(coef[1] / sqrt(vcov[1, 1])), df = length(ind.r) + sum(X == 1) - ncol(Z))
			})
    if (class(try.obj) == 'try-error') {
		stat <- NA
		pv <- NA
    }
	
	a0H = as.numeric(coef[1])
	a1H = as.numeric(coef[2])
	betaH = as.vector( coef[3:4])
	sigma1H = exp( coef[5]/2 )
	sigma2H = exp( coef[6]/2 )
	
	rhoH = (1-exp(coef[7]))/(1 + exp(coef[7]))
	a0Var = as.numeric(vcov[1, 1])
	Time_S1 = proc.time()[1] - start_S1
	#return(list(call = call, stat = stat, p.value = pv, optim.out = oout, optim.st = start ))
	return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "sigma1" = sigma1H, 
	            "sigma2" = sigma2H, "rho" = rhoH, "beta" = NULL, "objVec" = NULL, "Time" = Time_S1))
}

batch.correct.r1.bt <- function (Y, X, Z, ind.r, Y.r, start = NULL, method = 'BFGS', B = 50) {
	# Z include the intercept
	call <- match.call()
	
	obj <- batch.correct.r1(Y, X, Z, ind.r, Y.r, start, method)
	stat <- obj$stat
	stat.b <- NULL
	for (b in 1 : B) {
		ind <- sample(1 : length(Y.r), repl = TRUE)
		ind.r.b <- ind.r[ind]
		Y.r.b <- Y.r[ind]
		ind.nr.b <- sample(setdiff(1 : length(Y), ind.r), repl = TRUE)
		ind.b <- c(ind.r.b, ind.nr.b)
		Y.b <- Y[ind.b]
		X.b <- X[ind.b]
		Z.b <- Z[ind.b, , drop = FALSE]
		
		try.obj <- try(obj.b <- batch.correct.r1(Y.b, X.b, Z.b, 1 : length(Y.r), Y.r.b, start = obj$optim.out$par, method))
		if (class(try.obj) != 'try-error') {
			if (!is.na(obj.b$stat) ) {
				stat.b <- c(stat.b, obj.b$stat)
			}
		}
	}
	pv <- 2 * stats::pnorm(-abs(stat), sd = sd(stat.b))
	

	return(list(call = call, stat = stat, p.value = pv))
	
}

# batch.correct.r1.naive1 <- function (Y, X, Z, ind.r, Y.r) {
# 	# location/scale model
# 	# X <- as.numeric(factor(X)) - 1
# 	Y.o <- Y[ind.r]
# 	
# 	l <- mean(Y.r - Y.o)
# 	s <- sd(Y.r) / sd(Y.o)
# 	
# 	Y[X == 0] <- (Y[X == 0] + l) 
# 	m <- mean(Y[X == 0])
# 	Y[X == 0] <- (Y[X == 0]  - m) * s + m 
# 	
# 	lm.obj <- lm(Y ~ X + Z - 1)
# 	pv <- coef(summary(lm.obj))[1, 4]
# 
# 	return(list(p.value = pv))
# 	
# }

# batch.correct.r1.naive2 <- function (Y, X, Z, ind.r, Y.r) {
# 	#X <- as.numeric(factor(X)) - 1
# 	Y <- c(Y.r, Y[X == 1])
# 	Z <- rbind(Z[ind.r, ], Z[X == 1, ])
# 	X <- c(rep(0, length(ind.r)), X[X == 1])
# 	
# 	lm.obj <- lm(Y ~ X + Z - 1)
# 	pv <- coef(summary(lm.obj))[1, 4]
# 	
# 	return(list(p.value = pv))
# }

#batch.correct.r1.naive3 <- function (Y, X, Z, ind.r, Y.r, n = 10000) {
#	X <- as.numeric(factor(X)) - 1
#	
#	diff.obs <- Y.r -  Y[ind.r]
#
#	null.mean <- sapply(1 : n, function (i) mean(sample(diff.obs, repl  = TRUE)))
#	
#	X.mu <- model.matrix( ~ X + Z - 1)
#	X.sigma <- model.matrix(~ X)
#	lm.obj <- lmvar(y = as.vector(Y), X_mu = X.mu, X_sigma = X.sigma, intercept_mu = FALSE, intercept_sigma = FALSE)
#	
#	coeffs <- coef(summary(lm.obj))
#	estimate <- coeffs[1, 1]
#	stderr <- coeffs[1, 2]
#	
#	pv <- mean(2 * pnorm(-abs((estimate - null.mean) / stderr)))
#	
#	return(list(p.value = pv))
#
#
#}

find.ci.overlap.max <- function (est1, est1.se, est2, est2.se, est2.df) {
	
	fun.tmp <- function (p, est1, est1.se, est2, est2.se, est2.df) {
		if (est1 <= est2) {
			ci1.right <- est1 +  abs(qnorm((1 - p) / 2)) * est1.se
			ci2.left <- est2 - abs(qt((1 - p) / 2, df = est2.df)) * est2.se
			return(ci1.right - ci2.left)
		} else {
			ci1.left  <- est1 - abs(qnorm((1 - p) / 2)) * est1.se
			ci2.right <- est2  + abs(qt((1 - p) / 2, df = est2.df)) * est2.se
			return(ci1.left - ci2.right)
		}
	}
	
	return(uniroot(fun.tmp, c(0, 1), est1 = est1, est1.se = est1.se, 
					est2 = est2, est2.se = est2.se, est2.df = est2.df)$root)

	
}


is.ci.overlap <- function (est1, est1.se, est2, est2.se, est2.df, p = 0.95) {
	
	
	if (est1 <= est2) {
		ci1.right <- est1 +  abs(qnorm((1 - p) / 2)) * est1.se
		ci2.left <- est2 - abs(qt((1 - p) / 2, df = est2.df)) * est2.se
		return((ci1.right >= ci2.left))
	} else {
		ci1.left  <- est1 - abs(qnorm((1 - p) / 2)) * est1.se
		ci2.right <- est2  + abs(qt((1 - p) / 2, df = est2.df)) * est2.se
		return((ci2.right >=  ci1.left))
	}
	
	
}


batch.correct.r1.naive3 <- function (Y, X, Z, ind.r, Y.r) {
	
	X <- as.numeric(factor(X)) - 1

	X.mu <- model.matrix( ~ X + Z - 1)
	X.sigma <- model.matrix(~ X)
	lm.obj <- lmvar(y = as.vector(Y), X_mu = X.mu, X_sigma = X.sigma, intercept_mu = FALSE, intercept_sigma = FALSE)
	
	coeffs <- coef(summary(lm.obj))
	est1 <- coeffs[1, 1]
	est1.se <- coeffs[1, 2]
	
	t.obj  <- t.test(Y.r, Y[ind.r], pair = TRUE)
	
	est2 <- t.obj$estimate
	est2.se <- t.obj$stderr
	est2.df <- t.obj$parameter
	
#	p <- find.ci.overlap.max(est1, est1.se, est2, est2.se, est2.df) 
#	
#	return(list(p.value = 1 - p^2))
	
	if (is.ci.overlap(est1, est1.se, est2, est2.se, est2.df)) {
		return(list(p.value = 0.01))
	} else {
		return(list(p.value = 0.1))
	}
	
}


loglik.r2 <- function (paras, Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2) {
	n <- length(Y)
	n1 <- length(ind.r1)
	n2 <- length(ind.r2)
	
	p <- ncol(Z)
	X <- as.numeric(factor(X)) - 1
	XX <- cbind(X, X)
	a0 <- paras[1]
	a1 <- paras[2]
	a2 <- paras[3]
	b <- paras[4 : (p + 3)]
	v1 <- exp(paras[p + 4])
	v2 <- exp(paras[p + 5])	
	v3 <- exp(paras[p + 6])
	

	rho1 <- (1 - exp(paras[p + 7])) / (1 + exp(paras[p + 7]))
	rho2 <- sqrt(v1 / v2)
	
	mu <- Z %*% b + XX %*% c(a0, a1)
	s2 <- c(v1, v2)[X + 1]
	
	Z.r1 <- Z[ind.r1, ]
	Y.o1 <- Y[ind.r1]
	mu.r1 <- a2 + Z.r1 %*% b + rho1 * sqrt(v3 / v1) * (Y.o1 - Z.r1 %*% b)
	s2.r1 <- rep((1 - rho1^2) * v3, n1)
	
	Z.r2 <- Z[ind.r2, ]
	Y.o2 <- Y[ind.r12]
	mu.r2 <- a2 + Z.r2 %*% b + rho2 * sqrt(v3 / v2) * (Y.o2 - Z.r2 %*% b)
	s2.r2 <- rep((1 - rho2^2) * v3, n2)
	
	-sum(stats::dnorm(c(Y, Y.r), mean = c(mu, mu.r1, mu.r2), sd = sqrt(c(s2, s2.r1, s2.r2)), log = TRUE))
	
}

batch.correct.r2 <- function (Y, X, Z, ind.r1, ind.r2, Y.r1, Y.r2, start = NULL, method = 'BFGS') {
	# Z include the intercept
	call <- match.call()
	if (is.null(start)) {
		
		Y.o1 <- Y[ind.r1]
		c3 <-  mean(Y.r1 - Y.o1)
		
		Y.o2 <- Y[ind.r2]
		c23 <-  mean(Y.r2 - Y.o2)
		c2 <- c3 - c23
		
		lm.obj <- lm(Y ~ X + Z - 1)
		
		lm.coef <- coef(lm.obj)
		c1 <- lm.coef[1] - c2
		c4 <- lm.coef[2 : length(lm.coef)]
		c5 <- c6 <-  log(sigma(lm.obj)^2)
		c7 <- cor(Y.o1 - resid(lm.obj)[ind.r1], Y.r1 - resid(lm.obj)[ind.r1])
		c7 <- log((1 - c7) / (1 + c7))
		
		start <- c(c1, c2, c3, c4, c6, c6, c7)
	}
	
	oout <- stats::optim(start, loglik.r1, method = method, hessian = TRUE, 
			Y = Y, X = X, Z = Z, ind.r1 = ind.r1, ind.r2 = ind.r2,  Y.r1 = Y.r1, Y.r2 = Y.r2)
	coef <- oout$par
	try.obj <- try({	vcov <- solve(oout$hessian)
				pv <- 2 * stats::pnorm(-abs(coef[1] / sqrt(vcov[1, 1])))	
			})
	if (class(try.obj) == 'try-error') {
		pv <- NA
	}
	return(list(call = call, p.value = pv, optim.out = oout, optim.st = start ))
	
}


oneReplicate_Gen = function(seedJ){
  set.seed(seedJ + repID * 300)
  source("./oneReplicate/oneReplicate-New-S1.R")
  Estimate = batch.correct.r1(Y, X, Z, ind.r, Y.r)
  a0H = Estimate$a0 
  a0Var = Estimate$a0Var 
  a1H = Estimate$a1 
  betaH = Estimate$beta
  rhoH = Estimate$rho
  sigma1H = Estimate$sigma1
  sigma2H = Estimate$sigma2
  objVec = Estimate$objVec
  Time = Estimate$Time 
  return(list("a0" = a0H, "a0Var" = a0Var, "a1" = a1H, "sigma1" = sigma1H,
              "sigma2" = sigma2H, "rho" = rhoH, 
              "beta" = betaH,"objVec" = objVec, "Time" = Time))
}

oneReplicateWrap_Gen = function(seedJ) {
  eval =  try( oneReplicate_Gen(seedJ) )
  return(eval)
}

#require(lmvar); require(MASS)
#X = model.matrix(~ Bwt - 1, cats)
#
## Carry out the fit with the same model matrix for mu (the expected heart weight) and for log sigma (the standard deviation)
#fit = lmvar(cats$Hwt, X_mu = X, X_sigma = X)


require(utils) # for str

## some platforms hit zero exactly on the first step:
## if so the estimated precision is 2/3.
# f <- function (x, a) x - a
# str(xmin <- uniroot(f, c(0, 1), tol = 0.0001, a = 1/3))
