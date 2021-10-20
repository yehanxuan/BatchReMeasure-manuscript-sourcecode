p = 5
beta = c(1, 1, 1, 1, 1)


S1Data = Generate_S1(nc1, nt2, nc2, p, beta, a0, a1, sigma1, sigma2, rho)

Yc1 = S1Data$Yc1
Yt2 = S1Data$Yt2
Yc2 = S1Data$Yc2
Zc1 = S1Data$Zc1
Zt2 = S1Data$Zt2
Zc2 = S1Data$Zc2
Index = S1Data$Index
