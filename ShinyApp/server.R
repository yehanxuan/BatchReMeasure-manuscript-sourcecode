## server.R 
Calculate_Power_S1 = function(nc1, nt2, rho, a0, a1, nc2, alpha, seedJ = 2){
  set.seed(seedJ)
  n = nc1 + nt2 
  X =  as.numeric(gl(2, n / 2)) - 1
  Z <- cbind(rep(1, n), rnorm(n))
  b <- c(0, -0.5)
  v1 = 1
  v2 = 1
  Et <- rnorm(n, sd = ifelse (X == 0, sqrt(v1), sqrt(v2)))
  Y <- Z %*% b + cbind(X, X) %*% c(a0, a1) + Et
  Z.r.a <- Z[1 : (n / 2), ]
  Et.r.a <- Et[1 : (n / 2)]
  Y.r.a <- a1 + Z.r.a %*% b + rho * sqrt(v2) * Et.r.a/ sqrt(v1) +
    rnorm(n/2, sd = sqrt( (1 - rho^2) * v2 ) )
  
  ind.r <- 1:nc2 
  Y.r = Y.r.a[ind.r]
  ind0 <- X == 0
  ind1 <- X == 1
  Yc1 <- Y[ind0]
  Yt2 <- Y[ind1]
  Zc1 <- Z[ind0, , drop = F]
  Zt2 <- Z[ind1, , drop = F]
  Zc2 = Zc1[ind.r, , drop = F]
  Yc2 = Y.r
  
  a0Var = Oracle_Variance_a0(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, sqrt(v1), sqrt(v2), rho, ind.r)
  C = a0/sqrt(a0Var)
  C1 = qnorm(alpha/2, lower.tail = TRUE) - C
  C2 = qnorm(alpha/2, lower.tail = TRUE) + C
  Power = pnorm(C1) + pnorm(C2)
  return(Power)
}

PowerRatio = function(nc1, nt2, n2s, rho, a0, alpha, seedJ = 2) {
  a1 = 0.5
  r1 = 1
  Power.a = array(NA, c(length(a0), length(rho), length(n2s)), 
                  dimnames = list(TrueEffect = paste(a0), Cor = paste(rho), 
                                  RemeasureNo = paste(n2s)) ) 
  for (nc2 in n2s) {
    pw = Calculate_Power_S1(nc1, nt2, rho, a0, a1, nc2, alpha, seedJ)
    Power.a[as.character(a0), as.character(rho), as.character(nc2)] = pw
  }
  Power.df <- reshape::melt(Power.a)
  colnames(Power.df)[ncol(Power.df)] <- 'Power'
  pv.max = Power.df$Power[Power.df$RemeasureNo == n2s[length(n2s)]]
  ratio.df = Power.df 
  ratio.df$Power = Power.df$Power/pv.max
  colnames(ratio.df)[ncol(ratio.df)] <- 'Ratio'
  return(list("Power.df" = Power.df,
              "ratio.df" = ratio.df))
}

Oracle_Variance_a0 = function(Zc1, Zt2, Zc2, Yc1, Yt2, Yc2, sigma1, sigma2, rho, Index) {
  nc1 = nrow(Zc1)
  nt2 = nrow(Zt2)
  nc2 = nrow(Zc2)
  R = as.numeric(rho*sigma2/sigma1)
  
  Cov1 = t(Zc1)%*%Zc1/(sigma1^2) + t(Zc2)%*%Zc2*rho^2/(sigma1^2*(1-rho^2))
  Ztilde = t(t(Zc2) - (1-R)*colMeans(Zc2))
  Cov2 = (-2*rho*t(Zc2)%*%Ztilde/(sigma1*sigma2) + t(Ztilde)%*%Ztilde/(sigma2^2))/(1-rho^2)
  Zt2_ct = t(Zt2) - colMeans(Zt2)
  Cov3 = Zt2_ct%*%t(Zt2_ct)/(sigma2^2)
  S = Cov1 + Cov2 + Cov3
  
  A0 = (t(Zc2)/sigma1^2 - t(Ztilde)*rho/(sigma1*sigma2))/(1 - rho^2)
  B0 = (t(Zt2) - colMeans(Zt2)%*%t(rep(1, nt2)))/(sigma2^2)
  
  C0 = - (t(Zc2) - colMeans(Zc2)%*%t(rep(1, nc2)) )*rho/(sigma1*sigma2) +
    (t(Ztilde) - colMeans(Ztilde)%*%t(rep(1, nc2)))/(sigma2^2)
  C0 = C0/(1-rho^2)
  
  k = colMeans(Zt2) - (1 - R)*colMeans(Zc2)
  
  if (nc2 == nc1){
    D0 = 0
    d = 0
  } else {
    D0 = t(Zc1[-Index, , drop = F])/(sigma1^2)
    D = solve(S, D0)
    d = t(D)%*%k
  }
  
  A = solve(S, A0)
  B = solve(S, B0)
  C = solve(S, C0)
  a = t(A)%*%k + R*rep(1, nc2)/nc2
  b = t(B)%*%k + rep(1, nt2)/nt2
  c = t(C)%*%k - rep(1, nc2)/nc2
  Var = (sigma1^2)*( t(a)%*%a + t(d)%*%d ) +
    (sigma2^2)*(t(b)%*%b + t(c)%*%c) + 2*rho*sigma1*sigma2*t(a)%*%c
  Var = as.numeric(Var)
  return(Var)
}

shinyServer( function(input, output){
  
  data <- reactive({
    nc1 <- input$nc1
    nt2 <- input$nt2
    r1 <- 1
    r2 <- input$r2
    a0 <- input$a0 
    a1 <- 0.5
    alpha <- input$alpha
    n2s =  floor( seq(5, min(nc1, nt2), length = 20) )
    out <- PowerRatio(nc1, nt2, n2s, r2, a0, alpha)
    out
  })
  
  observeEvent(input$goButton, {
    # nc1 <- input$nc1
    # nt2 <- input$nt2
    # r1 <- 1
    # r2 <- input$r2
    a0 <- input$a0 
    # a1 <- 0.5
    # alpha <- input$alpha
    # n2s =  floor( seq(5, min(nc1, nt2), length = 20) )
    #Power.df = PowerRatio(nc1, nt2, n2s, r2, a0, alpha)$Power.df
    Power.df = data()$Power.df
    output$plot1 <- renderPlot({
      p1 <- ggplot(Power.df, aes(x = RemeasureNo, y=Power )) + 
        geom_line() + geom_point(size = 2) + xlab("No. of remeasured samples") +
        theme_bw(base_size = 22) + theme(legend.position="bottom") 
      if (a0 == '0') {
        p1 <- p1 + geom_hline(aes(yintercept=0.05), linetype=2) + 
          ylab('Type I error') 
      } else {
        p1 <- p1 + 
          ylab('Absolute Power') }
      p1
    })
  })
  
  observeEvent(input$goButton, {
    ratio.df = data()$ratio.df
    output$plot2 <- renderPlot({
      p2 <- ggplot(ratio.df,aes(x = RemeasureNo, y=Ratio) ) +
        geom_line() + geom_point(size = 2) + xlab("No. of remeasured samples") +
        theme_bw(base_size = 22) + theme(legend.position="bottom") +
        ylab("Relative Power")
      p2
    }) 
  }) 
  
  output$download <- downloadHandler(
    filename = function() {
      #paste0("Powerplot", ".pdf")
      "Powerplot.pdf"
    },
    content = function(file) {
      a0 <- input$a0 
      Power.df = data()$Power.df
      ratio.df = data()$ratio.df
      p1 <- ggplot(Power.df, aes(x = RemeasureNo, y=Power )) + 
        geom_line() + geom_point(size = 2) + xlab("No. of remeasured samples") +
        theme_bw(base_size = 22) + theme(legend.position="bottom") 
      if (a0 == '0') {
        p1 <- p1 + geom_hline(aes(yintercept=0.05), linetype=2) + 
          ylab('Type I error') 
      } else {
        p1 <- p1 + 
          ylab('Absolute Power') }
      
      p2 <- ggplot(ratio.df,aes(x = RemeasureNo, y=Ratio) ) +
        geom_line() + geom_point(size = 2) + xlab("No. of remeasured samples") +
        theme_bw(base_size = 22) + theme(legend.position="bottom") +
        ylab("Relative Power")
      pdf(file, width = 5, height = 4)
      print(p1)
      print(p2)
      dev.off()
    #  for (i in 1:2) {
    #    pdf(paste0(file, i, ".pdf"), width = 5, height = 4)
     #   print(get(paste0("p", i)) )
    #  }
      #ggsave(file, p1, width = 10, height = 8)
      #ggsave(file, p2, width = 10, height = 8)
    }
  )
})