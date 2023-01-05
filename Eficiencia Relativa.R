# eficiencia relativa
MSE <- function(theta, theta_estimado){
  mean((theta-theta_estimado)^2, na.rm=TRUE)
}

eficienciaRelativa<- function(N, a, b, rob=TRUE, alpha=0){
  
  MIEst  <- matrix(nrow = N , ncol = 2)
  MVEst  <- matrix(nrow = N , ncol = 2)
  MDPEst <- matrix(nrow = N,  ncol = 2)
  
  for(i in 1:N){
    skip_to_next <- FALSE
    tryCatch({
      mm <- rbeta(100, shape1 = a, shape2 = b)
      if(rob){
        xbar <- median(mm)
        s2   <- mad(mm)^2
        alpha0 <- -xbar*(-xbar+xbar^2+s2)/s2
        beta0 <- (alpha0-alpha0*xbar)/xbar
      }else{
        xbar <- mean(mm)
        s2   <- var(mm)
        alpha0 <- -xbar*(-xbar+xbar^2+s2)/s2
        beta0 <- (alpha0-alpha0*xbar)/xbar
      }
      
      MIEst[i,]  <- MI_estimador_beta(mm, rob)
      MDPEst[i,] <- MDPDE_estimador_beta(mm, alpha, rob)
      MVEst[i,] <- fitdistr(mm, densfun = "beta",
                            list(shape1=alpha0, shape2=beta0))$estimate
      
    }, error = function(e) {skip_to_next <<- TRUE },
    message = function(m) {skip_to_next <<- TRUE })
    
    if(skip_to_next) { next }
    
  }
  
  alphaMSE <- MSE(a,MVEst[,1])/MSE(a,MIEst[,1])
  betaMSE  <- MSE(b,MVEst[,2])/MSE(b,MIEst[,2])
  
  alphaMSE2 <- MSE(a,MVEst[,1])/MSE(a,MDPEst[,1])
  betaMSE2  <- MSE(b,MVEst[,2])/MSE(b,MDPEst[,2])

  retorno <- c(alphaMSE,betaMSE,alphaMSE2,betaMSE2)
}

set.seed(2022)
N <- 1000
ss <- eficienciaRelativa(N,
                         rob=FALSE,
                         a = 0.8,
                         b = 0.8,
                         alpha=0.4)
round(ss,4)












