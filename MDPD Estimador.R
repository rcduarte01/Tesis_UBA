#Estímador de  mínima divergencia de potencia de densidad
psibeta1_MDPDE <- function(x, muestra, alpha){
  (digamma(x[1] + x[2]) - digamma(x[1]) + log(muestra))*
    dbeta(muestra, shape1 = x[1], shape2 = x[2])^(alpha) - 
    
    (gamma(x[1] + x[2])/(gamma(x[1])*gamma(x[2])))^(alpha + 1)*
    (gamma(alpha*x[1] - alpha + x[1])*gamma(alpha*x[2] - alpha + x[2]))/
    gamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2])*
    (digamma(x[1] + x[2]) - digamma(x[1]) + digamma(alpha*x[1]-alpha + x[1]) -
       digamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2] ) )
}
psibeta2_MDPDE <- function(x, muestra, alpha){
  (digamma(x[1] + x[2]) - digamma(x[2]) + log(1-muestra))*
    dbeta(muestra, shape1 = x[1], shape2 = x[2])^(alpha) - 
    
    (gamma(x[1]+x[2])/(gamma(x[1])*gamma(x[2])))^(alpha + 1)*
    (gamma(alpha*x[1] - alpha + x[1])*gamma(alpha*x[2] - alpha + x[2]))/
    gamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2])*
    (digamma(x[1] + x[2]) - digamma(x[2]) + digamma(alpha*x[2] - alpha + x[2]) -
       digamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2] ) )
}

psii_MDPDE <- function(x,muestra, alpha){
  F1 <- psibeta1_MDPDE(x,muestra, alpha)
  F2 <- psibeta2_MDPDE(x,muestra, alpha)
  return(c(F1,F2))
}

psii_s_MDPDE <- function(x,muestra, alpha){
  F1 <- sum(psibeta1_MDPDE(x,muestra, alpha))
  F2 <- sum(psibeta2_MDPDE(x,muestra, alpha))
  return(F1^2+F2^2)
}

MDPDE_estimador_beta <- function(x, alpha=0){
  #ajuste <- rlm(x~1, method = "MM")
  #xbar <-  ajuste$coefficients
  #s2 <- summary(ajuste)$sigma^2
  
  
  xbar <- mean(x)
  s2 <- var(x)
  
  
  alpha0 <- -xbar*(-xbar+xbar^2+s2)/s2
  beta0 <- (alpha0-alpha0*xbar)/xbar
  #iniciales<- fitdist(x, "beta")$estimate
  
  #alpha0 <- iniciales[1]
  #beta0  <- iniciales[2]
  
  #alpha0 <- max(0.0001, (xbar*(1-xbar)/s2 - 1)*xbar)
  #beta0  <- max(0.0001, (xbar*(1-xbar)/s2 - 1)*(1-xbar))
  
  est <- stats::optim(par = c(alpha0,beta0),
                      fn=psii_s_MDPDE,
                      alpha=alpha,
                      muestra = x)$par
  
  return(est)
}


# varianza del estimador
J_beta <- function(x,alpha=0){
  J <- matrix(0, nrow = 2, ncol = 2)
  #---------------------------------------------------------------------------------------
  J [1,1] <- (gamma(x[1] + x[2])/(gamma(x[1])*gamma(x[2])))^(alpha + 1)*
    (gamma(alpha*x[1] - alpha + x[1])*gamma(alpha*x[2] - alpha + x[2]))/
    gamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2])*
    
    ( (psigamma(x[1] + x[2]) - psigamma(x[1]))^2 + 
        2*(psigamma(x[1]+x[2])-psigamma(x[1]))*
        (psigamma(alpha*x[1] - alpha + x[1]) -psigamma(alpha*x[1]+alpha*x[2]-2*alpha + x[1] + x[2]) ) + 
        (psigamma(alpha*x[1] - alpha + x[1]) - psigamma(alpha*x[1]+alpha*x[2]-2*alpha + x[1] + x[2]))^2+
        psigamma(alpha*x[1] - alpha + x[1], deriv = 1 ) - 
        psigamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2],deriv = 1))
  #---------------------------------------------------------------------------------------  
  J[1,2] <- (gamma(x[1] + x[2])/(gamma(x[1])*gamma(x[2])))^(alpha + 1)*
    (gamma(alpha*x[1] - alpha + x[1])*gamma(alpha*x[2] - alpha + x[2]))/
    gamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2])*
    
    ((psigamma(x[1] + x[2])-psigamma(x[1]))*(psigamma(x[1] + x[2])-psigamma(x[2])) +
       
       (psigamma(x[1] + x[2])-psigamma(x[1]))*(psigamma(alpha*x[2]- alpha + x[2]) -
                                                 psigamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1]+x[2])) +
       
       (psigamma(x[1] + x[2])-psigamma(x[2]))*(psigamma(alpha*x[1]- alpha + x[1]) -
                                                 psigamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1]+x[2])) + 
       
       (psigamma(alpha*x[1]-alpha+x[1]) - psigamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1] + x[2]) )*
       
       (psigamma(alpha*x[2]- alpha + x[2]) - psigamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1] + x[2]) ) -
       
       psigamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1] + x[2], deriv = 1))
  
  #---------------------------------------------------------------------------------------
  J[2,1] <- (gamma(x[1] + x[2])/(gamma(x[1])*gamma(x[2])))^(alpha + 1)*
    (gamma(alpha*x[1] - alpha + x[1])*gamma(alpha*x[2] - alpha + x[2]))/
    gamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2])*
    
    ((psigamma(x[1] + x[2])-psigamma(x[1]))*(psigamma(x[1] + x[2])-psigamma(x[2])) +
       
       (psigamma(x[1] + x[2])-psigamma(x[1]))*(psigamma(alpha*x[2]- alpha + x[2]) -
                                                 psigamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1]+x[2])) +
       
       (psigamma(x[1] + x[2])-psigamma(x[2]))*(psigamma(alpha*x[1]- alpha + x[1]) -
                                                 psigamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1]+x[2])) + 
       
       (psigamma(alpha*x[1]-alpha+x[1]) - psigamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1] + x[2]) )*
       
       (psigamma(alpha*x[2]- alpha + x[2]) - psigamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1] + x[2]) ) -
       
       psigamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1] + x[2], deriv = 1))
  
  
  J[2,2] <- (gamma(x[1] + x[2])/(gamma(x[1])*gamma(x[2])))^(alpha + 1)*
    (gamma(alpha*x[1] - alpha + x[1])*gamma(alpha*x[2] - alpha + x[2]))/
    gamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2])*
    
    ( (psigamma(x[1] + x[2]) - psigamma(x[2]))^2 + 
        2*(psigamma(x[1]+x[2])-psigamma(x[2]))*
        (psigamma(alpha*x[2] - alpha + x[2]) -psigamma(alpha*x[1]+alpha*x[2]-2*alpha + x[1] + x[2]) ) + 
        (psigamma(alpha*x[2] - alpha + x[2]) - psigamma(alpha*x[1]+alpha*x[2]-2*alpha + x[1] + x[2]))^2+
        psigamma(alpha*x[2] - alpha + x[2], deriv = 1 ) - 
        psigamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2],deriv = 1))
  
  return(J)
  
}

xi_beta <- function(x, alpha=0){
  xi1 <- (gamma(x[1] + x[2])/(gamma(x[1])*gamma(x[2])))^(alpha + 1)*
    (gamma(alpha*x[1] - alpha + x[1])*gamma(alpha*x[2] - alpha + x[2]))/
    gamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2])*
    (digamma(x[1] + x[2]) - digamma(x[1]) + digamma(alpha*x[1]-alpha + x[1]) -
       digamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2] ) )
  
  xi2 <- (gamma(x[1]+x[2])/(gamma(x[1])*gamma(x[2])))^(alpha + 1)*
    (gamma(alpha*x[1] - alpha + x[1])*gamma(alpha*x[2] - alpha + x[2]))/
    gamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2])*
    (digamma(x[1] + x[2]) - digamma(x[2]) + digamma(alpha*x[2] - alpha + x[2]) -
       digamma(alpha*x[1] + alpha*x[2] - 2*alpha + x[1] + x[2] ) )
  
  return(matrix(c(xi1,xi2), nrow = 2, ncol = 1))
  
}

K_beta <- function(x, alpha=0){
  xi <- xi_beta(x,alpha)
  K  <- J_beta(x,2*alpha) - xi%*%t(xi)
  return(K)
}

Ajuste_MDPE_estimador_beta<- function(muestra, alpha=0, conf=0.95){
  aa <- (1-conf)
  zz <- qnorm(aa/2, lower.tail = F)
  n <- length(muestra)
  
  # estimación puntual
  estimador <- MDPDE_estimador_beta(muestra, alpha)
  
  #matriz de varianza asintótica
  J <- J_beta(estimador, alpha)
  K <- K_beta(estimador, alpha)

  C <- solve(J)%*%K%*%solve(J)/n
  err_est <- sqrt(diag(C))
  
  LI_alpha <- estimador[1] - zz*err_est[1]
  LS_alpha <- estimador[1] + zz*err_est[1]
  
  LI_beta <- estimador[2] - zz*err_est[2]
  LS_beta <- estimador[2] + zz*err_est[2]
  
  CI <- rbind(c(LI_alpha,LS_alpha),
              c(LI_beta,LS_beta))
  
  ajuste<-cbind(estimador,err_est,CI)
  rownames(ajuste)<- c("a","b")
  colnames(ajuste)<-c("Estimaciones","Error Est Asintotico", "LI","LS")
  
  return(ajuste)
}

Finfluencia_MDPDE<- function(a,b, alpha){
  
  JJ <- J_beta(x=c(a,b), alpha)
  XI <- xi_beta(x=c(a,b), alpha)
  x0 <- seq(from=0,to=1, length=100)
  
  FF1 <- c()
  FF2 <- c()
  for(i in 1:100){
    UU <- matrix(c(psigamma(a+b) - psigamma(a) + log(x0[i]),
                   psigamma(a+b) - psigamma(b) + log(1-x0[i])), nrow = 2)
    
    FF <- solve(JJ)%*%(UU*dbeta(x0[i],shape1 = a, shape2 = b)^alpha-XI)
    FF1[i] <- FF[1,]
    FF2[i] <- FF[2,]
    
  }
  
  par(mfrow=c(1,2))
  plot(x0,FF1, type = "l", lwd=2,
       ylab=TeX("$IF_1(x_0)$"), xlab=TeX("$x_0$"),
       bty="n", axes = F,cex.lab=1)
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  plot(x0,FF2, type = "l", lwd=2,
       ylab=TeX("$IF_2(x_0)$"), xlab=TeX("$x_0$"),
       bty="n", axes = F,cex.lab=1)
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  par(mfrow=c(1,1))
  
  
}














