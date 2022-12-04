# este script contiene el código de las funciones
# relacionadas con el MI Estimador

library(MASS)
library(fitdistrplus)
library(latex2exp)
# funciones psi
psibeta1 <- function(x, muestra){
  pbeta(muestra, shape1=x[1], shape2=x[2]) - 1/2
}
psibeta2 <- function(x, muestra){
  pbeta(muestra, shape1=x[1], shape2=x[2])^2 - 1/3
}
psii <- function(x,muestra){
  F1 <- psibeta1(x,muestra)
  F2 <- psibeta2(x,muestra)
  return(c(F1,F2))
}
psiI <- function(x,muestra){
  F1 <- sum(pbeta(muestra, shape1=x[1], shape2=x[2]) -1/2)
  F2 <- sum(pbeta(muestra, shape1=x[1], shape2=x[2])^2 - 1/3)
  return(F1^2 + F2^2)
}

# estimador puntual
MI_estimador_beta <- function(x){
  ajuste <- rlm(x~1, method = "MM")
  xbar <-  ajuste$coefficients
  s2 <- summary(ajuste)$sigma^2
  
  alpha0 <- max(0.0001, (xbar*(1-xbar)/s2 - 1)*xbar)
  beta0  <- max(0.0001, (xbar*(1-xbar)/s2 - 1)*(1-xbar))
  
  est <- stats::optim(par = c(alpha0,beta0),
                      psiI,
                      muestra = x)$par
  return(est)
}

# varianza del estimador
Bfunc_MI <- function(muestra){
  n <- length(muestra)
  estimadores <- MI_estimador_beta(muestra)
  B <- matrix(0,nrow = 2, ncol = 2)
  for(i in 1:n){
    B <- B + rootSolve::gradient(psii,
                                 x=estimadores,
                                 muestra=muestra[i])
  }
  -B/n
}
Afunc_MI <- function(muestra){
  n <- length(muestra)
  estimadores <- MI_estimador_beta(muestra)
  A <- matrix(0,nrow = 2, ncol = 2)
  for(i in 1:n){
    A <- A + psii(x=estimadores, muestra=muestra[i]) %*% 
      t(psii(x=estimadores, muestra=muestra[i]))
  }
  A/n
}
Ajuste_MI_estimador_beta<- function(muestra, conf=0.95){
  alpha <- (1-conf)
  zz <- qnorm(alpha/2, lower.tail = F)
  n <- length(muestra)
  
  # estimación puntual
  estimador <- MI_estimador_beta(muestra)
  
  #matriz de varianza asintótica
  A <- Afunc_MI(muestra)
  B <- Bfunc_MI(muestra)
  
  C <- solve(B)%*%A%*%solve(t(B))/n
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

# Función de Influencia
Finfluencia_MI <- function(a,b){
  x0<- seq(from=0, to=1, length=100)
  
  muestra <- rbeta(1000, shape1 = a, shape2 = b)
  B  <- Bfunc_MI(muestra)
  FF <- solve(B)%*%sapply(x0, psii, x = c(a,b))
  
  par(mfrow=c(1,2))
  
  plot(x0,FF[1,], type = "l", lwd=2,
       ylab=TeX("$IF_1(x_0)$"), xlab=TeX("$x_0$"),
       bty="n", axes = F,cex.lab=1)
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  plot(x0,FF[2,], type = "l", lwd=2,
       ylab=TeX("$IF_2(x_0)$"), xlab=TeX("$x_0$"),
       bty="n", axes = F,cex.lab=1)
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  par(mfrow=c(1,1))
  
}

mm <- rbeta(100, 2,8)
ajuste1 <- Ajuste_MI_estimador_beta(mm)



