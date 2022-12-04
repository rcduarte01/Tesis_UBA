psibeta1_MLE <- function(x, muestra){
digamma(x[1] + x[2]) - digamma(x[1]) + log(muestra)
}

psibeta2_MLE <- function(x, muestra){
  digamma(x[1] + x[2]) - digamma(x[2]) + log(1-muestra)
}

psii_MLE <- function(x,muestra){
  F1 <- psibeta1_MLE(x,muestra)
  F2 <- psibeta2_MLE(x,muestra)
  return(c(F1,F2))
}

Jfun_MLE <- function(muestra, alfa, beta){
  n <- length(muestra)
  J <- matrix(0,nrow = 2, ncol = 2)
  
  for(i in 1:n){
    J <- J + rootSolve::gradient(psii_MLE,
                                 x=c(alfa,beta),
                                 muestra=muestra[i])
  }
  J/n
  
}

Finfluencia_MLE <- function(alfa,beta){
  x0 <- seq(from=0.01, to=0.99, length=100)
  muestra <- rbeta(100, shape1 = alfa, shape2 = beta)
  B  <- Jfun_MLE(muestra, alfa, beta)
  FF <- -solve(B)%*%sapply(x0, psii_MLE, x = c(alfa,beta))
  
  par(mfrow=c(1,2))
  plot(x0,FF[1,], type = "l", lwd=2, 
       ylab=TeX("$IF_1(x_0)$"), xlab=TeX("$x_0$"),
       axes=F)
  
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  
  plot(x0,FF[2,], type = "l", lwd=2, ylab=TeX("$IF_2(x_0)$"),
       xlab=TeX("$x_0$"),
       axes=F)
  
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  
  par(mfrow=c(1,1))

  
}
