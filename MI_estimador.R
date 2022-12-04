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