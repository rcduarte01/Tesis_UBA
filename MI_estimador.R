# este script contiene el código de las funciones
# relacionadas con el MI estimador

library(MASS)
library(fitdistrplus)
library(latex2exp)
library(patchwork)
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
  F1 <- sum(pbeta(muestra, shape1=x[1], shape2=x[2]) - 1/2)
  F2 <- sum(pbeta(muestra, shape1=x[1], shape2=x[2])^2 - 1/3)
  
  return(F1^2 + F2^2)
}

# estimador puntual
MI_estimador_beta <- function(x, rob=TRUE){

  if(rob){
    xbar <- median(x)
    s2   <- mad(x)^2
    alpha0 <- -xbar*(-xbar+xbar^2+s2)/s2
    beta0 <- (alpha0-alpha0*xbar)/xbar
  }else{
    xbar <- mean(x)
    s2   <- var(x)
    alpha0 <- -xbar*(-xbar+xbar^2+s2)/s2
    beta0 <- (alpha0-alpha0*xbar)/xbar
  }
  
  est <- stats::optim(par = c(alpha0,beta0),
                      fn=psiI,
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
  matrix(c(1/12,1/12,1/12,4/45), byrow = T, nrow = 2)
  #n <- length(muestra)
  #estimadores <- MI_estimador_beta(muestra)
  #A <- matrix(0,nrow = 2, ncol = 2)
  #for(i in 1:n){
   # A <- A + psii(x=estimadores, muestra=muestra[i]) %*% 
    #  t(psii(x=estimadores, muestra=muestra[i]))
  #}
  #A/n
}
Ajuste_MI_estimador_beta<- function(muestra, conf=0.95, rob=TRUE){
  alpha <- (1-conf)
  n <- length(muestra)
  zz <- qt(alpha/2, lower.tail = F, df=n-1)
  
  # estimación puntual
  estimador <- MI_estimador_beta(muestra, rob)
  
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
  x0 <- seq(from=0, to=1, length=100)
  
  muestra <- rbeta(100, shape1 = a, shape2 = b)
  #B <- BB2(muestra)
  B<- Bfunc_MI(muestra)
  FF1 <- c()
  FF2 <- c()
  for(i in 1:100){
    FF <- solve(B)%*%matrix(psii(x=c(a,b), muestra = x0[i] ),nrow = 2)
    FF1[i] <- FF[1]
    FF2[i] <- FF[2]
  }
  
  df <- data.frame(x0,FF1,FF2)
  grafico1 <- ggplot(df, aes(x=x0,y=FF1))+
    geom_line() +
    ylab(TeX("$IF_1(x_0)$") )+
    xlab(TeX("$x_0$") )+
    theme_minimal()+
    theme(axis.text = element_text(face="bold",size=15),
          legend.text = element_text(size=15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.position = c(0.8, 0.79))
  
  grafico2 <- ggplot(df, aes(x=x0,y=FF2))+
    geom_line() +
    ylab(TeX("$IF_2(x_0)$") )+
    xlab(TeX("$x_0$") )+
    theme_minimal()+
    theme(axis.text = element_text(face="bold",size=15),
          legend.text = element_text(size=15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.position = c(0.8, 0.79))
  
  grafico3 <- ggplot(data.frame(x = c(0,1)), aes(x = x)) +
    stat_function(fun = dbeta,n=1000,
                  args = list(shape1 = a, shape2=b)) +
    ylab(TeX("$f_{\\theta}(x)$") )+
    theme_minimal()+
    theme(axis.text = element_text(face="bold",size=15))
  grafico3+(grafico1/grafico2)
}




