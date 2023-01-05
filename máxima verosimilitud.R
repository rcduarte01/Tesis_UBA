library(ggplot2)
library(ggpubr)
library(gridExtra)
library(patchwork)
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
  
  
  df <- data.frame(x0,FF1=FF[1,],FF2 = FF[2,])
  
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
                  args = list(shape1 = alfa, shape2=beta)) +
    ylab(TeX("$f_{\\theta}(x)$") )+
    theme_minimal()+
    theme(axis.text = element_text(face="bold",size=15))
  grafico3 | (grafico1 / grafico2)
}

Finfluencia_MLE(7,7)
Finfluencia_MLE(0.7,0.9)

Finfluencia_MLE(3,0.8)
Finfluencia_MLE(2,6)
Finfluencia_MLE(0.8,0.8)





