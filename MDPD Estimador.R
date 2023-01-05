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

Halph <- function(x, muestra, alpha){
  (gamma(x[1]+x[2])/(gamma(x[1])*gamma(x[2])))^(1+alpha)*
    gamma(alpha*x[1]+x[1]-alpha)*gamma(alpha*x[2]+x[2]-alpha)/
    gamma(alpha*x[1]+alpha*x[2]-2*alpha+x[1]+x[2] ) - 
    (1+1/alpha)*mean(dbeta(muestra, shape1 = x[1], shape2 = x[2])^alpha)
}

MDPDE_estimador_beta <- function(x, alpha=0, rob=TRUE){
  
  if(rob){
    xbar <- median(x)
    s2   <- mad(x)^2
    alpha0 <- -xbar*(-xbar+xbar^2+s2)/s2
    beta0 <- (alpha0-alpha0*xbar)/xbar
    
    #ajuste <- robustbetareg(x~1, alpha = 0.5)
    #muhat  <- 1/(1+exp(-ajuste$coefficients$mean))
    #alpha0 <- as.numeric(muhat*ajuste$coefficients$precision)
    #beta0  <- as.numeric(ajuste$coefficients$precision - alpha0)
  }else{
    xbar <- mean(x)
    s2   <- var(x)
    alpha0 <- -xbar*(-xbar+xbar^2+s2)/s2
    beta0 <- (alpha0-alpha0*xbar)/xbar
  }

  est <- stats::optim(par = c(alpha0, beta0),
                      fn=Halph,
                      control = list(maxit=100),
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
 
  n <- length(muestra)
  zz <- qt(aa/2, lower.tail = F, df=n-1)
  
  # estimación puntual
  estimador <- MDPDE_estimador_beta(muestra, alpha, rob=FALSE)
  
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














