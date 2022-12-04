cub_empirico<- function(a,b, x0=0, cont=0){
  #x0: contaminción
  #cont: porcentaje de contaminación
  set.seed(2022)
  MVa <-  0
  MVb <-  0
  
  MIa <-  0
  MIb <-  0
  
  for(i in 1:1000){
    muestra <- rbeta(100, shape1 = a, shape2 = b)
    if(cont!=0){
      muestra[1:cont]<- rep(x0,cont)
    }
    
    ajusteMV <- confint(fitdist(muestra,"beta"))
    ajusteMI <- Ajuste_MI_estimador_beta(muestra)[,3:4]
    
    MVa <- MVa + ifelse(ajusteMV[1,1]<=a & ajusteMV[1,2]>=a,1,0)
    MVb <- MVb + ifelse(ajusteMV[2,1]<=b & ajusteMV[2,2]>=b,1,0)
    
    MIa <- MIa + ifelse(ajusteMI[1,1]<=a & ajusteMI[1,2]>=a,1,0)
    MIb <- MIb + ifelse(ajusteMI[2,1]<=b & ajusteMI[2,2]>=b,1,0)
  }
  
  return(c(MIa,MIb,MVa,MVb)/1000)
  
}

cub_empirico(a=2, b=6, x0=0.5, cont = 3)

