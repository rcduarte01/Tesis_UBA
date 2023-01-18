Finfluencia_MI(2,6)
Finfluencia_MDPDE(2,6, 0.3)
Finfluencia_MLE(0.7, 0.9)

tt <- rbeta(100, shape1 = 2, shape2 = 9)
curve(psibeta1_MDPDE(x, muestra = tt, alpha = 0.5))

curve(pbeta(x, shape1 = 0.7, shape2 = 0.9)^4-1/5)

curve(dbeta(x,shape1 = 0.7, shape2 = 0.9))
