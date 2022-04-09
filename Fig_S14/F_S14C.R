library(deSolve)
library(ggplot2)

#FUNCTIONS
  #LV model with Allee effect
  dy = c(0.0,0.0)
  LVAMonod <- function (t,y,parameters) {
    x <- y[1]
    y <- y[2]
    r1 <- parameters[1]
    r2 <- parameters[2]
    alpha12 <- parameters[3]
    alpha21 <- parameters[4]
    A <- parameters[5]
    t1 = t
    ta = parameters[8]
    IC51 = parameters[9]
    IC52 = parameters[10]
    Cm = parameters[11]
    
    if(t1 < ta | t1 > (ta+24) ) {
      dy[1] <- r1*x *((1-x)*x/(A+x) -alpha12*y)  
      dy[2] <- r2*y *(1 -y -alpha21*x)  
    } else {  #antibiotic shock
      dy[1] <- Rstatic(r1, Cm, IC51) *x * ((1-x)*x/(A+x) -alpha12*y) 
      dy[2] <- Rstatic(r2, Cm, IC52)*y *(1 -y -alpha21*x) 
    }
    list(dy)
  }
  
  # Serial dilutions with dispersal
  Dilute <- function(t, y, parms){
    y[1] <- y[1] / parms[7] + parms[6] 
    y[2] <- y[2] / parms[7] + parms[6]
    return(c(y[1], y[2]))
  }
  
  # Bacteriostatic model
  Rstatic <- function (R0, Bstat, I50) {
    Rstat <- R0 / (1 + Bstat / I50 ) 
    return (Rstat)
  }

#PARAMETER VALUES
  stepT = 1/100
  Tf = 240
  dilu = 30 #dilution rate
  rf= 0.52 #growth rate fast grower
  rs=0.41 #growth rate slow grower
  D = 0.02  # Dispersal applied during dilutions
  #AA=0.0
  AA=0.0085 #Allee effect
  alpha12 = 1.6 #interaction strength
  alpha21 = 1.6 #interaction strength
  Tshock = 96 #time at which antibiotic exposure begins
  IC50_1 = 3.1 #IC50 of fast grower
  IC50_2 = 8.3 #IC50 of slow grower
  Cm0 = 40 # Antibiotic concentration during shock
  
  #in order: r1,r2,alpha12,alpha21,A,D, dil, Tshock, IC50_1, IC50_2, Cm0 
  pars = c(rf,rs,alpha12,alpha21, AA, D, dilu, Tshock, IC50_1, IC50_2, Cm0)

  #simulation for Fig S14_C (left panel)
  yini = c(0.98, 0.01) 
  LVA.simulation <- ode(y = yini, func = LVAMonod,
                       times = seq(0,Tf,by = stepT), parms = pars, events = list(func = Dilute, time = c(24,48,72,96,120,144,168,192,216,240,264)))
  #plot it
  plot (LVA.simulation[,1], log10(LVA.simulation[,2]), xlim=c(0, Tf), ylim=c(-4,0), xlab='Time', ylab='Density', type = 'l', col='green', lwd = 2) 
  abline(h=log10(D), lty=2)
  lines (LVA.simulation[,1], log10(LVA.simulation[,3]), col = 'orange', lwd = 2)
  #export plot
  pdf('F_S14C_left.pdf', width=4, height=4)
  plot (LVA.simulation[,1], log10(LVA.simulation[,2]), xlim=c(0, Tf), ylim=c(-4,0), xlab='Time', ylab='Density', type = 'l', col='green', lwd = 2)
  lines (LVA.simulation[,1], log10(LVA.simulation[,3]), col = 'orange', lwd = 2)
  abline(h=log10(D), lty=2)
  dev.off()



#simulation for Fig S14_C (right panel)
  yini = c(0.01, 0.98) 
  LVA.simulation <- ode(y = yini, func = LVAMonod,
                      times = seq(0,Tf,by = stepT), parms = pars, events = list(func = Dilute, time = c(24,48,72,96,120,144,168,192,216,240,264)))
  #plot it
  plot (LVA.simulation[,1], log10(LVA.simulation[,2]), xlim=c(0, Tf), ylim=c(-4,0), xlab='Time', ylab='Density', type = 'l', col='green', lwd = 2) 
  abline(h=log10(D), lty=2)
  lines (LVA.simulation[,1], log10(LVA.simulation[,3]), col = 'orange', lwd = 2)
  #export plot
  pdf('F_S14C_right.pdf', width=4, height=4)
    plot (LVA.simulation[,1], log10(LVA.simulation[,2]), xlim=c(0, Tf), ylim=c(-4,0), xlab='Time', ylab='Density', type = 'l', col='green', lwd = 2)
    lines (LVA.simulation[,1], log10(LVA.simulation[,3]), col = 'orange', lwd = 2)
    abline(h=log10(D), lty=2)
  dev.off()
