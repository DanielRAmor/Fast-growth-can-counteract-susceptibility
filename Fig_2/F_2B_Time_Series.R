library(deSolve)
library(ggplot2)
library(phaseR)
library(rootSolve)

## The function to numerically solve the LV model:
dy = c(0.0,0.0)
LVAMonod <- function (t,y,parameters) {
  x <- y[1]
  y <- y[2]
  r1 <- parameters[1]
  r2 <- parameters[2]
  alpha12 <- parameters[3]
  alpha21 <- parameters[4]
  A <- parameters[5]
  M <- parameters[6]
  t1 = t
  ta = parameters[7]
  da1 = parameters[8]
  da2 = parameters[9]
  K1 = parameters[10]
  K2 = parameters[11]
  
  if(t1 < ta | t1 > (ta+24) ) {
    dy[1] <- r1*x *((x/(A+x)) * (1 -x/K1) -alpha12*y/K1) + M # dynamics of species F (fast grower, Ca)
    dy[2] <- r2*y *(1 -y/K2 -alpha21*x/K2) + M # dynamics of species S (slow grower, Lp)
  } else {  # temporary antibiotic shock
    dy[1] <- r1*x *((x/(A+x)) * (1 -x/K1) -alpha12*y/K1) + M - da1*x
    dy[2] <- r2*y *(1 -y/K2 -alpha21*x/K2) + M  - da2*y
  }
  list(dy)
}


#Parameter values
  #Change the 'rs' and 'Death2A' values to produce a time series for the
  #desired combination of growth rates ratio and death rates ratio
  Tf = 168  #Simulated hours
  stepT = 1/100 # time resolution of the output
  yini = c(0.9, 0.1) #initial normalized abundances of species R and S, respectively
  rf = 0.54 # growth rate of F
  rs = 0.43 # growth rate of S
  a12 = 1.25 #inhibition of F by S
  a21 = 1.25 #inhibition of S by F
  A = 0.0 #Allee effect
  Dispersal = 0.0001 #Dispersal rate
  Tshock = 72 # time at which antibiotic exposure begins
  Death1A = 0.9 #antibiotic-associated death rate for F (set to zero here so that there is not exposure when generating the phase space)
  Death2A = 0.9*0.4 #antibiotic-associated death rate for S (set to zero here so that there is not exposure when generating the phase space)
  Kf <- 1.0
  Ks <- 1.0
  pars = c(rf,rs,a12,a21, A, Dispersal, Tshock,  Death1A, Death2A, Kf, Ks) #vector of parameter values

###!!!

#Choose an initial species fraction
yini = c(0.9, 0.1) 
#yini = c(0.1, 0.9) 

#run a simulation and plot the time series of species abundances
LV.simulation <- ode(y = yini, func = LVAMonod,
                       times = seq(0,Tf,by = stepT), parms = pars)#, events= list(func = Shock, time = Tf/2))
#Note that vertical axis is in log scale:
plot (LV.simulation[,1], log10(LV.simulation[,2]), xlim=c(0, Tf), ylim=c(-4.5,0), xlab='Time', ylab='Abundance', type = 'l', col='green', lwd = 2)
lines (LV.simulation[,1], log10(LV.simulation[,3]), col = 'orange', lwd = 2)

#export pdf
pdf('rs0.43-delta0.4relativetoRf-InitiallyFastDominatedState.pdf', width=4, height=4)
  plot (LV.simulation[,1], log10(LV.simulation[,2]), xlim=c(0, Tf), ylim=c(-4.5,0), xlab='Time', ylab='Abundance', type = 'l', col='green', lwd = 2)
  lines (LV.simulation[,1], log10(LV.simulation[,3]), col = 'orange', lwd = 2)
dev.off()
