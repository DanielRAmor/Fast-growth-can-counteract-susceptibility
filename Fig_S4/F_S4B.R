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
Tf = 168  #Simulated hours
stepT = 1/30 # time resolution of the output
rf = 0.54 # growth rate of F
rs = 0.43 # growth rate of S
a12 = 1.25 #inhibition of F by S
a21 = 1.25 #inhibition of S by F
A = 0.0 #Allee effect (set to zero for No Allee effect)
Dispersal = 0.001 #Dispersal rate
Tshock = 72 # time at which antibiotic exposure begins
Death1A = 0.0 #no exposure to antibiotic
Death2A = 0.0  
Kf <- 1.0# Carrying capacity fast grower
Ks <- 1.0 # Carrying capacity slow grower

#set range of Kf and Ks to explore
KMax <- 3.0
KMin <- 0.3
Kf <- KMax
Ks <- KMax

#prepare data frame to store data
a = c(1.0,1.0)
phaseX <- data.frame(a,a,a,a)
colnames(phaseX) = c('Kf', 'Ks', 'xf','yf')


#begin simulations to plot the phase space
yini = c(0.1, 0.1)
Nsims = 30 #resolution of each axis of the phase space (number of simulations in each direction of the phase space)
#'for' loop to explore the range of alpha12 and Ks
for (i in 1:Nsims){
  print(i) #to report progress in the 'for' loop
  if (i > 1) {Kf = Kf - (KMax-KMin)/(Nsims-1)}  
  
  for (j in 1:Nsims) {
    if (j < 2) {Ks = KMax
    } else {
      Ks = Ks - (KMax-KMin)/(Nsims-1)  
    }
    
    pars = c(rf,rs,a12,a21, A, Dispersal, Tshock,  Death1A, Death2A, Kf, Ks)
    LV.simulation <- ode(y = yini, func = LVAMonod,
                         times = seq(0,Tf,by = stepT), parms = pars)#, events= list(func = Shock, time = Tf/2))
    
    phaseX[(i-1)*Nsims + j,1] = Kf
    phaseX[(i-1)*Nsims + j,2] = Ks
    phaseX[(i-1)*Nsims + j,3] = LV.simulation[Tf/stepT,2]
    phaseX[(i-1)*Nsims + j,4] = LV.simulation[Tf/stepT,3]
    
    
  }
}

phaseX[5] = phaseX[3] / (phaseX[3] + phaseX[4]) 
colnames(phaseX)[5] = c('Fraction1')

#plot the results
ploti <- ggplot(data = phaseX, aes(x= Kf, y= Ks ))
background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_rect (colour="black", size=1.5), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22))
titulos<- labs(title="", x=expression(paste("Carrying Capacity of F")), y=expression(paste("Carrying Capacity of S")))
escalaX <- scale_x_continuous(expand= c(0,0), limits=c(KMin,KMax))
escalaY <- scale_y_continuous(expand= c(0,0), limits=c(KMin,KMax)) 
couleur <- scale_fill_gradient(low = "orange", high = "chartreuse4", limits=c(0.0000, 1.0), guide="colorbar")
rast <- geom_tile(aes(fill = Fraction1))
plotii<- ploti + rast + background  + titulos + couleur + escalaX + escalaY 
plotii
#export plot
pdf(paste("F_S4B.pdf"), width=6, height=4.5)
print(plotii)
dev.off()
