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
yini = c(0.9, 0.1) #initial normalized abundances of species R and S, respectively
rf = 0.54 # growth rate of F
rs = rf # growth rate of S
a12 = 1.25 #inhibition of F by S
a21 = 1.25 #inhibition of S by F
A = 0.0 #Allee effect (set to zero for No Allee effect)
Dispersal = 0.0001 #Dispersal rate
Tshock = 72 # time at which antibiotic exposure begins
Death1A = 0.9 #antibiotic-associated death rate for F (set to 0.9 as the maximum death rate considered in the phase diagram in Fig. 2B)
Death2A = Death1A  #antibiotic-associated death rate for S (sset to 0.9 as the maximum death rate considered in the phase diagram in Fig. 2B)
Kf <- 1.0
Ks <- 1.0
pars = c(rf,rs,a12,a21, A, Dispersal, Tshock,  Death1A, Death2A, Kf, Ks) #vector of parameter values


#data frame to store the results
a = c(1.0,1.0)
phaseX <- data.frame(a,a,a,a)
colnames(phaseX) = c('g_ratio', 'd_ratio', 'xf','yf')

#Simulations to generate the phase space begin here
#First we compute the outcome when F initially dominates the community (first stable state)

#resolution of the phase space (Nsims * Nsims)
Nsims = 300
countSims <- 0 #just a counter for number of simulations
for (i in 1:Nsims){
  countSims <- countSims +1
  print(countSims)
  #explore the axis of growth rates ratio
  if (i > 1) {pars[2] = pars[2] - rs/Nsims}
  for (j in 1:Nsims) {
    #explore the axis of death rates ratio
    if (j < 2) {pars[9] = Death2A
    } else {
      pars[9] = pars[9] - Death2A/Nsims
    }
    #run a simulation
    LV.simulation <- ode(y = yini, func = LVAMonod,
                         times = seq(0,Tf,by = stepT), parms = pars)
    #export the result of the simulation at the final time
    phaseX[(i-1)*Nsims + j,1] = pars[2]/pars[1]
    phaseX[(i-1)*Nsims + j,2] = pars[9]/pars[8]
    phaseX[(i-1)*Nsims + j,3] = LV.simulation[Tf/stepT,2]
    phaseX[(i-1)*Nsims + j,4] = LV.simulation[Tf/stepT,3]
    
    
  }
}

#compute the fraction of species F in the outcome of the simulations
phaseX[5] = phaseX[3] / (phaseX[3] + phaseX[4]) 
colnames(phaseX)[5] = c('Fraction1')

#to see how the outcome looks when starting from the stable state dominated by F is exposed to the antibiotic:
ploti <- ggplot(data = phaseX, aes(x= g_ratio, y= d_ratio))
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
titulos<- labs(title="", x=expression(paste("Log(Growth rate ratio)")), y=expression(paste("Log(Death rate ratio)")))
escalaX <- scale_x_continuous(expand= c(0,0), limits=c(-2,0))
escalaY <- scale_y_continuous(expand= c(0,0), limits=c(-2,0))  
couleur <- scale_fill_gradient(low = "#FAA421", high = "#1AB24B", limits=c(0.00, 1.0), guide="colorbar")
rast <- geom_tile(aes(fill = Fraction1))
plotii<- ploti + rast + background  + titulos + couleur
plotii

#in case we want to export the results
#pdf(paste("F-dominated_state_exposed_to_antibiotics.pdf"), width=6, height=4.5)
#print(plotii)
#dev.off()

phaseY = phaseX   # store this information before we repeat the process starting from the S-dominated stable state

#Now run simulations starting from the S-dominated stable state
yini = c(0.1, 0.9) 
#reinitialize parameter values
pars = c(rf,rs,a12,a21, A, Dispersal, Tshock,  Death1A, Death2A, Kf, Ks) #vector of parameter values
#reinitialize data frame to store the results
a = c(1.0,1.0)
phaseX <- data.frame(a,a,a,a)
colnames(phaseX) = c('g_ratio', 'd_ratio', 'xf','yf')

#Simulations to generate the phase space begin here
#Now we compute the outcome when S initially dominates the community (first stable state)
countSims <- 0 #just a counter for number of simulations
for (i in 1:Nsims){
  countSims <- countSims +1
  print(countSims)
  #explore the axis of growth rates ratio
  if (i > 1) {pars[2] = pars[2] - rs/Nsims}
  for (j in 1:Nsims) {
    #explore the axis of death rates ratio
    if (j < 2) {pars[9] = Death2A
    } else {
      pars[9] = pars[9] - Death2A/Nsims
    }
    #run a simulation
    LV.simulation <- ode(y = yini, func = LVAMonod,
                         times = seq(0,Tf,by = stepT), parms = pars)
    #export the result of the simulation at the final time
    phaseX[(i-1)*Nsims + j,1] = pars[2]/pars[1]
    phaseX[(i-1)*Nsims + j,2] = pars[9]/pars[8]
    phaseX[(i-1)*Nsims + j,3] = LV.simulation[Tf/stepT,2]
    phaseX[(i-1)*Nsims + j,4] = LV.simulation[Tf/stepT,3]
    
    
  }
}

#compute the fraction of species F in the outcome of the simulations
phaseX[5] = phaseX[3] / (phaseX[3] + phaseX[4]) 
colnames(phaseX)[5] = c('Fraction1')

#to see how the outcome looks when starting from the stable state dominated by S is exposed to the antibiotic:
ploti <- ggplot(data = phaseX, aes(x= g_ratio, y= d_ratio))
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
titulos<- labs(title="", x=expression(paste("Log(Growth rate ratio)")), y=expression(paste("Log(Death rate ratio)")))
escalaX <- scale_x_continuous(expand= c(0,0), limits=c(-2,0))
escalaY <- scale_y_continuous(expand= c(0,0), limits=c(-2,0))  
couleur <- scale_fill_gradient(low = "#FAA421", high = "#1AB24B", limits=c(0.00, 1.0), guide="colorbar")
rast <- geom_tile(aes(fill = Fraction1))
plotii<- ploti + rast + background  + titulos + couleur
plotii

#in case we want to export the results
#pdf(paste("S-dominated_state_exposed_to_antibiotics.pdf"), width=6, height=4.5)
#print(plotii)
#dev.off()

#now combine the data from the two phase spaces of resilience (one for each stable state) into one phase diagram
#first convert the final fraction of species into a Boolean variable (to assign in which of the two stable states the system is at the end of each simulation)
for (i in 1:nrow(phaseX) ) {
  if (phaseX[i,5] < 0.5) { phaseX[i,6] = 0.0} 
  else {phaseX[i,6] = 1.0}
  if (phaseY[i,5] < 0.5) { phaseY[i,6] = 0.0} 
  else {phaseY[i,6] = 1.0}
}
phaseX[7] = phaseX[6] + phaseY[6] #combine information of the final state,
                                  #the outcome will be 0, 1, or 2 for each point in the phase space
colnames(phaseX)[7] <- c('Resilience')
#plot the phase diagram of resilience
ploti <- ggplot(data = phaseX, aes(x= g_ratio, y= d_ratio )) 
#some plot aesthetics
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22), legend.position = "none") #, axis.title.x=element_blank() )
titulos<- labs(title="", x=expression(paste("Log(Growth rate ratio)")), y=expression(paste("Log(Death rate ratio)")))
rast <- geom_tile(aes(fill = Resilience))
couleur <- scale_fill_gradient(low = "#FAA421", high = "#1AB24B", limits=c(0.00, 2.0), guide="colorbar")
plotii<- ploti + rast + background  + titulos + couleur 
plotii

pdf(paste("F_2B.pdf"), width=4.5, height=4.5)
print(plotii)
dev.off()

