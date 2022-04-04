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
rs = 0.43 # growth rate of S
a12 = 1.25 #inhibition of F by S
a21 = 1.25 #inhibition of S by F
A = 0.0 #Allee effect
Dispersal = 0.0001 #Dispersal rate
Tshock = 72 # time at which antibiotic exposure begins
Death1A = 0.0 #antibiotic-associated death rate for F (set to zero here so that there is not exposure when generating the phase space)
Death2A = 0.0 #antibiotic-associated death rate for S (set to zero here so that there is not exposure when generating the phase space)
Kf <- 1.0
Ks <- 1.0
pars = c(rf,rs,a12,a21, A, Dispersal, Tshock,  Death1A, Death2A, Kf, Ks) #vector of parameter values


#data frame to store the results
a = c(1.0,1.0)
phaseX <- data.frame(a,a,a,a)
colnames(phaseX) = c('xini', 'yini', 'xf','yf')

#Simulations to generate the phase space begin here
#resolution of the phase space (Nsims * Nsims)
Nsims = 200
yini = c(1.0, 1.0)
for (i in 1:Nsims){
  #explore the axis for initial abundance of F
  if (i > 1) {yini[1] = yini[1] /(10^(4/Nsims))} #where 4 is the order of magnitudes that will be covered
  
  for (j in 1:Nsims) {
    if (j < 2) {yini[2] = 1.0
    } else {
      #explore the axis for initial abundance of F
      yini[2] = yini[2] /(10^(4/Nsims)) #where 4 is the order of magnitudes that will be covered
    }
    #run a simulation
    LV.simulation <- ode(y = yini, func = LVAMonod,
                         times = seq(0,Tf,by = stepT), parms = pars)
    #export the result of the simulation at the final time
    phaseX[(i-1)*Nsims + j,1] = yini[1]
    phaseX[(i-1)*Nsims + j,2] = yini[2]
    phaseX[(i-1)*Nsims + j,3] = LV.simulation[Tf/stepT,2]
    phaseX[(i-1)*Nsims + j,4] = LV.simulation[Tf/stepT,3]
    
    
  }
}

#compute the fraction of species F in the outcome of the simulations
phaseX[5] = phaseX[3] / (phaseX[3] + phaseX[4]) 
colnames(phaseX)[5] = c('Fraction1')
#to plot with log axes
phaseX[1] = log10(phaseX[1]) 
phaseX[2] = log10(phaseX[2]) 

#Draw trajectory arrows
  Narrows = 28
  #prepare a data frame for the trajectory arrows
  a = c(0.1,0.1,0.1,0.1)
  arrows = data.frame(a,a,a,a)
  colnames(arrows) = c('x0','xf','y0','yf')
  #soft perturbation
  Death1A = 0.6
  Death2A = 0.6/1.3
  pars = c(rf,rs,a12,a21, A, Dispersal, Tshock,  Death1A, Death2A, Kf, Ks)
    #starting close to one stable state
    yini = c(0.0005,0.9995)
    LV.simulation <- ode(y = yini, func = LVAMonod,
                         times = seq(0,Tf,by = stepT), parms = pars)
    for (i in 1:Narrows) {
      arrows[i,1] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*(i-1)+1,2])
      arrows[i,2] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*i,2])
      arrows[i,3] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*(i-1)+1,3])
      arrows[i,4] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*i,3])
    }
    arrows1 <- arrows
    #starting close to the other stable state
    yini = c(0.9995,0.0005)
    LV.simulation <- ode(y = yini, func = LVAMonod,
                         times = seq(0,Tf,by = stepT), parms = pars)
    for (i in 1:Narrows) {
      arrows[i,1] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*(i-1)+1,2])
      arrows[i,2] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*i,2])
      arrows[i,3] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*(i-1)+1,3])
      arrows[i,4] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*i,3])
    }
    arrows2 <- arrows
    
  # harsh perturbation
  Death1A = 0.9
  Death2A = 0.9/1.3
  pars = c(rf,rs,a12,a21, A, Dispersal, Tshock,  Death1A, Death2A, Kf, Ks)
    #starting close to one stable state
    yini = c(0.0005,0.9995)
    LV.simulation <- ode(y = yini, func = LVAMonod,
                         times = seq(0,Tf,by = stepT), parms = pars)
    for (i in 1:Narrows) {
      arrows[i,1] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*(i-1)+1,2])
      arrows[i,2] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*i,2])
      arrows[i,3] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*(i-1)+1,3])
      arrows[i,4] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*i,3])
    }
    arrows3 <- arrows
    #starting close to the other stable state
    yini = c(0.9995,0.0005)
    LV.simulation <- ode(y = yini, func = LVAMonod,
                         times = seq(0,Tf,by = stepT), parms = pars)
    for (i in 1:Narrows) {
      arrows[i,1] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*(i-1)+1,2])
      arrows[i,2] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*i,2])
      arrows[i,3] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*(i-1)+1,3])
      arrows[i,4] = log10(LV.simulation[floor(Tf/(stepT*Narrows))*i,3])
    }
    arrows4 <- arrows

ploti <- ggplot(data = phaseX, aes(x= xini, y= yini ))
background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_rect (colour="black", size=1.5), axis.text=element_text(size=18), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )
titulos<- labs(title="", x=expression(paste("Normalized abundance F")), y=expression(paste("Normalized abundance S")))
escalaX <- scale_x_continuous(expand= c(0,0), limits=c(-4,0))
escalaY <- scale_y_continuous(expand= c(0,0), limits=c(-4,0)) 
couleur <- scale_fill_gradient(low = "orange", high = "chartreuse4", limits=c(0.0000, 1.0), guide="colorbar")

rast <- geom_tile(aes(fill = Fraction1))
trajectory1 <- geom_segment(data=arrows1[-c(12:16),], mapping= aes(x=x0,y=y0,xend=xf,yend=yf), arrow=arrow(length = unit(2,"mm")), size=0.8, color='white')
trajectory1a <- geom_segment(data=arrows1[12:16,], mapping= aes(x=x0,y=y0,xend=xf,yend=yf), arrow=arrow(length = unit(2,"mm")), size=0.8, color='white', linetype= 'longdash', alpha = 0.7)
trajectory2 <- geom_segment(data=arrows2[-c(12:16),], mapping= aes(x=x0,y=y0,xend=xf,yend=yf), arrow=arrow(length = unit(2,"mm")), size=0.8, color='white')
trajectory2a <- geom_segment(data=arrows2[12:16,], mapping= aes(x=x0,y=y0,xend=xf,yend=yf), arrow=arrow(length = unit(2,"mm")), size=0.8, color='white', linetype= 'longdash', alpha = 0.7)
trajectory3 <- geom_segment(data=arrows3[-c(12:16),], mapping= aes(x=x0,y=y0,xend=xf,yend=yf), arrow=arrow(length = unit(2,"mm")), size=0.8, color='gray2')
trajectory3a <- geom_segment(data=arrows3[12:16,], mapping= aes(x=x0,y=y0,xend=xf,yend=yf), arrow=arrow(length = unit(2,"mm")), size=0.8, color='gray2', linetype= 'longdash', alpha = 0.7)
trajectory4 <- geom_segment(data=arrows4[-c(12:16),], mapping= aes(x=x0,y=y0,xend=xf,yend=yf), arrow=arrow(length = unit(2,"mm")), size=0.8, color='gray2')
trajectory4a <- geom_segment(data=arrows4[12:16,], mapping= aes(x=x0,y=y0,xend=xf,yend=yf), arrow=arrow(length = unit(2,"mm")), size=0.8, color='gray2', linetype= 'longdash', alpha = 0.7)

plotii<- ploti + background + titulos + couleur + escalaX + escalaY  + rast + trajectory1 + trajectory1a + trajectory2 + trajectory2a + trajectory3 + trajectory3a + trajectory4 + trajectory4a
plotii

pdf(paste("F_2C.pdf"), width=6, height=4.5)
print(plotii)
dev.off()

