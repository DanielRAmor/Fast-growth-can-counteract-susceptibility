library(deSolve)
library (ggplot2)
library(viridis)
library(reshape2)

#FUNCTIONS
  Nspecies = 20 #number of species in the community. Used in the function definitions below to create matrices and vectors accordingly.
  
  #Species growth rates 
  GrowthRates <- function(gr){
    r = runif(Nspecies) #random growth rate for each species
    return (r)
  }
  
  #Antibiotic susceptibilities
  Susceptibilities <- function(dMax){
    Susc = dMax* runif(Nspecies)
    return (Susc)
  }
  
  #Carrying capacities
  Kapacities <- function(KMax, Nspecies){
    Kap = KMax* runif(Nspecies) #random carrying capacities (below there is an option to consider normalized carrying capacities)
    return (Kap)
  }
  
  #Species Interactions
  ComMatrix <- function (mu, Nspecies, Ks){
        A = matrix( 2*mu*runif(Nspecies^2), nrow = Nspecies, ncol = Nspecies, byrow = TRUE) #Interaction matrix from uniform distribution with mean interaction 'mu'
    # bias the interaction matrix to have two antagonistic groups
    for (i in 1:Nspecies){
      for (j in 1:Nspecies){
        if (i< (1+Nspecies/2) & j < (1+Nspecies/2)) {
          A[i,j] = mu*runif (1) / Ks[i]   # Weak interactions within Group A
        }else{
          if (i > (Nspecies/2) & j > (Nspecies/2)) {
            A[i,j] = mu*runif (1) / Ks[i]  # Weak interactions within Group B
          }
          else{
            A[i,j] = mu+mu*runif (1) / Ks[i]  # Strong competition between Group A and B
          }
          
        }
      }
    }
    
    for (i in 1:Nspecies) {A[i,i]= 1.0 / Ks[i]}
    return (A)
  }
  
  #Dispersal
  Dispersal <- function (dd){ 
    D = runif(Nspecies)
    D[] = dd #case in which dispersal is a constant equal for all species
    return(D)
  }
  
  # Initial Species Abundances
  init <- runif(Nspecies) #create a vector for initial species abundances
  #use InitA to begin near the stable state dominated by group A
  InitA <- function (Nspecies) {  #starting from high abundance of group A, looking for Stable State A
    for (i in 1:Nspecies){
      if (i< (1+Nspecies/2) ){
        init[i]= 0.1  
      }else{
        init[i]= 0.001
      }
    }
    return (init[]) 
  }
  #use InitB to begin near the stable state dominated by group B
  InitB <- function(Nspecies){
    for (i in 1:Nspecies){
      if (i> (Nspecies/2) ){ #starting from high abundance of group B, looking for Stable State B
        init[i]= 0.1  
      }else{
        init[i]= 0.001
      }
    }
    return (init[]) 
  }
  
  #To numerically solve the LV model with dispersal and temporary exposure to antibiotics
  multi <- function(t, n, parms) {
    with(as.list(c(n,parms)),{
      if (t >= ShockTi &  t < (ShockTi+24)){
        delta <- Susc[]
      }else{
        delta <- 0
      }
      dn <- r * n * (1 - (A %*% n)) + D - delta*n 
      return(list(c(dn)))
    })
  }


#PARAMETER VALUES
  times <- seq(from = 0, to = 2000, by = 1) #time range with time step for the output of the numerical integration
  init <- runif(Nspecies) #create a vector for the initial species abundances (actual abundance values modified below)
  deltaMax = 1.0 #maximum antibiotic associated death rate
  D = Dispersal(10^(-3))

#TO RUN A SIMULATION FOR A SINGLE COMMUNITY: 
  r = GrowthRates(1.0) #random growth rates
  Kmax <- 1.0
  Ks <- Kapacities(Kmax, Nspecies)
  Ks[] <- 1.0 #comment / uncomment this line depending on whether different / normalized carrying capacities are being considered, respectively:
  mu <- 0.2 #maximum interaction strength within groups (or minimum between groups)
  A = ComMatrix(mu, Nspecies,Ks) #generate the matrix of interactions
  Susc = Susceptibilities(deltaMax) #random species susceptibilities
  ShockT = 1000
  parameters <- list(r=r,A=A, D=D, Susc=Susc, ShockTi = ShockT)
  #Initial abundances close to Stable State A
  init <- InitA(Nspecies)
  # OR Initial abundances close to Stable State B
  init <- InitB(Nspecies) #comment this line to initiate the simulation in stable state dominated by group A
  out <- ode(y = init, times = times, func = multi, parms = parameters, method = 'ode45')

#PLOT THE TIME SERIES FOR THE SPECIES ABUNDANCES OF THE SIMULATED COMMUNITY  
  #prepare the data
  z1 <- as.data.frame(out)
  z2 <- melt(z1, id.vars='time')
  
  #Some aesthetics
  background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) #, axis.title.x=element_blank() )
  escalaY = scale_y_continuous( trans='log10', limits= c(0.001 , 1))#, ##,  labels = trans_format("log10", math_format(10^.x)))
  escalaX = scale_x_continuous(limits = c(900, 1200))
  lineas <- geom_line()
  colores <- scale_color_viridis(discrete=T, option = 'H')
  No_legend = theme(legend.position = "none")
  titulos <- labs(title=expression(paste("")), x=expression(paste("Time")), y=expression(paste("Abundance")))
  
  #Make the plot
  ploti <- ggplot(data=z2, aes(x= time, y = value, group= interaction(variable), color= interaction (variable)))#, fill = Replicate))
  plotii = ploti + background  + titulos  + lineas + No_legend + escalaY + colores + escalaX #+ virus
  plotii
  
  #Export the plot
  pdf(paste(  "CommunityTimeSeries", ".pdf"), width=8, height=6)
  print(plotii)
  dev.off()
  
