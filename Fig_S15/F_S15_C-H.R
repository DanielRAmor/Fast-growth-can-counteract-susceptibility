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
  deltaMax = 1.0 #maximum antibiotic associated death rate
  D = Dispersal(10^(-3))
  mu <- 0.2 #maximum interaction strength within groups (or minimum between groups)
  Kmax <- 1.0 #maximum value in the distribution of carrying capacities
  
#RUN MULTIPLE SIMULATIONS
  Nsims <- 100 #number of simulated communities
  #prepare data frames to store the data
  AltStateA <- data.frame(matrix(ncol=2, nrow= Nsims)) #reports if State A can stably dominate the community
  AltStateB <- data.frame(matrix(ncol=2, nrow= Nsims)) #reports if State B can stably dominate the community
  colnames(AltStateA)= c('InitA', 'InitB')
  colnames(AltStateB)= c('InitA', 'InitB')
  
  ShockToA <- data.frame(matrix(ncol=12, nrow= 1))  #reports if State A is resilient to antibiotic exposure
  ShockToB <- data.frame(matrix(ncol=12, nrow= 1))  #reports if State B is resilient to antibiotic exposure
  
  AltStateA[,]= 0
  AltStateB[,]= 0
  ShockToA[,]= 0
  ShockToB[,]= 0
  
  colnames(ShockToA)= c('ShockToA','ShockToB', 'r1', 'r11', 'Susc1', 'Susc11', 'MeanRGroupA','MeanRGroupB','MeanSuscA','MeanSuscB','MeanKA', 'MeanKB')
  colnames(ShockToB)= c('ShockToA','ShockToB', 'r1', 'r11', 'Susc1', 'Susc11', 'MeanRGroupA','MeanRGroupB','MeanSuscA','MeanSuscB','MeanKA', 'MeanKB')
  
  CountAltStates <- c() #to track which communities display bistability
  
  Communities = data.frame(matrix(ncol=(4+Nspecies), nrow= (Nspecies*Nsims))) #store community matrices, growth rates, susceptibilities and carrying capacities of all the simulated communities
  colnames(Communities)= c('Species', 1:Nspecies,'GrowthRate','Susceptibility','CarryingCapacity') 
  
  #'for' loop to explore bistability and resilience
  for (w in 1:Nsims){
    print(w) #just to report iteration number in the for loop
    r = GrowthRates(1.0) #set random growth rates for this community
    Ks <- Kapacities(Kmax, Nspecies) #set carrying capacities for these communities
    Ks[] <- 1.0 #comment this line if exploring the case for non-normalized carrying capacities
    A = ComMatrix(mu, Nspecies, Ks) #set the interaction matrix
    Susc = Susceptibilities(deltaMax) #set random susceptibilities
    
    Communities[(1+(w-1)*20):(20+(w-1)*20),1] = c(1:Nspecies) # Communities[] saves the data of each simulated interaction matrix, with associated growth rates and death rates
    Communities[(1+(w-1)*20):(20+(w-1)*20),2:21] = A
    Communities[(1+(w-1)*20):(20+(w-1)*20),22] = r[]
    Communities[(1+(w-1)*20):(20+(w-1)*20),23] = Susc[]
    Communities[(1+(w-1)*20):(20+(w-1)*20),24] = Ks[]
    #first we assess if this community is bistable (no antibitotic exposure)
    ShockT = 15000 #set the time to start of antibiotic exposure anytime after the end of the simulation so that the community is not exposed to the antibiotic  
    #list of parameter values
    parameters <- list(r=r,A=A, D=D, Susc=Susc, ShockTi = ShockT)
    
    #Evaluate stability of State A
    init <- InitA(Nspecies)
    out <- ode(y = init, times = times, func = multi, parms = parameters, method = 'ode45') #numerically solves the modified LV equations
    AbundanceGroupA = mean(c(out[nrow(out),2:11])) #mean final abundance of Group A species
    AbundanceGroupB = mean(c(out[nrow(out),12:21]))
    if (AbundanceGroupA > 0.1) {AltStateA[w,1] = 1} # records the (boolean) outcome for Group A 
    if (AbundanceGroupB > 0.1) {AltStateB[w,1] = 1} #records the (boolean) outcome for Group B
    
    #Run Stable State 2
    init <- InitB(Nspecies)
    out <- ode(y = init, times = times, func = multi, parms = parameters, method = 'ode45')
    AbundanceGroupA = mean(c(out[nrow(out),2:11])) #mean final abundance of Group A species
    AbundanceGroupB = mean(c(out[nrow(out),12:21]))
    if (AbundanceGroupA > 0.1) {AltStateA[w,2] = 1}  # records the (boolean) outcome for Group A 
    if (AbundanceGroupB > 0.1) {AltStateB[w,2] = 1}  # records the (boolean) outcome for Group A 
    
    
    # Apply perturbations ONLY TO COMMUNITIES EXHIBITING BISTABILITY
     if (AltStateA[w,1] != AltStateA[w,2]) {
       CountAltStates <- c(CountAltStates, w) #if the community is bistable, add it to this list
       ShockT = 1000 #antibiotic exposure begins at t=1000 hr
       parameters <- list(r=r,A=A, D=D, Susc=Susc, ShockTi = ShockT)
       #begin simulation in stable state A
       init <- InitA(Nspecies)
       out <- ode(y = init, times = times, func = multi, parms = parameters, method = 'ode45') #numerically solves the LV equations
       AbundanceGroupA = mean(c(out[nrow(out),2:11])) #mean final abundance of Group A species
       AbundanceGroupB = mean(c(out[nrow(out),12:21]))
       if (AbundanceGroupA > 0.1) {ShockToA[length(CountAltStates),1] = 1} # records the outcome for State A
       if (AbundanceGroupB > 0.1) {ShockToB[length(CountAltStates),1] = 1} #records the outcome for State B
       #begin simulation in stable state B
       init <- InitB(Nspecies)
       out <- ode(y = init, times = times, func = multi, parms = parameters, method = 'ode45')
       AbundanceGroupA = mean(c(out[nrow(out),2:11])) #mean final abundance of Group A species
       AbundanceGroupB = mean(c(out[nrow(out),12:21]))
       if (AbundanceGroupA > 0.1) {ShockToA[length(CountAltStates),2] = 1} 
       if (AbundanceGroupB > 0.1) {ShockToB[length(CountAltStates),2] = 1} 
       
       ShockToA[length(CountAltStates),3] = r[1] # stores the growth rate of species 1 (not needed, but useful if we would like to explore specific parameter values to this particular species in group A)
       ShockToA[length(CountAltStates),4] = r[11]# stores the growth rate of species 11 (not needed, but useful if we would like to explore specific parameter values to this particular species in group B)
       ShockToA[length(CountAltStates),5] = Susc[1] # stores the susceptibility of species 1
       ShockToA[length(CountAltStates),6] = Susc[11] # stores the susceptibility of species 11     
       ShockToA[length(CountAltStates),7] = mean(r[1:10])
       ShockToA[length(CountAltStates),8] = mean(r[11:20])
       ShockToA[length(CountAltStates),9] = mean(Susc[1:10])
       ShockToA[length(CountAltStates),10] = mean(Susc[11:20])
       }
   }
  
  ShockToA[is.na(ShockToA[,1]),1] = 0 # Change NA values to 0
  ShockToA[is.na(ShockToA[,2]),2] = 0 

  #export the results (ShockToA is the output that contains more information and will be used for plotting)
  write.csv(ShockToA,paste('mu-',mu, '-Kmax-Rmax-deltaMax-', deltaMax,'-D-', D[1],'-ShockToA.csv'), row.names = FALSE)
  write.csv(ShockToB,paste('mu-',mu, '-Kmax-Rmax-deltaMax-', deltaMax,'-D-', D[1],'-ShockToB.csv'), row.names = FALSE)
  write.csv(AltStateA,paste('mu-',mu, '-Kmax-Rmax-deltaMax-', deltaMax,'-D-', D[1],'-AltStateA.csv'), row.names = FALSE)
  write.csv(AltStateB,paste('mu-',mu, '-Kmax-Rmax-deltaMax-', deltaMax,'-D-', D[1],'-AltStateB.csv'), row.names = FALSE)
  write.csv(Communities,paste('mu-',mu, '-Kmax-Rmax-deltaMax-', deltaMax,'-D-', D[1],'-Communities.csv'), row.names = FALSE)
  
  ShockToA[11]= ShockToA[4]/ShockToA[3] # r11/r1 (growth rate ratio between species 11 and 1)
  ShockToA[12]= ShockToA[6]/ShockToA[5] # delta11/delta1
  ShockToA[13]= ShockToA[8]/ShockToA[7] # (mean growth rate in B) / (mean growth rate in A)
  ShockToA[14]= ShockToA[10]/ShockToA[9] # mean(susceptibilityB)/mean(susceptibilityA)
  
  # report result as: transition to A = '2', no transition (or transition in both cases) = '1', or transition to B '0'
  ShockToA[15]= ShockToA[1] + ShockToA[2] #if this quantity =1 it can mean either no transitions, or transitions in both directions (this case was not found in the data generating Fig S_15)
  for (i in 1:nrow(ShockToA)){
    if (ShockToA[i,15] == 1){
      if(ShockToA[i,1] == 0) {print(c('ACHTUNG: a transition in both directions!'))} #reports the strange case of a transition in both directions
    }     
  }
  ShockToA[,15] = as.integer(ShockToA[,15])
  colnames(ShockToA)[11:15] = c('GrowthRatioFocals','SuscRatioFocals','MeanGrowthBToA','MeanSuscBToA', 'Transitions')
  
  #Plot as function of the ratios for average growth rates and average susceptibilities of the two groups
  ploti <- ggplot(data=ShockToA, aes(x= MeanGrowthBToA, y = MeanSuscBToA, color=as.factor(Transitions), fill=as.factor(Transitions)))#, group= interaction(variable), color= interaction (variable)))#, fill = Replicate))
  background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_blank (), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) #, axis.title.x=element_blank() )
  colores <- scale_color_manual(values = c("red", "gray", "blue"))
  puntos = geom_point()
  rast <- geom_tile(aes(fill = Transitions))
  titulos <- labs(title=expression(paste("")), x=expression(paste("Time")), y=expression(paste("Abundance")))
  slope45 <- geom_segment(x = 0, y = 0, xend = 2, yend = 2, linetype="dashed", color = "gray")
  plotii = ploti + background + puntos + colores + slope45 #+ horz + vert# + xscale + yscale
  plotii
  
  #export the plot
  pdf(paste(  "Resilience of bistable communities", ".pdf"), width=8, height=6)
  print(plotii)
  dev.off()
  
  
  

  
