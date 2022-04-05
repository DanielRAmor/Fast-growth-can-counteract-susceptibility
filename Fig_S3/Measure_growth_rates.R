#README:
#This code reads the time series data for the growth of the 2 species in the presence of each of the 8 antibiotics at different concentrations
#The raw data contains the Optical Density (OD) of a 24hr incubation, with measurements every 15mins
#The algorithm prepares the data and measures the growth rate in each condition via fitting an exponential growth rate
#The growth rate fitting takes place within the OD range [0.01, 0.05], and within the first 12hr of growth
#The susceptibility is determined via the normalized OD (relative to the OD in the absence of antibiotic) at 24hr
#The output generates 2 output files (one per species) with the measured growth rates, as well as one plot for each measurement

library (ggplot2)
library (plyr)  
library(scales)
library(varhandle)

#read the data
G8all <- read.table ("Growth_Curves_in_Antibiotics.csv", header=T, sep = ",", dec = ".")
G8all <- G8all[, 2:ncol(G8all)]
G8all[G8all$OD<0, 5] = 0.00001 #remove potential negative values resulting from background OD correction (applied previously, right after reading the raw data)

# Use only the data of the first 12 hr
G8all = G8all[G8all[,3]< 12,]

# Reminder of antibiotic concentration in experimental plate column
#col_conc <- 1 # for 100 ug/ml (ACHTUNG! Divide by 4 for the case of Gentamicin and Kanamycin. Their maximum concentration was 25 ug/ml)
#col_conc <- 2 # for 50 ug/ml
#col_conc <- 3 # for 25 ug/ml
#col_conc <- 4 # for 12.5 ug/ml
#col_conc <- 5 # for 6.25 ug/ml
#col_conc <- 6 # for 3.125 ug/ml
#col_conc <- 7 # for 1.56 ug/ml
#col_conc <- 8 # for 0.78 ug/ml
#col_conc <- 9 # for 0.39 ug/ml
#col_conc <- 10 # for 0.195 ug/ml
#col_conc <- 11 # for 0.975 ug/ml
#col_conc <- 12 # for 0.0 ug/ml

#Reminder of antibiotic in experimental plate row
# A = Amp
# B = Car
# C = Cip
# D = Chl
# E = Ery
# F = Gen
# G = Kan
# H = Tet

#Plot aesthetics
background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_rect (colour="black", size=1.0), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12))
titulos <- labs(title="", x=expression(paste("Time (hr)")), y=expression(paste("OD")))
escalaY = scale_y_continuous(limits= c(0.0001 , 1), trans='log10',  labels = trans_format("log10", math_format(10^.x)))
escalaX = scale_x_continuous(limits = c(0, 23.8))
escalaX = scale_x_continuous(limits = c(0, 12))
lineas <- geom_line()
colores <- scale_colour_manual(values = c("#1AB24B", "#168936", "#24CE54", '#FAA51A', '#D88B1A', '#BA7D25', 'black', 'green', 'orange')) 
No_legend = theme(legend.position = "none")

## Prepare data frame for Lp
reportLp <- c()
count1 <- 0
outLp <- data.frame(c(100.1),c('D12'),c('R1'), c(0.0),c(0.0),c(0.0), stringsAsFactors=FALSE)
colnames(outLp) = c('Concentration','Well','Replicate','Slope', 'Intercept', 'Rsquare')
col_conc <- 0.0 # To store antibiotic concentration in experimental plate column
columna <- c(1,2,3,4,5,6,7,8,9,10,11,12) # Columns of experimental plate
fila <- c('A','B','C','D','E','F','G','H') # Rows of experimental plate
BioRep <- c('R1','R2','R3') # Replicate (different plates)
Especie <- c('Lp') 
#for loop to measure all growth rates for Lp and plot each result
for (m in 1:8){
  for (k in 1:12){
    for (j in 1:3){
      #choose one condition to analyze and plot
      Pozillo <- c(paste(fila[m], columna[k], sep = ""))
      G8 <- G8all[G8all$Well == Pozillo & G8all$Replicate == BioRep[j] & G8all$Species == Especie, ]
      count1 <- count1 + 1
      outLp[count1,] = c('col_conc', paste(G8[1,4]), G8[1,1], 0,0,0)
      ploti <- ggplot(data=G8, aes(x= Time_hr, y = OD))
      #fit an exponential growth rate
      if(max(G8$OD)>0.05){
        growth<-lm(log(G8[G8$OD>0.01 & G8$OD<0.05,5]) ~ G8[G8$OD>0.01 & G8$OD<0.05,3], data=G8)
        summary(growth)
        
        #store data of the fitting 
        reportLp <- summary(growth)
        outLp[count1,4] <-reportLp$coefficients[2, 'Estimate']
        outLp[count1,5] <-reportLp$coefficients[1, 'Estimate']
        outLp[count1,6] <-reportLp$r.squared
        alphaLp <-  growth$coefficients[[2]] 
        betaLp <- growth$coefficients[[1]] 
        #plot the slope of the measured (fitted) growth rate
        fitlineLp <- geom_segment(aes(x = 1.5, xend = 5.5, y = exp(betaLp + alphaLp*1.5 +0.3 ), yend = exp(betaLp + alphaLp*5.5 +0.3)))
        
        refer1 <- geom_hline(yintercept = 0.05) 
        refer2 <- geom_hline(yintercept = 0.01)    #these lines indicate the OD range considered in the fitting 
        
        plotii = ploti + background  + titulos + lineas + escalaY + escalaX + colores +fitlineLp + No_legend + refer1 + refer2
        plotii
        
      }else{
        outLp[count1,4] <- 0.0
        outLp[count1,5] <- 0.0
        outLp[count1,6] <- 0.0
        plotii = ploti + background  + titulos + lineas + escalaY + escalaX + colores + No_legend
        plotii
      }  
      
      pdf(paste(Especie,"-GrowthRate","-Replicate-",j,"-Well",G8[1,4], ".pdf"), width=4, height=3)
      print(plotii)
      dev.off()
      
    }
  }
}

#export the growth rates for Lp
write.csv(outLp,"Lp-GrowthRates-0.01to0.05-ODRange.csv", row.names = TRUE)

## Prepare data frame for Ca
reportCa <- c()
count1 <- 0
outCa <- data.frame(c(100.1),c('D12'),c('R1'), c(0.0),c(0.0),c(0.0), stringsAsFactors=FALSE)
colnames(outCa) = c('Concentration','Well','Replicate','Slope', 'Intercept', 'Rsquare')
col_conc <- 0.0 # To store antibiotic concentration in experimental plate column
columna <- c(1,2,3,4,5,6,7,8,9,10,11,12) # Columns of experimental plate
fila <- c('A','B','C','D','E','F','G','H') # Rows of experimental plate
BioRep <- c('R1','R2','R3') # Replicate (different plates)
Especie <- c('Ca') 
#for loop to measure all growth rates for Ca and plot each result
for (m in 1:8){
  for (k in 1:12){
    for (j in 1:3){
      #choose one condition to analyze and plot
      Pozillo <- c(paste(fila[m], columna[k], sep = ""))
      G8 <- G8all[G8all$Well == Pozillo & G8all$Replicate == BioRep[j] & G8all$Species == Especie, ]
      count1 <- count1 + 1
      outCa[count1,] = c('col_conc', paste(G8[1,4]), G8[1,1], 0,0,0)
      ploti <- ggplot(data=G8, aes(x= Time_hr, y = OD))
      #fit an exponential growth rate
      if(max(G8$OD)>0.05){
        growth<-lm(log(G8[G8$OD>0.01 & G8$OD<0.05,5]) ~ G8[G8$OD>0.01 & G8$OD<0.05,3], data=G8)
        summary(growth)
        #store data of the fitting 
        reportCa <- summary(growth)
        outCa[count1,4] <-reportCa$coefficients[2, 'Estimate']
        outCa[count1,5] <-reportCa$coefficients[1, 'Estimate']
        outCa[count1,6] <-reportCa$r.squared
        alphaCa <-  growth$coefficients[[2]] 
        betaCa <- growth$coefficients[[1]] 
        #plot the slope of the measured (fitted) growth rate
        fitlineCa <- geom_segment(aes(x = 1.5, xend = 5.5, y = exp(betaCa + alphaCa*1.5 +0.3 ), yend = exp(betaCa + alphaCa*5.5 +0.3)))
        
        refer1 <- geom_hline(yintercept = 0.05) 
        refer2 <- geom_hline(yintercept = 0.01)    #these lines indicate the OD range considered in the fitting 
        
        plotii = ploti + background  + titulos + lineas + escalaY + escalaX + colores +fitlineCa + No_legend + refer1 + refer2
        plotii
        
      }else{
        outCa[count1,4] <- 0.0
        outCa[count1,5] <- 0.0
        outCa[count1,6] <- 0.0
        plotii = ploti + background  + titulos + lineas + escalaY + escalaX + colores + No_legend
        plotii
      }  
      
      pdf(paste(  Especie,"-GrowthRate","-Replicate-",j,"-Well",G8[1,4], ".pdf"), width=4, height=3)
      print(plotii)
      dev.off()
      
    }
  }
}

#export the growth rates for Ca
write.csv(outCa,"Ca-GrowthRates-0.01to0.05-ODRange.csv", row.names = TRUE)

