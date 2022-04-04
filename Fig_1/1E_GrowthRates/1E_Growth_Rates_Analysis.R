library (ggplot2)
library (plyr)  
library(scales)
library(varhandle)

#This code reads the data from 6 experimental plates (3 for Ca, 3 for Lp)
#Column 12 in each experimental plate contains populations in the absence of antibiotics
#The raw data contains the Optical Density (OD) of a 24hr incubation, with measurements every 15mins

#first, it generates the 'NoAntibiotic_GrowthCurves.csv' 
#then it computes the individual growth rate for each replicate and exports the results as
#an individual plot for each replicate, a table of individual growth rates for each species,
#and the table ' Mean_Growth_Rates_and_Intercepts.csv' for the mean growth rate and mean intercept resulting from the fitting algorithm

#READ THE DATA OF THE 6 PLATES
#Ca plates
G1 <- read.table ("R1_Ca_growth_curves.dat", header=F, sep = "\t", dec = ".") 
G2 <- read.table ("R2_Ca_growth_curves.dat", header=F, sep = "\t", dec = ".")
G3 <- read.table ("R3_Ca_growth_curves.dat", header=F, sep = "\t", dec = ".")

#Lp plates
G5 <- read.table ("R1_Lp_growth_curves.dat", header=F, sep = "\t", dec = ".")
G6 <- read.table ("R2_Lp_growth_curves.dat", header=F, sep = "\t", dec = ".")
G7 <- read.table ("R3_Lp_growth_curves.dat", header=F, sep = "\t", dec = ".")

#prepare the format of the data
rownames(G1) <- G1[,1]
G1 <- G1[,2:ncol(G1)]
G1 <- as.data.frame(t(G1))  #t(G1) returns a matrix, which has only 1 data type
#as.data.frame allows to include multiple data types
G1 <- unfactor(G1)
G1[,3:ncol(G1)] <- lapply(G1[,3:ncol(G1)],as.numeric)#change data type from factor to numeric
typeof(G1[2,3]) #to check that elements are actually numeric
G1[3] <- G1[3] / 3600 # change time from seconds to hours
G1[,4:ncol(G1)] <- G1[,4:ncol(G1)] - 0.091 #correct for background OD
names(G1)[3] <- 'Time_hr'

rownames(G2) <- G2[,1]
G2 <- G2[,2:ncol(G2)]
G2 <- as.data.frame(t(G2))  
G2 <- unfactor(G2) 
G2[,3:ncol(G2)] <- lapply(G2[,3:ncol(G2)],as.numeric) 
typeof(G2[2,3]) 
G2[3] <- G2[3] / 3600
G2[,4:ncol(G2)] <- G2[,4:ncol(G2)] - 0.091
names(G2)[3] <- 'Time_hr'

rownames(G3) <- G3[,1]
G3 <- G3[,2:ncol(G3)]
G3 <- as.data.frame(t(G3))  
G3 <- unfactor(G3) 
G3[,3:ncol(G3)] <- lapply(G3[,3:ncol(G3)],as.numeric) 
typeof(G3[2,3]) 
G3[3] <- G3[3] / 3600
G3[,4:ncol(G3)] <- G3[,4:ncol(G3)] - 0.091 
names(G3)[3] <- 'Time_hr'

rownames(G5) <- G5[,1]
G5 <- G5[,2:ncol(G5)]
G5 <- as.data.frame(t(G5)) 
G5 <- unfactor(G5) 
G5[,3:ncol(G5)] <- lapply(G5[,3:ncol(G5)],as.numeric)
typeof(G5[2,3]) 
G5[3] <- G5[3] / 3600
G5[,4:ncol(G5)] <- G5[,4:ncol(G5)] - 0.091 
names(G5)[3] <- 'Time_hr'

rownames(G6) <- G6[,1]
G6 <- G6[,2:ncol(G6)]
G6 <- as.data.frame(t(G6))  
G6 <- unfactor(G6) 
G6[,3:ncol(G6)] <- lapply(G6[,3:ncol(G6)],as.numeric) 
typeof(G6[2,3]) 
G6[3] <- G6[3] / 3600
G6[,4:ncol(G6)] <- G6[,4:ncol(G6)] - 0.091 
names(G6)[3] <- 'Time_hr'

rownames(G7) <- G7[,1]
G7 <- G7[,2:ncol(G7)]
G7 <- as.data.frame(t(G7))  
G7 <- unfactor(G7) 
G7[,3:ncol(G7)] <- lapply(G7[,3:ncol(G7)],as.numeric) 
typeof(G7[2,3]) 
G7[3] <- G7[3] / 3600
G7[,4:ncol(G7)] <- G7[,4:ncol(G7)] - 0.091 
names(G7)[3] <- 'Time_hr'


# pepare data frame G4 format to store the relevant data (columns with 0.0 ug/mL of antibiotic)
a <- c(G1[1,1],G2[1,1], G3[1,1], G5[1,1],G6[1,1], G7[1,1])
b <- c(G1[1,2],G2[1,2], G3[1,2], G5[1,2],G6[1,2], G7[1,2] )
c <- c(G1[1,3],G2[1,3], G3[1,3], G5[1,3],G6[1,3], G7[1,3])
d <- c(colnames(G1[1:6]))
e <- c(G1[1,4],G2[1,4], G3[1,4], G5[1,4], G6[1,4], G7[1,4])
G4 <- data.frame(a,b,c,d,e, stringsAsFactors = FALSE)
names(G4)[1:3] <- colnames(G1[1:3])
names(G4)[4] <- 'Well'
names(G4)[5] <- 'OD'
G4 <- unfactor(G4)


# data for 1 concentration 
# CHOOSE A CONCENTRATION (here, 0.0ug/mL)
#col_conc <- 1 # for 100 ug/ml (take care with Gent and Kan, divide by 4!. Their maximum concentration was 25 ug/ml)
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
col_conc <- 12 # for 0.0 ug/ml

##  THIS LOOP GATHERS THE DATA OF ONE CONCENTRATION (here 0.0 ug/mL), FOR ALL ANTIBIOTICS (all rows of the experimental plate). 
for (i in 1:8){
  for (j in 1:nrow(G1)){
    G4[(i-1)*nrow(G1)+j,1:3] <- G1[j,1:3]
    G4[(i-1)*nrow(G1)+j,4] <- colnames(G1[12*i+3 - (12 - col_conc)])  #12*i+3 locates the no antibiotic condition
    # -(12 - col_conc) locates the specific concentration that we want, corresponding to the column in the plate
    G4[(i-1)*nrow(G1)+j,5] <- G1[j,12*i+3 -(12 - col_conc)]
    G4[96*8+ (i-1)*nrow(G2)+j,1:3] <- (G2[j,1:3])
    G4[96*8+(i-1)*nrow(G2)+j,4] <- colnames(G2[12*i+3- (12 - col_conc)])
    G4[96*8+(i-1)*nrow(G2)+j,5] <- G2[j,12*i+3- (12 - col_conc)]
    G4[96*16+(i-1)*nrow(G3)+j,1:3] <- G3[j,1:3]
    G4[96*16+(i-1)*nrow(G3)+j,4] <- colnames(G3[12*i+3- (12 - col_conc)])
    G4[96*16+(i-1)*nrow(G3)+j,5] <- G3[j,12*i+3- (12 - col_conc)]
    G4[96*24+(i-1)*nrow(G1)+j,1:3] <- G5[j,1:3]
    G4[96*24+(i-1)*nrow(G1)+j,4] <- colnames(G5[12*i+3- (12 - col_conc)])
    G4[96*24+(i-1)*nrow(G1)+j,5] <- G5[j,12*i+3- (12 - col_conc)]
    G4[96*32+ (i-1)*nrow(G2)+j,1:3] <- (G6[j,1:3])
    G4[96*32+(i-1)*nrow(G2)+j,4] <- colnames(G6[12*i+3- (12 - col_conc)])
    G4[96*32+(i-1)*nrow(G2)+j,5] <- G6[j,12*i+3- (12 - col_conc)]
    G4[96*40+(i-1)*nrow(G3)+j,1:3] <- G7[j,1:3]
    G4[96*40+(i-1)*nrow(G3)+j,4] <- colnames(G7[12*i+3- (12 - col_conc)])
    G4[96*40+(i-1)*nrow(G3)+j,5] <- G7[j,12*i+3- (12 - col_conc)]
  }  
}


write.csv(G4,"NoAntibiotic_GrowthCurves.csv", row.names = TRUE)

######


G8all <- read.table ("NoAntibiotic_GrowthCurves.csv", header=T, sep = ",", dec = ".")
G8all <- G8all[, 2:(ncol(G8all))]
G8all[G8all$OD<0, 5] = 0.00001 #remove negative values after background OD correction

#plot aesthetics
background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_rect (colour="black", size=1.0), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12), legend.position = 'none') 
titulos <- labs(title="", x=expression(paste("Time (hr)")), y=expression(paste("OD")))
escalaY = scale_y_continuous(limits= c(0.001 , 1), trans='log10',  labels = trans_format("log10", math_format(10^.x)))
escalaX = scale_x_continuous(limits = c(0, 20))
lineas <- geom_line()

## prepare data frame for the Ca species
reportCa <- c()
count1 <- 0
outCa <- data.frame(c(100.1),c('D12'),c('R1'), c(0.0),c(0.0),c(0.0), stringsAsFactors=FALSE)
colnames(outCa) = c('Concentration','Well','Replicate','Slope', 'Intercept', 'Rsquare')

col_conc <- 0.0 # Antibiotic concentration in experimental plate column


TechRep <- c('A12','B12','C12','D12','E12','F12','G12','H12') # to choose technical replicate
BioRep <- c('R1','R2','R3') #to choose biological replicate

## iterations start here
#growth rates for the Ca species
Especie <- c('Ca') 
for (i in 1:8){
  for (j in 1:3){
    G8 <- G8all[G8all$Well == TechRep[i] & G8all$Replicate == BioRep[j] & G8all$Species == Especie, ]
    count1 <- count1 + 1
    outCa[count1,] = c(col_conc, paste(G8[1,4]), G8[1,1], 0,0,0)
    
    ploti <- ggplot(data=G8, aes(x= Time_hr, y = OD))
    
    growth<-lm(log(G8[G8$OD>0.01 & G8$OD<0.05,5]) ~ G8[G8$OD>0.01 & G8$OD<0.05,3], data=G8)
    summary(growth)
    
    reportCa <- summary(growth)
    outCa[count1,4] <-reportCa$coefficients[2, 'Estimate']
    outCa[count1,5] <-reportCa$coefficients[1, 'Estimate']
    outCa[count1,6] <-reportCa$r.squared
    
    alphaCa <-  growth$coefficients[[2]] #/ log(10)
    betaCa <- growth$coefficients[[1]] #/ log(10)
    fitlineCa <- geom_segment(aes(x = 1.5, xend = 5.5, y = exp(betaCa + alphaCa*1.5 +0.3 ), yend = exp(betaCa + alphaCa*5.5 +0.3)))
    maltusian <- paste('mu[I]== ', alphaCa)
    notes <- annotate("text", size= 5, x = 16.5, y = 0.01, parse=T, label = as.character(maltusian))
    erre <- expression(R^2 ==0.9771)
    notes2 <- annotate("text", size= 5, x = 16.5, y = 0.001, parse=T, label = as.character(erre))
    
    refer1 <- geom_hline(yintercept = 0.05) 
    refer2 <- geom_hline(yintercept = 0.01)    #refer indicates the OD range considered in the fitting 
    
    plotii = ploti + background  + titulos + lineas + escalaY + escalaX +fitlineCa +  notes + notes2 + refer1 + refer2
    plotii
    
    pdf(paste(  Especie,"-GrowthRate-NoAntib-",G8[1,1],"-Well",G8[1,4], ".pdf"), width=4, height=3)
    print(plotii)
    dev.off()
    
  }
}

#Export Ca growth rates
write.csv(outCa,"GrowthRates-Ca-NoAntibiotic-ODRange-0.01to0.05.csv", row.names = TRUE)


## growth rates for the Lp species
reportLp <- c()
count1 <- 0
outLp <- data.frame(c(100.1),c('D12'),c('R1'), c(0.0),c(0.0),c(0.0), stringsAsFactors=FALSE)
colnames(outLp) = c('Concentration','Well','Replicate','Slope', 'Intercept', 'Rsquare')

col_conc <- 0.0 # Antibiotic concentration in experimental plate column

TechRep <- c('A12','B12','C12','D12','E12','F12','G12','H12') # to choose technical replicate
BioRep <- c('R1','R2','R3') #to choose biological replicate

Especie <- c('Lp') 
## iterations start here
for (i in 1:8){
  for (j in 1:3){
    G8 <- G8all[G8all$Well == TechRep[i] & G8all$Replicate == BioRep[j] & G8all$Species == Especie, ]
    count1 <- count1 + 1
    outLp[count1,] = c(col_conc, paste(G8[1,4]), G8[1,1], 0,0,0)
    
    ploti <- ggplot(data=G8, aes(x= Time_hr, y = OD))
    
    # tuned for Lp
    growth<-lm(log(G8[G8$OD>0.01 & G8$OD<0.05,5]) ~ G8[G8$OD>0.01 & G8$OD<0.05,3], data=G8)
    summary(growth)
    
    reportLp <- summary(growth)
    outLp[count1,4] <-reportLp$coefficients[2, 'Estimate']
    outLp[count1,5] <-reportLp$coefficients[1, 'Estimate']
    outLp[count1,6] <-reportLp$r.squared
    
    alphaLp <-  growth$coefficients[[2]] #/ log(10)
    betaLp <- growth$coefficients[[1]] #/ log(10)
    fitlineLp <- geom_segment(aes(x = 1.5, xend = 5.5, y = exp(betaLp + alphaLp*1.5 +0.3 ), yend = exp(betaLp + alphaLp*5.5 +0.3)))
    maltusian <- paste('mu[I]== ', alphaLp)
    notes <- annotate("text", size= 5, x = 16.5, y = 0.01, parse=T, label = as.character(maltusian))
    erre <- expression(R^2 ==0.9771)
    notes2 <- annotate("text", size= 5, x = 16.5, y = 0.001, parse=T, label = as.character(erre))
    
    refer1 <- geom_hline(yintercept = 0.05) 
    refer2 <- geom_hline(yintercept = 0.01)    #refer indicates the OD range considered in the fitting 
    
    plotii = ploti + background  + titulos + lineas + escalaY + escalaX + fitlineLp + notes + notes2 + refer1 + refer2
    plotii
    
    pdf(paste( Especie,"-GrowthRate-NoAntib-",G8[1,1],"-Well",G8[1,4],".pdf"), width=4, height=3)
    print(plotii)
    dev.off()
    
  }
}

#Export Lp Growth rates
write.csv(outLp,"GrowthRates-Lp-NoAntibiotic-ODRange-0.01to0.05.csv", row.names = TRUE)

#Export Mean Growth Rates
MeanGrowthRates <- data.frame(0, 0,0,0, stringsAsFactors = FALSE)
names(MeanGrowthRates)[1:4] <- c('MeanGrowthRateCa', 'MeanGrowthRateLp', 'MeanInterceptCa','MeanInterceptLp')
MeanGrowthRates[1,1] <- mean(sapply(outCa[4], as.numeric))
MeanGrowthRates[1,2] <- mean(sapply(outLp[4], as.numeric))
MeanGrowthRates[1,3] <- mean(sapply(outCa[5], as.numeric))
MeanGrowthRates[1,4] <- mean(sapply(outLp[5], as.numeric))

write.csv(MeanGrowthRates,"Mean_Growth_Rates_and_Intercepts.csv", row.names = FALSE)
print(MeanGrowthRates)
