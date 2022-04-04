library(ggplot2)
library (dplyr)  
library(scales)
library(varhandle)
library(reshape)
library(stringr)
library(plotrix)
#Load the raw data (Overnight culture in the plate reader with reads every 15mins)
G1 <- read.table ("Growth_Vs_Initial_Abundance.dat", header=F, sep = "\t", dec = ".")
#data processing
# use the first column as row labels and transpose the table
rownames(G1) <- G1[,1]  
G1 <- G1[,2:ncol(G1)]
G1 <- as.data.frame(t(G1))  
typeof(G1[2,3]) #to check that elements are actually numeric
G1[1] <- G1[1] / 3600 # Convert time units to hours
names(G1)[1] <- 'Time_hr'
G1[,4:ncol(G1)] <- G1[,4:ncol(G1)] - 0.08894 #correct for background Optical Density
#in case there is any negative value resulting from background correction
for (i in 1:nrow(G1)) {
  for (j in 1:ncol(G1)) {
      if (G1[i,j]<0) {G1[i,j]=0.0001}
  }
}

# LEGEND
# Ca was inoculated in experimental plate rows B, C, D
# Lp was inoculated in experimental plate rows E, F, G

# 1.2*10^6 cells inoculated in experimental plate column 11
# 3.6*10^5 cells inoculated in experimental plate column 10
# 1.2*10^5 cells inoculated in experimental plate column 9
# 3.6*10^4 cells inoculated in experimental plate column 8
# 1.2*10^4 cells inoculated in experimental plate column 7
# 3.6*10^3 cells inoculated in experimental plate column 6
# 1.2*10^3 cells inoculated in experimental plate column 5

#additional data formatting, for plotting purposes
G2 <- as.data.frame( melt(G1, id.vars=1, value.name="OD", variable.name="Well"))  ## data frame that will all OD values in one column only
G2[,4] <- str_sub(G2[,2], 1,1)  #plate row
G2[,5] <- str_sub(G2[,2], 2)  #plate column
G2[,5] <- sapply(G2[5], as.numeric)
G2[6] <- c('None')
# select columns with conditions of interest
G2 <- G2[G2[,5]>4 & G2[,5]!=12 ,]
#assign species to data
G2[G2[4]=='B'| G2[4]=='C' | G2[4]=='D',6] <- c('Ca')
G2[G2[4]=='F'| G2[4]=='G' | G2[4]=='H',6] <- c('Lp')
#select rows with cells
G2 <- G2[G2[,6]!= 'None' ,]
# names of the columns of the formatted data
colnames(G2) = c('Time_hr', 'Well', 'OD', 'row', 'column', 'species')
well <- distinct(G2[2])
well[] <- lapply(well, as.character)

#prepare a data frame to compute the fold growth
FoldG <- distinct(G2[2])
FoldG[1] <- lapply(FoldG[1], as.character)
FoldG[,2] <- str_sub(FoldG[,1], 1,1)  #plate row
FoldG[,3] <- str_sub(FoldG[,1], 2)  #plate column
FoldG[4] <-0
FoldG[5] <-0
FoldG[6] <-c('s')
FoldG[7] <-0
FoldG[8] <-0
names(FoldG) <- c('well', 'row', 'column', 'Tthreshold', 'ODthreshold',  'Species','Abundance0','EffGRate')

#Find time spend to reach the threshold OD!
ODThreshold <- 0.01
for (i in 1:nrow(FoldG)) {
  aa = FoldG[i,1]
  datos = G2[G2$Well ==aa,]
  Trow = which(abs(datos[3]-ODThreshold)==min(abs(datos[3]-ODThreshold))) #finds the data frame's row in which the abundance is closest to the threshold
  Trow = min(Trow)  #in case multiple OD values are equally distant to Threshold
  if (min(datos[3]-ODThreshold)<0) {
    if (max(datos[3]-ODThreshold)<0) {Trow = nrow(datos)}   #if it never crossed the threshold, we assign latest time in the data!
  } 
  FoldG[i,4:6] = c(datos[Trow,1],datos[Trow,3],datos[Trow,6]) #save the time it reached for the culture to reach the OD threshold
}

#inital OD is measured on the experimental plate (column 11) for an initial abundance of 10^6
ODiCa <- mean( c(G2[G2$Well == 'B11' & G2$Time_hr == 0.0 ,3],G2[G2$Well == 'C11' & G2$Time_hr == 0.0 ,3],G2[G2$Well == 'D11' & G2$Time_hr == 0.0 ,3]))
ODiLp <- mean( c(G2[G2$Well == 'F11' & G2$Time_hr == 0.0 ,3],G2[G2$Well == 'G11' & G2$Time_hr == 0.0 ,3],G2[G2$Well == 'H11' & G2$Time_hr == 0.0 ,3]))

#compute the effective Growth Rate through log(OD_NearestToThreshold / Initial_OD) / (Time_to_reach_OD_threshold), 
#where Initial_OD is the OD corresponding to and abundance of 10^6 CFU/mL divided by the corresponding dilution factor for each condition
ODi = ODiCa
FoldG[4] <- lapply(FoldG[4], as.numeric)
FoldG[5] <- lapply(FoldG[5], as.numeric)
for (i in 1:(nrow(FoldG))){
  if (FoldG[i,6] == 'Ca') {
    ODi = ODiCa
  }else{
    ODi = ODiLp
  }
  if (FoldG[i,3]==10) {FoldG[i,7:8] = c(3.6*10^5,log(FoldG[i,5] /(0.3*ODi))/FoldG[i,4])}
  if (FoldG[i,3]==9) {FoldG[i,7:8] = c(1.2*10^5,log(FoldG[i,5] /(0.1*ODi))/FoldG[i,4])}
  if (FoldG[i,3]==8) {FoldG[i,7:8] = c(3.6*10^4,log(FoldG[i,5] /(0.03*ODi))/FoldG[i,4])}
  if (FoldG[i,3]==7) {FoldG[i,7:8] = c(1.2*10^4,log(FoldG[i,5] /(0.01*ODi))/FoldG[i,4])}
  if (FoldG[i,3]==6) {FoldG[i,7:8] = c(3.6*10^3,log(FoldG[i,5] /(0.003*ODi))/FoldG[i,4])}
  if (FoldG[i,3]==5) {FoldG[i,7:8] = c(1.2*10^3,log(FoldG[i,5] /(0.001*ODi))/FoldG[i,4])}
}

#the condition with 10^6 initial cells is out of the region of interest in this analysis.
#we used this condition to measure the initial OD, which in this case it's too close to the threshold OD for this method to be accurate
FoldG <- FoldG[FoldG[,3]!=11,] 


#convert to numeric
FoldG[8] <- lapply(FoldG[8], as.numeric)

#compute mean and St Error
for (i in 1:nrow(FoldG)){
  FoldG[i,9] <- mean(FoldG[FoldG$Species==FoldG[i,6]&FoldG$Abundance0==FoldG[i,7],8]) 
  FoldG[i,10] <- std.error(FoldG[FoldG$Species==FoldG[i,6]&FoldG$Abundance0==FoldG[i,7],8]) 
}
names(FoldG)[9:10] <- c('MeanEffGRate','StErr')


#select the initial population densities of interest to make the plot
ploti <- ggplot(FoldG, aes(x=Abundance0, y= MeanEffGRate, fill= Species, group=Species, color=Species))


#plot aesthetics
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) 
err <- geom_errorbar(aes(ymin = MeanEffGRate - StErr, ymax =MeanEffGRate + StErr), size=0.8, width = 0.15)
puntos <- geom_point()
lineas <- geom_line()
escalaXlog <- scale_x_continuous(trans='log10', breaks = c(1.2*10^3,3.6*10^3,1.2*10^4,3.6*10^4,1.2*10^5,3.6*10^5))
escalaYLin <- scale_y_continuous(expand = c(0,0), limits=c(0.2,0.55))
couleur <- scale_colour_manual(values = c("#1AB24B","#FAA421"))
titulo <- ggtitle("")
xaxe <- xlab ("Initial abundance (CFU/ml)")
yaxe <- ylab("Effective per capita growth rate (1/hr)") 

plotii <- ploti + background + lineas + err + escalaXlog + escalaYLin + couleur + xaxe + yaxe

#plot it
plotii

#export the plot
pdf(paste("F_3B.pdf"), width=4, height=2.5)
print(plotii)
dev.off()


