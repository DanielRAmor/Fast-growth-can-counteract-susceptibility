#This code reads the data from 6 experimental plates (3 for Ca, 3 for Lp)
#The raw data contains the Optical Density (OD) of a 24hr incubation in different antibiotics and antibiotic concentrations, with measurements every 15mins
#The algorithm organizes the data of the 3 replicates per species into one output file 

library (plyr)  
library(scales)
library(varhandle)

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
G1 <- unfactor(G1) #avoids treating any variable under the data type 'factor'
G1[,3:ncol(G1)] <- lapply(G1[,3:ncol(G1)],as.numeric)#convert data to 'numeric' type
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

# Prepare the format of a data frame (that will contain the data for the specific antibiotic susceptibility that we want to plot)
a <- c(G1[1,1],G2[1,1], G3[1,1], G5[1,1],G6[1,1], G7[1,1])
b <- c(G1[1,2],G2[1,2], G3[1,2], G5[1,2],G6[1,2], G7[1,2] )
c <- c(G1[1,3],G2[1,3], G3[1,3], G5[1,3],G6[1,3], G7[1,3])
d <- c(colnames(G1[1:6]))
e <- c(G1[1,4],G2[1,4], G3[1,4], G5[1,4], G6[1,4], G7[1,4])
G4 <- data.frame(a,b,c,d,e, stringsAsFactors = FALSE)
names(G4)[1:3] <- colnames(G1[1:3])
names(G4)[4] <- 'Well'
names(G4)[5] <- 'OD'


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
#
#choose the antibiotic of interest! (leave uncommented only that one)
# A = Amp
# B = Car
# C = Cip
# D = Chl
# E = Ery
# F = Gen
# G = Kan
# H = Tet

#prepare data frame with convenient format
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

#organize data for all concentrations and antibiotics
col_conc <- 0
for (k in 1:12){
  col_conc <- col_conc + 1
  for (i in 1:8){
    for (j in 1:nrow(G1)){
      G4[(k-1)*96*48 + (i-1)*nrow(G1)+j,1:3] <- G1[j,1:3]
      G4[(k-1)*96*48 + (i-1)*nrow(G1)+j,4] <- colnames(G1[12*i+3 - (12 - col_conc)])  #12*i+3 locates the no antibiotic condition
                                                                                      # -(12 - col_conc) locates the specific concentration that we want, corresponding to the column in the plate
      G4[(k-1)*96*48 + (i-1)*nrow(G1)+j,5] <- G1[j,12*i+3 -(12 - col_conc)]
      G4[(k-1)*96*48 + 96*8+ (i-1)*nrow(G2)+j,1:3] <- (G2[j,1:3])
      G4[(k-1)*96*48 + 96*8+(i-1)*nrow(G2)+j,4] <- colnames(G2[12*i+3- (12 - col_conc)])
      G4[(k-1)*96*48 + 96*8+(i-1)*nrow(G2)+j,5] <- G2[j,12*i+3- (12 - col_conc)]
      G4[(k-1)*96*48 + 96*16+(i-1)*nrow(G3)+j,1:3] <- G3[j,1:3]
      G4[(k-1)*96*48 + 96*16+(i-1)*nrow(G3)+j,4] <- colnames(G3[12*i+3- (12 - col_conc)])
      G4[(k-1)*96*48 + 96*16+(i-1)*nrow(G3)+j,5] <- G3[j,12*i+3- (12 - col_conc)]
      G4[(k-1)*96*48 + 96*24+(i-1)*nrow(G1)+j,1:3] <- G5[j,1:3]
      G4[(k-1)*96*48 + 96*24+(i-1)*nrow(G1)+j,4] <- colnames(G5[12*i+3- (12 - col_conc)])
      G4[(k-1)*96*48 + 96*24+(i-1)*nrow(G1)+j,5] <- G5[j,12*i+3- (12 - col_conc)]
      G4[(k-1)*96*48 + 96*32+ (i-1)*nrow(G2)+j,1:3] <- (G6[j,1:3])
      G4[(k-1)*96*48 + 96*32+(i-1)*nrow(G2)+j,4] <- colnames(G6[12*i+3- (12 - col_conc)])
      G4[(k-1)*96*48 + 96*32+(i-1)*nrow(G2)+j,5] <- G6[j,12*i+3- (12 - col_conc)]
      G4[(k-1)*96*48 + 96*40+(i-1)*nrow(G3)+j,1:3] <- G7[j,1:3]
      G4[(k-1)*96*48 + 96*40+(i-1)*nrow(G3)+j,4] <- colnames(G7[12*i+3- (12 - col_conc)])
      G4[(k-1)*96*48 + 96*40+(i-1)*nrow(G3)+j,5] <- G7[j,12*i+3- (12 - col_conc)]
    }  
  }
}

#export data
write.csv(G4,"Growth_Curves_in_Antibiotics.csv", row.names = TRUE)

