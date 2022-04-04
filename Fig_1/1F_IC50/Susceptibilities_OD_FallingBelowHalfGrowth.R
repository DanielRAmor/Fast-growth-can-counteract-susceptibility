#this code reads the data on growth inhibition (via growth after 24hr) at different antibiotic concentrations
#and finds the IC50 for each of the 8 antibiotics

library (ggplot2)
library (plyr)  #for reordering the axis values
library(scales)
library(varhandle)

#read the data
M <- read.table ("Inhibition_vs_Concentration.dat", header=T, sep = "\t", dec = ".")   #read data
M <- unfactor(M)

#prepare format for data.frame
a <- c(M[1,16], M[13,16])
b <- c(M[1,3], M[2,3])
c <- c(M[1,7], M[2,7])
d <- c(M[1,7], M[2,7])
S <- data.frame(a,b,c,d,c,d,c,d,c,d, stringsAsFactors = FALSE)
names(S)[1:10] <- c('Antibiotic','Species', 'R1_MI50','R1_RelativeGrowthMI50', 'R2_MI50','R2_RelativeGrowthMI50', 'R3_MI50','R3_RelativeGrowthMI50', 'AverageIC50', 'StErrIC50')

#algorithm for finding the concentration at which growth falls below half maximal
  #Experimental Replicate 1
  count = 0
  for (aa in 1:16) {  #because we have 8 different antibiotics and 2 species
    count = 0
      #Replicate 1
    S[aa,4] <- 1.0
    for (cc in 1:12) {  # because we have 12 different concentrations
      if (count == 0){
        # to find the IC50 as the minimum concentration that inhibits growth by more than half the maximum growth:
        if ( M[cc + (aa-1)*12, 9] > 0.5) {   
          count = 1
          S[aa,1] <- M[cc+ (aa-1)*12, 16]  #Antibiotics
          S[aa,2] <- M[cc+ (aa-1)*12, 3]   #Species
          S[aa,3] <- M[cc+ (aa-1)*12 -1, 2]  #Concentration. the -1 is because we have found the first concentration that grows over 0.5. The previous concentration is the IC50.
                                              #But I am estimating IC50 as the previous concentration (first one that inhibits more of 50% growth)
          S[aa,4] <- M[cc+ (aa-1)*12 -1, 9]  #Relative growth 
        }
      }
    }
    
    #Same for experimental Replicate 2
    count = 0
    S[aa,6] <- 1.0
    for (cc in 1:12) {  
      if (count==0){
        if ( M[cc + (aa-1)*12, 12] > 0.5) {   
          count =1
          S[aa,5] <- M[cc+ (aa-1)*12 -1, 2]
          S[aa,6] <- M[cc+ (aa-1)*12 -1, 12]
        }
      }
    }
    
    #Same for experimental Replicate 3
    count = 0
    S[aa,8] <- 1.0
    for (cc in 1:12) {  
      if (count==0){
        if (M[cc + (aa-1)*12, 15] > 0.5) {   
          count = 1
          S[aa,7] <- M[cc+ (aa-1)*12 -1, 2]
          S[aa,8] <- M[cc+ (aa-1)*12 -1, 15]
        }
      }
    }
  }
  
  
  S[9] <- (S[3] + S[5] + S[7]) /3 #Average IC50
  for (i in 1:nrow(S)) {
  S[i,10] <- sd(c(S[i,3],S[i,5], S[i,7])) / sqrt(3) #Standard Error for the IC50
  }

write.csv (S, "F_1F_data.csv")   #export data


