library (ggplot2)
library (plyr)  #for reordering the axis values
library(scales)
library(varhandle)

M <- read.csv ("GrowthRateVsAntibioticConcentration-3Replicates.csv", header=T, dec = ".")   #read data
M <- unfactor(M)

#prepare a data frame with convenient format
a <- c(M[1,5], M[13,5])
b <- c(M[1,1], M[2,1])
c <- c(M[1,7], M[2,6])
d <- c(M[1,16], M[2,16])
S <- data.frame(a,b,c,d,c,d,c,d,c,d, stringsAsFactors = FALSE)
names(S)[1:10] <- c('Antibiotic','Species', 'R1_IC50','R1_RelativeGrowthRateIC50', 'R2_IC50','R2_RelativeGrowthRateIC50', 'R3_IC50','R3_RelativeGrowthRateIC50', 'AverageIC50', 'StErrIC50')

#'for' loop to look for the minimum concentration that inhibits the growth rate falls by more than half the maximum  
count = 0
for (aa in 1:16) {  #because we have 8 different antibiotics and 2 species
  count = 0
    #Replicate 1
  S[aa,4] <- 1.0
  for (cc in 1:12) {  # because we have 12 different concentrations
    if (count == 0){
      # in decreasing order of concentrations, we look for the first concentration at which the growth rate exceeds half the maximum
      if ( M[cc + (aa-1)*12, 16] > 0.5) {   
        count = 1
        S[aa,1] <- M[cc+ (aa-1)*12, 5]  #Antibiotics
        S[aa,2] <- M[cc+ (aa-1)*12, 1]   #Species
        S[aa,3] <- M[cc+ (aa-1)*12 -1, 6]  #This is the IC50. The -1 is because we found 
                                            #the minimum concentration for which the growth rate exceeds 1/2 maximum 
                                           #so the previous one is the IC50
        
        S[aa,4] <- M[cc+ (aa-1)*12 -1, 16]  #Relative growth at the measured IC50
      }
    }
  }
  
  #Replicate 2
  count = 0
  S[aa,6] <- 1.0
  for (cc in 1:12) {  
    if (count==0){
      if ( M[cc + (aa-1)*12, 17] > 0.5) {  
        count =1
        S[aa,5] <- M[cc+ (aa-1)*12 -1, 6]
        S[aa,6] <- M[cc+ (aa-1)*12 -1, 17]
      }
    }
  }
  
  #Replicate 3
  count = 0
  S[aa,8] <- 1.0
  for (cc in 1:12) {  
    if (count==0){
      if (M[cc + (aa-1)*12, 18] > 0.5) {  
        count = 1
        S[aa,7] <- M[cc+ (aa-1)*12 -1, 6]
        S[aa,8] <- M[cc+ (aa-1)*12 -1, 18]
      }
    }
  }
}


S[9] <- (S[3] + S[5] + S[7]) /3 #Average IC50
for (i in 1:nrow(S)) {
S[i,10] <- sd(c(S[i,3],S[i,5], S[i,7])) / sqrt(3) #Standard Error IC50
}

write.csv (S, "IC50-BasedOnGrowthRates.csv")   #export data

#Plot the IC50 of each species for each antibiotic
ploti <- ggplot(data = S, aes(x= Species, y = AverageIC50, color= Species, fill = Species))
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) 
titulos <- labs(title="", x=expression(paste("Species")), y=expression(paste("IC50 (Conc. that halves the max growth rate" )))
barras <- geom_bar(stat="identity")
err1 <- geom_errorbar(aes(ymin = AverageIC50 - StErrIC50, ymax =AverageIC50 + StErrIC50), size=0.75, color = 'gray2', width = 0.2) 
paneles <- facet_wrap(~ Antibiotic, scales = 'free_y')
colores <- scale_fill_manual(values=c('#1AB24B', '#FAA421')) 
colores2 <- scale_color_manual(values=c('#1AB24B', '#FAA421')) 
plotii = ploti + background  + titulos + barras +  paneles + colores + colores2 + err1
plotii

#export the plot
pdf(paste(  "S_3A_IC50_BasedOnGrowthRates", ".pdf"), width=4, height=4)
print(plotii)
dev.off()
