library (ggplot2)
library (plyr)  
library(scales)
library(varhandle)

#README:
#This code reads the data from 6 experimental plates (3 for Ca, 3 for Lp)
#The raw data contains the Optical Density (OD) of a 24hr incubation, with measurements every 15mins
#The algorithm analyzes and prepares the data for measuring and plotting the susceptibility to the antibiotic of interest
#The susceptibility is determined via the normalized OD (relative to the OD in the absence of antibiotic) at 24hr

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

#
#choose the antibiotic of interest! (leave uncommented only that one)
#a=0 # Amp
#a=12 # Car
#a=24 # Cip
#a=36 # Chl
#a=48 # Ery
#a=60 # Gen
a=72 # Kan
#a=84 # Tet

# the 'for' loop below transfers the experimental data for the antibiotic of interest to the dataframe G4
#and it also processes the data to obtain the Normalized OD
for (i in 1:12){
  for (j in 1:nrow(G1)){
    G4[(i-1)*nrow(G1)+j,1:3] <- G1[j,1:3]
    G4[(i-1)*nrow(G1)+j,4] <- colnames(G1[i+a+3])
    G4[(i-1)*nrow(G1)+j,5] <- G1[j,i+a+3]
    G4[96*12+ (i-1)*nrow(G2)+j,1:3] <- (G2[j,1:3])
    G4[96*12+(i-1)*nrow(G2)+j,4] <- colnames(G2[i+a+3])
    G4[96*12+(i-1)*nrow(G2)+j,5] <- G2[j,i+a+3]
    G4[96*24+(i-1)*nrow(G3)+j,1:3] <- G3[j,1:3]
    G4[96*24+(i-1)*nrow(G3)+j,4] <- colnames(G3[i+a+3])
    G4[96*24+(i-1)*nrow(G3)+j,5] <- G3[j,i+a+3]
    G4[96*36+(i-1)*nrow(G1)+j,1:3] <- G5[j,1:3]
    G4[96*36+(i-1)*nrow(G1)+j,4] <- colnames(G5[i+a+3])
    G4[96*36+(i-1)*nrow(G1)+j,5] <- G5[j,i+a+3]
    G4[96*48+ (i-1)*nrow(G2)+j,1:3] <- (G6[j,1:3])
    G4[96*48+(i-1)*nrow(G2)+j,4] <- colnames(G6[i+a+3])
    G4[96*48+(i-1)*nrow(G2)+j,5] <- G6[j,i+a+3]
    G4[96*60+(i-1)*nrow(G3)+j,1:3] <- G7[j,1:3]
    G4[96*60+(i-1)*nrow(G3)+j,4] <- colnames(G7[i+a+3])
    G4[96*60+(i-1)*nrow(G3)+j,5] <- G7[j,i+a+3]
  }  
}

#right after importing it above, the data was corrected for the background OD
#below we prevent any potential case of negative OD value after background correction
G4[G4$OD < 0, 5] = 0.001  

#just renaming some indexes for proper x-axis ordering when plotting
G4$Well[G4$Well=='A10'] <- 'AB10'  
G4$Well[G4$Well=='A11'] <- 'AB11'
G4$Well[G4$Well=='A12'] <- 'AB12'
G4$Well[G4$Well=='B10'] <- 'BB10'  
G4$Well[G4$Well=='B11'] <- 'BB11'
G4$Well[G4$Well=='B12'] <- 'BB12'
G4$Well[G4$Well=='C10'] <- 'CB10' 
G4$Well[G4$Well=='C11'] <- 'CB11'
G4$Well[G4$Well=='C12'] <- 'CB12'
G4$Well[G4$Well=='D10'] <- 'DB10'  
G4$Well[G4$Well=='D11'] <- 'DB11'
G4$Well[G4$Well=='D12'] <- 'DB12'
G4$Well[G4$Well=='E10'] <- 'EB10'  
G4$Well[G4$Well=='E11'] <- 'EB11'
G4$Well[G4$Well=='E12'] <- 'EB12'
G4$Well[G4$Well=='F10'] <- 'FB10'  
G4$Well[G4$Well=='F11'] <- 'FB11'
G4$Well[G4$Well=='F12'] <- 'FB12'
G4$Well[G4$Well=='G10'] <- 'GB10'  
G4$Well[G4$Well=='G11'] <- 'GB11'
G4$Well[G4$Well=='G12'] <- 'GB12'
G4$Well[G4$Well=='H10'] <- 'HB10'  
G4$Well[G4$Well=='H11'] <- 'HB11'
G4$Well[G4$Well=='H12'] <- 'HB12'

# G4 contains the OD time series under concentrations of the chosen antibiotic
# choose time point to plot OD and store the data in data frame G9
G9 <- G4
G9 = G9[G9$Time_hr > 23.7,]

#Make time point uniform, to correctly color the curves when plotting
G9[,3] <- 24.0

#plot aesthetics
background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_rect (colour="black", size=1.0), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) 
titulos <- labs(title="", x=expression(paste("Concentration (ug/mL)")), y=expression(paste("Normalized OD")))
escalaY = scale_y_continuous(trans='log10',  labels = trans_format("log10", math_format(10^.x)), limits = c(0.001, 1.3))
lineas <- geom_line()
leyenda <- theme(legend.position = "none")
colores <- scale_colour_manual(values = c("#1AB24B", '#FAA51A',"#168936", "#24CE54" , '#D88B1A', '#BA7D25')) 

#Normalization, CODE CHANGES WITH CONDITION! If editing this code, substitute below by the appropriate label 'XB12', where 'X' is the corresponding row in the experimental plate
maxCa24 <- mean (G9[G9$Species == 'Ca' & G9$Time_hr > 23 & G9$Well == 'GB12', 5])
maxLp24 <- mean (G9[G9$Species == 'Lp' & G9$Time_hr > 23 & G9$Well == 'GB12', 5])

G9[G9$Species == 'Ca' & G9$Time_hr > 23 , 6] <- G9[G9$Species == 'Ca' & G9$Time_hr > 23 , 5] / maxCa24
G9[G9$Species == 'Lp' & G9$Time_hr > 23 , 6] <- G9[G9$Species == 'Lp' & G9$Time_hr > 23 , 5] / maxLp24

names(G9)[6] <- 'Relative_OD'

#FOR ALL ANTIBIOTICS EXCEPT KAN & GEN, CONCENTRATIONS ARE:
#xlabels <- scale_x_discrete(limits = unique(rev(G9$Well)),labels=c("EB12" = "0.0", "EB11" = ".10","EB10" = "0.20", "E9" = ".39", "E8" = ".78", 
#                                                                   "E7" = "1.6","E6" = "3.1","E5" = "6.3","E4" = "12.5","E3" = "25","E2" = "50","E1" = "100"))

#FOR THE CASES OF KANAMYCIN AND GENTAMICIN, INSTEAD:
#xlabels <- scale_x_discrete(limits = unique(rev(G9$Well)),labels=c("FB12" = "0.0", "FB11" = ".02","FB10" = ".05", "F9" = ".10", "F8" = ".20", 
#                                                                   "F7" = ".39","F6" = ".78","F5" = "1.6","F4" = "3.1","F3" = "6.3","F2" = "12.5","F1" = "25"))
xlabels <- scale_x_discrete(limits = unique(rev(G9$Well)),labels=c("GB12" = "0.0", "GB11" = ".02","GB10" = ".05", "G9" = ".10", "G8" = ".20", 
                                                                   "G7" = ".39","G6" = ".78","G5" = "1.6","G4" = "3.1","G3" = "6.3","G2" = "12.5","G1" = "25"))


IC50line <- geom_hline(yintercept=0.5, linetype="dashed", color = "gray", size=1)
ploti <- ggplot(data = G9, aes(x= Well, y = Relative_OD , group= interaction(Replicate, Time_hr, Species), color= interaction (Time_hr, Species), fill = Replicate))
plotii = ploti + background  + titulos + lineas + colores  + leyenda + xlabels + IC50line
plotii

pdf(paste(  "Growth_in_Kanamycin", ".pdf"), width=6, height=4)
print(plotii)
dev.off()



