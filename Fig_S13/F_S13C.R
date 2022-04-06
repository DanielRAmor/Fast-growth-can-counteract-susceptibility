# Analisys on experimental data to assess cross-protective interactions via culturing the Lp species 
#in Ca's supernatant in the presence of antibiotics.

library(reshape)
library(stringr)
library(dplyr)
library(plotrix)
library(ggplot2)

#OD BACKGROUND OF MEDIA WITH AND WITHOUT ANTIBIOTICS
Day2Background <- read.table ('Day2_NoCells_AntibioticBackground.dat', header=F, sep = "\t", dec = ".")
Day2Background <- Day2Background[,2:13] 
#Columns 1-3 are for media + 16X antibiotic concentration. This one actually comes 50% from fresh media with antibiotics + 50% of media that spent the last 24hr at 30C, as in supernatant experiment 
#Columns 4-6 are for media + 1X antibiotic concentration. 100% fresh media with antibiotics, since I had not enough volume to mix it with media from the previous day controls (which I used to add cells to).
#Columns 7-9 are for a 1/8 dilution of the media in columns 4-6.
#Columns 10-12 are for Base Medium background (no antibiotics, no cells).
colnames(Day2Background) <- c('16Xa', '16Xb','16Xc','1Xa','1Xb','1Xc','1_8Xb','1_8Xb','1_8c','BMa','BMb', 'BMc')
rownames(Day2Background) <- c('Ery', 'Chl', 'Gen', 'Amp', 'Kan', 'Car', 'Tet', 'Cip')

Day2Background[,13] <-rowMeans(Day2Background[,1:3],na.rm=TRUE)  
Day2Background[,14] <-rowMeans(Day2Background[,4:6],na.rm=TRUE)  
Day2Background[,15] <-rowMeans(Day2Background[,7:9],na.rm=TRUE)  
Day2Background[,16] <-rowMeans(Day2Background[,10:12],na.rm=TRUE)  
colnames(Day2Background)[13:16] <- c('Bg16X','Bg1X','Bg1_8X','BgBM')
BMbackground<-mean(Day2Background[,16])

#From the output .csv in Fig. S12, measure the final OD of Ca in the absence of interactions and antibiotics
LpEco <- read.csv('LpEco.csv')
LpK = LpEco[11:13]# Carrying capacity as OD in absence of antibiotics
LpK[7:12,1] = LpK[1:6,2]
LpK[13:18,1] = LpK[1:6,3]
LpK = LpK[,-(2:3)]
LpMeanK = mean(LpK)
LpErrK = sd(LpK)/sqrt(18)


#OD Lp cells in 1X antibiotics, Ca Supernatant
Day2Lp_1X_in_CaS <- read.table ('Day2_Lp_in_1X_CaSupernatant.dat', header=F, sep = "\t", dec = ".")
colnames(Day2Lp_1X_in_CaS) <- c('Antibiotic','30X_1', '30X_2','30X_3','HighD_1','HighD_2','HighD_3','LowD_1','LowD_2','LowD_3','OnlyAnt_1','OnlyAnt_2', 'OnlyAnt_3')
Day2Lp_1X_in_CaS[1] <- c('Ery', 'Chl', 'Gen', 'Amp', 'Kan', 'Car', 'Tet', 'Cip')
#Subtract background OD for each antibiotic
for (i in 2:13){
  Day2Lp_1X_in_CaS[i] <- Day2Lp_1X_in_CaS[i] - Day2Background[14]
}

Lp_1X_in_CaS <-  as.data.frame( melt(Day2Lp_1X_in_CaS,  id.vars=c('Antibiotic'),  value.name="OD")) #formatting
#add replicate indexing
Lp_1X_in_CaS[,4] <- substr(Lp_1X_in_CaS[,2], nchar(as.character(Lp_1X_in_CaS[,2])),nchar(as.character(Lp_1X_in_CaS[,2])))
#remove replicate indexing from column 2
Lp_1X_in_CaS[,2] <- substr(Lp_1X_in_CaS[,2], 1,nchar(as.character(Lp_1X_in_CaS[,2]))-2)
colnames(Lp_1X_in_CaS) <- c('Antibiotic','Supernatant', 'OD', 'Replicate')

#Lp cells in 1X Ca supernatant
#Mean and StError for each Antibiotic&Supernatant condition
AntibClass <- unlist(Lp_1X_in_CaS[1] %>% unique)  # unique is creating a list, but we need a vector, so we also use unlist()
SupernClass <- unlist(dplyr::distinct(Lp_1X_in_CaS[2])) #alternative way of getting unique values
Stats1XLp_in_Ca <- data.frame(Antibiotic = c('Chl','Chl'),Supernatant=c('HighD','LowD'), MeanOD=c(0.1,0.5),StErrOD=c(0.01,0.02))
replicates <- 3  
for (anti in 1:length(AntibClass)){
  for (super in 1:length(SupernClass)){
    Stats1XLp_in_Ca[super+(anti-1)*length(SupernClass),1] <- AntibClass[anti]
    Stats1XLp_in_Ca[super+(anti-1)*length(SupernClass),2] <- SupernClass[super]
    Stats1XLp_in_Ca[super+(anti-1)*length(SupernClass),3] <- mean(Lp_1X_in_CaS[Lp_1X_in_CaS[,1]== AntibClass[anti] & Lp_1X_in_CaS[,2]== SupernClass[super],3])
    Stats1XLp_in_Ca[super+(anti-1)*length(SupernClass),4] <-   sd(Lp_1X_in_CaS[Lp_1X_in_CaS[,1]== AntibClass[anti] & Lp_1X_in_CaS[,2]== SupernClass[super],3])/sqrt(replicates)
  }
} 


#Plot
Plot1XLp <- ggplot(subset(Stats1XLp_in_Ca), aes(x= Supernatant, y= MeanOD/LpMeanK, group= Supernatant))
background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_rect (colour="black", size=1.5), axis.text=element_text(size=10), axis.title=element_text(size=12), title=element_text(size=22))
titulos<- labs(title="", x=expression(paste("Supernatant Type")), y=expression(paste("24hr abundance relative to maximum")))
escalaY <- scale_y_continuous(breaks = pretty(Stats1XLp_in_Ca$MeanOD/LpMeanK, n = 2), limits = c(0, NA))
escalaY <- scale_y_continuous(breaks = pretty(Stats1XLp_in_Ca$MeanOD/LpMeanK, n = 2))
couleur <- scale_colour_manual(values = c( 'gray'))
puntos <- geom_point()
barras <- geom_errorbar(aes(ymin= (MeanOD-StErrOD)/LpMeanK, ymax= (MeanOD+StErrOD)/LpMeanK, width=.15))
paneles <- facet_wrap(~ Antibiotic,  ncol=1)
plotii = Plot1XLp + background  + titulos + couleur + paneles + puntos + barras
plotii

pdf("F_S13C.pdf", width=2.5, height=10)
print(plotii)
dev.off()
