# Analisys on experimental data to assess cross-protective interactions via culturing the Ca species 
#in Lp's supernatant in the presence of antibiotics, diluted by a factor 8X,.

library(reshape)
library(stringr)
library(dplyr)
library(plotrix)
library(ggplot2)

#OD BACKGROUND OF MEDIA WITH AND WITHOUT ANTIBIOTICS
Day2Background <- read.table ('Day2_NoCells_AntibioticBackground.dat', header=F, sep = "\t", dec = ".")
Day2Background <- Day2Background[,2:13] 
#Columns 1-3 are for media + 16X antibiotic concentration (where it could be seen visually which antibiotics are less soluble, 
#potentially interfering more with background OD even at lower concentrations).
#Columns 4-6 are for media + 1X antibiotic concentration. 
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
CaEco <- read.csv('CaEco.csv')
CaK = CaEco[11:13]# Carrying capacity as OD in absence of antibiotics
CaK = CaK[-(4:5),]
CaK[7:12,1] = CaK[1:6,2]
CaK[13:18,1] = CaK[1:6,3]
CaK = CaK[,-(2:3)]
CaMeanK = mean(CaK)
CaErrK = sd(CaK)/sqrt(18)


#OD Ca cells in 1X antibiotics, Lp Supernatant
Day2Ca_1X_in_LpS <- read.table ('Day2_Ca_in_1over8_LpSupernatant.dat', header=F, sep = "\t", dec = ".")
colnames(Day2Ca_1X_in_LpS) <- c('Antibiotic','30X_1', '30X_2','30X_3','HighD_1','HighD_2','HighD_3','LowD_1','LowD_2','LowD_3','OnlyAnt_1','OnlyAnt_2', 'OnlyAnt_3')
Day2Ca_1X_in_LpS[1] <- c('NaN', 'NaN', 'NaN', 'Amp', 'NaN', 'Car', 'NaN', 'NaN')
#Subtract background OD for each antibiotic
for (i in 2:13){
  Day2Ca_1X_in_LpS[i] <- Day2Ca_1X_in_LpS[i] - Day2Background[14]
}

Ca_1X_in_LpS <-  as.data.frame( melt(Day2Ca_1X_in_LpS,  id.vars=c('Antibiotic'),  value.name="OD")) #formatting
#add replicate indexing
Ca_1X_in_LpS[,4] <- substr(Ca_1X_in_LpS[,2], nchar(as.character(Ca_1X_in_LpS[,2])),nchar(as.character(Ca_1X_in_LpS[,2])))
#remove replicate indexing from column 2
Ca_1X_in_LpS[,2] <- substr(Ca_1X_in_LpS[,2], 1,nchar(as.character(Ca_1X_in_LpS[,2]))-2)
colnames(Ca_1X_in_LpS) <- c('Antibiotic','Supernatant', 'OD', 'Replicate')

#Ca cells in 1X Lp supernatant
#Mean and StError for each Antibiotic&Supernatant condition
AntibClass <- unlist(Ca_1X_in_LpS[1] %>% unique)  # unique is creating a list, but we need a vector, so we also use unlist()
SupernClass <- unlist(dplyr::distinct(Ca_1X_in_LpS[2])) #alternative way of getting unique values
Stats1XCa_in_Lp <- data.frame(Antibiotic = c('Chl','Chl'),Supernatant=c('HighD','LowD'), MeanOD=c(0.1,0.5),StErrOD=c(0.01,0.02))
replicates <- 3  
for (anti in 1:length(AntibClass)){
  for (super in 1:length(SupernClass)){
    Stats1XCa_in_Lp[super+(anti-1)*length(SupernClass),1] <- AntibClass[anti]
    Stats1XCa_in_Lp[super+(anti-1)*length(SupernClass),2] <- SupernClass[super]
    Stats1XCa_in_Lp[super+(anti-1)*length(SupernClass),3] <- mean(Ca_1X_in_LpS[Ca_1X_in_LpS[,1]== AntibClass[anti] & Ca_1X_in_LpS[,2]== SupernClass[super],3])
    Stats1XCa_in_Lp[super+(anti-1)*length(SupernClass),4] <-   sd(Ca_1X_in_LpS[Ca_1X_in_LpS[,1]== AntibClass[anti] & Ca_1X_in_LpS[,2]== SupernClass[super],3])/sqrt(replicates)
  }
} 


#Plot
Plot1XCa <- ggplot(subset(Stats1XCa_in_Lp, Antibiotic != 'NaN'), aes(x= Supernatant, y= MeanOD/CaMeanK, group= Supernatant))
background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_rect (colour="black", size=1.5), axis.text=element_text(size=10), axis.title=element_text(size=12), title=element_text(size=22)) 
titulos<- labs(title="", x=expression(paste("Supernatant Type")), y=expression(paste("24hr abundance relative to maximum")))
escalaY <- scale_y_continuous(breaks = pretty(Stats1XCa_in_Lp$MeanOD/CaMeanK, n = 3), limits = c(0, NA))
couleur <- scale_colour_manual(values = c( 'gray'))
puntos <- geom_point()
barras <- geom_errorbar(aes(ymin= (MeanOD-StErrOD)/CaMeanK, ymax= (MeanOD+StErrOD)/CaMeanK, width=.15))
paneles <- facet_wrap(~ Antibiotic,  ncol=1)
plotii = Plot1XCa + background  + titulos + couleur + paneles + puntos + barras + escalaY
plotii

pdf("F_S13F.pdf", width=2.5, height=3)
print(plotii)
dev.off()
