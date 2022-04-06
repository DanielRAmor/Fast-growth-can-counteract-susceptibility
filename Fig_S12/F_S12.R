# Analisys on experimental data to assess ecological interactions via culturing each species 
#in the supernatant of either species in the absence of antibiotics.

library(reshape)
library(stringr)
library(dplyr)
library(plotrix)
library(ggplot2)
library(matrixStats)

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

#OD Ca cells in either Supernatant, No Antibiotics
Day2CaEco <- read.table ('Day2_Ca_in_BothSupernatants.dat', header=F, sep = "\t", dec = ".")
Day2CaEco <- Day2CaEco[,2:13] 
#Columns 1-3 are for 50% supernatant of cells from 30X dilution. 
#Columns 4-6 are for 50% supernatant of cells from high dispersal
#Columns 7-9 are for 50% supernatant of cells from low dispersal
#Columns 10-12 are for 50% supernatant of BM media (no antibiotics, no cells)
colnames(Day2CaEco) <- c('30Xa', '30Xb','30Xc','HighDa','HighDb','HighDc','LowDa','LowDb','LowDc','BMSa','BMSb', 'BMSc')
#Rows A-C correspond to supernatant from Ca cells
#Rows D & E contain no meaningful information (30X dilution of Ca cells into half (101.5uL) volume of fresh media)
#Rows F-H correspond to supernatant from Lp cells
rownames(Day2CaEco) <- c('CaSa', 'CaSb', 'CaSc', 'none1', 'none2', 'LpSa', 'LpSb', 'LpSc')
#export processed data (used in Fig S13 to get the Mean Abundance of Ca in the absence of antibiotics)
write.csv(Day2CaEco, 'CaEco.csv')

#OD Lp cells in either Supernatant, No Antibiotics
Day2LpEco <- read.table ('Day2_Lp_in_BothSupernatants.dat', header=F, sep = "\t", dec = ".")
Day2LpEco <- Day2LpEco[,2:13] 
#Eliminate non-meaningful rows
Day2LpEco <- Day2LpEco[-c(4,5),]
#Columns 1-3 are for 50% supernatant of cells from 30X dilution. 
#Columns 4-6 are for 50% supernatant of cells from high dispersal
#Columns 7-9 are for 50% supernatant of cells from low dispersal
#Columns 10-12 are for 50% supernatant of BM media (no antibiotics, no cells)
colnames(Day2LpEco) <- c('30Xa', '30Xb','30Xc','HighDa','HighDb','HighDc','LowDa','LowDb','LowDc','BMSa','BMSb', 'BMSc')
#Rows A-C correspond to supernatant from Ca cells
#Rows D & E contain no meaningful information (30X dilution of Lp cells into half (101.5uL) volume of fresh media)
#Rows F-H correspond to supernatant from Lp cells
Day2LpEco[13] = c('CaS1', 'CaS2', 'CaS3', 'LpS1', 'LpS2', 'LpS3')
colnames(Day2LpEco)[13] <- c('SSpecies')
#export processed data (used in Fig S13 to get the Mean Abundance of Lp in the absence of antibiotics)
write.csv(Day2LpEco, 'LpEco.csv')
LpEco <-  as.data.frame( melt(Day2LpEco,  id.vars=c('SSpecies'),  value.name="OD")) #formatting
#add Replicate indexing
LpEco[,4] <- substr(LpEco[,1], nchar(LpEco[,1]),nchar(LpEco[,1])) # from LpEco, read characters between the first and last number (reported by nchar)
LpEco[,5] <- substr(LpEco[,2], nchar(as.character(LpEco[,2])),nchar(as.character(LpEco[,2])))
#remove replicate indexing from column 1 and 2 (useful later when computing mean values)
LpEco[,1] <- substr(LpEco[,1], 1,nchar(LpEco[,1])-1)
LpEco[,2] <- substr(LpEco[,2], 1,nchar(as.character(LpEco[,2]))-1)
colnames(LpEco) <- c('SSpecies','SType','OD', 'ReplicateSSpecies', 'ReplicateSType')



MeanEcoCa <- data.frame(SType = c('30X','HighD'),SSpecies=c('Ca','Lp'), MeanOD=c(0.1,0.5),StErrOD=c(0.01,0.02))

MeanEcoCa[1,1:2] <- c('30X','Ca') 
MeanEcoCa[1,3] <- mean(as.matrix(Day2CaEco[1:3,1:3])) - BMbackground
MeanEcoCa[1,4] <- sd(c(as.matrix(Day2CaEco[1:3,1:3]))) / sqrt(9)

MeanEcoCa[2,1:2] <- c('HighD','Ca') 
MeanEcoCa[2,3] <- mean(as.matrix(Day2CaEco[1:3,4:6])) - BMbackground
MeanEcoCa[2,4] <- sd(c(as.matrix(Day2CaEco[1:3,4:6]))) / sqrt(9)

MeanEcoCa[3,1:2] <- c('LowD','Ca') 
MeanEcoCa[3,3] <- mean(as.matrix(Day2CaEco[1:3,7:9])) - BMbackground
MeanEcoCa[3,4] <- sd(c(as.matrix(Day2CaEco[1:3,7:9]))) / sqrt(9)

MeanEcoCa[4,1:2] <- c('NoCells','Ca')
MeanEcoCa[4,3] <- mean(as.matrix(Day2CaEco[1:3,10:12])) - BMbackground
MeanEcoCa[4,4] <- sd(c(as.matrix(Day2CaEco[1:3,10:12]))) / sqrt(9)

MeanEcoCa[5,1:2] <- c('30X','Lp')   # wells H3, H4 and H5 were not transferred (Viaflo failed the transfer). Note that I eliminate row H from the analysis
MeanEcoCa[5,3] <- mean(as.matrix(Day2CaEco[6:7,1:3])) - BMbackground
MeanEcoCa[5,4] <- sd(c(as.matrix(Day2CaEco[6:7,1:3]))) / sqrt(6)

MeanEcoCa[6,1:2] <- c('HighD','Lp')
MeanEcoCa[6,3] <- mean(as.matrix(Day2CaEco[6:7,4:6])) - BMbackground
MeanEcoCa[6,4] <- sd(c(as.matrix(Day2CaEco[6:7,4:6]))) / sqrt(6)

MeanEcoCa[7,1:2] <- c('LowD','Lp')
MeanEcoCa[7,3] <- mean(as.matrix(Day2CaEco[6:8,7:9])) - BMbackground
MeanEcoCa[7,4] <- sd(c(as.matrix(Day2CaEco[6:8,7:9]))) / sqrt(9)

MeanEcoCa[8,1:2] <- c('NoCells','Lp')
MeanEcoCa[8,3] <- mean(as.matrix(Day2CaEco[6:8,10:12])) - BMbackground
MeanEcoCa[8,4] <- sd(c(as.matrix(Day2CaEco[6:8,10:12]))) / sqrt(9)



#Barplot
PlotEcoCa <- ggplot(data = MeanEcoCa, aes(x= SType, y= MeanOD, group= SSpecies, color = SSpecies))

background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_rect (colour="black", size=1.5), axis.text=element_text(size=15), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )

titulos<- labs(title="", x=expression(paste("Supernatant Type")), y=expression(paste("OD24hr")))
#escalaX <- scale_x_continuous(expand= c(0,0), limits=c(0,0.88), breaks=seq(0, 0.8, 0.2)) # coord_cartesian(expand= c(0,0), xlim = c(0, 1.3)) 

couleur <- scale_colour_manual(values = c("#1AB24B", 'gray'))#"#FAA421"))

puntos <- geom_point()
barras <- geom_errorbar(aes(ymin= MeanOD-StErrOD, ymax= MeanOD+StErrOD), width=.15)

plotii = PlotEcoCa + background  + titulos + couleur + puntos + barras
plotii

pdf("F_S12A.pdf", width=6, height=4)
print(plotii)
dev.off()



MeanEcoLp <- data.frame(SType = c('30X','HighD'),SSpecies=c('Ca','Lp'), MeanOD=c(0.1,0.5),StErrOD=c(0.01,0.02))
#Note to myself: there has to be a way to make the next code lines more compact, maybe using factors
MeanEcoLp[1,1:2] <- c('30X','Ca') 
MeanEcoLp[1,3] <- mean(LpEco[LpEco[1]== 'CaS' & LpEco[2]== '30X',3]) - BMbackground
MeanEcoLp[1,4] <- sd(LpEco[LpEco[1]== 'CaS' & LpEco[2]== '30X',3]) / sqrt(9)

MeanEcoLp[2,1:2] <- c('HighD','Ca') 
MeanEcoLp[2,3] <- mean(LpEco[LpEco[1]== 'CaS' & LpEco[2]== 'HighD',3]) - BMbackground
MeanEcoLp[2,4] <- sd(LpEco[LpEco[1]== 'CaS' & LpEco[2]== 'HighD',3]) / sqrt(9)

MeanEcoLp[3,1:2] <- c('LowD','Ca') 
MeanEcoLp[3,3] <- mean(LpEco[LpEco[1]== 'CaS' & LpEco[2]== 'LowD',3]) - BMbackground
MeanEcoLp[3,4] <- sd(LpEco[LpEco[1]== 'CaS' & LpEco[2]== 'LowD',3]) / sqrt(9)

MeanEcoLp[4,1:2] <- c('NoCells','Ca') 
MeanEcoLp[4,3] <- mean(LpEco[LpEco[1]== 'CaS' & LpEco[2]== 'BMS',3]) - BMbackground
MeanEcoLp[4,4] <- sd(LpEco[LpEco[1]== 'CaS' & LpEco[2]== 'BMS',3]) / sqrt(9)

MeanEcoLp[5,1:2] <- c('30X','Lp') 
MeanEcoLp[5,3] <- mean(LpEco[LpEco[1]== 'LpS' & LpEco[2]== '30X',3]) - BMbackground
MeanEcoLp[5,4] <- sd(LpEco[LpEco[1]== 'LpS' & LpEco[2]== '30X',3]) / sqrt(9)

MeanEcoLp[6,1:2] <- c('HighD','Lp') 
MeanEcoLp[6,3] <- mean(LpEco[LpEco[1]== 'LpS' & LpEco[2]== 'HighD',3]) - BMbackground
MeanEcoLp[6,4] <- sd(LpEco[LpEco[1]== 'LpS' & LpEco[2]== 'HighD',3]) / sqrt(9)

MeanEcoLp[7,1:2] <- c('LowD','Lp') 
MeanEcoLp[7,3] <- mean(LpEco[LpEco[1]== 'LpS' & LpEco[2]== 'LowD',3]) - BMbackground
MeanEcoLp[7,4] <- sd(LpEco[LpEco[1]== 'LpS' & LpEco[2]== 'LowD',3]) / sqrt(9)

MeanEcoLp[8,1:2] <- c('NoCells','Lp') 
MeanEcoLp[8,3] <- mean(LpEco[LpEco[1]== 'LpS' & LpEco[2]== 'BMS',3]) - BMbackground
MeanEcoLp[8,4] <- sd(LpEco[LpEco[1]== 'LpS' & LpEco[2]== 'BMS',3]) / sqrt(9)


#Barplot
PlotEcoLp <- ggplot(data = MeanEcoLp, aes(x= SType, y= MeanOD, group= SSpecies, color = SSpecies))

background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_rect (colour="black", size=1.5), axis.text=element_text(size=15), axis.title=element_text(size=18), title=element_text(size=22)) #, axis.title.x=element_blank() )

titulos<- labs(title="", x=expression(paste("Supernatant Type")), y=expression(paste("OD24hr")))
#escalaX <- scale_x_continuous(expand= c(0,0), limits=c(0,0.88), breaks=seq(0, 0.8, 0.2)) # coord_cartesian(expand= c(0,0), xlim = c(0, 1.3)) 

couleur <- scale_colour_manual(values = c("#FAA421", 'gray'))#

puntos <- geom_point()
barras <- geom_errorbar(aes(ymin= MeanOD-StErrOD, ymax= MeanOD+StErrOD), width=.15)

plotii = PlotEcoLp + background  + titulos + couleur + puntos + barras
plotii

pdf("F_S12B.pdf", width=6, height=4)
print(plotii)
dev.off()
