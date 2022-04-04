library (ggplot2)
library (scales)
library(viridis)

#read the data
mydata <- read.csv ("TimeSeries_AntibioticExposure_1D_and_S2.csv", header=T)

#plot the time series for antibiotic exposure, community initially in the Ca-dominated state
dataplot <- subset(mydata, Initial =='Ca') #get the data from the initial state of interest

#NOTE: the data contains the colony counts for 10uL coming from a given dilution of the sample
      #to plot in CFU/mL, we need to multiply by 100 (from 10uL to 1mL), and by the dilution factor
      #this is done below when specifying what has to be plotted in the 'y' axis
plotReads <- ggplot(subset(dataplot, Replicate == '1' & Antibiotic == 'Chl'), aes(x= Day, y= Counts*10^(Dilution+2), group = Species, color = Species, fill = Species))

#other plot aesthetics
titulo <- ggtitle("")
xaxe <- xlab ("Day")
yaxe <- ylab("Density (CFU/ml") 
scale <- coord_cartesian(expand= c(0,0), xlim = c(0, 7), ylim = c(10^5,2*10^9))
yscale <- scale_y_log10()
couleur <- scale_colour_manual(values = c("#1AB24B","#FAA421"))#scale_color_viridis(discrete=TRUE, guide = FALSE)
panels <- facet_wrap(~ Antibiotic)
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) #, axis.title.x=element_blank() ) ## panel.background = element_rect(),# Background of plotting area

#curves for species populations in each replicate
Lineas1 <- geom_line(size = 0.5)
Lineas2 <- geom_line(data=subset(dataplot, Replicate == '2' & Antibiotic == 'Chl'), size = 0.5, alpha= 0.8)
Lineas3 <- geom_line(data=subset(dataplot, Replicate == '3' & Antibiotic == 'Chl'), size = 0.5, alpha= 0.6)

#for upper plot
ploti <- plotReads +  background  + Lineas1 + Lineas2 + Lineas3 +xaxe + yaxe + titulo + yscale + scale + couleur 
ploti

pdf(paste("AntibioticExposure_InitiallyCa-dominated.pdf"), width=6, height=4)
print(ploti)
dev.off()


#plot the time series for antibiotic exposure, community initially in the Lp-dominated state
dataplot <- subset(mydata, Initial =='Lp') #get the data from the initial state of interest

#NOTE: the data contains the colony counts for 10uL coming from a given dilution of the sample
#to plot in CFU/mL, we need to multiply by 100 (from 10uL to 1mL), and by the dilution factor
#this is done below when specifying what has to be plotted in the 'y' axis
plotReads <- ggplot(subset(dataplot, Replicate == '1' & Antibiotic == 'Chl'), aes(x= Day, y= Counts*10^(Dilution+2), group = Species, color = Species, fill = Species))

#other plot aesthetics
titulo <- ggtitle("")
xaxe <- xlab ("Day")
yaxe <- ylab("Density (CFU/ml") 
scale <- coord_cartesian(expand= c(0,0), xlim = c(0, 7), ylim = c(10^5,2*10^9))
yscale <- scale_y_log10()
couleur <- scale_colour_manual(values = c("#1AB24B","#FAA421"))#scale_color_viridis(discrete=TRUE, guide = FALSE)
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) #, axis.title.x=element_blank() ) ## panel.background = element_rect(),# Background of plotting area

#curves for species populations in each replicate
Lineas1 <- geom_line(size = 0.5)
Lineas2 <- geom_line(data=subset(dataplot, Replicate == '2' & Antibiotic == 'Chl'), size = 0.5, alpha= 0.8)
Lineas3 <- geom_line(data=subset(dataplot, Replicate == '3' & Antibiotic == 'Chl'), size = 0.5, alpha= 0.6)

#for upper plot
ploti <- plotReads +  background  + Lineas1 + Lineas2 + Lineas3 +xaxe + yaxe + titulo + yscale + scale + couleur 
ploti

pdf(paste("AntibioticExposure_InitiallyLp-dominated.pdf"), width=6, height=4)
print(ploti)
dev.off()
