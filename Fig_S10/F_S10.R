library (ggplot2)
library (scales)
library(viridis)

#Read the data
pares <- read.table ("Chl_lowDispersal_AdditionalReplicates.dat", header=T, sep = "\t", dec = ".")

#Plot the abundance time series for the community initially dominated by Lp
plotReads <- ggplot(data=subset(pares, Initial == 'Lp'), aes(x= Day, y= Counts1*10^(Dilution+2), group = Species, color = Species, fill = Species))
titulo <- ggtitle("")
xaxe <- xlab ("Time (Days)")
yaxe <- ylab("Abundance (CFU/mL)") 
yscale <- scale_y_log10()
scale <- coord_cartesian(expand= c(0,0), xlim = c(0, 7))
Lineas1 <- geom_line(aes(x=Day, y=Counts1*10^(Dilution+2)),size = 0.5) #The abundance is computed as the CFU counts per the dilution factor applied to the sample
#plus  extra 10^2 factor resulting from plating 10uL (and not 1mL), while we want to plot the abundance in CFU/mL units 
Lineas2 <- geom_line(aes(x=Day, y=Counts2*10^(Dilution+2)),size = 0.5, alpha= 0.8)
Lineas3 <- geom_line(aes(x=Day, y=Counts3*10^(Dilution+2)),size = 0.5, alpha= 0.6)
couleur <- scale_colour_manual(values = c("#1AB24B","#FAA421"))
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) #, axis.title.x=element_blank() ) ## panel.background = element_rect(),# Background of plotting area

#make F_S10, panel on top
ploti <- plotReads + background  +Lineas1 + Lineas2 + Lineas3 +xaxe + yaxe + titulo + yscale  + scale+ couleur 
ploti
#export it
pdf(paste("F_S10_top.pdf"), width=4, height=2.5)
print(ploti)
dev.off()

#Plot the abundance time series for the community initially dominated by Ca
plotReads <- ggplot(data=subset(pares, Initial == 'Ca'), aes(x= Day, y= Counts1*10^(Dilution+2), group = Species, color = Species, fill = Species))
titulo <- ggtitle("")
xaxe <- xlab ("Time (Days)")
yaxe <- ylab("Abundance (CFU/mL)") 
yscale <- scale_y_log10()
scale <- coord_cartesian(expand= c(0,0), xlim = c(0, 7))
Lineas1 <- geom_line(aes(x=Day, y=Counts1*10^(Dilution+2)),size = 0.5) #The abundance is computed as the CFU counts per the dilution factor applied to the sample
                                                                      #plus  extra 10^2 factor resulting from plating 10uL (and not 1mL), while we want to plot the abundance in CFU/mL units 
Lineas2 <- geom_line(aes(x=Day, y=Counts2*10^(Dilution+2)),size = 0.5, alpha= 0.8)
Lineas3 <- geom_line(aes(x=Day, y=Counts3*10^(Dilution+2)),size = 0.5, alpha= 0.6)
couleur <- scale_colour_manual(values = c("#1AB24B","#FAA421"))
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) #, axis.title.x=element_blank() ) ## panel.background = element_rect(),# Background of plotting area

#make F_S10, panel on bottom
ploti <- plotReads + background  +Lineas1 + Lineas2 + Lineas3 +xaxe + yaxe + titulo + yscale  + scale+ couleur 
ploti
#export it
pdf(paste("F_S10_bottom.pdf"), width=4, height=2.5)
print(ploti)
dev.off()


