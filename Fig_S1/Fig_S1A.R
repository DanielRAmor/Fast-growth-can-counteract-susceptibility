library (ggplot2)
library (scales)

#read the data
pHseries <- read.table ("Bistability-pH-timeSeries.dat", header=T, sep = "\t", dec = ".")

#plot pH time series
plotpH <- ggplot(data=pHseries, aes(x= Day, y= pH, color = InitialLp, fill = InitialLp))

#plot aesthetics
titulo <- ggtitle("")
xaxe <- xlab ("Day")
yaxe <- ylab("pH") 
yscale <- scale_y_log10()
scale <- coord_cartesian(expand= c(0,0), xlim = c(0, 5), ylim = c(10^4,2*10^9))
couleur <- scale_colour_manual(values = c("#D30B4E","#1C75BC"))
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) #, axis.title.x=element_blank() ) ## panel.background = element_rect(),# Background of plotting area

#curves for each replicate
Lineas1 <- geom_line(data=subset(pHseries,InitialLp == 'i10' & Replicate == '1'), size = 0.5)
Lineas2 <- geom_line(data=subset(pHseries,InitialLp == 'i10' & Replicate == '2'),  size = 0.5, alpha= 0.8)
Lineas3 <- geom_line(data=subset(pHseries, InitialLp == 'i10' & Replicate == '3'), size = 0.5, alpha= 0.6)
Lineas4 <- geom_line(data=subset(pHseries,InitialLp == 'i80' & Replicate == '1'), size = 0.5)
Lineas5 <- geom_line(data=subset(pHseries,InitialLp == 'i80' & Replicate == '2'),  size = 0.5, alpha= 0.8)
Lineas6 <- geom_line(data=subset(pHseries, InitialLp == 'i80' & Replicate == '3'), size = 0.5, alpha= 0.6)

#make the plot
ploti <- plotpH + background  +Lineas1 +Lineas2 +Lineas3 +Lineas4 +Lineas5 +Lineas6 +xaxe + yaxe + titulo + couleur 
ploti

#export it
pdf(paste("F_S1A.pdf"), width=4, height=2.5)
print(ploti)
dev.off()
