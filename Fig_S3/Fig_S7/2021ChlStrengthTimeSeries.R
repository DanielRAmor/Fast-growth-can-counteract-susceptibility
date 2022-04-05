library (ggplot2)
library (scales)

#read the data
mydata <- read.table ("TimeSeries-Chl-Strength.dat", header=T, sep = "\t", dec = ".")

#plot time series for each antibiotic concentration and initial (Day <2) stable state
plotReads <- ggplot(subset(mydata), aes(x= Day, y= Counts1*10^(Dilution+2), group = Species, color = Species, fill = Species))
titulo <- ggtitle("Panel labels: antibiotic concentration (ug/mL) and dominant species at Day 0")
xaxe <- xlab ("Day")
yaxe <- ylab("Abundance (CFU/mL)") 
yscale <- scale_y_log10()
scale <- coord_cartesian(expand= c(0,0), xlim = c(0, 7), ylim = c(10^5,2*10^9))
Lineas1 <- geom_line(size = 0.5)
Lineas2 <- geom_line(aes(x= Day, y= Counts2*10^(Dilution+2), group = Species, color = Species), size = 0.5, alpha= 0.8)
Lineas3 <- geom_line(aes(x= Day, y= Counts3*10^(Dilution+2), group = Species, color = Species), size = 0.5, alpha= 0.6)
couleur <- scale_colour_manual(values = c("#1AB24B","#FAA421"))
couleur2 <- scale_fill_manual(values = c("#1AB24B","#FAA421"))
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=6)) 
paneles <- facet_wrap( Strength ~ Initial, ncol=2)

#make the plot
ploti <- plotReads +  background  +Lineas1 + Lineas2 + Lineas3 +xaxe + yaxe + titulo + yscale + scale+ couleur2 + couleur + paneles 
ploti
#export the plot
pdf(paste("F_S7.pdf"), width=4.5, height=6)
print(ploti)
dev.off()

