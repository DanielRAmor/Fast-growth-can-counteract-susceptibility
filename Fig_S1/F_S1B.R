library (ggplot2)
library (scales)

#Read the data
pares <- read.table ("Chl_Exposure-pH-TimeSeries.dat", header=T, sep = "\t", dec = ".")

#Time series for the pH
plotpH <- ggplot(data=subset(pares), aes(x= Day, y= pH1, group = Initial, color = Initial, fill = Initial))

#plot aesthetics
titulo <- ggtitle("")
xaxe <- xlab ("Day")
yaxe <- ylab("pH") 
scale <- coord_cartesian(expand= c(0,0), xlim = c(0, 7))
couleur <- scale_colour_manual(values = c("#1C75BC","#D30B4E"))
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12))

#curves for each replicate
Lineas1 <- geom_line( aes(x=Day, y=pH1),size = 0.5)
Lineas2 <- geom_line(aes(x=Day, y=pH2),size = 0.5, alpha= 0.8)
Lineas3 <- geom_line(aes(x=Day, y=pH3),size = 0.5, alpha= 0.6)

#make the plot
ploti <- plotpH + background  +Lineas1 + Lineas2 + Lineas3 +xaxe + yaxe + titulo  + scale+ couleur 
ploti

pdf(paste("F_S1_B.pdf"), width=4, height=2.5)
print(ploti)
dev.off()
