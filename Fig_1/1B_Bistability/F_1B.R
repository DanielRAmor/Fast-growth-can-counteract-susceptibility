library (ggplot2)
library (scales)

#Read Data
pares <- read.table ("Bistability-3Replicates.dat", header=T, sep = "\t", dec = ".")

#Plot time series for the Ca-dominated stable state
plotReads <- ggplot(data=pares, aes(x= Day, y= Abundance, group = interaction(Species, Replicate), color = Species, fill = Species))
#some plotting aesthetics
titulo <- ggtitle("")
xaxe <- xlab ("Day")
yaxe <- ylab("Abundance (CFU/ml)") 
yscale <- scale_y_log10()
scale <- coord_cartesian(expand= c(0,0), xlim = c(0, 5), ylim = c(10^4,2*10^9))
couleur <- scale_colour_manual(values = c("#1AB24B","#FAA421"))
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) #, axis.title.x=element_blank() ) ## panel.background = element_rect(),# Background of plotting area

Lineas1 <- geom_line(data=subset(pares,Expected_Initial_Lp == '10' & Replicate == '1'), aes(group = Species), size = 0.5)
Lineas2 <- geom_line(data=subset(pares,Expected_Initial_Lp == '10' & Replicate == '2'), aes(group = Species), size = 0.5, alpha= 0.8)
Lineas3 <- geom_line(data=subset(pares, Expected_Initial_Lp == '10' & Replicate == '3'), aes(group = Species), size = 0.5, alpha= 0.6)

#make the plot
ploti <- plotReads + background  +Lineas1 +Lineas2 +Lineas3 +xaxe + yaxe + titulo + yscale  + scale+ couleur 
ploti

pdf(paste("Fig_1B_Top.pdf"), width=4, height=2.5)
print(ploti)
dev.off()

#Lp-dominated stable state
plotReads <- ggplot(data=pares, aes(x= Day, y= Abundance, group = interaction(Species, Replicate), color = Species, fill = Species))
titulo <- ggtitle("")
xaxe <- xlab ("Day")
yaxe <- ylab("Abundance (CFU/ml)") 
yscale <- scale_y_log10()
scale <- coord_cartesian(expand= c(0,0), xlim = c(0, 5), ylim = c(10^4,2*10^9))
couleur <- scale_colour_manual(values = c("#1AB24B","#FAA421"))
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) #, axis.title.x=element_blank() ) ## panel.background = element_rect(),# Background of plotting area

Lineas1 <- geom_line(data=subset(pares,Expected_Initial_Lp == '80' & Replicate == '1'), aes(group = Species), size = 0.5)
Lineas2 <- geom_line(data=subset(pares,Expected_Initial_Lp == '80' & Replicate == '2'), aes(group = Species), size = 0.5, alpha= 0.8)
Lineas3 <- geom_line(data=subset(pares, Expected_Initial_Lp == '80' & Replicate == '3'), aes(group = Species), size = 0.5, alpha= 0.6)

ploti <- plotReads + background  +Lineas1 +Lineas2 +Lineas3 +xaxe + yaxe + titulo + yscale  + scale+ couleur 
ploti

pdf(paste("Fig_1B_Bottom.pdf"), width=4, height=2.5)
print(ploti)
dev.off()

