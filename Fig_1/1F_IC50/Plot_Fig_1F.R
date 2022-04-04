#this code plots panel 1F
library(ggplot2)
#read the data
S1 <- read.csv('F_1F_data.csv',header = TRUE) 

# plot IC50 based on antibiotic concentration that inhibits half the maximal growth (final abundance)
ploti <- ggplot(data = S1, aes(x= Species, y = AverageIC50, color= Species, fill = Species))

#plot aesthetics
background <- theme_gray() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_blank(), panel.border = element_blank(), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) 
titulos <- labs(title="", x=expression(paste("Species")), y=expression(paste("Susceptibility (Conc. for ", 0.5, ' rel. growth)' )))
barras <- geom_bar(stat="identity")
err1 <- geom_errorbar(aes(ymin = AverageIC50 - StErrIC50, ymax =AverageIC50 + StErrIC50), size=1, color = 'gray2', width = 0.2) #position = position_dodge(width = 0.45),
paneles <- facet_wrap(~ Antibiotic, scales = 'free_y')
colores <- scale_fill_manual(values=c('#1AB24B', '#FAA421')) ##try this one!!
colores2 <- scale_color_manual(values=c('#1AB24B', '#FAA421')) ##try this one!!

#make the plot
plotii = ploti + background  + titulos + barras +  paneles + colores + colores2 + err1
plotii

#export it
pdf(paste(  "F_1F", ".pdf"), width=5, height=4)
print(plotii)
dev.off()

