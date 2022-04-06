library (ggplot2)
library(viridis)

#read the data
M <- read.table ("Growth-In-Buffered-Media.dat", header=T, sep = "\t", dec = ".")

#compute Fold Growth as ODf / ODi
M [5] = M[3]/0.0009  # 0.0009 is the initial OD computed as 
                     #OD of the prepared monocultures (adjusted to OD=0.063) multiplied by the inoculation factor 3uL in 210uL
                     # since the inial OD was lower than the detection limit
colnames(M)[5] <- c('Fold_Growth')

#plot the fold growth
ploti <- ggplot(data = M)
background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_rect (colour="black", size=1.0), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12)) 
titulos<- labs(title="", x=expression(paste("Initial pH")), y=expression(paste("Fold Growth (ODf/ODi)")))
escalaX <- scale_x_continuous(breaks=seq(2.5, 10.5, 1)) 
escalaY <- scale_y_continuous(breaks=seq(0, 1400, 400))
CaShade <- geom_ribbon(data=subset(M, Species == "Ca"), aes(x = pH, ymax = Fold_Growth, ymin = 0.0), alpha = .4, fill = "#22B14B")
LpShade <- geom_ribbon(data=subset(M, Species == "Lp"), aes(x = pH, ymax = Fold_Growth, ymin = 0.0), alpha = .37, fill = "#FAA421")
plotii = ploti + background  + titulos + CaShade + LpShade + escalaX + escalaY 
plotii

pdf(paste(  "F_S8", ".pdf"), width=2.7, height=2.5)
print(plotii)
dev.off()