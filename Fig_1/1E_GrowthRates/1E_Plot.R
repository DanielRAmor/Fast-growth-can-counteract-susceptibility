#This code reads the data and plots the time series for the growth of Ca and Lp in the absence of antibiotics
#and also visualizes their average maximum growth rates through the corresponding slopes

library (ggplot2)
library (plyr) 
library(scales)
library(varhandle)
library(viridis)

#read the data
G8all <- read.table ("NoAntibiotic_GrowthCurves.csv", header=T, sep = ",", dec = ".")
G8all <- G8all[, 2:(ncol(G8all))]
G8all[G8all$OD<0, 5] = 0.00001 #remove negative values after background OD correction

Means <- read.table ("Mean_Growth_Rates_and_Intercepts.csv", header=T, sep = ",", dec = ".")

#plot OD over time
ploti <- ggplot(data = G8all, aes(x= Time_hr, y = OD, group= interaction(Replicate, Species, Well), color= interaction (Replicate, Species), fill = Replicate))

#aesthetics
background <- theme_bw() + theme ( panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.line= element_line(colour="black"), panel.border = element_rect (colour="black", size=1.0), axis.text=element_text(size=9), axis.title=element_text(size=11), title=element_text(size=12), legend.position = 'none') 
titulos <- labs(title="", x=expression(paste("Time (hr)")), y=expression(paste("OD")))
escalaY = scale_y_continuous(limits= c(0.001 , 1), trans='log10',  labels = trans_format("log10", math_format(10^.x)))
escalaX = scale_x_continuous(limits = c(0, 20))
lineas <- geom_line()
colores <- scale_colour_manual(values = c("#1AB24B", "#168936", "#24CE54", '#FAA51A', '#D88B1A', '#BA7D25')) 

alphaCa <- Means[1,1]  #average growth rate over the 24 replicates (from the growth rate fitting analysis)
alphaLp <-  Means[1,2]
betaCa <- Means[1,3]  #average y-axis intercept (from the growth rate fitting analysis)
betaLp <- Means[1,4]
  
fitlineCa <- geom_segment(aes(x = 1.5, xend = 4.5, y = exp(betaCa + alphaCa*1.5 + 1.5), yend = exp(betaCa + alphaCa*4.5+ 1.5)), color = "#1AB24B")
fitlineLp <- geom_segment(aes(x = 4.5, xend = 7.5, y = exp(betaLp + alphaLp*4.5 - 1.5), yend = exp(betaLp + alphaLp*7.5 - 1.5)), color = '#FAA51A')

plotii = ploti + background  + titulos + lineas + escalaY + escalaX + fitlineCa + fitlineLp + colores# + meanCa + meanLp
plotii

pdf(paste(  "F_1E", ".pdf"), width=6, height=3)
print(plotii)
dev.off()

