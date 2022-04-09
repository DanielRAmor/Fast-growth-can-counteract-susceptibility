library(ggplot2)

#the data corresponds to the growth rates measured in Fig. S3
IMdata <- read.table ("GrowthRateVsChlConcentration.dat", header=T, sep = "\t", dec = ".")


## Bacteriostatics' model: Inverse Monod Growth
Rstatic <- function (R0, Bstat, I50) {
  Rstat <- R0 / (1 + Bstat / I50 ) 
  return (Rstat)
}

#parameter values for Ca
rf = 0.52
IC50_1 = 3.125   
# plot Ca's growth rate vs Chl concentration
Bstatic1<- seq(0,150, by= 0.01) 
plot (Bstatic1, Rstatic(rf, Bstatic1, IC50_1), xlab='Chl Concentration (ug/mL)', ylab='Growth rate (1/hr)', type = 'l', col='#1AB24B', lwd = 2, log='x', ylim =c(0.0,0.5), xlim = c(0.1, 100))
points(IMdata[1:12,1], IMdata[1:12,3], col = '#1AB24B', pch = 19)
arrows(x0=IMdata[1:12,1], y0=IMdata[1:12,3]-IMdata[1:12,4], x1=IMdata[1:12,1], y1=IMdata[1:12,3]+IMdata[1:12,4], code=3, angle=90, length=0.05, col="#1AB24B", lwd=0.4)
text(x=0.5, y=0.1, labels="IC50Ca = 3.125")
#export plot
pdf(paste("F_14A_right.pdf"), width=5, height=4.5)
plot (Bstatic1, Rstatic(rf, Bstatic1, IC50_1), xlab='Chl Concentration (ug/mL)', ylab='Growth rate (1/hr)', type = 'l', col='#1AB24B', lwd = 2, log='x', ylim =c(0.0,0.5), xlim = c(0.1, 100))
points(IMdata[1:12,1], IMdata[1:12,3], col = '#1AB24B', pch = 19)
arrows(x0=IMdata[1:12,1], y0=IMdata[1:12,3]-IMdata[1:12,4], x1=IMdata[1:12,1], y1=IMdata[1:12,3]+IMdata[1:12,4], code=3, angle=90, length=0.05, col="#1AB24B", lwd=0.4)
text(x=0.5, y=0.1, labels="IC50Ca = 3.125")
dev.off()

#parameter values for Lp
rs= 0.41
IC50_2 = 8.33 
# plot Lp's growth rate vs Chl concentrationrs 
Bstatic1<- seq(0,150, by= 0.01) 
plot (Bstatic1, Rstatic(rs, Bstatic1, IC50_2), xlab='Chl Concentration (ug/mL)', ylab='Growth rate (1/hr)', type = 'l', col='#FAA421', lwd = 2, log='x', ylim =c(0.0,0.5), xlim = c(0.1, 100))
points(IMdata[13:24,1], IMdata[13:24,3], col='#FAA421', pch = 19)
arrows(x0=IMdata[13:24,1], y0=IMdata[13:24,3]-IMdata[13:24,4], x1=IMdata[13:24,1], y1=IMdata[13:24,3]+IMdata[13:24,4], code=3, angle=90, length=0.05, col="#FAA421", lwd=0.4)
text(x=0.5, y=0.1, labels="IC50Lp = 8.33")
#export it
pdf(paste("F_14A_right.pdf"), width=5, height=4.5)
plot (Bstatic1, Rstatic(rs, Bstatic1, IC50_2), xlab='Chl Concentration (ug/mL)', ylab='Growth rate (1/hr)', type = 'l', col='#FAA421', lwd = 2, log='x', ylim =c(0.0,0.5), xlim = c(0.1, 100))
points(IMdata[13:24,1], IMdata[13:24,3], col='#FAA421', pch = 19)
arrows(x0=IMdata[13:24,1], y0=IMdata[13:24,3]-IMdata[13:24,4], x1=IMdata[13:24,1], y1=IMdata[13:24,3]+IMdata[13:24,4], code=3, angle=90, length=0.05, col="#FAA421", lwd=0.4)
text(x=0.5, y=0.1, labels="IC50Lp = 8.33")
dev.off()

