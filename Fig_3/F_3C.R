#This code generates Fig. 3C
#The plot shows the per capita growth rate of each species as a function of its normalized abundance.

# Values of the maximum growth rates for species F and S
Rf <- 0.54
Rs <- 0.43
#Create the x axis (for the range and precision of the numerical data)
x0 = seq(0.0,1.0, by=0.00001)
# the slow grower is subject to logistic growth
y_e = Rs*(1-x0) #y_e is the per capita growth rate of species S

#plot the per capita growth rate of S
plot (x0, y_e, xlab='Abundance', ylab='Per capita growth rate', col = 'orange', type='l', lwd = 2, xlim=c(0.0001,1.0), ylim=c(0,0.6), log='x')
#intensity of the Allee effect on species F
A= 0.0001
#plot the per capita growth rate of F
lines (x0, (Rff*x0*(1-x0)/(A+x0)), col = 'chartreuse4', lwd = 1.5)

#export the plot
pdf('F_3C.pdf', width=4.5, height=4.5)
  #plot the per capita growth rate of S
  plot (x0, y_e, xlab='Abundance', ylab='Per capita growth rate', col = 'orange', type='l', lwd = 2, xlim=c(0.0001,1.0), ylim=c(0,0.6), log='x')
  #intensity of the Allee effect on species F
  A= 0.0001
  #plot the per capita growth rate of F
  lines (x0, (Rff*x0*(1-x0)/(A+x0)), col = 'chartreuse4', lwd = 1.5)
dev.off()
