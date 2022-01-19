rm(list=ls())

# CODE TO PLOT PDF, CDF, SURVIVAL function
# library(ggplot2)
# library(survival)
library(Hmisc)

# load and format slr data
sl.2100.complete <- readRDS(file="complete_0118_slr_rcp26.rds")
sl.2100.complete <- sl.2100.complete[,ncol(sl.2100.complete)]

sl.2100.standard <- readRDS(file="standard_0118_slr_rcp26.rds")
sl.2100.standard <- sl.2100.standard[,ncol(sl.2100.standard)]

sl.2100.experts <- readRDS(file="experts_0118_slr_rcp26.rds")
sl.2100.experts <- sl.2100.experts[,ncol(sl.2100.experts)]

sl.2100.priors <- readRDS(file="priors_0118_slr_rcp26.rds")
sl.2100.priors <- sl.2100.priors[,ncol(sl.2100.priors)]

cdf.complete <- Ecdf(sl.2100.complete)
cdf.standard <- Ecdf(sl.2100.standard)
cdf.experts  <- Ecdf(sl.2100.experts)
cdf.priors   <- Ecdf(sl.2100.priors)
dev.off()

density.complete <- density(sl.2100.complete)
density.standard <- density(sl.2100.standard)
density.experts  <- density(sl.2100.experts)
density.priors   <- density(sl.2100.priors)

if(TRUE){
# dev.off()
par(mfcol= c(4,1),
    mar=c(0, 5, 0, 2) + 0.1
)
# par(mfrow=c(4,1))
mylwd = 2
mycexlab = 1.35

plot(density.standard$x, density.standard$y,
     xlim=c(0,5.5),
     axes = FALSE,
     type="l",
     lwd=mylwd, col="black",
     yaxt='n',
     ylab='',
     xlab=''
     )
axis(2,labels=FALSE)
axis(1,labels=FALSE)
title(ylab="Probability Density",cex.lab=mycexlab)
box()


lines(density.complete$x,density.complete$y,
      lwd=mylwd, col="red"
)

lines(density.experts$x, density.experts$y,
      lwd=mylwd, col="orange"
)

lines(density.priors$x, density.priors$y,
      lwd=1, col="black", lty="dashed")

## Plot cdf
standard.cdf.x <- cdf.standard$x
standard.cdf.y <- cdf.standard$y
plot(standard.cdf.x,standard.cdf.y, type="l",
     xlim=c(0,5.5),
     axes=FALSE,
     # main="",
     # xlab="Global mean sea level at 2100 [m]", ylab='Cumulative Probability', 
     xlab="",
     ylab="",
     lwd=mylwd, col="black"
     # , yaxt='n'
)
axis(2L,labels=c("0","0.2","0.4","0.6","0.8","1"),
     at=c(0,0.2,0.4,0.6,0.8,1))
axis(1,labels=FALSE)
title(ylab="Cumulative Probability",cex.lab=mycexlab)
box()

complete.cdf.x <- cdf.complete$x
complete.cdf.y <- cdf.complete$y
lines(complete.cdf.x, complete.cdf.y, type="l",
      lwd=mylwd, col="red")

experts.cdf.x <- cdf.experts$x
experts.cdf.y <- cdf.experts$y
lines(experts.cdf.x, experts.cdf.y, type="l",
      lwd=mylwd, col="orange")

priors.cdf.x <- cdf.priors$x
priors.cdf.y <- cdf.priors$y
lines(priors.cdf.x, priors.cdf.y, type="l",
      lwd=1, col="black", lty="dashed")

## Plot survival
plot(standard.cdf.x,1-standard.cdf.y, type="l",
     xlim=c(0,5.5),
     axes=FALSE,
     lwd=mylwd, col="black",
     # xlab="Global mean sea level at 2100 [m]", ylab='',
     xlab="",ylab=""
     )
axis(2L, labels=c("0","0.2","0.4","0.6","0.8","1"),
     at=c(0,0.2,0.4,0.6,0.8,1))
axis(1)
title(ylab="Survival Probability",cex.lab=mycexlab)
box()

lines(complete.cdf.x,1-complete.cdf.y,
      lwd=mylwd, col="red",
)

lines(experts.cdf.x,1-experts.cdf.y,
      lwd=mylwd,col="orange"
)

lines(priors.cdf.x,1-priors.cdf.y,
      lwd=1,col="black",lty="dashed")

mtext("Global mean sea level at 2100 [m]",side=1,line=3)



plot(1, type = "n", axes=FALSE, xlab="", ylab=""
     # ,main="Global mean sea level at 2100 [m]"
     )
legend(x = "center",inset = 0,
       legend = c("Wide priors combined with expert assessment and data" , 
                  "Wide priors combined with data",
                  "Wide priors combined with expert assessment",
                  "Wide priors"
       ),
       lwd=c(mylwd,mylwd,mylwd,1),
       lty=c("solid","solid","solid","dashed"),
       bty='n', cex=1.2,
       col=c("red","black","orange","black"), 
       horiz = FALSE,
       
       )
}

