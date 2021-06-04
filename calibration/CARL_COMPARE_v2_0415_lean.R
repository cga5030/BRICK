## PRODUCE TRACE PLOTS FOR CALIBRATION 

rm(list=ls())                        # Clear all previous variables
graphics.off()                       # clear plots

load("CARL_calib_MCMC_complete_v2_0415_1e715Apr2021.RData")

t.beg <- proc.time()

## Set up a filename for saving RData images along the way
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename.saveprogress <- paste('COMPARE_v2_0415_lean_',today,'.RData',sep='')

amcmc.out_experts  <- readRDS(file = "CALIB_amcmc_out_0415_4e7.rds",refhook = NULL)
amcmc.out_complete <- readRDS(file = "CALIB_amcmc_out_0415_1e7.rds",refhook = NULL)

## specify windows for marginal distributions
chain_experts <- amcmc.out_experts$samples[5e6:4e7,] # 5 million to 40 million
chain_complete <- amcmc.out_complete$samples[1e6:10e6,] # 1 million to 10 million

## get posterior densities
chain_experts_densities <- replicate(length(parnames), vector("list", 7), simplify = FALSE)
chain_complete_densities <- chain_experts_densities

for (pp in 1:length(parnames)){
  chain_experts_densities[[pp]] <- density(chain_experts[,pp])
}

for (pp in 1:length(parnames)){
  chain_complete_densities[[pp]] <- density(chain_complete[,pp])
}

## draw pdfs
sequence1 = c(26,23,24,22,20,14,15,17,19,21,16,18,25,27,28)
sequence2 = c(5,6,7,8)
sequence3 = c(9,10,11,12,13)
sequence4 = c(1,2,3,4)
sequence5 = c(29,30,31,32,33)

#added 6/1/2021
sequence6 = c(21,9,11)

print("begin Antarctic PDFs")
png(filename="!!COMPARE_0415_PDFs_1.png", width=1280, height=720, units ="px")
par(mfrow=c(3,5))
for (pp in sequence1){
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.label=3.0)
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=3,col="blue", cex.label=3.0)
}
dev.off()

print("begin rest of PDFs")
png(filename="!!COMPARE_0415_PDFs_1.png", width=1280, height=720, units ="px")
par(mfrow=c(3,6))
for (pp in sequence2){
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.label=3.0)
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=3,col="blue", cex.label=3.0)
}

for (pp in sequence3){
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.label=3.0)
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=3,col="blue", cex.label=3.0)
}

for (pp in sequence4){
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.label=3.0)
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=3,col="blue", cex.label=3.0)
}

for (pp in sequence5){
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.label=3.0)
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=3,col="blue", cex.label=3.0)
}
dev.off()

## create list of priors
my_priors <- data.frame(matrix(NA, nrow = length(parnames), ncol = 2))
names(my_priors) = c('lower_bound','upper_bound')
my_priors$lower_bound <- bound.lower
my_priors$upper_bound <- bound.upper

if(FALSE){
  if(TRUE){
    mylwd = 2
    graphics.off()  
    par(mfrow=c(1,2))
    pp=sequence6[1]
    plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
         type="l",
         xlab= expression(paste(kappa["DAIS "],"[\u00B0C"^-1,"]")), 
         ylab='Probability density', lwd=mylwd, col="red",
         yaxt='n',
         ylim=c(0,33),
         xlim=c(0.025,0.085)
        )
    lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,col="blue")
    x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
    lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
          col="black", lwd=1, lty=2)
    
    legend(x="topright", 
           legend=c("Wide physically reasonable ranges",
                    "Wide physical priors combined with expert assessments",
                    "Expert+data posteriors"),
           col=c("black","blue","red"), lwd=c(1,mylwd,mylwd), 
           lty=c(2,1,1))
    
    pp=sequence6[2]
    plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
         type="l",
         xlab=expression(paste(italic(c)["GIS "],"[m \u00B0C"^-1,"]")), 
         ylab='Probatility density', lwd=mylwd, col="red", 
         yaxt='n',
         ylim=c(0,0.7),
         xlim=c(-4,-0.001))
    
    
    lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,col="blue")
    x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
    lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
          col="black", lwd=1, lty=2)
    
    legend(x="topright", 
           legend=c("Wide physically reasonable ranges",
                    "Wide physical priors combined with expert assessments",
                    "Expert+data posteriors"),
           col=c("black","blue","red"), lwd=c(1,mylwd,mylwd), 
           lty=c(2,1,1))
    
    # pp=sequence6[3]
    # plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
    #      type="l",xlab=parnames[pp], ylab='PDF', lwd=mylwd, col="red", 
    #      yaxt='n',
    #      ylim=c(0,2400)
    #      )
    # lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,col="blue")
    # x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
    # lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
    #       col="black", lwd=1, lty=2)
    # 
    # legend(x="topright", legend=c("physical priors","expert priors","expert+data posteriors"),
    #        col=c("black","blue","red"), lwd=c(1,mylwd,mylwd), 
    #        lty=c(2,1,1))
  }
}

# save results in case you need to revisit later
save.image(file=filename.saveprogress)
t.end <- proc.time()
time.elapsed <- t.end - t.beg
print("line 28")
print(time.elapsed)
print("TEDDY IS A CUTE PUPPER")

