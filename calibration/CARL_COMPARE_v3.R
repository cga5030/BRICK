## PRODUCE TRACE PLOTS FOR CALIBRATION 

rm(list=ls())                        # Clear all previous variables
graphics.off()                       # clear plots

load("CARL_calib_MCMC_Complete_0614_2e7_14Jun2021.RData")

t.beg <- proc.time()

## Set up a filename for saving RData images along the way
configure <- 'COMPARE_v3_0715_'
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename.saveprogress <- paste(configure,today,'.RData',sep='')

amcmc.out_experts  <- readRDS(file = "CARL_calib_amcmc_out_Expert_0609_4e7_09Jun2021.rds",refhook = NULL)
amcmc.out_complete <- readRDS(file = "CARL_calib_amcmc_out_Complete_0614_2e7_14Jun2021.rds",refhook = NULL)

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

cexlab = 2

print("begin Antarctic PDFs")
filename.compare1<- paste('../output_calibration/',configure,'Set01','.png',sep='')
png(filename=filename.compare1, width=1920, height=1080, units ="px")
par(mfrow=c(3,5))
for (pp in sequence1){
  # experts
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.lab=cexlab)
  
  # complete
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=3,col="blue", cex.lab=cexlab)
  
  # priors:
  x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
  lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
        col="black", lwd=1, lty=2)
}
dev.off()

print("begin rest of PDFs")
filename.compare2<- paste('../output_calibration/',configure,'Set02','.png',sep='')
png(filename=filename.compare2, width=1920, height=1080, units ="px")
par(mfrow=c(3,6))
for (pp in sequence2){
  # experts
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.lab=cexlab)
  
  # complete
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=3,col="blue", cex.lab=cexlab)
  
  # priors:
  x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
  lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
        col="black", lwd=1, lty=2)
}

for (pp in sequence3){
  # experts
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.lab=cexlab)
  
  # complete
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=3,col="blue", cex.lab=cexlab)
  
  # priors:
  x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
  lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
        col="black", lwd=1, lty=2)
  }

for (pp in sequence4){
  # experts
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.lab=cexlab)
  
  # complete
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=3,col="blue", cex.lab=cexlab)
  
  # priors:
  x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
  lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
        col="black", lwd=1, lty=2)
}

for (pp in sequence5){
  # experts
  plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
       type="l",xlab=parnames[pp], ylab='PDF', lwd=3, col="red", yaxt='n', cex.lab=cexlab)
  
  # complete
  lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=3,col="blue", cex.lab=cexlab)
  
  # priors:
  x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
  lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
        col="black", lwd=1, lty=2)
}
dev.off()


# save results in case you need to revisit later
save.image(file=filename.saveprogress)
t.end <- proc.time()
time.elapsed <- t.end - t.beg
print("line 128")
print(time.elapsed)
print("TEDDY IS A CUTE PUPPER")

