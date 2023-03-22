rm(list=ls())
graphics.off()

t.beg <- proc.time()

## Load pipeline files
t.hind <- readRDS(file="priors_0426_t_hind.rds")
t.proj <- readRDS(file="priors_0426_t_proj.rds")
t.paleo <- readRDS(file="priors_0426_t_paleo.rds")

## NOTE: THERE IS ONLY ONE SCENARIO: "high temperature" from SEJ2018.
## here, "rcp26" refers to "high temperature" from SEJ2018.

#####################################################################
## Priors
priors.ais.paleo.05 <- readRDS(file="priors_0426_ais_paleo_05.rds")
priors.ais.paleo.50 <- readRDS(file="priors_0426_ais_paleo_50.rds")
priors.ais.paleo.95 <- readRDS(file="priors_0426_ais_paleo_95.rds")

priors.gsic.hind <- readRDS(file="priors_0426_gsic_hind.rds")
priors.te.hind <- readRDS(file="priors_0426_te_hind.rds")
priors.gis.hind <- readRDS(file="priors_0426_gis_hind.rds")
priors.ais.hind <- readRDS(file="priors_0426_ais_hind.rds")
priors.temp.hind <- readRDS(file="priors_0426_temp_hind.rds")
priors.ocheat.hind <- readRDS(file="priors_0426_ocheat_hind.rds")
priors.gsl.hind <- readRDS(file="priors_0426_gsl_hind.rds")

# for projections
priors.slr.rcp26 <- readRDS(file="priors_0426_slr_rcp26.rds")
priors.te.rcp26 <- readRDS(file="priors_0426_te_rcp26.rds") 
priors.gis.rcp26 <- readRDS(file="priors_0426_gis_rcp26.rds")
priors.gsic.rcp26 <- readRDS(file="priors_0426_gsic_rcp26.rds")
priors.ais.rcp26 <- readRDS(file="priors_0426_ais_rcp26.rds")
priors.temp.rcp26 <- readRDS(file="priors_0426_temp_rcp26.rds")
priors.ocheat.rcp26 <- readRDS(file="priors_0426_ocheat_rcp26.rds")
#####################################################################

#####################################################################
## Experts
experts.ais.paleo.05 <- readRDS(file="experts_0426_ais_paleo_05.rds")
experts.ais.paleo.50 <- readRDS(file="experts_0426_ais_paleo_50.rds")
experts.ais.paleo.95 <- readRDS(file="experts_0426_ais_paleo_95.rds")

experts.gsic.hind <- readRDS(file="experts_0426_gsic_hind.rds")
experts.te.hind <- readRDS(file="experts_0426_te_hind.rds")
experts.gis.hind <- readRDS(file="experts_0426_gis_hind.rds")
experts.ais.hind <- readRDS(file="experts_0426_ais_hind.rds")
experts.temp.hind <- readRDS(file="experts_0426_temp_hind.rds")
experts.ocheat.hind <- readRDS(file="experts_0426_ocheat_hind.rds")
experts.gsl.hind <- readRDS(file="experts_0426_gsl_hind.rds")

# for projections
experts.slr.rcp26 <- readRDS(file="experts_0426_slr_rcp26.rds")
experts.te.rcp26 <- readRDS(file="experts_0426_te_rcp26.rds") 
experts.gis.rcp26 <- readRDS(file="experts_0426_gis_rcp26.rds")
experts.gsic.rcp26 <- readRDS(file="experts_0426_gsic_rcp26.rds")
experts.ais.rcp26 <- readRDS(file="experts_0426_ais_rcp26.rds")
experts.temp.rcp26 <- readRDS(file="experts_0426_temp_rcp26.rds")
experts.ocheat.rcp26 <- readRDS(file="experts_0426_ocheat_rcp26.rds")
#####################################################################

#####################################################################
## Standard
standard.ais.paleo.05 <- readRDS(file="standard_0426_ais_paleo_05.rds")
standard.ais.paleo.50 <- readRDS(file="standard_0426_ais_paleo_50.rds")
standard.ais.paleo.95 <- readRDS(file="standard_0426_ais_paleo_95.rds")

standard.gsic.hind <- readRDS(file="standard_0426_gsic_hind.rds")
standard.te.hind <- readRDS(file="standard_0426_te_hind.rds")
standard.gis.hind <- readRDS(file="standard_0426_gis_hind.rds")
standard.ais.hind <- readRDS(file="standard_0426_ais_hind.rds")
standard.temp.hind <- readRDS(file="standard_0426_temp_hind.rds")
standard.ocheat.hind <- readRDS(file="standard_0426_ocheat_hind.rds")
standard.gsl.hind <- readRDS(file="standard_0426_gsl_hind.rds")

# for projections
standard.slr.rcp26 <- readRDS(file="standard_0426_slr_rcp26.rds")
standard.te.rcp26 <- readRDS(file="standard_0426_te_rcp26.rds") 
standard.gis.rcp26 <- readRDS(file="standard_0426_gis_rcp26.rds")
standard.gsic.rcp26 <- readRDS(file="standard_0426_gsic_rcp26.rds")
standard.ais.rcp26 <- readRDS(file="standard_0426_ais_rcp26.rds")
standard.temp.rcp26 <- readRDS(file="standard_0426_temp_rcp26.rds")
standard.ocheat.rcp26 <- readRDS(file="standard_0426_ocheat_rcp26.rds")
#####################################################################

#####################################################################
## Complete
complete.ais.paleo.05 <- readRDS(file="complete_0426_ais_paleo_05.rds")
complete.ais.paleo.50 <- readRDS(file="complete_0426_ais_paleo_50.rds")
complete.ais.paleo.95 <- readRDS(file="complete_0426_ais_paleo_95.rds")

complete.gsic.hind <- readRDS(file="complete_0426_gsic_hind.rds")
complete.te.hind <- readRDS(file="complete_0426_te_hind.rds")
complete.gis.hind <- readRDS(file="complete_0426_gis_hind.rds")
complete.ais.hind <- readRDS(file="complete_0426_ais_hind.rds")
complete.temp.hind <- readRDS(file="complete_0426_temp_hind.rds")
complete.ocheat.hind <- readRDS(file="complete_0426_ocheat_hind.rds")
complete.gsl.hind <- readRDS(file="complete_0426_gsl_hind.rds")

# for projections
complete.slr.rcp26 <- readRDS(file="complete_0426_slr_rcp26.rds")
complete.te.rcp26 <- readRDS(file="complete_0426_te_rcp26.rds") 
complete.gis.rcp26 <- readRDS(file="complete_0426_gis_rcp26.rds")
complete.gsic.rcp26 <- readRDS(file="complete_0426_gsic_rcp26.rds")
complete.ais.rcp26 <- readRDS(file="complete_0426_ais_rcp26.rds")
complete.temp.rcp26 <- readRDS(file="complete_0426_temp_rcp26.rds")
complete.ocheat.rcp26 <- readRDS(file="complete_0426_ocheat_rcp26.rds")
#####################################################################



## Other useful scripts
source('../Useful/colorblindPalette.R') # Get nice plotting colors: mycol array
source('../Useful/MultipleOutput.R')    # defines the useful ":=" operator

## And set the IPCC RCP colors
col26 <- c(0, 0, 255)/255
col45 <- c(121, 188, 255)/255
col60 <- c(255, 130, 45)/255
col85 <- c(255, 0, 0)/255

## Set colors to use for control model, and observations/experimental model,
## in all plots. This is indexed within "mycol", from "colorblindPalette.R".
colmod <- 2
colobs <- 11

## Where would you like to save the plots?
plotdir <- '../plots/'
# main_path=getwd()

##==============================================================================
logl.ar1 = function(r,sigma1,rho1,eps1=0) # default obs error is 0
{
  n = length(r) # r is the residuals
  if(length(eps1)==1) eps1 = rep(eps1,n)
  
  logl=0
  if(n>1) {
    w = r[2:n] - rho1*r[1:(n-1)] # this process whitens the residuals
    logl = logl + sum(dnorm(w,sd=sqrt((sigma1)^2+(eps1[c(-1)])^2),log=TRUE)) # add in the sum of
    # density of the whitened residuals with a standard deviation of the
    # variance and the obs. errors
  }
  return(logl)
}
##==============================================================================

##==============================================================================
##==============================================================================
## Supplemental Figure
## -- Hindcast, model vs observational data

## Initialize arrays for the output
# priors
priors.slr.05 = rep(NA,length(t.hind));		  priors.slr.50 = rep(NA,length(t.hind));	  	priors.slr.95 = rep(NA,length(t.hind))
priors.gsic.05 = rep(NA,length(t.hind));		priors.gsic.50 = rep(NA,length(t.hind));		priors.gsic.95 = rep(NA,length(t.hind))
priors.gis.05 = rep(NA,length(t.hind));		  priors.gis.50 = rep(NA,length(t.hind));		  priors.gis.95 = rep(NA,length(t.hind))
priors.te.05 = rep(NA,length(t.hind));			priors.te.50 = rep(NA,length(t.hind));			priors.te.95 = rep(NA,length(t.hind))
priors.temp.05 = rep(NA,length(t.hind));		priors.temp.50 = rep(NA,length(t.hind));		priors.temp.95 = rep(NA,length(t.hind))
priors.ocheat.05 = rep(NA,length(t.hind));  priors.ocheat.50 = rep(NA,length(t.hind));	priors.ocheat.95 = rep(NA,length(t.hind))

# experts
experts.slr.05 = rep(NA,length(t.hind));		experts.slr.50 = rep(NA,length(t.hind));	  experts.slr.95 = rep(NA,length(t.hind))
experts.gsic.05 = rep(NA,length(t.hind));		experts.gsic.50 = rep(NA,length(t.hind));		experts.gsic.95 = rep(NA,length(t.hind))
experts.gis.05 = rep(NA,length(t.hind));		experts.gis.50 = rep(NA,length(t.hind));		experts.gis.95 = rep(NA,length(t.hind))
experts.te.05 = rep(NA,length(t.hind));			experts.te.50 = rep(NA,length(t.hind));			experts.te.95 = rep(NA,length(t.hind))
experts.temp.05 = rep(NA,length(t.hind));		experts.temp.50 = rep(NA,length(t.hind));		experts.temp.95 = rep(NA,length(t.hind))
experts.ocheat.05 = rep(NA,length(t.hind)); experts.ocheat.50 = rep(NA,length(t.hind));	experts.ocheat.95 = rep(NA,length(t.hind))

# standard
standard.slr.05 = rep(NA,length(t.hind));		  standard.slr.50 = rep(NA,length(t.hind));	  	standard.slr.95 = rep(NA,length(t.hind))
standard.gsic.05 = rep(NA,length(t.hind));		standard.gsic.50 = rep(NA,length(t.hind));		standard.gsic.95 = rep(NA,length(t.hind))
standard.gis.05 = rep(NA,length(t.hind));		  standard.gis.50 = rep(NA,length(t.hind));		  standard.gis.95 = rep(NA,length(t.hind))
standard.te.05 = rep(NA,length(t.hind));			standard.te.50 = rep(NA,length(t.hind));			standard.te.95 = rep(NA,length(t.hind))
standard.temp.05 = rep(NA,length(t.hind));		standard.temp.50 = rep(NA,length(t.hind));		standard.temp.95 = rep(NA,length(t.hind))
standard.ocheat.05 = rep(NA,length(t.hind));  standard.ocheat.50 = rep(NA,length(t.hind));	standard.ocheat.95 = rep(NA,length(t.hind))

# complete
complete.slr.05 = rep(NA,length(t.hind));		  complete.slr.50 = rep(NA,length(t.hind));	  	complete.slr.95 = rep(NA,length(t.hind))
complete.gsic.05 = rep(NA,length(t.hind));		complete.gsic.50 = rep(NA,length(t.hind));		complete.gsic.95 = rep(NA,length(t.hind))
complete.gis.05 = rep(NA,length(t.hind));		  complete.gis.50 = rep(NA,length(t.hind));		  complete.gis.95 = rep(NA,length(t.hind))
complete.te.05 = rep(NA,length(t.hind));			complete.te.50 = rep(NA,length(t.hind));			complete.te.95 = rep(NA,length(t.hind))
complete.temp.05 = rep(NA,length(t.hind));		complete.temp.50 = rep(NA,length(t.hind));		complete.temp.95 = rep(NA,length(t.hind))
complete.ocheat.05 = rep(NA,length(t.hind));  complete.ocheat.50 = rep(NA,length(t.hind));	complete.ocheat.95 = rep(NA,length(t.hind))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.hind)){
  # priors  
  c(priors.slr.05[t],    priors.slr.50[t],    priors.slr.95[t])				:= quantile(priors.gsl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.gsic.05[t],   priors.gsic.50[t],   priors.gsic.95[t])			:= quantile(priors.gsic.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.gis.05[t],    priors.gis.50[t],    priors.gis.95[t])				:= quantile(priors.gis.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.te.05[t],     priors.te.50[t],     priors.te.95[t])				:= quantile(priors.te.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.temp.05[t],   priors.temp.50[t],   priors.temp.95[t])			:= quantile(priors.temp.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.ocheat.05[t], priors.ocheat.50[t], priors.ocheat.95[t])	  := quantile(priors.ocheat.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # experts 
  c(experts.slr.05[t],    experts.slr.50[t],    experts.slr.95[t])		:= quantile(experts.gsl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.gsic.05[t],   experts.gsic.50[t],   experts.gsic.95[t])		:= quantile(experts.gsic.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.gis.05[t],    experts.gis.50[t],    experts.gis.95[t])		:= quantile(experts.gis.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.te.05[t],     experts.te.50[t],     experts.te.95[t])			:= quantile(experts.te.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.temp.05[t],   experts.temp.50[t],   experts.temp.95[t])		:= quantile(experts.temp.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.ocheat.05[t], experts.ocheat.50[t], experts.ocheat.95[t])	:= quantile(experts.ocheat.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # standard  
  c(standard.slr.05[t],    standard.slr.50[t], standard.slr.95[t])				:= quantile(standard.gsl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.gsic.05[t],   standard.gsic.50[t], standard.gsic.95[t])			:= quantile(standard.gsic.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.gis.05[t],    standard.gis.50[t], standard.gis.95[t])				:= quantile(standard.gis.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.te.05[t],     standard.te.50[t], standard.te.95[t])					:= quantile(standard.te.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.temp.05[t],   standard.temp.50[t], standard.temp.95[t])			:= quantile(standard.temp.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.ocheat.05[t], standard.ocheat.50[t], standard.ocheat.95[t])	:= quantile(standard.ocheat.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # complete  
  c(complete.slr.05[t],    complete.slr.50[t], complete.slr.95[t])				:= quantile(complete.gsl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.gsic.05[t],   complete.gsic.50[t], complete.gsic.95[t])			:= quantile(complete.gsic.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.gis.05[t],    complete.gis.50[t], complete.gis.95[t])				:= quantile(complete.gis.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.te.05[t],     complete.te.50[t], complete.te.95[t])					:= quantile(complete.te.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.temp.05[t],   complete.temp.50[t], complete.temp.95[t])			:= quantile(complete.temp.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.ocheat.05[t], complete.ocheat.50[t], complete.ocheat.95[t])	:= quantile(complete.ocheat.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  
}

begyear = t.hind[1]
endyear = t.hind[length(t.hind)]
mod.time= begyear:endyear
begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time = length(mod.time)

source('../R/compute_indices.R')  # function to determine the model and
# data indices for comparisons

## Source the data for hindcast comparisons
source('../calibration/DOECLIM_readData.R')
source('../calibration/GSIC_readData.R')
source('../calibration/SIMPLE_readData.R')
source('../calibration/DAIS_readData.R')
source('../calibration/TE_readData.R')

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.temp,midx.ocheat,midx.gis,midx.gsic,midx.sl)
names(midx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )
oidx.all        = list(oidx.temp,oidx.ocheat,oidx.gis,oidx.gsic,oidx.sl)
names(oidx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.temp, obs.ocheat, obs.gis, obs.gsic, obs.sl)
names(obs.all) = c(    "temp"  , "ocheat"  , "gis"  , "gsic"  , "sl" )
obs.err.all        = list( obs.temp.err, obs.ocheat.err, obs.gis.err, obs.gsic.err, obs.sl.err)
names(obs.err.all) = c(    "temp"      , "ocheat"      , "gis"      , "gsic"      , "sl"      )

## Set the indices for normalization that are consistent with each data set
ind.norm.data = data.frame(
  c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
  c(which(mod.time==1850),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
  c(which(mod.time==1870),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

## Also normalize the available data to this time period.
ibeg=which(obs.temp.time==mod.time[1])
iend=which(obs.temp.time==mod.time[20])
obs.temp.norm = obs.temp - mean(obs.temp[ibeg:iend])

ibeg=which(obs.ocheat.time==mod.time[ind.norm[1]])
iend=which(obs.ocheat.time==mod.time[ind.norm[length(ind.norm)]])
obs.ocheat.norm = obs.ocheat - mean(obs.ocheat[ibeg:iend])

ibeg=which(obs.gsic.time==mod.time[ind.norm[1]])
iend=which(obs.gsic.time==mod.time[ind.norm[length(ind.norm)]])
obs.gsic.norm = obs.gsic #- mean(obs.gsic[ibeg:iend]) # GSIC does not need normalized - already is normalized to 1960

ibeg=which(obs.gis.time==mod.time[ind.norm[1]])
iend=which(obs.gis.time==mod.time[ind.norm[length(ind.norm)]])
obs.gis.norm = obs.gis - mean(obs.gis[ibeg:iend])

ibeg=which(obs.sl.time==mod.time[ind.norm[1]])
iend=which(obs.sl.time==mod.time[ind.norm[length(ind.norm)]])
obs.sl.norm = obs.sl - mean(obs.sl[ibeg:iend])

## TE trends
# read in TE_readData.R -- "trends.te"

## Precal windows 5-?:
## Last "precalibration window" is 1993-2010 mean trend, from the IPCC AR5 Ch13
## (Page 1151), for AIS SLR contribution: 0.27 +- 0.11 mm/year
## Note that model output is in meters SLE and these trends are mm, so a
## conversion is necessary.

trends.ais <- c(0.27 , 0.08 , 0.40 )/1000   # m/year (after the /1000)
trends.err <- c(0.11 , 0.185, 0.205)/1000   # m/year (after the /1000)
trends.2up <- trends.ais+2*trends.err
trends.2dn <- trends.ais-2*trends.err
ind.trends <- mat.or.vec( length(trends.ais), 2)
ind.trends[1,] <- c(which(date==-7) , which(date==10)) # 1993-2010
ind.trends[2,] <- c(which(date==-8) , which(date== 1)) # 1992-2001
ind.trends[3,] <- c(which(date== 2) , which(date==11)) # 2002-2011

## Precal window 4:
## Adding observational constraint
estimate.SLE.rate <- abs(-71/360)/1000
time.years <- 2002-1992      # using the midpoint of the 19-year interval
mid.cum.SLE_2002 <- estimate.SLE.rate*time.years
i1992 <- which(date==-8)

estimate.SLE.rate.error <- abs(-53/360)/1000     #1-sigma error
estimate.SLE.error <- sqrt(time.years)*estimate.SLE.rate.error #1-sigma error
# (*sqrt(10) because 10 years of potentially accumulated error:
#  total error^2 = year 1 error^2 + year 2 error^2 + ... year 10 error^2
#                = 10*year X error^2)
SE2_2002 <- estimate.SLE.error*2 #2-sigma error

positive_2SE <- mid.cum.SLE_2002 + SE2_2002 # Add the 2 standard error to the mean value
negative_2SE <- mid.cum.SLE_2002 - SE2_2002 # Subtract the 2 standard error to the mean value

## Precal windows 1-3:
## from Shaffer (2014). modified by Kelsey
upper.wind <- c(7.4, -6.9, -1.25, positive_2SE) # Windows 2-3 from Kelsey, Window 1 from DeConto and Pollard 2016
lower.wind <- c(3.6, -15.8, -4.0, negative_2SE)
#upper.wind <- c(6.0, -6.9, -1.25, positive_2SE) # Windows 1-3 from Kelsey
#lower.wind <- c(1.8, -15.8, -4.0, negative_2SE)
#upper.wind <- c(5.5, -8 , -2, positive_2SE) # Windows 1-3 fFrom Shaffer 2014, p 1809
#lower.wind <- c(2.5, -17, -4, negative_2SE)

windows <- matrix(c(lower.wind, upper.wind), nrow = length(upper.wind), ncol=2)
obs.targets <- (windows[,2]+windows[,1])*.5  # middle of window = obs to compare model to
obs.err.dais <- (windows[,2]-windows[,1])*.5      # half-width of window = uncertainty
obs.err.dais <- 0.5*obs.err.dais                       # assume all windows are 2*stdErr (last two actually are)

## Create a vector with each observation year
## 120kyr, 20Kyr, 6kyr (before present), 2002, and 1993 (first year of the IPCC trend)
obs.years <- c(120000, 220000, 234000, 240002)

##
## 5-95% CI of hindcasts, with obs in there
##

# pdf(paste(plotdir,'hindcasts_with_noise_and_normalization.pdf',sep=''),width=7,height=6,colormodel='cmyk')

# png(paste(plotdir,'standard_hindcasts.png',sep=''), width=958, height=1080, units ="px")

if(FALSE){
  
  n.sig = 2         # how many sigma to plot around the obs?
  # layout(cbind(c(1,4,7),c(2,5,7),c(3,6,7)))
  layout(cbind(c(1,3,5,6),c(2,4,5,6)))
  # par(mai=c(.5,.7,.2,.08)) #c(bottom, left, top, right)
  par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) # omi = c(0,0,0,2) for right side space
  plotall <- FALSE
  
  
  if(plotall){
    # >>> SURFACE TEMPERATURE <<<
    plot(mod.time[midx.temp], complete.temp.50[midx.temp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1850,2016), ylim=c(-.3,1.5), cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    # mtext(side=2, text='Surface temperature\n[deg C]', line=2.3, cex=.9);
    mtext(side=3, text=expression(bold('a) Surface temperature')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[midx.temp],rev(mod.time[midx.temp])), c(complete.temp.95,rev(complete.temp.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.temp.time[oidx.temp], obs.temp.norm[oidx.temp], type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.temp.time[oidx.temp],rev(obs.temp.time[oidx.temp])),
            c(obs.temp.norm[oidx.temp]+n.sig*obs.temp.err[oidx.temp],rev(obs.temp.norm[oidx.temp]-n.sig*obs.temp.err[oidx.temp])),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
    #legend(1839,1.65,c("5-95% range, modeled","2-sigma range, observations"),
    #       col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])), lwd=2, bty='n', cex=1.2)
    
    # >>> OCEAN HEAT <<<
    itmp=midx.ocheat[1]:nrow(complete.ocheat.hind)
    plot(mod.time[itmp], complete.ocheat.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1950,2016), ylim=c(-50,50), cex.lab=1.2, cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    # mtext(side=2, text='b) Ocean heat uptake', line=2.3, cex=.9);
    mtext(side=3, text=expression(bold('b) Ocean heat uptake')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(complete.ocheat.95[itmp],rev(complete.ocheat.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.ocheat.time, obs.ocheat.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.ocheat.time,rev(obs.ocheat.time)), c(obs.ocheat.norm+n.sig*obs.ocheat.err,rev(obs.ocheat.norm-n.sig*obs.ocheat.err)),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  }
  # >>> GSIC <<<
  itmp=midx.gsic[1]:nrow(complete.gsic.hind)
  plot(mod.time[itmp], complete.gsic.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.01,.04), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Glaciers and\nsmall ice caps [m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  # mtext(side=3, text=expression(bold(' c')), line=.25, cex=.9, adj=0);
  mtext(side=3, text=expression(bold(' a')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+21, .85 * par('usr')[4], labels = 'Glaciers and Ice Caps', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(complete.gsic.95[itmp],rev(complete.gsic.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gsic.time, obs.gsic.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  # polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  # >>> GIS <<<
  plot(mod.time, complete.gis.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+6, .85 * par('usr')[4], labels = 'GIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(complete.gis.95,rev(complete.gis.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gis.time, obs.gis.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  # polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  # >>> TE <<<
  x1971=seq(trends.te[1,4],trends.te[1,5])
  c1971=mean(x1971); yc1971=complete.te.50[which(mod.time==mean(x1971))]
  lo1971 = yc1971+(trends.te[1,2]/1000)*(x1971-c1971)
  hi1971 = yc1971+(trends.te[1,3]/1000)*(x1971-c1971)
  y1971 = yc1971+(trends.te[1,1]/1000)*(x1971-c1971)

  x1993=seq(trends.te[2,4],trends.te[2,5])
  c1993=mean(x1993); yc1993=complete.te.50[which(mod.time==mean(x1993))]
  lo1993 = yc1993+(trends.te[2,2]/1000)*(x1993-c1993)
  hi1993 = yc1993+(trends.te[2,3]/1000)*(x1993-c1993)
  y1993 = yc1993+(trends.te[2,1]/1000)*(x1993-c1993)

  plot(mod.time, complete.te.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.04,.06), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Thermal expansion\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+19, .8 * par('usr')[4], labels = 'Thermal expansion', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(complete.te.95,rev(complete.te.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(x1971,y1971, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1971,y1971, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  polygon(c(x1971,rev(x1971)), c(lo1971,rev(hi1971)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  lines(x1993,y1993, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1993,y1993, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(x1993,rev(x1993)), c(lo1993,rev(hi1993)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  
  # >>> TOTAL SLR <<<
  plot(mod.time, complete.slr.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+21, .75 * par('usr')[4], labels = 'GMSL', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(complete.slr.95,rev(complete.slr.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  
  # >>> AIS PALEO, SMOOTHED <<<
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], complete.ais.paleo.50[ipaleo], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=.25, cex=.9, adj=0);
  text(.95 * par('usr')[1], .8 * par('usr')[4], labels = 'AIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(complete.ais.paleo.95[ipaleo],rev(complete.ais.paleo.05[ipaleo])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # legend(-92000,12, c("5-95% range, model" , "2-sigma range, observations"), lwd=2, bty='n', cex=1.2,
  #        col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]) , rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])) )
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top",inset = 0,
         legend = c("median, model" ,
                    "5-95% range, model",
                    "observations"#,
                    # "2-sigma range, observations"
                    ),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,NA,8), bty='n', cex=1.2,
         pch=c(NA,NA,20,NA),
         col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
               rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
               rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
               ), 
               horiz = FALSE)
  
}
# dev.off()
##==============================================================================
##==============================================================================


##==============================================================================
##==============================================================================
## Supplemental Figure
## -- Sea-level projections to 2100

# ncdata <- nc_open(filename.brick.magicc)
# t.proj = ncvar_get(ncdata, 'time_proj')
# standard.slr.rcp26 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
# slr.rcp45 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
# slr.rcp85 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
# standard.te.rcp26 = ncvar_get(ncdata, 'TE_RCP26')
# te.rcp45 = ncvar_get(ncdata, 'TE_RCP45')
# te.rcp85 = ncvar_get(ncdata, 'TE_RCP85')
# standard.gis.rcp26 = ncvar_get(ncdata, 'GIS_RCP26')
# gis.rcp45 = ncvar_get(ncdata, 'GIS_RCP45')
# gis.rcp85 = ncvar_get(ncdata, 'GIS_RCP85')
# standard.gsic.rcp26 = ncvar_get(ncdata, 'GSIC_RCP26')
# gsic.rcp45 = ncvar_get(ncdata, 'GSIC_RCP45')
# gsic.rcp85 = ncvar_get(ncdata, 'GSIC_RCP85')
# standard.ais.rcp26 = ncvar_get(ncdata, 'AIS_RCP26')
# ais.rcp45 = ncvar_get(ncdata, 'AIS_RCP45')
# ais.rcp85 = ncvar_get(ncdata, 'AIS_RCP85')
# standard.temp.rcp26 = ncvar_get(ncdata, 'temp_RCP26')
# temp.rcp45 = ncvar_get(ncdata, 'temp_RCP45')
# temp.rcp85 = ncvar_get(ncdata, 'temp_RCP85')
# standard.ocheat.rcp26 = ncvar_get(ncdata, 'ocheat_RCP26')
# ocheat.rcp45 = ncvar_get(ncdata, 'ocheat_RCP45')
# ocheat.rcp85 = ncvar_get(ncdata, 'ocheat_RCP85')
# nc_close(ncdata)

#transpose:

# priors
priors.slr.rcp26 <- t(priors.slr.rcp26)
priors.te.rcp26 <- t(priors.te.rcp26)
priors.gis.rcp26 <- t(priors.gis.rcp26)
priors.gsic.rcp26 <- t(priors.gsic.rcp26)
priors.ais.rcp26 <- t(priors.ais.rcp26)
priors.temp.rcp26 <- t(priors.temp.rcp26)
priors.ocheat.rcp26 <- t(priors.ocheat.rcp26)

# experts
experts.slr.rcp26 <- t(experts.slr.rcp26)
experts.te.rcp26 <- t(experts.te.rcp26)
experts.gis.rcp26 <- t(experts.gis.rcp26)
experts.gsic.rcp26 <- t(experts.gsic.rcp26)
experts.ais.rcp26 <- t(experts.ais.rcp26)
experts.temp.rcp26 <- t(experts.temp.rcp26)
experts.ocheat.rcp26 <- t(experts.ocheat.rcp26)


# standard
standard.slr.rcp26 <- t(standard.slr.rcp26)
standard.te.rcp26 <- t(standard.te.rcp26)
standard.gis.rcp26 <- t(standard.gis.rcp26)
standard.gsic.rcp26 <- t(standard.gsic.rcp26)
standard.ais.rcp26 <- t(standard.ais.rcp26)
standard.temp.rcp26 <- t(standard.temp.rcp26)
standard.ocheat.rcp26 <- t(standard.ocheat.rcp26)

# complete
complete.slr.rcp26 <- t(complete.slr.rcp26)
complete.te.rcp26 <- t(complete.te.rcp26)
complete.gis.rcp26 <- t(complete.gis.rcp26)
complete.gsic.rcp26 <- t(complete.gsic.rcp26)
complete.ais.rcp26 <- t(complete.ais.rcp26)
complete.temp.rcp26 <- t(complete.temp.rcp26)
complete.ocheat.rcp26 <- t(complete.ocheat.rcp26)


## Initialize arrays for the output
# priors
priors.slr.rcp26.05 = rep(NA,length(t.proj)); priors.slr.rcp26.50 = rep(NA,length(t.proj)); priors.slr.rcp26.95 = rep(NA,length(t.proj))
priors.ais.rcp26.05 = rep(NA,length(t.proj)); priors.ais.rcp26.50 = rep(NA,length(t.proj)); priors.ais.rcp26.95 = rep(NA,length(t.proj))
priors.gis.rcp26.05 = rep(NA,length(t.proj)); priors.gis.rcp26.50 = rep(NA,length(t.proj)); priors.gis.rcp26.95 = rep(NA,length(t.proj))
priors.gsic.rcp26.05 = rep(NA,length(t.proj)); priors.gsic.rcp26.50 = rep(NA,length(t.proj)); priors.gsic.rcp26.95 = rep(NA,length(t.proj))
priors.te.rcp26.05 = rep(NA,length(t.proj)); priors.te.rcp26.50 = rep(NA,length(t.proj)); priors.te.rcp26.95 = rep(NA,length(t.proj))
priors.temp.rcp26.05 = rep(NA,length(t.proj)); priors.temp.rcp26.50 = rep(NA,length(t.proj)); priors.temp.rcp26.95 = rep(NA,length(t.proj))
priors.ocheat.rcp26.05 = rep(NA,length(t.proj)); priors.ocheat.rcp26.50 = rep(NA,length(t.proj)); priors.ocheat.rcp26.95 = rep(NA,length(t.proj))

# experts
experts.slr.rcp26.05 = rep(NA,length(t.proj)); experts.slr.rcp26.50 = rep(NA,length(t.proj)); experts.slr.rcp26.95 = rep(NA,length(t.proj))
experts.ais.rcp26.05 = rep(NA,length(t.proj)); experts.ais.rcp26.50 = rep(NA,length(t.proj)); experts.ais.rcp26.95 = rep(NA,length(t.proj))
experts.gis.rcp26.05 = rep(NA,length(t.proj)); experts.gis.rcp26.50 = rep(NA,length(t.proj)); experts.gis.rcp26.95 = rep(NA,length(t.proj))
experts.gsic.rcp26.05 = rep(NA,length(t.proj)); experts.gsic.rcp26.50 = rep(NA,length(t.proj)); experts.gsic.rcp26.95 = rep(NA,length(t.proj))
experts.te.rcp26.05 = rep(NA,length(t.proj)); experts.te.rcp26.50 = rep(NA,length(t.proj)); experts.te.rcp26.95 = rep(NA,length(t.proj))
experts.temp.rcp26.05 = rep(NA,length(t.proj)); experts.temp.rcp26.50 = rep(NA,length(t.proj)); experts.temp.rcp26.95 = rep(NA,length(t.proj))
experts.ocheat.rcp26.05 = rep(NA,length(t.proj)); experts.ocheat.rcp26.50 = rep(NA,length(t.proj)); experts.ocheat.rcp26.95 = rep(NA,length(t.proj))

# standard
standard.slr.rcp26.05 = rep(NA,length(t.proj)); standard.slr.rcp26.50 = rep(NA,length(t.proj)); standard.slr.rcp26.95 = rep(NA,length(t.proj))
standard.ais.rcp26.05 = rep(NA,length(t.proj)); standard.ais.rcp26.50 = rep(NA,length(t.proj)); standard.ais.rcp26.95 = rep(NA,length(t.proj))
standard.gis.rcp26.05 = rep(NA,length(t.proj)); standard.gis.rcp26.50 = rep(NA,length(t.proj)); standard.gis.rcp26.95 = rep(NA,length(t.proj))
standard.gsic.rcp26.05 = rep(NA,length(t.proj)); standard.gsic.rcp26.50 = rep(NA,length(t.proj)); standard.gsic.rcp26.95 = rep(NA,length(t.proj))
standard.te.rcp26.05 = rep(NA,length(t.proj)); standard.te.rcp26.50 = rep(NA,length(t.proj)); standard.te.rcp26.95 = rep(NA,length(t.proj))
standard.temp.rcp26.05 = rep(NA,length(t.proj)); standard.temp.rcp26.50 = rep(NA,length(t.proj)); standard.temp.rcp26.95 = rep(NA,length(t.proj))
standard.ocheat.rcp26.05 = rep(NA,length(t.proj)); standard.ocheat.rcp26.50 = rep(NA,length(t.proj)); standard.ocheat.rcp26.95 = rep(NA,length(t.proj))

# complete
complete.slr.rcp26.05 = rep(NA,length(t.proj)); complete.slr.rcp26.50 = rep(NA,length(t.proj)); complete.slr.rcp26.95 = rep(NA,length(t.proj))
complete.ais.rcp26.05 = rep(NA,length(t.proj)); complete.ais.rcp26.50 = rep(NA,length(t.proj)); complete.ais.rcp26.95 = rep(NA,length(t.proj))
complete.gis.rcp26.05 = rep(NA,length(t.proj)); complete.gis.rcp26.50 = rep(NA,length(t.proj)); complete.gis.rcp26.95 = rep(NA,length(t.proj))
complete.gsic.rcp26.05 = rep(NA,length(t.proj)); complete.gsic.rcp26.50 = rep(NA,length(t.proj)); complete.gsic.rcp26.95 = rep(NA,length(t.proj))
complete.te.rcp26.05 = rep(NA,length(t.proj)); complete.te.rcp26.50 = rep(NA,length(t.proj)); complete.te.rcp26.95 = rep(NA,length(t.proj))
complete.temp.rcp26.05 = rep(NA,length(t.proj)); complete.temp.rcp26.50 = rep(NA,length(t.proj)); complete.temp.rcp26.95 = rep(NA,length(t.proj))
complete.ocheat.rcp26.05 = rep(NA,length(t.proj)); complete.ocheat.rcp26.50 = rep(NA,length(t.proj)); complete.ocheat.rcp26.95 = rep(NA,length(t.proj))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.proj)){
  # priors
  c(priors.slr.rcp26.05[t], priors.slr.rcp26.50[t], priors.slr.rcp26.95[t]) := quantile(priors.slr.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.ais.rcp26.05[t], priors.ais.rcp26.50[t], priors.ais.rcp26.95[t]) := quantile(priors.ais.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.gis.rcp26.05[t], priors.gis.rcp26.50[t], priors.gis.rcp26.95[t]) := quantile(priors.gis.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.gsic.rcp26.05[t], priors.gsic.rcp26.50[t], priors.gsic.rcp26.95[t]) := quantile(priors.gsic.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.te.rcp26.05[t], priors.te.rcp26.50[t], priors.te.rcp26.95[t]) := quantile(priors.te.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.temp.rcp26.05[t], priors.temp.rcp26.50[t], priors.temp.rcp26.95[t]) := quantile(priors.temp.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(priors.ocheat.rcp26.05[t], priors.ocheat.rcp26.50[t], priors.ocheat.rcp26.95[t]) := quantile(priors.ocheat.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # experts
  c(experts.slr.rcp26.05[t], experts.slr.rcp26.50[t], experts.slr.rcp26.95[t]) := quantile(experts.slr.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.ais.rcp26.05[t], experts.ais.rcp26.50[t], experts.ais.rcp26.95[t]) := quantile(experts.ais.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.gis.rcp26.05[t], experts.gis.rcp26.50[t], experts.gis.rcp26.95[t]) := quantile(experts.gis.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.gsic.rcp26.05[t], experts.gsic.rcp26.50[t], experts.gsic.rcp26.95[t]) := quantile(experts.gsic.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.te.rcp26.05[t], experts.te.rcp26.50[t], experts.te.rcp26.95[t]) := quantile(experts.te.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.temp.rcp26.05[t], experts.temp.rcp26.50[t], experts.temp.rcp26.95[t]) := quantile(experts.temp.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(experts.ocheat.rcp26.05[t], experts.ocheat.rcp26.50[t], experts.ocheat.rcp26.95[t]) := quantile(experts.ocheat.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # standard
  c(standard.slr.rcp26.05[t], standard.slr.rcp26.50[t], standard.slr.rcp26.95[t]) := quantile(standard.slr.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.ais.rcp26.05[t], standard.ais.rcp26.50[t], standard.ais.rcp26.95[t]) := quantile(standard.ais.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.gis.rcp26.05[t], standard.gis.rcp26.50[t], standard.gis.rcp26.95[t]) := quantile(standard.gis.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.gsic.rcp26.05[t], standard.gsic.rcp26.50[t], standard.gsic.rcp26.95[t]) := quantile(standard.gsic.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.te.rcp26.05[t], standard.te.rcp26.50[t], standard.te.rcp26.95[t]) := quantile(standard.te.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.temp.rcp26.05[t], standard.temp.rcp26.50[t], standard.temp.rcp26.95[t]) := quantile(standard.temp.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(standard.ocheat.rcp26.05[t], standard.ocheat.rcp26.50[t], standard.ocheat.rcp26.95[t]) := quantile(standard.ocheat.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  
  # complete
  c(complete.slr.rcp26.05[t], complete.slr.rcp26.50[t], complete.slr.rcp26.95[t]) := quantile(complete.slr.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.ais.rcp26.05[t], complete.ais.rcp26.50[t], complete.ais.rcp26.95[t]) := quantile(complete.ais.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.gis.rcp26.05[t], complete.gis.rcp26.50[t], complete.gis.rcp26.95[t]) := quantile(complete.gis.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.gsic.rcp26.05[t], complete.gsic.rcp26.50[t], complete.gsic.rcp26.95[t]) := quantile(complete.gsic.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.te.rcp26.05[t], complete.te.rcp26.50[t], complete.te.rcp26.95[t]) := quantile(complete.te.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.temp.rcp26.05[t], complete.temp.rcp26.50[t], complete.temp.rcp26.95[t]) := quantile(complete.temp.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(complete.ocheat.rcp26.05[t], complete.ocheat.rcp26.50[t], complete.ocheat.rcp26.95[t]) := quantile(complete.ocheat.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
}

iproj = which(t.proj==2000):which(t.proj==2100)
i2100 = which(t.proj==2100)

# compute SLR densities
# GIS
priors.gis.2100 <- priors.gis.rcp26[nrow(priors.gis.rcp26),]
experts.gis.2100 <- experts.gis.rcp26[nrow(experts.gis.rcp26),]
standard.gis.2100 <- standard.gis.rcp26[nrow(standard.gis.rcp26),]
complete.gis.2100 <- complete.gis.rcp26[nrow(complete.gis.rcp26),]

priors.density.gis.2100 <- density(priors.gis.2100)
experts.density.gis.2100 <- density(experts.gis.2100)
standard.density.gis.2100 <- density(standard.gis.2100)
complete.density.gis.2100 <- density(complete.gis.2100)

# AIS
priors.ais.2100 <- priors.ais.rcp26[nrow(priors.ais.rcp26),]
experts.ais.2100 <- experts.ais.rcp26[nrow(experts.ais.rcp26),]
standard.ais.2100 <- standard.ais.rcp26[nrow(standard.ais.rcp26),]
complete.ais.2100 <- complete.ais.rcp26[nrow(complete.ais.rcp26),]

priors.density.ais.2100 <- density(priors.ais.2100)
experts.density.ais.2100 <- density(experts.ais.2100)
standard.density.ais.2100 <- density(standard.ais.2100)
complete.density.ais.2100 <- density(complete.ais.2100)

# TOTAL SLR
priors.slr.2100 <- priors.slr.rcp26[nrow(priors.slr.rcp26),]
experts.slr.2100 <- experts.slr.rcp26[nrow(experts.slr.rcp26),]
standard.slr.2100 <- standard.slr.rcp26[nrow(standard.slr.rcp26),]
complete.slr.2100 <- complete.slr.rcp26[nrow(complete.slr.rcp26),]

priors.density.2100 <- density(priors.slr.2100)
experts.density.2100 <- density(experts.slr.2100)
standard.density.2100 <- density(standard.slr.2100)
complete.density.2100 <- density(complete.slr.2100)

# TEMPERATURE
priors.temp.2100 <- priors.temp.rcp26[nrow(priors.temp.rcp26),]
experts.temp.2100 <- experts.temp.rcp26[nrow(experts.temp.rcp26),]
standard.temp.2100 <- standard.temp.rcp26[nrow(standard.temp.rcp26),]
complete.temp.2100 <- complete.temp.rcp26[nrow(complete.temp.rcp26),]

priors.density.temp.2100 <- density(priors.temp.2100)
experts.density.temp.2100 <- density(experts.temp.2100)
standard.density.temp.2100 <- density(standard.temp.2100)
complete.density.temp.2100 <- density(complete.temp.2100)

# TE
priors.te.2100 <- priors.te.rcp26[nrow(priors.te.rcp26),]
experts.te.2100 <- experts.te.rcp26[nrow(experts.te.rcp26),]
standard.te.2100 <- standard.te.rcp26[nrow(standard.te.rcp26),]
complete.te.2100 <- complete.te.rcp26[nrow(complete.te.rcp26),]

# GSIC
priors.gsic.2100 <- priors.gsic.rcp26[nrow(priors.gsic.rcp26),]
experts.gsic.2100 <- experts.gsic.rcp26[nrow(experts.gsic.rcp26),]
standard.gsic.2100 <- standard.gsic.rcp26[nrow(standard.gsic.rcp26),]
complete.gsic.2100 <- complete.gsic.rcp26[nrow(complete.gsic.rcp26),]


if(FALSE){

          
          print('==============================================================')
          print('min/5%/50%/95%/max of 2100 sea level relative to 1986-2005:')
          print(paste('RCP2.6: ',quantile(standard.slr.rcp26[251,],c(0,.05,.50,.95,1))))
          # print(paste('RCP4.5: ',quantile(slr.rcp45[251,],c(0,.05,.50,.95,1))))
          # print(paste('RCP8.5: ',quantile(slr.rcp85[251,],c(0,.05,.50,.95,1))))
          print('==============================================================')
          
          i2050 <- which(t.proj==2050)
          print('==============================================================')
          print('min/5%/50%/95%/max of 2050 sea level relative to 1986-2005:')
          print(paste('RCP2.6: ',quantile(standard.slr.rcp26[i2050,],c(0,.05,.50,.95,1))))
          # print(paste('RCP4.5: ',quantile(slr.rcp45[i2050,],c(0,.05,.50,.95,1))))
          # print(paste('RCP8.5: ',quantile(slr.rcp85[i2050,],c(0,.05,.50,.95,1))))
          print('==============================================================')
          
          ##==============================================================================
          
          ##
          ## 5-95% CI and median of 2100 SLR and all contributions under the RCP scenarios
          ##
          
          # slr.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
          #                               c(standard.slr.rcp26.50[i2100],slr.rcp45.50[i2100],slr.rcp85.50[i2100]),
          #                               c(standard.slr.rcp26.05[i2100],slr.rcp45.05[i2100],slr.rcp85.05[i2100]),
          #                               c(standard.slr.rcp26.95[i2100],slr.rcp45.95[i2100],slr.rcp85.95[i2100])
          # ))
          
          # gsic.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
          #                                c(standard.gsic.rcp26.50[i2100],gsic.rcp45.50[i2100],gsic.rcp85.50[i2100]),
          #                                c(standard.gsic.rcp26.05[i2100],gsic.rcp45.05[i2100],gsic.rcp85.05[i2100]),
          #                                c(standard.gsic.rcp26.95[i2100],gsic.rcp45.95[i2100],gsic.rcp85.95[i2100])
          # ))
          # gis.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
          #                               c(standard.gis.rcp26.50[i2100],gis.rcp45.50[i2100],gis.rcp85.50[i2100]),
          #                               c(standard.gis.rcp26.05[i2100],gis.rcp45.05[i2100],gis.rcp85.05[i2100]),
          #                               c(standard.gis.rcp26.95[i2100],gis.rcp45.95[i2100],gis.rcp85.95[i2100])
          # ))
          # te.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
          #                              c(standard.te.rcp26.50[i2100],te.rcp45.50[i2100],te.rcp85.50[i2100]),
          #                              c(standard.te.rcp26.05[i2100],te.rcp45.05[i2100],te.rcp85.05[i2100]),
          #                              c(standard.te.rcp26.95[i2100],te.rcp45.95[i2100],te.rcp85.95[i2100])
          # ))
          # ais.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
          #                               c(standard.ais.rcp26.50[i2100],ais.rcp45.50[i2100],ais.rcp85.50[i2100]),
          #                               c(standard.ais.rcp26.05[i2100],ais.rcp45.05[i2100],ais.rcp85.05[i2100]),
          #                               c(standard.ais.rcp26.95[i2100],ais.rcp45.95[i2100],ais.rcp85.95[i2100])
          # ))
          
          slr.ci90 = data.frame( cbind( c('RCP2.6'),
                                        c(standard.slr.rcp26.50[i2100]),
                                        c(standard.slr.rcp26.05[i2100]),
                                        c(standard.slr.rcp26.95[i2100])
          ))
          
          gsic.ci90 = data.frame( cbind( c('RCP2.6'),
                                         c(standard.gsic.rcp26.50[i2100]),
                                         c(standard.gsic.rcp26.05[i2100]),
                                         c(standard.gsic.rcp26.95[i2100])
          ))
          gis.ci90 = data.frame( cbind( c('RCP2.6'),
                                        c(standard.gis.rcp26.50[i2100]),
                                        c(standard.gis.rcp26.05[i2100]),
                                        c(standard.gis.rcp26.95[i2100])
          ))
          te.ci90 = data.frame( cbind( c('RCP2.6'),
                                       c(standard.te.rcp26.50[i2100]),
                                       c(standard.te.rcp26.05[i2100]),
                                       c(standard.te.rcp26.95[i2100])
          ))
          ais.ci90 = data.frame( cbind( c('RCP2.6'),
                                        c(standard.ais.rcp26.50[i2100]),
                                        c(standard.ais.rcp26.05[i2100]),
                                        c(standard.ais.rcp26.95[i2100])
          ))
          ## get into format for Latex table
          row.gsic = paste('Glaciers and small ice caps &',
                           1000*signif(standard.gsic.rcp26.50[i2100],4),' (',1000*signif(standard.gsic.rcp26.05[i2100],4),'-',1000*signif(standard.gsic.rcp26.95[i2100],4),') &',
                           # 1000*signif(gsic.rcp45.50[i2100],4),' (',1000*signif(gsic.rcp45.05[i2100],4),'-',1000*signif(gsic.rcp45.95[i2100],4),') &',
                           # 1000*signif(gsic.rcp85.50[i2100],4),' (',1000*signif(gsic.rcp85.05[i2100],4),'-',1000*signif(gsic.rcp85.95[i2100],4),') \\',
                           sep='')
          row.te = paste('Thermal expansion &',
                         1000*signif(standard.te.rcp26.50[i2100],4),' (',1000*signif(standard.te.rcp26.05[i2100],4),'-',1000*signif(standard.te.rcp26.95[i2100],4),') &',
                         # 1000*signif(te.rcp45.50[i2100],4),' (',1000*signif(te.rcp45.05[i2100],4),'-',1000*signif(te.rcp45.95[i2100],4),') &',
                         # 1000*signif(te.rcp85.50[i2100],4),' (',1000*signif(te.rcp85.05[i2100],4),'-',1000*signif(te.rcp85.95[i2100],4),') \\',
                         sep='')
          row.gis = paste('Greenland Ice Sheet &',
                          1000*signif(standard.gis.rcp26.50[i2100],4),' (',1000*signif(standard.gis.rcp26.05[i2100],4),'-',1000*signif(standard.gis.rcp26.95[i2100],4),') &',
                          # 1000*signif(gis.rcp45.50[i2100],4),' (',1000*signif(gis.rcp45.05[i2100],4),'-',1000*signif(gis.rcp45.95[i2100],4),') &',
                          # 1000*signif(gis.rcp85.50[i2100],4),' (',1000*signif(gis.rcp85.05[i2100],4),'-',1000*signif(gis.rcp85.95[i2100],4),') \\',
                          sep='')
          row.ais = paste('Antarctic Ice Sheet &',
                          1000*signif(standard.ais.rcp26.50[i2100],4),' (',1000*signif(standard.ais.rcp26.05[i2100],4),'-',1000*signif(standard.ais.rcp26.95[i2100],4),') &',
                          # 1000*signif(ais.rcp45.50[i2100],4),' (',1000*signif(ais.rcp45.05[i2100],4),'-',1000*signif(ais.rcp45.95[i2100],4),') &',
                          # 1000*signif(ais.rcp85.50[i2100],4),' (',1000*signif(ais.rcp85.05[i2100],4),'-',1000*signif(ais.rcp85.95[i2100],4),') \\',
                          sep='')
          row.slr = paste('Total sea level &',
                          1000*signif(standard.slr.rcp26.50[i2100],4),' (',1000*signif(standard.slr.rcp26.05[i2100],4),'-',1000*signif(standard.slr.rcp26.95[i2100],4),') &',
                          # 1000*signif(slr.rcp45.50[i2100],4),' (',1000*signif(slr.rcp45.05[i2100],4),'-',1000*signif(slr.rcp45.95[i2100],4),') &',
                          # 1000*signif(slr.rcp85.50[i2100],4),' (',1000*signif(slr.rcp85.05[i2100],4),'-',1000*signif(slr.rcp85.95[i2100],4),') \\',
                          sep='')
          
          # pdf(paste(plotdir,'projections_SLR_total.pdf',sep=''),width=3.5,height=2.45,colormodel='cmyk')
          # par(mfrow=c(1,1))
          # # RCP85
          # par(mai=c(.65,.65,.20,.2))
          # plot(t.proj[iproj],slr.rcp85.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann=FALSE,
          #      xlim=c(2000,2100), ylim=c(0,2), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
          # axis(1, seq(2000,2100,by=20)); axis(2, seq(0,2,by=.25), lab=c('0','','0.5','','1','','1.5','','2'));
          # mtext(side=2, text='Global mean sea level [m]', line=2.2, cex=1);
          # mtext(side=1, text='Year', line=2.2, cex=1);
          # polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp85.95[iproj],rev(slr.rcp85.05[iproj])),
          #         col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
          # # + UNIFORM RCP45
          # lines(t.proj[iproj],slr.rcp45.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
          # polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp45.95[iproj],rev(slr.rcp45.05[iproj])),
          #         col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
          # # + UNIFORM RCP26
          # lines(t.proj[iproj],standard.slr.rcp26.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
          # polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(standard.slr.rcp26.95[iproj],rev(standard.slr.rcp26.05[iproj])),
          #         col=rgb(col26[1],col26[2],col26[3],.5), border=NA);
          # # + legend
          # legend(t.proj[iproj[1]],2,c("5-95% range,","RCP2.6","RCP4.5","RCP8.5"),
          #        lty=c(NA,1,1,1), lwd=3, bty='n', cex=1,
          #        col=c(NA , rgb(col26[1],col26[2],col26[3]) ,
          #              rgb(col45[1],col45[2],col45[3]) , rgb(col85[1],col85[2],col85[3])))
          # 
          # dev.off()
}

if(FALSE){
# pdf(paste(plotdir,'projections_SLR_total_with_noise_and_normalization.pdf',sep=''),width=3.5,height=2.45,colormodel='cmyk')
png(paste(plotdir,'standard_projections.png',sep=''), width=866, height=516, units ="px")

par(mfrow=c(1,1))
# RCP85
par(mai=c(.65,.65,.20,.2))
plot(t.proj[iproj],standard.slr.rcp26.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2, #ann=FALSE,
     xlim=c(2000,2100), ylim=c(0,3), 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i'
     );
axis(1, seq(2000,2100,by=20)); 
axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
mtext(side=2, text='Global mean sea level [m]', line=2.2, cex=1);
mtext(side=1, text='Year', line=2.2, cex=1);
# polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp85.95[iproj],rev(slr.rcp85.05[iproj])),
#         col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + UNIFORM RCP45
# lines(t.proj[iproj],slr.rcp45.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
# polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp45.95[iproj],rev(slr.rcp45.05[iproj])),
#         col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
# + UNIFORM RCP26
# lines(t.proj[iproj],standard.slr.rcp26.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(standard.slr.rcp26.95[iproj],rev(standard.slr.rcp26.05[iproj])),
        col=rgb(col26[1],col26[2],col26[3],.5), border=NA);
# + legend
legend(
       # t.proj[iproj[1]],2.5,c("5-95% range","median"),
       "topleft", legend=c("5-95% range","median"),
       lty=c(1,1), lwd=c(9,3), bty='n', cex=1,
       # fill=c(rgb(col26[1],col26[2],col26[3],.5),NA),
       # pch=c(22,NA),
       col=c(rgb(col26[1],col26[2],col26[3],.5), rgb(col26[1],col26[2],col26[3])))

dev.off()
}
##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## Supplemental Figures
## Combined hindcasts and projections for each of Greenland, Antarctica, and Total SLR

## mycol colors, by row:
# 1 black
# 2 dark blue
# 3 aqua
# 4 pink
# 5 light orange
# 6 purple
# 7 blue
# 8 light purple
# 9 light blue
# 10 very light blue
# 11 dark red
# 12 brown
# 13 dark orange
# 14 neon green
# 15 neon yellow

col_priors <- 2
col_experts <- 2
col_standard <- 2
col_complete <- 2
col_obs <- 11
col_density <- 1
opacity <- 0.5

# dev.off()
if(TRUE){
# >>> GIS <<< ###################################################################################################
layout(cbind(c(1,4,7,10,13),c(2,5,8,11,13), c(3,6,9,12,13)))
par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) 
# omi = c(0,0,0,2) #for right side space

# priors #######################################

# plot hindcasts 
plot(mod.time, priors.gis.50, type='l', # median
     col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
mtext(side=1, text='Year', line=2.3, cex=.9);
# mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' a) Prior Hindcasts')), line=.25, cex=.9, adj=0);

polygon(c(mod.time,rev(mod.time)), c(priors.gis.95,rev(priors.gis.05)), # polygon
        col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3],opacity), border=NA);

# observations
# lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
points(obs.gis.time, obs.gis.norm, pch=20, 
       col=rgb(mycol[col_obs,1],mycol[col_obs,2],mycol[col_obs,3]));
# lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
abline(h=0,col="black",lty=2)

# plot priors projections
plot(t.proj[iproj],priors.gis.rcp26.50[iproj],type='l', ann=FALSE, # median
     col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3]),lwd=2,
     xlim=c(2000,2100), ylim=c(-1,3), 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i'
);
axis(1, seq(2000,2100,by=20)); 
# axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
axis(2, seq(-1,3,by=1), lab=c('-1','','1','','3'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Year', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' b) Prior Projections')), line=.25, cex=.9, adj=0);

polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(priors.gis.rcp26.95[iproj],rev(priors.gis.rcp26.05[iproj])),
        col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3],opacity), border=NA);
abline(h=0,col="black",lty=2)

# plot densities
plot(x=priors.density.gis.2100$y,
     y=priors.density.gis.2100$x,
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     ylim=c(-1,3),
     xlim=c(0,1.45),
     type="n",
     ann=FALSE,
     # xlab='Probability density',
     # ylab='[m SLE]', 
     lwd=2, col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3]), 
     xaxt='n')
axis(2, seq(-1,3,by=1), lab=c('-1','','1','','3'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Probability density', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' c) Prior SLR at 2100')), line=.25, cex=.9, adj=0);
polygon(x=c(rep(0,length(priors.density.gis.2100$x)),rev(priors.density.gis.2100$y)),
        y=c(priors.density.gis.2100$x,rev(priors.density.gis.2100$x)),
        col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
        # border=NA
)
abline(h=0,col="black",lty=2)
    
# experts #######################################
    
# plot hindcasts 
plot(mod.time, experts.gis.50, type='l', # median
     col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
mtext(side=1, text='Year', line=2.3, cex=.9);
# mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' d) Expert Hindcasts')), line=.25, cex=.9, adj=0);

polygon(c(mod.time,rev(mod.time)), c(experts.gis.95,rev(experts.gis.05)), # polygon
        col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3],opacity), border=NA);

# observations
# lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
points(obs.gis.time, obs.gis.norm, pch=20, 
       col=rgb(mycol[col_obs,1],mycol[col_obs,2],mycol[col_obs,3]));
# lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
abline(h=0,col="black",lty=2)

# plot experts projections
plot(t.proj[iproj],experts.gis.rcp26.50[iproj],type='l', ann=FALSE, # median
     col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3]),lwd=2,
     xlim=c(2000,2100), ylim=c(-1,3), 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i'
);
axis(1, seq(2000,2100,by=20)); 
# axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
axis(2, seq(-1,3,by=1), lab=c('-1','','1','','3'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Year', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' e) Expert Projections')), line=.25, cex=.9, adj=0);

polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(experts.gis.rcp26.95[iproj],rev(experts.gis.rcp26.05[iproj])),
        col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3],opacity), border=NA);
abline(h=0,col="black",lty=2)

# plot densities
plot(x=experts.density.gis.2100$y,
     y=experts.density.gis.2100$x,
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     ylim=c(-1,3),
     xlim=c(0,1.45),
     type="n",
     ann=FALSE,
     # xlab='Probability density',
     # ylab='[m SLE]', 
     lwd=2, col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3]), 
     xaxt='n')
axis(2, seq(-1,3,by=1), lab=c('-1','','1','','3'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Probability density', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' f) Expert SLR at 2100')), line=.25, cex=.9, adj=0);
polygon(x=c(rep(0,length(experts.density.gis.2100$x)),rev(experts.density.gis.2100$y)),
        y=c(experts.density.gis.2100$x,rev(experts.density.gis.2100$x)),
        col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
        # border=NA
)
abline(h=0,col="black",lty=2)

# standard #######################################

# plot hindcasts 
plot(mod.time, standard.gis.50, type='l', # median
     col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
mtext(side=1, text='Year', line=2.3, cex=.9);
# mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' g) Standard Hindcasts')), line=.25, cex=.9, adj=0);

polygon(c(mod.time,rev(mod.time)), c(standard.gis.95,rev(standard.gis.05)), # polygon
        col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3],opacity), border=NA);

# observations
# lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
points(obs.gis.time, obs.gis.norm, pch=20, 
       col=rgb(mycol[col_obs,1],mycol[col_obs,2],mycol[col_obs,3]));
# lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
abline(h=0,col="black",lty=2)

# plot standard projections
plot(t.proj[iproj],standard.gis.rcp26.50[iproj],type='l', ann=FALSE, # median
     col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]),lwd=2,
     xlim=c(2000,2100), ylim=c(-1,3), 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i'
);
axis(1, seq(2000,2100,by=20)); 
# axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
axis(2, seq(-1,3,by=1), lab=c('-1','','1','','3'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Year', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' h) Standard Projections')), line=.25, cex=.9, adj=0);

polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(standard.gis.rcp26.95[iproj],rev(standard.gis.rcp26.05[iproj])),
        col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3],opacity), border=NA);
abline(h=0,col="black",lty=2)

# plot densities
plot(x=standard.density.gis.2100$y,
     y=standard.density.gis.2100$x,
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     ylim=c(-1,3),
     xlim=c(0,1.45),
     type="n",
     ann=FALSE,
     # xlab='Probability density',
     # ylab='[m SLE]', 
     lwd=2, col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]), 
     xaxt='n')
axis(2, seq(-1,3,by=1), lab=c('-1','','1','','3'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Probability density', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' i) Standard SLR at 2100')), line=.25, cex=.9, adj=0);
polygon(x=c(rep(0,length(standard.density.gis.2100$x)),rev(standard.density.gis.2100$y)),
        y=c(standard.density.gis.2100$x,rev(standard.density.gis.2100$x)),
        col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
        # border=NA
)
abline(h=0,col="black",lty=2)

# complete #######################################

# plot hindcasts 
plot(mod.time, complete.gis.50, type='l', # median
     col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
mtext(side=1, text='Year', line=2.3, cex=.9);
# mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' j) Complete Hindcasts')), line=.25, cex=.9, adj=0);

polygon(c(mod.time,rev(mod.time)), c(complete.gis.95,rev(complete.gis.05)), # polygon
        col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);

# observations
# lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
points(obs.gis.time, obs.gis.norm, pch=20, 
       col=rgb(mycol[col_obs,1],mycol[col_obs,2],mycol[col_obs,3]));
# lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
abline(h=0,col="black",lty=2)

# plot complete projections
plot(t.proj[iproj],complete.gis.rcp26.50[iproj],type='l', ann=FALSE, # median
     col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]),lwd=2,
     xlim=c(2000,2100), ylim=c(-1,3), 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i'
);
axis(1, seq(2000,2100,by=20)); 
# axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
axis(2, seq(-1,3,by=1), lab=c('-1','','1','','3'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Year', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' k) Complete Projections')), line=.25, cex=.9, adj=0);

polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.gis.rcp26.95[iproj],rev(complete.gis.rcp26.05[iproj])),
        col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
abline(h=0,col="black",lty=2)

# plot densities
plot(x=complete.density.gis.2100$y,
     y=complete.density.gis.2100$x,
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     ylim=c(-1,3),
     xlim=c(0,1.45),
     type="n",
     ann=FALSE,
     # xlab='Probability density',
     # ylab='[m SLE]', 
     lwd=2, col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), 
     xaxt='n')
axis(2, seq(-1,3,by=1), lab=c('-1','','1','','3'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Probability density', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' l) Complete SLR at 2100')), line=.25, cex=.9, adj=0);
polygon(x=c(rep(0,length(complete.density.gis.2100$x)),rev(complete.density.gis.2100$y)),
        y=c(complete.density.gis.2100$x,rev(complete.density.gis.2100$x)),
        col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
        # border=NA
)
abline(h=0,col="black",lty=2)

# LEGEND
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("median, model" ,
                  "5-95% range, model",
                  "observations"#,
                  # "Probability density"
                  # "2-sigma range, observations"
                  ),
       # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
       lwd=c(2,8,NA#,
             #8
             ), bty='n', cex=1.2,
       pch=c(NA,NA,20,NA),
       col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
             rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
             rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
             # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])
             # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
       ), 
       horiz = FALSE)
}
##==============================================================================
##==============================================================================

if(TRUE){
  # >>> Antarctica <<< ###################################################################################################
  layout(cbind(c(1,4,7,10,13),c(2,5,8,11,13), c(3,6,9,12,13)))
  par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) 
  # omi = c(0,0,0,2) #for right side space
  
  # priors #######################################
  
  # plot hindcasts 
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], priors.ais.paleo.50[ipaleo], type='l', 
       col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' a) Prior Hindcasts')), line=.25, cex=.9, adj=0);
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(priors.ais.paleo.95[ipaleo],rev(priors.ais.paleo.05[ipaleo])), 
          col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  abline(h=0,col="black",lty=2)
  
  # plot priors projections
  plot(t.proj[iproj],priors.ais.rcp26.50[iproj],type='l', ann=FALSE, # median
       col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3]),lwd=2,
       xlim=c(2000,2100), ylim=c(-2,2), 
       xaxt='n', yaxt='n', xaxs='i', yaxs='i'
  );
  axis(1, seq(2000,2100,by=20)); 
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(-2,2,by=1), lab=c('-2','','0','','2'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b) Prior Projections')), line=.25, cex=.9, adj=0);
  
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(priors.ais.rcp26.95[iproj],rev(priors.ais.rcp26.05[iproj])),
          col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3],opacity), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot densities
  plot(x=priors.density.ais.2100$y,
       y=priors.density.ais.2100$x,
       xaxt='n', yaxt='n', xaxs='i', yaxs='i',
       ylim=c(-2,2),
       xlim=c(0,4.35),
       type="n",
       ann=FALSE,
       # xlab='Probability density',
       # ylab='[m SLE]', 
       lwd=2, col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3]), 
       xaxt='n')
  axis(2, seq(-2,2,by=1), lab=c('-2','','0','','2'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Probability density', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' c) Prior SLR at 2100')), line=.25, cex=.9, adj=0);
  polygon(x=c(rep(0,length(priors.density.ais.2100$x)),rev(priors.density.ais.2100$y)),
          y=c(priors.density.ais.2100$x,rev(priors.density.ais.2100$x)),
          col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
          # border=NA
  )
  abline(h=0,col="black",lty=2)
  
  # experts #######################################
  
  # plot hindcasts 
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], experts.ais.paleo.50[ipaleo], type='l', 
       col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d) Expert Hindcasts')), line=.25, cex=.9, adj=0);
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(experts.ais.paleo.95[ipaleo],rev(experts.ais.paleo.05[ipaleo])), 
          col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  abline(h=0,col="black",lty=2)
  
  # plot experts projections
  plot(t.proj[iproj],experts.ais.rcp26.50[iproj],type='l', ann=FALSE, # median
       col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3]),lwd=2,
       xlim=c(2000,2100), ylim=c(-2,2), 
       xaxt='n', yaxt='n', xaxs='i', yaxs='i'
  );
  axis(1, seq(2000,2100,by=20)); 
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(-2,2,by=1), lab=c('-2','','0','','2'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e) Expert Projections')), line=.25, cex=.9, adj=0);
  
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(experts.ais.rcp26.95[iproj],rev(experts.ais.rcp26.05[iproj])),
          col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3],opacity), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot densities
  plot(x=experts.density.ais.2100$y,
       y=experts.density.ais.2100$x,
       xaxt='n', yaxt='n', xaxs='i', yaxs='i',
       ylim=c(-2,2),
       xlim=c(0,4.35),
       type="n",
       ann=FALSE,
       # xlab='Probability density',
       # ylab='[m SLE]', 
       lwd=2, col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3]), 
       xaxt='n')
  axis(2, seq(-2,2,by=1), lab=c('-2','','0','','2'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Probability density', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' f) Expert SLR at 2100')), line=.25, cex=.9, adj=0);
  polygon(x=c(rep(0,length(experts.density.ais.2100$x)),rev(experts.density.ais.2100$y)),
          y=c(experts.density.ais.2100$x,rev(experts.density.ais.2100$x)),
          col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
          # border=NA
  )
  abline(h=0,col="black",lty=2)
  
  # standard #######################################
  
  # plot hindcasts 
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], standard.ais.paleo.50[ipaleo], type='l', 
       col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' g) Standard Hindcasts')), line=.25, cex=.9, adj=0);
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(standard.ais.paleo.95[ipaleo],rev(standard.ais.paleo.05[ipaleo])), 
          col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  abline(h=0,col="black",lty=2)
  
  # plot standard projections
  plot(t.proj[iproj],standard.ais.rcp26.50[iproj],type='l', ann=FALSE, # median
       col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]),lwd=2,
       xlim=c(2000,2100), ylim=c(-2,2), 
       xaxt='n', yaxt='n', xaxs='i', yaxs='i'
  );
  axis(1, seq(2000,2100,by=20)); 
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(-2,2,by=1), lab=c('-2','','0','','2'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' h) Standard Projections')), line=.25, cex=.9, adj=0);
  
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(standard.ais.rcp26.95[iproj],rev(standard.ais.rcp26.05[iproj])),
          col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3],opacity), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot densities
  plot(x=standard.density.ais.2100$y,
       y=standard.density.ais.2100$x,
       xaxt='n', yaxt='n', xaxs='i', yaxs='i',
       ylim=c(-2,2),
       xlim=c(0,4.35),
       type="n",
       ann=FALSE,
       # xlab='Probability density',
       # ylab='[m SLE]', 
       lwd=2, col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]), 
       xaxt='n')
  axis(2, seq(-2,2,by=1), lab=c('-2','','0','','2'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Probability density', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' i) Standard SLR at 2100')), line=.25, cex=.9, adj=0);
  polygon(x=c(rep(0,length(standard.density.ais.2100$x)),rev(standard.density.ais.2100$y)),
          y=c(standard.density.ais.2100$x,rev(standard.density.ais.2100$x)),
          col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
          # border=NA
  )
  abline(h=0,col="black",lty=2)
  
  # complete #######################################
  
  # plot hindcasts 
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], complete.ais.paleo.50[ipaleo], type='l', 
       col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' j) Complete Hindcasts')), line=.25, cex=.9, adj=0);
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(complete.ais.paleo.95[ipaleo],rev(complete.ais.paleo.05[ipaleo])), 
          col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  abline(h=0,col="black",lty=2)
  
  # plot complete projections
  plot(t.proj[iproj],complete.ais.rcp26.50[iproj],type='l', ann=FALSE, # median
       col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]),lwd=2,
       xlim=c(2000,2100), ylim=c(-2,2), 
       xaxt='n', yaxt='n', xaxs='i', yaxs='i'
  );
  axis(1, seq(2000,2100,by=20)); 
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(-2,2,by=1), lab=c('-2','','0','','2'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' k) Complete Projections')), line=.25, cex=.9, adj=0);
  
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.ais.rcp26.95[iproj],rev(complete.ais.rcp26.05[iproj])),
          col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot densities
  plot(x=complete.density.ais.2100$y,
       y=complete.density.ais.2100$x,
       xaxt='n', yaxt='n', xaxs='i', yaxs='i',
       ylim=c(-2,2),
       xlim=c(0,4.35),
       type="n",
       ann=FALSE,
       # xlab='Probability density',
       # ylab='[m SLE]', 
       lwd=2, col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), 
       xaxt='n')
  axis(2, seq(-2,2,by=1), lab=c('-2','','0','','2'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Probability density', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' l) Complete SLR at 2100')), line=.25, cex=.9, adj=0);
  polygon(x=c(rep(0,length(complete.density.ais.2100$x)),rev(complete.density.ais.2100$y)),
          y=c(complete.density.ais.2100$x,rev(complete.density.ais.2100$x)),
          col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
          # border=NA
  )
  abline(h=0,col="black",lty=2)
  
  
  
  # LEGEND
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top",inset = 0,
         legend = c("median, model" ,
                    "5-95% range, model",
                    "observations"#,
                    # "Probability density"
                    # "2-sigma range, observations"
         ),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,NA#,
               #8
         ), bty='n', cex=1.2,
         pch=c(NA,NA,20,NA),
         col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
               rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
               rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
         ), 
         horiz = FALSE)
}
##==============================================================================
##==============================================================================

if(TRUE){
  # >>> Total Sea Level <<< ###################################################################################################
  n.sig = 2         # how many sigma to plot around the obs?
  layout(cbind(c(1,4,7,10,13),c(2,5,8,11,13), c(3,6,9,12,13)))
  par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) 
  # omi = c(0,0,0,2) #for right side space
  
  # priors #######################################
  
  # plot hindcasts 
  plot(mod.time, priors.slr.50, type='l', 
       col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' a) Priors Hindcasts')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(priors.slr.95,rev(priors.slr.05)), 
          col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot priors projections
  plot(t.proj[iproj],priors.slr.rcp26.50[iproj],type='l', ann=FALSE, # median
       col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3]),lwd=2,
       xlim=c(2000,2100), ylim=c(0,4), 
       xaxt='n', yaxt='n', xaxs='i', yaxs='i'
  );
  axis(1, seq(2000,2100,by=20)); 
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(0,4,by=1), lab=c('0','','2','','4'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b) Prior Projections')), line=.25, cex=.9, adj=0);
  
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(priors.slr.rcp26.95[iproj],rev(priors.slr.rcp26.05[iproj])),
          col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3],opacity), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot densities
  plot(x=priors.density.2100$y,
       y=priors.density.2100$x,
       xaxt='n', yaxt='n', xaxs='i', yaxs='i',
       ylim=c(0,4),
       xlim=c(0,1.25),
       type="n",
       ann=FALSE,
       # xlab='Probability density',
       # ylab='[m SLE]', 
       lwd=2, col=rgb(mycol[col_priors,1],mycol[col_priors,2],mycol[col_priors,3]), 
       xaxt='n')
  axis(2, seq(0,4,by=1), lab=c('0','','2','','4'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Probability density', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' c) Prior SLR at 2100')), line=.25, cex=.9, adj=0);
  polygon(x=c(rep(0,length(priors.density.2100$x)),rev(priors.density.2100$y)),
          y=c(priors.density.2100$x,rev(priors.density.2100$x)),
          col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
          # border=NA
  )
  abline(h=0,col="black",lty=2)
  
  # experts #######################################
  
  # plot hindcasts 
  plot(mod.time, experts.slr.50, type='l', 
       col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d) Expert Hindcasts')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(experts.slr.95,rev(experts.slr.05)), 
          col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot experts projections
  plot(t.proj[iproj],experts.slr.rcp26.50[iproj],type='l', ann=FALSE, # median
       col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3]),lwd=2,
       xlim=c(2000,2100), ylim=c(0,4), 
       xaxt='n', yaxt='n', xaxs='i', yaxs='i'
  );
  axis(1, seq(2000,2100,by=20)); 
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(0,4,by=1), lab=c('0','','2','','4'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e) Expert Projections')), line=.25, cex=.9, adj=0);
  
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(experts.slr.rcp26.95[iproj],rev(experts.slr.rcp26.05[iproj])),
          col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3],opacity), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot densities
  plot(x=experts.density.2100$y,
       y=experts.density.2100$x,
       xaxt='n', yaxt='n', xaxs='i', yaxs='i',
       ylim=c(0,4),
       xlim=c(0,1.25),
       type="n",
       ann=FALSE,
       # xlab='Probability density',
       # ylab='[m SLE]', 
       lwd=2, col=rgb(mycol[col_experts,1],mycol[col_experts,2],mycol[col_experts,3]), 
       xaxt='n')
  axis(2, seq(0,4,by=1), lab=c('0','','2','','4'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Probability density', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' f) Expert SLR at 2100')), line=.25, cex=.9, adj=0);
  polygon(x=c(rep(0,length(experts.density.2100$x)),rev(experts.density.2100$y)),
          y=c(experts.density.2100$x,rev(experts.density.2100$x)),
          col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
          # border=NA
  )
  abline(h=0,col="black",lty=2)
  
  # standard #######################################
  
  # plot hindcasts 
  plot(mod.time, standard.slr.50, type='l', 
       col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' g) Standard Hindcasts')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(standard.slr.95,rev(standard.slr.05)), 
          col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot standard projections
  plot(t.proj[iproj],standard.slr.rcp26.50[iproj],type='l', ann=FALSE, # median
       col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]),lwd=2,
       xlim=c(2000,2100), ylim=c(0,4), 
       xaxt='n', yaxt='n', xaxs='i', yaxs='i'
  );
  axis(1, seq(2000,2100,by=20)); 
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(0,4,by=1), lab=c('0','','2','','4'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' h) Standard Projections')), line=.25, cex=.9, adj=0);
  
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(standard.slr.rcp26.95[iproj],rev(standard.slr.rcp26.05[iproj])),
          col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3],opacity), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot densities
  plot(x=standard.density.2100$y,
       y=standard.density.2100$x,
       xaxt='n', yaxt='n', xaxs='i', yaxs='i',
       ylim=c(0,4),
       xlim=c(0,1.25),
       type="n",
       ann=FALSE,
       # xlab='Probability density',
       # ylab='[m SLE]', 
       lwd=2, col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]), 
       xaxt='n')
  axis(2, seq(0,4,by=1), lab=c('0','','2','','4'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Probability density', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' i) Standard SLR at 2100')), line=.25, cex=.9, adj=0);
  polygon(x=c(rep(0,length(standard.density.2100$x)),rev(standard.density.2100$y)),
          y=c(standard.density.2100$x,rev(standard.density.2100$x)),
          col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
          # border=NA
  )
  abline(h=0,col="black",lty=2)
  
  # complete #######################################
  
  # plot hindcasts 
  plot(mod.time, complete.slr.50, type='l', 
       col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' j) Complete Hindcasts')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(complete.slr.95,rev(complete.slr.05)),
          col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  abline(h=0,col="black",lty=2) 
  
  
  # plot complete projections
  plot(t.proj[iproj],complete.slr.rcp26.50[iproj],type='l', ann=FALSE, # median
       col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]),lwd=2,
       xlim=c(2000,2100), ylim=c(0,4), 
       xaxt='n', yaxt='n', xaxs='i', yaxs='i'
  );
  axis(1, seq(2000,2100,by=20)); 
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(0,4,by=1), lab=c('0','','2','','4'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' k) Complete Projections')), line=.25, cex=.9, adj=0);
  
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.slr.rcp26.95[iproj],rev(complete.slr.rcp26.05[iproj])),
          col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
  abline(h=0,col="black",lty=2)
  
  # plot densities
  plot(x=complete.density.2100$y,
       y=complete.density.2100$x,
       xaxt='n', yaxt='n', xaxs='i', yaxs='i',
       ylim=c(0,4),
       xlim=c(0,1.25),
       type="n",
       ann=FALSE,
       # xlab='Probability density',
       # ylab='[m SLE]', 
       lwd=2, col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), 
       xaxt='n')
  axis(2, seq(0,4,by=1), lab=c('0','','2','','4'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Probability density', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' l) Complete SLR at 2100')), line=.25, cex=.9, adj=0);
  polygon(x=c(rep(0,length(complete.density.2100$x)),rev(complete.density.2100$y)),
          y=c(complete.density.2100$x,rev(complete.density.2100$x)),
          col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
          # border=NA
  )
  abline(h=0,col="black",lty=2)
  
  
  
  # LEGEND
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top",inset = 0,
         legend = c("median, model" ,
                    "5-95% range, model",
                    "observations"#,
                    # "Probability density"
                    # "2-sigma range, observations"
         ),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,NA#,
               #8
         ), bty='n', cex=1.2,
         pch=c(NA,NA,20,NA),
         col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
               rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
               rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
         ), 
         horiz = FALSE)
}
##==============================================================================
##==============================================================================


##########################################################################################################
## Hindcasts, Projections, 2100 Distribution (LEAN)


if(TRUE){
layout(cbind(c(1,4,6,8),c(2,4,6,8),c(3,5,7,8)))
par(mai=c(.5,.5,.3,.1)) #c(bottom, left, top, right) 

######################### ANTARCTIC ICE SHEET ################################


# complete #######################################

# plot hindcasts 
ipaleo=which(t.paleo==-149999):which(t.paleo==1)
plot(t.paleo[ipaleo], complete.ais.paleo.50[ipaleo], type='l', 
     col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), lwd=2, xlab='',
     ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
     # ylim=c(-20,10),
     ylim=c(-20,16),
     cex.lab=1.2, cex.axis=1.2);
mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
# mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' a) AIS Hindcasts')), line=.25, cex=.9, adj=0);
polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(complete.ais.paleo.95[ipaleo],rev(complete.ais.paleo.05[ipaleo])), 
        col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
for (i in 1:3) {
  polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
          c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
}
i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
# lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
abline(h=0,col="black",lty=2)

# plot complete projections
plot(t.proj[iproj],complete.ais.rcp26.50[iproj],type='l', ann=FALSE, # median
     col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]),lwd=2,
     xlim=c(2000,2100), ylim=c(0,1), 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i'
);
axis(1, seq(2000,2100,by=20)); 
# axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
axis(2, seq(0,1,by=.5), lab=c('0','0.5','1'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Year', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' b) AIS Projections')), line=.25, cex=.9, adj=0);

polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.ais.rcp26.95[iproj],rev(complete.ais.rcp26.05[iproj])),
        col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
abline(h=0,col="black",lty=2)

# plot densities
plot(x=complete.density.ais.2100$y,
     y=complete.density.ais.2100$x,
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     ylim=c(0,1),
     xlim=c(0,4.35),
     type="n",
     ann=FALSE,
     # xlab='Probability density',
     # ylab='[m SLE]', 
     lwd=2, col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), 
     xaxt='n')
axis(2, seq(0,1,by=.5), lab=c('0','0.5','1'));
# mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Probability density', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' c) AIS SLR at 2100')), line=.25, cex=.9, adj=0);
polygon(x=c(rep(0,length(complete.density.ais.2100$x)),rev(complete.density.ais.2100$y)),
        y=c(complete.density.ais.2100$x,rev(complete.density.ais.2100$x)),
        col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
        # border=NA
)
abline(h=0,col="black",lty=2)

######################### GREENLAND ICE SHEET ################################

# complete #######################################
ylim <- max(complete.gis.rcp26.95[iproj])
hind.time <- mod.time[1:151]
hind.complete.gis.50 <- complete.gis.50[1:151]
hind.complete.gis.95 <- complete.gis.95[1:151]
hind.complete.gis.05 <- complete.gis.05[1:151]

# plot hindcasts 
plot(hind.time, hind.complete.gis.50, type='l', # median
     col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2100), 
     # xaxt='n', 
     yaxt='n', xaxs='i', yaxs='i',
     # ylim=c(-.003,.01),
     ylim=c(-.05,2.2),
     cex.lab=1.2, cex.axis=1.2);
mtext(side=1, text='Year', line=2.3, cex=.9);
# mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' d) GIS hindcasts and projections')), line=.25, cex=.9, adj=0);

polygon(c(hind.time,rev(hind.time)), c(hind.complete.gis.95,rev(hind.complete.gis.05)), # polygon
        col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);

# plot complete projections
lines(t.proj[iproj],complete.gis.rcp26.50[iproj], ann=FALSE, # median
      col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]),lwd=2,
      # xlim=c(2000,2100), ylim=c(-1,3), 
      # xaxt='n', yaxt='n', xaxs='i', yaxs='i'
);

# observations
# lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
points(obs.gis.time, obs.gis.norm, pch=20, 
       col=rgb(mycol[col_obs,1],mycol[col_obs,2],mycol[col_obs,3]));
# lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
# abline(h=0,col="black",lty=2)

# axis(1, seq(2000,2100,by=20)); 
# axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
axis(2, seq(0,2,by=0.5), lab=c('0','','1','','2'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Year', line=2.3, cex=.9);
# mtext(side=3, text=expression(bold(' h) complete Projections')), line=.25, cex=.9, adj=0);

polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.gis.rcp26.95[iproj],rev(complete.gis.rcp26.05[iproj])),
        col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
# abline(h=0,col="black",lty=2)

# plot densities
plot(x=complete.density.gis.2100$y,
     y=complete.density.gis.2100$x,
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     ylim=c(-.05,2.2),
     xlim=c(0,1.45),
     type="n",
     ann=FALSE,
     # xlab='Probability density',
     # ylab='[m SLE]', 
     lwd=2, col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), 
     xaxt='n')
axis(2, seq(0,2,by=0.5), lab=c('0','','1','','2'));
# mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Probability density', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' e) GIS SLR at 2100')), line=.25, cex=.9, adj=0);
polygon(x=c(rep(0,length(complete.density.gis.2100$x)),rev(complete.density.gis.2100$y)),
        y=c(complete.density.gis.2100$x,rev(complete.density.gis.2100$x)),
        col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
        # border=NA
)

######################### TOTAL SEA LEVEL ################################

# complete #######################################

# plot hindcasts 
plot(mod.time, complete.slr.50, type='l', 
     col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), lwd=2, xlab='',
     ylab='', xlim=c(1850,2100), 
     ylim=c(-.3,3.25), 
     # xaxt='n',
     yaxt='n',
     xaxs='i',
     yaxs='i',
     cex.lab=1.2, cex.axis=1.2);
mtext(side=1, text='Year', line=2.3, cex=.9);
# mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
mtext(side=2, text='[m]', line=2.3, cex=.9);
# mtext(side=3, text=expression(bold(' j) Complete Hindcasts')), line=.25, cex=.9, adj=0);
polygon(c(mod.time,rev(mod.time)), c(complete.slr.95,rev(complete.slr.05)),
        col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);

abline(h=0,col="black",lty=2) 


# plot complete projections
lines(t.proj[iproj],complete.slr.rcp26.50[iproj],type='l', ann=FALSE, # median
     col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]),lwd=2,
     # xlim=c(2000,2100), ylim=c(0,4), 
     # xaxt='n', yaxt='n', xaxs='i', yaxs='i'
);


# axis(1, seq(2000,2100,by=20)); 
axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
# axis(2, seq(0,4,by=1), lab=c('0','','2','','4'));
# mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
# mtext(side=1, text='Year', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' f) Total Sea Level hindcasts and projections')), line=.25, cex=.9, adj=0);

polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.slr.rcp26.95[iproj],rev(complete.slr.rcp26.05[iproj])),
        col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
# abline(h=0,col="black",lty=2)

# plot observations
points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));


# plot densities
plot(x=complete.density.2100$y,
     y=complete.density.2100$x,
     xaxt='n', yaxt='n', xaxs='i', yaxs='i',
     ylim=c(-.3,3.25),
     xlim=c(0,1.25),
     type="n",
     ann=FALSE,
     # xlab='Probability density',
     # ylab='[m SLE]', 
     lwd=2, col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), 
     xaxt='n')
axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
# mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Probability density', line=2.3, cex=.9);
mtext(side=3, text=expression(bold(' g) Total SLR at 2100')), line=.25, cex=.9, adj=0);
polygon(x=c(rep(0,length(complete.density.2100$x)),rev(complete.density.2100$y)),
        y=c(complete.density.2100$x,rev(complete.density.2100$x)),
        col=rgb(mycol[col_density,1],mycol[col_density,2],mycol[col_density,3],opacity)#,
        # border=NA
)
abline(h=0,col="black",lty=2)

# LEGEND
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("median, model" ,
                  "5-95% range, model",
                  "observations"#,
                  # "Probability density"
                  # "2-sigma range, observations"
       ),
       # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
       lwd=c(2,8,NA#,
             #8
       ), bty='n', cex=1.2,
       pch=c(NA,NA,20,NA),
       col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
             rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
             rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
             # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])
             # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
       ), 
       horiz = TRUE)

}


##==============================================================================
##==============================================================================
## Figure 2
## Hindcasts for AIS, GIS and SLR

## mycol colors, by row:
# 1 black
# 2 dark blue
# 3 aqua
# 4 pink
# 5 light orange
# 6 purple
# 7 blue
# 8 light purple
# 9 light blue
# 10 very light blue
# 11 dark red
# 12 brown
# 13 dark orange
# 14 neon green
# 15 neon yellow

col_priors <- 2
col_experts <- 2
col_standard <- 2
col_complete <- 2
col_obs <- 11
col_density <- 1
opacity <- 0.5

if(TRUE){
  # layout(cbind(c(1,4,7),c(2,5,7), c(3,6,7))) # TWO ROWS OF THREE
  layout(cbind(c(1,3,5,7),c(2,4,6,7))) # THREE ROWS OF TWO
  par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) 
  # omi = c(0,0,0,2) #for right side space
  
  ######################### STANDARD ################################
  # a) Standard AIS Hindcasts
  # standard #######################################
  
  # plot hindcasts 
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], standard.ais.paleo.50[ipaleo], type='l', 
       col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold('a')), line=.25, cex=1.1, adj=0);
  text(.5* par('usr')[1], .8 * par('usr')[4], labels = 'AIS (BRICK)', cex=1.4 )
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(standard.ais.paleo.95[ipaleo],rev(standard.ais.paleo.05[ipaleo])), 
          col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  abline(h=0,col="black",lty=2)
  
  ######################### COMPLETE ################################
  # a) Complete AIS Hindcasts
  # complete #######################################
  
  # plot hindcasts 
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], complete.ais.paleo.50[ipaleo], type='l', 
       col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold('b')), line=.25, cex=1.1, adj=0);
  text(.5* par('usr')[1], .8 * par('usr')[4], labels = 'AIS (BRICK + SEJ)', cex=1.4 )
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(complete.ais.paleo.95[ipaleo],rev(complete.ais.paleo.05[ipaleo])), 
          col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  abline(h=0,col="black",lty=2)
  
  # b) Standard GIS Hindcasts
  # plot hindcasts 
  plot(mod.time, standard.gis.50, type='l', # median
       col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold('c')), line=.25, cex=1.1, adj=0);
  text(1* par('usr')[1]+35, .85 * par('usr')[4], labels = 'GIS (BRICK)', cex=1.4 )
  polygon(c(mod.time,rev(mod.time)), c(standard.gis.95,rev(standard.gis.05)), # polygon
          col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3],opacity), border=NA);
  
  # observations
  # lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gis.time, obs.gis.norm, pch=20, 
         col=rgb(mycol[col_obs,1],mycol[col_obs,2],mycol[col_obs,3]));
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  abline(h=0,col="black",lty=2)
  
  # b) complete GIS Hindcasts
  # plot hindcasts 
  plot(mod.time, complete.gis.50, type='l', # median
       col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold('d')), line=.25, cex=1.1, adj=0);
  text(1* par('usr')[1]+35, .85 * par('usr')[4], labels = 'GIS (BRICK + SEJ)', cex=1.4 )
  polygon(c(mod.time,rev(mod.time)), c(complete.gis.95,rev(complete.gis.05)), # polygon
          col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
  
  # observations
  # lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gis.time, obs.gis.norm, pch=20, 
         col=rgb(mycol[col_obs,1],mycol[col_obs,2],mycol[col_obs,3]));
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  abline(h=0,col="black",lty=2)
  
  # c) Standard SLR Hindcasts
  # plot hindcasts 
  plot(mod.time, standard.slr.50, type='l', 
       col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold('e')), line=.25, cex=1.1, adj=0);
  text(1* par('usr')[1]+90, .75 * par('usr')[4], labels = 'GMSL (BRICK)', cex=1.4 )
  polygon(c(mod.time,rev(mod.time)), c(standard.slr.95,rev(standard.slr.05)), 
          col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  abline(h=0,col="black",lty=2)
  
  # c) complete SLR Hindcasts
  # plot hindcasts 
  plot(mod.time, complete.slr.50, type='l', 
       col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold('f')), line=.25, cex=1.1, adj=0);
  text(1* par('usr')[1]+90, .75 * par('usr')[4], labels = 'GMSL (BRICK + SEJ)', cex=1.4 )
  polygon(c(mod.time,rev(mod.time)), c(complete.slr.95,rev(complete.slr.05)), 
          col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  # lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  abline(h=0,col="black",lty=2)



  # LEGEND
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top",inset = 0,
         legend = c("median, model" ,
                    "5-95% range, model",
                    "observations"#,
                    # "Probability density"
                    # "2-sigma range, observations"
         ),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,NA#,
               #8
         ), bty='n', cex=1.5,
         pch=c(NA,NA,20,NA),
         col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
               rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
               rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
         ), 
         horiz = FALSE)



}



##==============================================================================
##==============================================================================
## Figure 3 - Estimates of GMSL rise by 2100

# AIS Contributions to SLR
sl.2100.ais.complete <- complete.ais.rcp26[nrow(complete.ais.rcp26),]
sl.2100.ais.standard <- standard.ais.rcp26[nrow(standard.ais.rcp26),]
sl.2100.ais.experts <- experts.ais.rcp26[nrow(experts.ais.rcp26),]
sl.2100.ais.priors <- priors.ais.rcp26[nrow(priors.ais.rcp26),]

density.ais.complete <- density(sl.2100.ais.complete)
density.ais.standard <- density(sl.2100.ais.standard)
density.ais.experts  <- density(sl.2100.ais.experts)
density.ais.priors   <- density(sl.2100.ais.priors)

# AIS Contributions to SLR
sl.2100.gis.complete <- complete.gis.rcp26[nrow(complete.gis.rcp26),]
sl.2100.gis.standard <- standard.gis.rcp26[nrow(standard.gis.rcp26),]
sl.2100.gis.experts <- experts.gis.rcp26[nrow(experts.gis.rcp26),]
sl.2100.gis.priors <- priors.gis.rcp26[nrow(priors.gis.rcp26),]

density.gis.complete <- density(sl.2100.gis.complete)
density.gis.standard <- density(sl.2100.gis.standard)
density.gis.experts  <- density(sl.2100.gis.experts)
density.gis.priors   <- density(sl.2100.gis.priors)

# Total SLR
sl.2100.complete <- complete.slr.rcp26[nrow(complete.slr.rcp26),]
sl.2100.standard <- standard.slr.rcp26[nrow(standard.slr.rcp26),]
sl.2100.experts  <- experts.slr.rcp26[nrow(experts.slr.rcp26),]
sl.2100.priors   <- priors.slr.rcp26[nrow(priors.slr.rcp26),]

density.complete <- density(sl.2100.complete)
density.standard <- density(sl.2100.standard)
density.experts  <- density(sl.2100.experts)
density.priors   <- density(sl.2100.priors)

# Thermal Expansion
te.2100.complete <- complete.te.rcp26[nrow(complete.te.rcp26),]
te.2100.standard <- standard.te.rcp26[nrow(standard.te.rcp26),]
te.2100.experts  <- experts.te.rcp26[nrow(experts.te.rcp26),]
te.2100.priors   <- priors.te.rcp26[nrow(priors.te.rcp26),]

density.te.complete <- density(te.2100.complete)
density.te.standard <- density(te.2100.standard)
density.te.experts  <- density(te.2100.experts)
density.te.priors   <- density(te.2100.priors)

# Glaciers and Ice Caps
gsic.2100.complete <- complete.gsic.rcp26[nrow(complete.gsic.rcp26),]
gsic.2100.standard <- standard.gsic.rcp26[nrow(standard.gsic.rcp26),]
gsic.2100.experts  <- experts.gsic.rcp26[nrow(experts.gsic.rcp26),]
gsic.2100.priors   <- priors.gsic.rcp26[nrow(priors.gsic.rcp26),]

density.gsic.complete <- density(gsic.2100.complete)
density.gsic.standard <- density(gsic.2100.standard)
density.gsic.experts  <- density(gsic.2100.experts)
density.gsic.priors   <- density(gsic.2100.priors)


## mycol colors, by row:
# 1 black
# 2 dark blue
# 3 aqua
# 4 pink
# 5 light orange
# 6 purple
# 7 blue
# 8 light purple
# 9 light blue
# 10 very light blue
# 11 dark red
# 12 brown
# 13 dark orange
# 14 neon green
# 15 neon yellow
mycol_experts <- 5
mycol_standard <- 1
mycol_complete <- 9
mycol_priors <- 11
mycol_sej2018 <- 5
mycol_ar6 <- 11


if (TRUE){
# dev.off()
layout(cbind(c(1,1,4),c(1,1,4),c(1,1,4),c(2,3,4)))

# width by height: 658 by 693

par(mai=c(.15,.3,.3,.2), #c(bottom, left, top, right)
    oma=c(1,2,1,3) #c(bottom, left, top, right)
    )  
  

mylwd = 2
mycexlab = 1.2
# myylim = c(-1,4.5)

## TOTAL SLR
myxlim <- max(pmax(density.complete$y,density.standard$y,density.experts$y,density.priors$y))
plot(density.standard$y, density.standard$x, 
     xlim=c(0,1.4), # c(0,1.6 with AR6)
     ylim=c(.5,3.5),
     axes = FALSE,
     type="l",
     lwd=mylwd, col=mycol.rgb[mycol_standard],
     yaxt='n',
     ylab='',
     xlab=''
)
# abline(h=0, lty="dashed")
axis(1,labels=FALSE,tick=FALSE)
axis(2, seq(0.5,3.5,by=.5), lab=c('','1','','2','','3',''))
# title(ylab="Probability density",cex.lab=mycexlab)
box()
# mtext("Projected global mean sea level in 2100",side=2,line=3)
# mtext("relative to 1986-2005 average [m]",side=2,line=2)
mtext("[m]",side=2,line=2, cex=.9)
mtext(side=3, text=expression(bold('a')), line=.25, cex=mycexlab, adj=0);
text(.2 * par('usr')[1]+.1, .97 * par('usr')[4], labels = 'GMSL', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
mtext("Probability density",side=1,line=0.5)


# lines(density.experts$y, density.experts$x,
#       lwd=mylwd, col=mycol.rgb[mycol_experts], lty="dashed"
# )
# 
# lines(density.priors$y, density.priors$x,
#       lwd=mylwd, col=mycol.rgb[mycol_priors], lty="dashed")

lines(density.complete$y,density.complete$x,
      lwd=mylwd, col=mycol.rgb[mycol_complete]
)

text(x=1.1,y=2.2,labels=expression(bold("BRICK")), cex=1.4)
text(x=0.75,y=2.5,labels=expression(bold("BRICK+SEJ")), cex=1.4, col=mycol.rgb[mycol_complete])

# Bamber et al 2019
# 0.79 - 1.74 m above 2000 CE likely range
# 0.62 - 2.38 m above 2000 CE 5-95%
# abline(v=1.3)
# abline(v=1.4)

segments(x0=1.35, y0=0.62, x1 = 1.35, y1 = 2.38,
         col=mycol.rgb[mycol_sej2018],
         lwd=2
         )

polygon(x = c(1.3, 1.4, 1.4, 1.3),                             # X-Coordinates of polygon
        y = c(.79, 0.79, 1.74, 1.74),                         # Y-Coordinates of polygon
        col = mycol.rgb[mycol_sej2018],                            # Color of polygon
        border = mycol.rgb[mycol_sej2018],                                      # Color of polygon border
        lwd = 1,                                               # Thickness of border
)

text(x=1.2,y=1.3,labels=expression(bold("SEJ")), cex=1.4, col=mycol.rgb[mycol_sej2018])

# AR6 # .63 m to 1.01 m AR6 8.5 page from IPCC AR6 WGI Chapter 9 page 1299 (pdf page 89)
# .82 to 1.19 MICI
# since 1995 - 2014
# abline(v=1.5)
# abline(v=1.6)
# segments(x0=1.53, y0=0.82, x1 = 1.53, y1 = 1.19,
#          col=mycol.rgb[mycol_ar6],
#          lwd=2
# )
# polygon(x = c(1.48, 1.58, 1.58, 1.48),                             # X-Coordinates of polygon
#         y = c(0.63, 0.63, 1.01, 1.01),                         # Y-Coordinates of polygon
#         col = mycol.rgb[mycol_ar6],                            # Color of polygon
#         border = mycol.rgb[mycol_ar6],                                      # Color of polygon border
#         lwd = 1,                                               # Thickness of border
#         )                                             

## AIS ########################################################################
myxlim <- max(pmax(density.ais.complete$y,density.ais.standard$y,density.ais.experts$y,density.ais.priors$y)) # 4.3218
plot(density.ais.standard$y, density.ais.standard$x, 
     # xlim=c(0,myxlim),
     xlim=c(0,6), # c(0,7) with AR6
     ylim=c(-.15,1.5),
     axes = FALSE,
     type="l",
     lwd=mylwd, col=mycol.rgb[mycol_standard],
     yaxt='n',
     ylab='',
     xlab=''
)
# abline(h=0, lty="dashed")
axis(1,labels=FALSE,tick=FALSE)
axis(2, seq(0,1.5,by=0.25), lab=c('0','','0.5','','1','','1.5'))
box()
mtext("[m]",side=2,line=2, cex=.8)
mtext(side=3, text=expression(bold('b')), line=.25, cex=mycexlab, adj=0);
text(.2 * par('usr')[1]+1, .93 * par('usr')[4], labels = 'AIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
mtext("Probability density",side=1,line=0.2, cex=.8)

# lines(density.ais.experts$y, density.ais.experts$x,
#       lwd=mylwd, col=mycol.rgb[mycol_experts], lty="dashed"
# )
# 
# lines(density.ais.priors$y, density.ais.priors$x,
#       lwd=mylwd, col=mycol.rgb[mycol_priors], lty="dashed")

lines(density.ais.complete$y,density.ais.complete$x,
      lwd=mylwd, col=mycol.rgb[mycol_complete]
)

# AR6 # .03 m to .34 m AR6 8.5 likely page from IPCC AR6 WGI Chapter 9 page 89
# .19 to .53 m MICI
# since 1995 - 2014
# abline(v=5)
# abline(v=5.7)
# segments(x0=6.45, y0=.19, x1 = 6.45, y1 = .53,
#          col=mycol.rgb[mycol_ar6],
#          lwd=2
# )
# 
# polygon(x = c(6.1, 6.8, 6.8, 6.1),                             # X-Coordinates of polygon
#         y = c(0.03, .03, .34, .34),                         # Y-Coordinates of polygon
#         col = mycol.rgb[mycol_ar6],                            # Color of polygon
#         border = mycol.rgb[mycol_ar6],                         # Color of polygon border
#         lwd = 1,                                               # Thickness of border
# )

# Bamber et al 2019
# 0.02 - 0.57 m above 2000 CE likely range
# -.11 - 1.32 m above 2000 CE 5-95%
# abline(v=6.1)
# abline(v=6.8)

segments(x0=5.35, y0=-0.11, x1 = 5.35, y1 = 1.32,
         col=mycol.rgb[mycol_sej2018],
         lwd=2
)

polygon(x = c(5, 5.7, 5.7, 5),                             # X-Coordinates of polygon
        y = c(.02, .02, .57, .57),                         # Y-Coordinates of polygon
        col = mycol.rgb[mycol_sej2018],                            # Color of polygon
        border = mycol.rgb[mycol_sej2018],                                      # Color of polygon border
        lwd = 1,                                               # Thickness of border
)

## GIS ########################################################################
myxlim <- max(pmax(density.gis.complete$y,density.gis.standard$y,density.gis.experts$y,density.gis.priors$y)) #1.378
plot(density.gis.standard$y, density.gis.standard$x, 
     xlim=c(0,1.75), #c(0,2.05) with AR6
     ylim=c(0,2.5),
     axes = FALSE,
     type="l",
     lwd=mylwd, col=mycol.rgb[mycol_standard],
     yaxt='n',
     ylab='',
     xlab=''
)
# abline(h=0, lty="dashed")
axis(1,labels=FALSE,tick=FALSE)
axis(2, seq(0,2.5,by=0.5), lab=c('0','','1','','2',''))
box()
mtext("[m]",side=2,line=2, cex=.8)
mtext(side=3, text=expression(bold('c')), line=.25, cex=mycexlab, adj=0);
text(.2 * par('usr')[1]+.3, .93 * par('usr')[4], labels = 'GIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
mtext("Probability density",side=1,line=0.2, cex=.8)

# lines(density.gis.experts$y, density.gis.experts$x,
#       lwd=mylwd, col=mycol.rgb[mycol_experts], lty="dashed"
# )
# 
# lines(density.gis.priors$y, density.gis.priors$x,
#       lwd=mylwd, col=mycol.rgb[mycol_priors], lty="dashed")

lines(density.gis.complete$y,density.gis.complete$x,
      lwd=mylwd, col=mycol.rgb[mycol_complete]
)

# # AR6 # .09 m to .18 m AR6 8.5 page from IPCC AR6 WGI Chapter 9 page 89
# # since 1995 - 2014
# # abline(v=1.5)
# # abline(v=1.7)
# polygon(x = c(1.82, 2.02, 2.02, 1.82),                             # X-Coordinates of polygon
#         y = c(0.09, .09, .18, .18),                         # Y-Coordinates of polygon
#         col = mycol.rgb[mycol_ar6],                            # Color of polygon
#         border = mycol.rgb[mycol_ar6],                         # Color of polygon border
#         lwd = 1,                                               # Thickness of border
# )

# Bamber et al 2019
# 0.1 - 0.60 m above 2000 CE likely range
# .02 - .99 m above 2000 CE 5-95%
# abline(v=1.82)
# abline(v=2.02)

segments(x0=1.6, y0=.02, x1 = 1.6, y1 = .99,
         col=mycol.rgb[mycol_sej2018],
         lwd=2
)

polygon(x = c(1.5, 1.7, 1.7, 1.5),                             # X-Coordinates of polygon
        y = c(.1, .1, .60, .60),                         # Y-Coordinates of polygon
        col = mycol.rgb[mycol_sej2018],                            # Color of polygon
        border = mycol.rgb[mycol_sej2018],                                      # Color of polygon border
        lwd = 1,                                               # Thickness of border
)




## LEGEND #############################################################################
plot(1, type = "n", axes=FALSE, xlab="", ylab=""
     # ,main="Global mean sea level at 2100 [m]"
)
legend(x = "top",inset = 0,
       # legend = c("Wide priors combined with data (high-temperature scenario)" , 
       #            "Wide priors combined with data and expert assessment (high-temperature scenario)",
       #            "17th-83rd percentile - SEJ2018 expert assessment (high-temperature scenario)",
       #            "5-95th percentile - SEF2018 expert assessment (high-temperature scenario)"
       #            # "IPCC AR6 likely range (RCP 8.5)"
       # ),
       
       # legend = c("BRICK" , 
       #            "BRICK + SEJ2018",
       #            expression("5"^th~"-95"^th~" percentile, SEJ2018"),
       #            expression("17"^th~"-83"^rd~" percentile, SEJ2018")
       #            # "IPCC AR6 likely range (RCP 8.5)"
       # ), 
       
       legend = c("data-model calibration" , 
                  "expert-data-model calibration",
                  expression("5"^th~"-95"^th~" percentile, expert assessment"),
                  expression("17"^th~"-83"^rd~" percentile, expert assessment")
                  # "IPCC AR6 likely range (RCP 8.5)"
       ), 
       
       
       lwd=c(2,2,2,NA),
       pch=c(NA,NA,NA,15),
       lty=c("solid",
             "solid",
             "solid",
             NA
             ),
       bty='n', cex=1.5,
       col=c(mycol.rgb[mycol_standard],
             mycol.rgb[mycol_complete],
             mycol.rgb[mycol_experts],
             mycol.rgb[mycol_experts]
             # "#c77d7dff"
             ), 
       horiz = FALSE,
       
)




}




##==============================================================================
##==============================================================================
## Extended Data Figure 4 - Estimates of GMSL rise by 2100 with priors

# AIS Contributions to SLR
sl.2100.ais.complete <- complete.ais.rcp26[nrow(complete.ais.rcp26),]
sl.2100.ais.standard <- standard.ais.rcp26[nrow(standard.ais.rcp26),]
sl.2100.ais.experts <- experts.ais.rcp26[nrow(experts.ais.rcp26),]
sl.2100.ais.priors <- priors.ais.rcp26[nrow(priors.ais.rcp26),]

density.ais.complete <- density(sl.2100.ais.complete)
density.ais.standard <- density(sl.2100.ais.standard)
density.ais.experts  <- density(sl.2100.ais.experts)
density.ais.priors   <- density(sl.2100.ais.priors)

# AIS Contributions to SLR
sl.2100.gis.complete <- complete.gis.rcp26[nrow(complete.gis.rcp26),]
sl.2100.gis.standard <- standard.gis.rcp26[nrow(standard.gis.rcp26),]
sl.2100.gis.experts <- experts.gis.rcp26[nrow(experts.gis.rcp26),]
sl.2100.gis.priors <- priors.gis.rcp26[nrow(priors.gis.rcp26),]

density.gis.complete <- density(sl.2100.gis.complete)
density.gis.standard <- density(sl.2100.gis.standard)
density.gis.experts  <- density(sl.2100.gis.experts)
density.gis.priors   <- density(sl.2100.gis.priors)

# Total SLR
sl.2100.complete <- complete.slr.rcp26[nrow(complete.slr.rcp26),]
sl.2100.standard <- standard.slr.rcp26[nrow(standard.slr.rcp26),]
sl.2100.experts  <- experts.slr.rcp26[nrow(experts.slr.rcp26),]
sl.2100.priors   <- priors.slr.rcp26[nrow(priors.slr.rcp26),]

density.complete <- density(sl.2100.complete)
density.standard <- density(sl.2100.standard)
density.experts  <- density(sl.2100.experts)
density.priors   <- density(sl.2100.priors)

# Thermal Expansion
te.2100.complete <- complete.te.rcp26[nrow(complete.te.rcp26),]
te.2100.standard <- standard.te.rcp26[nrow(standard.te.rcp26),]
te.2100.experts  <- experts.te.rcp26[nrow(experts.te.rcp26),]
te.2100.priors   <- priors.te.rcp26[nrow(priors.te.rcp26),]

density.te.complete <- density(te.2100.complete)
density.te.standard <- density(te.2100.standard)
density.te.experts  <- density(te.2100.experts)
density.te.priors   <- density(te.2100.priors)

# Glaciers and Ice Caps
gsic.2100.complete <- complete.gsic.rcp26[nrow(complete.gsic.rcp26),]
gsic.2100.standard <- standard.gsic.rcp26[nrow(standard.gsic.rcp26),]
gsic.2100.experts  <- experts.gsic.rcp26[nrow(experts.gsic.rcp26),]
gsic.2100.priors   <- priors.gsic.rcp26[nrow(priors.gsic.rcp26),]

density.gsic.complete <- density(gsic.2100.complete)
density.gsic.standard <- density(gsic.2100.standard)
density.gsic.experts  <- density(gsic.2100.experts)
density.gsic.priors   <- density(gsic.2100.priors)


## mycol colors, by row:
# 1 black
# 2 dark blue
# 3 aqua
# 4 pink
# 5 light orange
# 6 purple
# 7 blue
# 8 light purple
# 9 light blue
# 10 very light blue
# 11 dark red
# 12 brown
# 13 dark orange
# 14 neon green
# 15 neon yellow
mycol_experts <- 9
mycol_standard <- 1
mycol_complete <- 9
mycol_priors <- 1
mycol_sej2018 <- 5
mycol_ar6 <- 11


if (TRUE){
  # dev.off()
  layout(cbind(c(1,1,4),c(1,1,4),c(1,1,4),c(2,3,4)))
  
  # width by height: 658 by 693
  
  par(mai=c(.15,.3,.3,.2), #c(bottom, left, top, right)
      oma=c(1,2,1,3) #c(bottom, left, top, right)
  )  
  
  
  mylwd = 2
  mylwd_priors = 1
  mycexlab = 1.2
  # myylim = c(-1,4.5)
  
  ## TOTAL SLR
  myxlim <- max(pmax(density.complete$y,density.standard$y,density.experts$y,density.priors$y))
  plot(density.standard$y, density.standard$x, 
       xlim=c(0,1.4), # c(0,1.6 with AR6)
       ylim=c(.5,3.5),
       axes = FALSE,
       type="l",
       lwd=mylwd, col=mycol.rgb[mycol_standard],
       yaxt='n',
       ylab='',
       xlab=''
  )
  # abline(h=0, lty="dashed")
  axis(1,labels=FALSE,tick=FALSE)
  axis(2, seq(0.5,3.5,by=.5), lab=c('','1','','2','','3',''))
  # title(ylab="Probability density",cex.lab=mycexlab)
  box()
  # mtext("Projected global mean sea level in 2100",side=2,line=3)
  # mtext("relative to 1986-2005 average [m]",side=2,line=2)
  mtext("[m]",side=2,line=2, cex=.9)
  mtext(side=3, text=expression(bold('a')), line=.25, cex=mycexlab, adj=0);
  text(.2 * par('usr')[1]+.15, .97 * par('usr')[4], labels = 'GMSL', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  mtext("Probability density",side=1,line=0.5)
  
  
  lines(density.experts$y, density.experts$x,
        lwd=mylwd_priors, col=mycol.rgb[mycol_experts], lty="dashed"
  )

  lines(density.priors$y, density.priors$x,
        lwd=mylwd_priors, col=mycol.rgb[mycol_priors], lty="dashed")
  
  lines(density.complete$y,density.complete$x,
        lwd=mylwd, col=mycol.rgb[mycol_complete]
  )
  
  text(x=1.1,y=2.2,labels=expression(bold("BRICK")), cex=1.4)
  text(x=0.75,y=2.5,labels=expression(bold("BRICK+SEJ")), cex=1.4, col=mycol.rgb[mycol_complete])
  text(x=0.45,y=0.5,labels="BRICK priors", cex=1.4)
  text(x=1,y=0.7,labels="Expert-informed priors", cex=1.4, col=mycol.rgb[mycol_complete])
  text(x=1.2,y=1.3,labels=expression(bold("SEJ")), cex=1.4, col=mycol.rgb[mycol_sej2018])
  
  # Bamber et al 2019
  # 0.79 - 1.74 m above 2000 CE likely range
  # 0.62 - 2.38 m above 2000 CE 5-95%
  # abline(v=1.3)
  # abline(v=1.4)
  
  segments(x0=1.35, y0=0.62, x1 = 1.35, y1 = 2.38,
           col=mycol.rgb[mycol_sej2018],
           lwd=2
  )
  
  polygon(x = c(1.3, 1.4, 1.4, 1.3),                             # X-Coordinates of polygon
          y = c(.79, 0.79, 1.74, 1.74),                         # Y-Coordinates of polygon
          col = mycol.rgb[mycol_sej2018],                            # Color of polygon
          border = mycol.rgb[mycol_sej2018],                                      # Color of polygon border
          lwd = 1,                                               # Thickness of border
  )
  
  
  # AR6 # .63 m to 1.01 m AR6 8.5 page from IPCC AR6 WGI Chapter 9 page 1299 (pdf page 89)
  # .82 to 1.19 MICI
  # since 1995 - 2014
  # abline(v=1.5)
  # abline(v=1.6)
  # segments(x0=1.53, y0=0.82, x1 = 1.53, y1 = 1.19,
  #          col=mycol.rgb[mycol_ar6],
  #          lwd=2
  # )
  # polygon(x = c(1.48, 1.58, 1.58, 1.48),                             # X-Coordinates of polygon
  #         y = c(0.63, 0.63, 1.01, 1.01),                         # Y-Coordinates of polygon
  #         col = mycol.rgb[mycol_ar6],                            # Color of polygon
  #         border = mycol.rgb[mycol_ar6],                                      # Color of polygon border
  #         lwd = 1,                                               # Thickness of border
  #         )                                             
  
  ## AIS ########################################################################
  myxlim <- max(pmax(density.ais.complete$y,density.ais.standard$y,density.ais.experts$y,density.ais.priors$y)) # 4.3218
  plot(density.ais.standard$y, density.ais.standard$x, 
       # xlim=c(0,myxlim),
       xlim=c(0,6), # c(0,7) with AR6
       ylim=c(-.15,1.5),
       axes = FALSE,
       type="l",
       lwd=mylwd, col=mycol.rgb[mycol_standard],
       yaxt='n',
       ylab='',
       xlab=''
  )
  # abline(h=0, lty="dashed")
  axis(1,labels=FALSE,tick=FALSE)
  axis(2, seq(0,1.5,by=0.25), lab=c('0','','0.5','','1','','1.5'))
  box()
  mtext("[m]",side=2,line=2, cex=.8)
  mtext(side=3, text=expression(bold('b')), line=.25, cex=mycexlab, adj=0);
  text(.2 * par('usr')[1]+1, .93 * par('usr')[4], labels = 'AIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  mtext("Probability density",side=1,line=0.2, cex=.8)
  
  lines(density.ais.experts$y, density.ais.experts$x,
        lwd=mylwd_priors, col=mycol.rgb[mycol_experts], lty="dashed"
  )

  lines(density.ais.priors$y, density.ais.priors$x,
        lwd=mylwd_priors, col=mycol.rgb[mycol_priors], lty="dashed")
  
  lines(density.ais.complete$y,density.ais.complete$x,
        lwd=mylwd, col=mycol.rgb[mycol_complete]
  )
  
  # AR6 # .03 m to .34 m AR6 8.5 likely page from IPCC AR6 WGI Chapter 9 page 89
  # .19 to .53 m MICI
  # since 1995 - 2014
  # abline(v=5)
  # abline(v=5.7)
  # segments(x0=6.45, y0=.19, x1 = 6.45, y1 = .53,
  #          col=mycol.rgb[mycol_ar6],
  #          lwd=2
  # )
  # 
  # polygon(x = c(6.1, 6.8, 6.8, 6.1),                             # X-Coordinates of polygon
  #         y = c(0.03, .03, .34, .34),                         # Y-Coordinates of polygon
  #         col = mycol.rgb[mycol_ar6],                            # Color of polygon
  #         border = mycol.rgb[mycol_ar6],                         # Color of polygon border
  #         lwd = 1,                                               # Thickness of border
  # )
  
  # Bamber et al 2019
  # 0.02 - 0.57 m above 2000 CE likely range
  # -.11 - 1.32 m above 2000 CE 5-95%
  # abline(v=6.1)
  # abline(v=6.8)
  
  segments(x0=5.35, y0=-0.11, x1 = 5.35, y1 = 1.32,
           col=mycol.rgb[mycol_sej2018],
           lwd=2
  )
  
  polygon(x = c(5, 5.7, 5.7, 5),                             # X-Coordinates of polygon
          y = c(.02, .02, .57, .57),                         # Y-Coordinates of polygon
          col = mycol.rgb[mycol_sej2018],                            # Color of polygon
          border = mycol.rgb[mycol_sej2018],                                      # Color of polygon border
          lwd = 1,                                               # Thickness of border
  )
  
  ## GIS ########################################################################
  myxlim <- max(pmax(density.gis.complete$y,density.gis.standard$y,density.gis.experts$y,density.gis.priors$y)) #1.378
  plot(density.gis.standard$y, density.gis.standard$x, 
       xlim=c(0,1.75), #c(0,2.05) with AR6
       ylim=c(0,2.5),
       axes = FALSE,
       type="l",
       lwd=mylwd, col=mycol.rgb[mycol_standard],
       yaxt='n',
       ylab='',
       xlab=''
  )
  # abline(h=0, lty="dashed")
  axis(1,labels=FALSE,tick=FALSE)
  axis(2, seq(0,2.5,by=0.5), lab=c('0','','1','','2',''))
  box()
  mtext("[m]",side=2,line=2, cex=.8)
  mtext(side=3, text=expression(bold('c')), line=.25, cex=mycexlab, adj=0);
  text(.2 * par('usr')[1]+.3, .93 * par('usr')[4], labels = 'GIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  mtext("Probability density",side=1,line=0.2, cex=.8)
  
  lines(density.gis.experts$y, density.gis.experts$x,
        lwd=mylwd_priors, col=mycol.rgb[mycol_experts], lty="dashed"
  )

  lines(density.gis.priors$y, density.gis.priors$x,
        lwd=mylwd_priors, col=mycol.rgb[mycol_priors], lty="dashed")
  
  lines(density.gis.complete$y,density.gis.complete$x,
        lwd=mylwd, col=mycol.rgb[mycol_complete]
  )
  
  # # AR6 # .09 m to .18 m AR6 8.5 page from IPCC AR6 WGI Chapter 9 page 89
  # # since 1995 - 2014
  # # abline(v=1.5)
  # # abline(v=1.7)
  # polygon(x = c(1.82, 2.02, 2.02, 1.82),                             # X-Coordinates of polygon
  #         y = c(0.09, .09, .18, .18),                         # Y-Coordinates of polygon
  #         col = mycol.rgb[mycol_ar6],                            # Color of polygon
  #         border = mycol.rgb[mycol_ar6],                         # Color of polygon border
  #         lwd = 1,                                               # Thickness of border
  # )
  
  # Bamber et al 2019
  # 0.1 - 0.60 m above 2000 CE likely range
  # .02 - .99 m above 2000 CE 5-95%
  # abline(v=1.82)
  # abline(v=2.02)
  
  segments(x0=1.6, y0=.02, x1 = 1.6, y1 = .99,
           col=mycol.rgb[mycol_sej2018],
           lwd=2
  )
  
  polygon(x = c(1.5, 1.7, 1.7, 1.5),                             # X-Coordinates of polygon
          y = c(.1, .1, .60, .60),                         # Y-Coordinates of polygon
          col = mycol.rgb[mycol_sej2018],                            # Color of polygon
          border = mycol.rgb[mycol_sej2018],                                      # Color of polygon border
          lwd = 1,                                               # Thickness of border
  )
  
  
  
  
  ## LEGEND #############################################################################
  plot(1, type = "n", axes=FALSE, xlab="", ylab=""
       # ,main="Global mean sea level at 2100 [m]"
  )
  legend(x = "top",inset = 0,
         # legend = c("Wide priors combined with data (high-temperature scenario)" , 
         #            "Wide priors combined with data and expert assessment (high-temperature scenario)",
         #            "17th-83rd percentile - SEJ2018 expert assessment (high-temperature scenario)",
         #            "5-95th percentile - SEF2018 expert assessment (high-temperature scenario)"
         #            # "IPCC AR6 likely range (RCP 8.5)"
         # ),
         
         # legend = c("BRICK" , 
         #            "BRICK + SEJ2018",
         #            expression("5"^th~"-95"^th~" percentile, SEJ2018"),
         #            expression("17"^th~"-83"^rd~" percentile, SEJ2018")
         #            # "IPCC AR6 likely range (RCP 8.5)"
         # ), 
         
         legend = c(
                    "standard BRICK priors",
                    "expert-informed priors",
                    "data-model calibration" , 
                    "expert-data-model calibration",
                    expression("5"^th~"-95"^th~" percentile, expert assessment"),
                    expression("17"^th~"-83"^rd~" percentile, expert assessment")
                    # "IPCC AR6 likely range (RCP 8.5)"
         ), 
         
         
         lwd=c(1,1,2,2,2,NA),
         pch=c(NA,NA,NA,NA,NA,15),
         lty=c(
               "dashed",
               "dashed",
               "solid",
               "solid",
               "solid",
               NA
         ),
         bty='n', cex=1.5,
         col=c(
               mycol.rgb[mycol_priors],
               mycol.rgb[mycol_experts],
               mycol.rgb[mycol_standard],
               mycol.rgb[mycol_complete],
               mycol.rgb[mycol_sej2018],
               mycol.rgb[mycol_sej2018]
               # "#c77d7dff"
         ), 
         horiz = FALSE,
         
  )
  
  
  
  
}











##########################################################################################################
# Extended Data Fig. 3 TE and GSIC
if(TRUE){
  
  ## mycol colors, by row:
  # 1 black
  # 2 dark blue
  # 3 aqua
  # 4 pink
  # 5 light orange
  # 6 purple
  # 7 blue
  # 8 light purple
  # 9 light blue
  # 10 very light blue
  # 11 dark red
  # 12 brown
  # 13 dark orange
  # 14 neon green
  # 15 neon yellow
  mycol_experts <- 9
  mycol_standard <- 1
  mycol_complete <- 9
  mycol_priors <- 1
  mycol_sej2018 <- 5
  mycol_ar6 <- 11
  
  mylwd = 2
  mycexlab = 1.35
  # myylim = c(-1,4.5)
  layout(cbind(c(1,3),c(2,3)))
  par(mai=c(.1,.5,.3,.3)) #c(bottom, left, top, right) 
  
  #### TE
  
  myxlim <- max(pmax(density.te.complete$y,density.te.standard$y,density.te.experts$y,density.te.priors$y))
  plot(density.te.standard$y, density.te.standard$x, 
       # xlim=c(0,myxlim),
       xlim=c(0,5),
       ylim=c(0,1),
       axes = FALSE,
       type="l",
       lwd=mylwd, col=mycol.rgb[mycol_standard],
       yaxt='n',
       ylab='',
       xlab=''
  )
  # abline(h=0, lty="dashed")
  axis(1,labels=FALSE,tick=FALSE)
  axis(2, seq(0,1.5,by=0.25), lab=c('0','','0.5','','1','','1.5'))
  box()
  mtext("[m]",side=2,line=2, cex=.8)
  mtext(side=3, text=expression(bold('a ')), line=.25, cex=1.2, adj=0);
  text(.2 * par('usr')[1]+2.2, .93 * par('usr')[4], labels = 'Thermal expansion', cex=1.4 ) 
  mtext("Probability density",side=1,line=0.2, cex=.8)
  
  lines(density.te.experts$y, density.te.experts$x,
        lwd=1, col=mycol.rgb[mycol_experts], lty="dashed"
  )

  lines(density.te.priors$y, density.te.priors$x,
        lwd=1, col=mycol.rgb[mycol_priors], lty="dashed")
  
  lines(density.te.complete$y,density.te.complete$x,
        lwd=mylwd, col=mycol.rgb[mycol_complete]
  )
  
  # text(x=3.3,y=0.33,labels=expression(bold("BRICK")), cex=1.3)
  # text(x=3.5,y=0.05,labels=expression(bold("BRICK+SEJ")), cex=1.3, col=mycol.rgb[mycol_complete])
  # text(x=0.45,y=0.5,labels="BRICK priors", cex=1.3)
  # text(x=1,y=0.7,labels="Expert-informed priors", cex=1.4, col=mycol.rgb[mycol_complete])


  
  #### GSIC
  
  myxlim <- max(pmax(density.gsic.complete$y,density.gsic.standard$y,density.gsic.experts$y,density.gsic.priors$y))
  myylim <- max(pmax(density.gsic.complete$x,density.gsic.standard$x,density.gsic.experts$x,density.gsic.priors$x))
  plot(density.gsic.standard$y, density.gsic.standard$x, 
       # xlim=c(0,myxlim),
       xlim=c(0,20),
       ylim=c(0,0.5),
       axes = FALSE,
       type="l",
       lwd=mylwd, col=mycol.rgb[mycol_standard],
       yaxt='n',
       ylab='',
       xlab=''
  )
  # abline(h=0, lty="dashed")
  axis(1,labels=FALSE,tick=FALSE)
  axis(2, seq(0,1.5,by=0.25), lab=c('0','0.25','0.5','','1','','1.5'))
  box()
  mtext("[m]",side=2,line=2, cex=.8)
  mtext(side=3, text=expression(bold('b')), line=.25, cex=1.2, adj=0);
  text(.2 * par('usr')[1]+10, .93 * par('usr')[4], labels = 'Glaciers and ice caps', cex=1.4 ) 
  mtext("Probability density",side=1,line=0.2, cex=.8)

  lines(density.gsic.experts$y, density.gsic.experts$x,
        lwd=1, col=mycol.rgb[mycol_experts], lty="dashed"
  )

  lines(density.gsic.priors$y, density.gsic.priors$x,
        lwd=1, col=mycol.rgb[mycol_priors], lty="dashed")
  
  lines(density.gsic.complete$y,density.gsic.complete$x,
        lwd=mylwd, col=mycol.rgb[mycol_complete]
  )
  
  # text(x=12,y=0.25,labels=expression(bold("BRICK")), cex=1.3)
  # text(x=15,y=0.13,labels=expression(bold("BRICK+SEJ")), cex=1.3, col=mycol.rgb[mycol_complete])
  
  
  ## LEGEND
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab=""
       # ,main="Global mean sea level at 2100 [m]"
  )
  
  legend(x = "top",inset = 0,
         # legend = c("Wide priors combined with data (high-temperature scenario)" , 
         #            "Wide priors combined with data and expert assessment (high-temperature scenario)",
         #            "17th-83rd percentile - SEJ2018 expert assessment (high-temperature scenario)",
         #            "5-95th percentile - SEF2018 expert assessment (high-temperature scenario)"
         #            # "IPCC AR6 likely range (RCP 8.5)"
         # ),
         
         # legend = c("BRICK" , 
         #            "BRICK + SEJ2018"
         #            # expression("5"^th~"-95"^th~" percentile, SEJ2018"),
         #            # # expression("17"^th~"-83"^rd~" percentile, SEJ2018")
         #            # "IPCC AR6 likely range (RCP 8.5)"
         # ), 

         
         legend = c("data-model calibration" , 
                    "expert-data-model calibration",
                    "wide priors",
                    "expert-informed priors"
                    # expression("5"^th~"-95"^th~" percentile, SEJ2018"),
                    # # expression("17"^th~"-83"^rd~" percentile, SEJ2018")
                    # "IPCC AR6 likely range (RCP 8.5)"
         ), 
         
         lwd=c(2,2,1,1),
         # pch=c(NA,
         #       NA
         #       # NA,
         #       # 15
         # ),
         lty=c("solid",
               "solid",
               "dashed",
               "dashed"
         ),
         bty='n', cex=mycexlab,
         col=c(mycol.rgb[mycol_standard],
               mycol.rgb[mycol_complete],
               mycol.rgb[mycol_priors],
               mycol.rgb[mycol_experts]
               # "#c77d7dff"
         ), 
         horiz = FALSE,
         
  )
  
}

##########################################################################################################



##########################################################################################################
# Figure 5 - GIS parameters

# start post rejection-sampled parameter sets
chain_priors <- readRDS("priors_0426_parameters_good.rds")
chain_experts <- readRDS("experts_0426_parameters_good.rds")
chain_standard <- readRDS("standard_0426_parameters_good.rds")
chain_complete <- readRDS("complete_0426_parameters_good.rds")

# how many to subsample to?
n.ensemble <- 2e6 # subsample to one million

##==============================================================================
## subsample

# priors
if (nrow(chain_priors)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_priors))
  sampleme <- seq(from=1, to=nrow(chain_priors),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_priors[k,]
  }
  
  chain_priors <- parameters
}

# experts
if (nrow(chain_experts)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_experts))
  sampleme <- seq(from=1, to=nrow(chain_experts),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_experts[k,]
  }
  
  chain_experts <- parameters
}

# standard
if (nrow(chain_standard)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_standard))
  sampleme <- seq(from=1, to=nrow(chain_standard),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_standard[k,]
  }
  
  chain_standard <- parameters
}

# complete
if (nrow(chain_complete)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_complete))
  sampleme <- seq(from=1, to=nrow(chain_complete),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_complete[k,]
  }
  
  chain_complete <- parameters
}

if (exists("parameters")){
  rm(parameters) 
}

##==============================================================================
##==============================================================================
## Which model(s) will you use?
## If you want to plug in your own model, insert the relevant "luse.XXX" line
## below, as well as into the "luse.brick = ..." command.
## Exactly one of te or tee must be TRUE.
luse.sneasy   = FALSE    # Simple Nonlinear EArth SYstem model (DOECLIM+CCM)
luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # thermosteric expansion contribution to SLR
luse.tee      = FALSE   # explicit thermosteric expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = TRUE    # Antarctic ice sheet model
luse.lws      = FALSE    # land water storage
luse.brick = cbind(luse.sneasy, luse.doeclim, luse.gsic, luse.te, luse.tee,
                   luse.simple, luse.dais, luse.lws)

## If you are using DAIS, include the fast dynamics emulator?
if (luse.dais) {
  l.aisfastdy <- TRUE
} else {
  l.aisfastdy <- FALSE  # force FALSE if not using DAIS
  slope.Ta2Tg <- NULL
  intercept.Ta2Tg <- NULL
}

if(luse.te & luse.tee) {
  luse.tee <- FALSE
  print('Only use 1 thermosteric expansion model; switching off explicit model.')
}
##==============================================================================
## Define parameters and their prior ranges
## -> Note: 'parnames' is defined here, which establishes how the parameters
##    are passed around into DEoptim, MCMC, likelihood functions, and the models

source('../calibration/CARL_BRICK_parameterSetup_v2.R')
##==============================================================================


## get posterior densities
chain_priors_densities <- replicate(length(parnames), vector("list", 7), simplify = FALSE)
chain_experts_densities <- chain_priors_densities
chain_standard_densities <- chain_priors_densities
chain_complete_densities <- chain_priors_densities


for (pp in 1:length(parnames)){
  chain_priors_densities[[pp]] <- density(chain_priors[,pp])
}

for (pp in 1:length(parnames)){
  chain_experts_densities[[pp]] <- density(chain_experts[,pp])
}

for (pp in 1:length(parnames)){
  chain_standard_densities[[pp]] <- density(chain_standard[,pp])
}

for (pp in 1:length(parnames)){
  chain_complete_densities[[pp]] <- density(chain_complete[,pp])
}

## get nice colors
source('../Useful/colorblindPalette.R') # Get nice plotting colors: mycol array

## mycol colors, by row:
# 1 black
# 2 dark blue
# 3 aqua
# 4 pink
# 5 light orange
# 6 purple
# 7 blue
# 8 light purple
# 9 light blue
# 10 very light blue
# 11 dark red
# 12 brown
# 13 dark orange
# 14 neon green
# 15 neon yellow

mycol_experts <- 9 # originally 5
mycol_standard <- 1
mycol_complete <- 9
mycol_priors <- 1

## draw pdfs
sequence_GIS = c(9,10,11,12,13,31,32)
sequence_GIS_select = c(9,10,11,12,13)

cexlab = 1.2
mylwd = 2

if(TRUE){
  layout(cbind(c(1,3,5,7),c(2,4,6,7)))
  par(mai=c(.6,.8,.3,.3) #c(bottom, left, top, right)
  )
  for (pp in sequence_GIS_select){
    
    myylim <- max(pmax(chain_experts_densities[[pp]]$y,chain_standard_densities[[pp]]$y,chain_complete_densities[[pp]]$y))
    
    # experts
    plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
         ylim = c(0,myylim),
         type="l",
         xlab = "", #xlab=parnames[pp],
         ylab='Probability density', lwd=1, lty="dashed", 
         col=mycol.rgb[mycol_experts],
         yaxt='n', cex.lab=cexlab)
    

    
    if(pp==9){ # a.simple
      mtext(side=3, text=expression(bold('a')), line=.25, cex=.9, adj=0);
      mtext("equilibrium volume sensitivity to",side=1,line=2, cex=.8)
      mtext(expression(paste("temperature changes (mSLE","\u00B0C"^-1,")" )),side=1,line=3.4, cex=.8)
    }
    
    if(pp==10){ # b.simple
      mtext(side=3, text=expression(bold('b')), line=.25, cex=.9, adj=0);
      mtext("equilibrium volume for zero",side=1,line=2, cex=.8)
      mtext(expression(paste("temperature anomaly (mSLE)")),side=1,line=3.4, cex=.8)
    }
    
    if(pp==11){ # alpha.simple
      mtext(side=3, text=expression(bold('c')), line=.25, cex=.9, adj=0);
      mtext("sensitivity of e-folding timescale of GIS volume response",side=1,line=2, cex=.8)
      mtext(expression(paste("to changes in temperature (","\u00B0C"^-1,"yr"^-1,")")),side=1,line=3.4, cex=.8)
    }
    
    if(pp==12){ # beta.simple
      mtext(side=3, text=expression(bold('d')), line=.25, cex=.9, adj=0);
      mtext("equilibrium timescale of GIS volume response",side=1,line=2, cex=.8)
      mtext(expression(paste("to changes in temperature (yr"^-1,")")),side=1,line=3.4, cex=.8)
    }
    
    if(pp==13){ # V0
      mtext(side=3, text=expression(bold('e')), line=.25, cex=.9, adj=0);
      mtext("initial condition for GIS contributions",side=1,line=2, cex=.8)
      mtext(expression(paste("to GMSL (mSLE)")),side=1,line=3.4, cex=.8)
    }
    
    # priors
    x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
    lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
          col="black", lwd=1, lty=2)
    
    # standard
    lines(chain_standard_densities[[pp]]$x,chain_standard_densities[[pp]]$y,lwd=mylwd,
          col=mycol.rgb[mycol_standard], 
          cex.lab=cexlab)
    
    # complete
    lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,
          col=mycol.rgb[mycol_complete],
          cex.lab=cexlab)
    
    if(pp==9){
      # 
      # text(x=-1.7,y=.39,labels=expression(bold("BRICK")), cex=1.4)
      # text(x=-0.8,y=.25,labels=expression(bold("BRICK+SEJ")), cex=1.4, col=mycol.rgb[mycol_complete])
      
    }
    
  }
  
  plot.new() # space
  
  ## LEGEND
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab=""
       # ,main="Global mean sea level at 2100 [m]"
  )
  
  legend(x = "top",inset = 0,

         legend = c("expert-data-model calibration" , 
                    "data-model calibration",
                    "expert-informed priors",
                    "standard priors"
         ), 
         
         lwd=c(2,2,1,1),
         
         lty=c("solid",
               "solid",
               "dashed",
               "dashed"
         ),
         bty='n', cex=mycexlab,
         col=c(mycol.rgb[mycol_complete],
               mycol.rgb[mycol_standard],
               mycol.rgb[mycol_experts],
               mycol.rgb[mycol_priors]
               # "#c77d7dff"
         ), 
         horiz = FALSE,
         
  )

}

##########################################################################################################



##########################################################################################################
# Supplemental figure - temperature and tcrit

# start post rejection-sampled parameter sets
chain_priors <- readRDS("priors_0426_parameters_good.rds")
chain_experts <- readRDS("experts_0426_parameters_good.rds")
chain_standard <- readRDS("standard_0426_parameters_good.rds")
chain_complete <- readRDS("complete_0426_parameters_good.rds")

# how many to subsample to?
n.ensemble <- 2e6 # subsample to one million

##==============================================================================
## subsample

# priors
if (nrow(chain_priors)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_priors))
  sampleme <- seq(from=1, to=nrow(chain_priors),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_priors[k,]
  }
  
  chain_priors <- parameters
}

# experts
if (nrow(chain_experts)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_experts))
  sampleme <- seq(from=1, to=nrow(chain_experts),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_experts[k,]
  }
  
  chain_experts <- parameters
}

# standard
if (nrow(chain_standard)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_standard))
  sampleme <- seq(from=1, to=nrow(chain_standard),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_standard[k,]
  }
  
  chain_standard <- parameters
}

# complete
if (nrow(chain_complete)>n.ensemble){
  parameters <- matrix(0, n.ensemble, ncol(chain_complete))
  sampleme <- seq(from=1, to=nrow(chain_complete),by=1)
  sample.index <- sample(sampleme, n.ensemble, replace=TRUE)
  
  for (i in 1:length(sample.index)){
    k = sample.index[i]
    parameters[i,] <- chain_complete[k,]
  }
  
  chain_complete <- parameters
}

if (exists("parameters")){
  rm(parameters) 
}

##==============================================================================
##==============================================================================
## Which model(s) will you use?
## If you want to plug in your own model, insert the relevant "luse.XXX" line
## below, as well as into the "luse.brick = ..." command.
## Exactly one of te or tee must be TRUE.
luse.sneasy   = FALSE    # Simple Nonlinear EArth SYstem model (DOECLIM+CCM)
luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # thermosteric expansion contribution to SLR
luse.tee      = FALSE   # explicit thermosteric expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = TRUE    # Antarctic ice sheet model
luse.lws      = FALSE    # land water storage
luse.brick = cbind(luse.sneasy, luse.doeclim, luse.gsic, luse.te, luse.tee,
                   luse.simple, luse.dais, luse.lws)

## If you are using DAIS, include the fast dynamics emulator?
if (luse.dais) {
  l.aisfastdy <- TRUE
} else {
  l.aisfastdy <- FALSE  # force FALSE if not using DAIS
  slope.Ta2Tg <- NULL
  intercept.Ta2Tg <- NULL
}

if(luse.te & luse.tee) {
  luse.tee <- FALSE
  print('Only use 1 thermosteric expansion model; switching off explicit model.')
}
##==============================================================================
## Define parameters and their prior ranges
## -> Note: 'parnames' is defined here, which establishes how the parameters
##    are passed around into DEoptim, MCMC, likelihood functions, and the models

source('../calibration/CARL_BRICK_parameterSetup_v2.R')
##==============================================================================


## get posterior densities
chain_priors_densities <- replicate(length(parnames), vector("list", 7), simplify = FALSE)
chain_experts_densities <- chain_priors_densities
chain_standard_densities <- chain_priors_densities
chain_complete_densities <- chain_priors_densities


for (pp in 1:length(parnames)){
  chain_priors_densities[[pp]] <- density(chain_priors[,pp])
}

for (pp in 1:length(parnames)){
  chain_experts_densities[[pp]] <- density(chain_experts[,pp])
}

for (pp in 1:length(parnames)){
  chain_standard_densities[[pp]] <- density(chain_standard[,pp])
}

for (pp in 1:length(parnames)){
  chain_complete_densities[[pp]] <- density(chain_complete[,pp])
}

## get nice colors
source('../Useful/colorblindPalette.R') # Get nice plotting colors: mycol array

## mycol colors, by row:
# 1 black
# 2 dark blue
# 3 aqua
# 4 pink
# 5 light orange
# 6 purple
# 7 blue
# 8 light purple
# 9 light blue
# 10 very light blue
# 11 dark red
# 12 brown
# 13 dark orange
# 14 neon green
# 15 neon yellow


mycol_experts <- 9 # originally 5
mycol_standard <- 1
mycol_complete <- 9
mycol_priors <- 1
mycol_tcrit <- 11
col_tcrit <- 11

cexlab = 1.2
mylwd = 2

myylim <- max(pmax(priors.density.temp.2100$y,standard.density.temp.2100$y,complete.density.temp.2100$y, experts.density.temp.2100$y))

# # priors
# plot(priors.density.temp.2100$x,priors.density.temp.2100$y,
#      ylim = c(0,myylim),
#      type="l",
#      xlab = "", #xlab=parnames[pp],
#      ylab='Probability density', lwd=1, lty="dashed", 
#      col=mycol.rgb[mycol_priors],
#      yaxt='n', cex.lab=cexlab)
# 
# # experts
# lines(experts.density.temp.2100$x,experts.density.temp.2100$y,
#       col=mycol.rgb[mycol_experts],
#       lty="dashed",
#       yaxt='n', cex.lab=cexlab)
# 
# # standard
# lines(standard.density.temp.2100$x,standard.density.temp.2100$y,lwd=mylwd,
#       col=mycol.rgb[mycol_standard], 
#       cex.lab=cexlab)
# 
# # complete
# lines(complete.density.temp.2100$x,complete.density.temp.2100$y,lwd=mylwd,
#       col=mycol.rgb[mycol_complete],
#       cex.lab=cexlab)
# 
# sequence_tcrit <- 28
# 
# for (pp in sequence_tcrit){
#   
#   myylim <- max(pmax(chain_experts_densities[[pp]]$y,chain_standard_densities[[pp]]$y,chain_complete_densities[[pp]]$y))
#   
#   # experts
#   plot(chain_experts_densities[[pp]]$x,chain_experts_densities[[pp]]$y,
#        ylim = c(0,myylim),
#        type="l",
#        xlab = "", #xlab=parnames[pp],
#        ylab='Probability density', lwd=1, lty="dashed", 
#        col=mycol.rgb[mycol_experts],
#        yaxt='n', cex.lab=cexlab)
#   
#   # priors
#   x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
#   lines(x,dunif(x, min = bound.lower[pp], max = bound.upper[pp]),
#         col="black", lwd=1, lty=2)
#   
#   # standard
#   lines(chain_standard_densities[[pp]]$x,chain_standard_densities[[pp]]$y,lwd=mylwd,
#         col=mycol.rgb[mycol_standard], 
#         cex.lab=cexlab)
#   
#   # complete
#   lines(chain_complete_densities[[pp]]$x,chain_complete_densities[[pp]]$y,lwd=mylwd,
#         col=mycol.rgb[mycol_complete],
#         cex.lab=cexlab)
#   
# }


#########################################################################
# Supplemental figure: Temperature and TCrit
# CONVERT Global mean temperature to Antarctic Ocean Temperature
anto <- function(a  =  0.26,
                 b  =  0.62,
                 Tf = -1.8,
                 Tg) {
  
  c  <- (Tf-b)/a
  Toc <-  Tf + (a*Tg + b-Tf) / (1 + exp(-Tg+c))
  
  return(Toc)
}


parameters_good_complete <- readRDS("complete_0426_parameters_good.rds")
tcrit_values <- parameters_good_complete[,28]

# Scale tcrit from AIS surface temperature to global mean temperature (from Wong et al 2017, Morice et al. 2012, shaffer 2014)
tgcrit <- tcrit_values 
tgcrit <- slope.Ta2Tg*tcrit_values + intercept.Ta2Tg
tgcrit_95 <- quantile(tgcrit, probs = c(0.95))
tgcrit_83 <- quantile(tgcrit, probs = c(0.83))
tgcrit_50 <- quantile(tgcrit, probs = c(0.5))
tgcrit_17 <- quantile(tgcrit, probs = c(0.17))
tgcrit_05 <- quantile(tgcrit, probs = c(0.05))

if(TRUE){
  
  # plot temp median
  plot(t.proj[iproj],complete.temp.rcp26.50[iproj],type='l', ann=FALSE, # median
       # col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]),
       col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]),
       lwd=2,
       xlim=c(2000,2100), ylim=c(0,6),
       xaxt='n', yaxt='n', xaxs='i', yaxs='i'
  );
  axis(1, seq(2000,2100,by=20));
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(0,6,by=1), lab=c('0','','2','','4','','6'));
  mtext(side=2, text='[deg C]', line=2.3, cex=.9);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=3, text=expression(bold(' Projected global mean temperature')), line=.25, cex=.9, adj=0);
  
  # plot tcrit median
  abline(h=tgcrit_50, 
         # lty="dashed", 
         lwd=2, 
         col=rgb(mycol[col_tcrit,1],mycol[col_tcrit,2],mycol[col_tcrit,3]))
  
  # plot tcrit 5-95 CI
  complete.tcrit.rcp26.95 <- complete.temp.rcp26.95
  complete.tcrit.rcp26.05 <- complete.temp.rcp26.05
  complete.tcrit.rcp26.95[] <-tgcrit_95
  complete.tcrit.rcp26.05[] <-tgcrit_05
  
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.tcrit.rcp26.95[iproj],rev(complete.tcrit.rcp26.05[iproj])),
          col=rgb(mycol[col_tcrit,1],mycol[col_tcrit,2],mycol[col_tcrit,3],0.3), border=NA);
  
  
  # plot temperature 5-95 CI
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.temp.rcp26.95[iproj],rev(complete.temp.rcp26.05[iproj])),
          col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3],opacity), border=NA);
  
  
  
  
  legend(x = "topleft",inset = 0,
         legend = c("global mean temperature, median" ,
                    "global mean temperature, 5-95% range",
                    "Tcrit, median",
                    "Tcrit, 5-95% range"
         ),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,2,8), bty='n', cex=.9,
         col=c(rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]),
               rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3],opacity),
               rgb(mycol[col_tcrit,1],mycol[col_tcrit,2],mycol[col_tcrit,3]),
               rgb(mycol[col_tcrit,1],mycol[col_tcrit,2],mycol[col_tcrit,3],0.3)
         ), 
         horiz = FALSE)
}



#########################################################################
# Extended Data Figure 1: Complete 

if(TRUE){
  
  n.sig = 2         # how many sigma to plot around the obs?
  # layout(cbind(c(1,4,7),c(2,5,7),c(3,6,7)))
  layout(cbind(c(1,3,5,6),c(2,4,5,6)))
  # par(mai=c(.5,.7,.2,.08)) #c(bottom, left, top, right)
  par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) # omi = c(0,0,0,2) for right side space
  plotall <- FALSE
  
  
  if(plotall){
    # >>> SURFACE TEMPERATURE <<<
    plot(mod.time[midx.temp], complete.temp.50[midx.temp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1850,2016), ylim=c(-.3,1.5), cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    # mtext(side=2, text='Surface temperature\n[deg C]', line=2.3, cex=.9);
    mtext(side=3, text=expression(bold('a) Surface temperature')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[midx.temp],rev(mod.time[midx.temp])), c(complete.temp.95,rev(complete.temp.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.temp.time[oidx.temp], obs.temp.norm[oidx.temp], type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.temp.time[oidx.temp],rev(obs.temp.time[oidx.temp])),
            c(obs.temp.norm[oidx.temp]+n.sig*obs.temp.err[oidx.temp],rev(obs.temp.norm[oidx.temp]-n.sig*obs.temp.err[oidx.temp])),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
    #legend(1839,1.65,c("5-95% range, modeled","2-sigma range, observations"),
    #       col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])), lwd=2, bty='n', cex=1.2)
    
    # >>> OCEAN HEAT <<<
    itmp=midx.ocheat[1]:nrow(complete.ocheat.hind)
    plot(mod.time[itmp], complete.ocheat.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1950,2016), ylim=c(-50,50), cex.lab=1.2, cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    # mtext(side=2, text='b) Ocean heat uptake', line=2.3, cex=.9);
    mtext(side=3, text=expression(bold('b) Ocean heat uptake')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(complete.ocheat.95[itmp],rev(complete.ocheat.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.ocheat.time, obs.ocheat.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.ocheat.time,rev(obs.ocheat.time)), c(obs.ocheat.norm+n.sig*obs.ocheat.err,rev(obs.ocheat.norm-n.sig*obs.ocheat.err)),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  }
  # >>> GSIC <<<
  itmp=midx.gsic[1]:nrow(complete.gsic.hind)
  plot(mod.time[itmp], complete.gsic.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.01,.04), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Glaciers and\nsmall ice caps [m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  # mtext(side=3, text=expression(bold(' c')), line=.25, cex=.9, adj=0);
  mtext(side=3, text=expression(bold(' a')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+21, .85 * par('usr')[4], labels = 'Glaciers and ice caps', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(complete.gsic.95[itmp],rev(complete.gsic.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gsic.time, obs.gsic.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  # polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  # >>> GIS <<<
  plot(mod.time, complete.gis.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+6, .85 * par('usr')[4], labels = 'GIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(complete.gis.95,rev(complete.gis.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gis.time, obs.gis.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  # polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  # >>> TE <<<
  x1971=seq(trends.te[1,4],trends.te[1,5])
  c1971=mean(x1971); yc1971=complete.te.50[which(mod.time==mean(x1971))]
  lo1971 = yc1971+(trends.te[1,2]/1000)*(x1971-c1971)
  hi1971 = yc1971+(trends.te[1,3]/1000)*(x1971-c1971)
  y1971 = yc1971+(trends.te[1,1]/1000)*(x1971-c1971)
  
  x1993=seq(trends.te[2,4],trends.te[2,5])
  c1993=mean(x1993); yc1993=complete.te.50[which(mod.time==mean(x1993))]
  lo1993 = yc1993+(trends.te[2,2]/1000)*(x1993-c1993)
  hi1993 = yc1993+(trends.te[2,3]/1000)*(x1993-c1993)
  y1993 = yc1993+(trends.te[2,1]/1000)*(x1993-c1993)
  
  plot(mod.time, complete.te.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.04,.06), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Thermal expansion\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' c')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+19, .8 * par('usr')[4], labels = 'Thermal expansion', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(complete.te.95,rev(complete.te.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(x1971,y1971, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1971,y1971, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  polygon(c(x1971,rev(x1971)), c(lo1971,rev(hi1971)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  lines(x1993,y1993, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1993,y1993, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(x1993,rev(x1993)), c(lo1993,rev(hi1993)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  
  # >>> TOTAL SLR <<<
  plot(mod.time, complete.slr.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+21, .75 * par('usr')[4], labels = 'GMSL', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(complete.slr.95,rev(complete.slr.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  
  # >>> AIS PALEO, SMOOTHED <<<
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], complete.ais.paleo.50[ipaleo], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=.25, cex=.9, adj=0);
  text(.95 * par('usr')[1], .8 * par('usr')[4], labels = 'AIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(complete.ais.paleo.95[ipaleo],rev(complete.ais.paleo.05[ipaleo])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # legend(-92000,12, c("5-95% range, model" , "2-sigma range, observations"), lwd=2, bty='n', cex=1.2,
  #        col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]) , rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])) )
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top",inset = 0,
         legend = c("median, model" ,
                    "5-95% range, model",
                    "observations"#,
                    # "2-sigma range, observations"
         ),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,NA,8), bty='n', cex=1.2,
         pch=c(NA,NA,20,NA),
         col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
               rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
               rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
         ), 
         horiz = FALSE)
  
}

#########################################################################
# Extended Data Figure 2: Standard 

if(TRUE){
  
  n.sig = 2         # how many sigma to plot around the obs?
  # layout(cbind(c(1,4,7),c(2,5,7),c(3,6,7)))
  layout(cbind(c(1,3,5,6),c(2,4,5,6)))
  # par(mai=c(.5,.7,.2,.08)) #c(bottom, left, top, right)
  par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) # omi = c(0,0,0,2) for right side space
  plotall <- FALSE
  
  
  if(plotall){
    # >>> SURFACE TEMPERATURE <<<
    plot(mod.time[midx.temp], standard.temp.50[midx.temp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1850,2016), ylim=c(-.3,1.5), cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    # mtext(side=2, text='Surface temperature\n[deg C]', line=2.3, cex=.9);
    mtext(side=3, text=expression(bold('a) Surface temperature')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[midx.temp],rev(mod.time[midx.temp])), c(standard.temp.95,rev(standard.temp.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.temp.time[oidx.temp], obs.temp.norm[oidx.temp], type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.temp.time[oidx.temp],rev(obs.temp.time[oidx.temp])),
            c(obs.temp.norm[oidx.temp]+n.sig*obs.temp.err[oidx.temp],rev(obs.temp.norm[oidx.temp]-n.sig*obs.temp.err[oidx.temp])),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
    #legend(1839,1.65,c("5-95% range, modeled","2-sigma range, observations"),
    #       col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])), lwd=2, bty='n', cex=1.2)
    
    # >>> OCEAN HEAT <<<
    itmp=midx.ocheat[1]:nrow(standard.ocheat.hind)
    plot(mod.time[itmp], standard.ocheat.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1950,2016), ylim=c(-50,50), cex.lab=1.2, cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    # mtext(side=2, text='b) Ocean heat uptake', line=2.3, cex=.9);
    mtext(side=3, text=expression(bold('b) Ocean heat uptake')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(standard.ocheat.95[itmp],rev(standard.ocheat.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.ocheat.time, obs.ocheat.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.ocheat.time,rev(obs.ocheat.time)), c(obs.ocheat.norm+n.sig*obs.ocheat.err,rev(obs.ocheat.norm-n.sig*obs.ocheat.err)),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  }
  # >>> GSIC <<<
  itmp=midx.gsic[1]:nrow(standard.gsic.hind)
  plot(mod.time[itmp], standard.gsic.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.01,.04), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Glaciers and\nsmall ice caps [m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  # mtext(side=3, text=expression(bold(' c')), line=.25, cex=.9, adj=0);
  mtext(side=3, text=expression(bold(' a')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+21, .85 * par('usr')[4], labels = 'Glaciers and ice caps', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(standard.gsic.95[itmp],rev(standard.gsic.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gsic.time, obs.gsic.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  # polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  # >>> GIS <<<
  plot(mod.time, standard.gis.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+6, .85 * par('usr')[4], labels = 'GIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(standard.gis.95,rev(standard.gis.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gis.time, obs.gis.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  # polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  # >>> TE <<<
  x1971=seq(trends.te[1,4],trends.te[1,5])
  c1971=mean(x1971); yc1971=standard.te.50[which(mod.time==mean(x1971))]
  lo1971 = yc1971+(trends.te[1,2]/1000)*(x1971-c1971)
  hi1971 = yc1971+(trends.te[1,3]/1000)*(x1971-c1971)
  y1971 = yc1971+(trends.te[1,1]/1000)*(x1971-c1971)
  
  x1993=seq(trends.te[2,4],trends.te[2,5])
  c1993=mean(x1993); yc1993=standard.te.50[which(mod.time==mean(x1993))]
  lo1993 = yc1993+(trends.te[2,2]/1000)*(x1993-c1993)
  hi1993 = yc1993+(trends.te[2,3]/1000)*(x1993-c1993)
  y1993 = yc1993+(trends.te[2,1]/1000)*(x1993-c1993)
  
  plot(mod.time, standard.te.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.04,.06), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Thermal expansion\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' c')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+19, .8 * par('usr')[4], labels = 'Thermal expansion', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(standard.te.95,rev(standard.te.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(x1971,y1971, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1971,y1971, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  polygon(c(x1971,rev(x1971)), c(lo1971,rev(hi1971)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  lines(x1993,y1993, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1993,y1993, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(x1993,rev(x1993)), c(lo1993,rev(hi1993)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  
  # >>> TOTAL SLR <<<
  plot(mod.time, standard.slr.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+21, .75 * par('usr')[4], labels = 'GMSL', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(standard.slr.95,rev(standard.slr.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  
  # >>> AIS PALEO, SMOOTHED <<<
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], standard.ais.paleo.50[ipaleo], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=.25, cex=.9, adj=0);
  text(.95 * par('usr')[1], .8 * par('usr')[4], labels = 'AIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(standard.ais.paleo.95[ipaleo],rev(standard.ais.paleo.05[ipaleo])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # legend(-92000,12, c("5-95% range, model" , "2-sigma range, observations"), lwd=2, bty='n', cex=1.2,
  #        col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]) , rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])) )
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top",inset = 0,
         legend = c("median, model" ,
                    "5-95% range, model",
                    "observations"#,
                    # "2-sigma range, observations"
         ),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,NA,8), bty='n', cex=1.2,
         pch=c(NA,NA,20,NA),
         col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
               rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
               rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
         ), 
         horiz = FALSE)
  
}

#########################################################################
# Extended Data Figure 2: Experts 

if(TRUE){
  
  n.sig = 2         # how many sigma to plot around the obs?
  # layout(cbind(c(1,4,7),c(2,5,7),c(3,6,7)))
  layout(cbind(c(1,3,5,6),c(2,4,5,6)))
  # par(mai=c(.5,.7,.2,.08)) #c(bottom, left, top, right)
  par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) # omi = c(0,0,0,2) for right side space
  plotall <- FALSE
  
  
  if(plotall){
    # >>> SURFACE TEMPERATURE <<<
    plot(mod.time[midx.temp], experts.temp.50[midx.temp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1850,2016), ylim=c(-.3,1.5), cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    # mtext(side=2, text='Surface temperature\n[deg C]', line=2.3, cex=.9);
    mtext(side=3, text=expression(bold('a) Surface temperature')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[midx.temp],rev(mod.time[midx.temp])), c(experts.temp.95,rev(experts.temp.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.temp.time[oidx.temp], obs.temp.norm[oidx.temp], type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.temp.time[oidx.temp],rev(obs.temp.time[oidx.temp])),
            c(obs.temp.norm[oidx.temp]+n.sig*obs.temp.err[oidx.temp],rev(obs.temp.norm[oidx.temp]-n.sig*obs.temp.err[oidx.temp])),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
    #legend(1839,1.65,c("5-95% range, modeled","2-sigma range, observations"),
    #       col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])), lwd=2, bty='n', cex=1.2)
    
    # >>> OCEAN HEAT <<<
    itmp=midx.ocheat[1]:nrow(experts.ocheat.hind)
    plot(mod.time[itmp], experts.ocheat.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1950,2016), ylim=c(-50,50), cex.lab=1.2, cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    # mtext(side=2, text='b) Ocean heat uptake', line=2.3, cex=.9);
    mtext(side=3, text=expression(bold('b) Ocean heat uptake')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(experts.ocheat.95[itmp],rev(experts.ocheat.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.ocheat.time, obs.ocheat.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.ocheat.time,rev(obs.ocheat.time)), c(obs.ocheat.norm+n.sig*obs.ocheat.err,rev(obs.ocheat.norm-n.sig*obs.ocheat.err)),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  }
  # >>> GSIC <<<
  itmp=midx.gsic[1]:nrow(experts.gsic.hind)
  plot(mod.time[itmp], experts.gsic.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.01,.28), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Glaciers and\nsmall ice caps [m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  # mtext(side=3, text=expression(bold(' c')), line=.25, cex=.9, adj=0);
  mtext(side=3, text=expression(bold(' a')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+21, .85 * par('usr')[4], labels = 'Glaciers and ice caps', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(experts.gsic.95[itmp],rev(experts.gsic.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gsic.time, obs.gsic.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  # polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  # >>> GIS <<<
  plot(mod.time, experts.gis.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+6, .85 * par('usr')[4], labels = 'GIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(experts.gis.95,rev(experts.gis.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gis.time, obs.gis.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  # polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  # >>> TE <<<
  x1971=seq(trends.te[1,4],trends.te[1,5])
  c1971=mean(x1971); yc1971=experts.te.50[which(mod.time==mean(x1971))]
  lo1971 = yc1971+(trends.te[1,2]/1000)*(x1971-c1971)
  hi1971 = yc1971+(trends.te[1,3]/1000)*(x1971-c1971)
  y1971 = yc1971+(trends.te[1,1]/1000)*(x1971-c1971)
  
  x1993=seq(trends.te[2,4],trends.te[2,5])
  c1993=mean(x1993); yc1993=experts.te.50[which(mod.time==mean(x1993))]
  lo1993 = yc1993+(trends.te[2,2]/1000)*(x1993-c1993)
  hi1993 = yc1993+(trends.te[2,3]/1000)*(x1993-c1993)
  y1993 = yc1993+(trends.te[2,1]/1000)*(x1993-c1993)
  
  # load TE observations
  # x1971 <- readRDS("TE_x1971.rds")
  # c1971 <- readRDS("TE_c1971.rds")
  # lo1971<- readRDS("TE_lo1971.rds")
  # hi1971<- readRDS("TE_hi1971.rds")
  # y1971 <- readRDS("TE_y1971.rds")
  # 
  # x1993 <- readRDS("TE_x1993.rds")
  # c1993 <- readRDS("TE_c1993.rds")
  # lo1993 <- readRDS("TE_lo1993.rds")
  # hi1993 <- readRDS("TE_hi1993.rds")
  # y1993 <- readRDS("TE_y1993.rds")
  
  plot(mod.time, experts.te.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.15,.16), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Thermal expansion\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+19, .8 * par('usr')[4], labels = 'Thermal expansion', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(experts.te.95,rev(experts.te.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(x1971,y1971, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1971,y1971, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  polygon(c(x1971,rev(x1971)), c(lo1971,rev(hi1971)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  lines(x1993,y1993, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1993,y1993, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(x1993,rev(x1993)), c(lo1993,rev(hi1993)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  
  # >>> TOTAL SLR <<<
  plot(mod.time, experts.slr.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+21, .75 * par('usr')[4], labels = 'GMSL', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(experts.slr.95,rev(experts.slr.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  
  # >>> AIS PALEO, SMOOTHED <<<
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], experts.ais.paleo.50[ipaleo], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=.25, cex=.9, adj=0);
  text(.95 * par('usr')[1], .8 * par('usr')[4], labels = 'AIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(experts.ais.paleo.95[ipaleo],rev(experts.ais.paleo.05[ipaleo])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # legend(-92000,12, c("5-95% range, model" , "2-sigma range, observations"), lwd=2, bty='n', cex=1.2,
  #        col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]) , rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])) )
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top",inset = 0,
         legend = c("median, model" ,
                    "5-95% range, model",
                    "observations"#,
                    # "2-sigma range, observations"
         ),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,NA,8), bty='n', cex=1.2,
         pch=c(NA,NA,20,NA),
         col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
               rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
               rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
         ), 
         horiz = FALSE)
  
}



#########################################################################
# Extended Data Figure 2: Priors

if(TRUE){
  
  n.sig = 2         # how many sigma to plot around the obs?
  # layout(cbind(c(1,4,7),c(2,5,7),c(3,6,7)))
  layout(cbind(c(1,3,5,6),c(2,4,5,6)))
  # par(mai=c(.5,.7,.2,.08)) #c(bottom, left, top, right)
  par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) # omi = c(0,0,0,2) for right side space
  plotall <- FALSE
  
  
  if(plotall){
    # >>> SURFACE TEMPERATURE <<<
    plot(mod.time[midx.temp], priors.temp.50[midx.temp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1850,2016), ylim=c(-.3,1.5), cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    # mtext(side=2, text='Surface temperature\n[deg C]', line=2.3, cex=.9);
    mtext(side=3, text=expression(bold('a) Surface temperature')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[midx.temp],rev(mod.time[midx.temp])), c(priors.temp.95,rev(priors.temp.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.temp.time[oidx.temp], obs.temp.norm[oidx.temp], type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.temp.time[oidx.temp],rev(obs.temp.time[oidx.temp])),
            c(obs.temp.norm[oidx.temp]+n.sig*obs.temp.err[oidx.temp],rev(obs.temp.norm[oidx.temp]-n.sig*obs.temp.err[oidx.temp])),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
    #legend(1839,1.65,c("5-95% range, modeled","2-sigma range, observations"),
    #       col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])), lwd=2, bty='n', cex=1.2)
    
    # >>> OCEAN HEAT <<<
    itmp=midx.ocheat[1]:nrow(priors.ocheat.hind)
    plot(mod.time[itmp], priors.ocheat.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1950,2016), ylim=c(-50,50), cex.lab=1.2, cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    # mtext(side=2, text='b) Ocean heat uptake', line=2.3, cex=.9);
    mtext(side=3, text=expression(bold('b) Ocean heat uptake')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(priors.ocheat.95[itmp],rev(priors.ocheat.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.ocheat.time, obs.ocheat.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.ocheat.time,rev(obs.ocheat.time)), c(obs.ocheat.norm+n.sig*obs.ocheat.err,rev(obs.ocheat.norm-n.sig*obs.ocheat.err)),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  }
  # >>> GSIC <<<
  itmp=midx.gsic[1]:nrow(priors.gsic.hind)
  plot(mod.time[itmp], priors.gsic.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.01,.28), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Glaciers and\nsmall ice caps [m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  # mtext(side=3, text=expression(bold(' c')), line=.25, cex=.9, adj=0);
  mtext(side=3, text=expression(bold(' a')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+21, .85 * par('usr')[4], labels = 'Glaciers and ice caps', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(priors.gsic.95[itmp],rev(priors.gsic.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gsic.time, obs.gsic.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  # polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  # >>> GIS <<<
  plot(mod.time, priors.gis.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+6, .85 * par('usr')[4], labels = 'GIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(priors.gis.95,rev(priors.gis.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gis.time, obs.gis.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
  #         col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  # polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  # >>> TE <<<
  x1971=seq(trends.te[1,4],trends.te[1,5])
  c1971=mean(x1971); yc1971=priors.te.50[which(mod.time==mean(x1971))]
  lo1971 = yc1971+(trends.te[1,2]/1000)*(x1971-c1971)
  hi1971 = yc1971+(trends.te[1,3]/1000)*(x1971-c1971)
  y1971 = yc1971+(trends.te[1,1]/1000)*(x1971-c1971)
  
  x1993=seq(trends.te[2,4],trends.te[2,5])
  c1993=mean(x1993); yc1993=priors.te.50[which(mod.time==mean(x1993))]
  lo1993 = yc1993+(trends.te[2,2]/1000)*(x1993-c1993)
  hi1993 = yc1993+(trends.te[2,3]/1000)*(x1993-c1993)
  y1993 = yc1993+(trends.te[2,1]/1000)*(x1993-c1993)
  
  # load TE observations
  # x1971 <- readRDS("TE_x1971.rds")
  # c1971 <- readRDS("TE_c1971.rds")
  # lo1971<- readRDS("TE_lo1971.rds")
  # hi1971<- readRDS("TE_hi1971.rds")
  # y1971 <- readRDS("TE_y1971.rds")
  # 
  # x1993 <- readRDS("TE_x1993.rds")
  # c1993 <- readRDS("TE_c1993.rds")
  # lo1993 <- readRDS("TE_lo1993.rds")
  # hi1993 <- readRDS("TE_hi1993.rds")
  # y1993 <- readRDS("TE_y1993.rds")
  
  plot(mod.time, priors.te.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.15,.16), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Thermal expansion\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+19, .8 * par('usr')[4], labels = 'Thermal expansion', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(priors.te.95,rev(priors.te.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(x1971,y1971, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1971,y1971, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  polygon(c(x1971,rev(x1971)), c(lo1971,rev(hi1971)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  lines(x1993,y1993, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1993,y1993, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(x1993,rev(x1993)), c(lo1993,rev(hi1993)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  
  # >>> TOTAL SLR <<<
  plot(mod.time, priors.slr.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d')), line=.25, cex=.9, adj=0);
  text(1 * par('usr')[1]+21, .75 * par('usr')[4], labels = 'GMSL', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(mod.time,rev(mod.time)), c(priors.slr.95,rev(priors.slr.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
  #         col=NA, 
  #         border=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5),
  #         lwd=2,
  #         lty="dashed");
  
  
  # >>> AIS PALEO, SMOOTHED <<<
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], priors.ais.paleo.50[ipaleo], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=.25, cex=.9, adj=0);
  text(.95 * par('usr')[1], .8 * par('usr')[4], labels = 'AIS', cex=1.4 ) # parusr xmin, xmax, ymin, ymax
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(priors.ais.paleo.95[ipaleo],rev(priors.ais.paleo.05[ipaleo])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  for (i in 1:3) {
    polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
            c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
  }
  i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  # legend(-92000,12, c("5-95% range, model" , "2-sigma range, observations"), lwd=2, bty='n', cex=1.2,
  #        col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]) , rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])) )
  
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "top",inset = 0,
         legend = c("median, model" ,
                    "5-95% range, model",
                    "observations"#,
                    # "2-sigma range, observations"
         ),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,NA,8), bty='n', cex=1.2,
         pch=c(NA,NA,20,NA),
         col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
               rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
               rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
               # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
         ), 
         horiz = FALSE)
  
}


############################################################################
# Extended Data Fig. 5 Time Series of AIS and GIS projections

if(TRUE){
# layout(cbind(c(1,3,5),c(2,4,5))) # for standard and complete
layout(cbind(c(1,2,3)))
# par(mai=c(.5,.7,.2,.08)) #c(bottom, left, top, right)
par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) # omi = c(0,0,0,2) for right side space

# AIS
# plot complete projections
plot(t.proj[iproj],complete.ais.rcp26.50[iproj],type='l', ann=FALSE, # median
     col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]),lwd=2,
     xlim=c(2000,2100), ylim=c(0,1), 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i'
);
axis(1, seq(2000,2100,by=20)); 
# axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
axis(2, seq(0,1,by=0.25), lab=c('0','','0.5','','1'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Year', line=2.3, cex=.9);
mtext(side=3, text=expression(bold('a')), line=.25, cex=.9, adj=0);
text(1* par('usr')[1]+3, .9 * par('usr')[4], labels = 'AIS', cex=1.1 )
polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.ais.rcp26.95[iproj],rev(complete.ais.rcp26.05[iproj])),
        col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
abline(h=0,col="black",lty=2)

# GIS
# plot complete projections
plot(t.proj[iproj],complete.gis.rcp26.50[iproj],type='l', ann=FALSE, # median
     col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3]),lwd=2,
     xlim=c(2000,2100), ylim=c(0,2), 
     xaxt='n', yaxt='n', xaxs='i', yaxs='i'
);
axis(1, seq(2000,2100,by=20)); 
# axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
axis(2, seq(0,2,by=0.5), lab=c('0','','1','','2'));
mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
mtext(side=1, text='Year', line=2.3, cex=.9);
mtext(side=3, text=expression(bold('b')), line=.25, cex=.9, adj=0);
text(1* par('usr')[1]+3, .9 * par('usr')[4], labels = 'GIS', cex=1.1 )
polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.gis.rcp26.95[iproj],rev(complete.gis.rcp26.05[iproj])),
        col=rgb(mycol[col_complete,1],mycol[col_complete,2],mycol[col_complete,3],opacity), border=NA);
abline(h=0,col="black",lty=2)
# 
# # GIS
# # plot standard projections
# plot(t.proj[iproj],standard.gis.rcp26.50[iproj],type='l', ann=FALSE, # median
#      col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]),lwd=2,
#      xlim=c(2000,2100), ylim=c(0,2), 
#      xaxt='n', yaxt='n', xaxs='i', yaxs='i'
# );
# axis(1, seq(2000,2100,by=20)); 
# # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
# axis(2, seq(0,2,by=0.5), lab=c('0','','1','','2'));
# mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
# mtext(side=1, text='Year', line=2.3, cex=.9);
# mtext(side=3, text=expression(bold('b')), line=.25, cex=.9, adj=0);
# text(1* par('usr')[1]+11, .9 * par('usr')[4], labels = 'AIS (BRICK+SEJ)', cex=1.1 )
# polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(standard.gis.rcp26.95[iproj],rev(standard.gis.rcp26.05[iproj])),
#         col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3],opacity), border=NA);
# abline(h=0,col="black",lty=2)



# # AIS
# # plot standard projections
# plot(t.proj[iproj],standard.ais.rcp26.50[iproj],type='l', ann=FALSE, # median
#      col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3]),lwd=2,
#      xlim=c(2000,2100), ylim=c(0,1), 
#      xaxt='n', yaxt='n', xaxs='i', yaxs='i'
# );
# axis(1, seq(2000,2100,by=20)); 
# # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
# axis(2, seq(0,1,by=0.25), lab=c('0','','0.5','','1'));
# mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
# mtext(side=1, text='Year', line=2.3, cex=.9);
# mtext(side=3, text=expression(bold('d')), line=.25, cex=.9, adj=0);
# text(1* par('usr')[1]+11, .9 * par('usr')[4], labels = 'AIS (BRICK+SEJ)', cex=1.1 )
# polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(standard.ais.rcp26.95[iproj],rev(standard.ais.rcp26.05[iproj])),
#         col=rgb(mycol[col_standard,1],mycol[col_standard,2],mycol[col_standard,3],opacity), border=NA);
# abline(h=0,col="black",lty=2)



# LEGEND
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("median, model" ,
                  "5-95% range, model"
                  # "observations"#,
                  # "Probability density"
                  # "2-sigma range, observations"
       ),
       # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
       lwd=c(2,8#,NA#,
             #8
       ), bty='n', cex=1.2,
       pch=c(NA,NA,20),
       col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
             rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5)#,
             #rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])#,
             # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])
             # rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
       ), 
       horiz = FALSE)

}


############################################################################
# Extended Data Fig. 5 v2 Time Series of AIS and GIS projections

if(TRUE){

  ## mycol colors, by row:
  # 1 black
  # 2 dark blue
  # 3 aqua
  # 4 pink
  # 5 light orange
  # 6 purple
  # 7 blue
  # 8 light purple
  # 9 light blue
  # 10 very light blue
  # 11 dark red
  # 12 brown
  # 13 dark orange
  # 14 neon green
  # 15 neon yellow
  

  col_GIS <- 9
  col_AIS <- 13

  
  # GIS
  # plot complete projections
  plot(t.proj[iproj],complete.gis.rcp26.50[iproj],type='l', ann=FALSE, # median
       col=rgb(mycol[col_GIS,1],mycol[col_GIS,2],mycol[col_GIS,3]),lwd=2,
       xlim=c(2000,2100), ylim=c(0,2), 
       xaxt='n', yaxt='n', xaxs='i', yaxs='i'
  );
  axis(1, seq(2000,2100,by=20)); 
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(0,2,by=0.5), lab=c('0','','1','','2'));
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=3, text=expression(bold('b')), line=.25, cex=.9, adj=0);
  # text(1* par('usr')[1]+3, .9 * par('usr')[4], labels = 'GIS', cex=1.1 )
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.gis.rcp26.95[iproj],rev(complete.gis.rcp26.05[iproj])),
          col=rgb(mycol[col_GIS,1],mycol[col_GIS,2],mycol[col_GIS,3],opacity), border=NA);
  abline(h=0,col="black",lty=2)
  
  # AIS
  # plot complete projections
  lines(t.proj[iproj],complete.ais.rcp26.50[iproj],type='l', ann=FALSE, # median
       col=rgb(mycol[col_AIS,1],mycol[col_AIS,2],mycol[col_AIS,3]),lwd=2
       
  );
  axis(1, seq(2000,2100,by=20)); 
  # axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'));
  axis(2, seq(0,1,by=0.25), lab=c('0','','0.5','','1'));
  
  # text(1* par('usr')[1]+3, .9 * par('usr')[4], labels = 'AIS', cex=1.1 )
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(complete.ais.rcp26.95[iproj],rev(complete.ais.rcp26.05[iproj])),
          col=rgb(mycol[col_AIS,1],mycol[col_AIS,2],mycol[col_AIS,3],opacity), border=NA);
  abline(h=0,col="black",lty=2)
  
  
  
  
  # LEGEND
  # plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "topleft",inset = 0,
         legend = c("median, GIS" ,
                    "5-95% range, GIS",
                    "median, AIS",
                    "5-95% range, AIS"
                    # "observations"#,
                    # "Probability density"
                    # "2-sigma range, observations"
         ),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,2,8
         ), bty='n', cex=.9,
         # pch=c(NA,NA,20),
         col=c(rgb(mycol[col_GIS,1],mycol[col_GIS,2],mycol[col_GIS,3]),
               rgb(mycol[col_GIS,1],mycol[col_GIS,2],mycol[col_GIS,3],opacity),
               rgb(mycol[col_AIS,1],mycol[col_AIS,2],mycol[col_AIS,3]),
               rgb(mycol[col_AIS,1],mycol[col_AIS,2],mycol[col_AIS,3],opacity)

         ), 
         horiz = FALSE)
  
}



## Save workspace image
t.end <- proc.time()
time.elapsed <- t.end - t.beg
# save.image(file=filename.saveprogress)
print("line 718")
print(time.elapsed)
print("TEDDY IS A CUTE PUPPER")
