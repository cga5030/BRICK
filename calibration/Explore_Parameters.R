# Explore Parameters <- maybe condense into MAP code later
rm(list=ls())                        # Clear all previous variables

t.beg <- proc.time() #start timer

# configure these
experts=TRUE # do you wish to invert expert assessment?
alldata=TRUE # do you wish to invert paleo and instrumental data?
configure <- '_EXPLORE_0715_' # configure filename for output files
deoptim.out <- readRDS("CARL_MAP_MAP_0713_2.4e3_13Jul2021.rds",refhook = NULL) # select MAP
amcmc.out <- readRDS("CARL_calib_amcmc_out_Complete_0614_2e7_14Jun2021.rds",refhook = NULL)
niter <- 100 # number of samples across prior range

## Set up a filename for saving RData images along the way
today <- Sys.Date(); today <- format(today,format="%d%b%Y")
filename.saveprogress <- paste('CARL_explore_parameters',configure,today,'.RData',sep='')

## Use the AIS fast dynamics emulator (Wong et al., 2017)
l.aisfastdy <- TRUE
fd.priors <- 'gamma'			# which prior distirbutions to use?

## Set the seed (for reproducibility)
set.seed(1234)
#set.seed(as.double(Sys.time())) # should yield same distributions... (a good test!)

## Needed libraries
library(DEoptim)
# library(ncdf4)
library(compiler)
library(invgamma)
enableJIT(3)

## Do you want to use RCP8.5 to make projections? (l.project=TRUE)
## Or do you want to use historical data to make hindcasts? (l.project=FALSE)
## Note -- l.project=FALSE => Kriegler (2005) data, same as Urban and Keller (2010)
l.project = TRUE
#begyear = 1765  # SNEASY start date
begyear = 1850  # DOECLIM start date
endyear = 2100
tstep   = 1
mod.time= seq(from=begyear, to=endyear, by=tstep)
begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)

## Source the models
source('../fortran/R/sneasyF.R')        # the SNEASY model (includes DOECLIM, and carbon cycle)
source('../fortran/R/doeclimF.R')       # the DOECLIM model
source('../fortran/R/GSIC_magiccF.R')   # the GSIC model
source('../fortran/R/brick_te_F.R')     # TE (thermosteric expansion) model
source('../fortran/R/brick_tee_F.R')    # TEE (explicit thermosteric expansion) model
source('../fortran/R/simpleF.R')        # GIS (Greenland Ice Sheet) model
source('../fortran/R/daisanto_fastdynF.R') # DAIS (Antarctic Ice Sheet) model
source('../R/brick_lws.R')              # LWS (land water storage)

## Source some useful functions for manipulating data
source('../R/forcing_total_carl.R')          # function to add up the total forcing
# source('../R/forcing_total.R')
source('../R/compute_indices.R')        # function to determine the model and
# data indices for comparisons

##==============================================================================
## BRICK:

## Read the model calibration data sets

## TODO
## TODO -- revise SNEASY_readData.R to match other components
## TODO

#source('../calibration/SNEASY_readData.R')    # read SNEASY calibration data
source('../calibration/DOECLIM_readData.R')   # read DOECLIM calibration data
source('../calibration/GSIC_readData.R')      # read GSIC calibration data
source('../calibration/TE_readData.R')        # read TE data
source('../calibration/SIMPLE_readData.R')    # GIS data, and trends in mass balance
#source('../calibration/DAIS_readData.R')     # DAIS forcing data (if at all uncoupled)


## TODO
## TODO -- add SNEASY calibration data to midx, oidx, obs, obs.err,
## TODO -- ind.norm.data, and i0
## TODO

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

## Set the indices of the initial condition for each sub-model
i0 = vector("list",nrow(ind.norm.data)); names(i0)=as.character(ind.norm.data[,1])

## GSIC initial conditions are actually relative to 1990 (Wigley and Raper 2005)
## Re-set these. The simulation is relative to 1990, but results and comparison
## to data is relative to 1960.
i0$gsic = which(mod.time==1990)

## GIS initial conditions are relative to 1961-1990
i0$gis = which(mod.time==1961)

## DAIS: 

## Read the data forcing for hindcasts and projections. Yields:
##  Ta (Antarctic temperature reduced to sea-level)
##  Toc (High latitude subsurface ocean temperature)
##  SL, obs.sl (Reconstructed sea level, Church and White (2011) modern-era sea level)
##  dSL (time rate of change of sea level)
##  date (240 Kyr before present to 2100 AD at 1-year intervals of forcings)

source('../calibration/DAIS_readData.R')

## For prior distributions:
alpha.var <- 2     # alpha parameter for inverse gamma for var (E[x]=beta/(alpha+1))
beta.var <- 1      # beta parameter for inverse gamma for var (uncertainty parameter)
# note that the upper bound on var.dais is not actually imposed; just for plotting
shape.lambda <- 8.1              # gives 5% quantile at lambda=0.005 and
rate.lambda <- 100*shape.lambda  # gives mean at 0.01 m/yr, DeConto and Pollard (2016)
rate.Tcrit <- 1.37               # gives 5% quantile at Tcrit = -10 deg C
shape.Tcrit <- 15*rate.Tcrit     # gives mean at -15 deg C (negative requires multiplication of Tcrit by -1)

source('../fortran/R/daisanto_fastdynF.R')

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
## Get the forcing data

if(luse.sneasy) {
  
  ## SNEASY
  setup.sneasy() # call this to initialize SNEASY
  
  ## TODO
  ## TODO -- (YG, TW) likely will need modified, to make sure forcing is as we want it
  ## TODO
  
  forcing <- vector('list',3)
  names(forcing) <- c('co2','aero','other')
  
  # load emissions data time series.
  rcp8.5.emis = read.csv("../data/sneasy_tmp/RCP85_EMISSIONS.csv")
  ibeg = which(rcp8.5.emis[,1]==begyear)
  iend = which(rcp8.5.emis[,1]==endyear)
  if(length(iend)*length(ibeg)==0) {print('ERROR - emissions data does not span mod.time')}
  emis = rcp8.5.emis[ibeg:iend,2] + rcp8.5.emis[ibeg:iend,3] # fossil + land use
  emisdata = data.frame(cbind(mod.time, emis))
  colnames(emisdata) = c("year","co2")
  forcing$co2 = emisdata$co2
  
  forcing.dat = read.table("../data/sneasy_tmp/forcing_rcp85.txt", header=TRUE)
  ibeg = which(forcing.dat[,1]==begyear)
  iend = which(forcing.dat[,1]==endyear)
  if(length(iend)*length(ibeg)==0) {print('ERROR - other radiative forcing data does not span mod.time')}
  aero.tmp  <- forcing.dat$aerosol.direct + forcing.dat$aerosol.indirect # aerosol radiative forcing
  other.tmp <- forcing.dat$ghg.nonco2 + forcing.dat$solar + forcing.dat$volcanic + forcing.dat$other
  forcing$aero  <- aero.tmp[ibeg:iend]
  forcing$other <- other.tmp[ibeg:iend]
  
} else {
  
  if(l.project) {
    forcing = read.csv( '../data/forcing_high6000.csv', header=TRUE )
    # forcing = read.csv( '../data/forcing_rcp85.csv', header=TRUE )
  } else {
    forcing = read.csv( '../data/forcing_hindcast.csv', header=TRUE )
  }
  
}

##==============================================================================
## Define the coupled model
## -> need it defined before DEoptim, so we can calculate the objective function

source('../R/CARL_BRICK_coupledModel.R')
# source('../R/BRICK_coupledModel.R')


##==============================================================================
## DAIS: Set up (pre-)calibration windows around the data
## The point: We run the model at many parameter values and see which ones
##            send the simulation through a window (width determined by the
##            observational errors) around the observational data.
## These windows are presented in Shaffer (2014) and Shepherd et al. (2012)
## 1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
## We want the cumulative sea-level equivalent in meters for the year 2002
## Note the conversion: 360Gt = 1mm SLE
## A fifth window is added to match IPCC AR5 Ch13 (page 1151) AIS SLR trend:
## 0.27 +/- 0.11 mm/year (convert to m/year here)

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

##==============================================================================

##==============================================================================
## Establish a gamma prior for drawing 1/tau.
## gamma distribution is the conjugate prior for the uncertain parameter beta
## in exponentially distributed random variable V(t):
##    dV/dt = -beta*V => V(t) = exp(-beta*t)
## For the TE model, this is not quite accurate, but similar enough (and
## previous testing with uniform priors on 1/tau confirm) that the gamma is an
## excellent approximation of the distribution of 1/tau (or tau)

if(luse.te) {
  # Original code did optimization to determine these parameters.
  # This, and future versions, will hand-pick.
  shape.invtau <- 1.81
  scale.invtau <- 0.00275
}

##==============================================================================
# begin exploring
MAP.deoptim <- deoptim.out$optim$bestmem
source('CARL_assimLikelihood.R') # load likelihood

my.logpost <- as.data.frame(matrix(data=NA,nrow=niter,ncol=length(parnames)))
my.changing.parameter <- my.logpost
names(my.logpost) <- parnames
names(my.changing.parameter) <- parnames

for (pp in 1:length(parnames)){
  
  my.parameters <- MAP.deoptim
  par.range <- seq(from=bound.lower[pp], to=bound.upper[pp], length = niter)
  
  if (pp==7)(
    par.range <- seq(from=0, to= 0.005, length=niter)
  )
  
  for (i in 1:length(par.range)){
    
    my.parameters[pp] <- par.range[i]
    my.changing.parameter[i,pp] <- my.parameters[pp]
    my.logpost[i,pp] <- log.post(parameters.in = my.parameters,
      
      # BRICK:
      parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
      slope.Ta2Tg.in=slope.Ta2Tg     , intercept.Ta2Tg.in=intercept.Ta2Tg                         ,
      ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
      oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
      obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
      bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
      luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       ,
      
      # DAIS:
      obs.in=obs.targets             , obs.err.in=obs.err.dais    , obs.step.in=obs.years         ,
      trends.ais.in=trends.ais       , trends.err.in=trends.err   , ind.trends.in=ind.trends      ,
      ind.norm.in=ind.relative       , alpha.var=alpha.var        , beta.var=beta.var             ,
      #slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , # already accounted for above 
      Tg.in=Tg.recon                 ,
      shape.lambda=shape.lambda      , rate.lambda=rate.lambda    , shape.Tcrit=shape.Tcrit       ,
      rate.Tcrit=rate.Tcrit          , 
      #l.aisfastdy=l.aisfastdy, # already accounted for above 
      #Ta.in=Ta  , Toc.in=Toc , # Tony left this commented
      SL.in=SL                       , dSL.in=dSL                 , 
      
      # Configure
      experts=experts                , alldata=alldata)
    
  }
  
}

## Begin plots
##############
# Make some trace plots
sequence1 = c(26,23,24,22,20,14,15,17,19,21,16,18,25,27,28)
sequence2 = c(5,6,7,8)
sequence3 = c(9,10,11,12,13)
sequence4 = c(1,2,3,4)
sequence5 = c(29,30,31,32,33)

cexlab = 2

## specify window for posterior marginal distributions
chain_complete <- amcmc.out$samples[5e6:20e6,] # 5 million to 20 million

## get posterior densities
chain_densities <- replicate(length(parnames), vector("list", 7), simplify = FALSE)
for (pp in 1:length(parnames)){
  chain_densities[[pp]] <- density(chain_complete[,pp])
}

# Trace plots for Antarctic parameters
if(TRUE){
  print("beginning antarctic plots")
  filename.trace1<- paste('../output_calibration/ExplorePlot',configure,'Set01','.png',sep='')
  png(filename=filename.trace1, width=1920, height=1080, units ="px")
  par(mfrow=c(3,5))
  for (pp in sequence1) {
    
    # prune posterior
    thisposterior = matrix(NA, nrow = length(chain_densities[[pp]]$x), ncol = 2)
    thisposterior[,1] = chain_densities[[pp]]$x
    thisposterior[,2] = chain_densities[[pp]]$y
    
    thisposterior_keep = as.data.frame(thisposterior)
    
    # prune posterior to 
    for (i in 1:length(chain_densities[[pp]]$x)){
      if (thisposterior[i,1]<min(my.changing.parameter[,pp]) | thisposterior[i,1]>max(my.changing.parameter[,pp])) {
        thisposterior_keep[i,2]=NA
      } else {
        thisposterior_keep[i,2]=thisposterior[i,2]
      }
    }
    
    thisposterior_keep <- thisposterior_keep[complete.cases(thisposterior_keep), ]
    
    #find matching length of logposts     
    this.logpost <- as.data.frame(matrix(NA, nrow = length(thisposterior_keep$V1), ncol = 1))
    this.changing.parameter <- this.logpost
    
    this.parameter <- MAP.deoptim
    par.range <- seq(from=bound.lower[pp], to=bound.upper[pp], length = length(thisposterior_keep$V1))
    
    if (pp==7)(
      par.range <- seq(from=0, to= 0.005, length=length(thisposterior_keep$V1))
    )
    
    for (i in 1:length(par.range)){
      this.parameter[pp] <- par.range[i]
      this.changing.parameter[i,1] <- this.parameter[pp]
      this.logpost[i,1] <- log.post(parameters.in = this.parameter,
                                    
                                    # BRICK:
                                    parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
                                    slope.Ta2Tg.in=slope.Ta2Tg     , intercept.Ta2Tg.in=intercept.Ta2Tg                         ,
                                    ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
                                    oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
                                    obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
                                    bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
                                    luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       ,
                                    
                                    # DAIS:
                                    obs.in=obs.targets             , obs.err.in=obs.err.dais    , obs.step.in=obs.years         ,
                                    trends.ais.in=trends.ais       , trends.err.in=trends.err   , ind.trends.in=ind.trends      ,
                                    ind.norm.in=ind.relative       , alpha.var=alpha.var        , beta.var=beta.var             ,
                                    #slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , # already accounted for above 
                                    Tg.in=Tg.recon                 ,
                                    shape.lambda=shape.lambda      , rate.lambda=rate.lambda    , shape.Tcrit=shape.Tcrit       ,
                                    rate.Tcrit=rate.Tcrit          , 
                                    #l.aisfastdy=l.aisfastdy, # already accounted for above 
                                    #Ta.in=Ta  , Toc.in=Toc , # Tony left this commented
                                    SL.in=SL                       , dSL.in=dSL                 , 
                                    
                                    # Configure
                                    experts=experts                , alldata=alldata)
    }
    
    # Scale posterior:
    log_max = max(this.logpost$V1)
    log_min = min(this.logpost$V1[which(this.logpost$V1 > 0)])
    span = log_max - log_min
    scale = span/(max(thisposterior_keep$V2)-min(thisposterior_keep$V2))
    
    thisposterior_scaled <- thisposterior_keep
    thisposterior_scaled[,2] <- scale*thisposterior_scaled[,2]+log_min
    names(thisposterior_scaled) <- c("x","y")
    
    # BEGIN PLOT!!!
    plot(x=thisposterior_scaled$x,y=thisposterior_scaled$y, type="l", lwd=3, 
         col="orange",xlab=parnames[pp],ylab="log post",
         cex.lab=cexlab)
    
    # MAP:
    abline(v=MAP.deoptim[pp], col="red")
    
    # marginal mode
    thisposterior_x <- thisposterior_scaled[which.max(thisposterior_scaled[,2]),1]
    abline(v=thisposterior_x,col="blue")
    
    # right axis
    axis(side = 4, 
         at=seq(from=min(thisposterior_scaled[,2]), to=max(thisposterior_scaled[,2]),length.out=5),
         col="orange",lwd=3, col.axis="orange")
    # mtext("log posterior", side = 4, line = 3, cex = 1.5)
    
    # priors:
    x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
    lines(x,scale*dunif(x, min = bound.lower[pp], max = bound.upper[pp])+log_min,
          col="black", lwd=1, lty=2)
    
    # log post
    lines(x=this.changing.parameter$V1,y=this.logpost$V1)
    
    # # legend
    # legend(x="topright",
    #        legend=c("Log posterior",
    #                 "Scaled marginal posterior",
    #                 "MAP estimate",
    #                 "Prior"),
    #        col=c("black","orange","red","black"), lwd=c(1,3,1,1),
    #        lty=c(1,1,1,2))
    
  }
  dev.off()
} # End of if true

# Trace Plots for the rest of the parameters.
if(TRUE){
  print("beginning trace plot 2")
  filename.trace2<- paste('../output_calibration/ExplorePlot',configure,'Set02','.png',sep='')
  png(filename=filename.trace2, width=1920, height=1080, units ="px")
  par(mfrow=c(3,6))
  for (pp in sequence2){
    
    # prune posterior
    thisposterior = matrix(NA, nrow = length(chain_densities[[pp]]$x), ncol = 2)
    thisposterior[,1] = chain_densities[[pp]]$x
    thisposterior[,2] = chain_densities[[pp]]$y
    
    thisposterior_keep = as.data.frame(thisposterior)
    
    # prune posterior to 
    for (i in 1:length(chain_densities[[pp]]$x)){
      if (thisposterior[i,1]<min(my.changing.parameter[,pp]) | thisposterior[i,1]>max(my.changing.parameter[,pp])) {
        thisposterior_keep[i,2]=NA
      } else {
        thisposterior_keep[i,2]=thisposterior[i,2]
      }
    }
    
    thisposterior_keep <- thisposterior_keep[complete.cases(thisposterior_keep), ]
    
    #find matching length of logposts     
    this.logpost <- as.data.frame(matrix(NA, nrow = length(thisposterior_keep$V1), ncol = 1))
    this.changing.parameter <- this.logpost
    
    this.parameter <- MAP.deoptim
    par.range <- seq(from=bound.lower[pp], to=bound.upper[pp], length = length(thisposterior_keep$V1))
    
    if (pp==7)(
      par.range <- seq(from=0, to= 0.005, length=length(thisposterior_keep$V1))
    )
    
    for (i in 1:length(par.range)){
      this.parameter[pp] <- par.range[i]
      this.changing.parameter[i,1] <- this.parameter[pp]
      this.logpost[i,1] <- log.post(parameters.in = this.parameter,
                                    
                                    # BRICK:
                                    parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
                                    slope.Ta2Tg.in=slope.Ta2Tg     , intercept.Ta2Tg.in=intercept.Ta2Tg                         ,
                                    ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
                                    oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
                                    obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
                                    bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
                                    luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       ,
                                    
                                    # DAIS:
                                    obs.in=obs.targets             , obs.err.in=obs.err.dais    , obs.step.in=obs.years         ,
                                    trends.ais.in=trends.ais       , trends.err.in=trends.err   , ind.trends.in=ind.trends      ,
                                    ind.norm.in=ind.relative       , alpha.var=alpha.var        , beta.var=beta.var             ,
                                    #slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , # already accounted for above 
                                    Tg.in=Tg.recon                 ,
                                    shape.lambda=shape.lambda      , rate.lambda=rate.lambda    , shape.Tcrit=shape.Tcrit       ,
                                    rate.Tcrit=rate.Tcrit          , 
                                    #l.aisfastdy=l.aisfastdy, # already accounted for above 
                                    #Ta.in=Ta  , Toc.in=Toc , # Tony left this commented
                                    SL.in=SL                       , dSL.in=dSL                 , 
                                    
                                    # Configure
                                    experts=experts                , alldata=alldata)
    }
    
    # Scale posterior:
    log_max = max(this.logpost$V1)
    log_min = min(this.logpost$V1[which(this.logpost$V1 > 0)])
    span = log_max - log_min
    scale = span/(max(thisposterior_keep$V2)-min(thisposterior_keep$V2))
    
    thisposterior_scaled <- thisposterior_keep
    thisposterior_scaled[,2] <- scale*thisposterior_scaled[,2]+log_min
    names(thisposterior_scaled) <- c("x","y")
    
    # BEGIN PLOT!!!
    plot(x=thisposterior_scaled$x,y=thisposterior_scaled$y, type="l", lwd=3, 
         col="orange",xlab=parnames[pp],ylab="log post",
         cex.lab=cexlab)
    
    # MAP:
    abline(v=MAP.deoptim[pp], col="red")
    
    # marginal mode
    thisposterior_x <- thisposterior_scaled[which.max(thisposterior_scaled[,2]),1]
    abline(v=thisposterior_x,col="blue")
    
    # right axis
    axis(side = 4, 
         at=seq(from=min(thisposterior_scaled[,2]), to=max(thisposterior_scaled[,2]),length.out=5),
         col="orange",lwd=3, col.axis="orange")
    # mtext("log posterior", side = 4, line = 3, cex = 1.5)
    
    # priors:
    x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
    lines(x,scale*dunif(x, min = bound.lower[pp], max = bound.upper[pp])+log_min,
          col="black", lwd=1, lty=2)
    
    # log post
    lines(x=this.changing.parameter$V1,y=this.logpost$V1)
    
    # # legend
    # legend(x="topright",
    #        legend=c("Log posterior",
    #                 "Scaled marginal posterior",
    #                 "MAP estimate",
    #                 "Prior"),
    #        col=c("black","orange","red","black"), lwd=c(1,3,1,1),
    #        lty=c(1,1,1,2))
    
  }
  
  for (pp in sequence3){
    
    # prune posterior
    thisposterior = matrix(NA, nrow = length(chain_densities[[pp]]$x), ncol = 2)
    thisposterior[,1] = chain_densities[[pp]]$x
    thisposterior[,2] = chain_densities[[pp]]$y
    
    thisposterior_keep = as.data.frame(thisposterior)
    
    # prune posterior to 
    for (i in 1:length(chain_densities[[pp]]$x)){
      if (thisposterior[i,1]<min(my.changing.parameter[,pp]) | thisposterior[i,1]>max(my.changing.parameter[,pp])) {
        thisposterior_keep[i,2]=NA
      } else {
        thisposterior_keep[i,2]=thisposterior[i,2]
      }
    }
    
    thisposterior_keep <- thisposterior_keep[complete.cases(thisposterior_keep), ]
    
    #find matching length of logposts     
    this.logpost <- as.data.frame(matrix(NA, nrow = length(thisposterior_keep$V1), ncol = 1))
    this.changing.parameter <- this.logpost
    
    this.parameter <- MAP.deoptim
    par.range <- seq(from=bound.lower[pp], to=bound.upper[pp], length = length(thisposterior_keep$V1))
    
    if (pp==7)(
      par.range <- seq(from=0, to= 0.005, length=length(thisposterior_keep$V1))
    )
    
    for (i in 1:length(par.range)){
      this.parameter[pp] <- par.range[i]
      this.changing.parameter[i,1] <- this.parameter[pp]
      this.logpost[i,1] <- log.post(parameters.in = this.parameter,
                                    
                                    # BRICK:
                                    parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
                                    slope.Ta2Tg.in=slope.Ta2Tg     , intercept.Ta2Tg.in=intercept.Ta2Tg                         ,
                                    ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
                                    oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
                                    obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
                                    bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
                                    luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       ,
                                    
                                    # DAIS:
                                    obs.in=obs.targets             , obs.err.in=obs.err.dais    , obs.step.in=obs.years         ,
                                    trends.ais.in=trends.ais       , trends.err.in=trends.err   , ind.trends.in=ind.trends      ,
                                    ind.norm.in=ind.relative       , alpha.var=alpha.var        , beta.var=beta.var             ,
                                    #slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , # already accounted for above 
                                    Tg.in=Tg.recon                 ,
                                    shape.lambda=shape.lambda      , rate.lambda=rate.lambda    , shape.Tcrit=shape.Tcrit       ,
                                    rate.Tcrit=rate.Tcrit          , 
                                    #l.aisfastdy=l.aisfastdy, # already accounted for above 
                                    #Ta.in=Ta  , Toc.in=Toc , # Tony left this commented
                                    SL.in=SL                       , dSL.in=dSL                 , 
                                    
                                    # Configure
                                    experts=experts                , alldata=alldata)
    }
    
    # Scale posterior:
    log_max = max(this.logpost$V1)
    log_min = min(this.logpost$V1[which(this.logpost$V1 > 0)])
    span = log_max - log_min
    scale = span/(max(thisposterior_keep$V2)-min(thisposterior_keep$V2))
    
    thisposterior_scaled <- thisposterior_keep
    thisposterior_scaled[,2] <- scale*thisposterior_scaled[,2]+log_min
    names(thisposterior_scaled) <- c("x","y")
    
    # BEGIN PLOT!!!
    plot(x=thisposterior_scaled$x,y=thisposterior_scaled$y, type="l", lwd=3, 
         col="orange",xlab=parnames[pp],ylab="log post",
         cex.lab=cexlab)
    
    # MAP:
    abline(v=MAP.deoptim[pp], col="red")
    
    # marginal mode
    thisposterior_x <- thisposterior_scaled[which.max(thisposterior_scaled[,2]),1]
    abline(v=thisposterior_x,col="blue")
    
    # right axis
    axis(side = 4, 
         at=seq(from=min(thisposterior_scaled[,2]), to=max(thisposterior_scaled[,2]),length.out=5),
         col="orange",lwd=3, col.axis="orange")
    # mtext("log posterior", side = 4, line = 3, cex = 1.5)
    
    # priors:
    x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
    lines(x,scale*dunif(x, min = bound.lower[pp], max = bound.upper[pp])+log_min,
          col="black", lwd=1, lty=2)
    
    # log post
    lines(x=this.changing.parameter$V1,y=this.logpost$V1)
    
    # # legend
    # legend(x="topright",
    #        legend=c("Log posterior",
    #                 "Scaled marginal posterior",
    #                 "MAP estimate",
    #                 "Prior"),
    #        col=c("black","orange","red","black"), lwd=c(1,3,1,1),
    #        lty=c(1,1,1,2))
    
  }
  
  for (pp in sequence4){
    

    
    # prune posterior
    thisposterior = matrix(NA, nrow = length(chain_densities[[pp]]$x), ncol = 2)
    thisposterior[,1] = chain_densities[[pp]]$x
    thisposterior[,2] = chain_densities[[pp]]$y
    
    thisposterior_keep = as.data.frame(thisposterior)
    
    # prune posterior to 
    for (i in 1:length(chain_densities[[pp]]$x)){
      if (thisposterior[i,1]<min(my.changing.parameter[,pp]) | thisposterior[i,1]>max(my.changing.parameter[,pp])) {
        thisposterior_keep[i,2]=NA
      } else {
        thisposterior_keep[i,2]=thisposterior[i,2]
      }
    }
    
    thisposterior_keep <- thisposterior_keep[complete.cases(thisposterior_keep), ]
    
    #find matching length of logposts     
    this.logpost <- as.data.frame(matrix(NA, nrow = length(thisposterior_keep$V1), ncol = 1))
    this.changing.parameter <- this.logpost
    
    this.parameter <- MAP.deoptim
    par.range <- seq(from=bound.lower[pp], to=bound.upper[pp], length = length(thisposterior_keep$V1))
    
    if (pp==7)(
      par.range <- seq(from=0, to= 0.005, length=length(thisposterior_keep$V1))
    )
    
    for (i in 1:length(par.range)){
      this.parameter[pp] <- par.range[i]
      this.changing.parameter[i,1] <- this.parameter[pp]
      this.logpost[i,1] <- log.post(parameters.in = this.parameter,
                                    
                                    # BRICK:
                                    parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
                                    slope.Ta2Tg.in=slope.Ta2Tg     , intercept.Ta2Tg.in=intercept.Ta2Tg                         ,
                                    ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
                                    oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
                                    obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
                                    bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
                                    luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       ,
                                    
                                    # DAIS:
                                    obs.in=obs.targets             , obs.err.in=obs.err.dais    , obs.step.in=obs.years         ,
                                    trends.ais.in=trends.ais       , trends.err.in=trends.err   , ind.trends.in=ind.trends      ,
                                    ind.norm.in=ind.relative       , alpha.var=alpha.var        , beta.var=beta.var             ,
                                    #slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , # already accounted for above 
                                    Tg.in=Tg.recon                 ,
                                    shape.lambda=shape.lambda      , rate.lambda=rate.lambda    , shape.Tcrit=shape.Tcrit       ,
                                    rate.Tcrit=rate.Tcrit          , 
                                    #l.aisfastdy=l.aisfastdy, # already accounted for above 
                                    #Ta.in=Ta  , Toc.in=Toc , # Tony left this commented
                                    SL.in=SL                       , dSL.in=dSL                 , 
                                    
                                    # Configure
                                    experts=experts                , alldata=alldata)
    }
    
    # Scale posterior:
    log_max = max(this.logpost$V1)
    log_min = min(this.logpost$V1[which(this.logpost$V1 > 0)])
    span = log_max - log_min
    scale = span/(max(thisposterior_keep$V2)-min(thisposterior_keep$V2))
    
    thisposterior_scaled <- thisposterior_keep
    thisposterior_scaled[,2] <- scale*thisposterior_scaled[,2]+log_min
    names(thisposterior_scaled) <- c("x","y")
    
    # BEGIN PLOT!!!
    plot(x=thisposterior_scaled$x,y=thisposterior_scaled$y, type="l", lwd=3, 
         col="orange",xlab=parnames[pp],ylab="log post",
         cex.lab=cexlab)
    
    # MAP:
    abline(v=MAP.deoptim[pp], col="red")
    
    # marginal mode
    thisposterior_x <- thisposterior_scaled[which.max(thisposterior_scaled[,2]),1]
    abline(v=thisposterior_x,col="blue")
    
    # right axis
    axis(side = 4, 
         at=seq(from=min(thisposterior_scaled[,2]), to=max(thisposterior_scaled[,2]),length.out=5),
         col="orange",lwd=3, col.axis="orange")
    # mtext("log posterior", side = 4, line = 3, cex = 1.5)
    
    # priors:
    x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
    lines(x,scale*dunif(x, min = bound.lower[pp], max = bound.upper[pp])+log_min,
          col="black", lwd=1, lty=2)
    
    # log post
    lines(x=this.changing.parameter$V1,y=this.logpost$V1)
    
    # # legend
    # legend(x="topright",
    #        legend=c("Log posterior",
    #                 "Scaled marginal posterior",
    #                 "MAP estimate",
    #                 "Prior"),
    #        col=c("black","orange","red","black"), lwd=c(1,3,1,1),
    #        lty=c(1,1,1,2))
  
  }
  
  for (pp in sequence5){
    
    # prune posterior
    thisposterior = matrix(NA, nrow = length(chain_densities[[pp]]$x), ncol = 2)
    thisposterior[,1] = chain_densities[[pp]]$x
    thisposterior[,2] = chain_densities[[pp]]$y
    
    thisposterior_keep = as.data.frame(thisposterior)
    
    # prune posterior to 
    for (i in 1:length(chain_densities[[pp]]$x)){
      if (thisposterior[i,1]<min(my.changing.parameter[,pp]) | thisposterior[i,1]>max(my.changing.parameter[,pp])) {
        thisposterior_keep[i,2]=NA
      } else {
        thisposterior_keep[i,2]=thisposterior[i,2]
      }
    }
    
    thisposterior_keep <- thisposterior_keep[complete.cases(thisposterior_keep), ]
    
    #find matching length of logposts     
    this.logpost <- as.data.frame(matrix(NA, nrow = length(thisposterior_keep$V1), ncol = 1))
    this.changing.parameter <- this.logpost
    
    this.parameter <- MAP.deoptim
    par.range <- seq(from=bound.lower[pp], to=bound.upper[pp], length = length(thisposterior_keep$V1))
    
    if (pp==7)(
      par.range <- seq(from=0, to= 0.005, length=length(thisposterior_keep$V1))
    )
    
    for (i in 1:length(par.range)){
      this.parameter[pp] <- par.range[i]
      this.changing.parameter[i,1] <- this.parameter[pp]
      this.logpost[i,1] <- log.post(parameters.in = this.parameter,
                                    
                                    # BRICK:
                                    parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
                                    slope.Ta2Tg.in=slope.Ta2Tg     , intercept.Ta2Tg.in=intercept.Ta2Tg                         ,
                                    ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
                                    oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
                                    obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
                                    bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
                                    luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       ,
                                    
                                    # DAIS:
                                    obs.in=obs.targets             , obs.err.in=obs.err.dais    , obs.step.in=obs.years         ,
                                    trends.ais.in=trends.ais       , trends.err.in=trends.err   , ind.trends.in=ind.trends      ,
                                    ind.norm.in=ind.relative       , alpha.var=alpha.var        , beta.var=beta.var             ,
                                    #slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , # already accounted for above 
                                    Tg.in=Tg.recon                 ,
                                    shape.lambda=shape.lambda      , rate.lambda=rate.lambda    , shape.Tcrit=shape.Tcrit       ,
                                    rate.Tcrit=rate.Tcrit          , 
                                    #l.aisfastdy=l.aisfastdy, # already accounted for above 
                                    #Ta.in=Ta  , Toc.in=Toc , # Tony left this commented
                                    SL.in=SL                       , dSL.in=dSL                 , 
                                    
                                    # Configure
                                    experts=experts                , alldata=alldata)
    }
    
    # Scale posterior:
    log_max = max(this.logpost$V1)
    log_min = min(this.logpost$V1[which(this.logpost$V1 > 0)])
    span = log_max - log_min
    scale = span/(max(thisposterior_keep$V2)-min(thisposterior_keep$V2))
    
    thisposterior_scaled <- thisposterior_keep
    thisposterior_scaled[,2] <- scale*thisposterior_scaled[,2]+log_min
    names(thisposterior_scaled) <- c("x","y")
    
    # BEGIN PLOT!!!
    plot(x=thisposterior_scaled$x,y=thisposterior_scaled$y, type="l", lwd=3, 
         col="orange",xlab=parnames[pp],ylab="log post",
         cex.lab=cexlab)
    
    # MAP:
    abline(v=MAP.deoptim[pp], col="red")
    
    # marginal mode
    thisposterior_x <- thisposterior_scaled[which.max(thisposterior_scaled[,2]),1]
    abline(v=thisposterior_x,col="blue")
    
    # right axis
    axis(side = 4, 
         at=seq(from=min(thisposterior_scaled[,2]), to=max(thisposterior_scaled[,2]),length.out=5),
         col="orange",lwd=3, col.axis="orange")
    # mtext("log posterior", side = 4, line = 3, cex = 1.5)
    
    # priors:
    x = seq(from=bound.lower[pp],to=bound.upper[pp],length.out=1000)
    lines(x,scale*dunif(x, min = bound.lower[pp], max = bound.upper[pp])+log_min,
          col="black", lwd=1, lty=2)
    
    # log post
    lines(x=this.changing.parameter$V1,y=this.logpost$V1)
    
    # # legend
    # legend(x="topright",
    #        legend=c("Log posterior",
    #                 "Scaled marginal posterior",
    #                 "MAP estimate",
    #                 "Prior"),
    #        col=c("black","orange","red","black"), lwd=c(1,3,1,1),
    #        lty=c(1,1,1,2))
    
  }
  
  
  dev.off()
} # End of if true

## Save workspace image
t.end <- proc.time()
time.elapsed <- t.end - t.beg
save.image(file=filename.saveprogress)
print("line 882")
print(time.elapsed)
print("TEDDY IS A CUTE PUPPER")
