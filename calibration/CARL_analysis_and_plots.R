rm(list=ls())
graphics.off()

t.beg <- proc.time()

## Initial set-up
# configure which pipeline you want to create figures for:
# for hindcasts
ais.paleo.05 <- readRDS(file="complete_0115_ais_paleo_05.rds")
ais.paleo.50 <- readRDS(file="complete_0115_ais_paleo_50.rds")
ais.paleo.95 <- readRDS(file="complete_0115_ais_paleo_95.rds")
t.paleo <- readRDS(file="complete_0115_t_paleo.rds")
t.hind <- readRDS(file="complete_0115_t_hind.rds")

gsic.hind <- readRDS(file="complete_0115_gsic_hind.rds")
te.hind <- readRDS(file="complete_0115_te_hind.rds")
gis.hind <- readRDS(file="complete_0115_gis_hind.rds")
ais.hind <- readRDS(file="complete_0115_ais_hind.rds")
temp.hind <- readRDS(file="complete_0115_temp_hind.rds")
ocheat.hind <- readRDS(file="complete_0115_ocheat_hind.rds")
gsl.hind <- readRDS(file="complete_0115_gsl_hind.rds")

# for projections
slr.rcp26 <- readRDS(file="complete_0115_slr_rcp26.rds")
te.rcp26 <- readRDS(file="complete_0115_te_rcp26.rds") 
gis.rcp26 <- readRDS(file="complete_0115_gis_rcp26.rds")
gsic.rcp26 <- readRDS(file="complete_0115_gsic_rcp26.rds")
ais.rcp26 <- readRDS(file="complete_0115_ais_rcp26.rds")
temp.rcp26 <- readRDS(file="complete_0115_temp_rcp26.rds")
ocheat.rcp26 <- readRDS(file="complete_0115_ocheat_rcp26.rds")
t.proj <- readRDS(file="complete_0115_t_proj.rds")

# ## File name for the BRICK physical model output (netCDF4)
# filename.brick.magicc = '../output_model/BRICK-model_physical_control_02Apr2017.nc'
# filename.brick.simple = '../output_model/BRICK-model_physical_SIMPLE-GSIC_02Apr2017.nc'
# filename.brick.gmsl   = '../output_model/BRICK-model_physical_R07_03Apr2017.nc'
# 
# ## File name for the Van Dantzig model output (netCDF4)
# filename.vandantzig = '../output_model/VanDantzig_RCP85_control_02Apr2017.nc'
# 
# ## File name for the BRICK post-calibrated parameters (csv) (the BRICK output came from these guys)
# filename.parameters.magicc  = '../output_calibration/BRICK-model_postcalibratedParameters_control_02Apr2017.nc'
# filename.parameters.simple  = '../output_calibration/BRICK-model_postcalibratedParameters_SIMPLE-GSIC_02Apr2017.nc'
# filename.parameters.gmsl    = '../output_calibration/BRICK-model_drawcalibratedParameters_R07_03Apr2017.nc'
# 
# ## Other files
# filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_01Nov2016.csv"
# filename.fingerprints = "../fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"

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

# ncdata <- nc_open(filename.brick.magicc)
# t.hind      = ncvar_get(ncdata, 'time_hind')
# gsl.hind    = ncvar_get(ncdata, 'GlobalSeaLevel_hind')
# temp.hind   = ncvar_get(ncdata, 'temp_hind')
# ocheat.hind = ncvar_get(ncdata, 'ocheat_hind')
# gis.hind    = ncvar_get(ncdata, 'GIS_hind')
# te.hind     = ncvar_get(ncdata, 'TE_hind')
# gsic.hind   = ncvar_get(ncdata, 'GSIC_hind')
# t.paleo     = ncvar_get(ncdata, 'time_paleo_avg')
# ais.paleo.05= ncvar_get(ncdata, 'AIS_paleo_avg_q05')
# ais.paleo.50= ncvar_get(ncdata, 'AIS_paleo_avg_q50')
# ais.paleo.95= ncvar_get(ncdata, 'AIS_paleo_avg_q95')
# nc_close(ncdata)


# ais.paleo.05 <- dais.paleo.05
# ais.paleo.50 <- dais.paleo.50
# ais.paleo.95 <- dais.paleo.95

if(FALSE) # if loading a local environment
{

  load("test_process_environment.RData")  
    
  ais.paleo.05 <- dais.paleo.05.avg
  ais.paleo.50 <- dais.paleo.50.avg
  ais.paleo.95 <- dais.paleo.95.avg
  
  t.paleo <- date.avg
}

## Initialize arrays for the output
slr.05 = rep(NA,length(t.hind));		slr.50 = rep(NA,length(t.hind));		slr.95 = rep(NA,length(t.hind))
gsic.05 = rep(NA,length(t.hind));		gsic.50 = rep(NA,length(t.hind));		gsic.95 = rep(NA,length(t.hind))
gis.05 = rep(NA,length(t.hind));		gis.50 = rep(NA,length(t.hind));		gis.95 = rep(NA,length(t.hind))
te.05 = rep(NA,length(t.hind));			te.50 = rep(NA,length(t.hind));			te.95 = rep(NA,length(t.hind))
temp.05 = rep(NA,length(t.hind));		temp.50 = rep(NA,length(t.hind));		temp.95 = rep(NA,length(t.hind))
ocheat.05 = rep(NA,length(t.hind)); ocheat.50 = rep(NA,length(t.hind));	ocheat.95 = rep(NA,length(t.hind))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.hind)){
  c(slr.05[t], slr.50[t], slr.95[t])					:= quantile(gsl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gsic.05[t], gsic.50[t], gsic.95[t])				:= quantile(gsic.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gis.05[t], gis.50[t], gis.95[t])					:= quantile(gis.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(te.05[t], te.50[t], te.95[t])							:= quantile(te.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(temp.05[t], temp.50[t], temp.95[t])				:= quantile(temp.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(ocheat.05[t], ocheat.50[t], ocheat.95[t])	:= quantile(ocheat.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
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
png(paste(plotdir,'complete_hindcasts.png',sep=''), width=958, height=1080, units ="px")

if(TRUE){
  n.sig = 2         # how many sigma to plot around the obs?
  # layout(cbind(c(1,4,7),c(2,5,7),c(3,6,7)))
  layout(cbind(c(1,3,5,6),c(2,4,5,6)))
  # par(mai=c(.5,.7,.2,.08)) #c(bottom, left, top, right)
  par(mai=c(.5,.7,.3,.2)) #c(bottom, left, top, right) # omi = c(0,0,0,2) for right side space
  plotall <- FALSE
  
  
  if(plotall){
    # >>> SURFACE TEMPERATURE <<<
    plot(mod.time[midx.temp], temp.50[midx.temp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1850,2016), ylim=c(-.3,1.5), cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    mtext(side=2, text='Surface temperature\n[deg C]', line=2.3, cex=.9);
    # mtext(side=3, text=expression(bold(' a')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[midx.temp],rev(mod.time[midx.temp])), c(temp.95,rev(temp.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.temp.time[oidx.temp], obs.temp.norm[oidx.temp], type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.temp.time[oidx.temp],rev(obs.temp.time[oidx.temp])),
            c(obs.temp.norm[oidx.temp]+n.sig*obs.temp.err[oidx.temp],rev(obs.temp.norm[oidx.temp]-n.sig*obs.temp.err[oidx.temp])),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
    #legend(1839,1.65,c("5-95% range, modeled","2-sigma range, observations"),
    #       col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])), lwd=2, bty='n', cex=1.2)
    
    # >>> OCEAN HEAT <<<
    itmp=midx.ocheat[1]:nrow(ocheat.hind)
    plot(mod.time[itmp], ocheat.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
         ylab='', xlim=c(1950,2016), ylim=c(-50,50), cex.lab=1.2, cex.axis=1.2);
    mtext(side=1, text='Year', line=2.3, cex=.9);
    mtext(side=2, text='Ocean heat uptake\n[10^22 J]', line=2.3, cex=.9);
    # mtext(side=3, text=expression(bold(' b')), line=.25, cex=.9, adj=0);
    polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(ocheat.95[itmp],rev(ocheat.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
    lines(obs.ocheat.time, obs.ocheat.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
    lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
    polygon(c(obs.ocheat.time,rev(obs.ocheat.time)), c(obs.ocheat.norm+n.sig*obs.ocheat.err,rev(obs.ocheat.norm-n.sig*obs.ocheat.err)),
            col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  }
  # >>> GSIC <<<
  itmp=midx.gsic[1]:nrow(gsic.hind)
  plot(mod.time[itmp], gsic.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.01,.04), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Glaciers and\nsmall ice caps [m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  # mtext(side=3, text=expression(bold(' c')), line=.25, cex=.9, adj=0);
  mtext(side=3, text=expression(bold(' a) Glaciers and Ice Caps')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gsic.95[itmp],rev(gsic.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gsic.time, obs.gsic.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  
  # >>> GIS <<<
  plot(mod.time, gis.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b) Greenland Ice Sheet')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(gis.95,rev(gis.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.gis.time, obs.gis.norm, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  
  # >>> TE <<<
  x1971=seq(trends.te[1,4],trends.te[1,5])
  c1971=mean(x1971); yc1971=te.50[which(mod.time==mean(x1971))]
  lo1971 = yc1971+(trends.te[1,2]/1000)*(x1971-c1971)
  hi1971 = yc1971+(trends.te[1,3]/1000)*(x1971-c1971)
  y1971 = yc1971+(trends.te[1,1]/1000)*(x1971-c1971)

  x1993=seq(trends.te[2,4],trends.te[2,5])
  c1993=mean(x1993); yc1993=te.50[which(mod.time==mean(x1993))]
  lo1993 = yc1993+(trends.te[2,2]/1000)*(x1993-c1993)
  hi1993 = yc1993+(trends.te[2,3]/1000)*(x1993-c1993)
  y1993 = yc1993+(trends.te[2,1]/1000)*(x1993-c1993)
  # 
  # # Save TE observations
  # saveRDS(x1971,"TE_x1971.rds")
  # saveRDS(c1971,"TE_c1971.rds")
  # saveRDS(lo1971,"TE_lo1971.rds")
  # saveRDS(hi1971,"TE_hi1971.rds")
  # saveRDS(y1971,"TE_y1971.rds")
  # 
  # saveRDS(x1993,"TE_x1993.rds")
  # saveRDS(c1993,"TE_c1993.rds")
  # saveRDS(lo1993,"TE_lo1993.rds")
  # saveRDS(hi1993,"TE_hi1993.rds")
  # saveRDS(y1993,"TE_y1993.rds")
  
  # # load TE observations
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
  
  plot(mod.time, te.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1950,2016), ylim=c(-.04,.06), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Thermal expansion\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' c) Thermal Expansion')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(te.95,rev(te.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(x1971,y1971, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1971,y1971, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  polygon(c(x1971,rev(x1971)), c(lo1971,rev(hi1971)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  lines(x1993,y1993, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  # points(x1993,y1993, pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(x1993,rev(x1993)), c(lo1993,rev(hi1993)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  
  # >>> TOTAL SLR <<<
  plot(mod.time, slr.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  # mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=2, text='[m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d) Total Sea Level')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(slr.95,rev(slr.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  # lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  points(obs.sl.time, obs.sl.norm,pch=20, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]));
  lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  
  # >>> AIS PALEO, SMOOTHED <<<
  ipaleo=which(t.paleo==-149999):which(t.paleo==1)
  plot(t.paleo[ipaleo], ais.paleo.50[ipaleo], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
       ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]),
       # ylim=c(-20,10),
       ylim=c(-20,16),
       cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  # mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=2, text='[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e) Antarctic Ice Sheet')), line=.25, cex=.9, adj=0);
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(ais.paleo.95[ipaleo],rev(ais.paleo.05[ipaleo])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
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
                    "observations",
                    "2-sigma range, observations"),
         # fill=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])),
         lwd=c(2,8,NA,8), bty='n', cex=1.2,
         pch=c(NA,NA,20,NA),
         col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),
               rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5),
               rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]),
               rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7)
               ), 
               horiz = FALSE)
  
}
dev.off()
##==============================================================================
##==============================================================================


##==============================================================================
##==============================================================================
## Supplemental Figure
## -- Sea-level projections to 2100

# ncdata <- nc_open(filename.brick.magicc)
# t.proj = ncvar_get(ncdata, 'time_proj')
# slr.rcp26 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
# slr.rcp45 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
# slr.rcp85 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
# te.rcp26 = ncvar_get(ncdata, 'TE_RCP26')
# te.rcp45 = ncvar_get(ncdata, 'TE_RCP45')
# te.rcp85 = ncvar_get(ncdata, 'TE_RCP85')
# gis.rcp26 = ncvar_get(ncdata, 'GIS_RCP26')
# gis.rcp45 = ncvar_get(ncdata, 'GIS_RCP45')
# gis.rcp85 = ncvar_get(ncdata, 'GIS_RCP85')
# gsic.rcp26 = ncvar_get(ncdata, 'GSIC_RCP26')
# gsic.rcp45 = ncvar_get(ncdata, 'GSIC_RCP45')
# gsic.rcp85 = ncvar_get(ncdata, 'GSIC_RCP85')
# ais.rcp26 = ncvar_get(ncdata, 'AIS_RCP26')
# ais.rcp45 = ncvar_get(ncdata, 'AIS_RCP45')
# ais.rcp85 = ncvar_get(ncdata, 'AIS_RCP85')
# temp.rcp26 = ncvar_get(ncdata, 'temp_RCP26')
# temp.rcp45 = ncvar_get(ncdata, 'temp_RCP45')
# temp.rcp85 = ncvar_get(ncdata, 'temp_RCP85')
# ocheat.rcp26 = ncvar_get(ncdata, 'ocheat_RCP26')
# ocheat.rcp45 = ncvar_get(ncdata, 'ocheat_RCP45')
# ocheat.rcp85 = ncvar_get(ncdata, 'ocheat_RCP85')
# nc_close(ncdata)

#transpose:
slr.rcp26 <- t(slr.rcp26)
te.rcp26 <- t(te.rcp26)
gis.rcp26 <- t(gis.rcp26)
gsic.rcp26 <- t(gsic.rcp26)
ais.rcp26 <- t(ais.rcp26)
temp.rcp26 <- t(temp.rcp26)
ocheat.rcp26 <- t(ocheat.rcp26)


## Initialize arrays for the output
slr.rcp26.05 = rep(NA,length(t.proj)); slr.rcp26.50 = rep(NA,length(t.proj)); slr.rcp26.95 = rep(NA,length(t.proj))
# slr.rcp45.05 = rep(NA,length(t.proj)); slr.rcp45.50 = rep(NA,length(t.proj)); slr.rcp45.95 = rep(NA,length(t.proj))
# slr.rcp85.05 = rep(NA,length(t.proj)); slr.rcp85.50 = rep(NA,length(t.proj)); slr.rcp85.95 = rep(NA,length(t.proj))
ais.rcp26.05 = rep(NA,length(t.proj)); ais.rcp26.50 = rep(NA,length(t.proj)); ais.rcp26.95 = rep(NA,length(t.proj))
# ais.rcp45.05 = rep(NA,length(t.proj)); ais.rcp45.50 = rep(NA,length(t.proj)); ais.rcp45.95 = rep(NA,length(t.proj))
# ais.rcp85.05 = rep(NA,length(t.proj)); ais.rcp85.50 = rep(NA,length(t.proj)); ais.rcp85.95 = rep(NA,length(t.proj))
gis.rcp26.05 = rep(NA,length(t.proj)); gis.rcp26.50 = rep(NA,length(t.proj)); gis.rcp26.95 = rep(NA,length(t.proj))
# gis.rcp45.05 = rep(NA,length(t.proj)); gis.rcp45.50 = rep(NA,length(t.proj)); gis.rcp45.95 = rep(NA,length(t.proj))
# gis.rcp85.05 = rep(NA,length(t.proj)); gis.rcp85.50 = rep(NA,length(t.proj)); gis.rcp85.95 = rep(NA,length(t.proj))
gsic.rcp26.05 = rep(NA,length(t.proj)); gsic.rcp26.50 = rep(NA,length(t.proj)); gsic.rcp26.95 = rep(NA,length(t.proj))
# gsic.rcp45.05 = rep(NA,length(t.proj)); gsic.rcp45.50 = rep(NA,length(t.proj)); gsic.rcp45.95 = rep(NA,length(t.proj))
# gsic.rcp85.05 = rep(NA,length(t.proj)); gsic.rcp85.50 = rep(NA,length(t.proj)); gsic.rcp85.95 = rep(NA,length(t.proj))
te.rcp26.05 = rep(NA,length(t.proj)); te.rcp26.50 = rep(NA,length(t.proj)); te.rcp26.95 = rep(NA,length(t.proj))
# te.rcp45.05 = rep(NA,length(t.proj)); te.rcp45.50 = rep(NA,length(t.proj)); te.rcp45.95 = rep(NA,length(t.proj))
# te.rcp85.05 = rep(NA,length(t.proj)); te.rcp85.50 = rep(NA,length(t.proj)); te.rcp85.95 = rep(NA,length(t.proj))
temp.rcp26.05 = rep(NA,length(t.proj)); temp.rcp26.50 = rep(NA,length(t.proj)); temp.rcp26.95 = rep(NA,length(t.proj))
# temp.rcp45.05 = rep(NA,length(t.proj)); temp.rcp45.50 = rep(NA,length(t.proj)); temp.rcp45.95 = rep(NA,length(t.proj))
# temp.rcp85.05 = rep(NA,length(t.proj)); temp.rcp85.50 = rep(NA,length(t.proj)); temp.rcp85.95 = rep(NA,length(t.proj))
ocheat.rcp26.05 = rep(NA,length(t.proj)); ocheat.rcp26.50 = rep(NA,length(t.proj)); ocheat.rcp26.95 = rep(NA,length(t.proj))
# ocheat.rcp45.05 = rep(NA,length(t.proj)); ocheat.rcp45.50 = rep(NA,length(t.proj)); ocheat.rcp45.95 = rep(NA,length(t.proj))
# ocheat.rcp85.05 = rep(NA,length(t.proj)); ocheat.rcp85.50 = rep(NA,length(t.proj)); ocheat.rcp85.95 = rep(NA,length(t.proj))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.proj)){
  c(slr.rcp26.05[t], slr.rcp26.50[t], slr.rcp26.95[t]) := quantile(slr.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(slr.rcp45.05[t], slr.rcp45.50[t], slr.rcp45.95[t]) := quantile(slr.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(slr.rcp85.05[t], slr.rcp85.50[t], slr.rcp85.95[t]) := quantile(slr.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(ais.rcp26.05[t], ais.rcp26.50[t], ais.rcp26.95[t]) := quantile(ais.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(ais.rcp45.05[t], ais.rcp45.50[t], ais.rcp45.95[t]) := quantile(ais.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(ais.rcp85.05[t], ais.rcp85.50[t], ais.rcp85.95[t]) := quantile(ais.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gis.rcp26.05[t], gis.rcp26.50[t], gis.rcp26.95[t]) := quantile(gis.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(gis.rcp45.05[t], gis.rcp45.50[t], gis.rcp45.95[t]) := quantile(gis.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(gis.rcp85.05[t], gis.rcp85.50[t], gis.rcp85.95[t]) := quantile(gis.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gsic.rcp26.05[t], gsic.rcp26.50[t], gsic.rcp26.95[t]) := quantile(gsic.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(gsic.rcp45.05[t], gsic.rcp45.50[t], gsic.rcp45.95[t]) := quantile(gsic.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(gsic.rcp85.05[t], gsic.rcp85.50[t], gsic.rcp85.95[t]) := quantile(gsic.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(te.rcp26.05[t], te.rcp26.50[t], te.rcp26.95[t]) := quantile(te.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(te.rcp45.05[t], te.rcp45.50[t], te.rcp45.95[t]) := quantile(te.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(te.rcp85.05[t], te.rcp85.50[t], te.rcp85.95[t]) := quantile(te.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(temp.rcp26.05[t], temp.rcp26.50[t], temp.rcp26.95[t]) := quantile(temp.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(temp.rcp45.05[t], temp.rcp45.50[t], temp.rcp45.95[t]) := quantile(temp.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(temp.rcp85.05[t], temp.rcp85.50[t], temp.rcp85.95[t]) := quantile(temp.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(ocheat.rcp26.05[t], ocheat.rcp26.50[t], ocheat.rcp26.95[t]) := quantile(ocheat.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(ocheat.rcp45.05[t], ocheat.rcp45.50[t], ocheat.rcp45.95[t]) := quantile(ocheat.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  # c(ocheat.rcp85.05[t], ocheat.rcp85.50[t], ocheat.rcp85.95[t]) := quantile(ocheat.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
}

iproj = which(t.proj==2000):which(t.proj==2100)
i2100 = which(t.proj==2100)

print('==============================================================')
print('min/5%/50%/95%/max of 2100 sea level relative to 1986-2005:')
print(paste('RCP2.6: ',quantile(slr.rcp26[251,],c(0,.05,.50,.95,1))))
# print(paste('RCP4.5: ',quantile(slr.rcp45[251,],c(0,.05,.50,.95,1))))
# print(paste('RCP8.5: ',quantile(slr.rcp85[251,],c(0,.05,.50,.95,1))))
print('==============================================================')

i2050 <- which(t.proj==2050)
print('==============================================================')
print('min/5%/50%/95%/max of 2050 sea level relative to 1986-2005:')
print(paste('RCP2.6: ',quantile(slr.rcp26[i2050,],c(0,.05,.50,.95,1))))
# print(paste('RCP4.5: ',quantile(slr.rcp45[i2050,],c(0,.05,.50,.95,1))))
# print(paste('RCP8.5: ',quantile(slr.rcp85[i2050,],c(0,.05,.50,.95,1))))
print('==============================================================')

##==============================================================================

##
## 5-95% CI and median of 2100 SLR and all contributions under the RCP scenarios
##

# slr.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
#                               c(slr.rcp26.50[i2100],slr.rcp45.50[i2100],slr.rcp85.50[i2100]),
#                               c(slr.rcp26.05[i2100],slr.rcp45.05[i2100],slr.rcp85.05[i2100]),
#                               c(slr.rcp26.95[i2100],slr.rcp45.95[i2100],slr.rcp85.95[i2100])
# ))

# gsic.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
#                                c(gsic.rcp26.50[i2100],gsic.rcp45.50[i2100],gsic.rcp85.50[i2100]),
#                                c(gsic.rcp26.05[i2100],gsic.rcp45.05[i2100],gsic.rcp85.05[i2100]),
#                                c(gsic.rcp26.95[i2100],gsic.rcp45.95[i2100],gsic.rcp85.95[i2100])
# ))
# gis.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
#                               c(gis.rcp26.50[i2100],gis.rcp45.50[i2100],gis.rcp85.50[i2100]),
#                               c(gis.rcp26.05[i2100],gis.rcp45.05[i2100],gis.rcp85.05[i2100]),
#                               c(gis.rcp26.95[i2100],gis.rcp45.95[i2100],gis.rcp85.95[i2100])
# ))
# te.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
#                              c(te.rcp26.50[i2100],te.rcp45.50[i2100],te.rcp85.50[i2100]),
#                              c(te.rcp26.05[i2100],te.rcp45.05[i2100],te.rcp85.05[i2100]),
#                              c(te.rcp26.95[i2100],te.rcp45.95[i2100],te.rcp85.95[i2100])
# ))
# ais.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
#                               c(ais.rcp26.50[i2100],ais.rcp45.50[i2100],ais.rcp85.50[i2100]),
#                               c(ais.rcp26.05[i2100],ais.rcp45.05[i2100],ais.rcp85.05[i2100]),
#                               c(ais.rcp26.95[i2100],ais.rcp45.95[i2100],ais.rcp85.95[i2100])
# ))

slr.ci90 = data.frame( cbind( c('RCP2.6'),
                              c(slr.rcp26.50[i2100]),
                              c(slr.rcp26.05[i2100]),
                              c(slr.rcp26.95[i2100])
))

gsic.ci90 = data.frame( cbind( c('RCP2.6'),
                               c(gsic.rcp26.50[i2100]),
                               c(gsic.rcp26.05[i2100]),
                               c(gsic.rcp26.95[i2100])
))
gis.ci90 = data.frame( cbind( c('RCP2.6'),
                              c(gis.rcp26.50[i2100]),
                              c(gis.rcp26.05[i2100]),
                              c(gis.rcp26.95[i2100])
))
te.ci90 = data.frame( cbind( c('RCP2.6'),
                             c(te.rcp26.50[i2100]),
                             c(te.rcp26.05[i2100]),
                             c(te.rcp26.95[i2100])
))
ais.ci90 = data.frame( cbind( c('RCP2.6'),
                              c(ais.rcp26.50[i2100]),
                              c(ais.rcp26.05[i2100]),
                              c(ais.rcp26.95[i2100])
))
## get into format for Latex table
row.gsic = paste('Glaciers and small ice caps &',
                 1000*signif(gsic.rcp26.50[i2100],4),' (',1000*signif(gsic.rcp26.05[i2100],4),'-',1000*signif(gsic.rcp26.95[i2100],4),') &',
                 # 1000*signif(gsic.rcp45.50[i2100],4),' (',1000*signif(gsic.rcp45.05[i2100],4),'-',1000*signif(gsic.rcp45.95[i2100],4),') &',
                 # 1000*signif(gsic.rcp85.50[i2100],4),' (',1000*signif(gsic.rcp85.05[i2100],4),'-',1000*signif(gsic.rcp85.95[i2100],4),') \\',
                 sep='')
row.te = paste('Thermal expansion &',
               1000*signif(te.rcp26.50[i2100],4),' (',1000*signif(te.rcp26.05[i2100],4),'-',1000*signif(te.rcp26.95[i2100],4),') &',
               # 1000*signif(te.rcp45.50[i2100],4),' (',1000*signif(te.rcp45.05[i2100],4),'-',1000*signif(te.rcp45.95[i2100],4),') &',
               # 1000*signif(te.rcp85.50[i2100],4),' (',1000*signif(te.rcp85.05[i2100],4),'-',1000*signif(te.rcp85.95[i2100],4),') \\',
               sep='')
row.gis = paste('Greenland Ice Sheet &',
                1000*signif(gis.rcp26.50[i2100],4),' (',1000*signif(gis.rcp26.05[i2100],4),'-',1000*signif(gis.rcp26.95[i2100],4),') &',
                # 1000*signif(gis.rcp45.50[i2100],4),' (',1000*signif(gis.rcp45.05[i2100],4),'-',1000*signif(gis.rcp45.95[i2100],4),') &',
                # 1000*signif(gis.rcp85.50[i2100],4),' (',1000*signif(gis.rcp85.05[i2100],4),'-',1000*signif(gis.rcp85.95[i2100],4),') \\',
                sep='')
row.ais = paste('Antarctic Ice Sheet &',
                1000*signif(ais.rcp26.50[i2100],4),' (',1000*signif(ais.rcp26.05[i2100],4),'-',1000*signif(ais.rcp26.95[i2100],4),') &',
                # 1000*signif(ais.rcp45.50[i2100],4),' (',1000*signif(ais.rcp45.05[i2100],4),'-',1000*signif(ais.rcp45.95[i2100],4),') &',
                # 1000*signif(ais.rcp85.50[i2100],4),' (',1000*signif(ais.rcp85.05[i2100],4),'-',1000*signif(ais.rcp85.95[i2100],4),') \\',
                sep='')
row.slr = paste('Total sea level &',
                1000*signif(slr.rcp26.50[i2100],4),' (',1000*signif(slr.rcp26.05[i2100],4),'-',1000*signif(slr.rcp26.95[i2100],4),') &',
                # 1000*signif(slr.rcp45.50[i2100],4),' (',1000*signif(slr.rcp45.05[i2100],4),'-',1000*signif(slr.rcp45.95[i2100],4),') &',
                # 1000*signif(slr.rcp85.50[i2100],4),' (',1000*signif(slr.rcp85.05[i2100],4),'-',1000*signif(slr.rcp85.95[i2100],4),') \\',
                sep='')

# pdf(paste(plotdir,'projections_SLR_total.pdf',sep=''),width=3.5,height=2.45,colormodel='cmyk')
# par(mfrow=c(1,1))
# # RCP85
# par(mai=c(.65,.65,.20,.2))
# plot(t.proj[iproj],slr.rcp85.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
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
# lines(t.proj[iproj],slr.rcp26.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
# polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp26.95[iproj],rev(slr.rcp26.05[iproj])),
#         col=rgb(col26[1],col26[2],col26[3],.5), border=NA);
# # + legend
# legend(t.proj[iproj[1]],2,c("5-95% range,","RCP2.6","RCP4.5","RCP8.5"),
#        lty=c(NA,1,1,1), lwd=3, bty='n', cex=1,
#        col=c(NA , rgb(col26[1],col26[2],col26[3]) ,
#              rgb(col45[1],col45[2],col45[3]) , rgb(col85[1],col85[2],col85[3])))
# 
# dev.off()

# pdf(paste(plotdir,'projections_SLR_total_with_noise_and_normalization.pdf',sep=''),width=3.5,height=2.45,colormodel='cmyk')
png(paste(plotdir,'complete_projections.png',sep=''), width=866, height=516, units ="px")

par(mfrow=c(1,1))
# RCP85
par(mai=c(.65,.65,.20,.2))
plot(t.proj[iproj],slr.rcp26.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2, #ann='',
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
# lines(t.proj[iproj],slr.rcp26.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp26.95[iproj],rev(slr.rcp26.05[iproj])),
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

##==============================================================================
##==============================================================================

##==============================================================================
##==============================================================================
# PDF of SLR at 2100
# sl.2100 <- proj.rcp26$slr[,251]
# sl.2100.density <- density(sl.2100)
# 
# plot(sl.2100.density$x,sl.2100.density$y,
#      type="l",xlab="Global mean sea level at 2100 [m]", ylab='PDF', lwd=2, col="black"
#      , yaxt='n'
#      )
# 
# legend(x="topright",legend="Standard priors combined with expert assessment",
#        lty=1, lwd=2, col="black"
#        )


## Save workspace image
t.end <- proc.time()
time.elapsed <- t.end - t.beg
# save.image(file=filename.saveprogress)
print("line 718")
print(time.elapsed)
print("TEDDY IS A CUTE PUPPER")
