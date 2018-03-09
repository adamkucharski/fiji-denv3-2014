# - - - - - - - - - - - - - - - - - - - - - - - 
# Transmission model characteristics
# Author: Adam Kucharski
# github.com/adamkucharski/fiji-denv3-2014
# - - - - - - - - - - - - - - - - - - - - - - - 

country.name="dengue_14"
epi.name="DENV"

# 2014 timeseries - split by week of chnge
timeser14_Sus <- read.csv(paste("data/Central_suspected.csv",sep=""), stringsAsFactors = F); timeser14_Sus$date = as.Date(timeser14_Sus$date)  # Load DLI dengue data
timeser14_All <- read.csv(paste("data/Central_lab_tested.csv",sep=""), stringsAsFactors = F); timeser14_All$date = as.Date(timeser14_All$date) # Load Lab tested cases
swap.date = as.Date("2014-03-08"); #Pick date of reporting change
cutt.date = as.Date("2014-08-01") # Pick max date to fit to - if smaller than swap.date, just fit case data  as.Date("2014-03-10")
start.date = as.Date("2013-11-04")
sample.collection.date = as.Date("2015-10-16")
dt =  7  #ode increment

# Set up parameters for control and seasonality
control.start = as.Date("2014-03-08"); control.end = as.Date("2014-03-22"); season.peak = as.Date("2014-02-15")
control.shift = as.numeric((control.start-start.date)/365); control.range = as.numeric(control.end - control.start + 14)/365 # Parameters for control measure range
season.shift = as.numeric((season.peak-start.date)/(2*pi*365))

# Set up timeseries

#timeser14 = timeser14_All
#timeser14[timeser14$date>swap.date,] = timeser14_Sus[timeser14_Sus$date>swap.date,]
timeser14_All = timeser14_All[timeser14_All$date < cutt.date,] 
timeser14_Sus = timeser14_Sus[timeser14_Sus$date < cutt.date,] 
timeser14_Sus$Central = timeser14_Sus$Central - timeser14_All$Central

# Set up weather fluctuation prior
weather.data <- read.csv(paste("data/data_Fiji_climate.csv",sep=""), stringsAsFactors = F) # Load climate data

temp.R0.data <- read.csv(paste("data/outputTemp_R0.csv",sep=""), stringsAsFactors = F) %>% data.frame() # Load relative R0 from Mordecai
weather2014 <- weather.data[weather.data$Date >= as.Date("2013-11-01") &weather.data$Date<= as.Date("2014-11-01"),] # select relevant dates
weather2014.relative.R0 <- c( max(temp.R0.data[min(weather2014$Av_temp)>=temp.R0.data$aegy.temps.DTR8,"R0.rel"]), temp.R0.data[max(weather2014$Av_temp)<=temp.R0.data$aegy.temps.DTR8,"R0.rel"][1] )
rel.temp <-  weather2014.relative.R0[1]/weather2014.relative.R0[2]

totalSS <- c(1,1,1) # number of sites in total
locationtab <- c("Central","Central17")
locationtabF <- c("Central","Western","Northern","Eastern") # labels


locnn <- 2 #length(locationtab) # Need this =2 for loop to work
iiH <- 1

# ELISA and MIA results - AGE GROUPS: 0<= x <20 and x>=20

n_ELISA_C_D = c(21,42,86)   # DENV-3 immunity ELISA
n_ELISA_A_D = c(133,154,176) # DENV-3 immunity ELISA

n_Luminex_C_D3 = c(6,36,85)   # DENV-3 MIA immunity
n_Luminex_A_D3 = c(81,104,176) # DENV-3 MIA immunity

# Choose whether to use
if(use.ELISA.data == T){n_Luminex_C_D3 = n_ELISA_C_D; n_Luminex_A_D3 = n_ELISA_A_D}


# Dengue priors
var_prior <- 0.1
# from Chan et al at 25C
prior_p_VEx <- c(15,var_prior); prior_p_Exp <- c(5.9,var_prior); prior_p_MuV <- c(8,var_prior); prior_p_Inf <- c(5,var_prior)
priorR0<-function(x){1+0*x} #dgamma(x,shape=prior_p_R0[1]/(prior_p_R0[2])+1, scale=prior_p_R0[2])}
priorExp<-function(x){dgamma(x,shape=prior_p_Exp[1]/(prior_p_Exp[2]), scale=prior_p_Exp[2])} #1+0*x} #
priorInf<-function(x){dgamma(x,shape=prior_p_Inf[1]/(prior_p_Inf[2]), scale=prior_p_Inf[2])} #1+0*x} #
priorVEx<-function(x){dgamma(x,shape=prior_p_VEx[1]/(prior_p_VEx[2]), scale=prior_p_VEx[2])} #1+0*x} #
priorMuV<-function(x){dgamma(x,shape=prior_p_MuV[1]/(prior_p_MuV[2]), scale=prior_p_MuV[2])} #1+0*x} #
priorAmplitude<-function(x){dgamma(x,shape=((1-rel.temp)/(1+rel.temp))/(0.1*var_prior), scale=(0.1*var_prior))} # Have strong prior on amplitude

itertab <- c(1); itertabM=c(1) # Iterate over locations in set up and MCMC
