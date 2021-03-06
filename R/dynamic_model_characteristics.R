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
cutt.date = as.Date("2014-06-01") # Pick max date to fit to - if smaller than swap.date, just fit case data  as.Date("2014-08-01")

start.date = as.Date("2013-11-04")
simulation.end.date = as.Date("2014-09-01") #  Note this is simulation range  # as.Date("2015-10-16")
dt =  7  #ode increment

# Set up parameters for control and seasonality
control.start = as.Date("2014-03-08"); control.end = as.Date("2014-03-22")
control.shift = as.numeric((control.start-start.date)/365); control.range = as.numeric(control.end - control.start + 14)/365 # Parameters for control measure range

# - - - - - - - - 
# Set up timeseries

timeser14_All = timeser14_All[timeser14_All$date < cutt.date,] 
timeser14_Sus = timeser14_Sus[timeser14_Sus$date < cutt.date,] 
timeser14_Sus$Central = timeser14_Sus$Central - timeser14_All$Central

# Omit initial wave of DLI cases, as before changeover:
timeser14_Sus[timeser14_Sus$date<"2014-02-01","Central"]=0

# Set up weather fluctuation prior
weather.data = read.csv(paste("data/data_Fiji_climate.csv",sep=""), stringsAsFactors = F) # Load climate data
weather2014 = weather.data[weather.data$Date >= as.Date("2013-10-01") &weather.data$Date<= as.Date("2014-10-01"),] # select relevant dates 

weather.data.daily = read.csv(paste("data/data_Fiji_climate_daily.csv",sep=""), stringsAsFactors = F) %>% data.frame() # Load daily temperature data
weather.data.daily$date = as.Date(weather.data.daily$date)

weather2014.daily = weather.data.daily[weather.data.daily$date >= as.Date("2013-10-01") &weather.data.daily$date<= as.Date("2016-11-01"),] # select relevant dates

tt_range = as.Date(weather2014.daily$date) # Date range
tt_actual = as.numeric(tt_range - start.date) + 7 # Convert to numeric and adjust for for dt in fitting

# Define temperature timeseries
yy_year1 = sapply(tt_range, function(x){ y = max(which(x>=as.Date(weather2014.daily$date))); weather2014.daily[y,"min_air_temp"] + (weather2014.daily[y,"max_air_temp"]-weather2014.daily[y,"min_air_temp"])/2 })

# Use actual data
tt_range2 = seq(start.date-7,start.date+800,1) # All dates
yy_yearT_mean = sapply(tt_range2, function(x){ y = yy_year1[which(x==(tt_range-7)):which(x==tt_range)]; mean(y) })
seasonaltemp <- function(x,theta){ yy_yearT_mean[round(x)+1]  } # pick out relative entry

# Temperature vals
tt_numeric = as.numeric(seq(start.date,start.date+300,1) - start.date)

# - - - - - 
# Fit moving average function to rainfall
tt_range = seq(start.date-30,start.date+820,1) # All dates
yy_yearR = sapply(tt_range, function(x){ y = max(which(x>=as.Date(weather.data$Date))); weather.data[y,"Rain_av"] })

tt_range2 = seq(start.date-7,start.date+800,1) # Create vector of times for moving average

# Interpolate between monthly averages
yy_yearR_mean = sapply(tt_range2, function(x){ 
                                  y = yy_yearR[which((x-15)==tt_range):which((x+15)==tt_range)] ; mean(y)
                                  })

tt_actual2 = as.numeric(tt_range2 - start.date) + 7

seasonalrain <- function(x,theta){ yy_yearR_mean[round(x)+1]  } # pick out relative entry

# - - - - - - - - 
# Load Mordecai prior data
EI_rate_temp <- function(temp){briere(temp, 6.111624e-05, 45.52661, 10.29733)/0.1017774} # scale to 1 at 25C baseline
mortality_rate_temp <- function(temp){ 1/(quad.2(temp, 9.023289, 37.66479, -0.143324)/29.00042 ) } # scale to 1 at 25C baseline

bite_temp <- function(temp){briere(temp, 0.0002016267, 40.04363, 13.76252)/0.2197016} # scale to 1 at 25C baseline
MD_rate_temp <- function(temp){briere(temp, 7.843086e-05, 39.10765, 11.56422)} 
eggs_per_female_temp<- function(temp){briere(temp, 0.008154744, 34.44309,  14.78163)} 
eggs_to_adult_temp <- function(temp){quad.2.trunc(temp, 13.58271, 38.288, -0.005985752)}
prob_to_h_temp <- function(temp){briere.trunc(temp, 0.0008328726, 35.77916, 17.22666)}
prob_to_v_temp <- function(temp){briere.trunc(temp, 0.0004880121, 37.38162, 12.67096)}

# Combined estimate of density change: eggs_per_female_temp * eggs_to_adult_temp * MD_rate_temp 
density_vary <- function(temp){eggs_per_female_temp(temp)*eggs_to_adult_temp(temp)*MD_rate_temp(temp)/0.5752381} # scale to 1 at 25C baseline

# - - - - - - - - 
# Set up model characterists

totalSS = c(1,1,1) # number of sites in total
locationtab = c("Central","Central17")
locationtabF = c("Central","Western","Northern","Eastern") # labels


locnn <- 2 # Need this =2 for MCMC loop to work (old multi-site functionality)
iiH <- 1

# - - - - - - - - 
# ELISA and MIA results - AGE GROUPS: 0<= x <20 and x>=20

n_ELISA_C_D = c(21,42,86)   # DENV-3 immunity ELISA
n_ELISA_A_D = c(133,154,177) # DENV-3 immunity ELISA

n_Luminex_C_D3 = c(6,36,86)   # DENV-3 MIA immunity
n_Luminex_A_D3 = c(81,104,177) # DENV-3 MIA immunity

if(use.ELISA.data == T){n_Luminex_C_D3 = n_ELISA_C_D; n_Luminex_A_D3 = n_ELISA_A_D} # Choose which to use


# - - - - - - - - 
# Dengue priors
var_prior <- 0.1

# from Chan et al at 25C
prior_p_Exp <- c(5.9,var_prior); 

# from Mordecai et al function above
prior_p_VEx <- c(10,var_prior); prior_p_MuV <- c(8,var_prior); prior_p_Inf <- c(5,var_prior)

prior_betaV <- c(0.25,var_prior)

priorR0<-function(x){1+0*x} # deprecated
priorExp<-function(x){dgamma(x,shape=prior_p_Exp[1]/(prior_p_Exp[2]), scale=prior_p_Exp[2])} #1+0*x} #
priorInf<-function(x){dgamma(x,shape=prior_p_Inf[1]/(prior_p_Inf[2]), scale=prior_p_Inf[2])} #1+0*x} #
priorVEx<-function(x){dgamma(x,shape=prior_p_VEx[1]/(prior_p_VEx[2]), scale=prior_p_VEx[2])} #1+0*x} #
priorMuV<-function(x){dgamma(x,shape=prior_p_MuV[1]/(prior_p_MuV[2]), scale=prior_p_MuV[2])} #1+0*x} #
priorDensity<-function(x){ifelse(x<50,1,0)} # Have flat prior on mosquito density

priorBeta_v<-function(x){dgamma(x,shape=prior_betaV[1]/prior_betaV[2], scale=prior_betaV[2])} # biting rate
priorBeta_h<-function(x){dgamma(x,shape=1/(var_prior), scale=(var_priorBeta))} # deprecated

priorAtRisk<-function(x){ifelse(x>0.8,1,0)} # Have strong prior on beta 1+ 0*x}#

itertab <- c(1); itertabM=c(1) # Iterate over locations in set up and MCMC


c(dbinom(30,size=100,prob=0.3,log=T), dbinom(70,size=100,prob=0.7,log=T)) %>% sum()

c(dbinom(10,size=100,prob=0.1,log=T), dbinom(50,size=100,prob=0.5,log=T)) %>% sum()




