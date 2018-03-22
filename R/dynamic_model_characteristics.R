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
control.start = as.Date("2014-03-08"); control.end = as.Date("2014-03-22")
control.shift = as.numeric((control.start-start.date)/365); control.range = as.numeric(control.end - control.start + 14)/365 # Parameters for control measure range

# - - - - - - - - 
# Set up timeseries

#timeser14 = timeser14_All
#timeser14[timeser14$date>swap.date,] = timeser14_Sus[timeser14_Sus$date>swap.date,]
timeser14_All = timeser14_All[timeser14_All$date < cutt.date,] 
timeser14_Sus = timeser14_Sus[timeser14_Sus$date < cutt.date,] 
timeser14_Sus$Central = timeser14_Sus$Central - timeser14_All$Central

# Set up weather fluctuation prior
weather.data.daily = read.csv(paste("data/suva_temp_data.csv",sep=""), stringsAsFactors = F); weather.data.daily$lsd = as.Date(weather.data.daily$lsd) # Load climate data

weather.data = read.csv(paste("data/data_Fiji_climate.csv",sep=""), stringsAsFactors = F) # Load climate data
temp.R0.data = read.csv(paste("data/outputTemp_R0.csv",sep=""), stringsAsFactors = F) %>% data.frame() # Load relative R0 from Mordecai
weather2014 = weather.data[weather.data$Date >= as.Date("2013-11-01") &weather.data$Date<= as.Date("2014-11-01"),] # select relevant dates
#weather2014.relative.R0 = c( max(temp.R0.data[min(weather2014$Av_temp)>=temp.R0.data$aegy.temps.DTR8,"R0.rel"]), temp.R0.data[max(weather2014$Av_temp)<=temp.R0.data$aegy.temps.DTR8,"R0.rel"][1] )
#rel.temp =  weather2014.relative.R0[1]/weather2014.relative.R0[2]

#weather2014.daily = weather.data.daily[weather.data.daily$lsd >= as.Date("2013-11-01") &weather.data.daily$lsd<= as.Date("2014-11-01"),] # select relevant dates

#plot(weather2014.daily$lsd,weather2014.daily$min_air_temp+ (weather2014.daily$max_air_temp-weather2014.daily$max_air_temp)/2,type="l")

# - - - - - - - - - - - - - - - 
# Fit sine wave to data
#tt_range = seq(start.date,start.date+200,1)
tt_range = as.Date(weather2014$Date) + 15
yy_year = sapply(tt_range, function(x){ y = max(which(x>=as.Date(weather2014$Date))); weather2014[y,"Av_temp"] })

#tt_range = as.Date(weather2014.daily$lsd) 
#yy_year = sapply(tt_range, function(x){ y = max(which(x>=as.Date(weather2014.daily$lsd))); weather2014.daily[y,"min_air_temp"] + (weather2014.daily[y,"max_air_temp"]-weather2014.daily[y,"min_air_temp"])/2 })

# Adjust for dt
tt_actual = as.numeric(tt_range - start.date) + dt

# - - - - - 
# Define temperature function
seasonaltemp <- function(x,theta){theta[["temp_base"]]*(1+ theta[["temp_amp"]]*sin((x - theta[["temp_shift"]])*2*pi/365) )  }
seasonaltemp_fit <- function(param,vals){ sum((seasonaltemp(tt_actual,param) - vals)^2)  }

# Fit sine function to temperature
param = c(temp_amp=0.1, temp_base=24, temp_shift=10)
optim_temp = optim(param, seasonaltemp_fit, method="L-BFGS-B",vals=yy_year, lower=c(rep(0,3)),upper=c(1,30,365), hessian=FALSE)$par

theta_fit = c(temp_amp=optim_temp[["temp_amp"]], temp_base=optim_temp[["temp_base"]], temp_shift=optim_temp[["temp_shift"]])

#theta_fit = c(temp_amp=0.2, temp_base=optim_temp[["temp_base"]], temp_shift=10)

# Plot temperature fit (to check:
tt_numeric = as.numeric(seq(start.date,start.date+300,1) - start.date)
plot(tt_actual,yy_year)
lines(tt_actual, seasonaltemp(tt_actual,theta_fit) )
lines(tt_numeric, seasonaltemp(tt_numeric,theta_fit) )

# - - - - - 
# Fit sine function to rainfall
yy_yearR = sapply(tt_range, function(x){ y = max(which(x>=as.Date(weather2014$Date))); weather2014[y,"Rain_av"] })

seasonalrain <- function(x,theta){theta[["rain_base"]]*(1+ theta[["rain_amp"]]*sin((x - theta[["rain_shift"]])*2*pi/365) )  }
seasonalrain_fit <- function(param,vals){ sum((seasonalrain(tt_actual,param) - vals)^2)  }

param = c(rain_amp=0.5, rain_base=210, rain_shift=10)

#tt_numeric = as.numeric(seq(start.date,start.date+400,1) - start.date)
#plot(tt_actual,yy_yearR)
#lines(tt_numeric, seasonalrain(tt_numeric,param) )

optim_rain = optim(param, seasonalrain_fit, method="L-BFGS-B",vals=yy_yearR, lower=c(rep(0,3)),upper=c(10,500,365), hessian=FALSE)$par
theta_fitRain = c(rain_amp=optim_rain[["rain_amp"]], rain_base=optim_rain[["rain_base"]], rain_shift=optim_rain[["rain_shift"]])


# Plot temperature fit (to check:
#tt_numeric = as.numeric(seq(start.date,start.date+400,1) - start.date)
# plot(tt_actual,yy_yearR)
# lines(tt_actual, seasonalrain(tt_actual,theta_fitRain) )
# lines(tt_numeric, seasonalrain(tt_numeric,theta_fitRain) )


# - - - - - - - - 
# Load Mordecai prior data

bite_temp <- function(temp){briere(temp, 0.0002016267, 40.04363, 13.76252)} # scale to 1 bite per day average at 26C baseline
EI_rate_temp <- function(temp){briere(temp, 6.111624e-05, 45.52661, 10.29733)/0.1103} # scale to 10 day value at 26C baseline
MD_rate_temp <- function(temp){briere(temp, 7.843086e-05, 39.10765, 11.56422)} 
eggs_per_female_temp<- function(temp){briere(temp, 0.008154744, 34.44309,  14.78163)} 
eggs_to_adult_temp <- function(temp){quad.2.trunc(temp, 13.58271, 38.288, -0.005985752)}
prob_to_h_temp <- function(temp){briere.trunc(temp, 0.0008328726, 35.77916, 17.22666)}
prob_to_v_temp <- function(temp){briere.trunc(temp, 0.0004880121, 37.38162, 12.67096)}
mortality_rate_temp <- function(temp){ 1/(quad.2(temp, 9.023289, 37.66479, -0.143324)/28.382 ) } # scale to value at 26C baseline

# Combined estimate of density change: eggs_per_female_temp * eggs_to_adult_temp * MD_rate_temp 
density_vary <- function(temp){eggs_per_female_temp(temp)*eggs_to_adult_temp(temp)*MD_rate_temp(temp)/0.67275} # scale to 1 at 26C baseline

# Plot values of different temperature variation

# xx_temp_plot = seq(22,26,0.1)
# plot(xx_temp_plot,bite_temp(xx_temp_plot),col="white",ylim=c(0,1.2))
# lines(xx_temp_plot,bite_temp(xx_temp_plot)*prob_to_h_temp(xx_temp_plot)*density_vary(xx_temp_plot)/(bite_temp(26)*prob_to_h_temp(26)*density_vary(26)),col="blue")
# lines(xx_temp_plot,bite_temp(xx_temp_plot)*prob_to_v_temp(xx_temp_plot)/(bite_temp(26)*prob_to_v_temp(26)),col="red")
# lines(xx_temp_plot,mortality_rate_temp(xx_temp_plot),col="green")
# lines(xx_temp_plot,EI_rate_temp(xx_temp_plot),col="grey")


# - - - - - - - - 
# Set up model characterists

totalSS = c(1,1,1) # number of sites in total
locationtab = c("Central","Central17")
locationtabF = c("Central","Western","Northern","Eastern") # labels


locnn <- 2 #length(locationtab) # Need this =2 for MCMC loop to work (old multi-site functionality)
iiH <- 1

# - - - - - - - - 
# ELISA and MIA results - AGE GROUPS: 0<= x <20 and x>=20

n_ELISA_C_D = c(21,42,86)   # DENV-3 immunity ELISA
n_ELISA_A_D = c(133,154,176) # DENV-3 immunity ELISA

n_Luminex_C_D3 = c(6,36,85)   # DENV-3 MIA immunity
n_Luminex_A_D3 = c(81,104,176) # DENV-3 MIA immunity

# Choose whether to use
if(use.ELISA.data == T){n_Luminex_C_D3 = n_ELISA_C_D; n_Luminex_A_D3 = n_ELISA_A_D}

# - - - - - - - - 
# Dengue priors
var_prior <- 0.1
var_priorBeta <- 0.1
# from Chan et al at 25C
prior_p_VEx <- c(10,var_prior); prior_p_Exp <- c(5.9,var_prior); prior_p_MuV <- c(14,var_prior); prior_p_Inf <- c(5,var_prior)

priorR0<-function(x){1+0*x} #dgamma(x,shape=prior_p_R0[1]/(prior_p_R0[2])+1, scale=prior_p_R0[2])}
priorExp<-function(x){dgamma(x,shape=prior_p_Exp[1]/(prior_p_Exp[2]), scale=prior_p_Exp[2])} #1+0*x} #
priorInf<-function(x){dgamma(x,shape=prior_p_Inf[1]/(prior_p_Inf[2]), scale=prior_p_Inf[2])} #1+0*x} #
priorVEx<-function(x){dgamma(x,shape=prior_p_VEx[1]/(prior_p_VEx[2]), scale=prior_p_VEx[2])} #1+0*x} #
priorMuV<-function(x){dgamma(x,shape=prior_p_MuV[1]/(prior_p_MuV[2]), scale=prior_p_MuV[2])} #1+0*x} #
priorAmplitude<-function(x){dgamma(x,shape=((1-rel.temp)/(1+rel.temp))/(0.1*var_prior), scale=(0.1*var_prior))} # Have strong prior on amplitude
priorDensity<-function(x){dgamma(x,shape=1/(var_prior), scale=(var_prior))} # Have weak prior on mosquito density
priorBetaH2M<-function(x){1 +0*x}#dgamma(x,shape=1/(var_priorBeta), scale=(var_priorBeta))} # Have strong prior on beta
priorBetaM2H<-function(x){1+0*x}#dgamma(x,shape=1/(var_priorBeta), scale=(var_priorBeta))} # Have strong prior on beta

itertab <- c(1); itertabM=c(1) # Iterate over locations in set up and MCMC
