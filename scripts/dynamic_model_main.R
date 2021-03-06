# - - - - - - - - - - - - - - - - - - - - - - - 
# Transmission modelling main code
# Author: Adam Kucharski
# github.com/adamkucharski/fiji-denv3-2014
# - - - - - - - - - - - - - - - - - - - - - - - 

library(foreach)
library(doMC)
library(deSolve)
library(mvtnorm)
library(MASS)
library(RColorBrewer)
library(magrittr)
library(lubridate)
library(colorspace)
library(coda)
library(IDPmisc)
library(VennDiagram)

rm(list=ls()) # Clear workspace

registerDoMC(4)  # Specify number of CPU cores

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Set up source functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

setwd("~/Documents/GitHub/fiji-denv3-2014/")
source("R/dynamic_model_analyse_posteriors.R")
source("R/dynamic_model_functions.R")
source("R/dynamic_model_mcmc.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

serological.data = c(F,T)


for(chainN in c(1)){ # Run MCMC across single or multiple chains

  for(s.type in serological.data){   # Fit models to two types of serological data (ELISA/MIA)
  
    use.ELISA.data = s.type # Fit to ELISA or MIA data
    
    # Fit using 4 model types (loop is inside function):
    # 1: SIR model, cases;  2: SIR model, serology and cases;  3: SIR + climate, serology and cases;  4: SIR + climate + control, serology and cases
    
    run_transmission_mcmc(MCMC.runs = 1e5) # set number of MCMC runs
    
  } # end serology loop
  
}# end chain loop


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot figures and outputs - Before plotting, need to define: use.ELISA.data = T or F
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Load relevant data
use.ELISA.data = F # Fit to ELISA or MIA data
chainN = 1 # pick MCMC chain
source("R/dynamic_model_characteristics.R",local=F)

# Compile the following:
# Figure 5 - this uses "Figure_5_FALSE3_4.pdf" output
# Table 5 - model parameters
# Supplementary model outputs: Figures S5-S10, Table S3-4 
for(p_pick in 1:4){
  plot_posteriors(p_pick)
  plot_figure_2014_dengue3(p_pick,long_time=F,DoubleFit=F)
}


# Table S6 - model comparison
model_comparison()

# Figure S4 - Illustration of parameters
plot_weather_and_control()


