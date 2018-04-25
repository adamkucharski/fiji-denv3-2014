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

aa = Sys.time()

serological.data = c(F,T)

# Run across multiple chains
for(chainN in c(1)){

  # Fit models to two types of serological data (ELISA/MIA)
  
  for(s.type in serological.data){
  #foreach(kk = c(1:2)) %dopar% {
  
    use.ELISA.data = s.type #serological.data[kk] # Fit to ELISA or MIA data  s.type #
    
    # Fit using 4 model types (loop inside function):
    # 1: SIR model, cases;  2: SIR model, serology and cases;  3: SIR + climate, serology and cases;  4: SIR + climate + control, serology and cases
    
    run_transmission_mcmc(MCMC.runs = 2e5) # set number of MCMC runs   
    
  } # end serology loop
  
}# end chain loop

Sys.time() - aa

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot figures and outputs - Before plotting, need to define: use.ELISA.data = T or F
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# Load relevant data
use.ELISA.data = F # Fit to ELISA or MIA data
chainN = 1 # pick MCMC chain
source("R/dynamic_model_characteristics.R",local=F)

# Compile the following:
# Figure 5 - this uses "Figure_5_FALSE3_4.pdf" output
# Supplementary model figures S7-S9
# Table 1 (in two parts). Main table uses "Table_5_params_part1_FALSE3_4" output
for(p_pick in 1:4){
  plot_posteriors(p_pick)
  plot_figure_2014_dengue3(p_pick,long_time=F,DoubleFit=F)
}


# Table S3 - model comparison
model_comparison()

# Figure S3 - illustration
plot_weather_and_control()


