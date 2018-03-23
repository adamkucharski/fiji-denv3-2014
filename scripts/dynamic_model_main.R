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

rm(list=ls()) # Clear workspace

clust1<-registerDoMC(4)  # Specify number of CPU cores

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

use.ELISA.data = F # Fit to ELISA or MIA data
exclude.p = 1 #  Subsequent weeks to skip after change from lab to DLI reporting (=2 implies skip 1 week) - DEPRECATED (BUT STILL USED)

# Load relevant data
source("R/dynamic_model_characteristics.R",local=F)


# Fit using 4 model types:
# 1: SIR model cases  2: SIR model serology and cases  3: SIR + climate  4: SIR + climate + control

run_transmission_mcmc(MCMC.runs = 1e2) # set number of MCMC runs   


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot figures and outputs - Before plotting, need to define: use.ELISA.data = T or F
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
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


