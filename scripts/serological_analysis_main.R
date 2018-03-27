# - - - - - - - - - - - - - - - - - - - - - - - 
# Serological analysis main code
# Author: Adam Kucharski
# github.com/adamkucharski/fiji-denv3-2014
# - - - - - - - - - - - - - - - - - - - - - - - 

library(ggplot2)
library(gridExtra)
library(reshape2)
library(fitdistrplus)
library(mgcv)
library(lhs)
library(mixtools)
library(magrittr)
library(geosphere)
library(irr)
library(RColorBrewer)
library(colorspace)

rm(list=ls(all=TRUE)) # Clear workspace

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Set up source functions

setwd("~/Documents/GitHub/fiji-denv3-2014/")
source("R/serology_analyse_functions.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot figures and outputs

# Compile the following:
# Table 2
# Table 3
# Table S1
fit_ELISA_univariable()

# Compile the following:
# Table 4
# Figure 2
# Figure 3
remove_noise()


# Compile the following:
# Figure S1
plot_surveillance_data()


