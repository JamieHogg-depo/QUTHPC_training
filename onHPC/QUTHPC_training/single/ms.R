## MASTER ## ------------------------------------------------------------------

# Set library
.libPaths('r_lib')

# set global variables
base_folder='QUTHPC_training/single'
loop_output_file='QUTHPC_training/single/outputs/r/asra1_Asthma_Persons'
model_spec='asra1'
sex='Persons'
condition='Asthma'
niter=800
nburnin=400 
thin=1
nchains=4

# load packages
library(tidyverse)
library(readxl)
library(sf)
library(spdep)
library(nimble)
library(Matrix)
library(HDInterval)
library(bayesplot)

# Load functions and maps
source(paste0(base_folder, "/r_src/functions.R"))
load(paste0(base_folder, "/data/mappopDATA.Rdata"))

## Run file ## ----------------------------------------------------------------
source(paste0(base_folder, "/r_src/", model_spec, ".R"))




