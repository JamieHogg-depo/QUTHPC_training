## MASTER ## ------------------------------------------------------------------

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




