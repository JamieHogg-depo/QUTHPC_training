## MASTER ## ------------------------------------------------------------------

# load packages
library(tidyverse)
library(readxl)
library(sf)
library(spdep)
library(nimble)
library(Matrix)
library(HDInterval)

# Load functions and maps
base_loc <- "WADOH/r_src/"
source(paste0(base_loc, "functions.R"))
load(paste0(base_loc, "mappopDATA.Rdata"))

## Run file ## ----------------------------------------------------------------
source(paste0(base_loc, Rfile, ".R"))




