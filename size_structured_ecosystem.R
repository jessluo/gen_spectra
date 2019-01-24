#!/usr/bin/env Rscript

# -------------------------------------------------
# Wrapper script for generating namelist input files
# Use for MARBL size-structured ecosystem
# 
# Receives 1 input: Casename
#
# -------------------------------------------------

library(stringr)

args <- commandArgs(TRUE)

if(length(args)==0){
     stop ("Need casename as command line input")
}

CASENAME <- args[1]

dir = str_c("../cases/", CASENAME)
dir.create(dir)

file.copy(from="input_parms.R", to=str_c(dir, "/input_parms.R"))
file.copy(from="gen_input_data.R", to=str_c(dir, "/gen_input_data.R"))
file.copy(from="gen_marbl_input_file.R", to=str_c(dir, "/gen_marbl_input_file.R"))

setwd(dir)
source("gen_input_data.R")
source("gen_marbl_input_file.R")


