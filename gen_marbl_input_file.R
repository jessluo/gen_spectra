#!/usr/bin/env Rscript
# --- load libraries ---
library('stringr')
library('dplyr')
library('reshape2')
`%ni%` = Negate(`%in%`) 

# --- read in data ---
ss <- read.csv("data/phytoplankton_input_data.csv", as.is=TRUE)
sz <- read.csv("data/zooplankton_input_data.csv", as.is=TRUE)
sg <- read.csv("data/grazing_input_data.csv", as.is=TRUE)

# --- namelist variables ---
namelist_vars_autotrophs <- c("sname", "lname", "Nfixer", "imp_calcifier", 
"exp_calcifier", "silicifier", "kFe", "kPO4", "kDOP", "kNO3", "kNH4", "kSiO3", 
"Qp_fixed", "gQfe_0", "gQfe_min", "alphaPI_per_day", "PCref_per_day", 
"thetaN_max", "loss_thres", "loss_thres2", "temp_thres", "mort_per_day", 
"mort2_per_day", "agg_rate_max", "agg_rate_min", "loss_poc", "Ea")

namelist_vars_zooplankton <- c("sname", "lname", "z_mort_0_per_day", 
"basal_metabolic_rate_per_day", "loss_thres", "z_mort2_0_per_day", "Ea")

namelist_vars_grazing <- c("sname", "lname", "auto_ind_cnt", "zoo_ind_cnt", 
"grazing_function", "z_umax_0_per_day", "z_grz", "graze_zoo", "graze_poc", 
"graze_doc", "f_zoo_detr", "auto_ind.1.", "zoo_ind.1.")

# --- convert to long form ---
# phyto
ss <- ss[,c("type", "index", namelist_vars_autotrophs)]
ss$sname <- str_c("\'", ss$sname, "\'")
ss$lname <- str_c("\'", ss$lname, "\'")

ssM <- melt(ss, id.vars=c("type", "index"), variable.name="varname")
ssM <- arrange(ssM, index, varname)

# zoo
sz <- sz[,c("type", "index", namelist_vars_zooplankton)]
sz$sname <- str_c("\'", sz$sname, "\'")
sz$lname <- str_c("\'", sz$lname, "\'")

szM <- melt(sz, id.vars=c("type", "index"), variable.name="varname")
szM <- arrange(szM, index, varname)

# grazing
sg <- sg[,c("type", "index1", "index2", namelist_vars_grazing)]
sg$sname <- str_c("\'", sg$sname, "\'")
sg$lname <- str_c("\'", sg$lname, "\'")
sgM <- melt(sg, id.vars=c("type", "index1", "index2"), variable.name="varname")
sgM <- sgM[complete.cases(sgM),]

sgM <- arrange(sgM, index2, index1, varname)
sgM$varname <- str_replace(sgM$varname, ".1.", "(1)")


# --- convert to namelist input form ---
dstr1 <- str_c(ssM$type, "_settings", "(", ssM$index, ")%", ssM$varname, " = ", ssM$value)
dstr2 <- str_c(szM$type, "_settings", "(", szM$index, ")%", szM$varname, " = ", szM$value)
dstr3 <- str_c(sgM$type, "_relationship_settings", "(", sgM$index1, ",", sgM$index2, ")%", sgM$varname, " = ", sgM$value)

d <- c(dstr1, "", "", dstr2, "", "", dstr3)

# just in case the TRUE/FALSE do not get interpreted in the same way
d <- str_replace(d, "TRUE", "T")
d <- str_replace(d, "FALSE", "F")

# --- add headers ---
nP <- nrow(ss)
nZ <- nrow(sz)
max_grazer_prey_cnt <- max(sg$index1)

dhead <- c("PFT_defaults = 'user-specified'", 
		str_c("autotroph_cnt = ", nP),
		str_c("zooplankton_cnt = ", nZ),
		str_c("max_grazer_prey_cnt = ", max_grazer_prey_cnt))

d <- c(dhead, "", "", d) 

write.table(d, "marbl_size_structured.input", col.names=FALSE, row.names=FALSE, quote=FALSE)

