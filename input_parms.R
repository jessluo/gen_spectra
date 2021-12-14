
####################################################
################### CONSTANTS ######################
####################################################
Tref <- 20
K <- 8.617330350E-5

####################################################
############### DEFAULT PARAMETERS #################
####################################################
# --- phytoplankton ---
median_diatom_vol_um3 <- 1.16E5
median_sp_mass_ugC <- 1.4E-4

# Kremer et al. 2017
Ea_pp <- 0.42 #0.42
Ea <- 0.32 # 0.317

PCref_beta <- -0.14
PCref_beta_pp_mp1 <- 0.05
PCref_beta_diat <- -0.14
PCref_const_pp <- 17
PCref_const_mp1 <- 13.3
PCref_const_diaz <- 11
PCref_const_mp <- 11.5
PCref_const_diat <- 12.4

# Edwards et al. 2015
# alphaPI
alphaPI_beta <- -0.15 # range (-0.17, -0.09)
alphaPI_const <- 10^(-0.25) # range (-1.5, -1.22)
alphaPI_diat_beta <- -0.12
alphaPI_diat_const <- 10^(-0.17)

#alphaPI_pp_const <- 10^(-0.6)

# thetaN_max
# derive a relationship based on alpha, Geider et al. 1997
# maximum is 8 mg Chla / mmol N (graph from Edwards et al. 2015)
# minimum is 0.8 mg Chla / mmol N
#thetaC_alpha_beta <- -0.4
#thetaC_alpha_const <- 0.47

#thetaC_alpha_beta_diat <- -0.5
#thetaC_alpha_const_diat <- 0.55

thetaC_beta_lg <- 0.026
thetaC_const_lg <- 0.285

thetaC_pp <- 0.12 # Pro ranges from 0.008-0.011 mg chl-a/mg C

# nutrients
# units of micromol N or P per liter, equivalent to mmol N or P per cubic meter
# Edwards et al. 2012
#VN_beta <- 0.82
#VN_const <- 10^(-8.1)
#VP_beta <- 1.0
#VP_const <- 10^(-9.1)
nut_beta_generic <- 0.3 # these were 0.3 before accoring to edwards et al. 2013
nut_beta_diatoms <- 0.3 # 0.3
diatom_nut_lim_factor = 1
diaz_Fe_lim_factor = 1.5 # N-fixation machinery?

kNO3_beta <- nut_beta_generic
kNO3_const <- 10^(-0.65)
kNO3_beta_diat <- nut_beta_diatoms
kNO3_const_diat <- kNO3_const * diatom_nut_lim_factor

kNO3_diaz = 2
kNH4_diaz = 0.2

kPO4_beta <- nut_beta_generic
kPO4_const <- 10^(-2.2)
kPO4_beta_diat <- nut_beta_diatoms
kPO4_const_diat <- kPO4_const * diatom_nut_lim_factor

kFe_const = 0.6e-5
kFe_beta = nut_beta_generic
kFe_const_diat = kFe_const * diatom_nut_lim_factor
kFe_beta_diat = nut_beta_diatoms

Qp_beta = 0.1
Qp_const = 0.022
Qp_fixed_diaz = 1/300
Qp_fixed_pp = 0.004651163
#Qp_fixed = 1/117

gQfe_min = 2.5E-6
gQfe_0 = 42E-6
gQfe_0_diaz = 60E-6

kSiO3_const <- 0.035

kNH4_const <- 0.02
kNH4_const_diat <- kNH4_const * diatom_nut_lim_factor

kDOP_marbl <- 0.5
kDOP_const <- 0.08
kDOP_const_diat <- kDOP_const * diatom_nut_lim_factor


# mortality & aggregation
#mort_beta <- -0.1
#mort_const <- 0.06
mort_PCref_factor <- 0.03
mort_PCref_diat_factor <-0.02
mort2_beta <- 0
mort2_const <- 0.035

## attmepting to do mortality for pp to be 0, and increase their mort per day
agg_min_beta <- 0
agg_min_const <- 0.01
agg_max_beta <- 0
agg_max_const <- 0.5

diat_agg_rate_min_multiplier = 1 

loss_poc_beta <- 0 #0.075
loss_poc_const <- 0 #0.02 / (median_sp_mass_ugC ^ loss_poc_beta)

phyto_PCref_modification = 2

# --- zooplankton ---
zoo_Ea = 0.55

# grazing parameters
# Hansen et al. 1997 converted to ESD_mm relationship
umax_beta = -0.25 # -0.21 +/- 0.02
umax_const = 4.5 # log a: -0.52+/- 0.29
km_beta <- 0#-0.03 # =-0.03 +/- 0.03
km_const <- 1.1 # log a: 0.35 +/- 0.58

opt_pred_prey_ratio_const = 12.5 # target 8:1 for ciliates, 18:1 for mesozooplankton
opt_pred_prey_ratio_beta = 0.15
#opt_pred_prey_sd_const = 9.5
#opt_pred_prey_sd_beta = 0.095

# assimilation and egestion parameters
AE = 0.65
graze_zoo_micro <- 0.3
graze_zoo_meso <- 0.3

graze_doc_frac <- 0.25

graze_poc_beta <- 0.3
graze_poc_const <- 0.2
graze_poc_diatom_scaling <- 1.2

# loss parameters
z_mort_beta <- -0.25 # McGurk et al. 1987 MEPS
z_mort_const <- 0.12 # it was 0.02 initially. 

#mesozoo_mort2_beta <- 0.2 # 
#mesozoo_mort2_const <- 0.005 # 
z_mort2_beta <- 0.2
z_mort2_const <- 0.2

#microzoo_mort2 = 0
lg_zoo_mort2_scaling = 4

zoo_loss_thres <- 0.075

# egestion 
zoo_detr_const = 0.1
zoo_detr_beta = 0.1
zoo_detr_diatom_scaling = 1

# grazing scalings
diat_grazing_scaling <- 0.65
med_diat_grazing_scaling <-1
lg_diat_grazing_scaling <-0.75
calcifiers_grazing_scaling <- 1.5 
zoo_grazing_scaling <- 0.45
microzoo_grazing_scaling = 0.85
microZ_diat_grazing = 0.45 # fraction relative to the feeding kernel-determined grazing rate
#mesoZ_diat_grazing = 2.3
microZ_pp_grazing_fixed = 0.35
#mesoZ_d2_grazing = 0.7

