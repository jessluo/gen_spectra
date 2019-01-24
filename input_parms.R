
####################################################
################### CONSTANTS ######################
####################################################
Tref <- 25
K <- 8.617330350E-5

####################################################
############### DEFAULT PARAMETERS #################
####################################################
# --- phytoplankton ---
median_diatom_vol_um3 <- 1.16E5
median_sp_mass_ugC <- 1.4E-4

# Kremer et al. 2017
Ea_pp <- 0.41722 #0.4
Ea <- 0.31664

PCref_beta <- -0.16
PCref_beta_pp_mp1 <- 0.08
PCref_const_pp <- 17.6
PCref_const_mp1 <- 13.85
PCref_const_diaz <- 10.85
PCref_const_mp <- 11.5
PCref_const_diat <- 11.65

# Edwards et al. 2015
# alphaPI
alphaPI_beta <- -0.13 # range (-0.17, -0.09)
alphaPI_const <- 10^(-0.25) # range (-1.5, -1.22)

# thetaN_max
# maximum is 8 mg Chla / mmol N (graph from Edwards et al. 2015)
# minimum is 0.8 mg Chla / mmol N
thetaC_beta_sm <- 0.0256
thetaC_coeff_sm <- 0.2659

thetaC_beta_lg <- -0.095
thetaC_coeff_lg <- 0.46

thetaC_beta_diat <- 0.132
thetaC_coeff_diat <- 0.055

# nutrients
# units of micromol N or P per liter, equivalent to mmol N or P per cubic meter
# Edwards et al. 2012
#VN_beta <- 0.82
#VN_const <- 10^(-8.1)
#VP_beta <- 1.0
#VP_const <- 10^(-9.1)
nut_beta_generic <- 0.3
nut_beta_diatoms <- 0.3
diatom_nut_lim_factor = 1/6

kNO3_beta <- nut_beta_generic
kNO3_const <- 10^(-0.84)
kNO3_beta_diat <- nut_beta_diatoms
kNO3_const_diat <- kNO3_const * diatom_nut_lim_factor

kNO3_diaz = 2
kNH4_diaz = 0.2

kPO4_beta <- nut_beta_generic
kPO4_const <- 10^(-2.5)
kPO4_beta_diat <- nut_beta_diatoms
kPO4_const_diat <- kPO4_const * diatom_nut_lim_factor

gQfe_const = exp(-11.2)
gQfe_beta = -0.15

gQfe_diat_scaling = 1.2

kFe_const = exp(-11)
kFe_beta = nut_beta_generic
kFe_const_diat = kFe_const / 2.1
kFe_beta_diat = 0.27

Qp_beta = 0.1
Qp_const = 0.022
Qp_fixed_diaz = 0.2735042735042735E-02
Qp_fixed_pp = 0.004651163
#Qp_fixed = 0.8547008547008548E-02

gQfe_min_diaz = 6E-6
gQfe_min = 3E-6

kSiO3_const <- 0.045

kNH4_const <- 0.006
kNH4_const_diat <- kNH4_const * diatom_nut_lim_factor

kDOP_marbl <- 0.5
kDOP_const <- 0.06
kDOP_const_diat <- kDOP_const * diatom_nut_lim_factor


# mortality & aggregation
mort_beta <- 0.1
mort_const <- exp(log(0.1)-mort_beta*(log(1E-2)))
mort2_beta <- 0.08
mort2_const <- exp(log(0.008)-mort2_beta*(log(1E-2)))

agg_min_beta <- -0.05
agg_min_const <- 0.02 / (median_sp_mass_ugC ^ agg_min_beta)
agg_max_beta <- -0.05
agg_max_const <- 0.58 / (median_sp_mass_ugC ^ agg_max_beta)

diat_agg_rate_min_multiplier = 1.25

loss_poc_beta <- 0.075
loss_poc_const <- 0.02 / (median_sp_mass_ugC ^ loss_poc_beta)

phyto_PCref_modification = 1.82

# --- zooplankton ---
zoo_Ea = 0.65

z_mort_beta <- -0.10
z_mort_const <- 0.10
z_mort2_beta <- -0.15
z_mort2_const <- 0.11

lg_zoo_mort2_scaling = 1

zoo_loss_thres <- 0.075

umax_beta = -0.45
umax_const = 6.2
km_beta <- -0.45
km_const <- 1.0

zoo_detr_const = 0.3
zoo_detr_beta = 0.2
zoo_detr_diatom_scaling = 1.5

graze_doc_beta <- 0.2
graze_doc_const <- 0.1

graze_poc_beta <- 0.2
graze_poc_const <- 0.29
diat_graze_poc_beta <- 0.03
diat_graze_poc_const <- 1.68

opt_pred_prey_ratio = 10
opt_pred_prey_sd_const = 9
opt_pred_prey_sd_beta = 0.1


diat_grazing_scaling <- 0.75
calcifiers_grazing_scaling <- 1 
zoo_grazing_scaling <- 1


