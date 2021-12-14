#!/usr/bin/env Rscript
# --- load libraries ---
library('stringr')
library('dplyr')
library('reshape2')
`%ni%` = Negate(`%in%`) 

source("input_parms.R")

####################################################
################## Autotrophs ######################
####################################################

# --- set up phytoplankton ---
nP <- 9

nMP <- ceiling((nP-2)/2)
nDIAT <- floor((nP-2)/2)
nMP_PP <- nMP + 1

mp_sizes <- 10^(seq(from=-0.05, to=2.2, len=nMP_PP))
diat_sizes <- 10^(seq(from=1.468, to=2.1345, len=nDIAT))

mp_sizes_log=seq(from=-0.05, to=2.2, len=nMP_PP)
dd=diff(mp_sizes_log)[1]
mp_sizes_min=10^(mp_sizes_log-dd/2)
mp_sizes_max=10^(mp_sizes_log+dd/2)
mp_size_range=paste0('[',round(mp_sizes_min,3),',',round(mp_sizes_max,3),')')

diaz_size=6.16667
diaz_size_min=10^(log10(diaz_size)-dd/2)
diaz_size_max=10^(log10(diaz_size)+dd/2)
diaz_size_range=paste0('[',round(diaz_size_min,3),',',round(diaz_size_max,3),')')

diat_sizes_log=seq(from=1.468, to=2.1345, len=nDIAT)
dd=diff(diat_sizes_log)[1]
diat_sizes_min=10^(diat_sizes_log-dd/2)
diat_sizes_max=10^(diat_sizes_log+dd/2)
diat_size_range=paste0('[',round(diat_sizes_min,3),',',round(diat_sizes_max,3),')')


mp_sname <- str_c("mp", seq_len(nMP))
diat_sname <- str_c("diat", seq_len(nDIAT))

mp_lname <- str_c("Mixed Phytoplankton ", seq_len(nMP))
diat_lname <- str_c("Diatoms ", seq_len(nDIAT))


ss <- data.frame(ESD_um=c(6.16667, mp_sizes, diat_sizes),
        sname=c("diaz", "pp", mp_sname, diat_sname),
        lname=c("Diazotroph", "Picoplankton", mp_lname, diat_lname))

ss$sname <- factor(ss$sname, levels=c("diaz", "pp", mp_sname, diat_sname))
ss$index <- as.numeric(ss$sname)
ss <- arrange(ss, sname)

ss$size_range_um=c(diaz_size_range, mp_size_range, diat_size_range)

# set nitrogen fixers, implicit calcifiers, silicifier flags
ss$Nfixer <- ifelse(ss$sname=="diaz", T, F)
ss$imp_calcifier <- ifelse(ss$ESD_um > 2.0 & ss$ESD_um < 25.0 & str_detect(ss$sname, "mp"), T, F)
ss$exp_calcifier <- F
ss$silicifier <- ifelse(str_detect(ss$sname, "diat"), T, F)


# --- calculate mass from ESD ---
ss$cell_rad_um <- ss$ESD_um / 2
ss$surf_area_um2 <- ss$cell_rad_um^2 * 4 * pi
ss$vol_um3 <- (4/3) * pi * ss$cell_rad_um^3


# mass <- 0.47 * vol_um3^0.99 * 10^-6 # reynolds 2006, for small nanoplankton
ss[ss$ESD_um < 10,"mass_ugC"] <- 0.47 * ss[ss$ESD_um < 10,"vol_um3"]^0.99 * 10^-6

# picoplankton (Bertilsson)
ss[ss$sname=="pp","mass_ugC"] = (ss[ss$sname=="pp","vol_um3"] * 250E-15) / 10^-6


# diatoms, Menden-Deuer and Lessard 2000
ss[str_detect(ss$sname, "diat"),"mass_ugC"] = 10^-6 * 
	10^(-0.933 + 0.881 * log10(ss[str_detect(ss$sname, "diat"),"vol_um3"]))
# small diatoms (< 3000 um3 vol)
ss[str_detect(ss$sname, "diat") & ss$vol_um3 < 3000,"mass_ugC"] = 10^-6 * 
	10^(-0.541 + 0.811 * log10(ss[str_detect(ss$sname, "diat") & ss$vol_um3 < 3000,"vol_um3"]))

# protists
ss[str_detect(ss$sname, "mp") & ss$ESD_um >= 10,"mass_ugC"] = 10^-6 * 
	10^(-0.665 + 0.939 * log10(ss[str_detect(ss$sname, "mp") & ss$ESD_um >= 10,"vol_um3"]))
ss[str_detect(ss$sname, "mp") & ss$ESD_um >= 10 & ss$vol_um3 < 3000,"mass_ugC"] = 10^-6 * 
	10^(-0.583 + 0.860 * log10(ss[str_detect(ss$sname, "mp") & ss$ESD_um >= 10 & ss$vol_um3 < 3000,"vol_um3"]))

ss$mass_umolC <- ss$mass_ugC / 12.011


# --- temperature and growth rate ---
PCref <- data.frame(sname=c("pp", "diaz", "mp1"), PCref_const=c(PCref_const_pp, PCref_const_diaz, PCref_const_mp1))
PCref <- rbind(PCref, data.frame(sname=mp_sname[2:length(mp_sname)], PCref_const=PCref_const_mp),
		data.frame(sname=diat_sname, PCref_const=PCref_const_diat))
PCref$sname <- factor(PCref$sname, levels=c("diaz", "pp", mp_sname, diat_sname))


ss$Ea <- ifelse(ss$sname=="pp", Ea_pp, Ea)
ss$PCref_beta <- ifelse(ss$sname %in% c("pp", "mp1"), PCref_beta_pp_mp1, PCref_beta)
ss[str_detect(ss$sname, "diat"),"PCref_beta"] <- PCref_beta_diat

ss <- left_join(ss, PCref, by="sname")

ss$PCref_per_day <- exp(ss$PCref_beta * log(ss$mass_ugC) - ss$Ea * (1/(K*(Tref+273.15))) + ss$PCref_const)


# --- photosynthesis parameters ---
ss$alphaPI_per_day <- alphaPI_const * ss$vol_um3^(alphaPI_beta)
# different scaling for diatoms
ss[str_detect(ss$sname,"diat"),"alphaPI_per_day"] <- alphaPI_diat_const * ss[str_detect(ss$sname,"diat"),]$vol_um3 ^ (alphaPI_diat_beta)

#ss[ss$sname=="pp","alphaPI_per_day"] <- alphaPI_pp_const


ss$thetaC <- thetaC_const_lg*(ss$vol_um3)^thetaC_beta_lg #units of mg Chla / mmol C

#ss[ss$ESD_um < 5, "thetaC"] = thetaC_coeff_sm*(ss[ss$ESD_um < 5, "vol_um3"]^thetaC_beta_sm)
#ss[str_detect(ss$sname, "diat"),"thetaC"] <- thetaC_coeff_diat * ss[str_detect(ss$sname, "diat"),"vol_um3"]^thetaC_beta_diat

#ss$thetaC <- ss$alphaPI_per_day*thetaC_alpha_beta + thetaC_alpha_const
#ss[str_detect(ss$sname, "diat"),"thetaC"] <- ss[str_detect(ss$sname, "diat"),"alphaPI_per_day"]*thetaC_alpha_beta_diat + thetaC_alpha_const_diat


ss[ss$sname == "pp","thetaC"] = thetaC_pp

ss$thetaN_max <- ss$thetaC * (117/16) # convert to units of mg Chla / mmol N

# override
#ss$thetaN_max[str_detect(ss$sname, 'p')] = 1.25 # pp/mp = sp
#ss$thetaN_max[str_detect(ss$sname, 'diat')] = 2.22 # diat
#ss$thetaN_max[str_detect(ss$sname, 'diaz')] = 1.9 # diaz

tmp = ss[,c("sname", "mass_ugC", "thetaC", "thetaN_max")]
tmp$thetaC_gg = tmp$thetaC / 12.011
print(tmp)

# --- nutrient utilization traits ---
ss$kNO3 <- kNO3_const*ss$vol_um3^kNO3_beta
ss$kPO4 <- kPO4_const*ss$vol_um3^kPO4_beta

ss[str_detect(ss$sname, "diat"),]$kNO3 = kNO3_const_diat*ss[str_detect(ss$sname, "diat"),"vol_um3"]^kNO3_beta_diat
ss[str_detect(ss$sname, "diat"),]$kPO4 = kPO4_const_diat*ss[str_detect(ss$sname, "diat"),"vol_um3"]^kPO4_beta_diat

ss[str_detect(ss$sname, "diaz"),]$kNO3 = kNO3_diaz

# nutrient half-sat constant scalings
ss$kSiO3 <- 0
ss[str_detect(ss$sname, "diat"),]$kSiO3 <- kSiO3_const * ss[str_detect(ss$sname, "diat"),"vol_um3"] ^ 0.3

ss$kNH4 <- kNH4_const * ss$vol_um3 ^ nut_beta_generic
ss[str_detect(ss$sname, "diat"),]$kNH4 = kNH4_const_diat * ss[str_detect(ss$sname, "diat"),"vol_um3"] ^ nut_beta_diatoms

ss[str_detect(ss$sname, "diaz"),]$kNH4 = kNH4_diaz

ss$kDOP <- kDOP_const * ss$vol_um3 ^ nut_beta_generic
ss[str_detect(ss$sname, "diat"),]$kDOP = kDOP_const_diat * ss[str_detect(ss$sname, "diat"),"vol_um3"] ^ nut_beta_diatoms


# --- stochiometry ---
ss$Qp_fixed <- Qp_const * ss$mass_ugC ^ Qp_beta
#ss$Qp_fixed = 1/117.0 # override

ss[ss$sname=='diaz', "Qp_fixed"] <- Qp_fixed_diaz
ss[ss$sname=='pp', "Qp_fixed"] <- Qp_fixed_pp


ss$gQfe_0 <- ifelse(ss$sname=="diaz", gQfe_0_diaz, gQfe_0)
ss$gQfe_min <- gQfe_min

ss$kFe <- kFe_const * ss$vol_um3 ^ kFe_beta
ss[str_detect(ss$sname, "diat"),]$kFe <- kFe_const_diat * ss[str_detect(ss$sname, "diat"),"vol_um3"] ^ kFe_beta_diat
ss[str_detect(ss$sname, "diaz"),]$kFe <- ss[str_detect(ss$sname, "diaz"),]$kFe * diaz_Fe_lim_factor


# loss scalings
ss$loss_thres <- ifelse(ss$sname=="diaz" | str_detect(ss$sname, "diat"), 0.02, 0.01)
ss$loss_thres2 <- ifelse(ss$sname=="diaz" , 0.001, 0)

ss$temp_thres <- ifelse(ss$sname=="diaz", 15, -10)

# mortality exponent
#ss$mort_per_day <- mort_const * ss$mass_ugC ^ mort_beta
ss$mort_per_day <- ss$PCref_per_day * mort_PCref_factor
ss$mort2_per_day <- mort2_const * ss$mass_ugC ^ mort2_beta

#adding a differnt mortality scaling for diatioms
ss[str_detect(ss$sname, "diat"),"mort_per_day"] <- ss[str_detect(ss$sname, "diat"),"PCref_per_day"] * mort_PCref_diat_factor


## attmepting to do mortality for pp to be 0, and increase their mort per day
#ss[ss$sname=='pp', "mort2_per_day"] <- mort2_pp_const
#ss[ss$sname=='pp', "mort_per_day"] <- mort_pp_const
#ss[ss$sname=='diat1', "mort_per_day"] <- mort_diat1_const
#ss[ss$sname=='diat2', "mort_per_day"] <- mort_diat2_const
#ss[ss$sname=='diat3', "mort_per_day"] <- mort_diat3_const

ss$agg_rate_min <- agg_min_const * ss$mass_ugC ^ agg_min_beta
ss$agg_rate_max <- agg_max_const * ss$mass_ugC ^ agg_max_beta

ss[str_detect(ss$sname, "diat"),"agg_rate_min"] <- ss[str_detect(ss$sname, "diat"),"agg_rate_min"] * diat_agg_rate_min_multiplier

ss$loss_poc <- loss_poc_const * ss$mass_ugC ^ loss_poc_beta
ss[str_detect(ss$sname, "diat"),"loss_poc"] <- ss[str_detect(ss$sname, "diat"),"loss_poc"] * 3
ss[ss$imp_calcifier,"loss_poc"] <- ss[ss$imp_calcifier,"loss_poc"]* 2


namelist_vars_autotrophs <- c("sname", "lname", "Nfixer", "imp_calcifier", "exp_calcifier", 
"silicifier", "kFe", "kPO4", "kDOP", "kNO3", "kNH4", "kSiO3", 
"Qp_fixed", "gQfe_0", "gQfe_min", "alphaPI_per_day", "PCref_per_day", 
"thetaN_max", "loss_thres", "loss_thres2", "temp_thres", "mort_per_day", 
"mort2_per_day", "agg_rate_max", "agg_rate_min", "loss_poc", 
"Ea")

numeric_vars <- c("kFe", "kPO4", "kDOP", "kNO3", "kNH4", "kSiO3", "gQfe_0", 
"alphaPI_per_day", "PCref_per_day", "thetaN_max", "mort_per_day", 
"mort2_per_day", "agg_rate_max", "agg_rate_min", "loss_poc")

all(namelist_vars_autotrophs %in% names(ss))
# print(str(ss[,namelist_vars_autotrophs]))


####################################################
################## Zooplankton #####################
####################################################

# --- set up zooplankton ---
nZ = 6
#zoo_sizes <- 10^(seq(from=-1.69897, to=1.30103, len=nZ)) # mm ESD
#zoo_sizes <- 10^(seq(from=-1.82391, to=1.30103, len=nZ)) # mm ESD
zoo_sizes_log = seq(from=-1.449, to=1.0510273, len=nZ)
zoo_sizes = 10^(zoo_sizes_log) # mm ESD

# add in min and max zoo sizes within each bin
dz = diff(zoo_sizes_log)[1]
zoo_sizes_min = 10^(zoo_sizes_log-dz/2)
zoo_sizes_max = 10^(zoo_sizes_log+dz/2)

zoo_size_range = paste0('[',round(zoo_sizes_min,5),',',round(zoo_sizes_max,5),')')

sz = data.frame(ESD_mm=zoo_sizes, 
	sname=str_c('zoo', seq_len(nZ)),
        size_range=zoo_size_range, stringsAsFactors=FALSE)

nMicro = nrow(sz[sz$ESD_mm < 0.2,])
nMeso = nrow(sz[sz$ESD_mm >= 0.2,])

sz$lname = c(str_c("Microzooplankton ", seq_len(nMicro)), str_c("Mesozooplankton ", seq_len(nMeso)))

# biovolume
sz$vol_mm3 <- 4/3 * pi * (sz$ESD_mm / 2) ^ 3
vol_to_ESD <- 2*(3/4 * 1/pi)^(1/3)

sz$mass_mgC <- ifelse(sz$ESD_mm > 0.1,
0.06281 * sz$ESD_mm ^ 3.0, # Pitt et al. 2013, comparative physiology - for mesozooplankton
(0.22 * (sz$vol_mm3*1E9) ^ 0.939) * 1E-9) # Menden-Deuer and Lessard 2000 - ciliates and other protists

sz$Ea <- zoo_Ea

sz$z_mort_0_per_day <- z_mort_const * sz$ESD_mm ^ z_mort_beta
## this is equivalent to:
# sz$z_mort_0_per_day <- z_mort_const * (vol_to_ESD^z_mort_beta) * (sz$vol_mm3 ^ (z_mort_beta*1/3))

#quadratic mortality
sz$z_mort2_0_per_day <- z_mort2_const * sz$ESD_mm ^ z_mort2_beta
## this is equivalent to:
# sz$z_mort2_0_per_day <- z_mort2_const * (vol_to_ESD^z_mort2_beta) * (sz$vol_mm3 ^ (z_mort2_beta*1/3))
sz[nZ,"z_mort2_0_per_day"] <- sz[nZ,"z_mort2_0_per_day"] * lg_zoo_mort2_scaling


sz$loss_thres <- zoo_loss_thres

sz$z_umax_base <- umax_const * sz$ESD_mm ^ umax_beta
## this is equivalent to:
# sz$z_umax_base <- umax_const * (vol_to_ESD^umax_beta) * (sz$vol_mm3 ^ (umax_beta*1/3))

sz$z_grz <- km_const * sz$ESD_mm ^ km_beta

# routing of grazing to POC (equivalent to egestion)
# make sure it is less than assimilation efficiency
sz$graze_poc <- graze_poc_const * sz$ESD_mm^graze_poc_beta
sz[sz$graze_poc > (1-AE),"graze_poc"] = 1-AE

# change to fixed fraction of non-POC detritus
sz$graze_doc <- graze_doc_frac * (1-AE-sz$graze_poc)

# POC routing
sz$f_zoo_detr <- zoo_detr_const * sz$ESD_mm ^ zoo_detr_beta

# --- set up grazing ---
g = array(NA, dim=c(nP + nZ, nZ))

pred_size <- sz$ESD_mm
prey_size <- c(ss$ESD_um / 1000, sz$ESD_mm)

pred_names <- as.character(sz$sname)
prey_names <- c(as.character(ss$sname), as.character(sz$sname))

# set grazing pred-prey size ratio matrix
for (n in 1:length(pred_size)){
	for (m in 1:length(prey_size)){
		pred_prey_ratio = pred_size[n] / prey_size[m]
		g[m,n] <- pred_prey_ratio
	}
}

# do not allow for self grazing 
# or grazing when the predator-prey size ratio < 0.4 or > 30
for (m in 1:(nP+nZ)){
	if (m >= nP+1){
		g[m, m-nP] <- NA
		#g[m,] <- g[m,] / 2
	}
	tmp <- g[m,]
	g[m,][which(tmp < 1 | tmp > 50)] <- NA	
}

# --- grazing preference ---
opt_pred_prey_ratio = opt_pred_prey_ratio_const * sz$ESD_mm ^ opt_pred_prey_ratio_beta
pred_prey_sd = opt_pred_prey_ratio / 2

gf <- array(NA, dim=c(nP + nZ, nZ))

for (z in 1:nZ){
	gf[,z] = pnorm(abs(g[,z]-opt_pred_prey_ratio[z]), 
	mean=0, sd=pred_prey_sd[z], lower=FALSE) * 2
}

### Add in exceptions to the default size-based grazing preference
# Diatoms scaling (highest grazing for small diatoms, then scaling down for larger diatoms)
# trying to do this in a way that makes it not sensitive to nDIAT=3
idx = which(str_detect(ss$sname, 'diat'))
# first third
gf[idx[1:floor(nDIAT/3)],] = gf[idx[1:floor(nDIAT/3)],] * diat_grazing_scaling
# second third
gf[idx[ceiling(nDIAT/3):floor(nDIAT*2/3)],] = gf[idx[ceiling(nDIAT/3):floor(nDIAT*2/3)],] * med_diat_grazing_scaling
# final third
gf[idx[ceiling(nDIAT*2/3):nDIAT],] = gf[idx[ceiling(nDIAT*2/3):nDIAT],] *lg_diat_grazing_scaling

# Diatom grazing by microzooplankton
idy = seq_len(nMicro)
gf[idx,idy] = gf[idx,idy] * microZ_diat_grazing

# picoplankton grazing by microzooplankton
# fixed value, not scaling
idx = which(ss$sname == "pp")
gf[idx,idy] = (gf[idx,idy] / gf[idx,idy] ) * microZ_pp_grazing_fixed

# scaling of implicit calcifiers
gf[which(ss$imp_calcifier),] = gf[which(ss$imp_calcifier),]  * calcifiers_grazing_scaling

# scaling of zooplankton feeding on zooplankton
gf[(nP+1):(nP+nZ),] = gf[(nP+1):(nP+nZ),] * zoo_grazing_scaling

# scaling of predation on microzooplankton
gf[(nP+1):(nP+nMicro),] = gf[(nP+1):(nP+nMicro),] * microzoo_grazing_scaling

# set zeros to be NAs to avoid making excess computation
gf[which(gf == 0)] = NA

### normalize grazing preference by column
gf <- apply(gf, 2, function(x){ 
	x <- x / sum(x, na.rm=TRUE)
	return(x)
})	


# scale grazing umax by grazing function
g_z_umax <- rep(sz$z_umax_base, each=(nP+nZ)) * gf


# --- grazing array ---
max_grazer_prey_cnt = max(apply(g, 2, function(x) length(which(!is.na(x)))))

sg = data.frame()

for (n in 1:dim(g)[2]){
	k=0
	
	for (m in 1:dim(g)[1]){
		if(!(is.na(g[m,n]))){
			k=k+1
			
			if (m <= nP) {	auto_ind_cnt=1; zoo_ind_cnt=0
							auto_ind=m; zoo_ind=NA}
			if (m > nP) {	auto_ind_cnt=0; zoo_ind_cnt=1
							auto_ind=NA; zoo_ind=(m-nP)}
			
			tmp <- data.frame(type="grazing", index1 = k, index2 = n, 
						sname1 = prey_names[m], sname2=pred_names[n], 
						sname=str_c("grz_", prey_names[m], "_", pred_names[n]),
						lname=str_c("Grazing of ", prey_names[m], " by ", pred_names[n]),
						auto_ind_cnt, zoo_ind_cnt, auto_ind, zoo_ind,
						z_umax_0_per_day = g_z_umax[m,n], stringsAsFactors=FALSE)
			
			sg <- rbind(sg, tmp)
		}
		
	}
}


# set grazing function
#sg$grazing_function <- ifelse(sg$sname1 %in% ss$sname[1:4], 1, 2)
sg$grazing_function[str_detect(sg$sname1,"zoo")] <- 1

# join together max grazing rates from zooplankton df
sg <- left_join(sg, data.frame(sname2=sz$sname, 
	z_grz=sz$z_grz, graze_doc=sz$graze_doc,
	graze_poc=sz$graze_poc,
    f_zoo_detr=sz$f_zoo_detr, stringsAsFactors=FALSE), by="sname2")

# set routing
# change graze_poc amount for diatom prey and picoplankton prey
sg[str_detect(sg$sname1, "diat"), "graze_poc"] = sg[str_detect(sg$sname1, "diat"), "graze_poc"] * graze_poc_diatom_scaling 
sg[sg$sname1=="pp", "graze_poc"] = 0.0


sg$graze_zoo <- ifelse(sg$sname2 %in% sz[sz$ESD_mm < 0.2,"sname"], graze_zoo_micro, graze_zoo_meso)
#sg$f_zoo_detr <- ifelse(sg$sname2=="zoo1", 0.12, 0.24)

# check to make sure grazing routing terms are less than 1
idx = which(sg$graze_doc + sg$graze_poc + sg$graze_zoo > 1)
sg[idx,c('graze_doc','graze_poc')] <- sg[idx,c('graze_doc','graze_poc')] / (rowSums(sg[idx,c('graze_doc','graze_poc')])/(1-1e-10-sg$graze_zoo[idx]))

# increase fraction of f_zoo_detr when prey are diatoms
sg[str_detect(sg$sname1, "diat"), "f_zoo_detr"] <- sg[str_detect(sg$sname1, "diat"), "f_zoo_detr"] * zoo_detr_diatom_scaling

sg <- dplyr::rename(sg, auto_ind.1.=auto_ind, zoo_ind.1.=zoo_ind)


# --- fill in the grazing array ---
tmp_g <- data.frame(type="grazing", index1=rep(1:max_grazer_prey_cnt, times=nZ),
                                        index2=rep(1:nZ, each=max_grazer_prey_cnt), stringsAsFactors=FALSE)

sg <- left_join(tmp_g, sg, by=c("type", "index1", "index2"))


# full in NAs with nulls and zeros
sg[which(is.na(sg$sname1)),c("sname1", "sname2", "sname", "lname")] <- "null"

sg[which(is.na(sg$auto_ind_cnt)),
        c("auto_ind_cnt", "zoo_ind_cnt", "z_umax_0_per_day", "z_grz", "graze_zoo",
        "graze_poc", "graze_doc", "f_zoo_detr")] <- 0

# grazing function must be either 1 or 2, can't be 0
sg[which(is.na(sg$grazing_function)), "grazing_function"] <- 1



####################################################
############# TUNING & WRITING FILES ###############
####################################################

# --- phytoplankton ---
ss$PCref_per_day <- ss$PCref_per_day * phyto_PCref_modification
#ss$gQfe_0 = 3.5E-5


# --- save sizes ---
tmp_s = ss[,c("sname", "mass_ugC", "vol_um3", "ESD_um", "Qp_fixed")]
tmp_s$ESD_mm = tmp_s$ESD_um / 1000
tmp_s$type = "phyto"

tmp_z = sz[,c("sname", "mass_mgC", "vol_mm3", "ESD_mm")]
tmp_z$Qp_fixed = 1/117.0
tmp_z$mass_ugC = tmp_z$mass_mgC * 1000
tmp_z$vol_um3 = tmp_z$vol_mm3 * 1e9
tmp_z$type = "zoo"

cols = c("type", "sname", "mass_ugC", "vol_um3", "ESD_mm", "Qp_fixed")
sizes = rbind(tmp_s[,cols], tmp_z[,cols])
sizes$mmolC = sizes$mass_ugC / 12.011 / 1000


# --- save sizes ---
tmp_s = ss[,c("sname", "mass_ugC", "vol_um3", "ESD_um", "Qp_fixed")]
tmp_s$ESD_mm = tmp_s$ESD_um / 1000
tmp_s$type = "phyto"

tmp_z = sz[,c("sname", "mass_mgC", "vol_mm3", "ESD_mm")]
tmp_z$Qp_fixed = 1/117.0
tmp_z$mass_ugC = tmp_z$mass_mgC * 1000
tmp_z$vol_um3 = tmp_z$vol_mm3 * 1e9
tmp_z$type = "zoo"

cols = c("type", "sname", "mass_ugC", "vol_um3", "ESD_mm", "Qp_fixed")
sizes = rbind(tmp_s[,cols], tmp_z[,cols])
sizes$mmolC = sizes$mass_ugC / 12.011 / 1000


# --- write files ---
ss <- cbind(type="autotroph", ss) 

sz <- cbind(type="zooplankton", index=1:nZ, sz)
sz$sname <- as.character(sz$sname)

# set namelist vars order
namelist_vars_grazing <- c("sname", "lname", "auto_ind_cnt", "zoo_ind_cnt", "grazing_function", 
"z_umax_0_per_day", "z_grz", "graze_zoo", "graze_poc", "graze_doc", 
"f_zoo_detr", "auto_ind.1.", "zoo_ind.1.")

all(namelist_vars_grazing %in% names(sg))


dir.create("data", showWarnings=FALSE)
write.csv(sizes, "data/plankton_sizes.csv", row.names=FALSE)
write.csv(ss, "data/phytoplankton_input_data.csv", row.names=FALSE)
write.csv(sz, "data/zooplankton_input_data.csv", row.names=FALSE)
write.csv(sg, "data/grazing_input_data.csv", row.names=FALSE)
write.csv(gf, "data/grazing_FK_value.csv", row.names=FALSE)



