#!/usr/bin/env Rscript

#setwd("/Users/jluo/Dropbox/Research/NCAR/models/MARBL-gen_inputs/gen-spectra")

library("ggplot2")
library("plyr")
library("reshape2")
library("stringr")
library("RColorBrewer")
library("scales")

# load data
ss <- read.csv("data/phytoplankton_input_data.csv", stringsAsFactors=FALSE)
sz <- read.csv("data/zooplankton_input_data.csv", stringsAsFactors=FALSE)
sg <- read.csv("data/grazing_input_data.csv", stringsAsFactors=FALSE)
gf <- read.csv("data/grazing_FK_value.csv")

dir.create("plots", showWarnings=FALSE)

# some initial plots
ss$type <- c("diaz", rep("mp",5), rep("diat",3))
sz$type <- 'zoo'

p <- ggplot() + 
	geom_point(aes(y=log10(mass_ugC / 1000), x=log10(ESD_um / 1000), color=type, shape=type), data=ss) +
  geom_line(aes(y=log10(mass_ugC / 1000), x=log10(ESD_um / 1000), linetype=type), data=ss) +
	geom_text(aes(y=log10(mass_ugC / 1000)-0.2, x=log10(ESD_um / 1000) + 0.2, label=sname, color=type), alpha=0.5, data=ss, show.legend=FALSE) + 
	geom_point(aes(y=log10(mass_mgC), x=log10(ESD_mm), color=type, shape=type), data=sz) + 
  geom_line(aes(y=log10(mass_mgC), x=log10(ESD_mm), linetype=type), data=sz) + 
	geom_text(aes(y=log10(mass_mgC)+0.7, x=log10(ESD_mm), label=sname, color=type), alpha=0.5, data=sz, show.legend=FALSE) + 
  scale_linetype_manual('PFT Type', values=c(4,0,1,3), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto", "Zooplankton")) +
  scale_shape_manual('PFT Type', values=c(15,17,19,9), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto", "Zooplankton")) +
	scale_color_manual('PFT Type', values=c("#984ea3", "#4daf4a", "#377eb8", "#e41a1c"), 
	labels=c("Diatoms", "Diazotrophs", "Mixed Phyto", "Zooplankton")) + 
	labs(x=expression(paste(log[10], " Equivalent Spherical Diameter (ESD, mm)")), y=expression(paste(log[10]," Mass (mg C)"))) + theme_classic()

ggsave("plots/PFT_mass_v_ESD.png", p, height=6, width=6)

# stochiometry & growth 
d1 = melt(ss[,c("type", "mass_ugC", "ESD_um", "Qp_fixed", "PCref_per_day", "alphaPI_per_day", "thetaN_max")], id.vars=c("type", "mass_ugC", "ESD_um"))

p1 <- ggplot(d1, aes(x=log10(ESD_um/1000), y=log10(value))) + 
  geom_point(aes(color=type, shape=type)) + 
  geom_line(aes(linetype=type)) + 
  theme_classic() + 
  facet_grid(variable~., scales="free_y") +
  scale_shape_manual('PFT Type', values=c(15,17,19), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) +
	scale_linetype_manual('PFT Type', values=c(4,0,1), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) + 
  scale_color_manual("PFT Type", values=c("#984ea3", "#4daf4a", "#377eb8"), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) +
  labs(x=expression(paste(log[10], " Equivalent Spherical Diameter (ESD, mm)")), y=expression(paste(log[10]," Parameter value"))) 

ggsave("plots/PFT_growth_parameters.png", p1, height=6, width=6)

d1$value2 = 1/d1$value
p1 = ggplot(mapping=aes(x=log10(ESD_um/1000), y=log10(value2))) + 
  geom_point(aes(shape=type), data=d1[d1$variable=="Qp_fixed",]) + 
	geom_line(linetype=2, data=d1[d1$variable=="Qp_fixed" & d1$type != "diaz",]) + theme_bw() + 
	scale_shape_manual('PFT Type', values=c(19,8,17)) + labs(x="Log10 Equiv. Spherical Diameter (mm)", y="Log10 P:C ratio", colour="Parameter")
ggsave("plots/PFT_stochiometry.png", p1, height=2.5, width=4.5)


# nutrient uptake
d2 = melt(ss[,c("type", "mass_ugC", "ESD_um", "kNO3", "kNH4", "kPO4", "kDOP", "kSiO3", "kFe", "gQfe_0")], id.vars=c("type", "mass_ugC", "ESD_um"))

p2 <- ggplot(d2[d2$value != 0,], aes(x=log10(ESD_um/1000), y=log10(value))) + 
  geom_point(aes(color=variable, shape=type)) + 
  geom_line(aes(linetype=type, color=variable)) + 
  theme_classic() + 
  scale_color_brewer("Parameter", palette="Paired") + 
  scale_shape_manual('PFT Type', values=c(15,17,19), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) +
  scale_linetype_manual('PFT Type', values=c(4,0,1), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) + 
  labs(x=expression(paste(log[10], " Equivalent Spherical Diameter (ESD, mm)")), y=expression(paste(log[10]," Parameter value"))) 

ggsave("plots/PFT_nutrient_parameters.png", p2, height=6, width=6)

p2b <- ggplot(d2[d2$value != 0,], aes(x=log10(ESD_um/1000), y=log10(value))) + 
  geom_point(aes(color=type, shape=type)) + 
  geom_line(aes(linetype=type)) + 
  theme_classic() + 
  facet_grid(variable~., scales="free_y") +
  scale_color_manual("PFT Type", values=c("#984ea3", "#4daf4a", "#377eb8"), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) +
  scale_shape_manual('PFT Type', values=c(15,17,19), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) +
  scale_linetype_manual('PFT Type', values=c(4,0,1), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) + 
  labs(x=expression(paste(log[10], " Equivalent Spherical Diameter (ESD, mm)")), y=expression(paste(log[10]," Parameter value"))) 

ggsave("plots/PFT_nutrient_parameters_b.png", p2b, height=6, width=6)

# loss terms
d3 = melt(ss[,c("type", "mass_ugC", "ESD_um", "mort_per_day", "mort2_per_day", "agg_rate_min", "agg_rate_max", "loss_poc")], id.vars=c("type", "mass_ugC", "ESD_um"))

p3 <- ggplot(d3, aes(x=log10(ESD_um/1000), y=log10(value))) + 
  geom_point(aes(color=type, shape=type)) + 
  geom_line(aes(linetype=type)) + 
  theme_classic() + 
  facet_grid(variable~., scales="free_y") +
  scale_color_manual("PFT Type", values=c("#984ea3", "#4daf4a", "#377eb8"), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) +
  scale_shape_manual('PFT Type', values=c(15,17,19), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) +
  scale_linetype_manual('PFT Type', values=c(4,0,1), labels=c("Diatoms", "Diazotrophs", "Mixed Phyto")) + 
  labs(x=expression(paste(log[10], " Equivalent Spherical Diameter (ESD, mm)")), y=expression(paste(log[10]," Parameter value"))) 

ggsave("plots/PFT_loss_parameters.png", p3, height=6, width=6)


# grazing
d4 = melt(sz[,c("ESD_mm", "mass_mgC", "z_umax_base", "z_grz", "z_mort_0_per_day", "z_mort2_0_per_day", "graze_doc", "f_zoo_detr")], id.vars=c("ESD_mm", "mass_mgC"))

p4 <- ggplot(d4, aes(x=log10(ESD_mm), y=log10(value), color=variable)) + 
  geom_point(shape=9) + 
  geom_line() + 
  theme_classic() + 
  scale_color_brewer("Parameter", palette="Paired") + 
  labs(x=expression(paste(log[10], " Equivalent Spherical Diameter (ESD, mm)")), y=expression(paste(log[10]," Parameter value"))) 

ggsave("plots/Zoo_parameters.png", p4, height=6, width=6)


sz$ESD_um = sz$ESD_mm *1000
s = rbind(ss[,c("sname", "ESD_um")], sz[,c("sname", "ESD_um")])
names(gf) <- sz$sname
s <- arrange(s, ESD_um)


gf <- cbind(prey=c(as.character(ss$sname), sz$sname), gf)
gf_melt <- melt(gf, id.vars="prey",variable.name ='pred')
gf_melt$prey <- factor(gf_melt$prey, levels=rev(s$sname))

p5 <- ggplot(gf_melt) + geom_tile(aes(y=prey, x=pred, fill=value)) + 
	scale_fill_gradient2(low="#2166ac", high="#b2182b", mid="#f7f7f7", midpoint=0.3, limits=c(0,0.6)) +
	labs(x="Predators", y="Prey", fill="FK value") + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))
	
ggsave("plots/Zoo_FK_grid.pdf", p5, height=8, width=6)

# grazing functions
P = seq(0,40,by=0.5)
dg = ddply(sg[sg$sname!="null",], ~sname1+sname2, function(x){
  umax=x$z_umax_0_per_day
  km=x$z_grz
  g=(umax * P^2) / (P^2 + km^2)
  return(data.frame(umax,km,P,g))
})
dg$grazing = str_c('graze_',dg$sname1,'_',dg$sname2)
ggplot(dg) + geom_path(aes(x=P,y=g,color=sname1)) + facet_wrap(~sname2) + theme_classic()

gf_umax = sg[sg$sname!="null", c("sname1", "sname2", "z_umax_0_per_day")]
names(gf_umax) = c("prey", "pred", "umax")
gf_nas = gf_melt[is.na(gf_melt$value),]
gf_umax <- rbind(gf_umax, data.frame(cbind(gf_nas[,c("prey","pred")], umax=gf_nas$value)))
gf_umax$prey = factor(gf_umax$prey, levels=rev(s$sname))
gf_umax$pred = factor(gf_umax$pred, levels=sz$sname)


p6 <- ggplot(gf_umax) + geom_tile(aes(y=prey, x=pred, fill=umax)) + 
  geom_text(aes(y=prey, x=pred, label=ifelse(gf_umax$umax < 0.1, sprintf('%0.1e', gf_umax$umax), sprintf('%0.1f', gf_umax$umax)))) +
  scale_fill_gradientn(colors=rev(brewer.pal(11,"Spectral")), na.value = 'grey50', limits=c(0,0.6)) + 
  labs(x="Predators", y="Prey", fill="Maximum\ngrazing rate") + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))

ggsave("plots/Zoo_umax_grid.pdf", p6, height=7, width=6)

# feeding kernel
opt_pred_prey_ratio = 10
opt_pred_prey_sd_const = 9
opt_pred_prey_sd_beta = 0.1

sz$pred_prey_sd = opt_pred_prey_sd_const * sz$ESD_mm ^ opt_pred_prey_sd_beta

x=seq(0.1,30, by=0.1)

#plot(pnorm(x, mean=10, sd=sz$pred_prey_sd[1])~x, type='l')
y = pnorm(abs(x-opt_pred_prey_ratio), mean=0, sd=sz$pred_prey_sd[1], lower=FALSE) * 2
pdf("plots/dist_spectra.pdf",height=4,width=6)
plot(y~x, type='l', ylab="SPECTRA pdf", xlab="Predator-prey size ratio")
dev.off()


y = dnorm(x, mean=opt_pred_prey_ratio, sd=sz$pred_prey_sd[1]) 
pdf("plots/dist_gaussian.pdf",height=4,width=6)
plot(y~x, type='l', ylab="Gaussian pdf", xlab="Predator-prey size ratio")
dev.off()

#require(rmutil)
#pdf("plots/dist_laplace.pdf",height=4,width=6)
#plot(dlaplace(x,m=opt_pred_prey_ratio,s=sz$pred_prey_sd[1])~x, type='l', ylab = "Laplace pdf", xlab="Predator-prey size ratio")
#dev.off()

### -- Code to generate the feeding kernel (before adjustments) -- ### 
fk = data.frame(x= rep(x, times=nrow(sz)))
fk$sname = rep(sz$sname, each=length(x))
#fk$sname = factor(fk$sname)
fk = join(fk, sz[,c("sname","pred_prey_sd")])

fk = ddply(fk, ~sname, function(x){
  x$fk_val = pnorm(abs(x$x-opt_pred_prey_ratio), mean=0, sd=x$pred_prey_sd, lower=FALSE) * 2
  return(x)
})

p7 <- ggplot(fk) + 
  geom_line(aes(x=x, y=fk_val, color=sname)) +
  theme_classic() +
  scale_color_manual("Zooplankton", values=brewer.pal(9,"OrRd")[4:9]) +
  labs(x="Predator-prey size ratio", y="Feeding Kernel Value")

ggsave("plots/Zoo_FK_distribution.pdf", p7, height=5, width=6)

