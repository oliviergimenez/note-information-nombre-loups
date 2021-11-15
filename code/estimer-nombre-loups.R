#-------------------------------------------------------------------------------
# Olivier Gimenez (CNRS) et Chistophe Duchamp (OFB) - novembre 2021
# Estimation du nombre de loups en France par capture-recapture

# Citation :
# Cubaynes, S. Pradel, R. Choquet, R. Duchamp, C. Gaillard, J-M., 
# Lebreton, J-D., Marboutin, E., Miquel, C., Reboulet, A-M., Poillot, C., 
# Taberlet, P. and O. Gimenez. (2010). Importance of accounting for detection 
# heterogeneity when estimating abundance: the case of French wolves. 
# Conservation Biology 24: 621-626.
#-------------------------------------------------------------------------------

#----------- 1. Preliminaires : packages R, donnees, et fonctions utiles

# packages R requis
library(tidyverse)
theme_set(theme_light(base_size = 35))
library(janitor)
library(R2ucare)
library(broom)

# donnees capture-recapture hiver 1997/1998 -> hiver 2018/2019
dat <- readRDS("dat/cmrlouphiver.rds")

# fonctions utiles
source('code/loupfns.R')

#----------- 2. Qualite ajustement modele Cormack-Jolly-Seber

# ignore reprises d animaux morts
dataloup_gof <- as.data.frame(dat)
dataloup_gof[dataloup_gof == 2] <- 0
eff_gof <- rep(1,nrow(dataloup_gof))

# composantes du test d adequation du modele CJS (trap-dep et transient)
test3sr(dataloup_gof,eff_gof)
test3sm(dataloup_gof,eff_gof)
test2ct(dataloup_gof,eff_gof)
test2cl(dataloup_gof,eff_gof)

# test pour modele avec effet transient pris en compte
stat_cjsage <- overall_CJS(dataloup_gof,eff_gof)$chi2 - test3sr(dataloup_gof,eff_gof)$test3sr['stat']
df_cjsage <- overall_CJS(dataloup_gof,eff_gof)$degree_of_freedom - test3sr(dataloup_gof,eff_gof)$test3sr['df']
stat_cjsage 
df_cjsage
1 - pchisq(stat_cjsage, df_cjsage)
# ajustement satisfaisant une fois effet age sur la survie pris en compte
# on considere que effet transient est un melange de transient et de jeunes

#----------- 3. Selection meilleur modele

# Modeles consideres
# - pi: proportion constante
# - phi: survie avec 2 classes d’âge, hom ou het, effet temps ou periode
# - psi: transition avec ou pas
# - p: detection het ou pas
# - r: recovery het ou pas
# 1: pi, phi(a2), p(mix), r(mix), psi
# 2: pi, phi(a2), p, r(mix), psi
# 3: pi, phi(a2), p(mix), r, psi
# 4: pi, phi(a2), p, r, psi
# 5: pi, phi(a2*mix), p(mix), r(mix), psi
# 6: pi, phi(a2*mix), p, r(mix), psi
# 7: pi, phi(a2*mix), p(mix), r, psi
# 8: pi, phi(a2*mix), p, r, psi
# 9: pi, phi(a2), p(mix), r(mix)
# 10: pi, phi(a2), p, r(mix)
# 11: pi, phi(a2), p(mix), r
# 12: pi, phi(a2), p, r
# 13: pi, phi(a2*mix), p(mix), r(mix)
# 14: pi, phi(a2*mix), p, r(mix)
# 15: pi, phi(a2*mix), p(mix), r
# 16: pi, phi(a2*mix), p, r
# 17: pi, phi(a2*t), p(mix), r(mix), psi
# 18: pi, phi(a2*mix*t), p(mix), r(mix), psi
# 19: pi, phi(a2*period), p(mix), r(mix), psi
# nb param: 9, 8, 8, 7, 11, 10, 10, 9, 7, 6, 6, 5, 9, 8, 8, 7, 2*(k-1)+6, 4*(k-1)+5, 11 (k = nb occ capture)

# quantites misc
dataloup_hiver <- as.data.frame(dat)
s <- 5 # nb etats (vivant classe 1, vivant classe 2, juste mort classe 1, juste mort classe 2, mort)
m <- 3 # nb obs (pas detecte, detecte vivant, detecte mort)
k <- dim(dataloup_hiver)[2]
km1 <- k-1
nh <- dim(dataloup_hiver)[1]

# effectifs
eff <- rep(1,nh)

# premiere capture
fc <- NULL
init.state <- NULL
for (jj in 1:nh){
  temp <- 1:k
  fc <- c(fc,min(temp[dataloup_hiver[jj,]==1]))
  init.state <- c(init.state, dataloup_hiver[jj,fc[jj]])
}
	
# derniere capture
lc <- readRDS(file = "dat/lc.rds")

# transpose donnees
dataloup_hiver <- t(dataloup_hiver)
	
# nb parametres
nb_param <- c(9, 8, 8, 7, 11, 10, 10, 9, 7, 6, 6, 5, 9, 8, 8, 7, 2*km1+6, 4*km1+5, 11)
nb_modl <- length(nb_param)

# ajuste 19 modeles avec 3 jeux differents de valeurs initiales 
# pour tenir compte de min locaux eventuels
# les temps de calcul sont longs
res_mdl <- vector("list", nb_modl)
nb_inits <- 3
for (i in 1:nb_modl){
  # considere modele courant
  mod <- paste('devCJS',i,sep='')
  res_tpmin <- vector("list", nb_inits)
  for (j in 1:nb_inits){
    binit <- runif(nb_param[i])
    tpmin <- optim(par = binit,
                   fn = eval(parse(text = mod)),
                   gr = NULL,
                   hessian = FALSE, 
                   dataloup_hiver,eff,fc,lc,init.state,nh,km1,
                   method = "BFGS",
                   control = list(trace = 1, REPORT = 1, maxit = 500))
    res_tpmin[[j]] <- tpmin
  }
  dev_values <- NULL
  for (kk in 1: nb_inits) dev_values <- c(dev_values,res_tpmin[[kk]]$value)
  global.dev <- which.min(dev_values)
  res_mdl[[i]] <- res_tpmin[[global.dev]]
}

# calcule AIC
AIC <- NULL
for (kk in 1: nb_modl){
  AIC <- c(AIC, 2 * res_mdl[[kk]]$value + 2 * nb_param[kk])
} 
res <- data.frame(AIC = AIC,
                  model.name = c('1: pi, phi(a2), p(mix), r(mix), psi',
                                 '2: pi, phi(a2), p, r(mix), psi',
                                 '3: pi, phi(a2), p(mix), r, psi',
                                 '4: pi, phi(a2), p, r, psi',
                                 '5: pi, phi(a2*mix), p(mix), r(mix), psi',
                                 '6: pi, phi(a2*mix), p, r(mix), psi',
                                 '7: pi, phi(a2*mix), p(mix), r, psi',
                                 '8: pi, phi(a2*mix), p, r, psi',
                                 '9: pi, phi(a2), p(mix), r(mix)',
                                 '10: pi, phi(a2), p, r(mix)',
                                 '11: pi, phi(a2), p(mix), r',
                                 '12: pi, phi(a2), p, r',
                                 '13: pi, phi(a2*mix), p(mix), r(mix)',
                                 '14: pi, phi(a2*mix), p, r(mix)',
                                 '15: pi, phi(a2*mix), p(mix), r',
                                 '16: pi, phi(a2*mix), p, r',
                                 '17: pi, phi(a2*t), p(mix), r(mix), psi',
                                 '18: pi, phi(a2*mix*t), p(mix), r(mix), psi',
                                 '19: pi, phi(a2*period), p(mix), r(mix), psi'))

# reordonne selon valeurs croissantes AIC
res.ordered <- res[order(AIC),]
res.ordered

# AIC                                  model.name
# 17 3829.721      17: pi, phi(a2*t), p(mix), r(mix), psi
# 19 3831.733      19: pi, phi(a2*period), p(mix), r(mix), psi
# 5  3848.414       5: pi, phi(a2*mix), p(mix), r(mix), psi
# 18 3853.162      18: pi, phi(a2*mix*t), p(mix), r(mix), psi
# 13 3856.115      13: pi, phi(a2*mix), p(mix), r(mix)
# 1  3856.479       1: pi, phi(a2), p(mix), r(mix), psi
# 9  3869.122       9: pi, phi(a2), p(mix), r(mix)
# 3  3871.769       3: pi, phi(a2), p(mix), r, psi
# 7  3875.308       7: pi, phi(a2*mix), p(mix), r, psi
# 11 3876.958      11: pi, phi(a2), p(mix), r
# 15 3880.387      15: pi, phi(a2*mix), p(mix), r
# 14 3894.766      14: pi, phi(a2*mix), p, r(mix)
# 2  3895.738       2: pi, phi(a2), p, r(mix), psi
# 6  3898.750       6: pi, phi(a2*mix), p, r(mix), psi
# 12 3912.349      12: pi, phi(a2), p, r
# 10 3914.349      10: pi, phi(a2), p, r(mix)
# 4  3916.349       4: pi, phi(a2), p, r, psi
# 16 3916.349      16: pi, phi(a2*mix), p, r
# 8  3919.938       8: pi, phi(a2*mix), p, r, psi

#----------- 4. Parametres estimes

# effet periode sur survie, et heterogeneite sur parametres detection
# on se base sur modele 19 pour estimer les effectifs
# les temps de calcul sont un tout petit peu longs
nb_inits <- 1
label_best_model <- 19
mod <- paste('devCJS', label_best_model,sep='')
res_tpmin <- vector("list", nb_inits)
for (j in 1: nb_inits){
	binit <- runif(nb_param[label_best_model])
	tpmin = optim(binit,eval(parse(text=mod)),NULL, hessian=T, dataloup_hiver,eff,
	fc, lc, init.state,nh,km1,method="BFGS",control=list(trace=1, REPORT=1))
	res_tpmin[[j]] <- tpmin
	}
dev_values <- NULL
for (kk in 1: nb_inits) dev_values <- c(dev_values,res_tpmin[[kk]]$value)
global.dev <- which.min(dev_values)
res_mdl <- res_tpmin[[global.dev]]
x <- res_mdl$par
fisher_info <- solve(0.5 * res_mdl$hessian)
prop_sigma <- sqrt(diag(fisher_info))

# les parametres estimes et intervalles de confiance
estim <- function(bb, prop_sigma){
param <- plogis(bb)
IClower <- plogis(bb - 1.96 * prop_sigma)
ICupper <- plogis(bb + 1.96 * prop_sigma)
list(mle = param, CI = c(IClower, ICupper))
}

# prop d individus de chaque classe
(prop <- estim(x[5], prop_sigma[5]))

# pr detection
(pp1 <- estim(x[1], prop_sigma[1])) 
(pp2 <- estim(x[2], prop_sigma[2]))

# pr reprise
(ll1 <- estim(x[3], prop_sigma[3]))
(ll2 <- estim(x[4], prop_sigma[4])) 

# pr transitions entre classes
(psi12 <- estim(x[6], prop_sigma[6]))
(psi21 <- estim(x[7], prop_sigma[7]))

# pr de survie
(phi1_period1 <- estim(x[8], prop_sigma[8]))
(phi1_period2 <- estim(x[9], prop_sigma[9]))
(phi2_period1 <- estim(x[10], prop_sigma[10]))
(phi2_period2 <- estim(x[11], prop_sigma[11]))

#----------- 5. effectif estime et int de confiance par bootstrap

# les temps de calcul sont tres tres longs
nbMC <- 500 # nb bootstraps
xx <- matrix(NA,nbMC,nb_param[label_best_model]) # store parameters
Nhet <- NA
dataloup_hiver <- t(dataloup_hiver)
for (i in 1:nbMC){
  
  # genere pseudo jeu de donnees
  mask <- sample.int(nh,replace=T)
  pseudodata <- dataloup_hiver[mask,]
  
  # date premier capture
  pseudo_fc <- fc[mask]
  pseudo_init.state <- init.state[mask]
	pseudo_lc <- lc[mask]
	
	# transpose donnees
	pseudodata <- t(pseudodata)

	mod <- paste('devCJS', label_best_model,sep='')
	res_tpmin <- vector("list", nb_inits)
	for (j in 1: nb_inits){
		binit <- x
		tpmin <- optim(par = binit,
		               fn = eval(parse(text=mod)),
		               gr = NULL, 
		               hessian = FALSE,
		               pseudodata, eff, pseudo_fc, pseudo_lc, pseudo_init.state, nh, km1,
		               method = "BFGS",
		               control = list(trace = 1, REPORT = 1))
		res_tpmin[[j]] <- tpmin
	}
	dev_values <- NULL
	for (kk in 1: nb_inits) dev_values <- c(dev_values,res_tpmin[[kk]]$value)
	global.dev <- which.min(dev_values)
	res_mdl <- res_tpmin[[global.dev]]
	xx[i,] <- res_mdl$par

	prop <- plogis(xx[i,5])
	p1 <- plogis(xx[i,1])
	p2 <- plogis(xx[i,2])

	# re-ordonne estimations (label switching)
	if (p1>p2){
		temp <- p2 
		p2 <- p1
		p1 <- temp
		prop <- 1-prop
	}

	phi_period1 <- plogis(xx[i,10])
	phi_period2 <- plogis(xx[i,11])

	# calcule effectifs estimes	
	h <- t(pseudodata)
	
	# remplace les morts, les 2, par des 0 pour ne pas les comptabiliser dans les effectifs
	h[h==2] <- 0
	
	# detectes par occasion
	cijdata <- apply(h,2,sum)
	
	# date de premiere detection
	d <- apply(h,1,premieredetection)
	udata <- NULL
	for(kk in 1:ncol(h)) {udata[kk] <- sum(d==kk)}
	
	# deja marques par occasion
	mdata <- cijdata - udata
	
	# on enleve la premiere occasion
	cij <- cijdata[-1]
	m <- mdata[-1]
	u <- udata[-1]
	
	#--- calcul des nouveaux marques attendus
	bigU <- matrix(0, nrow = length(p1), ncol = length(u))
	for(zz in 1:ncol(bigU)) {
		bigU[,zz] <- (1-prop) * u[zz] / p2 + prop * u[zz] / p1
	}
	# apply(bigU,2,mean)

	#--- calcul des deja marques attendus
	#M2 = u1 (pi phi1 + (1-pi) phi1)
	#M3 = u1 (pi phi1 phi2 + (1-pi) phi1 phi2) + u2 (pi phi2 + (1-pi) phi2)
	#M4 = u1 (pi phi1 phi2 phi3 + (1-pi) phi1 phi2 phi3) + u2 (pi phi2 phi3 + (1-pi) phi2 phi3) + u3 (pi phi3 + (1-pi) phi3)
	#...
	
	survie <- t(as.matrix(c(rep(phi_period1,18), rep(phi_period2, length(19:km1))))) # replique phi par int de temps
	bigM <- matrix(0,nrow = nrow(survie),ncol = ncol(survie))

	for(ii in 1:nrow(bigM)) {
		for(t in 1:ncol(bigM)) {
			temp <- rep(NA,t)
			for(j in 1:t) {
				temp[j] <- u[j] * ((1-prop[ii]) * prod(survie[ii,1:j]) + prop[ii] * prod(survie[ii,1:j]))
			}
		bigM[ii,t] <- sum(temp) 
		}
	}
	# apply(bigM,2,mean)

	#--- calcul des effectifs estimes
	Nhet <- rbind(Nhet,bigU + bigM)

}

# calcule effectifs annuels et intervalles de confiance
Nhet <- Nhet[-1,] # effectif premiere annee non-estimable puisqu'on estime des probs de recapture
Nhet2 <- Nhet
# on prend la mediane de la distribution bootstrap (la moyenne donne des valeurs proches)
apply(Nhet2, 2, quantile, c(2.5,50,97.5)/100, na.rm = TRUE)

#----------- 6.  prediction des effectifs pour 2019/2020 et 2020/2021 
#-----------     via calibration capture-recapture - effectifs minimums retenus

#-- les donnees de capture-recapture brutes sont consolidees 
#-- au fil de l'eau au regard des analyses génétiques en continu (confirmation de génotype 
#-- ou information a posteriori par exemple au fur et a mesure de la consolidation 
#-- des empreintes de génotypages grâce aux recaptures) 
#-- les estimations communiquees a la DREAL AURA en 2021 sont les suivantes :

dat <- matrix(c(1995,   12.00,  17.1,
                1996,   15.00,  35.4,
                1997,   20.00,  47.7,
                1998,   24.00,  24.7,
                1999,   29.00,  58.4,
                2000,   26.00,  45.7,
                2001,   29.00,  76.1,
                2002,   36.00,  103.1,
                2003,   39.00,  97.2,
                2004,   43.50,  122.5,
                2005,   61.50,  125.2,
                2006,   52.50,  101.1,
                2007,   62.50,  118.7,
                2008,   76.50,  135.9,
                2009,   68.00,  134.3,
                2010,   78.00,  164.7,
                2011,   91.50,  201.2,
                2012,   89.50,  167.4,
                2013,   139.00, 337.1,
                2014,   129.00, 275.7,
                2015,   132.50, 341.4,
                2016,   185.50, 529.9,
                2017,   224.00, 540.7,
                2018,   273.00, 563.6), byrow = T, ncol = 3)
colnames(dat) <- c("an", "EMR", "CMR")
dat <- as.data.frame(dat)

# regression lineaire - estimation - variance - valeurs predites
mod <- lm(CMR ~ EMR, data = dat)
mod %>% tidy()
mod %>% glance()
#mod %>% augment()

# recherche point de rupture
b1 <- function(x, bp) ifelse(x < bp, bp - x, 0)
b2 <- function(x, bp) ifelse(x < bp, 0, x - bp)

# fonction qui calcule la deviance de la regression par morceaux
foo <- function(bp){
  mod <- lm(CMR ~ b1(EMR, bp) + b2(EMR, bp), data = dat)
  deviance(mod)
}

# cherche point de rupture
search.range <- c(min(dat$EMR) + 10, max(dat$EMR) - 10)
foo.opt <- optimize(foo, interval = search.range)
bp <- foo.opt$minimum
bp

# intervalle de confiance pour la valeur au point de rupture
foo.root <- function(bp, tgt){ foo(bp) - tgt }
tgt <- foo.opt$objective + qchisq(0.95,1)
lb95 <- uniroot(foo.root, lower = search.range[1], upper = bp, tgt = tgt)
ub95 <- uniroot(foo.root, lower = bp, upper = search.range[2], tgt = tgt)
lb95$root
ub95$root

# ajuste regression par morceaux avec cette valeur
mod2 <- lm(CMR ~ b1(EMR, bp) + b2(EMR, bp), data = dat)

# comparaison ajustement des 2 modeles mod = lineaire, et mod2 = reg par morceaux
# la regression par morceaux recoit plus de support selon AIC
AIC(mod, mod2)

# graphiquement
dat %>% 
  add_column(segreg = predict(mod2),
             linreg = predict(mod)) %>%
  ggplot() + 
  geom_line(aes(x = EMR, y = segreg, color = "green"), size = 1.3) + 
  geom_line(aes(x = EMR, y = linreg, color = "blue"), size = 1.3) + 
  scale_color_discrete(name = "Méthodes de régression", labels = c("linéaire", "segmentée")) + 
  geom_point(aes(x = EMR, y = CMR), size = 5, color = "white") +
  geom_point(aes(x = EMR, y = CMR), size = 3, color = "black") +
  labs(x = 'EMR', y = 'CMR')

#-- prediction du nombre de loups 
#-- hivers 2019-2020 et 2020-2021
#-- sur base calibration capture-recapture EMR
new <- data.frame(an = c(2019, 2020), 
                  EMR = c(301,403))
pred_segreg <- predict(mod2, new, interval = "prediction")
pred_segreg

#----------- 7. figure des effectifs estimes par cmr, et prediction via la calibration emr - cmr

dat <- read_csv2("dat/nbloupsestime.csv") %>%
  mutate(an = as.integer(an),
         hiver = paste0(1995:2020, "-", 1996:2021))
dat_cmr <- dat %>%
  filter(an < 2019)
dat_calib <- dat %>%
  filter(an > 2018)

tendance_nb_loups <- dat_cmr %>%
  ggplot(aes(x = an, y = cmr)) + 
  geom_line() + 
  geom_point(size = 5, color = "white") +
  geom_point(size = 3, aes(color = "black")) +
  geom_ribbon(aes(ymin = ic_bas, ymax = ic_haut), alpha = 0.3) + 
  geom_line(data = dat_calib, aes(x = an, y = cmr), color = "orange") + 
  geom_point(data = dat_calib, size = 5, color = "white") +
  geom_point(data = dat_calib, size = 3, aes(color = "orange")) +
  geom_ribbon(data = dat_calib, aes(ymin = ic_bas, ymax = ic_haut), alpha = 0.3, fill = "orange") +  
  labs(x = "hivers", y = "nombre estimé de loups") +
  scale_colour_manual("", values = c("black", "orange"), 
                      labels = c("estimation par capture-recapture", 
                                 "prédiction via calibration")) + 
  theme(legend.position = c(0.3, 0.8),
        legend.background = element_rect(fill = alpha('white', 0.4)))
tendance_nb_loups

ggsave(plot = tendance_nb_loups, 
       filename = "fig/nbloupestim.png", 
       dpi = 600, 
       width = 16, 
       height = 10)

#----------- 8.  Info sur la session R

sessionInfo()

# R version 4.1.0 (2021-05-18)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] broom_0.7.9      R2ucare_1.0.0    janitor_2.1.0    forcats_0.5.1    stringr_1.4.0   
# [6] dplyr_1.0.7      purrr_0.3.4.9000 readr_2.0.0      tidyr_1.1.3      tibble_3.1.4    
# [11] ggplot2_3.3.5    tidyverse_1.3.1


