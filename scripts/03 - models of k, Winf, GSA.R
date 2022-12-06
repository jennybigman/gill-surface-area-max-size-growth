# 06 k, Winf, GSA 

	# transfomrations
	RawGSA8_phylo$Log_k <- log10(RawGSA8_phylo$mean_K)
	RawGSA8_phylo$Log_Winf <- log10(RawGSA8_phylo$mean_Winfinity) 
	RawGSA8_phylo$Winf_std <- ((RawGSA8_phylo$Log_Winf) - mean(RawGSA8_phylo$Log_Winf)) / sd(RawGSA8_phylo$Log_Winf)
	
	Log_k <- log10(RawGSA8_phylo$mean_K)
	Log_Winf <- log10(RawGSA8_phylo$mean_Winfinity) 
	Winf_std <- ((RawGSA8_phylo$Log_Winf) - mean(RawGSA8_phylo$Log_Winf)) / sd(RawGSA8_phylo$Log_Winf)

	#### k ~ Winf + GSA intercept ####
	RawGSA8_phylo$id = as.numeric(factor(RawGSA8_phylo$Binomial))
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo)
	
	dat <- list(
							N=nrow(RawGSA8_phylo),
							J=len(RawGSA8_phylo$Binomial),
							sp=RawGSA8_phylo$id,
							LogGSAcm2=RawGSA8_phylo$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo$LogCenterMassG,
							k = unique(RawGSA8_phylo$Log_k),
							Winf = unique(RawGSA8_phylo$Winf_std))


	k_Winf_INT <- stan(file = here("./stan models/K_WINF_INT.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "a_k_int",
	              		  				 "b_GSA_int",
	              		  				 "b_Winf_int",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		           "sigma",
	              		  				 "sigma_k_int",
	              		  				 "beta_ints_std",
	              		  				 "log_lik"))

	#save(k_Winf_INT, file = here("./output/k_Winf_INT.rds"))
	#load(file = here("./output/k_Winf_INT.rds"))

	post_k_Winf_INT <- as.data.frame(k_Winf_INT)
	post_k_Winf_INT_mean <- stack(apply(post_k_Winf_INT, 2, mean)) %>% rename(mean = values)
	post_k_Winf_INT_LCIs <- stack(apply(post_k_Winf_INT, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_k_Winf_INT_HCIs <- stack(apply(post_k_Winf_INT, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_k_Winf_INT_sum <- merge(post_k_Winf_INT_mean, post_k_Winf_INT_LCIs, by = "ind") %>%
												 merge(post_k_Winf_INT_HCIs) 
	post_k_Winf_INT_sum <- round_df(post_k_Winf_INT_sum, 2) %>% round_df(2)

	head(post_k_Winf_INT_sum)
	
	# proportion > 0
	k_Winf_INT_GSA_prob0 <- (length(which(post_k_Winf_INT$b_GSA_int >0)))/length(post_k_Winf_INT$b_GSA_int)*100 
	k_Winf_INT_GSA_prob0 <- round(k_Winf_INT_GSA_prob0, 1)
	k_Winf_INT_Winf_prob0 <- (length(which(post_k_Winf_INT$b_Winf_int >0)))/length(post_k_Winf_INT$b_Winf_int)*100
	k_Winf_INT_Winf_prob0 <- round(k_Winf_INT_Winf_prob0, 1)
	
	# calculate VIFs
	ints <- post_k_Winf_INT[72:103]
	sp_names <- sort(unique(RawGSA8_phylo$Binomial))
	names(ints) <- sp_names 
	ints_sum <- stack(apply(ints, 2, mean)) %>%
							rename(Binomial = ind,
										 int = values)
	
	kWinf_int_dat <- distinct(RawGSA8_phylo, Log_k, Log_Winf, Winf_std) 
	
	ints_sum <- merge(ints_sum, kWinf_int_dat, by = "Binomial")

	# fit lm
	
	lm_kWinf_int <- lm(Log_k ~ int + Winf_std, data = ints_sum)
	summary(lm_kWinf_int)
	confint(lm_kWinf_int)
	
	vif(lm_kWinf_int)
	
	dat_int <- ints_sum %>% dplyr::select(Winf_std, int)
	cor(dat_int)


	####  k ~ Winf + GSA slope ####

	RawGSA8_phylo$id = as.numeric(factor(RawGSA8_phylo$Binomial))
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo)
	
	dat <- list(
							N=nrow(RawGSA8_phylo),
							J=len(RawGSA8_phylo$Binomial),
							sp=RawGSA8_phylo$id,
							LogGSAcm2=RawGSA8_phylo$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo$LogCenterMassG,
							k = unique(RawGSA8_phylo$Log_k),
							Winf = unique(RawGSA8_phylo$Winf_std))


	k_Winf_slope <- stan(file = here("./stan models/K_WINF_SLOPE.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "a_k_slope",
	              		  				 "b_GSA_slope",
	              		  				 "b_Winf_slope",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		           "sigma",
	              		  				 "sigma_k_slope",
	              		  				 "beta_slopes_std",
	              		  				 "log_lik"))
	
	#save(k_Winf_slope, file = here("./output/k_Winf_slope.rds"))
	#load(file = here("./output/k_Winf_slope.rds"))
			 
	post_k_Winf_slope <- as.data.frame(k_Winf_slope)
	post_k_Winf_slope_mean <- stack(apply(post_k_Winf_slope, 2, mean)) %>% rename(mean = values)
	post_k_Winf_slope_LCIs <- stack(apply(post_k_Winf_slope, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_k_Winf_slope_HCIs <- stack(apply(post_k_Winf_slope, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_k_Winf_slope_sum <- merge(post_k_Winf_slope_mean, post_k_Winf_slope_LCIs, by = "ind") %>%
												   merge(post_k_Winf_slope_HCIs) 
	post_k_Winf_slope_sum <- round_df(post_k_Winf_slope_sum, 2)
	
	head(post_k_Winf_slope_sum)

	# proportion > 0
	k_Winf_slope_GSA_prob0 <- (length(which(post_k_Winf_slope$b_GSA_slope >0)))/length(post_k_Winf_slope$b_GSA_slope)
	k_Winf_slope_Winf_prob0 <- (length(which(post_k_Winf_slope$b_Winf_slope >0)))/length(post_k_Winf_slope$b_GSA_slope)

	# calculate VIFs
	slopes <- post_k_Winf_slope[72:103]
	sp_names <- sort(unique(RawGSA8_phylo$Binomial))
	names(slopes) <- sp_names 
	slopes_sum <- stack(apply(slopes, 2, mean)) %>%
							rename(Binomial = ind,
										 slopes = values)
	kWinf_slope_dat <- distinct(RawGSA8_phylo, Log_k, Log_Winf, Winf_std) 
	slopes_sum <- merge(slopes_sum, kWinf_slope_dat, by = "Binomial")

	# fit lm
	lm_kWinf_slope <- lm(Log_k ~ slopes + Winf_std, data = slopes_sum)
	summary(lm_kWinf_slope)
	confint(lm_kWinf_slope)
	vif(lm_kWinf_slope)
	dat_slope <- slopes_sum %>% dplyr::select(Winf_std, slopes)
	cor(dat_slope)


	#### run models with phylogeny ####

	# identity matrix for stan
	d_mat <- diag(1, 32, 32)
	
	# prepare tree for stan
	inv.tree <- inverseA(tree_pruned, nodes = "TIPS", scale = TRUE) #works
	A <- solve(inv.tree$Ainv) # is this reversing the inverse so the matrix is in the form of a covariance matrix and not its inverse (required by brms)
	rownames(A) <- rownames(inv.tree$Ainv) # assigning rownames of species ( = to the tips of the tree) to the covariance matrix
	
	vcov_mat <- as.matrix(A)

	#### k ~ Winf + GSA intercept ####
	RawGSA8_phylo$id = as.numeric(factor(RawGSA8_phylo$Binomial))
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo)
		
	dat <- list(
							N=nrow(RawGSA8_phylo),
							J=len(RawGSA8_phylo$Binomial),
							sp=RawGSA8_phylo$id,
							LogGSAcm2=RawGSA8_phylo$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo$LogCenterMassG,
							k = unique(RawGSA8_phylo$Log_k),
							Winf = unique(RawGSA8_phylo$Winf_std),
							d_mat = d_mat,
							vcov_mat = vcov_mat)

	k_Winf_INT_PHYLO <- stan(file = here("./stan models/K_WINF_INT_PHYLO.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "a_k_int",
	              		  				 "b_GSA_int",
	              		  				 "b_Winf_int",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		           "sigma",
	              		  				 "sigma_k_int",
	              		  				 "beta_ints_std",
	              		  				 "log_lik",
	              		  				 "lambda_k"))

	#save(k_Winf_INT_PHYLO, file = here("./output/k_Winf_INT_PHYLO.rds"))
	#load(file = here("./output/k_Winf_INT_PHYLO.rds"))
	
	post_k_Winf_INT_PHYLO <- as.data.frame(k_Winf_INT_PHYLO)
	post_k_Winf_INT_PHYLO_mean <- stack(apply(post_k_Winf_INT_PHYLO, 2, mean)) %>% 
																rename(mean = values)
	post_k_Winf_INT_PHYLO_LCIs <- stack(apply(post_k_Winf_INT_PHYLO, 2, quantile, 
																						prob = (0.025))) %>% rename(LCI = values)
	post_k_Winf_INT_PHYLO_HCIs <- stack(apply(post_k_Winf_INT_PHYLO, 2, quantile, 
																						prob = (0.975))) %>% rename(HCI = values)
	post_k_Winf_INT_PHYLO_sum <- merge(post_k_Winf_INT_PHYLO_mean, post_k_Winf_INT_PHYLO_LCIs, 
																		 by = "ind") %>%
												       merge(post_k_Winf_INT_PHYLO_HCIs) 
	
	post_k_Winf_INT_PHYLO_sum <- round_df(post_k_Winf_INT_PHYLO_sum, 2)

	head(post_k_Winf_INT_PHYLO_sum)
	
	# proportion > 0
	k_Winf_int_phylo_GSAint_prob0 <- (length(which(post_k_Winf_INT_PHYLO$b_GSA_int >0)))/length(post_k_Winf_INT_PHYLO$b_GSA_int)
	k_Winf_int_phylo_Winfint_prob0 <- (length(which(post_k_Winf_INT_PHYLO$b_Winf_int >0)))/length(post_k_Winf_INT_PHYLO$b_Winf_int)

	####  k ~ Winf + GSA slope ####
	RawGSA8_phylo$id = as.numeric(factor(RawGSA8_phylo$Binomial))
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo)
		
	dat <- list(
							N=nrow(RawGSA8_phylo),
							J=len(RawGSA8_phylo$Binomial),
							sp=RawGSA8_phylo$id,
							LogGSAcm2=RawGSA8_phylo$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo$LogCenterMassG,
							k = unique(RawGSA8_phylo$Log_k),
							Winf = unique(RawGSA8_phylo$Winf_std),
							d_mat = d_mat,
							vcov_mat = vcov_mat)
	
	k_Winf_slope_PHYLO <- stan(file = here("./stan models/K_WINF_SLOPE_PHYLO.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "a_k_slope",
	              		  				 "b_GSA_slope",
	              		  				 "b_Winf_slope",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		           "sigma",
	              		  				 "sigma_k_slope",
	              		  				 "beta_slopes_std",
	              		  				 "log_lik",
	              		  				 "lambda_k"))

	#save(k_Winf_slope_PHYLO, file = here("./output/k_Winf_slope_PHYLO.rds"))
	#load(file = here("./output/k_Winf_slope_PHYLO.rds"))
			 
	post_k_Winf_slope_PHYLO <- as.data.frame(k_Winf_slope_PHYLO)
	post_k_Winf_slope_PHYLO_mean <- stack(apply(post_k_Winf_slope_PHYLO, 2, mean)) %>% 
																	rename(mean = values)
	post_k_Winf_slope_PHYLO_LCIs <- stack(apply(post_k_Winf_slope_PHYLO, 2, quantile, 
																							prob = (0.025))) %>% rename(LCI = values)
	post_k_Winf_slope_PHYLO_HCIs <- stack(apply(post_k_Winf_slope_PHYLO, 2, quantile, 
																							prob = (0.975))) %>% rename(HCI = values)
	post_k_Winf_slope_PHYLO_sum <- merge(post_k_Winf_slope_PHYLO_mean, post_k_Winf_slope_PHYLO_LCIs, 
																			 by = "ind") %>%
												         merge(post_k_Winf_slope_PHYLO_HCIs) 
	post_k_Winf_slope_PHYLO_sum <- round_df(post_k_Winf_slope_PHYLO_sum, 2)
	
	# proportion > 0
	k_Winf_slope_phylo_GSA_slope_prob0 <- (length(which(post_k_Winf_slope_PHYLO$b_GSA_slope >0)))/length(post_k_Winf_slope_PHYLO$b_GSA_slope)
	k_Winf_slope_phylo_Win_slope_prob0 <- (length(which(post_k_Winf_slope_PHYLO$b_Winf_slope >0)))/length(post_k_Winf_slope_PHYLO$b_Winf_slope)


#### no aquaculture species ####

	RawGSA8_phylo_no_ac <- filter(RawGSA8_phylo, Binomial %notin% aquaculture_sp_mean)
	
	#### k ~ Winf + GSA intercept ####
	
	RawGSA8_phylo_no_ac$id = as.numeric(factor(RawGSA8_phylo_no_ac$Binomial))
	
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo_no_ac)
		
	dat <- list(
							N=nrow(RawGSA8_phylo_no_ac),
							J=len(RawGSA8_phylo_no_ac$Binomial),
							sp=RawGSA8_phylo_no_ac$id,
							LogGSAcm2=RawGSA8_phylo_no_ac$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo_no_ac$LogCenterMassG,
							k = unique(RawGSA8_phylo_no_ac$Log_k),
							Winf = unique(RawGSA8_phylo_no_ac$Winf_std))
	
	
	k_Winf_INT_no_ac <- stan(file = here("./stan models/K_WINF_INT.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "a_k_int",
	              		  				 "b_GSA_int",
	              		  				 "b_Winf_int",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		           "sigma",
	              		  				 "sigma_k_int",
	              		  				 "beta_ints_std",
	              		  				 "log_lik"))
	
	save(k_Winf_INT_no_ac, file = ("./output/k_Winf_INT_no_ac.rds"))
	
	load(file = here("./output/k_Winf_INT_no_ac.rds"))
	
	post_k_Winf_INT_no_ac <- as.data.frame(k_Winf_INT_no_ac)
	
	post_k_Winf_INT_no_ac_mean <- stack(apply(post_k_Winf_INT_no_ac, 2, mean)) %>% rename(mean = values)
	
	post_k_Winf_INT_no_ac_LCIs <- stack(apply(post_k_Winf_INT_no_ac, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_k_Winf_INT_no_ac_HCIs <- stack(apply(post_k_Winf_INT_no_ac, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_k_Winf_INT_no_ac_sum <- merge(post_k_Winf_INT_no_ac_mean, post_k_Winf_INT_no_ac_LCIs, 
																		 by = "ind") %>%
												 merge(post_k_Winf_INT_no_ac_HCIs) 
	
	post_k_Winf_INT_no_ac_sum <- round_df(post_k_Winf_INT_no_ac_sum, 2) %>% round_df(2)
	
	head(post_k_Winf_INT_no_ac_sum)
	
	####  k ~ Winf + GSA slope ####
	
	RawGSA8_phylo_no_ac$id = as.numeric(factor(RawGSA8_phylo_no_ac$Binomial))
	
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo_no_ac)
		
	dat <- list(
							N=nrow(RawGSA8_phylo_no_ac),
							J=len(RawGSA8_phylo_no_ac$Binomial),
							sp=RawGSA8_phylo_no_ac$id,
							LogGSAcm2=RawGSA8_phylo_no_ac$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo_no_ac$LogCenterMassG,
							k = unique(RawGSA8_phylo_no_ac$Log_k),
							Winf = unique(RawGSA8_phylo_no_ac$Winf_std))
	
	
	k_Winf_slope_no_ac <- stan(file = ("./stan models/K_WINF_SLOPE.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "a_k_slope",
	              		  				 "b_GSA_slope",
	              		  				 "b_Winf_slope",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		           "sigma",
	              		  				 "sigma_k_slope",
	              		  				 "beta_slopes_std",
	              		  				 "log_lik"))
	
	save(k_Winf_slope_no_ac, file = ("./output/k_Winf_slope_no_ac.rds"))
	
	load(file = here("./output/k_Winf_slope_no_ac.rds"))
			 
	post_k_Winf_slope_no_ac <- as.data.frame(k_Winf_slope_no_ac)
	
	post_k_Winf_slope_no_ac_mean <- stack(apply(post_k_Winf_slope_no_ac, 2, mean)) %>% rename(mean = values)
	
	post_k_Winf_slope_no_ac_LCIs <- stack(apply(post_k_Winf_slope_no_ac, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_k_Winf_slope_no_ac_HCIs <- stack(apply(post_k_Winf_slope_no_ac, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_k_Winf_slope_no_ac_sum <- merge(post_k_Winf_slope_no_ac_mean, 
																			 post_k_Winf_slope_no_ac_LCIs, by = "ind") %>%
												         merge(post_k_Winf_slope_no_ac_HCIs) 
	
	post_k_Winf_slope_no_ac_sum <- round_df(post_k_Winf_slope_no_ac_sum, 2)
	
	head(post_k_Winf_slope_no_ac_sum)
	

	#### without airbreathing species ####
	
	RawGSA8_phylo_no_abr <- filter(RawGSA8_phylo, Binomial %notin% airbreathing_mean)
	
	#### k ~ Winf + GSA intercept ####
	
	RawGSA8_phylo_no_abr$id = as.numeric(factor(RawGSA8_phylo_no_abr$Binomial))
	
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo_no_abr)
		
	dat <- list(
							N=nrow(RawGSA8_phylo_no_abr),
							J=len(RawGSA8_phylo_no_abr$Binomial),
							sp=RawGSA8_phylo_no_abr$id,
							LogGSAcm2=RawGSA8_phylo_no_abr$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo_no_abr$LogCenterMassG,
							k = unique(RawGSA8_phylo_no_abr$Log_k),
							Winf = unique(RawGSA8_phylo_no_abr$Winf_std))
	
	
	k_Winf_INT_no_abr <- stan(file = ("./stan models/K_WINF_INT.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "a_k_int",
	              		  				 "b_GSA_int",
	              		  				 "b_Winf_int",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		           "sigma",
	              		  				 "sigma_k_int",
	              		  				 "beta_ints_std",
	              		  				 "log_lik"))
	
	save(k_Winf_INT_no_abr, file = ("./output/k_Winf_INT_no_abr.rds"))
	
	load(file = here("./output/k_Winf_INT_no_abr.rds"))
	
	post_k_Winf_INT_no_abr <- as.data.frame(k_Winf_INT_no_abr)
	
	post_k_Winf_INT_no_abr_mean <- stack(apply(post_k_Winf_INT_no_abr, 2, mean)) %>% rename(mean = values)
	
	post_k_Winf_INT_no_abr_LCIs <- stack(apply(post_k_Winf_INT_no_abr, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_k_Winf_INT_no_abr_HCIs <- stack(apply(post_k_Winf_INT_no_abr, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_k_Winf_INT_no_abr_sum <- merge(post_k_Winf_INT_no_abr_mean, post_k_Winf_INT_no_abr_LCIs, 
																		 by = "ind") %>%
												 merge(post_k_Winf_INT_no_abr_HCIs) 
	
	post_k_Winf_INT_no_abr_sum <- round_df(post_k_Winf_INT_no_abr_sum, 2) %>% round_df(2)
	
	head(post_k_Winf_INT_no_abr_sum)
	
	####  k ~ Winf + GSA slope ####
	
	RawGSA8_phylo_no_abr$id = as.numeric(factor(RawGSA8_phylo_no_abr$Binomial))
	
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo_no_abr)
		
	dat <- list(
							N=nrow(RawGSA8_phylo_no_abr),
							J=len(RawGSA8_phylo_no_abr$Binomial),
							sp=RawGSA8_phylo_no_abr$id,
							LogGSAcm2=RawGSA8_phylo_no_abr$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo_no_abr$LogCenterMassG,
							k = unique(RawGSA8_phylo_no_abr$Log_k),
							Winf = unique(RawGSA8_phylo_no_abr$Winf_std))
	
	
	k_Winf_slope_no_abr <- stan(file = ("./stan models/K_WINF_SLOPE.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "a_k_slope",
	              		  				 "b_GSA_slope",
	              		  				 "b_Winf_slope",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		           "sigma",
	              		  				 "sigma_k_slope",
	              		  				 "beta_slopes_std",
	              		  				 "log_lik"))
	
	save(k_Winf_slope_no_abr, file = ("./output/k_Winf_slope_no_abr.rds"))
	
	load(file = here("./output/k_Winf_slope_no_abr.rds"))
			 
	post_k_Winf_slope_no_abr <- as.data.frame(k_Winf_slope_no_abr)
	
	post_k_Winf_slope_no_abr_mean <- stack(apply(post_k_Winf_slope_no_abr, 2, mean)) %>% rename(mean = values)
	
	post_k_Winf_slope_no_abr_LCIs <- stack(apply(post_k_Winf_slope_no_abr, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_k_Winf_slope_no_abr_HCIs <- stack(apply(post_k_Winf_slope_no_abr, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_k_Winf_slope_no_abr_sum <- merge(post_k_Winf_slope_no_abr_mean, 
																			 post_k_Winf_slope_no_abr_LCIs, by = "ind") %>%
												         merge(post_k_Winf_slope_no_abr_HCIs) 
	
	post_k_Winf_slope_no_abr_sum <- round_df(post_k_Winf_slope_no_abr_sum, 2)

	head(post_k_Winf_slope_no_abr_sum)


