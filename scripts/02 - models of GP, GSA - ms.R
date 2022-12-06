# 05a - models of GP ~ intercept or slope of GSA allometry

	

	## without phylo ##
	
	# GP ~ intercept #

	RawGSA8_phylo$id = as.numeric(factor(RawGSA8_phylo$Binomial))

	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo)

	dat <- list(
							N=nrow(RawGSA8_phylo),
							J=len(RawGSA8_phylo$Binomial),
							sp=RawGSA8_phylo$id,
							LogGSAcm2=RawGSA8_phylo$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo$LogCenterMassG,
							GP_int = unique(RawGSA8_phylo$mean_gperf_simple))


	MLM_GP_int_std <- stan(file = here("./stan models/correct_MLM_GP_INT_STD.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "aGP_int",
	              		  				 "bGP_int",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		           "sigma",
	              		  				 "sigma_GP_int",
	              		  				 "beta_ints_std",
	              		  				 "log_lik"))

	#save(MLM_GP_int_std, file = here("./output/MLM_GP_int_std.rds"))
	#load(file = here("./output/MLM_GP_int_std.rds"))

	post_MLM_GP_int_std <- as.data.frame(MLM_GP_int_std)
	names_post_MLM_GP_int_std <- names(post_MLM_GP_int_std)

	# summary of output: mean effect sizes, sds of effect sizes, & 95% Bayesian Credible Intervals
	post_MLM_GP_int__std_mean <- stack(apply(post_MLM_GP_int_std, 2, mean)) %>% rename(mean = values)
	post_MLM_GP_int_std_LCIs <- stack(apply(post_MLM_GP_int_std, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_MLM_GP_int_std_HCIs <- stack(apply(post_MLM_GP_int_std, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_MLM_GP_int_std_sum <- merge(post_MLM_GP_int__std_mean, post_MLM_GP_int_std_LCIs, by = "ind") %>%
												     merge(post_MLM_GP_int_std_HCIs) %>% round_df(2)
	post_MLM_GP_int_std_sd <- stack(apply(post_MLM_GP_int_std, 2, sd)) %>% rename(mean = values)

	# proportion of the distribution that is > 3
	GP_int_slope_prob0 <- (length(which(post_MLM_GP_int_std$bGP_int >0)))/length(post_MLM_GP_int_std$bGP_int)


	# GP ~ slope #

	RawGSA8_phylo$id = as.numeric(factor(RawGSA8_phylo$Binomial))
	
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo)
	
	dat <- list(
							N=nrow(RawGSA8_phylo),
							J=len(RawGSA8_phylo$Binomial),
							sp=RawGSA8_phylo$id,
							LogGSAcm2=RawGSA8_phylo$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo$LogCenterMassG,
							GP_slope = unique(RawGSA8_phylo$mean_gperf_simple))


	MLM_GP_slope_std <- stan(file = here("./stan models/correct_MLM_GP_slope_STD.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "aGP_slope",
	              		  				 "bGP_slope",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		  				 "beta_slopes_std",
	              		           "sigma",
	              		  				 "sigma_GP_slope",
	              		  				 "log_lik"))

	save(MLM_GP_slope_std, file = here("./output/MLM_GP_slope_std.rds"))
	load(file = here("./output/MLM_GP_slope_std.rds"))

	post_MLM_GP_slope_std <- as.data.frame(MLM_GP_slope_std)
	names_post_MLM_GP_slope_std <- names(post_MLM_GP_slope_std)

	# summary of output: mean effect sizes, sds of effect sizes, & 95% Bayesian Credible Intervals
	post_MLM_GP_slope_std_mean <- stack(apply(post_MLM_GP_slope_std, 2, mean)) %>% rename(mean = values)
	post_MLM_GP_slope_std_LCIs <- stack(apply(post_MLM_GP_slope_std, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_MLM_GP_slope_std_HCIs <- stack(apply(post_MLM_GP_slope_std, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_MLM_GP_slope_std_sum <- merge(post_MLM_GP_slope_std_mean, post_MLM_GP_slope_std_LCIs, by = "ind") %>%
												       merge(post_MLM_GP_slope_std_HCIs) %>% round_df(2)
	post_MLM_GP_slope_std_sd <- stack(apply(post_MLM_GP_slope_std, 2, sd)) %>% rename(mean = values)

	# proportion of the distribution that is > 3
	GP_slope_slope_prob0 <- (length(which(post_MLM_GP_slope_std$bGP_slope >0)))/length(post_MLM_GP_slope_std$bGP_slope)

	## with phylo ####

	# load phylo and wrangle into correct format

	# identity matrix for stan
	d_mat <- diag(1, 32, 32)

	# prepare tree for stan
	inv.tree <- inverseA(tree_pruned, nodes = "TIPS", scale = TRUE) #works
	A <- solve(inv.tree$Ainv) # is this reversing the inverse so the matrix is in the form of a covariance matrix and not its inverse (required by brms)
	rownames(A) <- rownames(inv.tree$Ainv) # assigning rownames of species ( = to the tips of the tree) to the covariance matrix
	
	vcov_mat <- as.matrix(A)

	# ontogenetic intercept #

	RawGSA8_phylo$id = as.numeric(factor(RawGSA8_phylo$Binomial))

	dat <- list(
							N=nrow(RawGSA8_phylo),
							J=len(RawGSA8_phylo$Binomial),
							sp=RawGSA8_phylo$id,
							LogGSAcm2=RawGSA8_phylo$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo$LogCenterMassG,
							GP_int = unique(RawGSA8_phylo$mean_gperf_simple),
							d_mat = d_mat,
							vcov_mat = vcov_mat)
	
	
	MLM_GP_int_phylo_std <- stan(
	                file = ("./stan models/correct_MLM_stan_phylo_GP_int_STD.stan"),
	                data = dat,
	                iter = 5000,
	                warmup = 1000,
	                chains = 4,
	                control = list("adapt_delta" = 0.99,
	                							 "max_treedepth" = 18),
	                pars = c(
	                				 "global_slope",
	              		  		 "global_int",
	              		  		 "aGP_int",
	              		  		 "bGP_int",
	              		       "beta_int",
	              		  		 "beta_slope",
	                				 "log_lik",
	              		       "sigma",
	              		  		 "sigma_GP_int",
	                				 "sigma_total_GP_int",
	                				 "lambda_GP_int"))

	#save(MLM_GP_int_phylo_std, file = here("./output/MLM_GP_int_phylo_std.rds"))
	#load(file = here("./output/MLM_GP_int_phylo_std.rds"))

	post_MLM_GP_int_std_phylo <- as.data.frame(MLM_GP_int_phylo_std)
	post_MLM_GP_int_phylo_std_mean <- stack(apply(post_MLM_GP_int_std_phylo, 2, mean))
	post_MLM_GP_int_phylo_std_LCIs <- stack(apply(post_MLM_GP_int_std_phylo, 2, quantile, 
																						prob = (0.025))) %>% rename(LCI = values)
	post_MLM_GP_int_phylo_std_HCIs <- stack(apply(post_MLM_GP_int_std_phylo, 2, quantile, 
																						prob = (0.975))) %>% rename(HCI = values)
	
	post_MLM_GP_int_phylo_std_SD <- stack(apply(post_MLM_GP_int_std_phylo, 2, sd))
	
	post_MLM_GP_int_phylo_std_sum <- merge(post_MLM_GP_int_phylo_std_mean, post_MLM_GP_int_phylo_std_LCIs, by = "ind") %>%
												           merge(post_MLM_GP_int_phylo_std_HCIs) %>% round_df(2)
	
	head(post_MLM_GP_int_phylo_std_sum)
	
	fortable <- post_MLM_GP_int_phylo_std[grep("lambda", post_MLM_GP_int_phylo_std$ind), ]
	fortable <- post_MLM_GP_int_phylo_std[grep("bGP", post_MLM_GP_int_phylo_std$ind), ]

	GP_int_slope_phylo_prob0 <- (length(which(post_MLM_GP_int_std_phylo$bGP_int >0)))/length(post_MLM_GP_int_std_phylo$bGP_int)


	# ontogenetic slope #

	RawGSA8_phylo$id = as.numeric(factor(RawGSA8_phylo$Binomial))
	
	dat <- list(
							N=nrow(RawGSA8_phylo),
							J=len(RawGSA8_phylo$Binomial),
							sp=RawGSA8_phylo$id,
							LogGSAcm2=RawGSA8_phylo$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo$LogCenterMassG,
							GP_slope = unique(RawGSA8_phylo$mean_gperf_simple),
							d_mat = d_mat,
							vcov_mat = vcov_mat)
	
	MLM_GP_slope_phylo_std <- stan(
	                file = here("./stna models/correct_MLM_stan_phylo_GP_slope_std.stan"),
	                data = dat,
	                iter = 5000,
	                warmup = 1000,
	                chains = 4,
	                control = list("adapt_delta" = 0.99,
	                							 "max_treedepth" = 18),
	                pars = c(
	                				 "global_slope",
	              		  		 "global_int",
	              		  		 "aGP_slope",
	              		  		 "bGP_slope",
	              		       "beta_int",
	              		  		 "beta_slope",
	                				 "log_lik",
	              		       "sigma",
	              		  		 "sigma_GP_slope",
	                				 "sigma_total_GP_slope",
	                				 "lambda_GP_slope"))
	
	#save(MLM_GP_slope_phylo_std, file = here("./output/MLM_GP_slope_phylo_std.rds"))
	#load(file = here("./output/MLM_GP_slope_phylo_std.rds"))
	
	post_MLM_GP_slope_phylo_std <- as.data.frame(MLM_GP_slope_phylo_std)
	
	post_MLM_GP_slope_phylo_std_mean <- stack(apply(post_MLM_GP_slope_phylo_std, 2, mean))
	post_MLM_GP_slope_phylo_std_LCIs <- stack(apply(post_MLM_GP_slope_phylo_std, 2, quantile, 
																						prob = (0.025)))  %>% rename(LCI = values)
	post_MLM_GP_slope_phylo_std_HCIs <- stack(apply(post_MLM_GP_slope_phylo_std, 2, quantile, 
																						prob = (0.975)))  %>% rename(HCI = values)
	
	post_MLM_GP_slope_phylo_std_SD <- stack(apply(post_MLM_GP_slope_phylo_std, 2, sd))
	
	
	post_MLM_GP_slope_phylo_std_sum <- merge(post_MLM_GP_slope_phylo_std_mean, post_MLM_GP_slope_phylo_std_LCIs, 
																					 by = "ind") %>%
												             merge(post_MLM_GP_slope_phylo_std_HCIs) %>% round_df(2)
	
	head(post_MLM_GP_slope_phylo_std_sum)

	GP_slope_slope_phylo_prob0 <- (length(which(post_MLM_GP_slope_phylo_std$bGP_slope >0)))/length(post_MLM_GP_slope_phylo_std$bGP_slope)

	#### without aquaculture species ####
	
	RawGSA8_phylo_no_ac <- RawGSA8_phylo %>%
		filter(Binomial %notin% aquaculture_sp_mean)
	
	#### intercept model ####

	RawGSA8_phylo_no_ac$id = as.numeric(factor(RawGSA8_phylo_no_ac$Binomial))
	
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo_no_ac)
	
	dat <- list(
							N=nrow(RawGSA8_phylo_no_ac),
							J=len(RawGSA8_phylo_no_ac$Binomial),
							sp=RawGSA8_phylo_no_ac$id,
							LogGSAcm2=RawGSA8_phylo_no_ac$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo_no_ac$LogCenterMassG,
							GP_int = unique(RawGSA8_phylo_no_ac$mean_gperf_simple))
	

	GP_INT_NO_AC <- stan(file = "./stan models/MLM_GP_INT_EX.stan",
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "aGP_int",
	              		  				 "bGP_int",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		  				 "beta_ints_std",
	              		           "sigma",
	              		  				 "sigma_GP_int",
	              		  				 "log_lik"))
	
	#save(GP_INT_NO_AC, file = ("./analyses/manuscript code/without aquaculture species/GP_INT_NO_AC.rds"))
	#load(file = ("./analyses/manuscript code/without aquaculture species/GP_INT_NO_AC.rds"))


	post_GP_INT_NO_AC <- as.data.frame(GP_INT_NO_AC)
	names_post_GP_INT_NO_AC <- names(post_GP_INT_NO_AC)
	
	post_GP_INT_NO_AC_mean <- stack(apply(post_GP_INT_NO_AC, 2, mean)) %>% rename(mean = values)
	
	post_GP_INT_NO_AC_LCIs <- stack(apply(post_GP_INT_NO_AC, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_INT_NO_AC_HCIs <- stack(apply(post_GP_INT_NO_AC, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_GP_INT_NO_AC_sum <- merge(post_GP_INT_NO_AC_mean, post_GP_INT_NO_AC_LCIs, by = "ind") %>%
												   merge(post_GP_INT_NO_AC_HCIs) 
	
	post_GP_INT_NO_AC_sd <- stack(apply(post_GP_INT_NO_AC, 2, sd)) %>% rename(mean = values)
	
	fortable <- post_GP_INT_NO_AC_sum[grep("GP_int", post_GP_INT_NO_AC_sum$ind), ] %>% round_df(2)


	#### slope model ####
	
	RawGSA8_phylo_no_ac$id = as.numeric(factor(RawGSA8_phylo_no_ac$Binomial))
	
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo_no_ac)
	
	dat <- list(
							N=nrow(RawGSA8_phylo_no_ac),
							J=len(RawGSA8_phylo_no_ac$Binomial),
							sp=RawGSA8_phylo_no_ac$id,
							LogGSAcm2=RawGSA8_phylo_no_ac$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo_no_ac$LogCenterMassG,
							GP_slope = unique(RawGSA8_phylo_no_ac$mean_gperf_simple))
	
	
	MLM_GP_slope_std_ex <- stan(file = "./stan models/MLM_GP_SLOPE_EX.stan",
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "aGP_slope",
	              		  				 "bGP_slope",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		  				 "beta_slopes_std",
	              		           "sigma",
	              		  				 "sigma_GP_slope",
	              		  				 "log_lik"))
	
	#save(GP_SLOPE_NO_AC, file = ("./analyses/manuscript code/without aquaculture species/GP_SLOPE_NO_AC.rds"))
	#load(file = ("./stan models/MLM_GP_SLOPE_EX.rds"))
	
	
	post_GP_SLOPE_NO_AC <- as.data.frame(MLM_GP_slope_std_ex)
	names_post_GP_SLOPE_NO_AC <- names(post_GP_SLOPE_NO_AC)
	
	post_GP_SLOPE_NO_AC_mean <- stack(apply(post_GP_SLOPE_NO_AC, 2, mean)) %>% rename(mean = values)
	
	post_GP_SLOPE_NO_AC_LCIs <- stack(apply(post_GP_SLOPE_NO_AC, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_SLOPE_NO_AC_HCIs <- stack(apply(post_GP_SLOPE_NO_AC, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_GP_SLOPE_NO_AC_sum <- merge(post_GP_SLOPE_NO_AC_mean, post_GP_SLOPE_NO_AC_LCIs, by = "ind") %>%
												     merge(post_GP_SLOPE_NO_AC_HCIs) 
	
	post_GP_SLOPE_NO_AC_sd <- stack(apply(post_GP_SLOPE_NO_AC, 2, sd)) %>% rename(mean = values)
	
	fortable <- post_GP_SLOPE_NO_AC_sum[grep("GP_slope", post_GP_SLOPE_NO_AC_sum$ind), ] %>% round_df(2)
	

	#### without air-breathers ####
	
	#### intercept model ####
	
	RawGSA8_phylo_no_abr$id = as.numeric(factor(RawGSA8_phylo_no_abr$Binomial))
	
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo_no_abr)
	
	dat <- list(
							N=nrow(RawGSA8_phylo_no_abr),
							J=len(RawGSA8_phylo_no_abr$Binomial),
							sp=RawGSA8_phylo_no_abr$id,
							LogGSAcm2=RawGSA8_phylo_no_abr$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo_no_abr$LogCenterMassG,
							GP_int = unique(RawGSA8_phylo_no_abr$mean_gperf_simple))
	
	
	GP_INT_NO_ABR <- stan(file = here("./stan models/correct_MLM_GP_INT_STD.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "aGP_int",
	              		  				 "bGP_int",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		  				 "beta_ints_std",
	              		           "sigma",
	              		  				 "sigma_GP_int",
	              		  				 "log_lik"))
	
	
	#save(GP_INT_NO_ABR, file = ("./analyses/manuscript code/without aquaculture species/GP_INT_NO_ABR.rds"))
	#load(file = ("./analyses/manuscript code/without aquaculture species/GP_INT_NO_ABR.rds"))

		post_GP_INT_NO_ABR <- as.data.frame(GP_INT_NO_ABR)
	
	post_GP_INT_NO_ABR_mean <- stack(apply(post_GP_INT_NO_ABR, 2, mean)) %>% rename(mean = values)
	
	post_GP_INT_NO_ABR_LCIs <- stack(apply(post_GP_INT_NO_ABR, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_INT_NO_ABR_HCIs <- stack(apply(post_GP_INT_NO_ABR, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_GP_INT_NO_ABR_sum <- merge(post_GP_INT_NO_ABR_mean, post_GP_INT_NO_ABR_LCIs, by = "ind") %>%
												       merge(post_GP_INT_NO_ABR_HCIs) 
	
	fortable <- post_GP_INT_NO_ABR_sum[grep("GP_int", post_GP_INT_NO_ABR_sum$ind), ] %>% round_df(2)
	
	
	#### slope ####
	
	RawGSA8_phylo_no_abr$id = as.numeric(factor(RawGSA8_phylo_no_abr$Binomial))
	
	data_matrix <- model.matrix(LogGSAcm2 ~ LogCenterMassG, data = RawGSA8_phylo_no_abr)
	
	dat <- list(
							N=nrow(RawGSA8_phylo_no_abr),
							J=len(RawGSA8_phylo_no_abr$Binomial),
							sp=RawGSA8_phylo_no_abr$id,
							LogGSAcm2=RawGSA8_phylo_no_abr$LogGSAcm2,
							LogCenterMassG = RawGSA8_phylo_no_abr$LogCenterMassG,
							GP_slope = unique(RawGSA8_phylo_no_abr$mean_gperf_simple))
	
	
	GP_SLOPE_NO_ABR <- stan(file = here("./stan models/correct_MLM_GP_slope_STD.stan"),
	              		  data = dat,
	              		  iter = 5000,
	              		  warmup = 1000,
	              		  chains = 4,
	              		  control = list("adapt_delta" = 0.9999,
	              		  							 "max_treedepth" = 18),
	              		  pars = c("global_slope",
	              		  				 "global_int",
	              		  				 "aGP_slope",
	              		  				 "bGP_slope",
	              		           "beta_int",
	              		  				 "beta_slope",
	              		  				 "beta_slopes_std",
	              		           "sigma",
	              		  				 "sigma_GP_slope",
	              		  				 "log_lik"))
	
	#save(GP_SLOPE_NO_ABR, file = here("./output/GP_SLOPE_NO_ABR.rds"))
	#load(file = here("./output/GP_SLOPE_NO_ABR.rds"))
	
	post_GP_SLOPE_NO_ABR <- as.data.frame(GP_SLOPE_NO_ABR)
	
	post_GP_SLOPE_NO_ABR_mean <- stack(apply(post_GP_SLOPE_NO_ABR, 2, mean)) %>% rename(mean = values)
	
	post_GP_SLOPE_NO_ABR_LCIs <- stack(apply(post_GP_SLOPE_NO_ABR, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_SLOPE_NO_ABR_HCIs <- stack(apply(post_GP_SLOPE_NO_ABR, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	
	post_GP_SLOPE_NO_ABR_sum <- merge(post_GP_SLOPE_NO_ABR_mean, post_GP_SLOPE_NO_ABR_LCIs, by = "ind") %>%
												       merge(post_GP_SLOPE_NO_ABR_HCIs) 
	
	fortable <- post_GP_SLOPE_NO_ABR_sum[grep("GP_slope", post_GP_SLOPE_NO_ABR_sum$ind), ] %>% round_df(2)
	