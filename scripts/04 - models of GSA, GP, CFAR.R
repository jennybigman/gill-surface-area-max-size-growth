# 07 - models of GSA, growth, activity


	# read in activity data and estimate CFAR

	CFAR_dat <- read.csv(here("./data/Mahan_CFAR_datasheet.csv"),
	                           header = TRUE, stringsAsFactors = FALSE,
	                           strip.white = TRUE)  %>%
	                           dplyr::select(-X.1, -X, - X.2, -X.3, -X.4, -X.5)

	CFAR_dat <- CFAR_dat %>%
							dplyr::select(Species, height, surface_area, image_source) %>%
						  drop_na() %>%
							mutate(Binomial = Species,
										 CFAR = (height^2)/surface_area)
	
	len(CFAR_dat$Species)

	# filter CFAR dat to only the species you have raw data for
	
	CFAR_dat_raw <- CFAR_dat %>%
							    filter(., Binomial %in% RawGSA8_phylo$Binomial)
	
	len(CFAR_dat_raw$Species)

	# which species do not have CFAR estimated from FAO?
	
	CFAR_dat_raw_FAO <- CFAR_dat_raw %>% filter(., image_source == "FAO")
	
	CFAR_dat_raw_no_FAO <- CFAR_dat_raw %>% filter(., image_source != "FAO")
	
	CFAR_dat_raw_no_FAO <- CFAR_dat_raw_no_FAO %>% filter(., Binomial %notin% CFAR_dat_raw_FAO$Binomial)
	
	# if a species has an FAO image, use this preferentially
	
	CFAR_dat_raw_sum <- CFAR_dat_raw %>%
											group_by(Binomial) %>%
										  summarise(mean_CFAR = if(any('image_source' == "FAO")) CFAR
										  											else mean(CFAR))
	
	len(CFAR_dat_raw_sum$Binomial)

	# standardize data
	
	CFAR_dat_raw_sum <- CFAR_dat_raw_sum %>%
											mutate(meanCFARstd = (mean_CFAR - mean(mean_CFAR)/sd(mean_CFAR)))
	
	# include only those species that have CFAR data
	
	RawGSA8_CFAR <- RawGSA8_phylo %>%
									filter(., Binomial %in% CFAR_dat_raw_sum$Binomial)
	
	len(RawGSA8_CFAR$Binomial)
	
	setdiff((unique(RawGSA8_phylo$Binomial)), (unique(RawGSA8_CFAR$Binomial)))

	# one eel, one ray

	## make sure data line up
	setdiff(unique(CFAR_dat_raw_sum$Binomial), unique(RawGSA8_CFAR$Binomial))
	RawGSA8_CFAR <- merge(RawGSA8_CFAR, CFAR_dat_raw_sum, by = "Binomial")


	## run stan models

	#### full raw dataset, n = 30 ####

	#### GP ~ intercept + CFAR ####

	RawGSA8_CFAR$id = as.numeric(factor(RawGSA8_CFAR$Binomial))
	
	dat <- list(
							N=nrow(RawGSA8_CFAR),
							J=len(RawGSA8_CFAR$Binomial),
							sp=RawGSA8_CFAR$id,
							LogGSAcm2=RawGSA8_CFAR$LogGSAcm2,
							LogCenterMassG = RawGSA8_CFAR$LogCenterMassG,
							CFAR = unique(RawGSA8_CFAR$meanCFARstd),
							GP = unique(RawGSA8_CFAR$mean_gperf_simple))
	
	GP_CFAR_int_std <- stan(
	                file = here("./stan models/MLM_CFAR_INT.stan"),
	                data = dat,
	                iter = 5000,
	                warmup = 2000,
	                chains = 4,
	                control = list("adapt_delta" = 0.99,
	                							 "max_treedepth" = 18),
	                pars = c( "global_slope",
	              		  			"global_int",
	              		  			"aGP",
	              		  			"bGP_int",
	                					"bGP_CFAR",
	              		        "beta_int",
	              		  			"beta_slope",
	              		  			"beta_ints_std",
	              		        "sigma",
	              		  			"sigma_GP",
	              		  			"log_lik"))

	#save(GP_CFAR_int_std, file = here("./output/GP_CFAR_int_std.rds"))
	#load(file = here("./output/GP_CFAR_int_std.rds"))
	
	post_GP_CFAR_int_std <- as.data.frame(GP_CFAR_int_std)
	post_GP_CFAR_int_std_mean <- stack(apply(post_GP_CFAR_int_std, 2, mean))
	post_GP_CFAR_int_std_sd <- stack(apply(post_GP_CFAR_int_std, 2, sd))
	post_GP_CFAR_int_std_LCI <- stack(lapply(post_GP_CFAR_int_std, quantile, 
																	prob = (0.025))) %>%
									    						rename(low_CI = values)
	post_GP_CFAR_int_std_HCI <- stack(lapply(post_GP_CFAR_int_std, quantile, 
																	prob = (0.975))) %>%
									    					  rename(high_CI = values)
	post_GP_CFAR_int_std_sum <- merge(post_GP_CFAR_int_std_mean, post_GP_CFAR_int_std_LCI) %>% 
														  merge(post_GP_CFAR_int_std_HCI) %>% round_df(2)
	
	fortable <- post_GP_CFAR_int_std_sum[grep("bGP", post_GP_CFAR_int_std_sum$ind), ]
	fortable <- post_GP_CFAR_int_std_sum[grep("aGP", post_GP_CFAR_int_std_sum$ind), ]

	# proportion > 0
	GSAint_CFAR_GSAint_prob0 <- (length(which(post_GP_CFAR_int_std$bGP_int >0)))/
		length(post_GP_CFAR_int_std$bGP_int)

	GSAint_CFAR_CFAR_prob0 <- (length(which(post_GP_CFAR_int_std$bGP_CFAR >0)))/
		length(post_GP_CFAR_int_std$bGP_CFAR)

	#### GP ~ slope + CFAR ####

	RawGSA8_CFAR$id = as.numeric(factor(RawGSA8_CFAR$Binomial))

	dat <- list(
							N=nrow(RawGSA8_CFAR),
							J=len(RawGSA8_CFAR$Binomial),
							sp=RawGSA8_CFAR$id,
							LogGSAcm2=RawGSA8_CFAR$LogGSAcm2,
							LogCenterMassG = RawGSA8_CFAR$LogCenterMassG,
							CFAR = unique(RawGSA8_CFAR$meanCFARstd),
							GP = unique(RawGSA8_CFAR$mean_gperf_simple))
	
	GP_CFAR_slope_std <- stan(
	                file = here("./stan models/MLM_CFAR_slope.stan"),
	                data = dat,
	                iter = 5000,
	                warmup = 2000,
	                chains = 4,
	                control = list("adapt_delta" = 0.99,
	                							 "max_treedepth" = 18),
	                pars = c( "global_slope",
	              		  			"global_int",
	              		  			"aGP",
	              		  			"bGP_slope",
	                					"bGP_CFAR",
	              		        "beta_int",
	              		  			"beta_slope",
	              		  			"beta_slopes_std",
	              		        "sigma",
	              		  			"sigma_GP",
	              		  			"log_lik"))
	
	#save(GP_CFAR_slope_std, file = ("./output/GP_CFAR_slope_std.rds"))
	#load(file = here("./output/GP_CFAR_slope_std.rds")) 
	
	post_GP_CFAR_slope_std <- as.data.frame(GP_CFAR_slope_std)
	post_GP_CFAR_slope_std_mean <- stack(apply(post_GP_CFAR_slope_std, 2, mean))
	post_GP_CFAR_slope_std_sd <- stack(apply(post_GP_CFAR_slope_std, 2, sd))
	post_GP_CFAR_slope_std_LCI <- stack(lapply(post_GP_CFAR_slope_std, quantile, 
																	prob = (0.025))) %>%
									    						rename(low_CI = values)
	post_GP_CFAR_slope_std_HCI <- stack(lapply(post_GP_CFAR_slope_std, quantile, 
																	prob = (0.975))) %>%
									    					  rename(high_CI = values)
	post_GP_CFAR_slope_std_sum <- merge(post_GP_CFAR_slope_std_mean, post_GP_CFAR_slope_std_HCI) %>% 
														    merge(post_GP_CFAR_slope_std_LCI) %>% round_df(2)
	
	fortable <- post_GP_CFAR_slope_std_sum[grep("bGP", post_GP_CFAR_slope_std_sum$ind), ]
	fortable <- post_GP_CFAR_slope_std_sum[grep("aGP", post_GP_CFAR_slope_std_sum$ind), ]

	# proportion > 0
	GSAslope_CFAR_GSAslope_prob0 <- (length(which(post_GP_CFAR_slope_std$bGP_slope >0)))/
		length(post_GP_CFAR_slope_std$bGP_slope)

	GSAslope_CFAR_CFAR_prob0 <- (length(which(post_GP_CFAR_slope_std$bGP_CFAR >0)))/
		length(post_GP_CFAR_slope_std$bGP_CFAR)

	#### with phylo ####
	
	### which species on the tree?
	sp_CFAR <- unique(RawGSA8_CFAR$phylo)
	sp_drop <- setdiff(tree_pruned$tip.label, sp_CFAR)
	tree_pruned_ex <- drop.tip(tree_pruned, sp_drop)

	## prep tree
	library(here)
	library(MCMCglmm)
	library(phytools)
	library(phangorn)

	# identity matrix for stan
	d_mat <- diag(1, 30, 30)
	
	# prepare tree for stan
	inv.tree <- inverseA(tree_pruned_ex, nodes = "TIPS", scale = TRUE) #works
	A <- solve(inv.tree$Ainv) # is this reversing the inverse so the matrix is in the form of a covariance matrix and not its inverse (required by brms)
	rownames(A) <- rownames(inv.tree$Ainv) # assigning rownames of species ( = to the tips of the tree) to the covariance matrix
	
	vcov_mat <- as.matrix(A)

	### run models with phylogeny ####
	
	## run stan models
	
	#### GP ~ intercept + CFAR 
	
	RawGSA8_CFAR$id = as.numeric(factor(RawGSA8_CFAR$Binomial))
	
	dat <- list(
							N=nrow(RawGSA8_CFAR),
							J=len(RawGSA8_CFAR$Binomial),
							sp=RawGSA8_CFAR$id,
							LogGSAcm2=RawGSA8_CFAR$LogGSAcm2,
							LogCenterMassG = RawGSA8_CFAR$LogCenterMassG,
							CFAR = unique(RawGSA8_CFAR$meanCFARstd),
							d_mat = d_mat,
							vcov_mat = vcov_mat,
							GP = unique(RawGSA8_CFAR$mean_gperf_simple))
	
	GP_CFAR_int_std_phylo <- stan(
	                file = here("./stan models/MLM_CFAR_INT_PHYLO.stan"),
	                data = dat,
	                iter = 5000,
	                warmup = 2000,
	                chains = 4,
	                control = list("adapt_delta" = 0.99,
	                							 "max_treedepth" = 18),
	                pars = c( "global_slope",
	              		  			"global_int",
	              		  			"aGP",
	              		  			"bGP_int",
	                					"bGP_CFAR",
	              		        "beta_int",
	              		  			"beta_slope",
	              		  			"beta_ints_std",
	              		        "sigma",
	              		  			"sigma_GP",
	              		  			"log_lik",
	                					"lambda_GP"))

	#save(GP_CFAR_int_std_phylo, file = here("./output/GP_CFAR_int_std_phylo.rds"))
	#load(file = here("./output/GP_CFAR_int_std_phylo.rds"))
	
	post_GP_CFAR_int_std_phylo <- as.data.frame(GP_CFAR_int_std_phylo)
	post_GP_CFAR_int_std_phylo_mean <- stack(apply(post_GP_CFAR_int_std_phylo, 2, mean))
	post_GP_CFAR_int_std_phylo_sd <- stack(apply(post_GP_CFAR_int_std_phylo, 2, sd))
	post_GP_CFAR_int_std_phylo_LCI <- stack(lapply(post_GP_CFAR_int_std_phylo, quantile, 
																	prob = (0.025))) %>%
									    						rename(low_CI = values)
	post_GP_CFAR_int_std_phylo_HCI <- stack(lapply(post_GP_CFAR_int_std_phylo, quantile, 
																	prob = (0.975))) %>%
									    					  rename(high_CI = values)
	post_GP_CFAR_int_std_phylo_sum <- merge(post_GP_CFAR_int_std_phylo_mean, 
																					post_GP_CFAR_int_std_phylo_LCI) %>% 
														  merge(post_GP_CFAR_int_std_phylo_HCI) %>% round_df(2)



	#### GP ~ slope + CFAR ####
	
	RawGSA8_CFAR$id = as.numeric(factor(RawGSA8_CFAR$Binomial))
	
	dat <- list(
							N=nrow(RawGSA8_CFAR),
							J=len(RawGSA8_CFAR$Binomial),
							sp=RawGSA8_CFAR$id,
							LogGSAcm2=RawGSA8_CFAR$LogGSAcm2,
							LogCenterMassG = RawGSA8_CFAR$LogCenterMassG,
							CFAR = unique(RawGSA8_CFAR$meanCFARstd),
							vcov_mat = vcov_mat,
							d_mat = d_mat,
							GP = unique(RawGSA8_CFAR$mean_gperf_simple))
	
	GP_CFAR_slope_std_phylo <- stan(
	                file = here("./stan models/MLM_CFAR_SLOPE_PHYLO.stan"),
	                data = dat,
	                iter = 5000,
	                warmup = 2000,
	                chains = 4,
	                control = list("adapt_delta" = 0.99,
	                							 "max_treedepth" = 18),
	                pars = c( "global_slope",
	              		  			"global_int",
	              		  			"aGP",
	              		  			"bGP_slope",
	                					"bGP_CFAR",
	              		        "beta_int",
	              		  			"beta_slope",
	              		  			"beta_slopes_std",
	              		        "sigma",
	              		  			"sigma_GP",
	                					"lambda_GP",
	              		  			"log_lik"))
	
	#save(GP_CFAR_slope_std_phylo, file = ("./output/GP_CFAR_slope_std_phylo.rds"))
	#load(file = here("./output/GP_CFAR_slope_std_phylo.rds")) 
	
	post_GP_CFAR_slope_std_phylo <- as.data.frame(GP_CFAR_slope_std_phylo)
	post_GP_CFAR_slope_std_phylo_mean <- stack(apply(post_GP_CFAR_slope_std_phylo, 2, mean))
	post_GP_CFAR_slope_std_phylo_sd <- stack(apply(post_GP_CFAR_slope_std_phylo, 2, sd))
	post_GP_CFAR_slope_std_phylo_LCI <- stack(lapply(post_GP_CFAR_slope_std_phylo, quantile, 
																	prob = (0.025))) %>%
									    						rename(low_CI = values)
	post_GP_CFAR_slope_std_phylo_HCI <- stack(lapply(post_GP_CFAR_slope_std_phylo, quantile, 
																	prob = (0.975))) %>%
									    					  rename(high_CI = values)
	post_GP_CFAR_slope_std_phylo_sum <- merge(post_GP_CFAR_slope_std_phylo_mean,
																						post_GP_CFAR_slope_std_phylo_HCI) %>% 
														    merge(post_GP_CFAR_slope_std_phylo_LCI) %>% round_df(2)
	


	#### without aquaculture species ####

	# which species in the raw dataset of 30 CFAR species are used in aquaculture?
	RawGSA8_CFAR_no_ac <- filter(RawGSA8_CFAR, Binomial %notin% aquaculture_sp_mean)
	len(RawGSA8_CFAR_no_ac$Binomial) # n = 27 species

	#### GP ~ intercept + CFAR ####
	
	RawGSA8_CFAR$id = as.numeric(factor(RawGSA8_CFAR$Binomial))
	
	dat <- list(
							N=nrow(RawGSA8_CFAR),
							J=len(RawGSA8_CFAR$Binomial),
							sp=RawGSA8_CFAR$id,
							LogGSAcm2=RawGSA8_CFAR$LogGSAcm2,
							LogCenterMassG = RawGSA8_CFAR$LogCenterMassG,
							CFAR = unique(RawGSA8_CFAR$meanCFARstd),
							GP = unique(RawGSA8_CFAR$mean_gperf_simple))
	
	GP_CFAR_int_std_ex <- stan(
	                file = ("./stan models/MLM_CFAR_INT.stan"),
	                data = dat,
	                iter = 5000,
	                warmup = 2000,
	                chains = 4,
	                control = list("adapt_delta" = 0.99,
	                							 "max_treedepth" = 18),
	                pars = c( "global_slope",
	              		  			"global_int",
	              		  			"aGP",
	              		  			"bGP_int",
	                					"bGP_CFAR",
	              		        "beta_int",
	              		  			"beta_slope",
	              		  			"beta_ints_std",
	              		        "sigma",
	              		  			"sigma_GP",
	              		  			"log_lik"))
	
	#save(GP_CFAR_int_std_ex, file = ("./outputs/GP_CFAR_int_std_ex.rds"))
	#load(file = here("./output/GP_CFAR_int_std_ex.rds"))
	
	post_GP_CFAR_int_std_no_ac <- as.data.frame(GP_CFAR_int_std_ex)
	post_GP_CFAR_int_std_no_ac_mean <- stack(apply(post_GP_CFAR_int_std_no_ac, 2, mean))
	post_GP_CFAR_int_std_no_ac_sd <- stack(apply(post_GP_CFAR_int_std_no_ac, 2, sd))
	post_GP_CFAR_int_std_no_ac_LCI <- stack(lapply(post_GP_CFAR_int_std_no_ac, quantile, 
																	prob = (0.025))) %>%
									    						rename(low_CI = values)
	post_GP_CFAR_int_std_no_ac_HCI <- stack(lapply(post_GP_CFAR_int_std_no_ac, quantile, 
																	prob = (0.975))) %>%
									    					  rename(high_CI = values)
	post_GP_CFAR_int_std_no_ac_sum <- merge(post_GP_CFAR_int_std_no_ac_mean, 
																					post_GP_CFAR_int_std_no_ac_LCI) %>% 
														  merge(post_GP_CFAR_int_std_no_ac_HCI) %>% round_df(2)
	fortable <- post_GP_CFAR_int_std_sum[grep("bGP", post_GP_CFAR_int_std_sum$ind), ]
	fortable <- post_GP_CFAR_int_std_sum[grep("aGP", post_GP_CFAR_int_std_sum$ind), ]


	#### GP ~ slope + CFAR ####

	RawGSA8_CFAR$id = as.numeric(factor(RawGSA8_CFAR$Binomial))
	
	dat <- list(
							N=nrow(RawGSA8_CFAR),
							J=len(RawGSA8_CFAR$Binomial),
							sp=RawGSA8_CFAR$id,
							LogGSAcm2=RawGSA8_CFAR$LogGSAcm2,
							LogCenterMassG = RawGSA8_CFAR$LogCenterMassG,
							CFAR = unique(RawGSA8_CFAR$meanCFARstd),
							GP = unique(RawGSA8_CFAR$mean_gperf_simple))
	
	GP_CFAR_slope_std_ex <- stan(
	                file = ("./stan models/MLM_CFAR_slope.stan"),
	                data = dat,
	                iter = 5000,
	                warmup = 2000,
	                chains = 4,
	                control = list("adapt_delta" = 0.99,
	                							 "max_treedepth" = 18),
	                pars = c( "global_slope",
	              		  			"global_int",
	              		  			"aGP",
	              		  			"bGP_slope",
	                					"bGP_CFAR",
	              		        "beta_int",
	              		  			"beta_slope",
	              		  			"beta_slopes_std",
	              		        "sigma",
	              		  			"sigma_GP",
	              		  			"log_lik"))
	
	#save(GP_CFAR_slope_std_ex, file = ("./output/GP_CFAR_slope_std_ex.rds"))
	#load(file = here("./output/GP_CFAR_slope_std_ex.rds"))
	
	post_GP_CFAR_slope_std_no_ac <- as.data.frame(GP_CFAR_slope_std_ex)
	post_GP_CFAR_slope_std_no_ac_mean <- stack(apply(post_GP_CFAR_slope_std_no_ac, 2, mean))
	post_GP_CFAR_slope_std_no_ac_sd <- stack(apply(post_GP_CFAR_slope_std_no_ac, 2, sd))
	post_GP_CFAR_slope_std_no_ac_LCI <- stack(lapply(post_GP_CFAR_slope_std_no_ac, quantile, 
																	prob = (0.025))) %>%
									    						rename(low_CI = values)
	post_GP_CFAR_slope_std_no_ac_HCI <- stack(lapply(post_GP_CFAR_slope_std_no_ac, quantile, 
																	prob = (0.975))) %>%
									    					  rename(high_CI = values)
	post_GP_CFAR_slope_std_no_ac_sum <- merge(post_GP_CFAR_slope_std_no_ac_mean, 
																						post_GP_CFAR_slope_std_no_ac_HCI) %>% 
														    merge(post_GP_CFAR_slope_std_no_ac_LCI) %>% round_df(2)
	


	#### without airbreathing species ####

	# which species in the raw dataset of 30 CFAR species are used in aquaculture?
	RawGSA8_CFAR_no_abr <- filter(RawGSA8_CFAR, Binomial %notin% airbreathing_mean)
	len(RawGSA8_CFAR_no_abr$Binomial) # n = 25 species

	#### GP ~ intercept + CFAR 
	
	RawGSA8_CFAR_no_abr$id = as.numeric(factor(RawGSA8_CFAR_no_abr$Binomial))
	
	dat <- list(
							N=nrow(RawGSA8_CFAR_no_abr),
							J=len(RawGSA8_CFAR_no_abr$Binomial),
							sp=RawGSA8_CFAR_no_abr$id,
							LogGSAcm2=RawGSA8_CFAR_no_abr$LogGSAcm2,
							LogCenterMassG = RawGSA8_CFAR_no_abr$LogCenterMassG,
							CFAR = unique(RawGSA8_CFAR_no_abr$meanCFARstd),
							GP = unique(RawGSA8_CFAR_no_abr$mean_gperf_simple))
	
	GP_CFAR_int_std_no_abr <- stan(
	                file = ("./stan models/MLM_CFAR_INT.stan"),
	                data = dat,
	                iter = 5000,
	                warmup = 2000,
	                chains = 4,
	                control = list("adapt_delta" = 0.99,
	                							 "max_treedepth" = 18),
	                pars = c( "global_slope",
	              		  			"global_int",
	              		  			"aGP",
	              		  			"bGP_int",
	                					"bGP_CFAR",
	              		        "beta_int",
	              		  			"beta_slope",
	              		  			"beta_ints_std",
	              		        "sigma",
	              		  			"sigma_GP",
	              		  			"log_lik"))

	#save(GP_CFAR_int_std_no_abr, file = ("./output/GP_CFAR_int_std_no_abr.rds"))
	#load(file = here("./output/GP_CFAR_int_std_no_abr.rds"))
	
	post_GP_CFAR_int_std_no_abr <- as.data.frame(GP_CFAR_int_std_no_abr)
	post_GP_CFAR_int_std_no_abr_mean <- stack(apply(post_GP_CFAR_int_std_no_abr, 2, mean))
	post_GP_CFAR_int_std_no_abr_sd <- stack(apply(post_GP_CFAR_int_std_no_abr, 2, sd))
	post_GP_CFAR_int_std_no_abr_LCI <- stack(lapply(post_GP_CFAR_int_std_no_abr, quantile, 
																	prob = (0.025))) %>%
									    						rename(low_CI = values)
	post_GP_CFAR_int_std_no_abr_HCI <- stack(lapply(post_GP_CFAR_int_std_no_abr, quantile, 
																	prob = (0.975))) %>%
									    					  rename(high_CI = values)
	post_GP_CFAR_int_std_no_abr_sum <- merge(post_GP_CFAR_int_std_no_abr_mean, 
																					post_GP_CFAR_int_std_no_abr_LCI) %>% 
														  merge(post_GP_CFAR_int_std_no_abr_HCI) %>% round_df(2)


	#### GP ~ slope + CFAR ####
	
	RawGSA8_CFAR_no_abr$id = as.numeric(factor(RawGSA8_CFAR_no_abr$Binomial))
	
	dat <- list(
							N=nrow(RawGSA8_CFAR_no_abr),
							J=len(RawGSA8_CFAR_no_abr$Binomial),
							sp=RawGSA8_CFAR_no_abr$id,
							LogGSAcm2=RawGSA8_CFAR_no_abr$LogGSAcm2,
							LogCenterMassG = RawGSA8_CFAR_no_abr$LogCenterMassG,
							CFAR = unique(RawGSA8_CFAR_no_abr$meanCFARstd),
							GP = unique(RawGSA8_CFAR_no_abr$mean_gperf_simple))
	
	GP_CFAR_slope_std_no_abr <- stan(
	                file = ("./stan models/MLM_CFAR_slope.stan"),
	                data = dat,
	                iter = 5000,
	                warmup = 2000,
	                chains = 4,
	                control = list("adapt_delta" = 0.99,
	                							 "max_treedepth" = 18),
	                pars = c( "global_slope",
	              		  			"global_int",
	              		  			"aGP",
	              		  			"bGP_slope",
	                					"bGP_CFAR",
	              		        "beta_int",
	              		  			"beta_slope",
	              		  			"beta_slopes_std",
	              		        "sigma",
	              		  			"sigma_GP",
	              		  			"log_lik"))

	#save(GP_CFAR_slope_std_no_abr, file = ("./output/GP_CFAR_slope_std_no_abr.rds"))
	#load(file = here("./output/GP_CFAR_slope_std_no_abr.rds"))
	
	post_GP_CFAR_slope_std_no_abr <- as.data.frame(GP_CFAR_slope_std_no_abr)
	post_GP_CFAR_slope_std_no_abr_mean <- stack(apply(post_GP_CFAR_slope_std_no_abr, 2, mean))
	post_GP_CFAR_slope_std_no_abr_sd <- stack(apply(post_GP_CFAR_slope_std_no_abr, 2, sd))
	post_GP_CFAR_slope_std_no_abr_LCI <- stack(lapply(post_GP_CFAR_slope_std_no_abr, quantile, 
																	prob = (0.025))) %>%
									    						rename(low_CI = values)
	post_GP_CFAR_slope_std_no_abr_HCI <- stack(lapply(post_GP_CFAR_slope_std_no_abr, quantile, 
																	prob = (0.975))) %>%
									    					  rename(high_CI = values)
	post_GP_CFAR_slope_std_no_abr_sum <- merge(post_GP_CFAR_slope_std_no_abr_mean, 
																						post_GP_CFAR_slope_std_no_abr_HCI) %>% 
														    merge(post_GP_CFAR_slope_std_no_abr_LCI) %>% round_df(2)
	
	
