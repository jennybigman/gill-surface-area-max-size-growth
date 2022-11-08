# figure 1 -- plots of GP & slope and int WITHOUT PHYLO

	## plot prep ####
	
	#### GP ~ intercept no phylo ####

	# load model that has GP ~ intercept (intercepts were standardized after first level)
	load(file = here("./output/MLM_GP_int_std.rds"))

	# extract posterior samples and summarize
	post_MLM_GP_int_std <- as.data.frame(MLM_GP_int_std)

	# summary of output: mean effect sizes, sds of effect sizes, & 95% Bayesian Credible Intervals
	post_MLM_GP_int_std_mean <- stack(apply(post_MLM_GP_int_std, 2, mean)) %>% rename(mean = values)
	post_MLM_GP_int_std_LCIs <- stack(apply(post_MLM_GP_int_std, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_MLM_GP_int_std_HCIs <- stack(apply(post_MLM_GP_int_std, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_MLM_GP_int_std_sum <- merge(post_MLM_GP_int_std_mean, post_MLM_GP_int_std_LCIs, by = "ind") %>%
												     merge(post_MLM_GP_int_std_HCIs) 
	
	# all intercepts and slopes are reported in the output in terms of how different
	# the value is to the reference category, so need to add those to the reference
	# categories
	
	#names(post_MLM_GP_int_std)
	
	beta_ints <- post_MLM_GP_int_std[[2]] + post_MLM_GP_int_std[5:36]
	
	beta_ints_std_mod <- post_MLM_GP_int_std[69:100]
	# standardize beta ints because that is what were

	## for plotting

	# create vectors of the correct name for each variable
	species_list <- sort(unique(RawGSA8_phylo$Binomial))
	
	# add a "_" between genus and species
	Species <- lapply(strsplit(as.character(species_list), "\\ "), "[", 2)
	Genus <- lapply(strsplit(as.character(species_list), "\\ "), "[", 1)
	species_list <- tibble(species_list, Genus, Species)
	species_list$Species2 <- paste(species_list$Genus, species_list$Species,
	                                sep = "_")
	species_list <- species_list$Species2
	
	# change the column names to the species name
	names(beta_ints_std_mod) <- species_list

	# put coefficients in a dataframe together for plotting
	
	# intercepts
	std_beta_ints_mean <- stack(apply(beta_ints_std_mod, 2, mean)) %>%
											  rename(species = ind,
											 			   intercept = values)
	

	# get growth performance and add into dataset for plotting
	GP_species <- distinct(RawGSA8_phylo, mean_gperf_simple, .keep_all = TRUE) %>%
								dplyr::select(Binomial, mean_gperf_simple) %>% ungroup()
	
	GP_species$Species <- lapply(strsplit(as.character(GP_species$Binomial), "\\ "), "[", 2)
	GP_species$Genus <- lapply(strsplit(as.character(GP_species$Binomial), "\\ "), "[", 1)
	GP_species$species <- paste(GP_species$Genus, GP_species$Species,
	                                sep = "_")
	
	GP_species <- GP_species %>% dplyr::select(species, mean_gperf_simple)

	std_beta_ints_mean <- merge(std_beta_ints_mean, GP_species)

	ggplot() +
	geom_point(data = std_beta_ints_mean, 
		aes(x = intercept, y = mean_gperf_simple))

	# fit and prediction lines #

	# fitted line
	intercept_obs <- sort(std_beta_ints_mean$intercept)
	
	fitted_matrix_intercept <- t(sapply(1:nrow( as.data.frame(MLM_GP_int_std)), function(i) {
	                             post_MLM_GP_int_std[i, 'aGP_int'] +
	                             post_MLM_GP_int_std[i, 'bGP_int'] * intercept_obs 
		}))
	
	mean_fitted_matrix_intercept <- colMeans(fitted_matrix_intercept)
	
	BCI_int <- t(apply(fitted_matrix_intercept, 2, quantile,
	                    prob = c(0.025, 0.975))) %>%
	            as.data.frame() %>%
	            rename(lower_BCI = "2.5%",
	                   upper_BCI = "97.5%") %>%
	            mutate(intercept_obs = intercept_obs,
	                   fitted_GP = mean_fitted_matrix_intercept)


	#### GP ~ slope no phylo ####
	
	load(file = here("./output/MLM_GP_slope_std.rds"))
	post_MLM_GP_slope_std <- as.data.frame(MLM_GP_slope_std)
	names(post_MLM_GP_slope_std)
	
	# summary of output: mean effect sizes, sds of effect sizes, & 95% Bayesian Credible Intervals
	post_MLM_GP_slope_std_mean <- stack(apply(post_MLM_GP_slope_std, 2, mean)) %>% rename(mean = values)
	post_MLM_GP_slope_std_LCIs <- stack(apply(post_MLM_GP_slope_std, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_MLM_GP_slope_std_HCIs <- stack(apply(post_MLM_GP_slope_std, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_MLM_GP_slope_std_sum <- merge(post_MLM_GP_slope_std_mean, post_MLM_GP_slope_std_LCIs, by = "ind") %>%
												    	 merge(post_MLM_GP_slope_std_HCIs) 
	
	beta_slopes_std <- post_MLM_GP_slope_std[69:100]
	
	## for plotting
	
	# change the column names to the species name
	names(beta_slopes_std) <- species_list
	
	# put coefficients in a dataframe together for plotting
	
	# slopes
	beta_slopes_std_mean <- stack(apply(beta_slopes_std, 2, mean)) %>%
											  rename(species = ind,
											 			   slope = values)
	
	# get growth performance and add into dataset for plotting
	beta_slopes_std_mean <- merge(beta_slopes_std_mean, GP_species)
	
	ggplot(data = beta_slopes_std_mean, 
		aes(x = slope, y = mean_gperf_simple)) +
		geom_point()
	
	# fit and prediction lines #
	
	# fitted line
	slope_obs <- sort(beta_slopes_std_mean$slope)
	
	fitted_matrix_slope <- t(sapply(1:nrow( as.data.frame(MLM_GP_slope_std)), function(i) {
	                             post_MLM_GP_slope_std[i, 'aGP_slope'] +
	                             post_MLM_GP_slope_std[i, 'bGP_slope'] * slope_obs 
		}))
	
	mean_fitted_matrix_slope <- colMeans(fitted_matrix_slope)
	
	BCI_slope <- t(apply(fitted_matrix_slope, 2, quantile,
	                    prob = c(0.025, 0.975))) %>%
	            as.data.frame() %>%
	            rename(lower_BCI = "2.5%",
	                   upper_BCI = "97.5%") %>%
	            mutate(slope_obs = slope_obs,
	                   fitted_GP = mean_fitted_matrix_slope)
	

	#### scaling plots ####
	
	# intercept no phylo #
	
	GP_int_nophylo <- 
		ggplot(BCI_int) +
		geom_ribbon(data = BCI_int, 
			aes(x = intercept_obs,	ymin = lower_BCI,
      ymax = upper_BCI),	fill = "grey", alpha = 0.3) +
		geom_line(data = BCI_int, aes(x = intercept_obs, y = fitted_GP),
     color = "darkgrey", size = 1.25, alpha = 0.6) +
    geom_point(data = std_beta_ints_mean,
     aes(x = intercept,y = mean_gperf_simple),
     alpha = 0.6, size = 4) +
		annotate("text", x = -2.7 , y = 5.4,
     label = "(a)",size = 5) +
		scale_y_continuous(
     name = "Growth performance",
     breaks = c(0, 1, 2, 3, 4, 5),
     labels = c(0, 1, 2, 3, 4, 5),
     limits = c(0, 5.4)) +
		scale_x_continuous(
     name = "Ontogenetic intercept",
     breaks = c(-2, 0, 1),
     labels = c(-2, 0, 1),
     limits = c(-2.70, 1.53)) +
		theme(legend.position = "none") +
		theme_bw() +
  	theme(
  	  axis.text = element_text(size= 10, colour = "black"),
  	  axis.title = element_text(size=12),
  	  panel.grid.major = element_blank(),
  	  panel.grid.minor = element_blank(),
  	  panel.border = element_rect(linetype = 1, size= 1, fill = NA))


	# slope no phylo #

	GP_slope_nophylo <- 
		ggplot(BCI_slope) +
		geom_ribbon(data = BCI_slope, 
			aes(x = slope_obs, ymin = lower_BCI,
      ymax = upper_BCI), fill = "grey", alpha = 0.3) +
		geom_line(data = BCI_slope,
      aes(x = slope_obs, y = fitted_GP),
      color = "darkgrey", size = 1.25, alpha = 0.6) +
		geom_point(data = beta_slopes_std_mean,
       aes(x = slope, y = mean_gperf_simple),
       alpha = 0.6, size = 4) +
		annotate("text", x = -2 , y = 5.4,
      label = "(b)", size = 5) +
		scale_y_continuous(
      name = "Growth performance",
      breaks = c(0, 1, 2, 3, 4, 5),
      labels = c(0, 1, 2, 3, 4, 5),
      limits = c(0, 5.4)) +
		scale_x_continuous(
      name = "Ontogenetic slope",
      breaks = c(-2, 0, 2),
      labels = c(-2, 0, 2),
      limits = c(-2, 2)) +
    theme_bw() +
  	theme(
  	legend.position = "none",
    axis.text.x = element_text(size= 10, colour = "black"),
  	axis.title.x = element_text(size=12),
  	axis.title.y = element_blank(),
  	axis.text.y = element_blank(),
  	axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linetype = 1, size= 1, fill = NA))


	
	#### posterior distribution plots ####
	
	# intercept no phylo #

	post_MLM_GP_int_std_sum_slope <- post_MLM_GP_int_std_sum %>% 
		filter(ind == "bGP_int")

	GP_int_dist <- 
		ggplot(data = post_MLM_GP_int_std, aes(bGP_int)) +
		geom_density(color = "black", fill = "grey50", alpha = 0.5) +
		scale_x_continuous(
			breaks = c(0, 1),
			labels = c(0, 1)) +
		geom_vline(aes(xintercept = 0), color = "black") +
		theme_bw() +
		theme(
			axis.text.x = element_text(size = 6),
			axis.ticks.length.x = unit(-0.2, "cm"),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			axis.title = element_blank(),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank(),
		  panel.border = element_rect(color = "black",
		  	size = 0.5, fill = NA),
			plot.margin = unit(c(0, 0, 0, 0), "null"))
	
	# shaded by proportion > 0 
	post_MLM_GP_int_std <- as.data.frame(MLM_GP_int_std)
	GP_int_slope_prob0 <- (length(which(post_MLM_GP_int_std$bGP_int >0)))/length(post_MLM_GP_int_std$bGP_int) * 100

  # create df of density of posterior to plot by hand
	GP_int_dist <- 
		ggplot(data = post_MLM_GP_int_std, aes(bGP_int)) +
		geom_density(color = "black", fill = "grey50", alpha = 0.5)
	
	GP_int_dist_df <- ggplot_build(GP_int_dist)$data[[1]] 
	
	GSA_int_dist_base <- 
		ggplot() + 
		geom_area(data = subset(GP_int_dist_df, x < 0),
                         aes(x=x,y=y),
                         fill = "#cccccc",
                         color = "black", 
												 alpha = 0.5) +
		geom_area(data = subset(GP_int_dist_df, x > 0),
                         aes(x=x,y=y),
                         fill = "#4c4c4c",
                         color = "black", 
												 alpha = 0.5) +
	geom_segment(aes(x = -0.792, xend = 1.27, y = 0, yend = 0), color = "black", size = 0.5) +
		scale_x_continuous(
			breaks = c(0, 1),
			labels = c(0, 1)) +
		annotate("text", x = 0.39, y = 0.15,
		          label = paste0(GP_int_slope_prob0, "%"),
		          size = 3) +
		theme_bw() +
		theme(
			axis.text.x = element_text(size = 6),
			axis.ticks.length.x = unit(-0.2, "cm"),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			axis.title = element_blank(),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank(),
		  panel.border = element_rect(color = "black",
		  	size = 0.5, fill = NA),
			plot.margin = unit(c(0, 0, 0, 0), "null"))

		#ggplot_build(GSA_int_dist_base)$layout$panel_scales_x[[1]]$range$range


	# slope no phylo
	
	post_MLM_GP_slope_std_sum_slope <- post_MLM_GP_slope_std_sum %>%
		filter(ind == "bGP_slope")

	GP_slope_dist <- 
		ggplot(data = post_MLM_GP_slope_std, aes(bGP_slope)) +
    geom_density(color = "black", fill = "grey50",  alpha = 0.5) +
    geom_segment(data = post_MLM_GP_slope_std_sum_slope,
    	aes(x = LCI, xend = HCI, y = 0, yend = 0),
      size = 0.5, color = "black") +
    geom_point(data = post_MLM_GP_slope_std_sum_slope,
  		aes(x = mean, y = 0), size = 3) +
    annotate("text", x = post_MLM_GP_slope_std_sum_slope$mean, y = 0.35,
      label = post_MLM_GP_slope_std_sum_slope$mean, size = 3) +
    scale_x_continuous(
			breaks = c(0, 1),
			labels = c(0, 1)) +
		theme_bw() +
		theme(
			axis.text.x = element_text(size = 6),
			axis.ticks.length.x = unit(-0.2, "cm"),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			axis.title = element_blank(),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank(),
		  panel.border = element_rect(color = "black",
		                 size = 0.5, fill = NA),
			 plot.margin = unit(c(0, 0, 0, 0), "null"))
	
	# shaded by proportion > 0 
	 
	# create df of density of posterior to plot by hand
	post_MLM_GP_slope_std <- as.data.frame(MLM_GP_slope_std)
	GP_slope_slope_prob0 <- (length(which(post_MLM_GP_slope_std$bGP_slope >0)))/length(post_MLM_GP_slope_std$bGP_slope) * 100

	GP_slope_dist <- 
		ggplot(data = post_MLM_GP_slope_std, aes(bGP_slope)) +
		geom_density(color = "black", fill = "grey50", alpha = 0.5)
	
	GP_slope_dist_df <- ggplot_build(GP_slope_dist)$data[[1]] 
	
	GSA_slope_dist_base <- 
		ggplot() + 
		geom_area(data = subset(GP_slope_dist_df, x < 0),
                         aes(x=x,y=y),
                         fill = "#cccccc",
                         color = "black", alpha = 0.5) +
		geom_area(data = subset(GP_slope_dist_df, x > 0),
                         aes(x=x,y=y),
                         fill = "#4c4c4c",
                         color = "black", alpha = 0.5) +
		scale_x_continuous(
			breaks = c(0, 1),
			labels = c(0, 1)) +
		annotate("text", x = 0.38, y = 0.15,
		          label = paste0(GP_slope_slope_prob0, "%"),
		          size = 3) +
		geom_segment(aes(x = -0.726, xend = 1.230, y = 0, yend = 0), color = "black") +
		theme_bw() +
		theme(
			axis.text.x = element_text(size = 6),
			axis.ticks.length.x = unit(-0.2, "cm"),
			axis.ticks.y = element_blank(),
			axis.text.y = element_blank(),
			axis.title = element_blank(),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank(),
		  panel.border = element_rect(color = "black",
		  	size = 0.5, fill = NA),
			plot.margin = unit(c(0, 0, 0, 0), "null"))
	
		#ggplot_build(GP_slope_dist)$layout$panel_scales_x[[1]]$range$range


	## plot together ####
	
	plot1 <- GP_int_nophylo + 				 
		       ggtitle("Intercept") + theme(plot.title = element_text(size = 12, hjust = 0.5,  face = "bold")) +
					 theme(plot.margin = unit(c(0.2, -0.05, 0, 0.2), "in"))
	
	plot2 <- GP_slope_nophylo + 
			     ggtitle("Slope") + theme(plot.title = element_text(size = 12, hjust = 0.5,  face = "bold")) +
					 theme(plot.margin = unit(c(0.2, 0.2, 0, -0.05), "in"))

	plot3 <- GSA_int_dist_base + theme(
		panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
		panel.border = element_blank())

	plot4 <- GSA_slope_dist_base  + theme(
		panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
		panel.border = element_blank())
	
	plot1_inset <- plot1 + inset_element(plot3, 0.7, 0.005, 1, 0.3, clip = TRUE)

	plot2_inset <- plot2 + inset_element(plot4,0.7, 0.005, 1, 0.3, clip = TRUE)
	
	
	Fig_gsa_growth <- plot1_inset + plot2_inset + plot_layout(ncol = 2, widths = c(1, 1))

	ggsave(Fig_gsa_growth, 
			 file = here("./output/plots/Fig_gsa_growth.png"),
									 height = 6, width = 12, units = "in")


