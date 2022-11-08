# FIGURE -- GSA,GP, CFAR

	# read in activity data and estimate CFAR

	CFAR_dat <- read.csv(here("./data/activity/Mahan_CFAR_datasheet.csv"),
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
	mean_dat <- distinct(RawGSA8_CFAR, Binomial, .keep_all = TRUE) %>%
							dplyr::select(Binomial, mean_gperf_simple, mean_CFAR, meanCFARstd)
	
	
	#### GP ~ GSA intercept, CFAR ####
	load(file = here("./output/GP_CFAR_int_std.rds"))
	post_GP_CFAR_int <- as.data.frame(GP_CFAR_int_std)
	ints <- post_GP_CFAR_int %>% dplyr::select(contains("beta_ints_std"))
	mean_dat$ints <- colMeans(ints)

	# prediction lines 

	# first quartile
	ints_mean <- (mean_dat$ints)
	fq_CFAR <- unname(quantile((mean_dat$meanCFARstd), probs = c(0.25)))
	
	fitted_matrix_fq_int <- t(sapply(1:nrow(post_GP_CFAR_int), function(i) {
	                    post_GP_CFAR_int[i, 'aGP'] +
	                    post_GP_CFAR_int[i, 'bGP_int'] * ints_mean +
					            post_GP_CFAR_int[i, 'bGP_CFAR'] * fq_CFAR 
		}))
	
	mean_fitted_matrix_fq_int <- colMeans(fitted_matrix_fq_int)
	
	BCI_int_fq <- t(apply(fitted_matrix_fq_int, 2, quantile,
	                    prob = c(0.025, 0.975))) %>%
	            as.data.frame() %>%
	            rename(lower_BCI = "2.5%",
	                   upper_BCI = "97.5%") %>%
	            mutate(ints_mean = ints_mean,
	                   fitted_GP_fq = mean_fitted_matrix_fq_int)
	
	# 3rd quartile
	tq_CFAR <- unname(quantile((mean_dat$meanCFARstd), probs = c(0.75)))
	fitted_matrix_tq_int <- t(sapply(1:nrow(post_GP_CFAR_int), function(i) {
	                    post_GP_CFAR_int[i, 'aGP'] +
	                    post_GP_CFAR_int[i, 'bGP_int'] * ints_mean +
					            post_GP_CFAR_int[i, 'bGP_CFAR'] * tq_CFAR 
		}))
	
	mean_fitted_matrix_tq_int <- colMeans(fitted_matrix_tq_int)
	
	BCI_int_tq <- t(apply(fitted_matrix_tq_int, 2, quantile,
	                    prob = c(0.025, 0.975))) %>%
	            as.data.frame() %>%
	            rename(lower_BCI = "2.5%",
	                   upper_BCI = "97.5%") %>%
	            mutate(ints_mean = ints_mean,
	                   fitted_GP_tq = mean_fitted_matrix_tq_int)
	
	# plot ####
	CFAR_int <- 
		ggplot(mean_dat) +
		geom_line(data = BCI_int_fq,
							aes(x = ints_mean, y = fitted_GP_fq),
							color = "grey", size = 1.25, alpha = 0.6) +
		geom_line(data = BCI_int_tq,
							aes(x = ints_mean, y = fitted_GP_tq),
							color = "red", size = 1.25, alpha = 0.6) +
		geom_point(data = mean_dat,
							 aes(x = ints, y = mean_gperf_simple,
							 colour = meanCFARstd), alpha = 0.6, size = 4) +
		scale_y_continuous(
			   	name = "Growth performance",
	        breaks = c(1, 3, 5),
	        labels = c(1, 3, 5),
	        limits = c(0, 5.4)) +
			annotate("text", x = -2.55 , y = 5.4,
	                         label = "(a)",
	                         size = 5) +
		scale_color_gradient(
	    low = "lightgrey",
	    high = "red") +
		xlab("Ontogenetic intercept") +
		pos1_theme() +
		theme(legend.position = "none")
	
	# GP ~ GSA slope, CFAR ####
	load(file = here("./output/GP_CFAR_slope_std.rds"))
	post_GP_CFAR_slope <- as.data.frame(GP_CFAR_slope_std)
	slopes <- post_GP_CFAR_slope %>% dplyr::select(contains("beta_slopes_std"))
	mean_dat$slopes <- colMeans(slopes)
	
	# prediction lines 
	
	# first quartile
	slopes_mean <- (mean_dat$slopes)
	fq_CFAR <- unname(quantile((mean_dat$meanCFARstd), probs = c(0.25)))
	
	fitted_matrix_fq_slope <- t(sapply(1:nrow(post_GP_CFAR_slope), function(i) {
	                    post_GP_CFAR_slope[i, 'aGP'] +
	                    post_GP_CFAR_slope[i, 'bGP_slope'] * slopes_mean +
					            post_GP_CFAR_slope[i, 'bGP_CFAR'] * fq_CFAR 
		}))
	
	mean_fitted_matrix_fq_slope <- colMeans(fitted_matrix_fq_slope)
	
	BCI_slope_fq <- t(apply(fitted_matrix_fq_slope, 2, quantile,
	                    prob = c(0.025, 0.975))) %>%
	            as.data.frame() %>%
	            rename(lower_BCI = "2.5%",
	                   upper_BCI = "97.5%") %>%
	            mutate(slopes_mean = slopes_mean,
	                   fitted_GP_fq = mean_fitted_matrix_fq_slope)
	
	# 3rd quartile
	tq_CFAR <- unname(quantile((mean_dat$meanCFARstd), probs = c(0.75)))
	fitted_matrix_tq_slope <- t(sapply(1:nrow(post_GP_CFAR_slope), function(i) {
	                    post_GP_CFAR_slope[i, 'aGP'] +
	                    post_GP_CFAR_slope[i, 'bGP_slope'] * slopes_mean +
					            post_GP_CFAR_slope[i, 'bGP_CFAR'] * tq_CFAR 
		}))
	
	mean_fitted_matrix_tq_slope <- colMeans(fitted_matrix_tq_slope)
	
	BCI_slope_tq <- t(apply(fitted_matrix_tq_slope, 2, quantile,
	                    prob = c(0.025, 0.975))) %>%
	            as.data.frame() %>%
	            rename(lower_BCI = "2.5%",
	                   upper_BCI = "97.5%") %>%
	            mutate(slopes_mean = slopes_mean,
	                   fitted_GP_tq = mean_fitted_matrix_tq_slope)
	
	# plot ####
	CFAR_slope <- 
		ggplot(mean_dat) +
		geom_line(data = BCI_slope_fq,
							aes(x = slopes_mean, y = fitted_GP_fq),
							color = "grey", size = 1.25, alpha = 0.6) +
		geom_line(data = BCI_slope_tq,
							aes(x = slopes_mean, y = fitted_GP_tq),
							color = "red", size = 1.25, alpha = 0.6) +
		geom_point(data = mean_dat,
							 aes(x = slopes, y = mean_gperf_simple,
							 colour = meanCFARstd), alpha = 0.6, size = 4) +
			scale_color_gradient(
	    low = "lightgrey",
	    high = "red") +
		annotate("text", x = -1.6 , y = 5.4,
	                         label = "(b)",
	                         size = 5) +
		xlab("Ontogenetic slope") +
			scale_y_continuous(
			   	name = "Growth performance",
	        breaks = c(1, 3, 5),
	        labels = c(1, 3, 5),
	        limits = c(0, 5.4)) +
		pos2n_theme() +
		theme(legend.position = "none") 
	
	#### coefficient plot ####
	post_int_sum <- post_GP_CFAR_int %>%
		dplyr::select(bGP_int, bGP_CFAR) %>%
		rename(bGP_int_CFAR = bGP_CFAR)
	
	post_slope_sum <- post_GP_CFAR_slope %>%
		dplyr::select(bGP_slope, bGP_CFAR) %>%
		rename(bGP_slope_CFAR = bGP_CFAR)

	posts_sum <- bind_cols(post_int_sum, post_slope_sum)
	
	vars <- colnames(posts_sum)

	# plot function 
	plot_func <- function(var){
		
		new_dat <- posts_sum %>% dplyr::select(var)
		
		base_plot <- 
			ggplot() +
			geom_density(data = new_dat, aes(x = new_dat[, 1])) 
		
		dens_df <- ggplot_build(base_plot)$data[[1]]

		p <- 
			ggplot() + 
			geom_area(data = subset(dens_df, x < 0),
  	                       aes(x=x,y=y),
  	                       fill = "#cccccc",
  	                       color = "white") +
			geom_area(data = subset(dens_df, x > 0),
  	                       aes(x=x,y=y),
  	                       fill = "#4c4c4c",
  	                       color = "white") +
			scale_x_continuous(
				name = "effect size",
				breaks = c(-1, 0, 1),
				labels = c(-1, 0, 1),
				limits = c(-1, 1.5)) +
			scale_y_continuous(
				sec.axis=sec_axis(~., breaks = NULL)) + 
			geom_segment(aes(x = -1, xend = 1.5, y = 0, yend = 0), color = "black") +
			theme_bw() +
			theme(
				panel.background = element_rect(fill = "transparent"),
    		plot.background = element_rect(fill = "transparent", color = NA),
				panel.grid   = element_blank(),
				legend.position = "none",
				axis.title = element_blank(),
				axis.text = element_blank(),
				axis.ticks = element_blank(),
				axis.line.y = element_line(),
				panel.border = element_blank())
	}
	
	vars <- colnames(posts_sum)
	plot_list <- map(vars, plot_func)	
	
	# make plot top
	
	plot_top <- 
		plot_list[[1]] +
			scale_x_continuous(
				limits = c(-1, 1.5),
				position = "top") +
			scale_y_continuous(
				sec.axis=sec_axis(~., breaks = NULL)) + 
			annotate("text", x = -1 , y = 1.6,
      	label = "(c)", size = 5) +
			theme(axis.line = element_line())
	
	# make bottom plot

	plot_bot <- 
		plot_list[[4]] +
			scale_x_continuous(
				name = "effect size",
				breaks = c(-1, 0, 1),
				labels = c(-1, 0, 1),
				limits = c(-1, 1.5)) +
			scale_y_continuous(
				sec.axis=sec_axis(~., breaks = NULL)) + 
			theme(
				axis.line = element_line(),
				axis.title.y = element_blank(),
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
				axis.text.x = element_text(size = 12),
				axis.title.x = element_text(size = 14))
	
	### plot together ####
	
	plot1 <- plot_top + 
		geom_vline(aes(xintercept = 0), color = "black", size = 0.5) +
		annotate("text", label = "68.3%",
						 x = 1.34, y = 0.1, size = 4) +
		annotate("text", label = "GSA\nintercept",
						 x = -0.8, y = 0.2, size = 4) 
	
	plot2 <- plot_list[[2]] + 
		geom_vline(aes(xintercept = 0), color = "black", size = 0.5) +
		annotate("text", label = "99.8%",
						 x = 1.34, y = 0.1, size = 4) +
		annotate("text", label = "caudal fin aspect ratio\nin intercept model",
						 x = -0.57, y = 0.23, size = 4) 
	
	plot3 <- plot_list[[3]] +	
		geom_vline(aes(xintercept = 0), color = "black", size = 0.5) +
		annotate("text", label = "79.6%",
						 x = 1.34, y = 0.1, size = 4) +
		annotate("text", label = "GSA slope",
						 x = -0.78, y = 0.12, size = 4) 
	
	plot4 <- plot_bot + 
		geom_vline(aes(xintercept = 0), color = "black", size = 0.5) +
		annotate("text", label = "99.9%",
						 x = 1.35, y = 0.1, size = 4) +
		annotate("text", label = "caudal fin aspect ratio\nin slope model",
						 x = -0.52, y = 0.24, size = 4)
	
	coef_plot <- 
		plot1 + plot_spacer() + 
		plot2 + plot_spacer() + 
		plot3 + plot_spacer() + 
		plot4 + 
		plot_layout(ncol = 1,  
		heights = c(1, -0.5, 1, -0.5, 1, -0.6, 1)) 
		

	#### plot all together ###
	
	plot5 <- CFAR_int + 
		       ggtitle("Intercept") + 
		       theme(plot.title = element_text(size = 12, hjust = 0.5,  face = "bold")) +
					 theme(plot.margin = unit(c(0.2, 0, 0.2, 0.2), "in"))
	
	plot6 <- CFAR_slope + 
		       ggtitle("Slope") + 
		       theme(plot.title = element_text(size = 12, hjust = 0.5,  face = "bold")) +
					 theme(plot.margin = unit(c(0.2, 0, 0.2, 0.06), "in"))

	Fig_CFAR <- plot5 + plot6  + coef_plot + plot_layout(ncol = 3, widths = c(1, 1, 1)) 
	
	ggsave(path = here("./output/plots/"), filename = "Fig_CFAR.png",
				 plot = Fig_CFAR,
	       height = 7, width = 14.5, units = "in")

