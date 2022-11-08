# FIGURE - k, Winf, GSA

	# transformations
	RawGSA8_phylo$Log_k <- log10(RawGSA8_phylo$mean_K)
	RawGSA8_phylo$Log_Winf <- log10(RawGSA8_phylo$mean_Winfinity) 
	RawGSA8_phylo$Winf_std <- ((RawGSA8_phylo$Log_Winf) - mean(RawGSA8_phylo$Log_Winf)) / sd(RawGSA8_phylo$Log_Winf)
	
	## extract intercepts, slopes, and GAI and put into dataset for plotting ##
	mean_dat <- distinct(RawGSA8_phylo, Binomial, .keep_all = TRUE) %>%
							dplyr::select(Binomial, mean_K, mean_Winfinity, Log_k, Log_Winf, Winf_std)

	# intercepts
	load(file = here("./output/k_Winf_INT.rds"))
	post_k_Winf_INT <- as.data.frame(k_Winf_INT)
	ints <- post_k_Winf_INT %>% dplyr::select(contains("beta_ints_std"))
	mean_dat$ints <- colMeans(ints)
	
	# slopes
	load(file = here("./output/k_Winf_slope.rds"))
	post_k_Winf_slope <- as.data.frame(k_Winf_slope)
	slopes <- post_k_Winf_slope %>% dplyr::select(contains("beta_slopes_std"))
	mean_dat$slopes <- colMeans(slopes)

	#### plot with no lines ####

	mid <- 0

	colfunc <- colorRampPalette(c( "slategray", "navy"))
	colfunc(3)

	# k, Winf, GSA intercept 
	KWINF_INT <- 
		ggplot(mean_dat) +
		geom_point(data = mean_dat,
			aes(x = mean_dat$Winf_std, y = mean_dat$Log_k,
			colour = ints), alpha = 0.6, size = 4) +
		scale_color_gradient(low = "steelblue1",
    	high = "darkblue") +
		annotate("text", x = -1.92 , y = 0.64,
    	label = "(a)", size = 4) +
		scale_x_continuous(
      name = (expression(italic(W['\U221E']))),
      breaks=c(-2, 0, 2),
      labels=c(-2, 0, 2),
      limits = c(-1.92, 2.23))  +
    scale_y_continuous(
      name = expression(italic("k")),
      breaks = c(-1.5, -0.5, 0.5),
      labels = c(-1.5, -0.5, 0.5),
      limits = c(-1.47, 0.64)) +
    theme_bw() + 
    theme(
  		legend.position = "none",
  		axis.title = element_text(size=12),
  		panel.grid.major = element_blank(),
  		panel.grid.minor = element_blank(),
  		panel.border = element_rect(linetype = 1, size= 1, fill = NA),
  		plot.margin = unit(c(0.2, 0, 0.2, 0.2), "in"))
	
	# slope
	KWINF_slope <- 
		ggplot(mean_dat) +
		geom_point(data = mean_dat,
			aes(x = mean_dat$Winf_std, y = mean_dat$Log_k,
			colour = slopes), alpha = 0.6, size = 4) +
		scale_color_gradient(
			low = "steelblue1",
      high = "darkblue") +
		annotate("text", x = -1.92 , y = 0.64,
      label = "(b)", size = 4) +
		scale_x_continuous(
      name = (expression(italic(W['\U221E']))),
      breaks=c(-2, 0, 2),
      labels=c(-2, 0, 2),
      limits = c(-1.92, 2.23))  +
    scale_y_continuous(
      breaks = c(-1.5, -0.5, 0.5),
      labels = c(-1.5, -0.5, 0.5),
      limits = c(-1.47, 0.64)) +
		theme_bw() + 
    theme(
  		legend.position = "none",
  		axis.text.y = element_blank(),
  		axis.ticks.y = element_blank(),
  		axis.title.x = element_text(size=12),
  		axis.title.y = element_blank(),
  		panel.grid.major = element_blank(),
  		panel.grid.minor = element_blank(),
  		panel.border = element_rect(linetype = 1, size= 1, fill = NA),
  		plot.margin = unit(c(0.2, 0, 0.2, -0.05), "in"))

	# make the coef plot by hand
	
	post_slope_sum <- post_k_Winf_slope %>%
		dplyr::select(b_GSA_slope, b_Winf_slope) 
		
	post_int_sum <- post_k_Winf_INT %>%
		dplyr::select(b_GSA_int, b_Winf_int) 
	
	posts_sum <- bind_cols(post_slope_sum, post_int_sum)
	
	vars <- colnames(posts_sum)

	## top plot 
	top_plot <- 
			ggplot() +
			geom_density(data = posts_sum, aes(x = b_GSA_int)) 
		
	top_dens_df <- ggplot_build(top_plot)$data[[1]]

	plot_top <- 
			ggplot() + 
			geom_area(data = subset(top_dens_df, x < 0),
  	                       aes(x=x,y=y),
  	                       fill = "#cccccc",
  	                       color = "white") +
			geom_area(data = subset(top_dens_df, x > 0),
  	                       aes(x=x,y=y),
  	                       fill = "#4c4c4c",
  	                       color = "white") +
			scale_x_continuous(
				limits = c(-0.7, 0.6),
				position = "top") +
			scale_y_continuous(
				sec.axis=sec_axis(~., breaks = NULL)) + 
			geom_segment(aes(x = -0.7, xend = 0.6, y = 0, yend = 0), color = "black") +
			annotate("text", x = -0.68 , y = 4.8,
      	label = "(c)", size = 4) +
			theme_bw() +
			theme(
				panel.background = element_rect(fill = "transparent"),
    		plot.background = element_rect(fill = "transparent", color = NA),
				panel.grid   = element_blank(),
				axis.line = element_line(),
				legend.position = "none",
				axis.title = element_blank(),
				axis.text = element_blank(),
				axis.ticks = element_blank(),
			  panel.border = element_blank())
	
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
				breaks = c(-0.5, 0, 0.5),
				labels = c(-0.5, 0, 0.5),
				limits = c(-0.7, 0.6)) +
			scale_y_continuous(
				sec.axis=sec_axis(~., breaks = NULL)) + 
			geom_segment(aes(x = -0.7, xend = 0.6, y = 0, yend = 0), color = "black") +
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
	
	# make last plot with axis 
	axis_var <- vars[4]

	bot_plot <- 
			ggplot() +
			geom_density(data = posts_sum, aes(x = b_Winf_slope)) 
		
	bot_dens_df <- ggplot_build(bot_plot)$data[[1]]

	plot_bot <- 
			ggplot() + 
			geom_area(data = subset(bot_dens_df, x < 0),
  	                       aes(x=x,y=y),
  	                       fill = "#cccccc",
  	                       color = "white") +
			geom_area(data = subset(bot_dens_df, x > 0),
  	                       aes(x=x,y=y),
  	                       fill = "#4c4c4c",
  	                       color = "white") +
			scale_x_continuous(
				name = "effect size",
				breaks = c(-0.5, 0, 0.5),
				labels = c(-0.5, 0, 0.5),
				limits = c(-0.7, 0.6)) +
			scale_y_continuous(
				sec.axis=sec_axis(~., breaks = NULL)) + 
			geom_segment(aes(x = -0.7, xend = 0.6, y = 0, yend = 0), color = "black") +
			theme_bw() +
			theme(
				panel.background = element_rect(fill = "transparent"),
    		plot.background = element_rect(fill = "transparent", color = NA),
				panel.grid   = element_blank(),
				axis.line = element_line(),
				legend.position = "none",
				axis.title.y = element_blank(),
				axis.text.y = element_blank(),
				axis.ticks.y = element_blank(),
			  panel.border = element_blank(),
				axis.text.x = element_text(size = 10),
				axis.title.x = element_text(size = 12))
	
	plot1 <- plot_top + 
		geom_vline(aes(xintercept = 0), color = "black", size = 0.5) +
		annotate("text", label = "13.0%",
						 x = 0.53, y = 0.5, size = 4) +
		annotate("text", label = "GSA\nintercept",
						 x = -0.6, y = 0.75, size = 4) +
		theme(plot.margin = unit(c(0.2, 0.2, 0.2, -0.5), "in"))

	plot2 <- plot_list[[2]] + 
		geom_vline(aes(xintercept = 0), color = "black", size = 0.5) +
		annotate("text", label = "0%",
						 x = 0.55, y = 0.6, size = 4) +
		annotate("text", label = "maximum\nsize in\nintercept\nmodel",
						 x = -0.58, y = 1.8, size = 4) +
		theme(plot.margin = unit(c(0.2, 0.2, 0.2, -0.5), "in"))

	
	plot3 <- plot_list[[1]] +	
		geom_vline(aes(xintercept = 0), color = "black", size = 0.5) +
		annotate("text", label = "25.0%",
						 x = 0.53, y = 0.5, size = 4) +
		annotate("text", label = "GSA slope",
						 x = -0.58, y = 0.5, size = 4) +
		theme(plot.margin = unit(c(0.2, 0.2, 0.2, -0.5), "in"))

	
	plot4 <- plot_bot + 
		geom_vline(aes(xintercept = 0), color = "black", size = 0.5) +
		annotate("text", label = "0%",
						 x = 0.55, y = 0.6, size = 4) +
		annotate("text", label = "maximum\nsize in\nslope\nmodel",
						 x = -0.58, y = 1.8, size = 4) +
		theme(plot.margin = unit(c(0, 0, 0, 0), "pt"))

	
	coef_plot <- 
		plot1 + plot_spacer() + 
		plot2 + plot_spacer() + 
		plot3 + plot_spacer() + 
		plot4 + 
		plot_layout(ncol = 1,  
		heights = c(1, -0.5, 1, -0.5, 1, -0.5, 1)) 
	

	plot5 <- KWINF_INT + 
		       ggtitle("Intercept") + 
		       theme(plot.title = element_text(size = 12, hjust = 0.5,  face = "bold")) 
	plot6 <- KWINF_slope + 
		       ggtitle("Slope") + 
		       theme(plot.title = element_text(size = 12, hjust = 0.5,  face = "bold")) 

	FigK_WINF_NOLINES <- plot5 + plot_spacer() +
											 plot6 + plot_spacer() +
											 coef_plot +
		plot_layout(ncol = 5, widths = c(1, -0.02, 1, -0.035, 1)) 
	
	ggsave(path = here("./output/plots/"), filename = "FigK_WINF_NOLINES_2Colors.png",
				 plot = FigK_WINF_NOLINES,
	       height = 7, width = 14.5, units = "in")


