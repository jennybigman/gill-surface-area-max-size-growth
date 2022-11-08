	# 08 -- estimating the relationship of growth performance and gill area index only for the raw dataset 
	
	# here, we are working with the raw dataset, but to estimate gill area index, we need to take a mean of 
	# gill surface area and body mass 

	mean_dat <- RawGSA8_phylo %>% 
						  group_by(Binomial, mean_gperf_simple) %>%
						  summarise(meanGSAcm2 = mean(GSAcm2),
						  					meanMassG = mean(MassG))
	
	### estimate gill area index using d = 0.8
	
	GAI_formula <- function(G, W){
	  GAI = G / (W^0.8)
	  return(GAI)
	}
	
	mean_dat <- mean_dat %>%
	            rowwise() %>%
	            mutate(GAI_setd = GAI_formula(G = meanGSAcm2, W = meanMassG))
	

	## estimate gill area index using Pauly's relationship between d & max size
	
	GSA_growth_dat_raw <- filter(GSA_growth_dat, Binomial %in% mean_dat$Binomial)
	
	
	max_size_dat <- distinct(GSA_growth_dat_raw, Binomial, .keep_all = TRUE) %>%
									dplyr::select(Binomial, mean_Winfinity)
	
	mean_dat <- merge(mean_dat, max_size_dat, by = "Binomial")
	
	mean_dat$LogWinfinity <- log10(mean_dat$mean_Winfinity)


	# Pauly's relationship (Pauly 1981)
	
	d_formula <- function(x){
	  d = 0.6742 + 0.03574 * x
	  return(d)
	}
	
	GAI_formula_est <- function(G, W, x){
	  GAI = G / (W^x)
	  return(GAI)
	}
	
	mean_dat <- mean_dat %>%
	            rowwise() %>%
	            mutate(est_d = d_formula(x = LogWinfinity),
	            			 GAI_estd = GAI_formula_est(G = meanGSAcm2, W = meanMassG, x = est_d),
	            			 LogGAI_estd = log10(GAI_estd),
	            			 LogGAI_setd = log10(GAI_setd))
	
	mean_setd <- mean(mean_dat$LogGAI_setd)
	sd_setd <- sd(mean_dat$LogGAI_setd)
	
	mean_estd <- mean(mean_dat$LogGAI_estd)
	sd_estd <- sd(mean_dat$LogGAI_estd)
	            			 
	mean_dat <- mean_dat %>%
	            mutate(LogGAI_setd_std = ((LogGAI_setd - mean_setd) / sd_setd),
	            	     LogGAI_estd_std = ((LogGAI_estd - mean_estd) / sd_estd))
	

#### relationship between GAI_setd & growth performance ##

	#get_prior(mean_gperf_simple ~ LogGAI_setd_std,
	#          data = mean_dat,
	#          family = gaussian())
	# 
	#GP_GAI_setd_raw <-  brm(data = mean_dat, family = gaussian,
	#      						      mean_gperf_simple ~ LogGAI_setd_std,
	#      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	#      						           	    prior(student_t(3, 5, 10), class = b),
	#      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	#      						           	chains = 4,
	#      						           	iter = 5000, warmup = 1000)
	#
	#save(GP_GAI_setd_raw, file = here("./output/GP_GAI_setd_raw.rds"))
	load(file = here("./output/GP_GAI_setd_raw.rds"))

	post_GP_GAI_setd_raw_INT <- as.data.frame(GP_GAI_setd_raw)
	post_GP_GAI_setd_raw_INT_mean <- stack(apply(post_GP_GAI_setd_raw_INT, 2, mean)) %>% rename(mean = values)
	post_GP_GAI_setd_raw_INT_LCIs <- stack(apply(post_GP_GAI_setd_raw_INT, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_GAI_setd_raw_INT_HCIs <- stack(apply(post_GP_GAI_setd_raw_INT, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_GP_GAI_setd_raw_INT_sum <- merge(post_GP_GAI_setd_raw_INT_mean, post_GP_GAI_setd_raw_INT_LCIs, by = "ind") %>%
														      merge(post_GP_GAI_setd_raw_INT_HCIs) 
	post_GP_GAI_setd_raw_INT_sum <- round_df(post_GP_GAI_setd_raw_INT_sum, 2) %>% round_df(2)

	head(post_GP_GAI_setd_raw_INT_sum)
	
	# proportion > 0
	GP_GAI_setd_raw_INT_sum_prob0 <- (length(which(post_GP_GAI_setd_raw_INT$b_LogGAI_setd_std >0)))/length(post_GP_GAI_setd_raw_INT$b_LogGAI_setd_std)*100 
	GP_GAI_setd_raw_INT_sum_prob0 <- round(GP_GAI_setd_raw_INT_sum_prob0, 1) 
	
	# slope = 0.44 (0 to 0.88)

	loo(GP_GAI_setd_raw)

## relationship between GAI_estd & growth performance 

	#get_prior(mean_gperf_simple ~ LogGAI_estd_std,
	#          data = mean_dat,
	#          family = gaussian())
	# 
	#GP_GAI_estd_raw <-  brm(data = mean_dat, family = gaussian,
	#      						      mean_gperf_simple ~ LogGAI_estd_std, 
	#      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
	#      						           	    prior(student_t(3, 5, 10), class = b),
	#      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
	#      						           	chains = 4,
	#      						           	iter = 5000, warmup = 1000)
	#
	#
	#save(GP_GAI_estd_raw, file = here("./output/GP_GAI_estd_raw.rds"))
	load(here("./output/GP_GAI_estd_raw.rds"))
	
	summary(GP_GAI_estd_raw)
	# slope = 0.06 (-0.43 to 0.54)

	loo(GP_GAI_estd_raw)
	
	post_GP_GAI_estd_raw <- as.data.frame(GP_GAI_estd_raw)
	post_GP_GAI_estd_raw_mean <- stack(apply(post_GP_GAI_estd_raw, 2, mean)) %>% rename(mean = values)
	post_GP_GAI_estd_raw_LCIs <- stack(apply(post_GP_GAI_estd_raw, 2, quantile, prob = (0.025))) %>% rename(LCI = values)
	post_GP_GAI_estd_raw_HCIs <- stack(apply(post_GP_GAI_estd_raw, 2, quantile, prob = (0.975))) %>% rename(HCI = values)
	post_GP_GAI_estd_raw_sum <- merge(post_GP_GAI_estd_raw_mean, post_GP_GAI_estd_raw_LCIs, by = "ind") %>%
														  merge(post_GP_GAI_estd_raw_HCIs) 
	post_GP_GAI_estd_raw_sum <- round_df(post_GP_GAI_estd_raw_sum, 2) %>% round_df(2)

	head(post_GP_GAI_estd_raw_sum)
	
	# proportion > 0
	GP_GAI_estd_raw_sum_prob0 <- (length(which(post_GP_GAI_estd_raw$b_LogGAI_estd_std >0)))/length(post_GP_GAI_estd_raw$b_LogGAI_estd_std)*100 
	GP_GAI_estd_raw_sum_prob0 <- round(GP_GAI_estd_raw_sum_prob0, 1) 
	


#### predictions and plots ####

## for GAI estimated with set d 

# create new dataframe for fitted values
LogGAI_setd_std <- mean_dat$LogGAI_setd_std
newdat_setd_std <- data.frame(LogGAI_setd_std)

# extract fitted values
LogGAI_setd_fitted <-
  fitted(GP_GAI_setd_raw, 
         newdata = newdat_setd_std) %>%
  as_tibble() %>%
  bind_cols(newdat_setd_std) %>%
	rename(lowerCI = "Q2.5",
				 higherCI = "Q97.5")


GP_GAI_setd_plot <- 		ggplot(LogGAI_setd_fitted) +
	                      geom_ribbon(data = LogGAI_setd_fitted, 
                  								  aes(x = LogGAI_setd_std, 
                  								  		ymin = lowerCI,
                  								  		ymax = higherCI), 
                  									fill = "grey", alpha = 0.3) +
          						  geom_line(data = LogGAI_setd_fitted,
          						           aes(x = LogGAI_setd_std, y = Estimate),
          						               color = "darkgrey", size = 1.25, alpha = 0.6) +
          						  geom_point(data = mean_dat,
          						            aes(x = LogGAI_setd_std,
          						                y = mean_gperf_simple),
          						                alpha = 0.6, size = 4) +
											  ylab(expression(paste("Growth performance", ' '(log(k * W[inf]), g/yr)))) +
											  xlab(expression(paste("Gill area index", ' '(G/W^{d}, cm^{2}/g)))) +
											  theme(legend.position = "none") +
          						  white_theme_noleg()

ggsave(GP_GAI_setd_plot, 
			 file = here("./output/plots/GP_GAI_setd_plot.jpg"),
									 height = 7.5, width = 10, units = "in")

## for GAI estimated with est d 

# create new dataframe for fitted values
LogGAI_estd_std <- mean_dat$LogGAI_estd_std
newdat_estd_std <- data.frame(LogGAI_estd_std)

# extract fitted values
LogGAI_estd_fitted <-
  fitted(GP_GAI_estd_raw, 
         newdata = newdat_estd_std) %>%
  as_tibble() %>%
  bind_cols(newdat_estd_std) %>%
	rename(lowerCI = "Q2.5",
				 higherCI = "Q97.5")


GP_GAI_estd_plot <- 		ggplot(LogGAI_estd_fitted) +
	                      geom_ribbon(data = LogGAI_estd_fitted, 
                  								  aes(x = LogGAI_estd_std, 
                  								  		ymin = lowerCI,
                  								  		ymax = higherCI), 
                  									fill = "grey", alpha = 0.3) +
          						  geom_line(data = LogGAI_estd_fitted,
          						           aes(x = LogGAI_estd_std, y = Estimate),
          						               color = "darkgrey", size = 1.25, alpha = 0.6) +
          						  geom_point(data = mean_dat,
          						            aes(x = LogGAI_estd_std,
          						                y = mean_gperf_simple),
          						                alpha = 0.6, size = 4) +
											  ylab(expression(paste("Growth performance", ' '(log(k * W[inf]), g/yr)))) +
											  xlab(expression(paste("Gill area index", ' '(G/W^{d}, cm^{2}/g)))) +
											  theme(legend.position = "none") +
          						  white_theme_noleg()

ggsave(GP_GAI_estd_plot, 
			 file = here("./output/plots/GP_GAI_estd_plot.jpg"),
									 height = 7.5, width = 10, units = "in")


### with phylogeny
#
## setd 
#
#get_prior(mean_gperf_simple ~ LogGAI_setd + (1|gr(phylo, cov = A)),
#          data = mean_dat,
#					data2 = list(A = A),
#          family = gaussian())
# 
#GP_GAI_setd_raw_phylo <-  brm(data = mean_dat, family = gaussian,
#      						      mean_gperf_simple ~ LogGAI_setd +
#												(1|gr(phylo, cov = A)), 
#												data2 = list(A = A),
#      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
#      						           	    prior(student_t(3, 0, 5), class = b),
#      						      					prior(student_t(3, 0, 2.5), class = sd),
#      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
#      						           	chains = 4,
#												      control = list("adapt_delta" = 0.9999,
#              		  							           "max_treedepth" = 20),
#      						           	iter = 15000, warmup = 2500)
#
#save(GP_GAI_setd_raw_phylo, file = here("./output/GP_GAI_setd_raw_phylo.rds"))
#
## slope = 

#loo(GP_GAI_setd_raw_phylo)

## relationship between GAI_estd & growth performance 

#get_prior(mean_gperf_simple ~ LogGAI_estd + (1|gr(phylo, cov = A)),
#          data = mean_dat,
#					data2 = list(A = A),
#          family = gaussian())
# 
#GP_GAI_estd_raw_phylo <-  brm(data = mean_dat, family = gaussian,
#      						      mean_gperf_simple ~ LogGAI_estd + 
#												(1|gr(phylo, cov = A)), 
#												data2 = list(A = A),
#      						      prior = c(prior(student_t(3, 2.9, 2.5), class = Intercept),
#      						           	    prior(student_t(3, 5, 10), class = b),
#      						      					prior(student_t(3, 0, 2.5), class = sd),
#      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
#      						           	chains = 4,
#      						           	iter = 5000, warmup = 1000)
#
## slope = 0.08 (-0.88 to 1.05)

#save(GP_GAI_estd_raw_phylo, file = here("./output/GP_GAI_estd_raw_phylo.rds"))

#loo(GP_GAI_estd_raw_phylo)



### plot all three together (the plot from 5d)

# GAI estimated empirically

load(file = here("./output/MLM_GAI_prog.rds"))

load(file = here("./output/MLM_GAI_CORRECT.rds"))

post_GAI1 <- as.data.frame(MLM_GAI)

mod_GAI <- post_GAI1 %>% dplyr::select(., contains("LogGAI"))

species_list <- sort(unique(RawGSA8_phylo$Binomial))

# add a "_" between genus and species
Species <- lapply(strsplit(as.character(species_list), "\\ "), "[", 2)
Genus <- lapply(strsplit(as.character(species_list), "\\ "), "[", 1)
species_list <- tibble(species_list, Genus, Species)
species_list$Species2 <- paste(species_list$Genus, species_list$Species,
                                sep = "_")
species_list <- species_list$Species2

names(mod_GAI) <- species_list


GAI_mean <- stack(apply(mod_GAI, 2, mean)) %>%
										  rename(species = ind,
										 			   GAI = values)


# get growth performance and add into dataset for plotting
GP_species <- distinct(RawGSA8_phylo, mean_gperf_simple, .keep_all = TRUE) %>%
							dplyr::select(Binomial, mean_gperf_simple) %>% ungroup()

GP_species$Species <- lapply(strsplit(as.character(GP_species$Binomial), "\\ "), "[", 2)
GP_species$Genus <- lapply(strsplit(as.character(GP_species$Binomial), "\\ "), "[", 1)
GP_species$species <- paste(GP_species$Genus, GP_species$Species,
                                sep = "_")

GP_species <- GP_species %>% dplyr::select(species, mean_gperf_simple)

GAI_mean <- merge(GAI_mean, GP_species)

ggplot(data = GAI_mean, 
											aes(x = GAI, y = mean_gperf_simple)) +
							 geom_point()

# fitted line
GAI_obs <- sort(GAI_mean$GAI)

fitted_matrix_GAI <- t(sapply(1:nrow( as.data.frame(MLM_GAI_prog)), function(i) {
                             post_GAI1[i, 'aGP_GAI'] +
                             post_GAI1[i, 'bGP_GAI'] * GAI_obs 
	}))

mean_fitted_matrix_GAI <- colMeans(fitted_matrix_GAI)

BCI_GAI <- t(apply(fitted_matrix_GAI, 2, quantile,
                    prob = c(0.025, 0.975))) %>%
            as.data.frame() %>%
            rename(lower_BCI = "2.5%",
                   upper_BCI = "97.5%") %>%
            mutate(GAI_obs = GAI_obs,
                   fitted_GP = mean_fitted_matrix_GAI)

GP_GAI_mlm_plot <- ggplot(BCI_GAI) +
	                      geom_ribbon(data = BCI_GAI, 
                  								  aes(x = GAI_obs, 
                  								  		ymin = lower_BCI,
                  								  		ymax = upper_BCI), 
                  									fill = "grey", alpha = 0.3) +
          						  geom_line(data = BCI_GAI,
          						           aes(x = GAI_obs, y = fitted_GP),
          						               color = "darkgrey", size = 1.25, alpha = 0.6) +
          						  geom_point(data = GAI_mean,
          						            aes(x = GAI,
          						                y = mean_gperf_simple),
          						                alpha = 0.6, size = 4) +
												xlim(-3.47, 2.18) +
												scale_y_continuous(breaks=c(1, 3, 5)) +
											  theme(legend.position = "none") +
												annotate("text", x = 2.11 , y = 5.4,
                         label = "(a)",
                         size = 3) +
												annotate("text", x = -2.73 , y = 5.4,
                         label = "d empirically estimated",
                         size = 3) +
												annotate("text", x = -2.55, y = 5.0,
                         label = "slope = 0.13 (-0.40 to 0.65)",
                         size = 3) +
											  theme_bw() +
											  theme(
											    axis.text.y =element_text(size= 12, colour = "black"),
											    axis.title= element_blank(),
											    axis.text.x = element_blank(),
											    axis.ticks.x = element_blank(),
											    panel.grid.major = element_blank(),
											    panel.grid.minor = element_blank(),
											    panel.border = element_rect(linetype = 1, size= 1, fill = NA))
								
ggplot_build(GP_GAI_estd_plot)$layout$panel_scales_y[[1]]$range$range
ggplot_build(GP_GAI_estd_plot)$layout$panel_scales_x[[1]]$range$range

# GAI estimated with set d

GP_GAI_setd_plot <- 		ggplot(LogGAI_setd_fitted) +
	                      geom_ribbon(data = LogGAI_setd_fitted, 
                  								  aes(x = LogGAI_setd_std, 
                  								  		ymin = lowerCI,
                  								  		ymax = higherCI), 
                  									fill = "grey", alpha = 0.3) +
          						  geom_line(data = LogGAI_setd_fitted,
          						           aes(x = LogGAI_setd_std, y = Estimate),
          						               color = "darkgrey", size = 1.25, alpha = 0.6) +
          						  geom_point(data = mean_dat,
          						            aes(x = LogGAI_setd_std,
          						                y = mean_gperf_simple),
          						                alpha = 0.6, size = 4) +
												annotate("text", x = 2.11 , y = 5.4,
                         label = "(b)",
                         size = 3) +
												annotate("text", x = -3.35 , y = 5.4,
                         label = "d = 0.8",
                         size = 3) +
												annotate("text", x = -2.55, y = 5.0,
                         label = "slope = 0.44 (0.00 to 0.88)",
                         size = 3) +
											  theme(legend.position = "none") +
												xlim(-3.47, 2.18) +
												scale_y_continuous(breaks=c(1, 3, 5)) +
											  theme_bw() +
											  theme(
											    axis.text.y=element_text(size= 12, colour = "black"),
											    axis.title= element_blank(),
											    axis.text.x = element_blank(),
											    panel.grid.major = element_blank(),
											    panel.grid.minor = element_blank(),
											    axis.ticks.x = element_blank(),
											    panel.border = element_rect(linetype = 1, size= 1, fill = NA))

# GAI estimated with est d from Pauly's relationship

GP_GAI_estd_plot <- 		ggplot(LogGAI_estd_fitted) +
	                      geom_ribbon(data = LogGAI_estd_fitted, 
                  								  aes(x = LogGAI_estd_std, 
                  								  		ymin = lowerCI,
                  								  		ymax = higherCI), 
                  									fill = "grey", alpha = 0.3) +
          						  geom_line(data = LogGAI_estd_fitted,
          						           aes(x = LogGAI_estd_std, y = Estimate),
          						               color = "darkgrey", size = 1.25, alpha = 0.6) +
          						  geom_point(data = mean_dat,
          						            aes(x = LogGAI_estd_std,
          						                y = mean_gperf_simple),
          						                alpha = 0.6, size = 4) +
												annotate("text", x = 2.11 , y = 5.4,
                         label = "(c)",
                         size = 3) +
												annotate("text", x = -3.2 , y = 5.4,
                         label = "d predicted",
                         size = 3) +
												annotate("text", x = -2.55, y = 5.0,
                         label = "slope = 0.06 (-0.43 to 0.54)",
                         size = 3) +
											  xlab(expression(paste("Gill area index", ' '(G/W^{d}, cm^{2}/g)))) +
												theme(legend.position = "none") +
												xlim(-3.47, 2.18) +
												scale_y_continuous(breaks=c(1, 3, 5)) +
											  theme_bw() +
											  theme(
											    axis.text=element_text(size= 14, colour = "black"),
											    axis.title.y= element_blank(),
											    panel.grid.major = element_blank(),
											    panel.grid.minor = element_blank(),
											    panel.border = element_rect(linetype = 1, size= 1, fill = NA))
# plot together

plot1 <- GP_GAI_mlm_plot + theme(plot.margin = unit(c(0.2, 0.2, 0, 0.2), "in"))
plot2 <- GP_GAI_setd_plot + theme(plot.margin = unit(c(-0.05, 0.2, 0, 0.2), "in"))
plot3 <- GP_GAI_estd_plot + theme(plot.margin = unit(c(-0.05, 0.2, 0.2, 0.2), "in"))

Fig_GAI <- (plot1 / plot2 / plot3) + plot_layout(heights  = c(1, 1, 1.05)) 

Fig_GAI <- add_global_label((Fig_GAI), Ylab = "Growth performance [log(k * Winf)]", size = 4, Ygap = 0.04)

ggsave(Fig_GAI, file = here("./output/plots/Comparison GP and GAI Pauly 3C.jpg"),
			 height = 7, width = 5, units = "in" )













#### without aquaculture species ####

## set d ## 

#estimate gill area index using d = 0.8

GAI_formula <- function(G, W){
  GAI = G / (W^0.8)
  return(GAI)
}

mean_dat_no_ac <- mean_dat_no_ac %>%
                   rowwise() %>%
                   mutate(GAI_setd = GAI_formula(G = meanGSAcm2, W = meanMassG))

## estimate gill area index using Pauly's relationship between d & max size

GSA_growth_dat_raw_no_ac <- filter(GSA_growth_dat_no_ac, Binomial %in% mean_dat_no_ac$Binomial)


max_size_dat_no_ac <- distinct(GSA_growth_dat_raw_no_ac, Binomial, .keep_all = TRUE) %>%
							      	select(Binomial, mean_Winfinity)

mean_dat_no_ac <- merge(mean_dat_no_ac, max_size_dat_no_ac, by = "Binomial")

mean_dat_no_ac$LogWinfinity <- log10(mean_dat_no_ac2$mean_Winfinity)

## est d ##

# d predicted from Pauly's relationship (Pauly 1981)

d_formula <- function(x){
  d = 0.6742 + 0.03574 * x
  return(d)
}

GAI_formula_est <- function(G, W, x){
  GAI = G / (W^x)
  return(GAI)
}

mean_dat_no_ac <- mean_dat_no_ac %>%
            rowwise() %>%
            mutate(est_d = d_formula(x = LogWinfinity),
            			 GAI_estd = GAI_formula_est(G = meanGSAcm2, W = meanMassG, x = est_d),
            			 LogGAI_estd = log10(GAI_estd),
            			 LogGAI_setd = log10(GAI_setd))

mean_setd_no_ac <- mean(mean_dat_no_ac$LogGAI_setd)
sd_setd_no_ac <- sd(mean_dat_no_ac$LogGAI_setd)

mean_estd_no_ac <- mean(mean_dat_no_ac$LogGAI_estd)
sd_estd_no_ac <- sd(mean_dat_no_ac$LogGAI_estd)
            			 
mean_dat_no_ac <- mean_dat_no_ac %>%
            mutate(LogGAI_setd_std = ((LogGAI_setd - mean_setd_no_ac) / sd_setd_no_ac),
            	     LogGAI_estd_std = ((LogGAI_estd - mean_estd_no_ac) / sd_estd_no_ac))


#### relationship between GAI_setd & growth performance ##

get_prior(mean_gperf_simple ~ LogGAI_setd_std,
          data = mean_dat_no_ac,
          family = gaussian())
 
GP_GAI_setd_raw_no_ac <-  brm(data = mean_dat_no_ac, family = gaussian,
      						      mean_gperf_simple ~ LogGAI_setd_std,
      						      prior = c(prior(student_t(3, 3, 2.5), class = Intercept),
      						           	    prior(student_t(3, 5, 10), class = b),
      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
      						           	chains = 4,
      						           	iter = 5000, warmup = 1000)

save(GP_GAI_setd_raw_no_ac, file = here("./output/GP_GAI_setd_raw_no_ac.rds"))
load(file = here("./output/GP_GAI_setd_raw_no_ac.rds"))

summary(GP_GAI_setd_raw_no_ac)

post_GAI_setd_ac <- data.frame(GP_GAI_setd_raw_no_ac)

post_GAI_setd_ac_sum <- stack(apply(post_GAI_setd_ac, 2, mean)) %>% round_df(2)

# slope = 0.51 (0.04 to 0.99)

loo(GP_GAI_setd_raw_no_ac)

## relationship between GAI_estd & growth performance 

get_prior(mean_gperf_simple ~ LogGAI_estd_std,
          data = mean_dat_no_ac,
          family = gaussian())
 
GP_GAI_estd_raw_no_ac <-  brm(data = mean_dat_no_ac, family = gaussian,
      						      mean_gperf_simple ~ LogGAI_estd_std, 
      						      prior = c(prior(student_t(3, 3, 2.5), class = Intercept),
      						           	    prior(student_t(3, 5, 10), class = b),
      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
      						           	chains = 4,
      						           	iter = 5000, warmup = 1000)


save(GP_GAI_estd_raw_no_ac, file = here("./output/GP_GAI_estd_raw_no_ac.rds"))
load(here("./output/GP_GAI_estd_raw_no_ac.rds"))

summary(GP_GAI_estd_raw_no_ac)
# slope = 0.10 (-0.41 to 0.62)

loo(GP_GAI_estd_raw_no_ac)


#### no air-breathers ####

## set d ##

# estimate gill area index using d = 0.8

GAI_formula <- function(G, W){
  GAI = G / (W^0.8)
  return(GAI)
}

mean_dat_no_abr <- mean_dat_no_abr %>%
                    rowwise() %>%
                    mutate(GAI_setd = GAI_formula(G = meanGSAcm2, W = meanMassG))


## estimate gill area index using Pauly's relationship between d & max size

GSA_growth_dat_raw_no_abr <- filter(GSA_growth_dat_no_abr, Binomial %in% 
																			mean_dat_no_abr$Binomial)


max_size_dat_no_abr <- distinct(GSA_growth_dat_raw_no_abr, Binomial, .keep_all = TRUE) %>%
								       select(Binomial, mean_Winfinity)

mean_dat_no_abr <- merge(mean_dat_no_abr, max_size_dat_no_abr, by = "Binomial")

mean_dat_no_abr$LogWinfinity <- log10(mean_dat_no_abr$mean_Winfinity)

## est d ##
 
## est d from Pauly's relationship (Pauly 1981)

d_formula <- function(x){
  d = 0.6742 + 0.03574 * x
  return(d)
}

GAI_formula_est <- function(G, W, x){
  GAI = G / (W^x)
  return(GAI)
}

mean_dat_no_abr <- mean_dat_no_abr %>%
            rowwise() %>%
            mutate(est_d = d_formula(x = LogWinfinity),
            			 GAI_estd = GAI_formula_est(G = meanGSAcm2, W = meanMassG, x = est_d),
            			 LogGAI_estd = log10(GAI_estd),
            			 LogGAI_setd = log10(GAI_setd))

mean_setd_no_abr <- mean(mean_dat_no_abr$LogGAI_setd)
sd_setd_no_abr <- sd(mean_dat_no_abr$LogGAI_setd)

mean_estd_no_abr <- mean(mean_dat_no_abr$LogGAI_estd)
sd_estd_no_abr <- sd(mean_dat_no_abr$LogGAI_estd)
            			 
mean_dat_no_abr <- mean_dat_no_abr %>%
            mutate(LogGAI_setd_std = ((LogGAI_setd - mean_setd_no_abr) / sd_setd_no_abr),
            	     LogGAI_estd_std = ((LogGAI_estd - mean_estd_no_abr) / sd_estd_no_abr))


#### relationship between GAI_setd & growth performance ##

get_prior(mean_gperf_simple ~ LogGAI_setd_std,
          data = mean_dat_no_abr,
          family = gaussian())
 
GP_GAI_setd_raw_no_abr <-  brm(data = mean_dat_no_abr, family = gaussian,
      						      mean_gperf_simple ~ LogGAI_setd_std,
      						      prior = c(prior(student_t(3, 3.5, 2.5), class = Intercept),
      						           	    prior(student_t(3, 5, 10), class = b),
      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
      						           	chains = 4,
      						           	iter = 5000, warmup = 1000)

save(GP_GAI_setd_raw_no_abr, file = here("./output/GP_GAI_setd_raw_no_abr.rds"))
load(file = here("./output/GP_GAI_setd_raw_no_abr.rds"))

summary(GP_GAI_setd_raw_no_abr)

post_GAI_setd_abr <- data.frame(GP_GAI_setd_raw_no_abr)

post_GAI_setd_abr_sum <- stack(apply(post_GAI_setd_abr, 2, mean)) %>% round_df(2)

loo(GP_GAI_setd_raw_no_abr)

## relationship between GAI_estd & growth performance 

get_prior(mean_gperf_simple ~ LogGAI_estd_std,
          data = mean_dat_no_abr,
          family = gaussian())
 
GP_GAI_estd_raw_no_abr <-  brm(data = mean_dat_no_abr, family = gaussian,
      						      mean_gperf_simple ~ LogGAI_estd_std, 
      						      prior = c(prior(student_t(3, 3.5, 2.5), class = Intercept),
      						           	    prior(student_t(3, 5, 10), class = b),
      						           	    prior(student_t(3, 0, 2.5),  class = sigma)),
      						           	chains = 4,
      						           	iter = 5000, warmup = 1000)


save(GP_GAI_estd_raw_no_abr, file = here("./output/GP_GAI_estd_raw_no_abr.rds"))
load(here("./output/GP_GAI_estd_raw_no_abr.rds"))

summary(GP_GAI_estd_raw_no_abr)

loo(GP_GAI_estd_raw_no_abr)