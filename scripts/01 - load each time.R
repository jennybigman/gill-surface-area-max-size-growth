# 01-load all data

	# load libraries
	library(ape)
	library(beepr)
	library(broom)
	library(brms)
	library(car)
	library(cowplot)
	library(data.table)
	library(dplyr)
	library(forcats)
	library(forestmangr)
	library(geiger)
	library(ggeffects)
	library(ggplot2)
	library(here)
	library(lme4)
	library(lmodel2)
	library(loo)
	library(magick)
	library(MASS)
	library(patchwork)
	library(purrr)
	library(readr)
	library(reshape2)
	library(rfishbase)
	library(rstan)
	library(stringr)
	library(tidyr)
	library(tidyverse)
	library(quantreg)
	
	`%notin%` <- Negate(`%in%`)

	# set options for stan
	rstan_options(auto_write = TRUE)
	options(mc.cores = parallel::detectCores())

	# helpful functions
	len <- function(y) {
	  z <-length(unique(y))
	  return(z)
	}
	
		# load gill surface area & growth data species with raw data
	
	# all raw data (for later filtering by body size range)
	RawGSA <- read_csv(file = "./data/all_raw_data.csv")

	# for those species with at least 8 individuals and on phylogeny
	RawGSA8_phylo <- read_csv(file = here("./data/RawGSA8_phylo.csv"))

	# load phylogeny
	fish_elasmo_supertree <- read.tree(file = here("./data/fish_elasmo_supertree.tre"))
	len(fish_elasmo_supertree$tip.label)

	## which species are on the tree?
	Species <- lapply(strsplit(as.character(RawGSA8_phylo$Binomial), "\\ "), "[", 2)
	Genus <- lapply(strsplit(as.character(RawGSA8_phylo$Binomial), "\\ "), "[", 1)
	species_list <- tibble(Genus, Species)
	RawGSA8_phylo$phylo <- paste(species_list$Genus, species_list$Species,
	                                sep = "_")

len(RawGSA8_phylo$Binomial)

sp_drop <- setdiff(fish_elasmo_supertree$tip.label, RawGSA8$phylo)

tree_pruned <- drop.tip(fish_elasmo_supertree, sp_drop)

len(tree_pruned$tip.label)
phylo_list <- tree_pruned$tip.label

len(unique(RawGSA8$Binomial))
sp_df <- unique(RawGSA8$Binomial)

RawGSA8_phylo <- filter(RawGSA8, phylo %in% phylo_list)

raw_phylo <- sort(unique(RawGSA8_phylo$Binomial))

len(RawGSA8_phylo$Binomial)

setdiff(RawGSA8_phylo$phylo, tree_pruned$tip.label)
setdiff(tree_pruned$tip.label, RawGSA8_phylo$phylo)

setdiff(raw8, raw_phylo)
