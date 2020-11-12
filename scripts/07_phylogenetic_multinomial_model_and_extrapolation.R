library(caret)
library(fishualize)
library(tidyverse)
library(reshape)
library(lattice)
require(gridExtra) 
library(viridis)
library(drake)
library(parallel)
library(bipartite)
library(mcclust)
library(brms)
library(fishtree)
library(fishflux)
library(tidyverse)
library(parallel)
library(rfishbase)
library(drake)
library(MCMCglmm)
library(ape)
library(picante)
library(parallel)
library(rstan)
library(tidybayes)
library(ggnewscale)
library(tidytree)
library(ggtree)

Sys.setlocale("LC_MESSAGES", "en_US.utf8")


    # phylogenetic analysis
    diet_cat = load_and_clean(modules$mod_predators),
    tree = fishtree_complete_phylogeny(diet_cat$species),
    dietcat = modify_dietcat(diet_cat, tree),
    delta_statistic = delta_stat(dietcat),  #25
    fit_diet = run_dietmodel(dietcat, tree),
    r_phylo = extract_r_phylo(fit_diet),
    b = extract_b(fit_diet),
    extrap_test = test_extrapolate(b, r_phylo, dietcat, tree),
    extrapolation = extrapolate(b, r_phylo),
    delta_null = delta_nullm(dietcat, n = 200),
    adapt_extrapolation_cats(extrapolation)
  

