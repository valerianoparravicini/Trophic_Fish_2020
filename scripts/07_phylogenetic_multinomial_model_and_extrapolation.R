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

plan <- 
  drake_plan(
    
    ## agreements
    classif = read.csv("data/experts_classification.csv", sep=",", dec=".", row.names = 1, na.strings = "NA"),
    figure1 = make_plot1(classif),
    figure2 = make_plot2,
    
    ## cluster analysis
    # load diet data
    diets = read.csv("data/data_ISfull_grp_final_no_sp_no_wormy.csv", sep=",", dec=".", row.names=1, na = "NA"),
    # standardize diet data
    diets_p = get_dietp(diets),
    # run custer analysis
    modules = cluster_analysis(diets_p),
    # save results
    result_mod_predators = write.csv(modules$mod_predators, "results/mod_predators.csv", row.names = FALSE),
    result_mod_prey_items = write.csv(modules$mod_prey_items, "results/mod_prey_items.csv", row.names = FALSE),
    
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
  )

