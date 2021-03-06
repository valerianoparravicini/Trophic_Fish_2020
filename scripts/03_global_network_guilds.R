
# --- >> dataset load and standardization

diets <- read.csv("data/data_guts.csv", sep=",", dec=".", row.names=1, na = "NA") 


diets$fish_sp <- dplyr::recode(diets$fish_sp, 
                               "Apogon ellioti" = "Jaydia ellioti",               
                               "Neomyxus chaptalii" = "Neomyxus leuciscus",
                               "Zebrasoma veliferum" = "Zebrasoma velifer",
                               "Epinephelus macrosplos" = "Epinephelus macrospilos",
                               "Labracinus melanotaenia" = "Labracinus cyclophthalmus",
                               "Sargocentron spniferum" = "Sargocentron spiniferum",
                               "Parapercis cephalopunctata" = "Parapercis millepunctata",
                               "Tylosurus acus" = "Tylosurus acus acus",
                               "Archamia fucata" = "Taeniamia fucata",
                               "Chromis caerulea" = "Chromis ternatensis",
                               "Arothron hispdus" = "Arothron hispidus",
                               "Dasyatis americana" = "Hypanus americanus", # this is a stingray btw
                               "Apogon hyalosoma" ="Yarica hyalosoma"
)


diets_p <- reshape::cast(diets, fish_sp ~ grp, value = "quantity", fun.aggregate=sum)

diets_p[is.na(diets_p)] = 0

diets_samplesize <- sapply(as.character(unique(diets_p$fish_sp)), function(x) {sum(diets[diets$fish_sp == x,]$nb_guts)})

min_samplesize <- 5

diets_p <- diets_p[diets_p$fish_sp%in% names(diets_samplesize[diets_samplesize >=min_samplesize]),]
species <- diets_p[,1]
diets_p <- diets_p[,2:ncol(diets_p)]

sum_row <- rowSums(diets_p)

diets_p <- lapply(1:length(sum_row), function(r) {diets_p[r,]/sum_row[r]})
diets_p <- do.call(rbind, diets_p)

rownames(diets_p) <- species 

diets_p <- na.omit(diets_p)


# --- >> create 500 classifications according to modularity maximization

clusters <- parallel::mclapply(1:500, function(i) {
  
  cat(i, "\n")
  set.seed(i)
  
  bipartite::computeModules(t(diets_p))
  
}, mc.cores = 40)#eo mclapply clusters


# numper of modules 

n_mod <- sapply(clusters, function(x) {
  
  length(bipartite::listModuleInformation(x)[[2]])
  
})#eo sapply n_mod


# --- >> find module ID for predators

pred_modules <- parallel::mclapply(clusters, function(x) {
  
  mat_memb <- x@modules[2:nrow(x@modules), 2:ncol(x@modules)]
  
  sapply(1:ncol(mat_memb), function (i) {
    
    which(mat_memb [,i] == max(mat_memb [,i]))
    
  })#eo sapply 
  
}, mc.cores = 40)#eo mclapply pred_modules


# --- >> define the medoid solutions according to VII statistics of Meila (2007)


combinations <- expand.grid(1:500, 1:500)


vii <- parallel::mclapply(1:nrow(combinations), function (x) {
  
  mcclust::vi.dist(pred_modules[[combinations[x,1]]], pred_modules[[combinations[x,2]]], parts = FALSE, base = 2)
  
}, mc.cores = 40)#eo mclapply vii



vii_per_cluster <- data.frame(x = combinations[,1], y = combinations[,2], vii = do.call(rbind,vii))

cast_vii <- reshape::cast(vii_per_cluster, x ~ y, value = "vii")

medoid <- which(colMeans(cast_vii) == min(colMeans(cast_vii)))


medoid_guilds = clusters[[medoid]]

resList <- bipartite::listModuleInformation(medoid_guilds)



# --- >> extract and save modules information


# print results
bipartite::printoutModuleInformation(medoid_guilds)


getModules_predators <- lapply(1:length(resList[[2]]), function(x) {
  
  cbind(resList[[2]][[x]][[2]], rep(x, length(resList[[2]][[x]][[2]])))
  
  
})

getModules_prey <- lapply(1:length(resList[[2]]), function(x) {
  
  cbind(resList[[2]][[x]][[1]], rep(x, length(resList[[2]][[x]][[1]])))
  
  
})


mod_predators <- do.call(rbind, getModules_predators)
mod_prey_items <- do.call(rbind, getModules_prey)

colnames(mod_predators) <- c("species", "cluster")
colnames(mod_prey_items) <- c("prey_item", "cluster")

# Standardize category numbers 
# motivation: numbers of clusters are random each time analysis is run
library(tidyverse)

ref <- # reference dataframe with species
  data.frame(
    species = c(
      "Cephalopholis argus",
      "Parupeneus multifasciatus",
      "Acanthurus lineatus",
      "Chaetodon ornatissimus",
      "Abudefduf sexfasciatus",
      "Balistapus undulatus",
      "Ostracion meleagris",
      "Stethojulis bandanensis"),
    ref = c(4, 7, 2, 3, 8, 6, 1, 5)
  ) %>% 
  dplyr::mutate(species = as.character(species))

ref <- mod_predators %>%
  dplyr::inner_join(ref) %>%
  dplyr::select(-species)
  
mod_predators <- as.data.frame(mod_predators) %>%
  dplyr::left_join(ref) %>%
  dplyr::select(-cluster) %>%
  dplyr::rename(cluster = ref)

mod_prey_items <- as.data.frame(mod_prey_items) %>%
  dplyr::left_join(ref) %>%
  dplyr::select(-cluster) %>%
  dplyr::rename(cluster = ref)


# Save results
save(mod_predators, mod_prey_items, file = "output/results/trophic_guilds_medoid.RData")

write.csv(mod_predators, "output/results/mod_predators.csv", row.names = FALSE)
write.csv(mod_prey_items, "output/results/mod_prey_items.csv", row.names = FALSE)

