# ---------------------------------------------------------------------->> grp6 - multi-site


diets <- read.csv("data/data_ISfull_grp_final_no_sp_no_wormy.csv", sep=",", dec=".", row.names=1, na = "NA") 

#check_taxonomy <- rfishbase::validate_names(unique(diets$fish_sp))
#unique(subset(diets$fish_sp, !(diets$fish_sp %in% check_taxonomy)))

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


diets$sp_reg <- paste(diets$fish_sp, diets$site_code, sep="-")

diets_p <- reshape::cast(diets, sp_reg ~ grp6, value = "quantity", fun.aggregate=sum)
diets_p[is.na(diets_p)] = 0

diets_samplesize <- sapply(as.character(unique(diets_p$sp_reg)), function(x) {sum(diets[diets$sp_reg == x,]$nb_guts)})

min_samplesize <- 10


diets_p <- diets_p[diets_p$sp_reg %in% names(diets_samplesize[diets_samplesize >=min_samplesize]),]
species <- diets_p[,1]
diets_p <- diets_p[,2:ncol(diets_p)]

sum_row <- rowSums(diets_p)

diets_p <- lapply(1:length(sum_row), function(r) {diets_p[r,]/sum_row[r]})
diets_p <- do.call(rbind, diets_p)

rownames(diets_p) <- species 

diets_p <- na.omit(diets_p)


fish_sp_retained <- data.frame(do.call(rbind,strsplit(rownames(diets_p), "[-]")))

#number of locations

n_locations <- sapply(as.character(unique(fish_sp_retained[,1])), function(x) {
  
  length(unique(diets[diets$fish_sp == x,]$site_code))
  
})

#hist(n_locations)
table(n_locations)


## perform clustering

clusters <- parallel::mclapply(1:1000, function(i) {
  
  cat(i, "\n")
  set.seed(i)
  bipartite::computeModules(t(diets_p))
  
}, mc.cores = 50)


pred_modules <- parallel::mclapply(clusters, function(x) {
  
  mat_memb <- x@modules[2:nrow(x@modules), 2:ncol(x@modules)]
  
  sapply(1:ncol(mat_memb), function (i) {
    
    which(mat_memb [,i] == max(mat_memb [,i]))
  })
  
}, mc.cores = 50)


combin <- expand.grid(1:1000, 1:1000)

vii <- parallel::mclapply(1:nrow(combin), function (x) {
  
  mcclust::vi.dist(pred_modules[[combin[x,1]]], pred_modules[[combin[x,2]]], parts = FALSE, base = 2)
  
}, mc.cores = 50)



vii_clust <- data.frame(x = combin[,1], y = combin[,2], vii = do.call(rbind,vii))


cast_vii <- reshape::cast(vii_clust, x ~ y, value = "vii")

medoid <- which(colMeans(cast_vii) == min(colMeans(cast_vii)))

n_cast_vii = (cast_vii-min(cast_vii))/(max(cast_vii)-min(cast_vii))

medoid_guilds = clusters[[medoid]]

resList <- bipartite::listModuleInformation(medoid_guilds)
#printoutModuleInformation(model)


getModules_predators <- lapply(1:length(resList[[2]]), function(x) {
  
  cbind(resList[[2]][[x]][[2]], rep(x, length(resList[[2]][[x]][[2]])))
  
  
})

getModules_prey <- lapply(1:length(resList[[2]]), function(x) {
  
  cbind(resList[[2]][[x]][[1]], rep(x, length(resList[[2]][[x]][[1]])))
  
  
})


clusters <- do.call(rbind, getModules_predators)

clust <- t(sapply(1:dim(clusters)[1], function(x) {
  
  
  unlist(strsplit(clusters[,1][x], "\\-"))
}))

clust = cbind(clust, clusters[,2])

#number of multiple sites species

length(clust[,1][duplicated(clust[,1])])

test_multiple_clusters <- lapply(clust[,1][duplicated(clust[,1])], function(x) {
  
  data.frame(species = x, identical = dim(table(clust[which(clust[,1] == x), 3])), n.sites = length(clust[which(clust[,1] == x), 3]))
  
})


sum(test_multiple_clusters[,2]) - sum(test_multiple_clusters[,3])

test_multiple_clusters <- do.call(rbind, test_multiple_clusters)









