
diets <- read.csv("data/data_guts.csv", sep=",", dec=".", row.names=1, na = "NA") 

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
                               "Dasyatis americana" = "Hypanus americanus", 
                               "Apogon hyalosoma" ="Yarica hyalosoma"
)


# ----> create standardized interaction matrix diets_p with a sample size of at least 10 guts

diets$sp <- paste(diets$fish_sp)
diets_p <- reshape::cast(diets, sp ~ grp, value = "quantity", fun.aggregate=sum)
diets_p[is.na(diets_p)] = 0

diets_samplesize <- sapply(as.character(unique(diets_p$sp)), function(x) {sum(diets[diets$sp == x,]$nb_guts)})
min_samplesize <- 10

diets_p <- diets_p[diets_p$sp %in% names(diets_samplesize[diets_samplesize >=min_samplesize]),]
species <- diets_p[,1]
diets_p <- diets_p[,2:ncol(diets_p)]

sum_row <- rowSums(diets_p)
diets_p <- lapply(1:length(sum_row), function(r) {diets_p[r,]/sum_row[r]})
diets_p <- do.call(rbind, diets_p)

rownames(diets_p) <- species 
diets_p <- na.omit(diets_p)

# go back to long format 

diets_p$sp <- rownames(diets_p)                       
d_long <- reshape::melt(diets_p, id=c("sp"))
colnames(d_long) <- c("predator", "prey", "quantity")
d_long$predator <- gsub(" ", "_", d_long$predator)

# get phylogenetic data for species whose trophic guilds have been extrapolated

sp_extra <- read.csv("results/extrapolation_trophic_guilds.csv", sep=",", dec=".",  na = "NA")
ft <- fishtree::fishtree_complete_phylogeny(sp_extra$species)

#create maximum credibility tree out of 100 phylogenetic trees
mct <- phangorn::maxCladeCred(ft)

#get phylogenetic eigenvector maps according to Guénard et al. (2013) 
eig <- MPSEM::PEM.build(MPSEM::Phylo2DirectedGraph(mct))

save(eig, file="results/phylo_eig.RData")

all_species <- read.csv("data/all_reef_fish.csv", sep=",", dec=".",  na = "NA")
our_sp <- unique(d_long$predator)

#sizemax <- rfishbase::species(gsub("_" , " ", sp_extra$species))
#sizemax <- pull(species(gsub("_", " ", (diet_cat$species)), fields = "Length"), Length)
#size_our_sp <- sizemax[sizemax$Species %in% gsub("_" , " ", our_sp ),]$Length

# add missing ones manually
#diet_cat[diet_cat$species == "Antennablennius velifer", "sizemax"] <- 7
#diet_cat[diet_cat$species == "Ostorhinchus aroubiensis", "sizemax"] <- 12
#diet_cat[diet_cat$species == "Acanthurus bahianus", "sizemax"] <- 38

#pcoa_our_sp <- pcoa$points[rownames(pcoa$points) %in% our_sp,] 

em_e <- eig$u
em_our_sp <- em_e[rownames(em_e) %in% our_sp,] 

#modify here fr the species list

#d_long_matching <- d_long[d_long$predator %in% rownames(pcoa_our_sp),]
d_long_matching <- d_long[d_long$predator %in% rownames(em_e),]


#d_long_phylo <- mclapply(d_long_matching$predator, function(x) {
#  
#            pcoa_our_sp[rownames(pcoa_our_sp) == x]
#  
#}, mc.cores = 10)

d_long_phylo <- mclapply(d_long_matching$predator, function(x) {
  
  em_our_sp[rownames(em_our_sp) == x]
  
}, mc.cores = 10)

maxsize <- mclapply(d_long_matching$predator, function(x) {
  
  all_species[all_species$species == x,]$sizemax
  
}, mc.cores = 10)

maxsize <- do.call(rbind,maxsize)

d_long_matching$size <- maxsize

cal_data <- cbind(d_long_matching, do.call(rbind,d_long_phylo))
cal_data <- cal_data[,-1]

weights <- cal_data$quantity
cal_data$quantity <- ifelse(cal_data$quantity>0,1,0)

#select 500 eigenvectors
cal_data <- cal_data[,1:503]

data.table::fwrite(cal_data, "data/cal_data_h2o_mod.csv")
