load_and_clean <- function(diet_cat){
  #all taxa fishbase
  tax <- load_taxa()
  
  diet_cat <- diet_cat %>% as.data.frame()
  
  #load data on diet clusters
  colnames(diet_cat) <- c("species", "cat") 
  
  # check for name errors
  #error <- fishflux::name_errors(diet_cat$species)
  
  # fix wrong species names
  diet_cat$species <- 
    dplyr::recode(diet_cat$species, 
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
  
  diet_cat <- diet_cat %>% mutate(Species = species) %>% left_join(tax)
  return(diet_cat)
}

modify_dietcat <- function(diet_cat, tree){
  
  # some species not in tree (among which the ray) --> 7
  diet_cat <- diet_cat[diet_cat$species %in% gsub("_", " ", tree[[1]]$tip.label), ] %>% unique()
  # Remove fish that are not in tree 
  
  # add lengths
  # Get max length for all fish through fishbase
  diet_cat$sizemax <- pull(species(gsub("_", " ", (diet_cat$species)), fields = "Length"), Length)
  # add missing ones manually
  diet_cat[diet_cat$species == "Antennablennius velifer", "sizemax"] <- 7
  diet_cat[diet_cat$species == "Ostorhinchus aroubiensis", "sizemax"] <- 12
  diet_cat[diet_cat$species == "Acanthurus bahianus", "sizemax"] <- 38
  
  
  diet_cat$phylo <- gsub(" ", "_", diet_cat$species)
  return(diet_cat)
}

delta_stat <- function(dietcat){
  #SOME PARAMETERS... 
  lambda0 <- 0.1   #rate parameter of the proposal 
  se      <- 0.5   #standard deviation of the proposal
  sim     <- 10000 #number of iterations
  thin    <- 10    #we kept only each 10th iterate 
  burn    <- 100   #100 iterates are burned-in
  #subset of species for which phylotree is sure --> 535
  tree <- fishtree_phylogeny(dietcat$species)
  dietcat <- filter(dietcat, species %in% gsub("_"," ", tree$tip.label)) %>% unique()
  treedata <- data.frame(phylo = tree$tip.label) %>%
    left_join(select(dietcat, phylo, cat))
  trait <- as.numeric(as.character(treedata$cat))
  #CALCULATE DELTA A
  d <- delta(trait,tree,lambda0,se,sim,thin,burn)
  
  return(d)
}
##########
nullm <- function(dietcat){
  lambda0 <- 0.1   #rate parameter of the proposal 
  se      <- 0.5   #standard deviation of the proposal
  sim     <- 1000 #number of iterations
  thin    <- 10    #we kept only each 10th iterate 
  burn    <- 100   #100 iterates are burned-in
  tree <- fishtree_phylogeny(dietcat$species)
  dietcat <- filter(dietcat, species %in% gsub("_"," ", tree$tip.label)) %>% unique()
  trait <- sample(as.character(dietcat$cat)) # random shuffle
  deltaA <- delta(trait,tree,lambda0,se,sim,thin,burn)
  return(deltaA)
}

nullm_safe <- possibly(nullm, otherwise = NA)

##### run nullm 100 times #### 

delta_nullm <- function(dietcat, n = 100){
  list <- mclapply(1:n, function(x){
    nullm_safe(dietcat)
  }, mc.cores = 50)
  return(unlist(list))
}



##### regression #####
run_dietmodel <- function(dietcat, tree){
  # Prerequesite for phylo brms is the covariance matrix 
  # --> Do this for each tree
  corphy_list <- parallel::mclapply(tree, function(x){
    phy2 <- x %>% ape::chronoMPL()
    inv.phylo <- MCMCglmm::inverseA(phy2, nodes = "TIPS", scale = TRUE)
    A <- base::solve(inv.phylo$Ainv)
    rownames(A) <- rownames(inv.phylo$Ainv)
    A <- A[sort(rownames(A)), sort(colnames(A))]
    return(A)
  }, mc.cores = 50)
  
  # Transform list of covariance matrices to 3d array
  corphy_array <- simplify2array(corphy_list, higher = TRUE) 
  dim(corphy_array)
  
  # Summarize these matrices 
  corphy_m <- apply(corphy_array, 1:2, mean)
  corphy_sd <- apply(corphy_array, 1:2, sd)
  
  fit <- brm(
    formula= cat ~ log(sizemax) + (1|phylo),
    family = categorical (link = logit ),
    data = dietcat,
    cores = 3, chains = 3, warmup = 3000, iter = 6000,     
    cov_ranef = list(phylo = corphy_m), control = list(adapt_delta = 0.99, max_treedepth = 15)
  )
  return(fit)
}

extract_r_phylo <- function(fit_diet){
  
  r_phylo <- tidybayes::spread_draws(fit_diet, r_phylo__mu2[species,Intercept], 
                                     r_phylo__mu3[species,Intercept], 
                                     r_phylo__mu4[species,Intercept],
                                     r_phylo__mu5[species,Intercept],
                                     r_phylo__mu6[species,Intercept],
                                     r_phylo__mu7[species,Intercept],
                                     r_phylo__mu8[species,Intercept])
  return(r_phylo)
}

extract_b <- function(fit_diet){
  
  b <- tidybayes::spread_draws(fit_diet, b_mu2_Intercept, b_mu3_Intercept,
                               b_mu4_Intercept, b_mu5_Intercept, 
                               b_mu6_Intercept, b_mu7_Intercept,
                               b_mu8_Intercept,
                               b_mu2_logsizemax, b_mu3_logsizemax,
                               b_mu4_logsizemax, b_mu5_logsizemax,
                               b_mu6_logsizemax, b_mu7_logsizemax,
                               b_mu8_logsizemax)
  return(b)
}

### test extrapolation method
test_extrapolate <- function(b, r_phylo, dietcat, tree){
  
  #set.seed(seed = 98765)
  sdraws <- sample(1:9000, size = 2000)
  
  extrap <- lapply(1:nrow(dietcat), function(x){ #nrow(dietcat)
    print(x)
    
    diet_prune <- dietcat[-x,]
    r_phylo <- filter(r_phylo, species %in% diet_prune$phylo)
    
    all_draw_diet <- parallel::mclapply(sdraws, function(x){
      
      print(paste("draw", x))
      draw = x
      
      set.seed(seed = x)
      i <- sample(1:100, 1)
      
      tree1 <- tree[[i]]
      
      b_sub <- filter(b, .draw == draw)
      b_int2 <- b_sub$b_mu2_Intercept
      b_int3 <- b_sub$b_mu3_Intercept
      b_int4 <- b_sub$b_mu4_Intercept
      b_int5 <- b_sub$b_mu5_Intercept
      b_int6 <- b_sub$b_mu6_Intercept
      b_int7 <- b_sub$b_mu7_Intercept
      b_int8 <- b_sub$b_mu8_Intercept
      
      b_logsize2 <- b_sub$b_mu2_logsizemax
      b_logsize3 <- b_sub$b_mu3_logsizemax
      b_logsize4 <- b_sub$b_mu4_logsizemax
      b_logsize5 <- b_sub$b_mu5_logsizemax
      b_logsize6 <- b_sub$b_mu6_logsizemax
      b_logsize7 <- b_sub$b_mu7_logsizemax
      b_logsize8 <- b_sub$b_mu8_logsizemax
      
      r_spec <- dplyr::filter(r_phylo, .draw == draw)
      
      # group 2
      print("group 2")
      
      st2 <- r_spec$r_phylo__mu2
      names(st2) <- r_spec$species
      
      r_phy <- picante::phyEstimate(tree1, st2)
      
      #extrapolated
      rphy_m2 <- data.frame(
        species = rownames(r_phy),
        r_phylo2 = r_phy[,1])
      #known from model
      rphy_2 <- data.frame(
        species = names(st2),
        r_phylo2 = st2
      )
      rphy2 <- rbind(rphy_m2, rphy_2)
      
      # group 3
      print("group 3")
      
      st3 <- r_spec$r_phylo__mu3
      names(st3) <- r_spec$species
      
      r_phy <- picante::phyEstimate(tree1, st3)
      
      #extrapolated
      rphy_m3 <- data.frame(
        species = rownames(r_phy),
        r_phylo3 = r_phy[,1])
      #known from model
      rphy_3 <- data.frame(
        species = names(st3),
        r_phylo3 = st3
      )
      rphy3 <- rbind(rphy_m3, rphy_3)
      
      # group 4
      print("group 4")
      
      st4 <- r_spec$r_phylo__mu4
      names(st4) <- r_spec$species
      
      r_phy <- picante::phyEstimate(tree1, st4)
      
      #extrapolated
      rphy_m4 <- data.frame(
        species = rownames(r_phy),
        r_phylo4 = r_phy[,1])
      #known from model
      rphy_4 <- data.frame(
        species = names(st4),
        r_phylo4 = st4
      )
      rphy4 <- rbind(rphy_m4, rphy_4)
      
      # group 5
      print("group 5")
      
      st5 <- r_spec$r_phylo__mu5
      names(st5) <- r_spec$species
      
      r_phy <- picante::phyEstimate(tree1, st5)
      
      #extrapolated
      rphy_m5 <- data.frame(
        species = rownames(r_phy),
        r_phylo5 = r_phy[,1])
      #known from model
      rphy_5 <- data.frame(
        species = names(st5),
        r_phylo5 = st5
      )
      rphy5 <- rbind(rphy_m5, rphy_5)
      
      # group 6
      print("group 6")
      
      st6 <- r_spec$r_phylo__mu6
      names(st6) <- r_spec$species
      
      r_phy <- picante::phyEstimate(tree1, st6)
      
      #extrapolated
      rphy_m6 <- data.frame(
        species = rownames(r_phy),
        r_phylo6 = r_phy[,1])
      #known from model
      rphy_6 <- data.frame(
        species = names(st6),
        r_phylo6 = st6
      )
      rphy6 <- rbind(rphy_m6, rphy_6)
      
      
      # group 7
      print("group 7")
      
      st7 <- r_spec$r_phylo__mu7
      names(st7) <- r_spec$species
      
      r_phy <- picante::phyEstimate(tree1, st7)
      
      #extrapolated
      rphy_m7 <- data.frame(
        species = rownames(r_phy),
        r_phylo7 = r_phy[,1])
      #known from model
      rphy_7 <- data.frame(
        species = names(st7),
        r_phylo7 = st7
      )
      rphy7 <- rbind(rphy_m7, rphy_7)
      
      # group 8
      print("group 8")
      
      st8 <- r_spec$r_phylo__mu8
      names(st8) <- r_spec$species
      
      r_phy <- picante::phyEstimate(tree1, st8)
      
      #extrapolated
      rphy_m8 <- data.frame(
        species = rownames(r_phy),
        r_phylo8 = r_phy[,1])
      #known from model
      rphy_8 <- data.frame(
        species = names(st8),
        r_phylo8 = st8
      )
      rphy8 <- rbind(rphy_m8, rphy_8)
      
      
      diet_ex <- left_join(rphy_m2, rphy_m3) %>% left_join(rphy_m4) %>% left_join(rphy_m5) %>% left_join(rphy_m6) %>% 
        left_join(rphy_m7) %>% left_join(rphy_m8)  %>% unique() %>%
        left_join(dplyr::select(dietcat, Family, species = phylo, sizemax)) %>%
        mutate(b_logsize2 = b_logsize2, b_logsize3 = b_logsize3,
               b_logsize4 = b_logsize4, b_logsize5 = b_logsize5, 
               b_logsize6 = b_logsize6, b_logsize7 = b_logsize7, 
               b_logsize8 = b_logsize8, 
               
               b_int2 = b_int2, b_int3 = b_int3, b_int4 = b_int4, b_int5 = b_int5,
               b_int6 = b_int6, b_int7 = b_int7, b_int8 = b_int8) %>%
        mutate(mu1 = (0),
               mu2 = (r_phylo2 + b_int2 + (b_logsize2 * log(sizemax))), 
               mu3 = (r_phylo3 + b_int3 + (b_logsize3 * log(sizemax))), 
               mu4 = (r_phylo4 + b_int4 + (b_logsize4 * log(sizemax))), 
               mu5 = (r_phylo5 + b_int5 + (b_logsize5 * log(sizemax))),
               mu6 = (r_phylo6 + b_int6 + (b_logsize6 * log(sizemax))),
               mu7 = (r_phylo7 + b_int7 + (b_logsize7 * log(sizemax))),
               mu8 = (r_phylo8 + b_int8 + (b_logsize8 * log(sizemax)))) %>% 
        mutate(draw = draw)
      
      return(diet_ex)
      
    }, mc.cores = 45) %>% plyr::ldply()
  }) %>% bind_rows()
  
  
  
  probs <-
    parallel::mclapply(1:nrow(extrap), function(x){
      sub <- extrap[x,]
      sm <- softmax(c(sub$mu1, sub$mu2, sub$mu3, sub$mu4, sub$mu5, sub$mu6, sub$mu7, sub$mu8))
      result <- data.frame(species = sub$species, family = sub$Family, p1 = sm[1], p2 = sm[2], p3 = sm[3], p4 = sm[4], p5 = sm[5],
                           p6 = sm[6], p7 = sm[7], p8 = sm[8])  
      result$nentropy <- 1 - nentropy(result[,3:10])
      return(result)
    }, mc.cores = 50) %>% plyr::ldply()
  
  diet_predict <- probs  %>%
    group_by(family, species) %>% 
    dplyr::summarise(p1_m = mean(p1, na.rm = TRUE), p1_sd = sd(p1, na.rm = TRUE), p2_m = mean(p2, na.rm = TRUE), p2_sd = sd(p2, na.rm = TRUE),
                     p3_m = mean(p3, na.rm = TRUE), p3_sd = sd(p3, na.rm = TRUE), p4_m = mean(p4, na.rm = TRUE), p4_sd = sd(p4, na.rm = TRUE),
                     p5_m = mean(p5, na.rm = TRUE), p5_sd = sd(p5, na.rm = TRUE), p6_m = mean(p6, na.rm = TRUE), p6_sd = sd(p6, na.rm = TRUE),
                     p7_m = mean(p7, na.rm = TRUE), p7_sd = sd(p7, na.rm = TRUE),  p8_m = mean(p8, na.rm = TRUE), p8_sd = sd(p8, na.rm = TRUE),
                     nentropy_m = mean(nentropy), nentropy_sd = sd(nentropy)) %>% ungroup() %>%
    left_join(select(dietcat, species = phylo, cat)) %>%
    mutate(sd_tot = sqrt(p1_sd^2 + p2_sd^2 + p3_sd^2 + p4_sd^2 + p5_sd^2 + p6_sd^2 + p7_sd^2 + p8_sd^2))
  
  
  cats <- as.factor(as.character(c(1:8)))
  diet_predict$pred_cat <- apply(select(diet_predict, p1_m, p2_m, p3_m, p4_m, p5_m, p6_m, p7_m, p8_m), 1, function(x){cats[which(x == max(x))]})
  diet_predict$match <- diet_predict$cat == diet_predict$pred_cat
  perc_match <- sum(diet_predict$match)/nrow(diet_predict)
  
  #correct family names
  diet_predict[diet_predict$family == "Scaridae", "family"] <- "Labridae"
  diet_predict[diet_predict$family == "Microdesmidae", "family"] <- "Gobiidae"
  
  cont_table <- table(diet_predict$pred_cat, diet_predict$cat)
  
  
  return(list(diet_predict = diet_predict, perc_match = perc_match, cont_table = cont_table))
}


### extrapolate ####

extrapolate <- function(b, r_phylo){
  
  
  allsp <- read.csv("data/all_reef_fish.csv")
  
  alltree <- fishtree_complete_phylogeny(allsp$Species)
  
  
  #set.seed(seed = 98765)
  sdraws <- sample(1:9000, size = 2000)
  
  
  extrap <- parallel::mclapply(sdraws, function(x){
    
    print(paste("draw", x))
    draw = x
    
    #set.seed(seed = x)
    i <- sample(1:100, 1)
    
    tree1 <- alltree[[i]]
    
    b_sub <- filter(b, .draw == draw)
    b_int2 <- b_sub$b_mu2_Intercept
    b_int3 <- b_sub$b_mu3_Intercept
    b_int4 <- b_sub$b_mu4_Intercept
    b_int5 <- b_sub$b_mu5_Intercept
    b_int6 <- b_sub$b_mu6_Intercept
    b_int7 <- b_sub$b_mu7_Intercept
    b_int8 <- b_sub$b_mu8_Intercept
    
    b_logsize2 <- b_sub$b_mu2_logsizemax
    b_logsize3 <- b_sub$b_mu3_logsizemax
    b_logsize4 <- b_sub$b_mu4_logsizemax
    b_logsize5 <- b_sub$b_mu5_logsizemax
    b_logsize6 <- b_sub$b_mu6_logsizemax
    b_logsize7 <- b_sub$b_mu7_logsizemax
    b_logsize8 <- b_sub$b_mu8_logsizemax
    
    r_spec <- dplyr::filter(r_phylo, .draw == draw)
    
    # group 2
    print("group 2")
    
    st2 <- r_spec$r_phylo__mu2
    names(st2) <- r_spec$species
    
    r_phy <- picante::phyEstimate(tree1, st2)
    
    #extrapolated
    rphy_m2 <- data.frame(
      species = rownames(r_phy),
      r_phylo2 = r_phy[,1])
    #known from model
    rphy_2 <- data.frame(
      species = names(st2),
      r_phylo2 = st2
    )
    rphy2 <- rbind(rphy_m2, rphy_2)
    
    # group 3
    print("group 3")
    
    st3 <- r_spec$r_phylo__mu3
    names(st3) <- r_spec$species
    
    r_phy <- picante::phyEstimate(tree1, st3)
    
    #extrapolated
    rphy_m3 <- data.frame(
      species = rownames(r_phy),
      r_phylo3 = r_phy[,1])
    #known from model
    rphy_3 <- data.frame(
      species = names(st3),
      r_phylo3 = st3
    )
    rphy3 <- rbind(rphy_m3, rphy_3)
    
    # group 4
    print("group 4")
    
    st4 <- r_spec$r_phylo__mu4
    names(st4) <- r_spec$species
    
    r_phy <- picante::phyEstimate(tree1, st4)
    
    #extrapolated
    rphy_m4 <- data.frame(
      species = rownames(r_phy),
      r_phylo4 = r_phy[,1])
    #known from model
    rphy_4 <- data.frame(
      species = names(st4),
      r_phylo4 = st4
    )
    rphy4 <- rbind(rphy_m4, rphy_4)
    
    # group 5
    print("group 5")
    
    st5 <- r_spec$r_phylo__mu5
    names(st5) <- r_spec$species
    
    r_phy <- picante::phyEstimate(tree1, st5)
    
    #extrapolated
    rphy_m5 <- data.frame(
      species = rownames(r_phy),
      r_phylo5 = r_phy[,1])
    #known from model
    rphy_5 <- data.frame(
      species = names(st5),
      r_phylo5 = st5
    )
    rphy5 <- rbind(rphy_m5, rphy_5)
    
    # group 6
    print("group 6")
    
    st6 <- r_spec$r_phylo__mu6
    names(st6) <- r_spec$species
    
    r_phy <- picante::phyEstimate(tree1, st6)
    
    #extrapolated
    rphy_m6 <- data.frame(
      species = rownames(r_phy),
      r_phylo6 = r_phy[,1])
    #known from model
    rphy_6 <- data.frame(
      species = names(st6),
      r_phylo6 = st6
    )
    rphy6 <- rbind(rphy_m6, rphy_6)
    
    
    # group 7
    print("group 7")
    
    st7 <- r_spec$r_phylo__mu7
    names(st7) <- r_spec$species
    
    r_phy <- picante::phyEstimate(tree1, st7)
    
    #extrapolated
    rphy_m7 <- data.frame(
      species = rownames(r_phy),
      r_phylo7 = r_phy[,1])
    #known from model
    rphy_7 <- data.frame(
      species = names(st7),
      r_phylo7 = st7
    )
    rphy7 <- rbind(rphy_m7, rphy_7)
    
    # group 8
    print("group 8")
    
    st8 <- r_spec$r_phylo__mu8
    names(st8) <- r_spec$species
    
    r_phy <- picante::phyEstimate(tree1, st8)
    
    #extrapolated
    rphy_m8 <- data.frame(
      species = rownames(r_phy),
      r_phylo8 = r_phy[,1])
    #known from model
    rphy_8 <- data.frame(
      species = names(st8),
      r_phylo8 = st8
    )
    rphy8 <- rbind(rphy_m8, rphy_8)
    
    
    diet_ex <- left_join(rphy2, rphy3) %>% left_join(rphy4) %>% left_join(rphy5) %>% left_join(rphy6) %>% 
      left_join(rphy7) %>% left_join(rphy8)  %>% unique() %>%
      inner_join(dplyr::select(allsp, Family, species, sizemax)) %>%
      mutate(b_logsize2 = b_logsize2, b_logsize3 = b_logsize3,
             b_logsize4 = b_logsize4, b_logsize5 = b_logsize5, 
             b_logsize6 = b_logsize6, b_logsize7 = b_logsize7, 
             b_logsize8 = b_logsize8, 
             
             b_int2 = b_int2, b_int3 = b_int3, b_int4 = b_int4, b_int5 = b_int5,
             b_int6 = b_int6, b_int7 = b_int7, b_int8 = b_int8) %>%
      mutate(mu1 = (0),
             mu2 = (r_phylo2 + b_int2 + (b_logsize2 * log(sizemax))), 
             mu3 = (r_phylo3 + b_int3 + (b_logsize3 * log(sizemax))), 
             mu4 = (r_phylo4 + b_int4 + (b_logsize4 * log(sizemax))), 
             mu5 = (r_phylo5 + b_int5 + (b_logsize5 * log(sizemax))),
             mu6 = (r_phylo6 + b_int6 + (b_logsize6 * log(sizemax))),
             mu7 = (r_phylo7 + b_int7 + (b_logsize7 * log(sizemax))),
             mu8 = (r_phylo8 + b_int8 + (b_logsize8 * log(sizemax)))) %>% 
      mutate(draw = draw)
    
    return(diet_ex)
    
  }, mc.cores = 6) %>% plyr::ldply()
  
  
  
  probs <-
    parallel::mclapply(1:nrow(extrap), function(x){
      sub <- extrap[x,]
      sm <- softmax(c(sub$mu1, sub$mu2, sub$mu3, sub$mu4, sub$mu5, sub$mu6, sub$mu7, sub$mu8))
      result <- data.frame(species = sub$species, family = sub$Family, p1 = sm[1], p2 = sm[2], p3 = sm[3], p4 = sm[4], p5 = sm[5],
                           p6 = sm[6], p7 = sm[7], p8 = sm[8])  
      result$nentropy <- 1 - nentropy(result[,3:10])
      return(result)
    }, mc.cores = 10) %>% plyr::ldply()
  
  
  write.csv(extrap, "results/dietalldraws_reg.csv", row.names = FALSE)
  write.csv(probs, "results/dietalldraws_probs.csv", row.names = FALSE)
  
  diet_predict <- probs  %>%
    group_by(family, species) %>% 
    dplyr::summarise(p1_m = mean(p1, na.rm = TRUE), p1_sd = sd(p1, na.rm = TRUE), p2_m = mean(p2, na.rm = TRUE), p2_sd = sd(p2, na.rm = TRUE),
                     p3_m = mean(p3, na.rm = TRUE), p3_sd = sd(p3, na.rm = TRUE), p4_m = mean(p4, na.rm = TRUE), p4_sd = sd(p4, na.rm = TRUE),
                     p5_m = mean(p5, na.rm = TRUE), p5_sd = sd(p5, na.rm = TRUE), p6_m = mean(p6, na.rm = TRUE), p6_sd = sd(p6, na.rm = TRUE),
                     p7_m = mean(p7, na.rm = TRUE), p7_sd = sd(p7, na.rm = TRUE),  p8_m = mean(p8, na.rm = TRUE), p8_sd = sd(p8, na.rm = TRUE),
                     nentropy_m = mean(nentropy), nentropy_sd = sd(nentropy)) %>% ungroup() %>%
    mutate(sd_tot = sqrt(p1_sd^2 + p2_sd^2 + p3_sd^2 + p4_sd^2 + p5_sd^2 + p6_sd^2 + p7_sd^2 + p8_sd^2))
  
  
  cats <- as.factor(as.character(c(1:8)))
  diet_predict$pred_cat <- apply(select(diet_predict, p1_m, p2_m, p3_m, p4_m, p5_m, p6_m, p7_m, p8_m), 1, function(x){cats[which(x == max(x))]})
  
  #correct family names
  diet_predict[diet_predict$family == "Scaridae", "family"] <- "Labridae"
  diet_predict[diet_predict$family == "Microdesmidae", "family"] <- "Gobiidae"
  
  write.csv(diet_predict, "results/dietall_summary.csv", row.names = FALSE)
  
  return(diet_predict)
}

## softmax function
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}

softmax <- function (x) {
  exp(x - logsumexp(x))
}

#entropy
nentropy <- function(prob) {
  
  k              <- ncol(prob)                       #number of states
  prob[prob>1/k] <- prob[prob>1/k]/(1-k) - 1/(1-k)   #state entropies
  tent           <- apply(prob,1,sum)                #node entropy
  
  #correct absolute 0/1
  tent[tent == 0] <- tent[tent == 0] + runif(1,0,1)/10000
  tent[tent == 1] <- tent[tent == 1] - runif(1,0,1)/10000
  
  return(tent)
}


summarize_extrapolation <- function(dietall){
  
  probs <-
    parallel::mclapply(1:nrow(dietall), function(x){
      sub <- dietall[x,]
      sm <- softmax(c(sub$mu1, sub$mu2, sub$mu3, sub$mu4, sub$mu5))
      return(data.frame(species = sub$species, family = sub$Family, p1 = sm[1], p2 = sm[2], p3 = sm[3], p4 = sm[4], p5 = sm[5]))
    }, mc.cores = 30) %>% plyr::ldply()
  
  
  probs$nentropy <- 1 - nentropy(probs[,3:7])
  
  diet_predict <- probs  %>%
    group_by(family, species) %>% 
    dplyr::summarise(p1_m = mean(p1, na.rm = TRUE), p1_sd = sd(p1, na.rm = TRUE), p2_m = mean(p2, na.rm = TRUE), p2_sd = sd(p2, na.rm = TRUE),
                     p3_m = mean(p3, na.rm = TRUE), p3_sd = sd(p3, na.rm = TRUE), p4_m = mean(p4, na.rm = TRUE), p4_sd = sd(p4, na.rm = TRUE),
                     p5_m = mean(p5, na.rm = TRUE), p5_sd = sd(p5, na.rm = TRUE),
                     nentropy_m = mean(nentropy), nentropy_sd = sd(nentropy)) %>% ungroup()
  
  cats <- as.factor(as.character(c(1:5)))
  diet_predict$pred_cat <- apply(select(diet_predict, p1_m, p2_m, p3_m, p4_m, p5_m), 1, function(x){cats[which(x == max(x))]})
  #correct family names
  diet_predict[diet_predict$family == "Scaridae", "family"] <- "Labridae"
  diet_predict[diet_predict$family == "Microdesmidae", "family"] <- "Gobiidae"
  return(diet_predict)
}

##### plot #####

make_treeplot <-
  function(fit_diet, dietcat){
  
  tax <- rfishbase::load_taxa()
  
  tax <- read.csv("data/PFC_taxonomy.csv") %>%
    select(Family = family, Species = genus.species)
  
  
  ppred <- predict(fit_diet) %>% as.data.frame()

  colnames(ppred) <- c("p1","p2","p3", "p4", "p5", "p6", "p7", "p8")

  ppred <- ppred %>% mutate(species = fit_diet$data$phylo)

  rownames(ppred) <- ppred$species

tree1 <- fishtree::fishtree_phylogeny(species = ppred$species)

tree1t <- as_tibble(tree1)
tax$label <- gsub(" ", "_", tax$Species)
tree1t <- left_join(tree1t, tax)
#tree1t[!is.na(tree1t$Family) & tree1t$Family == "Scaridae", "Family"] <- "Labridae"

tree1d <- as.treedata(tree1t)


data <- filter(dietcat, phylo%in% tree1$tip.label)
ppred <- as.data.frame(ppred) %>% filter(species%in% tree1$tip.label)
rownames(ppred) <- ppred$species


p <- ggtree(tree1d, layout = "circular")
p

p <- 
  gheatmap(p, ppred[,7, drop = FALSE], offset = 0.2, width=0.1, 
           colnames_position = "top", font.size=2, color="white") + 
  scale_fill_gradient(low = "white", high = "#BC5380", limits = c(0,1), breaks = c(0, 1),
                      guide = guide_colorbar(ticks = FALSE)) 
p
p <- p + new_scale_fill()


p2 <- gheatmap(p, ppred[,3, drop = FALSE], offset = 15, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#5CD400", limits = c(0,1), breaks = c(0, 1), 
                      guide = guide_colorbar(ticks = FALSE)) 

p2 <- p2 + new_scale_fill()

p3 <- gheatmap(p2, ppred[,4, drop = FALSE], offset = 30, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#1B229B", limits = c(0,1), breaks = c(0, 1), 
                      guide = guide_colorbar(ticks = FALSE))

p3 <- p3 + new_scale_fill() 

p4 <- gheatmap(p3, ppred[,1, drop = FALSE], offset = 45, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#EE3B00", limits = c(0,1), breaks = c(0, 1))
p4 <- p4 + new_scale_fill() 

p5 <- gheatmap(p4, ppred[,8, drop = FALSE], offset = 60, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#46B3B6", limits = c(0,1), breaks = c(0, 1), 
                      guide = guide_colorbar(ticks = FALSE)) +
  theme(legend.position = "none") +
  guides(fill = guide_colourbar(ticks = FALSE))

p5 <- p5 + new_scale_fill()

p6 <- gheatmap(p5, ppred[,6, drop = FALSE], offset = 75, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#E2AF79", limits = c(0,1), breaks = c(0, 1), 
                      guide = guide_colorbar(ticks = FALSE)) +
  theme(legend.position = "none") 

p6 <- p6 + new_scale_fill()

p7 <- gheatmap(p6, ppred[,2, drop = FALSE], offset = 90, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#855EB0", limits = c(0,1), breaks = c(0, 1), 
                      guide = guide_colorbar(ticks = FALSE)) +
  theme(legend.position = "none") 

p7 <- p7 + new_scale_fill()

p8 <- gheatmap(p7, ppred[,5, drop = FALSE], offset = 105, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#F3DC00", limits = c(0,1), breaks = c(0, 1), 
                      guide = guide_colorbar(ticks = FALSE))  +
  theme(legend.position = "none") 

p8 <- p8 + new_scale_fill()

##entropy
nentropy <- function(prob) {
  
  k              <- ncol(prob)                       #number of states
  prob[prob>1/k] <- prob[prob>1/k]/(1-k) - 1/(1-k)   #state entropies
  tent           <- apply(prob,1,sum)                #node entropy
  
  #correct absolute 0/1
  tent[tent == 0] <- tent[tent == 0] + runif(1,0,1)/10000
  tent[tent == 1] <- tent[tent == 1] - runif(1,0,1)/10000
  
  return(tent)
}

data$entropy <- unlist(lapply(1:nrow(ppred), function(x){
  ent <- nentropy(ppred[x, 1:8])}))
data$Species <- gsub("_", " ", data$phylo)
data <- left_join(data, tax)
data$label <- data$phylo
data[data$Family == "Scaridae", "Family"] <- "Labridae"
data[data$Family == "Microdesmidae", "Family"] <- "Gobiidae"
unique(data$Family)

ent_sum <- summarise(group_by(data, Family), ent = mean(entropy))
data <- left_join(data, ent_sum)

tree1t <- as_tibble(tree1)
tree1t <- left_join(tree1t, data)
tree1d <- as.treedata(tree1t)
# tree1t$entropy

ppred[,10] <- (1 - data$ent)

p9 <- gheatmap(p8, ppred[, 10, drop = FALSE], offset = 125, width = 0.05,
               colnames_position = "top", font.size = 2, color = NULL) +
  scale_fill_gradient(low = "lightgrey", high = "black" , limits = c(0,1), breaks = c(0, 1), 
                      guide = guide_colorbar(ticks = FALSE)) +
  theme(legend.position = "none") 



ggsave("figures/diet_phyloplot.pdf", p9, height = 12, width = 12)


p10 <- gheatmap(p8, ppred[, 10, drop = FALSE], offset = 125, width = 0.05,
                colnames_position = "top", font.size = 2, color = NULL) +
  scale_fill_gradient(low = "lightgrey", high = "black", limits = c(0,1), breaks = c(0, 1)) +
  guides(fill = guide_colourbar(ticks = FALSE))
ggsave("figures/diet_phyloplot_legend.pdf",p10, height = 18, width = 12)


###### species label #####

p <- ggtree(tree1d, layout = "circular")
p + geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 5)
ggsave('figures/treefamily.pdf', height = 12, width = 12)

p9 + geom_tiplab(aes(angle = angle), size = 2, offset = 145)
ggsave("figures/diet_phyloplot_species.pdf", height = 12, width = 12)

}


plot_extrap_test <- function(fit_diet, extrap_test){
  
  cont <- extrap_test$cont
  sums <- colSums(cont)
  sums <- data.frame(cat = names(sums), sum = sums)
  cont <- cont %>% as.data.frame()
  colnames(cont) <- c("predict", "cat", "freq")
  cont <- left_join(cont, sums) %>% mutate(prop = round(freq/sum, 2)) %>%
    mutate(predict = 
    dplyr::recode(predict,  "7" = "1", "3" = "2", "4" = "3", "1" = "4",
                  "8" = "5", "6" = "6", "2" = "7", "5" = "8"),
    cat = 
    dplyr::recode(cat,  "7" = "1", "3" = "2", "4" = "3", "1" = "4",
                  "8" = "5", "6" = "6", "2" = "7", "5" = "8"))
  
  
  cont$cat <- factor(cont$cat, levels = as.character(1:8))
  cont$predict <- factor(cont$predict, levels = as.character(1:8))
  
  ggplot(cont) +
    geom_raster(aes(x = cat, y = fct_reorder(predict, -as.numeric(predict)), fill = prop), alpha = 0.8) +
    geom_text(aes(x = cat, y = fct_reorder(predict, -as.numeric(predict)), label = prop), size = 3) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "trophic guild", y = "predicted trophic guild", fill = "proportion per trophic guild") +
    scale_fill_viridis(begin = 0.1) +
    #scale_fill_fish(option = "Ostracion_whitleyi", direction = -1, begin = 0.5, end = 0.9) +
    coord_equal()
  
  ggsave("figures/confusion_matrix.pdf", width = 6, height = 8)
  ggsave("figures/confusion_matrix.tiff", width = 6, height = 8)
  
}


adapt_extrapolation_cats <- function(extrapolation){
  
  result <- 
    extrapolation %>%
  dplyr::mutate(pred_cat = 
                  dplyr::recode(pred_cat,  "7" = "1", "3" = "2", "4" = "3", "1" = "4",
                              "8" = "5", "6" = "6", "2" = "7", "5" = "8")) %>%
    rename(trophic_guild_predicted = pred_cat,
           p1_m = p7_m,
           p1_sd = p7_sd,
           p2_m = p3_m,
           p2_sd = p3_sd,
           p3_m = p4_m,
           p3_sd = p4_sd,
           p4_m = p1_m,
           p4_sd = p1_sd,
           p5_m = p8_m,
           p5_sd = p8_sd,
           p7_m = p2_m,
           p7_sd = p2_sd,
           p8_m = p5_m,
           p8_sd = p5_sd)
write_csv(result, "results/extrapolation_trophic_guilds.csv")
return(result)
}

