library(h2o)
library(pbmcapply)

#set up H20 Cluster
h2o.clust <- tryCatch(h2o.init(startH2O = FALSE),error=function(e){
  
  h2o.init(ip = "localhost", port = 54321, startH2O = TRUE,
           forceDL = FALSE, enable_assertions = TRUE, license = NULL,
           nthreads = 10, max_mem_size = "150g", min_mem_size = "10g",
           ice_root = tempdir(), strict_version_check = TRUE,
           proxy = NA_character_, https = FALSE, insecure = FALSE,
           cookies = NA_character_)
  
  #h2o.removeAll() # (Optional) Remove all objects in H2O cluster
  
})


###### H2o import

filename_ens <- list.files(paste0(getwd(),"/results/ensemble"))
mod <- h2o.loadModel(paste0(getwd(),"/results/ensemble", filename_ens))

da <- read.csv2("data/cal_data_h2o_mod.csv", header=TRUE, sep=",", dec=".")

dat <- as.h2o(da)

y <- "quantity"


pred <- as.data.frame(h2o.predict(mod, dat))


dt <- data.frame(a=1:nrow(dat), b=da$quantity, c=pred$predict)
ot <- PresenceAbsence::optimal.thresholds(dt)
a <- PresenceAbsence::auc(dt)

optim <- lapply(seq(0,1,0.01), function(x) {
  
  c <- PresenceAbsence::cmx(dt, threshold = x)
  tss <- PresenceAbsence::sensitivity(c)+PresenceAbsence::specificity(c)-1
  c(x,tss[1])
}) 

optim <- do.call(rbind, optim)   

optim[which(as.numeric(optim[,2]) == max(as.numeric(optim[,2])))]

c <- PresenceAbsence::cmx(dt, threshold = 0.31)
PresenceAbsence::sensitivity(c)+PresenceAbsence::specificity(c)-1


################ extrapolation on the same species used for trophic guilds

sp_extra <- read.csv("results/extrapolation_trophic_guilds.csv", sep=",", dec=".",  na = "NA")

load("results/phylo_eig.RData")

ext_da <- lapply(sp_extra$species, function(x) {
  
  eig$u[rownames(eig$u) == x, ][1:500]
  
})

ext_da <- do.call(rbind, ext_da)

all_species <- read.csv("data/all_reef_fish.csv", sep=",", dec=".",  na = "NA")

maxsize <- sapply(sp_extra$species, function(x) {
  
  all_species[all_species$species == x,]$sizemax
  
})

extra_data <- cbind(as.data.frame(sp_extra$species),maxsize, ext_da)

cal_data <- read.csv2("data/cal_data_h2o_mod.csv", sep=",", dec=".")

colnames(extra_data) <- c("predator",colnames(cal_data)[3:503])

prey <- unique(cal_data$prey)


extra_data_multiprey <- pbmclapply(1:nrow(extra_data), function(x) {
  
  df <- extra_data[x,]
  
  df_temp <- lapply(1:length(prey), function(y){
    
         temp_temp <- data.frame(c(y, df))
         colnames(temp_temp)[1] <- "prey" 
         temp_temp 
  })
  
 do.call(rbind, df_temp)
  
}, mc.cores=30)

extra_data_multiprey <- do.call(rbind, extra_data_multiprey )

extra_data_multiprey$prey <- sapply(extra_data_multiprey$prey, function(x) {
  
                  prey[x]
})

extra_df <- extra_data_multiprey 

save(extra_df, "data/extrapolation_dataset.RData")

extra_data_multiprey  <- as.h2o(extra_data_multiprey)

extra <- as.data.frame(h2o.predict(mod, extra_data_multiprey))

extra_df$interaction <- round(extra$predict, digits = 2) 
extra_df$interaction[extra_df$interaction<0] <- 0

extrapolation <- data.frame(prey = extra_df$prey, predator = extra_df$predator, interaction = extra_df$interaction)

test_extra <- reshape2::dcast(extrapolation, predator ~ prey, value.var = "interaction")

write.csv2(test_extra, file="results/extrapolation_trophic_interactions.csv")

