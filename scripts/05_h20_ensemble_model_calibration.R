
#### ---------> set up cluster with h20

library(h2o)

h2o.clust <- tryCatch(h2o.init(startH2O = FALSE),error=function(e){
  
  h2o.init(ip = "localhost", port = 54321, startH2O = TRUE,
           forceDL = FALSE, enable_assertions = TRUE, license = NULL,
           nthreads = 40, max_mem_size = "200g", min_mem_size = "10g",
           ice_root = tempdir(), strict_version_check = TRUE,
           proxy = NA_character_, https = FALSE, insecure = FALSE,
           cookies = NA_character_)
  
  #h2o.removeAll() # (Optional) Remove all objects in H2O cluster
  
})


## H2o data import
da <- read.csv2("data/cal_data_h2o_mod.csv", header=TRUE, sep=",", dec=".")

dat <- as.h2o(da)

#h2o.describe(dat)
y <- "quantity"

#extreme gradient boosting
dir.create(paste0(getwd(),"/results/xgb"))
m.xgb <- h2o.xgboost(y = y, training_frame = dat,  nfolds = 10, seed = 1, ntrees = 2000, booster = "dart",
                       keep_cross_validation_predictions=TRUE)
xgb_path <- h2o.saveModel(object = m.xgb, path = paste0(getwd(),"/results/xgb"), force = TRUE)

#gradient boosting
dir.create(paste0(getwd(),"/results/gbm"))
m.gbm <- h2o.gbm(y = y, training_frame = dat,  nfolds = 10, seed = 1, ntrees = 2000,
                 keep_cross_validation_predictions=TRUE)
gbm_path <- h2o.saveModel(object = m.gbm, path = paste0(getwd(),"/results/gbm"), force = TRUE)

#random forest
dir.create(paste0(getwd(),"/results/rf"))
m.rf <- h2o.randomForest(y = y, training_frame = dat,  nfolds = 10, seed = 1, ntrees = 2000,
                         keep_cross_validation_predictions=TRUE)
rf_path <- h2o.saveModel(object = m.rf, path = paste0(getwd(),"/results/rf"), force = TRUE)

print(xgb_path)
print(gbm_path)
print(rf_path)

ensemble <- h2o.stackedEnsemble(y = y, 
                                training_frame = dat ,
                                base_models = list (m.xgb, m.gbm, m.rf), seed=1)


dir.create(paste0(getwd(),"/results/ensemble"))
ens_path <- h2o.saveModel(object =ensemble, path = paste0(getwd(),"/results/ensemble"), force = TRUE)


print(ens_path)

h2o.shutdown()
