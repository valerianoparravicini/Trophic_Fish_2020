

expert <- read.csv("data/converted_experts_classification.csv")


# -------------------> select minimum of species for comparison 


expert <- expert[,9:ncol(expert)]

#generate experts pairs
pairs <- combn(colnames(expert), 2)

#n species in common between pairs
common_sp <- sapply(1:ncol(pairs), function (x) {
  
  combination <- pairs[,x] 
  common <- which(!is.na(expert[,combination[1]])) %in%  which(!is.na(expert[,combination[2]]))
  table(common)["TRUE"]
})#eo sapply 

min_n_sp <- 50
pairs_s <- pairs[,which(common_sp > min_n_sp)]



# ------------------------> calculate agreement

nsp_compared <-  sapply(1:ncol(pairs_s), function(x) {     
  
df <- na.omit(data.frame( x = expert[,pairs_s[1,x]], y = expert[,pairs_s[2,x]]))

common_guilds <- intersect(unique(df[,1]),unique(df[,2]))

df <- df[which(df$x %in% common_guilds),]
df <- df[which(df$y %in% common_guilds),]

dim(df)[1]

})

overall_agreement <- lapply(1:ncol(pairs_s), function(x) {
  

    df <- na.omit(data.frame( x = expert[,pairs_s[1,x]], y = expert[,pairs_s[2,x]]))
    
    common_guilds <- intersect(unique(df[,1]),unique(df[,2]))
    
    df <- df[which(df$x %in% common_guilds),]
    df <- df[which(df$y %in% common_guilds),]
    
    xtab <- table(df[,1], df[,2])
    caret::confusionMatrix(xtab)
    
})#eo lapply


#calculate overall accuracy
overall_accuracy <- sapply(1:length(overall_agreement), function(x) {
  
  overall_agreement [[x]]$overall[1]
  
})#eo lapply


weighted_agreement <- sum(overall_accuracy * nsp_compared) / sum(nsp_compared)

weighted_agreement


#accuracy per guild

accuracies_class <- lapply(1:length(overall_agreement), function(x) {

  cat(x,"\n")  
overall_agreement[[x]]$byClass[,11]
  
})#eo lapply

accuracies_class_df <- lapply(accuracies_class, function(x) {
  
  as.data.frame(t(x))
  
})#eo lapply

accuracies_class_values <- plyr:::rbind.fill(accuracies_class_df)

accuracies_class_values_wm <- apply(accuracies_class_values,2, function(x)
  
                          x * nsp_compared 
  )

accuracies_class_values_wm <- apply(accuracies_class_values_wm,2, function(x) {
  
        sum(x[!is.na(x)])/sum(nsp_compared[!is.na(x)])  
  
})



# -------------> plots figure 1

hist(overall_accuracy)
boxplot(accuracies_class_values )

