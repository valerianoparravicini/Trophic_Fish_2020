
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
  
  df$x <- as.factor(as.character(df$x))
  df$y <- as.factor(as.character(df$y))
  
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



# -------------> plots figure 1: overall accuracy and accuracy per guild

h = hist(overall_accuracy, plot=FALSE)

library(ggplot2)
library(tidyverse)
library(fishualize)

df <- data.frame(
  overall_accuracy = overall_accuracy,
  groups = cut(overall_accuracy, breaks = h$breaks)
) %>% group_by(groups) %>%
  dplyr::summarize(dens = n()/23)

gg_accuracies <- reshape::melt(accuracies_class_values)

p1 <- ggplot(df) +
  geom_bar(aes(x = groups, y = dens), stat = "identity", fill = "grey", color = "black",
           width = 1,
           position = position_nudge(x = 0.5)) +
  labs(y = "Experts' proportion", x = "Agreement") +
  scale_x_discrete(labels = c("0.6", "", "0.7", "", "0.8", "", "0.9", ""), expand = c(0.01, 0.1)) +
  scale_y_continuous(expand = c(0, 0.005) )+
  theme_bw() + 
  geom_vline(xintercept = 4.6, linetype = 2, color = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none", panel.border = element_blank(), axis.line = element_line()) 


p2 <- 
  ggplot(gg_accuracies) +
  geom_boxplot(aes(x = variable, y = value, fill = variable), alpha = 0.7, outlier.alpha = 0) +
  geom_jitter(aes(x = variable, y = value), width = 0.1, size = 0.8, alpha = 0.5) +
  labs(x = "Trophic guild", y = "Agreement") +
  scale_fill_fish_d(option = "Melichthys_vidua") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none", axis.line = element_line(), panel.border = element_blank()) +
  scale_x_discrete(labels = c("H", "I", "P", "PK", "O"))

p <- cowplot::plot_grid(p1,p2, nrow = 1, labels = c("a)", "b)"), label_fontface = "plain", hjust = 0)


ggsave(plot = p, "output/plots/Fig_1_Agreement.pdf", width = 10, height = 5)
#ggsave(plot = p, "output/plots/Fig_1_Agreement.tiff", width = 10, height = 5)


