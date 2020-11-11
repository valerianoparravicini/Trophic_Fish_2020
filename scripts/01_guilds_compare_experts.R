# ----------------------->>  Comparing Experts' classifications


make_plot1 <- function(classif){
#count experts per species
classif$count <- sapply(1:nrow(classif), function(x) {
  
     8 - length(which(is.na(classif[x,18:25])==T))
  
})#eo sapply


#generate experts pairs
pairs <- combn(colnames(classif)[18:25], 2)

#n species in common between pairs
common_sp <- sapply(1:ncol(pairs), function (x) {
  
            combination <- pairs[,x] 
            common <- which(!is.na(classif[,combination[1]])) %in%  which(!is.na(classif[,combination[2]]))
            table(common)["TRUE"]
})#eo sapply 

min_n_sp <- 30
pairs_s <- pairs[,which(common_sp > min_n_sp)]


overall_agreement <- lapply(1:ncol(pairs_s), function(x) {
  
                            if (pairs_s[1,x] == "rnickwio") {
                              
                             df <- na.omit(data.frame( x = classif[,pairs_s[1,x]], y = classif[,pairs_s[2,x]]))
                             xtab <- table(df[,1], df[,2])
                             xtab <- xtab[,-which(colnames(xtab) == "O")]
                             caret::confusionMatrix(xtab)
                            
                            } else if (pairs_s[2,x] == "rnickwio") { 
                              
                              df <- na.omit(data.frame( x = classif[,pairs_s[1,x]], y = classif[,pairs_s[2,x]]))
                              xtab <- table(df[,1], df[,2])
                              xtab <- xtab[-which(rownames(xtab) == "O"),]
                              caret::confusionMatrix(xtab)
                              
                            } else {
                              
                              df <- na.omit(data.frame( x = classif[,pairs_s[1,x]], y = classif[,pairs_s[2,x]]))
                              xtab <- table(df[,1], df[,2])
                              caret::confusionMatrix(xtab)
                             
                              }

  
})#eo lapply


overall_accuracy <- sapply(1:length(overall_agreement), function(x) {
  
                      overall_agreement [[x]]$overall[1]

})#eo lapply


#accuracy per guild

accuracies_class <- lapply(1:length(overall_agreement), function(x) {
  
  if(dim(overall_agreement[[x]]$byClass)[1] == 4) {
    
    data.frame(H = overall_agreement[[x]]$byClass[1,11], I = overall_agreement[[x]]$byClass[2,11], O = NA, P = overall_agreement[[x]]$byClass[3,11], PK = overall_agreement[[x]]$byClass[4,11])
    
    
  } else if (dim(overall_agreement[[x]]$byClass)[2] == 4) {
    
    data.frame(H = overall_agreement[[x]]$byClass[1,11], I = overall_agreement[[x]]$byClass[2,11], O = NA, P = overall_agreement[[x]]$byClass[3,11], PK = overall_agreement[[x]]$byClass[4,11])
    
  } else {
    
    data.frame(H = overall_agreement[[x]]$byClass[1,11], I = overall_agreement[[x]]$byClass[2,11], O = overall_agreement[[x]]$byClass[3,11], P = overall_agreement[[x]]$byClass[4,11], PK = overall_agreement[[x]]$byClass[5,11])
    
  }
})#eo lapply



# ------------------------>>  plot overall accuracy and accuracy per guild

h = hist(overall_accuracy, plot=FALSE)

df <- data.frame(
  overall_accuracy = overall_accuracy,
  groups = cut(overall_accuracy, breaks = h$breaks)
) %>% group_by(groups) %>%
  dplyr::summarize(dens = n()/23)

gg_accuracies <- reshape::melt(accuracies_class)

p1 <- ggplot(df) +
  geom_bar(aes(x = groups, y = dens), stat = "identity", fill = "grey", color = "black",
           width = 1,
           position = position_nudge(x = 0.5)) +
  labs(y = "Experts' proportion", x = "Agreement") +
  scale_x_discrete(labels = c("0.6", "", "", "0.75", "", "", "0.9"), expand = c(0.01, 0.1)) +
  scale_y_continuous(expand = c(0, 0.005) )+
  theme_bw() + 
  geom_vline(xintercept = 4.6, linetype = 2, color = "red", size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none", panel.border = element_blank(), axis.line = element_line()) 
p1

p2 <- 
  ggplot(gg_accuracies) +
  geom_boxplot(aes(x = variable, y = value, fill = variable), alpha = 0.7, outlier.alpha = 0) +
  geom_jitter(aes(x = variable, y = value), width = 0.1, size = 0.8, alpha = 0.5) +
  labs(x = "Trophic guild", y = "Agreement") +
  scale_fill_fish_d(option = "Melichthys_vidua") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none", axis.line = element_line(), panel.border = element_blank()) 
p2

p <- cowplot::plot_grid(p1,p2, nrow = 1, labels = c("a)", "b)"), label_fontface = "plain", hjust = 0)
p 

ggsave(plot = p, "figures/Fig_1_Agreement.pdf", width = 10, height = 5)

}

## ----------------------->> levelplot of individual experts' pairs

make_plot2 <- function(classif){
  
  
  #generate experts pairs
  pairs <- combn(colnames(classif)[18:25], 2)
  
  #n species in common between pairs
  common_sp <- sapply(1:ncol(pairs), function (x) {
    
    combination <- pairs[,x] 
    common <- which(!is.na(classif[,combination[1]])) %in%  which(!is.na(classif[,combination[2]]))
    table(common)["TRUE"]
  })#eo sapply 
  

#only experts with more than 100 species in common
min_n_sp <- 200
pairs_s <- pairs[,which(common_sp > min_n_sp)]



dfs <- lapply(1:ncol(pairs_s), function (i) {
  
  pair  <- pairs_s[,i]
  
   na.omit(data.frame( x = classif[, pair[1]], y = classif[,pair[2]]))
  
}) #eo lapply

tabs <-lapply(1:ncol(pairs_s), function (i) {
  
  pair  <- pairs_s[,i]
  
  xtab <- table(dfs[[i]]$x, dfs[[i]]$y)
  t(apply(xtab, 1, function(x) x/colSums(xtab)))
  
}) #eo lapply 




p1 <- levelplot(tabs[[1]], col.regions = viridis(100), xlab=list("Mouillot et al. (2014)", cex=0.9), ylab = list("Morais & Bellwood (2018)", cex=0.75), scales=list(x=list(rot=90)))
p2 <- levelplot(tabs[[2]], col.regions = viridis(100), xlab=list("Mouillot et al. (2014)", cex=0.9), ylab = list("Brandl et al. (2016)", cex=0.75), scales=list(x=list(rot=90)))
p3 <- levelplot(tabs[[3]], col.regions = viridis(100), xlab=list("Mouillot et al. (2014)", cex=0.9), ylab = list("Halpern & Floeter (2008)", cex=0.75), scales=list(x=list(rot=90)))
p4 <- levelplot(tabs[[4]], col.regions = viridis(100), xlab=list("Mouillot et al. (2014)", cex=0.9), ylab = list("Yeager et al. (2017)", cex=0.75), scales=list(x=list(rot=90)))
p5 <- levelplot(tabs[[5]], col.regions = viridis(100), xlab=list("Mouillot et al. (2014)", cex=0.9), ylab = list("Stuart-Smith et al. (2018)", cex=0.75), scales=list(x=list(rot=90)))
p6 <- levelplot(tabs[[6]], col.regions = viridis(100), xlab=list("Morais & Bellwood (2018)", cex=0.9), ylab = list("Stuart-Smith et al. (2018)", cex=0.75), scales=list(x=list(rot=90)))
p7 <- levelplot(tabs[[7]], col.regions = viridis(100), xlab=list("Brandl et al. (2016)", cex=0.9), ylab = list("Stuart-Smith et al. (2018)", cex=0.75), scales=list(x=list(rot=90)))
p8 <- levelplot(tabs[[8]], col.regions = viridis(100), xlab=list("Halpern & Floeter (2008)", cex=0.9), ylab = list("Stuart-Smith et al. (2018)", cex=0.75), scales=list(x=list(rot=90)))
p9 <- levelplot(tabs[[9]], col.regions = viridis(100), xlab=list("Yeager et al. (2017)", cex=0.9), ylab = list("Stuart-Smith et al. (2018)", cex=0.75), scales=list(x=list(rot=90)))

plot_pairs <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9, ncol=3) 

postscript("figures/Fig_2_Agreement.eps", horizontal=T, width=10, height=10)

plot(plot_pairs)

dev.off()
}

