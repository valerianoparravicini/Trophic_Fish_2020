#### This script recreates the figure of the phylogenetic tree with probabilities of all trophic guilds

library(ape)
library(adephylo)
library(phylobase)
library(phylosignal)
library(ggtree)
library(phytools)
library(fishtree)
library(tidytree)
library(tidyverse)
library(ggnewscale)

tax <- rfishbase::load_taxa()
dietcat <- read.csv("data/dietcat.csv")

tree <- fishtree::fishtree_complete_phylogeny(dietcat$species)
tree1 <- tree[[1]]

ppred <- read.csv("data/ppred.csv")
rownames(ppred) <- ppred$species
ppred

tree1 <- fishtree::fishtree_phylogeny(species = ppred$species)

tree1t <- as_tibble(tree1)
tax$label <- gsub(" ", "_", tax$Species)
tree1t <- left_join(tree1t, tax)
tree1t[!is.na(tree1t$Family) & tree1t$Family == "Scaridae", "Family"] <- "Labridae"

tree1d <- as.treedata(tree1t)


data <- filter(dietcat, phylo%in% tree1$tip.label)
ppred <- filter(ppred, species%in% tree1$tip.label)
rownames(ppred) <- ppred$species


p <- ggtree(tree1d, layout = "circular")


p <- 
  gheatmap(p, ppred[,1, drop = FALSE], offset = 0.2, width=0.1, 
           colnames_position = "top", font.size=2, color="white") + 
  scale_fill_gradient(low = "white", high = "#BC5380")

p <- p + new_scale_fill()


p2 <- gheatmap(p, ppred[,2, drop = FALSE], offset = 15, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#5CD400")

p2 <- p2 + new_scale_fill()

p3 <- gheatmap(p2, ppred[,3, drop = FALSE], offset = 30, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#1B229B")
p3 <- p3 + new_scale_fill()

p4 <- gheatmap(p3, ppred[,4, drop = FALSE], offset = 45, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#EE3B00")
p4 <- p4 + new_scale_fill()

p5 <- gheatmap(p4, ppred[,5, drop = FALSE], offset = 60, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#46B3B6") +
  theme(legend.position = "none")

p5 <- p5 + new_scale_fill()

p6 <- gheatmap(p5, ppred[,6, drop = FALSE], offset = 75, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#E2AF79") +
  theme(legend.position = "none")

p6 <- p6 + new_scale_fill()

p7 <- gheatmap(p6, ppred[,7, drop = FALSE], offset = 90, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#855EB0") +
  theme(legend.position = "none")

p7 <- p7 + new_scale_fill()

p8 <- gheatmap(p7, ppred[,8, drop = FALSE], offset = 105, width = 0.1,
               colnames_position = "top", font.size = 2, color = "white") +
  scale_fill_gradient(low = "white", high = "#F3DC00") +
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
  scale_fill_gradient(low = "lightgrey", high = "black") +
  theme(legend.position = "none")


ggsave("output/plots/diet_phyloplot.pdf", p9, height = 12, width = 12)


p10 <- gheatmap(p8, ppred[, 10, drop = FALSE], offset = 125, width = 0.05,
                colnames_position = "top", font.size = 2, color = NULL) +
  scale_fill_gradient(low = "lightgrey", high = "black") 
ggsave("output/plots/diet_phyloplot_legend.pdf",p10, height = 18, width = 12)


###### species label #####

p + geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 5)
ggsave('output/plots/treefamily.pdf', height = 12, width = 12)

p9 + geom_tiplab(aes(angle = angle), size = 2, offset = 145)
ggsave("output/plots/diet_phyloplot_species.pdf", height = 12, width = 12)

