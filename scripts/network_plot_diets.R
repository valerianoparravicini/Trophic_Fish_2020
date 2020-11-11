library(plyr)
library(dplyr)
library(tidyverse)
library(vegan)
library(nlme)
library(lme4)
library(ggpubr)
library(ggrepel)
library(colorRamps)
library(tidybayes)
library(brms)
library(fishualize)
library(modelr)
library(beepr)
library(png)
library(grid)
library(scales)
library(iNEXT)
library(geomnet)
library(GGally)
library(extrafont)


theme_sjb <- function (base_size = 8, base_family = "") {
  theme_bw() %+replace% 
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 8, family = "Arial"),
          axis.title = element_text(color = "black", size = 8, family = "Arial"),
          legend.title=element_text(size=8), 
          legend.text=element_text(size=8))
}


### FROM VALE'S SCRIPT ###########################################################################
diets <- read.csv("data/data_ISfull_grp_final_no_sp_no_wormy.csv", sep=",", dec=".", row.names=1, na = "NA")

diets_p <- reshape::cast(diets, fish_sp ~ grp6, value = "quantity", fun.aggregate=sum)
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
#############################################################################################
mod_predators <- read.csv(file = "results/mod_predators.csv")
mod_prey_items <- read.csv(file = "results/mod_prey_items.csv")

fish.clusters <- as.data.frame(mod_predators) %>%
  rename(genspe = V1, fishcluster = V2)

prey.clusters <- as.data.frame(mod_prey_items) %>%
  rename(prey = V1, preycluster = V2)


diet.matrix <- diets_p %>%
  mutate(genspe = rownames(.)) %>%
  select(genspe, everything())


# meta.23s <- sub23s.loc[c(1:3)]
diet.network <- diet.matrix %>%
  pivot_longer(cols = 2:37, names_to = "prey") %>%
  filter(value > 0) %>%
  left_join(fish.clusters, by = "genspe") %>%
  full_join(prey.clusters, by = "prey") %>%
  group_by(fishcluster, prey, preycluster) %>%
  summarize(sum.prop = mean(value)) %>%
  ungroup() %>%
  mutate(fishcluster = as.factor(fishcluster)) %>%
  filter(sum.prop > 0.0)



trophic.network <- ggplot(data = diet.network) +
  geom_net(layout.alg = "circle", aes(from_id = prey, to_id = fishcluster,  linewidth = sum.prop*4), 
           labelon = T, repel = TRUE, curvature = 0,directed = F, arrowgap = 0.001, arrowsize = 0.01) +
  theme_net() +
  theme(legend.position = "none")
trophic.network

ggsave("plots/networkplot_raw.pdf", trophic.network, width = 15, height = 15, useDingbats = F)



families <- read.csv(file = "data/diet_families.csv") %>%
  right_join(fish.clusters) %>%
  group_by(family, fishcluster) %>%
  summarize(total_n = n()) %>%
  ungroup() %>%
  group_by(fishcluster) %>%
  mutate(percent_fam = total_n/sum(total_n)) %>%
  filter(percent_fam > 0.02) %>%
  drop_na()


c = do.call(colorRampPalette(c("darkblue", "cyan")),list(10))


clus1 <- families %>%
  filter(fishcluster == 1) %>%
  mutate(colorvector = do.call(colorRampPalette(c("deeppink4", "lightpink")),list(length(family)))) 
clus1plot <- ggplot(clus1, aes(x = "", y = percent_fam, fill = family)) +
  geom_bar(width = 1, stat = "identity", color = "white", size = 0.2) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = clus1$colorvector) +
  theme_sjb() + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.line = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.ticks = element_blank(),
                     legend.position = "right")
clus1plot

clus2 <- families %>%
  filter(fishcluster == 2) %>%
  mutate(colorvector = do.call(colorRampPalette(c("darkgreen", "chartreuse")),list(length(family))))

 

clus2plot <-  ggplot(clus2, aes(x = "", y = percent_fam, fill = family)) +
  geom_bar(width = 1, stat = "identity", color = "white", size = 0.2) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = clus2$colorvector) +
  theme_sjb() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.line = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.ticks = element_blank(),
                      legend.position = "right")
clus2plot


clus3 <- families %>%
  filter(fishcluster == 3) %>%
  mutate(colorvector = do.call(colorRampPalette(c("darkblue", "lightblue1")),list(length(family))))

clus3plot <-  ggplot(clus3, aes(x = "", y = percent_fam, fill = family)) +
  geom_bar(width = 1, stat = "identity", color = "white", size = 0.2) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = clus3$colorvector) +
  theme_sjb() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.line = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.ticks = element_blank(),
                      legend.position = "right")
clus3plot

clus4 <- families %>%
  filter(fishcluster == 4) %>%
  mutate(colorvector = do.call(colorRampPalette(c("darkred", "salmon1")),list(length(family))))

  

clus4plot <-  ggplot(clus4, aes(x = "", y = percent_fam, fill = family)) +
  geom_bar(width = 1, stat = "identity", color = "white", size = 0.2) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = clus4$colorvector) +
  theme_sjb() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.line = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.ticks = element_blank(),
                      legend.position = "right")
clus4plot

clus5 <- families %>%
  filter(fishcluster == 5) %>%
  mutate(colorvector = do.call(colorRampPalette(c("turquoise4", "paleturquoise1")),list(length(family))))

clus5plot <-  ggplot(clus5, aes(x = "", y = percent_fam, fill = family)) +
  geom_bar(width = 1, stat = "identity", color = "white", size = 0.2) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = clus5$colorvector) +
  theme_sjb() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.line = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.ticks = element_blank(),
                      legend.position = "right")
clus5plot


clus6 <- families %>%
  filter(fishcluster == 6) %>%
  mutate(colorvector = do.call(colorRampPalette(c("saddlebrown", "burlywood1")),list(length(family))))

clus6plot <- ggplot(clus6, aes(x = "", y = percent_fam, fill = family)) +
  geom_bar(width = 1, stat = "identity", color = "white", size = 0.2) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = clus6$colorvector) +
  theme_sjb() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.line = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.ticks = element_blank(),
                      legend.position = "right")
clus6plot

clus7 <- families %>%
  filter(fishcluster == 7) %>%
  mutate(colorvector = do.call(colorRampPalette(c("purple4", "lavender")),list(length(family))))

  

  clus7plot <-  ggplot(clus7, aes(x = "", y = percent_fam, fill = family)) +
  geom_bar(width = 1, stat = "identity", color = "seashell", size = 0.2) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = clus7$colorvector) +
  theme_sjb() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.line = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.ticks = element_blank(),
                      legend.position = "right")
clus7plot


clus8 <- families %>%
  filter(fishcluster == 8) %>%
  mutate(colorvector = do.call(colorRampPalette(c("darkorange3", "yellow1")),list(length(family))))

clus8plot <-  ggplot(clus8, aes(x = "", y = percent_fam, fill = family)) +
  geom_bar(width = 1, stat = "identity", color = "white", size = 0.2) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = clus8$colorvector) +
  theme_sjb() + theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.line = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.ticks = element_blank(),
                      legend.position = "right")
clus8plot


piecharts = ggarrange(clus1plot, clus2plot, clus3plot, clus4plot, 
                      clus5plot, clus6plot, clus7plot, clus8plot,
                      nrow = 2, ncol = 4, legend = "none")

piecharts

ggsave("plots/piecharts.pdf", piecharts, width = 10, height = 10, useDingbats = F)

piechartlegends = ggarrange(clus1plot, clus2plot, clus3plot, clus4plot, 
                      clus5plot, clus6plot, clus7plot, clus8plot,
                      nrow = 2, ncol = 4, legend = "right")

ggsave("plots/piechartslegends.pdf", piechartlegends, width = 10, height = 10, useDingbats = F)
