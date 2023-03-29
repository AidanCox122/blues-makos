
# setup -------------------------------------------------------------------

library(tidyverse)
library(tags2etuff)
library(NbClust)
library(dendextend)
library(hms)
library(vegan)
library(ggdendro)

high_res <- read_csv('data/clean/high_resolution_summaries.csv')


# clustering --------------------------------------------------------------

# select only clustering variables (depth bins and st. deviations of vertical movement)
Ctad2 <- 
  high_res %>% 
  dplyr::select(d.b1:n.sd)

## Heirarchical Clustering

# Step 1: make the clusters
# create a distance matrix
combo_dist2 <- dist(Ctad2, method = "manhattan")
# perform cluster assignment
Clust2 <- hclust(combo_dist2, method = "average")

# Step 2: partition the data into its rightful cluster
combo_sub_grp2 <- dendextend::cutree(Clust2,
                                     k = 12,
                                     order_clusters_as_data = TRUE)
table(combo_sub_grp2)

# assign clusters to the tad input data
Ctad2$cluster <- combo_sub_grp2

# assign clusters to high_resolution data
high_res <- high_res %>% 
  mutate(cluster = combo_sub_grp2)

clust_stamp2 <- high_res %>%
  dplyr::select(kode, species, cluster)


# cluster summaries -------------------------------------------------------

# how much time do sharks spend across depth bins within each cluster?
bins <- c(0, 10, 50, 100, 200, 300, 400, 500, 2000)


high_res %>% 
  group_by(cluster) %>% 
  summarise(
    # daytime bins
    'day.0-10' = mean(d.b1),
    'day.10-50' = mean(d.b2),
    'day.50-100' = mean(d.b3),
    'day.100-200' = mean(d.b4),
    'day.200-300' = mean(d.b5),
    'day.300-400' = mean(d.b6),
    'day.400-500' = mean(d.b7),
    'day.500-2000' = mean(d.b8),
    d.sd = mean(d.sd),
    
    # nighttime bins
    'night.0-10' = mean(n.b1),
    'night.10-50' = mean(n.b2),
    'night.50-100' = mean(n.b3),
    'night.100-200' = mean(n.b4),
    'night.200-300' = mean(n.b5),
    'night.300-400' = mean(n.b6),
    'night.400-500' = mean(n.b7),
    'night.500-2000' = mean(n.b8),
    n.sd = mean(n.sd)) %>% View()

  # how much time did sharks spend in the mesopelagic across clusters?
  group_by(cluster) %>% 
  summarise(
    # daytime bins
    'd.epi' = sum(c(`day.0-10`, `day.10-50`, `day.50-100`, `day.100-200`)),
    'd.meso' = sum(c(`day.200-300`, `day.300-400`, `day.400-500`, `day.500-2000`)),
    
    # nighttime bins
    'n.epi' = sum(c(`night.0-10`, `night.10-50`, `night.50-100`, `night.100-200`)),
    'n.meso' = sum(c(`night.200-300`, `night.300-400`, `night.400-500`, `night.500-2000`)))

  
# what was the frequency of clusters between species?

high_res %>% 
  group_by(species, cluster) %>% 
  summarize(n = n())

# environmental conditions within each cluster
high_res %>% 
  filter(species == 'I.oxyrinchus' & cluster == 5) %>% 
  pull(ssh) %>% 
  summary()

high_res %>% 
  filter(species == 'P.glauca' & cluster == 5) %>% 
  pull(ssh) %>% 
  summary()

