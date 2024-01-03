
# setup -------------------------------------------------------------------

library(tidyverse)
library(tags2etuff)
library(NbClust)
library(dendextend)
library(hms)
library(vegan)
library(ggdendro)
# mapping
library(rnaturalearthdata)
library(rgeos)
library(ggspatial)
library(sf)

high_res <- read_csv('data/clean/high_resolution_summaries.csv')

combo_series <- 
  read_csv('data/clean/Series/combo_series.csv')

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

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


# removing outlier clusters -----------------------------------------------

high_res_prime_clusts <- high_res %>% filter(cluster <= 5)

Tad_prime <- 
  high_res_prime_clusts %>% 
  dplyr::select(d.b1:n.sd)

## Heirarchical Clustering

# Step 1: make the clusters
# create a distance matrix
man.dist.prime <- dist(Tad_prime, method = "manhattan")
# perform cluster assignment
clusters.prime <- hclust(man.dist.prime, method = "average")

# Step 2: partition the data into its rightful cluster
cut.clusters.prime <- dendextend::cutree(clusters.prime,
                                     k = 8,
                                     order_clusters_as_data = TRUE)
table(cut.clusters.prime)

# testing most parsimonious cluster number ---
## CH INDEX:
NbClust(data = Tad_prime, diss = man.dist.prime, distance = NULL, min.nc = 2,
        max.nc = 20, method = "average", index = "ch")
# Raw: 2

## DUDA TEST: Smallest number of clusters which produces a critical value greater than the Duda score
NbClust(data = Tad_prime, diss = man.dist.prime, distance = NULL, min.nc = 2,
        max.nc = 20, method = "average", index = "duda")
# Raw: 8 clusters 

## C INDEX
NbClust(data = Tad_prime, diss = man.dist.prime, distance = NULL, min.nc = 2,
        max.nc = 20, method = "average", index = "cindex")
# Raw: 6 clusters 

## BEALE INDEX
NbClust(data = Tad_prime, diss = man.dist.prime, distance = NULL, min.nc = 2,
        max.nc = 20, method = "average", index = "beale")
# Raw: 8 clusters 

# cluster summaries -------------------------------------------------------

# how much time do sharks spend across depth bins within each cluster?
bins <- c(0, 10, 50, 100, 200, 300, 400, 500, 2000)

## 1. WHO: what was the frequency of clusters between species? ----

high_res %>% 
  filter(cluster <= 5) %>% 
  group_by(species, cluster) %>% 
  summarize(n = n())

# what was the frequency of clusters among individuals
high_res %>% 
  filter(cluster <=5) %>% 
  group_by(species, cluster) %>% 
  summarize(uniqueIndividuals = n_distinct(ptt))

# are a few individuals driving specific clusters
high_res %>% 
  filter(cluster <=5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'DVM 1',
    cluster == 2 ~ 'Epipelagic',
    cluster == 3 ~ 'DVM 2',
    cluster == 4 ~ 'DVM 3',
    cluster == 5 ~ 'DVM 4')) %>%
  group_by(species, ptt, cluster) %>% 
  summarize(count = n()) %>% #pivot_wider(names_from = cluster, values_from = count)
  left_join(
    high_res %>% 
              filter(cluster <=5) %>% 
              mutate(cluster = case_when(
                cluster == 1 ~ 'DVM 1',
                cluster == 2 ~ 'Epipelagic',
                cluster == 3 ~ 'DVM 2',
                cluster == 4 ~ 'DVM 3',
                cluster == 5 ~ 'DVM 4')) %>%
              group_by(species, ptt) %>% 
              summarize(total = n())) %>% # group_by(species, ptt) %>% summarize(sum = sum(count), total = mean(total))
  mutate(percObs = count / total) %>% 
  group_by(species, cluster) %>% 
  summarize(avgPercObs = mean(percObs), sd = sd(percObs))

## 2. WHAT: what did each cluster look like ----

combo_series %>%
  left_join(clust_stamp2, by = 'kode') %>% 
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'DVM 1',
    cluster == 2 ~ 'Epipelagic',
    cluster == 3 ~ 'DVM 2',
    cluster == 4 ~ 'DVM 3',
    cluster == 5 ~ 'DVM 4')) %>%
  group_by(cluster, dn) %>%
  summarize(
    min = min(depth),
    q1 = quantile(depth, 0.25),
    median = median(depth),
    mean = mean(depth),
    q3 = quantile(depth, 0.75),
    max = max(depth)) %>% 
  filter(dn == 'd')

# TAD summaries

high_res %>% 
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'DVM 1',
    cluster == 2 ~ 'Epipelagic',
    cluster == 3 ~ 'DVM 2',
    cluster == 4 ~ 'DVM 3',
    cluster == 5 ~ 'DVM 4')) %>%
  group_by(ptt, cluster) %>% 
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
    'd.epi' = sum(c(`day.0-10`, `day.10-50`, `day.50-100`, `day.100-200`)),
    'd.meso' = sum(c(`day.200-300`, `day.300-400`, `day.400-500`, `day.500-2000`)),
    
    # nighttime bins
    'night.0-10' = mean(n.b1),
    'night.10-50' = mean(n.b2),
    'night.50-100' = mean(n.b3),
    'night.100-200' = mean(n.b4),
    'night.200-300' = mean(n.b5),
    'night.300-400' = mean(n.b6),
    'night.400-500' = mean(n.b7),
    'night.500-2000' = mean(n.b8),
    n.sd = mean(n.sd),
    'n.epi' = sum(c(`night.0-10`, `night.10-50`, `night.50-100`, `night.100-200`)),
    'n.meso' = sum(c(`night.200-300`, `night.300-400`, `night.400-500`, `night.500-2000`))) %>% #View() 
  # average mesopelagic and epipelagic occupancy by cluster
  group_by(cluster) %>% 
  summarise(day.epi = mean(d.epi),
            d.epi.sd = sd(d.epi),
            day.meso = mean(d.meso),
            d.meso.sd = sd(d.meso),
            night.epi = mean(n.epi),
            n.epi.sd = sd(n.epi),
            night.meso = mean(n.meso),
            n.meso.sd = sd(n.meso))

  # dive profiles

# environmental conditions within each cluster
high_res %>%
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'DVM 1',
    cluster == 2 ~ 'Epipelagic',
    cluster == 3 ~ 'DVM 2',
    cluster == 4 ~ 'DVM 3',
    cluster == 5 ~ 'DVM 4'),
    cluster = factor(cluster, levels = c(
      'Epipelagic', 
      'DVM 1', 
      'DVM 2', 
      'DVM 3',
      'DVM 4'),
      ordered = T
    )) %>% 
  group_by(cluster) %>% 
  # summary of ild.5
  summarize(
    min = min(ild.5, na.rm = T),
    q1 = quantile(ild.5, 0.25),
    median = median(ild.5, na.rm = T),
    mean = mean(ild.5, na.rm = T),
    q3 = quantile(ild.5, 0.75),
    max = max(ild.5, na.rm = T)
  )

## 3. WHERE: average position of clusters ----

# plot cluster locations

r <- raster::raster('data/raw/global_bathy_0.01.nc')
## this is a global grid with pacific-centered coordinates (longitudes 0 to 360)
## raster::rotate converts from 0-360 longitudes to 180 longitudes (atlantic-centered)
r <- raster::rotate(r)
r <- raster::aggregate(r, fact = 5)
r_df <- raster::as.data.frame(r, xy = TRUE)
r_df <- r_df %>% filter(z < 0)

ggplot(data = r_df) +
  geom_sf(data = world, fill = 'black') + 
  geom_point(data = high_res %>% 
               filter(cluster <= 5) %>% 
               mutate(cluster = case_when(
                 cluster == 1 ~ 'DVM 1',
                 cluster == 2 ~ 'Epipelagic',
                 cluster == 3 ~ 'DVM 2',
                 cluster == 4 ~ 'DVM 3',
                 cluster == 5 ~ 'DVM 4')) %>% 
               mutate(
                 cluster = factor(cluster,
                                   levels = c('Epipelagic',
                                              'DVM 1',
                                              'DVM 2',
                                              'DVM 3' ,
                                              'DVM 4',
                                              'DVM 5'),
                                   ordered = T)),
             aes(x = x, y = y, color = cluster)) +
  scale_color_manual(values = c("#FFFF5CFF",
                                "#78CEA3FF",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  scale_x_continuous(breaks = c(-80, -70, -60, -50, -40)) +
  scale_y_continuous(breaks = c(10, 20, 30, 40)) +
  coord_sf(xlim = c(-80, -35), ylim = c(9, 45)) +
  labs(x = "Longitude",
       y = 'Latitude',
       color = 'Cluster') 

# latitude distribution of observations from each cluster
high_res %>%
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'DVM 1',
    cluster == 2 ~ 'Epipelagic',
    cluster == 3 ~ 'DVM 2',
    cluster == 4 ~ 'DVM 3',
    cluster == 5 ~ 'DVM 4'),
    cluster = factor(cluster, levels = c(
      'Epipelagic', 
      'DVM 1', 
      'DVM 2', 
      'DVM 3',
      'DVM 4'),
      ordered = T
      )) %>% 
  ggplot() +
  geom_boxplot(aes(x = cluster, y = y, fill = cluster)) +
  scale_fill_manual(values = c("#FFFF5CFF",
                               "#78CEA3FF",
                               "#488E9EFF",
                               "#404C8BFF",
                               "#281A2CFF")) + 
  labs(x = NULL, y = 'Latitude') + 
  guides(fill = 'none')+
  theme_classic()

## 4. WHEN: when do the clusters occur in time ----

high_res %>%
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'DVM 1',
    cluster == 2 ~ 'Epipelagic',
    cluster == 3 ~ 'DVM 2',
    cluster == 4 ~ 'DVM 3',
    cluster == 5 ~ 'DVM 4')) %>% 
  rbind((high_res %>% 
           filter(cluster <= 5) %>% 
           mutate(cluster = 0))) %>% 
  mutate(yday = lubridate::yday(Date),
         # add a correction to center yday on tag date
         yday = if_else(
           yday >=238,
           (yday - 238),
           (yday + 127))) %>%  # pull(yday) %>% summary()
  mutate(season = 
           case_when(
             yday < 27 ~ 'Summer',
             yday >= 27 & yday < 117 ~ 'Fall',
             yday >= 117 & yday < 205 ~ 'Winter',
             yday >= 205 & yday < 298 ~ 'Spring',
             yday >= 298 ~ 'Summer')) %>% 
  group_by(cluster, season) %>% 
  summarize(count = n()) # %>%  View()
  summarize(
    min = min(yday),
    q1 = quantile(yday, 0.25),
    median = median(yday),
    mean = mean(yday),
    q3 = quantile(yday, 0.75),
    max = max(yday))
