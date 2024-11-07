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

require(vegan)
source("scripts/analysis/biostats.R")

high_res <- read_csv('data/clean/high_resolution_summaries.csv')

combo_series <- 
  read_csv('data/clean/Series/combo_series.csv')

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# depth bins used for analysis (top level included)
bins <- c(0, 10, 50, 100, 200, 300, 400, 500, 2000)

# clustering --------------------------------------------------------------

# select only clustering variables (depth bins and st. deviations of vertical movement)
Ctad2 <- 
  high_res %>% 
  dplyr::select(d.b1:n.sd)

# Data transformations (Arcsine for proportions | Square-root for standard deviations)
Ctad2[,c(1:8,10:17)] <- Ctad2[,c(1:8,10:17)] / 100 # change to proportions
Ctad2[,c(1:8,10:17)] <- data.trans(Ctad2[,c(1:8,10:17)],method="asin",plot=F)
Ctad2[,c(9,18)] <- data.trans(Ctad2[,c(9,18)],method="power",exp=0.5,plot=F)

### PCA ###
#PCA with correlation matrix if variables are in different units or different scales
#  "" variance-covariance matrix if same units or data types
divecomp <- Ctad2
#Conduct PCA (w/ scaled and centered data)
dive.pca<-prcomp(divecomp, center = T, scale=T) 

#Monte Carlo randomization test for eigenvalue stat significance
#ordi.monte(divecomp,ord='pca') 

#Scores
dive.scores<-dive.pca$x 
dive.scores<-dive.scores[,1:5] # only use sig principal components as determined by Monte Carlo randomization test

#Calculate Euclidean distance matrix on centered/scaled PCA 1-3 scores
dive.scores.adj <- scale(dive.scores,center=T,scale=T)
site.eucd<-vegdist(dive.scores.adj,method='manhattan')

#Ward clustering (Murtagh & Legendre 2014 Journal of Classification 31: 274-295)
sitecl.ave<-hclust(site.eucd,method='ward.D2') 

#Mode of clusters was 5 (testing between 2 and 20, with options 'kl' thru 'gap' [20 indices])
#NbClust(dive.scores.adj,distance="manhattan",min.nc=2, max.nc=20, method="ward.D2",index="gap") 

#2:  3
#3:  4
#4:  2
#5:  7
#7:  1
#8:  2
#20: 1

Clust2 <- sitecl.ave

# Step 2: partition the data into its rightful cluster
combo_sub_grp2 <- dendextend::cutree(Clust2,
                                     k = 5,
                                     order_clusters_as_data = TRUE)
table(combo_sub_grp2)

# assign clusters to the tad input data
# select only clustering variables (depth bins and st. deviations of vertical movement)
Ctad2 <- 
  high_res %>% 
  dplyr::select(d.b1:n.sd)
Ctad2$cluster <- combo_sub_grp2

# assign clusters to high_resolution data
high_res <- high_res %>% 
  mutate(cluster = combo_sub_grp2)

clust_stamp2 <- high_res %>%
  dplyr::select(kode, species, cluster)



# cluster summaries -------------------------------------------------------

# 1. WHAT: what did each cluster look like --------------------------------------------------------------------

## percent time spent in each TAD bin ----
high_res %>% 
  # rename the clusters for easy reading
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>%
# calculate the average % time spent in each depth bin
  group_by(cluster) %>% 
  summarise(
    # daytime bins
    'day.0-10' = mean(d.b1),
    # d.b1.sd = sd(d.b1),
    'day.10-50' = mean(d.b2),
    # d.b2.sd = sd(d.b2),
    'day.50-100' = mean(d.b3),
    # d.b3.sd = sd(d.b3),
    'day.100-200' = mean(d.b4),
    # d.b4.sd = sd(d.b4),
    'day.200-300' = mean(d.b5),
    # d.b5.sd = sd(d.b5),
    'day.300-400' = mean(d.b6),
    # d.b6.sd = sd(d.b6),
    'day.400-500' = mean(d.b7),
    # d.b7.sd = sd(d.b7),
    'day.500-2000' = mean(d.b8),
    d.b8.sd = sd(d.b8),
    d.sd = mean(d.sd),
    
    # nighttime bins
    'night.0-10' = mean(n.b1),
    # n.b1.sd = sd(n.b1),
    'night.10-50' = mean(n.b2),
    # n.b2.sd = sd(n.b2),
    'night.50-100' = mean(n.b3),
    # n.b3.sd = sd(n.b3),
    'night.100-200' = mean(n.b4),
    n.b4.sd = sd(n.b4),
    'night.200-300' = mean(n.b5),
    # n.b5.sd = sd(n.b5),
    'night.300-400' = mean(n.b6),
    # n.b6.sd = sd(n.b6),
    'night.400-500' = mean(n.b7),
    # n.b7.sd = sd(n.b7),
    'night.500-2000' = mean(n.b8),
    # n.b8.sd = sd(n.b8),
    n.sd = mean(n.sd)) %>% View() 

## calculate the percentage of time spent in the epi/mesopelagic ----

high_res %>% 
  # rename the clusters for easy reading
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>%
  # create columns for % time in each zone
  group_by(kode, cluster) %>%
  summarize(
    # daytime depth use
    'day_epi' = sum(c(d.b1, d.b2, d.b3, d.b4)),
    'day_meso' = sum(c(d.b5, d.b6, d.b7, d.b8)),
    'night_epi' = sum(c(n.b1, n.b2, n.b3, n.b4)),
    'night_meso' = sum(c(n.b5, n.b6, n.b7, n.b8))) %>% 
  # calculate the average % time spent in each depth bin
  group_by(cluster) %>% 
  summarise(
    # daytime habitat
    # 'day_mean_epi' = mean(day_epi),
    # 'day_sd_epi' = sd(day_epi),
    'day_mean_meso' = mean(day_meso),
    'day_sd_meso' = sd(day_meso),
    # nightly habitat
    # 'night_mean_epi' = mean(night_epi),
    # 'night_sd_epi' = sd(night_epi),
    'night_mean_meso' = mean(night_meso),
    'night_sd_meso' = sd(night_meso),
    'overall_mean_meso' = mean(c(day_meso, night_meso)),
    'overall_sd_meso' = sd(c(day_meso, night_meso))) %>% View()
    
## exact average depth targeted ----

### read in fine-scale data ----

combo_series %>%
  left_join(clust_stamp2) %>% 
  mutate(
    cluster = case_when(
      cluster == 1 ~ 'EPI 2',
      cluster == 2 ~ 'DVM 1',
      cluster == 3 ~ 'EPI 1',
      cluster == 4 ~ 'DVM 2',
      cluster == 5 ~ 'DVM 3')) %>%
  # remove 1207 incomplete / continental shelf observations 
  filter(!is.na(cluster)) %>% # pull(kode) %>% unique() %>% length()
  # remove any observations not classified as d or n (281 time recs.)
  filter(!is.na(dn)) %>% 
  group_by(cluster, dn) %>% 
  summarize(
    Med.depth = median(depth))

## 2. WHEN: when do the clusters occur in time ----

high_res %>%
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>% 
  rbind((high_res %>% 
           filter(cluster <= 5) %>% 
           mutate(cluster = 0))) %>% 
  mutate(yday = lubridate::yday(Date)) %>% 
  group_by(cluster) %>% 
  summarise(
    Avg.date. = mean(yday),
    SD.yday = sd(yday))
         
         
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
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>%
  group_by(species, ptt, cluster) %>% 
  summarize(count = n()) %>% #pivot_wider(names_from = cluster, values_from = count)
  left_join(
    high_res %>% 
      filter(cluster <=5) %>% 
      mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>%
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
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>%
  group_by(cluster, dn) %>%
  summarize(
    min = min(depth),
    q1 = quantile(depth, 0.25),
    median = median(depth),
    mean = mean(depth),
    q3 = quantile(depth, 0.75),
    max = max(depth)) %>% 
  filter(dn == 'd')

combo_series %>%
  left_join(clust_stamp2, by = 'kode') %>% 
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>%
  group_by(cluster, dn) %>%
  summarize(
    min = min(depth),
    q1 = quantile(depth, 0.25),
    median = median(depth),
    mean = mean(depth),
    q3 = quantile(depth, 0.75),
    max = max(depth)) %>% 
  filter(dn == 'n')

# TAD summaries

high_res %>% 
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>%
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
            d.sd = mean(d.sd),
            night.epi = mean(n.epi),
            n.epi.sd = sd(n.epi),
            night.meso = mean(n.meso),
            n.meso.sd = sd(n.meso),
            n.sd = mean(n.sd),)

# dive profiles

# environmental conditions within each cluster
high_res %>%
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3'),
    cluster = factor(cluster, levels = c(
      'EPI 1', 
      'EPI 2', 
      'DVM 1', 
      'DVM 2',
      'DVM 3'),
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
  geom_raster(aes(x,y,fill = z)) +
  geom_sf(data = world, fill = 'black') + 
  geom_contour(
    aes(x = x,
        y = y,
        z = z),
    color = "black",
    breaks = c(-1000)) +
  geom_point(data = high_res %>% 
               filter(cluster <= 5) %>% 
               mutate(cluster = case_when(
                 cluster == 1 ~ 'EPI 2',
                 cluster == 2 ~ 'DVM 1',
                 cluster == 3 ~ 'EPI 1',
                 cluster == 4 ~ 'DVM 2',
                 cluster == 5 ~ 'DVM 3')) %>%
               mutate(
                 cluster = factor(cluster,
                                  levels = c(
                                    'EPI 1', 
                                    'EPI 2', 
                                    'DVM 1', 
                                    'DVM 2',
                                    'DVM 3'),
                                  ordered = T)),
             aes(x = x, y = y, color = cluster), size = 0.5) +
  scale_color_manual(values = c("yellow",
                                "yellow3",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  scale_x_continuous(breaks = c(-80, -70, -60, -50, -40)) +
  scale_y_continuous(breaks = c(10, 20, 30, 40)) +
  # facet_wrap(~cluster, nrow = 2) +
  coord_sf(xlim = c(-80, -35), ylim = c(9, 45)) +
  labs(x = "Longitude",
       y = 'Latitude',
       color = 'Cluster') 

# latitude distribution of observations from each cluster
high_res %>%
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3'),
    cluster = factor(cluster, levels = c(
      'EPI 1', 
      'EPI 2', 
      'DVM 1', 
      'DVM 2',
      'DVM 3'),
      ordered = T)) %>% 
  group_by(cluster
           #ptt
           ) %>%
  summarize(
    Q1.Lat. = quantile(latitude, 0.05),
    Avg.Lat = mean(latitude),
    SE.Lat = std.error(latitude),
    Q3.Lat = quantile(latitude, 0.95)
    # Q1.Lon. = quantile(longitude, 0.05),
    # Avg.Long = mean(longitude),
    # SD.Long = sd(longitude)
    # Q3.Lon = quantile(longitude, 0.90)
    )
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
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>% 
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
  summarize(count = n())  #%>%  View()



# figure 4 ----------------------------------------------------------------

cluster_series <- 
  combo_series %>% 
  inner_join(clust_stamp2,
             by = 'kode') %>% 
  filter(cluster <= 5) %>% 
  mutate(ogCluster = cluster,
         cluster = case_when(
           cluster == 1 ~ 'EPI 2',
           cluster == 2 ~ 'DVM 1',
           cluster == 3 ~ 'EPI 1',
           cluster == 4 ~ 'DVM 2',
           cluster == 5 ~ 'DVM 3'),
         cluster = factor(cluster, levels = c('EPI 1', 'EPI 2', 'DVM 1', 'DVM 2', 'DVM 3'), ordered = T))
require(viridis)
ggplot(data = cluster_series, aes(x = Local_Time, y = depth) ) +
  stat_bin_2d(aes(fill = after_stat(ndensity)), geom = "tile", bins = c(100, 100)) +
  scale_fill_gradientn(colors = viridis(100), limits = c(0,0.5), oob = scales::squish) + 
  labs(fill = "Relative Density (log)") +
  scale_y_reverse() +
  scale_x_continuous(breaks = c(21600, 64800), labels = c("06:00", "18:00")) +
  xlab("Local Time") +
  ylab("Depth (m)") +
  facet_grid(species.x~cluster) +
  theme_linedraw() +
  theme(legend.position = 'bottom')
