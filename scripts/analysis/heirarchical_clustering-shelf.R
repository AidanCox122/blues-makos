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

high_res_all <- read_csv('/Users/aidansmacpro/Desktop/blues-makos/data/clean/high_resolution_all_obs.csv')

combo_series <- 
  read_csv('data/clean/Series/combo_series.csv')

combo_track <- 
  read_csv('data/clean/Tracks/combo_track.csv')

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# depth bins used for analysis (top level included)
bins <- c(0, 10, 50, 100, 200, 300, 400, 500, 2000)

# cluster with all observations -------------------------------------------

Ctad2_test <- 
  high_res_all %>% 
  dplyr::select(d.b1:d.sd, n.b1:n.sd)

# Data transformations (Arcsine for proportions | Square-root for standard deviations)
Ctad2_test[,c(1:8,10:17)] <- Ctad2_test[,c(1:8,10:17)] / 100 # change to proportions
Ctad2_test[,c(1:8,10:17)] <- data.trans(Ctad2_test[,c(1:8,10:17)],method="asin",plot=F)
Ctad2_test[,c(9,18)] <- data.trans(Ctad2_test[,c(9,18)],method="power",exp=0.5,plot=F)

### PCA ###
#PCA with correlation matrix if variables are in different units or different scales
#  "" variance-covariance matrix if same units or data types
divecomp_test <- Ctad2_test
#Conduct PCA (w/ scaled and centered data)
dive.pca_test<-prcomp(divecomp_test, center = T, scale=T) 

#Monte Carlo randomization test for eigenvalue stat significance
#ordi.monte(divecomp,ord='pca') 

#Scores
dive.scores_test<-dive.pca_test$x 
dive.scores_test<-dive.scores_test[,1:5] # only use sig principal components as determined by Monte Carlo randomization test

#Calculate Euclidean distance matrix on centered/scaled PCA 1-3 scores
dive.scores.adj_test <- scale(dive.scores_test,center=T,scale=T)
site.eucd_test<-vegdist(dive.scores.adj_test,method='manhattan')

#Ward clustering (Murtagh & Legendre 2014 Journal of Classification 31: 274-295)
sitecl.ave_test<-hclust(site.eucd_test,method='ward.D2') 

#Mode of clusters was 5 (testing between 2 and 20, with options 'kl' thru 'gap' [20 indices])
# NbClust(dive.scores.adj_test,distance="manhattan",min.nc=2, max.nc=20, method="ward.D2",index="ptbiserial") 

Clust2_test <- sitecl.ave_test

# Step 2: partition the data into its rightful cluster
combo_sub_grp2_test <- dendextend::cutree(Clust2_test,
                                     k = 6,
                                     order_clusters_as_data = TRUE)
table(combo_sub_grp2_test)

high_res_all <- 
  high_res_all %>% 
  mutate(cluster = combo_sub_grp2_test)

clust_stamp_test <- high_res_all %>%
  dplyr::select(kode, species, cluster)

## Visualize vertical movement in clusters ----
cluster_series_test <- 
  combo_series %>% 
  inner_join(clust_stamp_test,
             by = 'kode') %>% # summary()
  mutate(ogCluster = cluster,
         cluster = case_when(
           cluster == 1 ~ 'EPI 2',
           cluster == 2 ~ 'EPI 1',
           cluster == 3 ~ 'DVM 1', 
           cluster == 4 ~ 'DVM 2', 
           cluster == 5 ~ 'DVM 3',
           cluster == 6 ~ 'RDVM'), 
           cluster = factor(cluster, levels = c('RDVM', 'EPI 1', 'EPI 2', 'DVM 1', 'DVM 2', 'DVM 3'), ordered = T))
require(viridis)
ggplot(data = cluster_series_test, aes(x = Local_Time, y = depth) ) +
  stat_bin_2d(aes(fill = after_stat(ndensity)), geom = "tile", bins = c(100, 100)) +
  scale_fill_gradientn(colors = viridis(100), limits = c(0.05,1), oob = scales::censor) + 
  labs(fill = "Relative Density (log)") +
  scale_y_reverse() +
  scale_x_continuous(breaks = c(21600, 64800), labels = c("06:00", "18:00")) +
  xlab("Local Time") +
  ylab("Depth (m)") +
  facet_grid(species.x~cluster) +
  theme_linedraw() +
  theme(legend.position = 'bottom')

## Spatial segregation of clusters ----
# add location to high-res-all
# create a stamp with coordinates for each tracking day
position.stmp <- combo_track %>% dplyr::select(latitude, longitude, kode)

# remove duplicates in days w. more than one location
position.stmp <- position.stmp[which(!duplicated(position.stmp$kode)),] 

# apply position data to high_res data frame
high_res_all <- left_join(high_res_all, position.stmp, by = c("kode"))

# load in bathymetry raster
r <- raster::raster('data/raw/global_bathy_0.01.nc')
## this is a global grid with pacific-centered coordinates (longitudes 0 to 360)
## raster::rotate converts from 0-360 longitudes to 180 longitudes (atlantic-centered)
r <- raster::rotate(r)
r <- raster::aggregate(r, fact = 5)
r_df <- raster::as.data.frame(r, xy = TRUE)
r_df <- r_df %>% filter(z < 0)

ggplot(data = r_df) +
  geom_sf(data = world, fill = 'black') +
  geom_contour(
    aes(x = x,
        y = y,
        z = z),
    color = "black",
    breaks = c(-1000)) +
  geom_point(data = high_res_all %>% 
               mutate(cluster = case_when(
                 cluster == 1 ~ 'EPI 2',
                 cluster == 2 ~ 'EPI 1',
                 cluster == 3 ~ 'DVM 1', 
                 cluster == 4 ~ 'DVM 2', 
                 cluster == 5 ~ 'DVM 3',
                 cluster == 6 ~ 'RDVM')) %>%
               mutate(
                 cluster = factor(cluster,
                                  levels = c(
                                    'RDVM',
                                    'EPI 1', 
                                    'EPI 2', 
                                    'DVM 1', 
                                    'DVM 2',
                                    'DVM 3'),
                                  ordered = T)),
             aes(x = longitude, y = latitude, color = cluster), size = 0.5) +
  scale_color_manual(values = c('goldenrod',
                                "yellow",
                                "yellow3",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  scale_x_continuous(breaks = c(-80, -70, -60, -50, -40)) +
  scale_y_continuous(breaks = c(10, 20, 30, 40)) +
  facet_wrap(~cluster, nrow = 2) +
  coord_sf(xlim = c(-80, -35), ylim = c(9, 45)) +
  labs(x = "Longitude",
       y = 'Latitude',
       color = 'Cluster') 

# Compare Shelf Behavior to __
shelf_series <- 
  combo_series %>%
  left_join(clust_stamp2) %>% 
  filter(
    kode %in% shelf$kode |
      cluster %in% c(1,3)) %>% 
  mutate(ogcluster = 0,
         cluster = case_when(
           is.na(cluster) ~ 'Shelf',
           cluster == 1 ~ 'EPI 2',
           cluster == 3 ~ 'EPI 1'),
         cluster = factor(cluster, levels = c('EPI 1', 'Shelf', 'EPI 2'), ordered = T))

require(viridis)
ggplot(data = shelf_series, aes(x = Local_Time, y = depth) ) +
  stat_bin_2d(aes(fill = after_stat(ndensity)), geom = "tile", bins = c(100, 100)) +
  scale_fill_gradientn(colors = viridis(100), limits = c(0,0.5), oob = scales::squish) + 
  labs(fill = "Relative Density (log)") +
  scale_y_reverse() +
  scale_x_continuous(breaks = c(21600, 64800), labels = c("06:00", "18:00")) +
  xlab("Local Time") +
  ylab("Depth (m)") +
  facet_grid(species~cluster) +
  theme_linedraw() +
  theme(legend.position = 'bottom')
