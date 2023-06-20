
# setup -------------------------------------------------------------------

library(tidyverse)
library(tags2etuff)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(ggspatial)
library(sf)
library(raster)
library(cmocean)
library(HMMoce)
library(scales)
devtools::load_all('analyzePSAT')

# create an object with continent boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")

# create an object with bin bounds
bins <- c(0, 10, 50, 100, 200, 300, 400, 500, 2000)

source('scripts/analysis/general_analysis.R')

combo_hdr <- read_csv('data/clean/combo_hdr.csv')


# read in track information, combined and individual
combo_track <- read_csv('data/clean/Tracks/combo_track.csv')
read_rds('data/clean/Tracks/all_sharks.rds') %>% 
  list2env(.GlobalEnv)

position.stmp <- combo_track %>% dplyr::select(latitude, longitude, kode)

# remove duplicates in days w. more than one location
position.stmp <- position.stmp[which(!duplicated(position.stmp$kode)),] 


# read in series information, combined and individual
combo_series <- read_csv('data/clean/Series/combo_series.csv')
all_series <- read_rds('data/clean/Series/all_sharks.rds')

# read in high resolution time-at-depth summaries
high_res <- read_csv('data/clean/high_resolution_summaries.csv')

# figure 1 ----------------------------------------------------------------
r <- raster::raster('data/raw/global_bathy_0.01.nc')
## this is a global grid with pacific-centered coordinates (longitudes 0 to 360)
## raster::rotate converts from 0-360 longitudes to 180 longitudes (atlantic-centered)
r <- raster::rotate(r)
r <- raster::aggregate(r, fact = 5)
r_df <- raster::as.data.frame(r, xy = TRUE)
r_df <- r_df %>% filter(z < 0)

start_pts <- mako %>% 
  # remove greg's tags
  filter(ptt != 78680 &
           ptt != 78682 &
           ptt != 78683) %>% 
  group_by(ptt) %>% 
  filter(DateTime == min(DateTime)) %>% 
  rbind(
    blue %>% 
      group_by(ptt) %>% 
      filter(DateTime == min(DateTime)))

end_pts <- 
  mako %>% 
  # remove greg's tags
  filter(ptt != 78680 &
           ptt != 78682 &
           ptt != 78683) %>%
  group_by(ptt) %>% 
  filter(DateTime == max(DateTime)) %>% 
  rbind(
    blue %>% 
      group_by(ptt) %>% 
      filter(DateTime == max(DateTime)))


ggplot(data = r_df) +
  geom_sf(data = world) + 
  geom_contour(
    aes(x = x,
        y = y,
        z = z),
    color = "black",
    breaks = c(-1000)) +
  scale_x_continuous(breaks = c(-80, -70, -60, -50, -40)) +
  scale_y_continuous(breaks = c(10, 20, 30, 40)) +
  coord_sf(xlim = c(-80, -35), ylim = c(9, 45)) +
  geom_path(data = (combo_track %>% 
              # remove greg's tags
              filter(ptt != 78680 &
                       ptt != 78682 &
                       ptt != 78683)),
            aes(longitude,
                latitude,
                color = species,
                group = ptt),
            lineend = 'butt') +
  # add tagging location as open green circle
  geom_point(data = start_pts, aes(longitude, latitude),
             color = 'green',
             shape = 16,
             size = 1.5) +
  
  # add pop-up locatoin as black X
  geom_point(data = end_pts, aes(longitude, latitude),
             fill = 'black',
             shape = 4,
             stroke = 0.75,
             size = 1.75) + 
  xlab("Longitude") +
  ylab("Latitude") +
  theme_linedraw() +
  theme(legend.position = 'none',
        element_blank())

## mesopelagic inset ----
# making map of % mesopelagic occupancy in a gridded raster
# See this youtube video for a helpful tutorial: https://www.youtube.com/watch?v=_bzqAuBnOag 

# set the extent of the region 
e <- extent(c(-80, -30, 5, 50))
# the reference system is:
crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
# now put them together in 50x50 raster, can back-calculate the area of each cell if necessary based on extent
r <- raster(nrow = 50, ncol = 50, ext = e, crs = crs)
p <- as(r@extent, 'SpatialPolygons')

# create a binary variable, in mesopelagic or not
combo_series$meso <- cut(combo_series$depth, breaks = c(0, 200, 2000), labels = c(1,2), include.lowest = TRUE, right = FALSE)

# count the observations in the mesopelagic for each tracking day
meso_freq <- 
  combo_series %>% 
  # filter out observations from over the shelf and with incomplete records
  # remove data from over the continental shelf
  filter(bathy <= -1000) %>% 
  group_by(ptt) %>% 
  count(Date) %>% # there are 1634 total tracking days with series data
  filter(n >= 288) %>% # there are 1249 tracking days with more than 50% of their series data
  right_join(combo_series, by = c("Date", "ptt")) %>% 
  filter(!is.na(n)) %>% 
  group_by(kode) %>% 
  count(meso, .drop = FALSE) %>% 
  mutate(tot = sum(n), perc = (n/tot) * 100) %>% 
  filter(meso == 2) %>% 
  dplyr::select("kode", "tot", "perc") %>% 
  rename(totalObservations = tot, percentMesopelagic = perc)

# Add % time in mesopelagic to track data for each shark
meso_tracks <- 
  combo_track %>% 
  left_join(meso_freq, by = c("kode"))

# separate by species
# an object for mako sharks
mako_tracks <- 
  meso_tracks %>% 
  filter(species == "I.oxyrinchus")
# and an object for blue sharks
blue_tracks <- 
  meso_tracks %>% 
  filter(species == "P.glauca")

# now convert each object to a grided raster with the same dimensions as r, where the value of each cell is the mean of percent time in the mesopelagic
makor <- raster::rasterize(mako_tracks[,c('longitude','latitude')], r, mako_tracks$percentMesopelagic, fun = mean)
bluer <- raster::rasterize(blue_tracks[,c('longitude','latitude')], r, blue_tracks$percentMesopelagic, fun = mean)

# convert these rasters back into data frames for plotting
makor_df <- 
  as.data.frame(makor, xy = TRUE) %>% 
  filter(!is.na(layer))

bluer_df <- 
  as.data.frame(bluer, xy = TRUE) %>% 
  filter(!is.na(layer))

# plot the data for mako sharks
ggplot(data = makor_df) +
  geom_tile(aes(x = x, y = y, fill = layer), color = "black") +
  geom_sf(data = world) +
  scale_fill_cmocean(name = "deep") +
  #scale_fill_cmocean(name = "deep", limits = c(0, 100), breaks = c(20, 40, 60, 80), labels = c(20, 40, 60, 80)) +
  coord_sf(xlim = c(-80, -35), ylim = c(10, 45)) +
  #guides(fill=guide_colorbar(title="Mesopelagic (%)", direction = "horizontal")) +
  theme_linedraw() +
  xlab("")+
  ylab("") +
  theme(legend.position = "none")
#theme(legend.position = c(0.80, 0.085), legend.box = "horizontal", legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid", color = "black"))

ggplot(data = bluer_df) +
  geom_tile(aes(x = x, y = y, fill = layer), color = "black") +
  geom_sf(data = world) +
  scale_fill_cmocean(name = "deep") +
  #scale_fill_cmocean(name = "deep", limits = c(0, 100), breaks = c(20, 40, 60, 80), labels = c(20, 40, 60, 80)) +
  coord_sf(xlim = c(-80, -35), ylim = c(10, 45)) +
  guides(fill=guide_colorbar(title="Mesopelagic (%)", direction = "horizontal")) +
  theme_linedraw() +
  xlab("")+
  ylab("") +
  theme(legend.position = "none") # +
  # theme(legend.position = c(0.50, 0.085),
  #       legend.box = "horizontal",
  #       legend.background = element_rect(fill = "white",
  #                                        size = 0.5,
  #                                        linetype = "solid",
  #                                        color = "black"))

# figure 2 ----------------------------------------------------------------
# blue shark series
blue_series <-
  combo_series %>% 
  filter(species == 'P.glauca')

# mako shark series
mako_series <-
  combo_series %>% 
  filter(species == 'I.oxyrinchus')

## blue shark -------------------------------------------------------------

btuff <- 
  read_csv('data/clean/Series/b133018_fullSeries.csv')
  
# series plot
main.b <-
  ggplot(data = btuff) +
  geom_ribbon(aes(x = DateTime_local,
                  y = bathy,
                  ymin = bathy,
                  ymax = 1000)) +
  geom_point(aes(x = DateTime_local,
                 y = depth,
                 color = temperature)) +
  scale_color_cmocean(name = "thermal", limits = c(5,30)) +
  geom_line(aes(x = DateTime_local,
                y = ild.5),
            color = "grey22",
            linewidth = 1.5) +
  geom_hline(yintercept = 200,
             color = "black",
             alpha = 0.4,
             linetype = 'dashed') +
  ylab("Depth (m)") +
  xlab("Month") +
  scale_y_reverse(limits = c(1000, 0), breaks = c(seq(0,1000,by = 200))) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))
#guides(color=guide_colorbar(title="Temp (ºC)", direction = "horizontal")) +
#theme(legend.position = c(0.18, 0.095), legend.box = "horizontal", legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid", color = "black"))

# save the plot
ggsave(plot = main.b,
       filename = 'products/figures/figure2/BlueVerticalHabitat.png',
       width = 300,
       height = 240,
       units = 'mm')

# Time-at-depth Histogram sub-plot 
# step 1: calculate the percentage of time in each depth bin for each day for each individual
blue.bin.depth.indiv <- 
  blue_series %>% 
  group_by(ptt, Date, dn) %>% 
  # place each depth measurement in its appropriate bin
  mutate(
    bin = cut(depth, breaks = bins, labels = c(seq(1:8)), include.lowest = TRUE)) %>% 
  # count the number of observations in each bin
  count(bin, .drop = FALSE) %>% 
  # find the total number of observations from each day from day and night
  mutate(tot = sum(n)) %>% 
  # filter out observations which are less than 50% complete
  filter(tot >= 144) %>% 
  # calculate the percent of time in each bin
  mutate(perc = (n / tot)*100) %>% 
  ungroup() %>%
  # summarize at the individual level
  group_by(ptt, dn, bin) %>%
  summarize(m = mean(perc))

# step 2: calculate the mean and standard error for each bin at the species level
blue.bin.depth <- 
  blue.bin.depth.indiv %>% 
  ungroup() %>%
  group_by(dn, bin) %>%
  summarize(perc = mean(m), se = sd(m)) %>% 
  mutate(bin = as.numeric(bin))

# step 3: graph it
sub.b <- ggplot() +
  geom_col(data = filter(blue.bin.depth, dn == "d"), aes(x = bin, y = -perc), width = 1, fill = NA, color = "black") +
  geom_errorbar(data = filter(blue.bin.depth, dn == "d"), aes(x = bin, ymin = (-perc), ymax = (-perc-se)), width = 0.5, color = "black", alpha = 0.6) +
  geom_col(data = filter(blue.bin.depth, dn == "n"), aes(x = bin, y = perc), width = 1, fill = "grey", color= "black") +
  geom_errorbar(data = filter(blue.bin.depth, dn == "n"), aes(x = bin, ymin = (perc), ymax = (perc+se)), width = 0.5, color = "black", alpha = 0.6) +
  scale_x_reverse(breaks = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5), labels = c("0", "10", "50", "100","200","300", "400", "500", "2000")) +
  scale_y_continuous(limits = c(-60,60), breaks = c(-60,-40,-20, 0,20,40,60), labels = c(60,40,20,0,20,40,60), position = "right") +
  xlab("Depth (m)") +
  ylab("Time at Depth (%)") +
  coord_flip() +
  theme_classic() + 
  theme(plot.margin = unit(rep(0, 4), "cm"),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))


#A viewport taking up a fraction of the plot area
ggsave(plot = sub.b,
       filename = 'products/figures/figure2/BlueAllTAD.png',
       width = 300,
       height = 240,
       units = 'mm')

## mako shark -------------------------------------------------------------
mtuff <- 
  read_csv('data/clean/Series/m163096_fullSeries.csv')

# series graph
breaks <- as.POSIXct(c("2017-11-01 00:00:00", "2017-12-01 00:00:00", "2018-01-01 00:00:00"))

main.m <-
  ggplot(data = mtuff) +
  geom_ribbon(aes(x = DateTime_local, y = bathy, ymin = bathy, ymax = 1000)) +
  geom_point(aes(x = DateTime_local, y = depth, color = temperature)) +
  scale_color_cmocean(name = "thermal", limits = c(5,30)) +
  geom_line(aes(x = DateTime_local, y = ild.5), color = "grey22", linewidth = 1.5) +
  geom_hline(yintercept = 200, color = "black", alpha = 0.4, linetype = 'dashed') +
  ylab("Depth (m)") +
  xlab("Month") +
  scale_y_reverse(limits = c(1000, 0), breaks = c(seq(0,1000,by=200))) +
  scale_x_continuous(breaks = breaks, labels = c("Nov", "Dec", "Jan")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))  # guides(color=guide_colorbar(title="Temp (ºC)", direction = "horizontal")) +
  # theme(legend.position = c(0.18, 0.095), legend.box = "horizontal", legend.background = element_rect(fill = "white", linewidth = 0.5, linetype = "solid", color = "black"))

# save the plot
ggsave(plot = main.m,
       filename = 'products/figures/figure2/MakoVerticalHabitat.png',
       width = 300,
       height = 240,
       units = 'mm')

# Time-at-depth Histogram sub-plot (2)
# step 1: calculate the percentage of time in each depth bin for each day for each individual
mako.bin.depth <- mako_series %>% group_by(ptt, Date, dn)
mako.bin.depth$bin <- cut(mako.bin.depth$depth, breaks = bins, labels = c(seq(1:8)), include.lowest = TRUE)

mako.bin.depth <- mako.bin.depth %>% count(bin, .drop = FALSE)
mako.bin.depth <- mako.bin.depth %>% mutate(tot = sum(n))
mako.bin.depth <- mako.bin.depth %>% filter(tot >= 144) # 144 because half of 288 (half of the whole day)
mako.bin.depth <- mako.bin.depth %>% mutate(perc = (n / tot)*100)

mako.bin.depth <- ungroup(mako.bin.depth) %>% group_by(ptt, dn, bin) %>% summarize(m = mean(perc))

# step 2: calculate the mean and standard error for each bin
mako.bin.depth <- ungroup(mako.bin.depth) %>% group_by(dn, bin) %>% summarize(perc = mean(m), se = sd(m))

# step 3: graph it
mako.bin.depth$bin <- as.numeric(mako.bin.depth$bin)
sub.m <- ggplot() +
  geom_col(data = filter(mako.bin.depth, dn == "d"), aes(x = bin, y = -perc), width = 1, fill = NA, color = "black") +
  geom_errorbar(data = filter(mako.bin.depth, dn == "d"), aes(x = bin, ymin = (-perc), ymax = (-perc-se)), width = 0.5, color = "black", alpha = 0.6) +
  geom_col(data = filter(mako.bin.depth, dn == "n"), aes(x = bin, y = perc), width = 1, fill = "grey", color= "black") +
  geom_errorbar(data = filter(mako.bin.depth, dn == "n"), aes(x = bin, ymin = (perc), ymax = (perc+se)), width = 0.5, color = "black", alpha = 0.6) +
  scale_x_reverse(breaks = c(0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5), labels = c("0", "10", "50", "100","200","300", "400", "500", "2000")) +
  scale_y_continuous(limits = c(-60,60), breaks = c(-60,-40,-20, 0,20,40,60), labels = c(60,40,20,0,20,40,60), position = "right") +
  xlab("Depth (m)") +
  ylab("Time at Depth (%)") +
  coord_flip() +
  theme_classic() +
  theme(plot.margin = unit(rep(0,4), 'cm'),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))


#Just save the plot
ggsave(plot = sub.m,
       filename = 'products/figures/figure2/MakoAllTAD.png',
       width = 300,
       height = 240,
       units = 'mm')
# figure 3 ----------------------------------------------------------------

source('scripts/analysis/heirarchical_clustering.R')

# create the dendrogram object
cluster_dendro <- 
  as.dendrogram(Clust2) 

# make a key for navigating the dendrogram

# create a list of species ordered by leaves on the dendrogram
leaf.species <- 
  high_res$species 
leaf.species <- 
  leaf.species[order.dendrogram(cluster_dendro)] 

# vector of leaf orders (index) relative to original input (high_res) 
leaf.index <- 
  order.dendrogram(cluster_dendro)

# create a list of clusters, ordered by leaves on the dendrogram
leaf.clusters <- 
  high_res$cluster
leaf.clusters <- 
  leaf.clusters[order.dendrogram(cluster_dendro)]

leaf.key <- data.frame(
  species = leaf.species,
  index = leaf.index,
  cluster = leaf.clusters) %>% 
  mutate(clusterName = 
           case_when(
             cluster == 1 ~ 'DVM 1',
             cluster == 2 ~ 'Epipelagic',
             cluster == 3 ~ 'DVM 2',
             cluster == 4 ~ 'DVM 3',
             cluster == 5 ~ 'DVM 5',
             cluster >= 6 ~ as.character(cluster)),
         clusterName = factor(clusterName, 
                          levels = c('Epipelagic',
                                     'DVM 1',
                                     'DVM 2',
                                     'DVM 3',
                                     'DVM 4',
                                     'DVM 5',
                                     '6',
                                     '7',
                                     '8',
                                     '9',
                                     '10',
                                     '11',
                                     '12'),
                          ordered = T)) 

rm(leaf.species,
   leaf.index,
   leaf.clusters)

# prune leaves from rare clusters

# list the leaves which belong to rare clusters
rare <- list(
  six = c(113,696,737),
  seven = c(117, 208, 118, 121, 122, 123, 116, 119),
  eight = c(124),
  nine = c(786, 557, 564, 184, 725, 615, 785, 236, 309, 442, 619),
  ten = c(336),
  eleven = c(413),
  twelve = c(642))

# create a key for all clusters
color.leafs.all <- 
  leaf.key %>%
  # order of the leaves on the pruned tree
  mutate(order = seq(1:786)) %>% 
  arrange(clusterName)

# create a key for common clusters
color.leafs <- 
  leaf.key %>%
  filter(cluster <= 5) %>% 
  # order of the leaves on the pruned tree
  mutate(order = seq(1:760)) %>% 
  arrange(clusterName)

# prune the leaves of rare clusters
pruned_dendro <- 
  cluster_dendro %>% 
  prune(
    (rare %>% unlist)) 
  
# reorder the leaves so clusters progress from shallowest to deepest
ordered_dendro <-
  pruned_dendro %>% 
  rotate(color.leafs$order)

# for ALL clusters
# ordered_dendro <-
#   cluster_dendro %>% 
#   rotate(color.leafs.all$order)

# assign cluster labels to pruned dendrogram
labels_colors(ordered_dendro) <- 
  color.leafs$cluster

# # for ALL clusters
# labels_colors(ordered_dendro) <- 
#   color.leafs.all$cluster

# color the labels by cluster 
labels_colors(ordered_dendro)
plot(ordered_dendro)

# make a list of ggplot colors
cmocean("deep")(5)
show_col(cmocean("deep")(5))

# color the branches of the dendrogram by cluster: 
color_dendro <- 
  color_branches(
    ordered_dendro,
    clusters = color.leafs$cluster,
    col = c("#FFFF5CFF", "#78CEA3FF", "#488E9EFF", "#404C8BFF", "#281A2CFF")) 

# for ALL cclusters
color_dendro <- 
  color_branches(
    ordered_dendro,
    clusters = color.leafs.all$cluster,
    col = c("#FFFF5CFF", "#78CEA3FF", "#488E9EFF", "#404C8BFF", "#281A2CFF",
            'red', 'red', 'red', 'red', 'red', 'red', 'red' )) 

# create a pdf object for the dendrogram plot
pdf(file = 'products/figures/figure3/color_dendrogram.pdf', width = 10, height = 7)

# save the plot to the above pdf object
plot(color_dendro)

# unlink pdf object
dev.off()


## fig. 3 subplots -------------------------------------------------------------

# Step 1: recreate daily TAD summaries from series data
## can't use those transmitted by tags due to lower temporal resolution
blue.bin.depth <- 
  combo_series %>% 
  filter(species == 'P.glauca') %>% 
  mutate(
    bin = 
      cut(depth,
          breaks = bins,
          labels = c(seq(1:8)),
          right = FALSE,
          include.lowest = TRUE)) %>% 
  group_by(ptt, Date) %>%
  count(bin, .drop = FALSE) %>% 
  mutate(tot = sum(n)) %>% 
  filter(tot >= 280) %>% 
  mutate(perc = (n / tot)*100) %>% 
  ungroup()

mako.bin.depth <- 
  combo_series %>% 
  filter(species == 'I.oxyrinchus') %>% 
  mutate(
    bin = 
      cut(depth,
          breaks = bins,
          labels = c(seq(1:8)),
          right = FALSE,
          include.lowest = TRUE)) %>% 
  group_by(ptt, Date) %>%
  count(bin, .drop = FALSE) %>% 
  mutate(tot = sum(n)) %>% 
  filter(tot >= 280) %>% 
  mutate(perc = (n / tot)*100) %>% 
  ungroup()

tad_combo <- 
  rbind(mako.bin.depth, blue.bin.depth) %>% 
  mutate(kode = paste(Date, ptt, sep = "_")) %>% 
  inner_join(clust_stamp2, by = c("kode")) %>% 
  mutate(bin = as.numeric(bin))

# Step 2: use tad_combo to make heatmaps for each cluster
library(scales) 
library(viridis)

add_rows <- function(x) {
  c1 <- filter(tad_combo, cluster == x)
  
  date_stmp <- c1 %>% dplyr::select(c(Date, kode))
  date_stmp <- date_stmp[which(!duplicated(date_stmp$kode)),] %>%
    mutate(Row = c(seq(1:nrow(.))))
  c1 <- inner_join(c1, date_stmp, by = c("Date", "kode"))
  return(c1)
}
  
heatmap_tad <- 
  c(1:12) %>% 
  map(
    ~ add_rows(.)
  ) %>% 
  set_names(paste('c', c(1:12), sep = ''))

rm(tad_combo, mako.bin.depth, blue.bin.depth)

# plot and save the heatmaps
for(x in names(heatmap_tad)) {
  heatmap <- 
    ggplot(data = heatmap_tad[[x]]) +
    geom_tile(aes(x = Row, y = bin, fill = perc), linejoin = "round") +
    scale_fill_gradientn(colors = viridis(100), limits = c(0,50), oob = scales::squish) +
    labs(x = 'Day', y = 'Depth (m)', title = x) +
    theme_minimal() +
    theme(legend.position = 'top')+
    scale_y_reverse(breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5), labels = c("0", "10", "50", "100", "200", "300", "400", "500", "2000")) +
    theme(aspect.ratio = 20/10)
  
  ggsave(plot = heatmap,
         filename = paste('products/figures/figure3/', x, '.png', sep = ''),
         width = 300,
         height = 240,
         units = 'mm')
}




# figure 4 ----------------------------------------------------------------

cluster_series <- 
  combo_series %>% 
  inner_join(clust_stamp2,
             by = 'kode') %>% 
  filter(cluster <= 5) %>% 
  mutate(ogCluster = cluster,
         cluster = case_when(
           cluster == 1 ~ 'DVM 1',
           cluster == 2 ~ 'Epipelagic',
           cluster == 3 ~ 'DVM 2',
           cluster == 4 ~ 'DVM 3',
           cluster == 5 ~ 'DVM 4'),
         cluster = factor(cluster, levels = c('Epipelagic', 'DVM 1', 'DVM 2', 'DVM 3', 'DVM 4'), ordered = T))
  
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

# figure 5 ----------------------------------------------------------------

library(mclogit)

# get cluster information
source('scripts/analysis/heirarchical_clustering.R')


## formatting --------------------------------------------------------------

combo_mod <- 
  high_res %>% 
  # remove rare clusters
  filter(cluster <= 5 ) %>% 
  # reset cluster names
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
      'DVM 4')
    )) %>%
  # set the reference level to the most common cluster
  mutate(
    clus2 = relevel(
      as.factor(cluster),
      ref = 'DVM 1'),
    # reference as blue sharks because we have more tags from them
    species = relevel(
      as.factor(species),
      ref = 'P.glauca'),
    ptt = factor(ptt),
    depth = bathy,
    l2 = lunar) %>% 
  # select response and predictors
  dplyr::select(c(Date,
                  kode,
                  species,
                  ptt,
                  cluster,
                  clus2,
                  latitude,
                  longitude,
                  depth,
                  sst,
                  sst_sd,
                  ssh,
                  ssh_sd,
                  eke,
                  ild.5,
                  n2,
                  chloro,
                  l2,
                  lunar)) %>% 
  # convert lunar to factor variable based on phase
  mutate(
    lunar = case_when(
      lunar >= 0.75 ~ 2,
      lunar < 0.75 ~ 1),
    lunar = relevel(
      # setting lunar to a factor with reference 1
      as.factor(lunar), ref = 1))

# train the best fitting model
m.mod <- 
  mblogit(
    formula = clus2 ~ ssh +
      ssh_sd +
      lunar +
      n2 +
      ssh:n2 +
      ssh:lunar +
      ssh_sd:lunar +
      ssh_sd:species + 
      species:n2,
    random = ~1|ptt,
    data = combo_mod)

## Effect of SSH -----------------------------------------------------------

# vary ssh while holding other predictors constant
p_ssh <- 
  tibble(
    ssh = rep(seq(-1,1, 0.025),times = 34), 
    lunar = as.factor(rep(rep(c(0,1),each = 81), times = 17)),
    # hold ssh_sd at the mean value
    ssh_sd = rep(c(0.02605838), times = 2754),
    n2 = rep(c(1.472e-05), times = 2754),
    ptt = as.factor(rep(c("106754",
                          "133016",
                          "133017",
                          "133018",
                          "133021",
                          "141247",
                          "141254",
                          "141255",
                          "141256",
                          "141257",
                          "141258",
                          "141259",
                          "154096",
                          "163096",
                          "163097",
                          "163098",
                          "206771"),
                        each = 162))) %>% 
  # bind to another tibble with high values of n2
  rbind(
    tibble(
      ssh = rep(seq(-1,1, 0.025),times = 34), 
      lunar = as.factor(rep(rep(c(0,1),each = 81), times = 17)),
      # hold ssh_sd at the mean value
      ssh_sd = rep(c(0.02605838), times = 2754),
      n2 = rep(c(7.457e-05), times = 2754),
      ptt = as.factor(rep(c("106754",
                            "133016",
                            "133017",
                            "133018",
                            "133021",
                            "141247",
                            "141254",
                            "141255",
                            "141256",
                            "141257",
                            "141258",
                            "141259",
                            "154096",
                            "163096",
                            "163097",
                            "163098",
                            "206771"),
                          each = 162)))) %>% 
  left_join(high_res %>% 
              dplyr::select(ptt, species) %>%
              distinct() %>% 
              mutate(ptt = factor(ptt)),
            by = 'ptt')
  
# add predictions based on the best performing model
p_ssh <- 
  cbind(
    p_ssh,
    predict(m.mod,
            newdata = p_ssh,
            type = "response",
            se.fit = FALSE))
# condense predictions to a single column (value) dependent on cluster (variable)
 lssh <- 
  reshape2::melt(p_ssh,
       id.vars = c("ssh",
                   "lunar",
                   "ssh_sd",
                   'n2',
                   "ptt",
                   "species"),
       measure.vars = c("DVM 1",
                        "Epipelagic",
                        "DVM 2",
                        "DVM 3",
                        "DVM 4")) %>% 
  mutate(
    lunar = factor(lunar))

rm(p_ssh)

# create a mean line by averaging all individuals 
opall <- 
  lssh %>% 
  group_by(species, 
           variable,
           ssh,
           lunar,
           n2) %>%
  summarize(
    # lower limit of predictions
    LL = min(value),
    # upper limit of predictions
    UL = max(value),
    # average predictions 
    value = mean(value)) %>% 
  ungroup()

# plot relationships for blues
opall %>% 
  filter(species == 'I.oxyrinchus') %>%
  ggplot(aes(x = ssh,
             y = value,
             color = variable)) + 
  geom_ribbon(aes(x = ssh,
                  ymin = LL,
                  ymax = UL,
                  fill = variable),
              color = 'white',
              linewidth = 0.25,
              alpha = 0.25) +
  geom_line(linewidth = 1.5, alpha = 1.2) +
  scale_x_continuous(expand = c(0,0.01)) +
  scale_y_continuous(expand = c(0,0.01)) +
  # use this code to generate the below colors (with some edits to yellow): show_col(cmocean(name = 'deep')(5))
  scale_color_manual(values = c("#78CEA3FF",
                                "#FFFF5CFF",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  scale_fill_manual(values = c("#78CEA3FF",
                                "#FFFF5CFF",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  facet_wrap(n2~lunar) + 
  labs(fill = "Cluster") +
  guides(color = 'none') +
  ylab("Probability") +
  theme_linedraw()

## Effect of SSH_SD -----------------------------------------------------------

# vary ssh_sd while holding other predictors constant
p_ssh_sd <- 
  data.frame(
    ssh_sd = rep(seq(0,0.13, 0.001625),times = 34), 
    lunar = as.factor(rep(rep(c(0,1),each = 81),times = 17)),
    # hold n2 at median value
    n2 = rep(c(3.096e-05), times = 2754),
    # hold ssh at the median value
    ssh = rep(c(-0.0500), times = 2754),
    ptt = as.factor(rep(c("106754",
                          "133016",
                          "133017",
                          "133018",
                          "133021",
                          "141247",
                          "141254",
                          "141255",
                          "141256",
                          "141257",
                          "141258",
                          "141259",
                          "154096",
                          "163096",
                          "163097",
                          "163098",
                          "206771"),
                        each = 162))) %>% 
  left_join(
    (high_res %>% 
       select(ptt, species) %>% 
       distinct() %>% 
       mutate(ptt = factor(ptt))),
    by = 'ptt') 

# add predictions based on the best performing model
p_ssh_sd <- 
  cbind(
    p_ssh_sd,
    predict(m.mod,
            newdata = p_ssh_sd,
            type = "response",
            se.fit = FALSE))
# condense predictions to a single column (value) dependent on cluster (variable)
lssh_sd <- 
  reshape2::melt(p_ssh_sd,
       id.vars = c("ssh",
                   "lunar",
                   'n2',
                   "ssh_sd",
                   "ptt",
                   "species"),
       measure.vars = c("DVM 1",
                        "Epipelagic",
                        "DVM 2",
                        "DVM 3",
                        "DVM 4"))

rm(p_ssh_sd)

# create a mean line by averaging all individuals 
op <- 
  lssh_sd %>% 
  group_by(ssh_sd,
           lunar,
           variable,
           species) %>%
  summarize(
    # lower limit of predictions
    LL = min(value),
    # upper limit of predictions
    UL = max(value),
    # average predictions 
    value = mean(value)) %>% 
  ungroup()

# Set the depth order of the clusters
op$variable <- 
  factor(op$variable,
         levels = c(2, 1, 3, 4, 5),
         ordered = TRUE) 

# plot relationships for makos
ggplot(data = lssh_sd) + 
  geom_path(aes(x = ssh_sd, y = value, color = (factor(variable,
                                                    levels = c(2, 1, 3, 4, 5),
                                                    ordered = TRUE) )), alpha = 0.25) +
  geom_line(data = op,aes(x = ssh_sd, y = value, color = variable), size = 1.2, alpha = 1.2) +
  scale_color_manual(values = c("#FFFF5CFF",
                                "#78CEA3FF",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  facet_grid(species~.) + 
  labs(color = "Cluster") +
  ylab("Probability") +
  theme_minimal()

op %>% 
  ggplot(aes(x = ssh_sd,
             y = value,
             color = variable)) + 
  geom_ribbon(aes(x = ssh_sd,
                  ymin = LL,
                  ymax = UL,
                  fill = variable),
              color = 'white',
              linewidth = 0.25,
              alpha = 0.25) +
  geom_line(linewidth = 1.5, alpha = 1.2) +
  scale_x_continuous(expand = c(0,0.01)) +
  scale_y_continuous(expand = c(0,0.01)) +
  # use this code to generate the below colors (with some edits to yellow): show_col(cmocean(name = 'deep')(5))
  scale_color_manual(values = c("#78CEA3FF",
                                "#FFFF5CFF",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  scale_fill_manual(values = c("#78CEA3FF",
                               "#FFFF5CFF",
                               "#488E9EFF",
                               "#404C8BFF",
                               "#281A2CFF")) +
  facet_wrap(species~lunar) + 
  labs(x = 'Sea-surface Height Standard Deviation (m)', 
       y = 'Probability',
       fill = "Cluster") +
  guides(color = 'none') +
  theme_linedraw()

# supplementary -----------------------------------------------------------


## industrial series plots -------------------------------------------------

# create an object with series data for the cluster of interest
c1 <- 
  cluster_series %>%
  # update with the cluster of interest
  filter(cluster == 'DVM 4')

# get unique dates 
series_dates <- unique(c1$Date) 

for(i in 1:length(series_dates)) {
  # find all observations from on each unique day
  day <- c1 %>%
    filter(Date == series_dates[i])
  
  # filter out days w. <60% coverage
  sharks <- 
    day %>% 
    count(ptt) %>%
    filter(n >= 260)
  
  # create a list of tag ids
  ptt_n <- c(sharks$ptt)
  full_series <- 
    day %>% 
    filter(ptt %in% ptt_n)
    
  if (nrow(full_series) != 0) {
    ggplot() + 
      geom_line(data = filter(full_series, species.x == "P.glauca"), aes(x = DateTime, y = depth), color = "#00BFC4") +
      geom_line(data = filter(full_series, species.x == "I.oxyrinchus"), aes(x = DateTime, y = depth), color = "#F8766D") +
      scale_y_reverse() +
      facet_grid(ptt~.) +
      theme_classic()
    
    cluster <- 
      c1 %>% 
      pull(cluster) %>% 
      unique() %>% 
      as.character()
    
    ggsave(paste("products/series/", cluster, '/', full_series[1,1], ".png", sep = ""),
           width = 200,
           height = 150,
           units = 'mm')
    
    rm(day, sharks, ptt_n, full_series)
  }
  else {
    rm(day, sharks, ptt_n, full_series)
  }
}


## all timeseries ---------------------------------------------------------

cluster_series %>%
  # update with the cluster of interest
  filter(cluster == 'DVM 3') %>% 
  ggplot() +
  geom_path(aes(x = Local_Time, y = depth, color = cluster), alpha = 0.5) + 
  scale_y_reverse() +
  facet_grid(cluster ~ .)

## when do clusters occur ----
library(ggridges)

high_res %>%
  filter(cluster <= 5) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'DVM 1',
    cluster == 2 ~ 'Epipelagic',
    cluster == 3 ~ 'DVM 2',
    cluster == 4 ~ 'DVM 3',
    cluster == 5 ~ 'DVM 4'),
    cluster = factor(cluster, levels = c(
      'DVM 4',
      'DVM 3',
      'DVM 2',
      'DVM 1',
      'Epipelagic'))) %>% 
  mutate(yday = lubridate::yday(Date)) %>%  # pull(yday) %>% summary()
  ggplot(aes(x = yday, y = cluster, fill = factor(cluster))) +
  geom_density_ridges(scale = 2.5, rel_min_height = 0.01, alpha = 0.65) +
  scale_x_continuous(limits = c(0, 366), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c("#281A2CFF",
                               "#404C8BFF",
                               "#488E9EFF",
                               "#78CEA3FF",
                               "#FFFF5CFF"),
                    name = 'Cluster') +
  labs(x = 'Day of the Year', y = 'Cluster') +
  theme_classic()

## corrected yday boxplot

high_res %>%
  filter(cluster <= 5) %>% 
  rbind((high_res %>% 
           filter(cluster <= 5) %>% 
           mutate(cluster = 0))) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'DVM 1',
    cluster == 2 ~ 'Epipelagic',
    cluster == 3 ~ 'DVM 2',
    cluster == 4 ~ 'DVM 3',
    cluster == 5 ~ 'DVM 4',
    cluster == 0 ~ 'All Data'),
    cluster = factor(cluster, levels = c(
      'DVM 4',
      'DVM 3',
      'DVM 2',
      'DVM 1',
      'Epipelagic',
      'All Data'))) %>% 
  mutate(yday = lubridate::yday(Date)) %>% 
  # add a correction to center yday on tag date
  mutate(yday = if_else(
           yday >=238,
           (yday - 238),
           (yday + 127))) %>%  # pull(yday) %>% summary()
  ggplot(aes(x = cluster, y = yday, fill = cluster)) +
  geom_jitter(aes(color = cluster), alpha = 0.4) +
  # add extent of the seasons
  geom_rect(aes(ymin = 0,
                ymax = 27,
                xmin = 0,
                xmax = 7),
            fill = 'grey80',
            color = 'grey75',
            alpha = 0.01) +
  geom_rect(aes(ymin = 27,
                ymax = 117,
                xmin = 0,
                xmax = 7),
            fill = NA,
            color = 'grey75',
            alpha = 0.02) +
  geom_rect(aes(ymin = 117,
                ymax = 205,
                xmin = 0,
                xmax = 7),
            fill = 'grey80',
            color = 'grey75',
            alpha = 0.01) +
  geom_boxplot(alpha = 0.75) + 
  scale_fill_manual(values = c("#281A2CFF",
                               "#404C8BFF",
                               "#488E9EFF",
                               "#78CEA3FF",
                               "#FFFF5CFF",
                               '#AAAAAAAA'),
                    name = 'Cluster') +
  scale_color_manual(values = c("#281A2CFF",
                               "#404C8BFF",
                               "#488E9EFF",
                               "#78CEA3FF",
                               "#FFFF5CFF",
                               '#AAAAAAAA')) +
  scale_y_continuous(limits = c(0, 275)) +
  annotate("text", x = 0.25, y = 65, label = "Fall", size = 3) +
  annotate("text", x = 0.25, y = 160, label = "Winter", size = 3) +
  annotate("text", x = 0.25, y = 255, label = "Spring", size = 3) +
  labs(x = 'Cluster', y = 'Days since Tagging') +
  guides(color = 'none') +
  coord_flip() +
  facet_wrap(~species) +
  theme_classic()


## Cluster Location Map ----------------------------------------------------

high_res %>% 
  filter(cluster <= 5) %>% 
  group_by(cluster) %>% 
  summarize(lon.min = min(x),
            lon.max = max(x),
            lat.min = min(y),
            lat.max = max(y),
            med.lon = median(x),
            med.lat = median(y)) %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'DVM 1',
    cluster == 2 ~ 'Epipelagic',
    cluster == 3 ~ 'DVM 2',
    cluster == 4 ~ 'DVM 3',
    cluster == 5 ~ 'DVM 4'),
    cluster = factor(cluster, levels = c(
      'DVM 1',
      'DVM 2',
      'DVM 3',
      'DVM 4',
      'Epipelagic'))) %>% 
  ggplot() +
  geom_sf(data = world) +
  geom_rect(aes(xmin = lon.min,
                xmax = lon.max,
                ymin = lat.min,
                ymax = lat.max,
                fill = factor(cluster)),
            alpha = 0.5) +
  geom_point(data = (high_res %>% 
               filter(cluster <= 5) %>% 
               mutate(cluster = case_when(
                 cluster == 1 ~ 'DVM 1',
                 cluster == 2 ~ 'Epipelagic',
                 cluster == 3 ~ 'DVM 2',
                 cluster == 4 ~ 'DVM 3',
                 cluster == 5 ~ 'DVM 4'),
                 cluster = factor(cluster, levels = c(
                   'DVM 1',
                   'DVM 2',
                   'DVM 3',
                   'DVM 4',
                   'Epipelagic')))),
             aes(x = x, y = y),
             color = 'grey22',
             alpha = 0.6,
             shape = 20) +
  geom_hline(aes(yintercept = med.lat),
             linetype = 'dashed',
             alpha = 0.8) +
  geom_vline(aes(xintercept = med.lon),
             linetype = 'dashed',
             alpha = 0.8) +
  geom_point(aes(x = med.lon,
                 y = med.lat),
             color = 'red',
             size = 3) +
  scale_fill_manual(values = c("#78CEA3FF",
                               "#488E9EFF",
                               "#404C8BFF",
                               "#281A2CFF",
                               "#FFFF5CFF"),
                    name = 'Cluster') +
  scale_color_manual(values = c("#78CEA3FF",
                               "#488E9EFF",
                               "#404C8BFF",
                               "#281A2CFF",
                               "#FFFF5CFF"),
                    name = 'Cluster') +
  coord_sf(xlim = c(-80, -35), ylim = c(9, 45)) +
  guides(color = 'none') +
  labs(x = 'Longitude', y = 'Latitude') +
  facet_wrap(~cluster, nrow = 3)
