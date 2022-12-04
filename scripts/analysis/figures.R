
# setup -------------------------------------------------------------------

library(tidyverse)
library(tags2etuff)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(ggspatial)
library(sf)

world <- ne_countries(scale = "medium", returnclass = "sf")

source('scripts/analysis/general_analysis.R')

combo_hdr <- read_csv('data/clean/combo_hdr.csv')


# read in track information, combined and individual
combo_track <- read_csv('data/clean/Tracks/combo_track.csv')
read_rds('data/clean/Tracks/all_sharks.rds') %>% 
  list2env(.GlobalEnv)

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
  scale_x_continuous(minor_breaks = c(-70, -50, -30)) +
  scale_y_continuous(minor_breaks = c(10, 20, 30, 40)) +
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
             shape = 1,
             size = 3) +
  
  # add pop-up locatoin as black X
  geom_point(data = end_pts, aes(longitude, latitude),
             color = 'black',
             shape = 4,
             size = 3) + 
  xlab("Longitude") +
  ylab("Latitude") +
  theme_linedraw() +
  theme(element_blank())

# figure 2 ----------------------------------------------------------------

## blue shark -------------------------------------------------------------

# select the meta data for an individual with good records
b_hdr <- 
  combo_hdr %>% 
  filter(ptt == 133018)

btuff <-
  combo_series %>% 
  filter(ptt == 133018) %>% 
  filter(!is.na(temperature)) %>% 
  mutate(bathy = bathy * -1) %>% 
  left_join((high_res %>% 
              select(kode, ild.5)),
            by = 'kode') %>% 
  left_join((omega_combo %>% 
               select(kode, ild.5)))

# series plot
#main.b <- 
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
            alpha = 0.4,
            size = 1.5) +
  geom_ribbon(aes(x = DateTime_local,
                  y = bathy,
                  ymin = 200,
                  ymax = 1000),
              color = "black",
              fill = "NA",
              alpha = 0.8,
              linetype = 'dashed') +
  ylab("Depth (m)") +
  xlab("Month") +
  scale_y_reverse(limits = c(1000, 0)) +
  theme_classic() +
  theme(legend.position = "none")
#guides(color=guide_colorbar(title="Temp (ºC)", direction = "horizontal")) +
#theme(legend.position = c(0.18, 0.095), legend.box = "horizontal", legend.background = element_rect(fill = "white", size = 0.5, linetype = "solid", color = "black"))

# Time-at-depth Histogram sub-plot 
# step 1: calculate the percentage of time in each depth bin for each day for each individual
blue.bin.depth <- blue_series %>% group_by(ptt, Date, dn)
blue.bin.depth$bin <- cut(blue.bin.depth$depth, breaks = bins, labels = c(seq(1:8)), include.lowest = TRUE)

blue.bin.depth <- blue.bin.depth %>% count(bin, .drop = FALSE)
blue.bin.depth <- blue.bin.depth %>% mutate(tot = sum(n))
blue.bin.depth <- blue.bin.depth %>% filter(tot >= 144)
blue.bin.depth <- blue.bin.depth %>% mutate(perc = (n / tot)*100)

blue.bin.depth <- ungroup(blue.bin.depth) %>% group_by(ptt, dn, bin) %>% summarize(m = mean(perc))

# step 2: calculate the mean and standard error for each bin
blue.bin.depth <- ungroup(blue.bin.depth) %>% group_by(dn, bin) %>% summarize(perc = mean(m), , se = sd(m))

# step 3: graph it
blue.bin.depth$bin <- as.numeric(blue.bin.depth$bin)
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
  theme(plot.margin = unit(rep(0, 4), "cm"))


#A viewport taking up a fraction of the plot area
library(grid)
vp <- viewport(width = 0.4, height = 0.4, x = 0.8, y = 0.3)
#Just draw the plot twice
png("BlueVertHabitat.png")
main.b
print(sub.b, vp = vp)
dev.off()

## mako shark -------------------------------------------------------------



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
  cluster = leaf.clusters
)

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

# prune the leaves of rare clusters
pruned_dendro <- 
  cluster_dendro %>% 
  prune(
    (rare %>% unlist))

# create a key for common clusters
color.leafs <- 
  leaf.key %>%
  filter(cluster <= 5)

# assign cluster labels to pruned dendrogram
labels_colors(pruned_dendro) <- 
  color.leafs$cluster

# color the labels by cluster 
labels_colors(pruned_dendro)
plot(pruned_dendro)

# make a list of ggplot colors
cmocean("deep")(5)
show_col(cmocean("deep")(5))

# color the branches of the dendrogram by cluster: 
color_dendro <- 
  color_branches(
    pruned_dendro,
    clusters = color.leafs$cluster,
    col = c("#FFFF5CFF", "#78CEA3FF", "#488E9EFF", "#404C8BFF", "#281A2CFF")) 
plot(color_dendro)

## fig. 3 subplots -------------------------------------------------------------

bins <- c(0, 10, 50, 100, 200, 300, 400, 500, 2000)

# Step 1: recreate daily TAD summaries from series data
## can't use those transmitted by tags due to lower coverage
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
  c(1:8) %>% 
  map(
    ~ add_rows(.)
  ) %>% 
  set_names(paste('c', c(1:8), sep = ''))

rm(tad_combo, mako.bin.depth, blue.bin.depth)

ggplot(data = heatmap_tad$c5) +
  geom_tile(aes(x = Row, y = bin, fill = perc), linejoin = "round") +
  scale_fill_viridis(limits = c(0,75), oob = scales::squish) +
  theme_minimal() +
  scale_y_reverse(breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5), labels = c("0", "10", "50", "100", "200", "300", "400", "500", "2000")) +
  theme(aspect.ratio = 20/10)


# figure 4 ----------------------------------------------------------------

cluster_series <- 
  combo_series %>% 
  inner_join(clust_stamp2,
             by = 'kode') %>% 
  filter(cluster <= 5) %>% 
  mutate(ogCluster = cluster,
         cluster = factor(cluster, levels = c(2,1,3,4,5)))

ggplot(data = cluster_series, aes(x = Local_Time, y = depth) ) +
  stat_bin_2d(aes(fill = (log10((..ndensity..)+0.01))), geom = "tile", bins = c(100, 100)) +
  scale_fill_cmocean(name = "dense") + 
  labs(fill = "Relative Density (log)") +
  scale_y_reverse(limits = c(1000, 0), breaks = c(0, 500, 1000), labels = c("0", "500", "1000")) +
  scale_x_continuous(breaks = c(21600, 64800), labels = c("06:00", "18:00")) +
  xlab("Local Time") +
  ylab("Depth (m)") +
  facet_grid(species.x~cluster) +
  theme_linedraw()


# figure 5 ----------------------------------------------------------------

library(mclogit)

# get cluster information
source('scripts/analysis/heirarchical_clustering.R')


## formatting --------------------------------------------------------------

combo_mod <- 
  high_res %>% 
  # remove rare clusters
  filter(cluster <= 5 ) %>% 
  # set the reference level to the most common cluster
  mutate(
    clus2 = relevel(
      as.factor(cluster),
      ref = 1),
    # reference as blue sharks because we have more tags from them
    species = relevel(
      as.factor(cluster),
      ref = 1),
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

m.mod <- 
  mblogit(
    formula = clus2 ~ ssh + ssh_sd + lunar + species + n2 + ssh:n2,
    random = ~1|ptt,
    data = combo_mod)

## Effect of SSH -----------------------------------------------------------

# vary ssh while holding other predictors constant
p_ssh <- 
  data.frame(
    ssh = rep(seq(-1,1, 0.025),times = 34), 
    lunar = as.factor(rep(rep(c(0,1),each = 81),times = 17)),
    # hold n2 at first quartile value (contrast low and high values to see interaction)
    n2 = rep(c(1.492e-05), times = 2754),
    # hold ssh_sd at the mean value
    ssh_sd = rep(c(0.02605838), times = 2754),
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
  rbind(
    data.frame(
      ssh = rep(seq(-1,1, 0.025), times = 34), 
      lunar = as.factor(rep(rep(c(0,1), each = 81), times = 17)),
      # same thing with 3rd quartile value
      n2 = rep(c(7.461e-05), times = 2754),
      ssh_sd = rep(c(0.02605838), times = 2754),
      
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
  left_join(
    (high_res %>% 
      select(ptt, species) %>% 
       distinct() %>% 
       mutate(ptt = factor(ptt))),
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
  melt(p_ssh,
       id.vars = c("ssh",
                   "lunar",
                   "n2",
                   "ssh_sd",
                   "ptt",
                   "species"),
       measure.vars = c("1",
                        "2",
                        "3",
                        "4",
                        "5")) %>% 
  mutate(
    lunar = factor(lunar),
    n2 = factor(n2))

rm(p_ssh)

# create a mean line by averaging all individuals 
op1 <- 
  lssh %>% 
  filter(n2 == 1.492e-05) %>% 
  group_by(n2,
           ssh,
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

op2 <- 
  lssh %>% 
  filter(n2 == 7.461e-05) %>% 
  group_by(n2,
           ssh,
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
  
# combine 
opall <- 
  rbind(op1, op2)

rm(op1, op2)

# Set the depth order of the clusters
opall$variable <- 
  factor(opall$variable,
         levels = c(2, 1, 3, 4, 5),
         ordered = TRUE) 

# plot relationships for makos
ggplot(data = lssh) + 
  geom_path(aes(x = ssh, y = value, color = (factor(variable,
                                                    levels = c(2, 1, 3, 4, 5),
                                                    ordered = TRUE) )), alpha = 0.25) +
  geom_line(data = opall,aes(x = ssh, y = value, color = variable), size = 1.2, alpha = 1.2) +
  scale_color_manual(values = c("#FFFF5CFF",
                                "#78CEA3FF",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  facet_grid(species~n2) + 
  labs(color = "Cluster") +
  ylab("Probability") +
  theme_minimal()

## Effect of SSH_SD -----------------------------------------------------------

# vary ssh_sd while holding other predictors constant
p_ssh_sd <- 
  data.frame(
    ssh_sd = rep(seq(0,0.13, 0.001625),times = 34), 
    lunar = as.factor(rep(rep(c(0,1),each = 81),times = 17)),
    # hold n2 at mean
    n2 = rep(c(5.949e-05), times = 2754),
    # hold ssh at the mean value
    ssh = rep(c(-0.1304), times = 2754),
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
lssh <- 
  melt(p_ssh_sd,
       id.vars = c("ssh",
                   "lunar",
                   "n2",
                   "ssh_sd",
                   "ptt",
                   "species"),
       measure.vars = c("1",
                        "2",
                        "3",
                        "4",
                        "5")) %>% 
  mutate(
    lunar = factor(lunar),
    n2 = factor(n2))

rm(p_ssh_sd)

# create a mean line by averaging all individuals 
op <- 
  lssh %>% 
  group_by(ssh_sd,
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
ggplot(data = lssh) + 
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



# industrial series plots -------------------------------------------------

# create an object with series data for the cluster of interest
c1 <- 
  cluster_series %>%
  # update with the cluster of interest
  filter(cluster == 2)

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
    
    ggsave(paste("products/series/cluster2/", full_series[1,1], ".png", sep = ""))
    
    rm(day, sharks, ptt_n, full_series)
  }
  else {
    rm(day, sharks, ptt_n, full_series)
  }
}

