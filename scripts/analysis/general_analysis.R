
# notes -------------------------------------------------------------------

# for instructions on running an ANOVA test in R, see: https://statsandr.com/blog/anova-in-r/
# instructions on repeated measure ANOVA IN R: https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#google_vignette
# instructions on kruskal-walis test in R: https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/#google_vignette
# instructions on friedman test in R: https://www.datanovia.com/en/lessons/friedman-test-in-r/
# instructions on welch's t-test in R: https://whitlockschluter3e.zoology.ubc.ca/RLabs/R_tutorial_Comparing_means_of_2_groups.html#welch’s_t-test

# setup -------------------------------------------------------------------
library(tidyverse)
library(tags2etuff)
library(multcomp)
library(sf)
library(ggpubr)
library(rstatix) # for ANOVA functions
# library(car) # leveneTest

fList <- list.files("data/raw/etuff", full.names = T)
blues <- fList[grep('160424', fList)]
makos <- fList[grep('159924', fList)]

combo_series <- read_csv('data/clean/Series/combo_series.csv')

combo_track <- read_csv('data/clean/Tracks/combo_track.csv')

bathy_stmp <- read_rds('data/raw/bathymetry.rds')

source('scripts/analysis/hierarchical_clustering-MCA.R')

#define standard error of mean function
std.error <- function(x) sd(x)/sqrt(length(x))

# general movement summaries ----------------------------------------------

# find number of days occuring over the continental shelf
filter(combo_series, bathy > -1000) %>%
  pull(kode) %>% 
  unique() # 358 tracking days occured over the shelf (this is how many we lose) 

# filter out tracking days with less than 50% of their data
## 86400 daily seconds / 150 seconds per sample = 576

## create an object with only 1/2 complete series records-----
complete_series_0.5 <- 
  combo_series %>% 
  # remove data from over the continental shelf
  filter(bathy <= -1000) %>% 
  group_by(ptt) %>% 
  count(Date) %>% # there are 1634 total tracking days with series data
  filter(n >= 288) %>% # there are 1249 tracking days with more than 50% of their series data
  right_join(combo_series, by = c("Date", "ptt")) %>% 
  filter(!is.na(n)) %>% 
  ungroup() %>% 
  # add latitude and longitude data
  left_join(
    combo_track %>% 
      group_by(kode) %>% 
      filter(DateTime == min(DateTime)),
    by = 'kode') %>% 
  # add cluster 
  left_join(
    clust_stamp2,
    by = 'kode') %>% 
  # convert cluster to a factor
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3'),
    cluster = factor(cluster, levels = c(
      'DVM 1',
      'DVM 2',
      'DVM 3',
      'EPI 1',
      'EPI 2'))) %>% 
  # remove any observations without clusters
  filter(!is.na(cluster)) #85633 observations removed

## average maximum depth by species----
complete_series_0.5 %>%
  group_by(species, ptt, kode) %>% # maximum across individual days
  summarize(max = max(depth)) %>%
  group_by(species, ptt) %>% 
  summarize(all_max = max(max), mean = mean(max)) %>% 
  group_by(species) %>%
  summarize(all_max = max(all_max), grandMean = mean(mean), StErr = std.error(mean))

## Differences in median daytime depth ----

# identify extreme outliers
complete_series_0.5 %>% 
  group_by(cluster, kode, dn) %>%
  summarize(med.depth = median(depth)) %>% 
  ungroup() %>% 
  filter(dn == 'd') %>% 
  identify_outliers(med.depth) # there are nine extreme outliers

# test homoscedasticity
complete_series_0.5 %>% 
  group_by(cluster, kode, dn) %>%
  summarize(med.depth = median(depth)) %>% 
  ungroup() %>% 
  filter(dn == 'd') %>% 
  car::leveneTest(med.depth~cluster, data = .) # p > 0.05            


oneW_depth_aov <- 
  complete_series_0.5 %>% 
  group_by(cluster, kode, dn) %>%
  summarize(med.depth = median(depth)) %>% 
  ungroup() %>% 
  filter(dn == 'd') %>%
  aov(med.depth ~ cluster,
      data = .)

summary(oneW_depth_aov)

one.way.depth <- TukeyHSD(oneW_depth_aov)

# run the non-parametric alternative
complete_series_0.5 %>% 
  group_by(cluster, kode, dn) %>%
  summarize(med.depth = median(depth)) %>% 
  ungroup() %>% 
  filter(dn == 'd') %>% 
  rstatix::kruskal_test(med.depth ~ cluster,
                        data = .) # there is a significant difference

# difference in daytime median depth between species 

complete_series_0.5 %>% 
  group_by(species, kode, dn) %>%
  summarize(med.depth = median(depth)) %>% 
  ungroup() %>% 
  filter(dn == 'd') %>% 
  t.test(med.depth ~ species, data = ., var.equal = FALSE) # p = 0.65


## Differences in standard deviation of depth use --------------------------

# identify any extreme outliers
complete_series_0.5 %>% 
  # recreate d.sd from time-series data
  group_by(species, cluster, kode, dn) %>%
  summarize(sd.depth = sd(depth)) %>% 
  ungroup() %>% 
  # select only daytime values
  filter(dn == 'd') %>% 
  group_by(species) %>% 
  identify_outliers(sd.depth) # there are no extreme outliers

# test for homoscedasticity
complete_series_0.5 %>% 
  # recreate d.sd from time-series data
  group_by(species, cluster, kode, dn) %>%
  summarize(sd.depth = sd(depth)) %>% 
  ungroup() %>% 
  # select only daytime values
  filter(dn == 'd') %>%
  leveneTest(sd.depth~species, data = .) # p < 0.005 so variance equal

complete_series_0.5 %>% 
  # recreate d.sd from time-series data
  group_by(species, kode, dn) %>%
  summarize(sd.depth = sd(depth)) %>% 
  ungroup() %>% 
  # select only daytime values
  filter(dn == 'd') %>% 
  t.test(sd.depth ~ species, data = ., var.equal = FALSE) # p = 0.05579


## Difference in latitude distribution----
# where do most of the track locations fall?
# what is the range of locations in each direction?
complete_series_0.5 %>%
  dplyr::select(latitude, longitude) %>% 
  summary()

# identify any extreme outliers
complete_series_0.5 %>% 
  dplyr::select(species, kode, latitude, cluster) %>% 
  group_by(cluster) %>% 
  identify_outliers(latitude) %>% View()

# test for homoscedasticity
complete_series_0.5 %>% 
  dplyr::select(species, kode, latitude, cluster) %>% 
  leveneTest(latitude~species, data = .) # p < 0.001 there is unequal variance between species

lat_aov <- aov(latitude ~ cluster,
                 data = complete_series_0.5 %>% 
                 group_by(species, cluster, kode) %>%
                 summarize(latitude = mean(latitude)))

summary(lat_aov)

# non-parametric version
complete_series_0.5 %>% 
  group_by(cluster, kode) %>%
  summarize(latitude = mean(latitude)) %>% 
  ungroup() %>% 
  rstatix::kruskal_test(latitude ~ cluster,
                        data = .) # there is a highly signif. difference btwn. clusters

lat.tukey <- TukeyHSD(lat_aov)


## Difference in distance from shelf by cluster ----------------------------

# create a crs object for WGS84 Azimuth Equidistant projection
my_crs <- 
st_crs('PROJCRS["World_Azimuthal_Equidistant",
    BASEGEOGCRS["WGS 84",
        DATUM["World Geodetic System 1984",
            ELLIPSOID["WGS 84",6378137,298.257223563,
                LENGTHUNIT["metre",1]]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["Degree",0.0174532925199433]]],
    CONVERSION["World_Azimuthal_Equidistant",
        METHOD["Modified Azimuthal Equidistant",
            ID["EPSG",9832]],
        PARAMETER["Latitude of natural origin",0,
            ANGLEUNIT["Degree",0.0174532925199433],
            ID["EPSG",8801]],
        PARAMETER["Longitude of natural origin",0,
            ANGLEUNIT["Degree",0.0174532925199433],
            ID["EPSG",8802]],
        PARAMETER["False easting",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8806]],
        PARAMETER["False northing",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8807]]],
    CS[Cartesian,2],
        AXIS["(E)",east,
            ORDER[1],
            LENGTHUNIT["metre",1]],
        AXIS["(N)",north,
            ORDER[2],
            LENGTHUNIT["metre",1]],
    USAGE[
        SCOPE["Not known."],
        AREA["World."],
        BBOX[-90,-180,90,180]],
    ID["ESRI",54032]]')
# read in WGS84 shapefiles and reproject to NAD83 Azimuth equidistant
## isobath 
isobath1000 <-
  st_read('/Users/aidansmacpro/Documents/GIS/blues-makos/isobath1000.shp')

projected_isobath1000 <- 
  isobath1000 %>% 
  st_transform(crs = my_crs) #26919 (NAD83 UTM Z19N does not cover full extent of study area)

# st_write(projected_isobath1000, '/Users/aidansmacpro/Documents/GIS/blues-makos/projected_isobath1000.shp')

## cluster locations
cluster_locs <- 
  read_csv('/Users/aidansmacpro/Desktop/blues-makos/data/processed/clustered_high_res_tad.csv') %>%
  dplyr::select("ptt",       "Date",      "species",   "kode",      "bathy",    
                "latitude",  "longitude",  "cluster") %>% 
  mutate(distance = NA)

## convert locs to SF object
locs_sf <- st_as_sf(cluster_locs,coords = c("longitude","latitude"),
                    # points were colleted using WGS84
                    crs =4326)  

projected_cluster_locs <- 
  locs_sf %>% 
  st_transform(crs = my_crs)

# st_write(projected_cluster_locs, '/Users/aidansmacpro/Documents/GIS/blues-makos/projected_cluster_locs.shp')

# load distance from 1000m isobath for each cluster obs.

dist_1000 <- 
  # unprojected
  # read_csv('/Users/aidansmacpro/Desktop/blues-makos/data/processed/distance_to_1000.csv') %>% 
  # projected
  read_csv('/Users/aidansmacpro/Desktop/blues-makos/data/processed/projected_dist_to_1000.csv') %>% 
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
      ordered = T))

# load distance from shelf for each cluster obs.
dist_shelf <- 
  # unprojected
  # read_csv('/Users/aidansmacpro/Desktop/blues-makos/data/processed/distance_to_shelf.csv') %>% 
  # projected
  read_csv('/Users/aidansmacpro/Desktop/blues-makos/data/processed/projected_dist_to_shelf.csv') %>%
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
      ordered = T))

# test for difference in distance to 1000m isobath by clusters
aov_1000 <- aov(distance ~ cluster,
                    data = dist_1000)

summary(aov_1000) # p < 0.001

# test if data met anova assumptions
## equal variance
dist_1000 %>% 
  car::leveneTest(distance~cluster, data = .) # p < 0.001, unequal var.

## normal distribution of variance
car::qqPlot(aov_1000$residuals,
            id = FALSE # id = FALSE to remove point identification
) # our data seem to have a right skew towards farther distances

## no extreme outliers
dist_1000 %>% 
  group_by(cluster) %>% 
  identify_outliers(distance) %>% View()

# non-parametric version
kw.1000 <- 
  dist_1000 %>% 
  rstatix::kruskal_test(distance ~ cluster,
                        data = .) # there is a highly signif. difference btwn. clusters

summary(kw.1000)

# post hoc tests
tukey.1000 <- glht(aov_1000, linfct = mcp(cluster = "Tukey"))

summary(tukey.1000)

# test for difference in distance to shelf by cluster
aov_shelf <- aov(distance ~ cluster,
                data = dist_shelf)

summary(aov_shelf) # p < 0.001

# test if data met anova assumptions
## equal variance
dist_shelf %>% 
  car::leveneTest(distance~cluster, data = .) # p < 0.001, unequal var.

## normal distribution of variance
car::qqPlot(aov_shelf$residuals,
            id = FALSE # id = FALSE to remove point identification
) # our data seem to have a right skew towards farther distances

## no extreme outliers
dist_shelf %>% 
  group_by(cluster) %>% 
  identify_outliers(distance) %>% View()

# non-parametric version
dist_shelf %>% 
  rstatix::kruskal_test(distance ~ cluster,
                        data = .) # there is a highly signif. difference btwn. clusters

summary(kw.1000)

# post hoc tests
tukey.1000 <- glht(aov_1000, linfct = mcp(cluster = "Tukey"))

summary(tukey.1000)

tukey.shelf <- glht(aov_shelf, linfct = mcp(cluster = "Tukey"))

summary(tukey.shelf) # All clusters p < 0.001 except EPI2 - DVM 1

# visualize the differences
ggplot() +
  # closest to farthest order of distance to 1000: DVM1, EPI2, EPI1, DVM3, DVM2
  geom_boxplot(data = dist_1000, aes(x = cluster, y = distance, fill = cluster)) +
  # same pattern but more exaggerated 
  # geom_boxplot(data = dist_shelf, aes(x = cluster, y = distance, fill = cluster)) +
  scale_color_manual(values = c("yellow",
                                "yellow3",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF")) +
  coord_flip() +
  theme_minimal()

## temperature distribution ----
complete_series_0.5 %>%
  filter(!is.na(temperature) & species == "I.oxyrinchus") %>%
  pull(temperature) %>%
  summary()

complete_series_0.5 %>%
  filter(!is.na(temperature) & species == "P.glauca") %>%
  pull(temperature) %>%
  summary()


# Difference in Temperature @ Median/Q3 Depth --------------------------------
# find temperature at MEDIAN
temp_at_median <- 
  complete_series_0.5 %>%
  dplyr::select(species, cluster, kode, dn, depth, temperature, latitude, longitude) %>% 
  filter(!is.na(temperature) & dn == 'd') %>% 
  left_join(
    # find the median depth for each species on each day
    complete_series_0.5 %>% 
      group_by(cluster, kode, dn) %>%
      summarize(med.depth = median(depth)) %>% 
      ungroup() %>% 
      filter(dn == 'd'),
    by = 'kode') %>% 
  # identify median depth
  filter(depth == med.depth) %>%  # 443 unique kodes
  # there are 14,000 entries because sharks may have visited median depth as many as 256 times
  group_by(kode, cluster.x) %>% 
  summarize(mean.temp = mean(temperature)) %>% 
  ungroup()

# test assumptions
## equal variance
temp_at_median %>% 
  car::leveneTest(mean.temp~cluster.x, data = .) # p < 0.001, unequal var.

## no extreme outliers
temp_at_median %>%
  group_by(cluster.x) %>% 
  identify_outliers(mean.temp) # no extreme outliers

temp_at_median %>%
  kruskal_test(mean.temp ~ cluster.x, data = .) # p < 0.001

temp_at_median %>% 
  ggplot() +
  geom_boxplot(aes(x = cluster.x, y = mean.temp, color = cluster.x)) +
  scale_color_manual(values = c("#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF",
                                "yellow",
                                "yellow3")) +
  theme_classic()

## Find temperature at Q3
temp_at_Q3 <- 
  complete_series_0.5 %>%
  dplyr::select(species, cluster, kode, dn, depth, temperature, latitude, longitude) %>% 
  filter(!is.na(temperature) & dn == 'd') %>% 
  left_join(
    # find the median depth for each species on each day
    complete_series_0.5 %>% 
      group_by(cluster, kode, dn) %>%
      summarize(q3.depth = quantile(depth, 0.75)) %>% 
      ungroup() %>% 
      filter(dn == 'd'),
    by = 'kode') %>% 
  # identify median depth
  filter(depth == q3.depth) %>%  # 443 unique kodes
  # there are 9,972 entries because sharks may have visited median depth as many as 256 times
  group_by(kode, cluster.x) %>% 
  summarize(mean.temp = mean(temperature)) %>% 
  ungroup()

temp_at_Q3 %>%
  kruskal_test(mean.temp ~ cluster.x, data = .) # p < 0.001

temp_at_Q3 %>% 
  ggplot() +
  geom_boxplot(aes(x = cluster.x, y = mean.temp, color = cluster.x)) +
  scale_color_manual(values = c("#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF",
                                "yellow",
                                "yellow3")) +
  theme_classic()

## diel patterns ----

# calculate the percentage of time in the epipelagic
complete_series_0.5 %>% 
  mutate(epipelagic = if_else(
    depth <= 200,
    1,
    0)) %>% # View()
  # specify during day 'd' or night 'n' periods
  filter(dn == 'n') %>% 
  group_by(kode, ptt, species) %>% 
  summarize(n = n(),
            propAbove200 = sum(epipelagic)/n,
            .groups = 'drop') %>% 
  group_by(ptt, species) %>% 
  summarize(meanEpiProp = mean(propAbove200)) %>% 
  group_by(species) %>% 
  summarize(propEpi = mean(meanEpiProp),
            StDev = sd(meanEpiProp),
            minEpi = min(meanEpiProp),
            maxEpi = max(meanEpiProp))

# calculate the percentage of time in the mesopelagic
complete_series_0.5 %>% 
  mutate(mesopelagic = if_else(
    depth > 200,
    1,
    0)) %>% # View()
  # specify during day 'd' or night 'n' periods
  filter(dn == 'd') %>% 
  group_by(kode, ptt, species) %>% 
  summarize(n = n(),
            propBelow200 = sum(mesopelagic)/n,
            .groups = 'drop') %>% 
  group_by(ptt, species) %>% 
  summarize(meanMesoProp = mean(propBelow200)) %>% 
  group_by(species) %>% 
  summarize(propMeso = mean(meanMesoProp),
            StDev = sd(meanMesoProp),
            minMeso = min(meanMesoProp),
            maxMeso = max(meanMesoProp))

# binned time-at-depth breakdown 
complete_series_0.5 %>% 
  filter(dn == 'd') %>% 
  group_by(kode, species) %>% 
  summarise(
    b1 = sum((depth < 10)),
    b2 = sum((depth >= 10 & depth <50)),
    b3 = sum((depth >= 50 & depth<100)),
    b4 = sum((depth >= 100 & depth<200)),
    b5 = sum((depth >= 200 & depth<300)),
    b6 = sum((depth >= 300 & depth<400)),
    b7 = sum((depth >= 400 & depth<500)),
    b8 = sum((depth >= 500)),
    sd = sd(depth),
    total = sum(b1,b2,b3,b4,b5,b6,b7,b8),
    .groups = 'drop') %>% # View()
  transmute(
    kode = kode,
    species = species,
    b1 = (b1/total) * 100,
    b2 = (b2/total) * 100,
    b3 = (b3/total) * 100,
    b4 = (b4/total) * 100,
    b5 = (b5/total) * 100,
    b6 = (b6/total) * 100,
    b7 = (b7/total) * 100,
    b8 = (b8/total) * 100
  ) %>% 
  group_by(species) %>% 
  summarize(
    m.b1 = mean(b1),
    m.b2 = mean(b2),
    m.b3 = mean(b3),
    m.b4 = mean(b4),
    m.b5 = mean(b5),
    m.b6 = mean(b6),
    m.b7 = mean(b7),
    m.b8 = mean(b8),
    sd.b1 = sd(b1),
    sd.b2 = sd(b2),
    sd.b3 = sd(b3),
    sd.b4 = sd(b4),
    sd.b5 = sd(b5),
    sd.b6 = sd(b6),
    sd.b7 = sd(b7),
    sd.b8 = sd(b8)) # %>% View


## Differences in bathymetry -----------------------------------------------


