
# setup -------------------------------------------------------------------
library(tidyverse)
library(tags2etuff)

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

## differences in median daytime depth ----

oneW_depth_aov <- aov(med.depth ~ species,
               data = complete_series_0.5 %>% 
                 group_by(species, cluster, kode) %>%
                 summarize(med.depth = median(depth)))

summary(oneW_depth_aov)

one.way.depth <- TukeyHSD(oneW_depth_aov)

# difference in daytime depth between species by cluster
twoW_depth_aov <- aov(med.depth ~ species + cluster,
                 data = complete_series_0.5 %>% 
                   group_by(species, cluster, kode) %>%
                   summarize(med.depth = median(depth)))

summary(twoW_depth_aov)

depth.tukey <- TukeyHSD(twoWdepth_aov)


## Differences in standard deviation of depth use --------------------------

oneW_dSD_aov <- aov(sd.depth ~ species,
                      data = complete_series_0.5 %>% 
                        group_by(species, cluster, kode, dn) %>%
                        summarize(sd.depth = sd(depth)) %>% 
                      # select only daytime values
                      filter(dn == 'd'))

summary(oneW_dSD_aov)

one.way.dSD <- TukeyHSD(oneW_dSD_aov)

## location distribution----
# where do most of the track locations fall?
# what is the range of locations in each direction?
complete_series_0.5 %>%
  dplyr::select(latitude, longitude) %>% 
  summary()

lat_aov <- aov(latitude ~ cluster,
                 data = complete_series_0.5 %>% 
                 group_by(species, cluster, kode) %>%
                 summarize(latitude = mean(latitude)))

summary(lat_aov)

lat.tukey <- TukeyHSD(lat_aov)

## temperature distribution ----
complete_series_0.5 %>%
  filter(!is.na(temperature) & species == "I.oxyrinchus") %>%
  pull(temperature) %>%
  summary()

complete_series_0.5 %>%
  filter(!is.na(temperature) & species == "P.glauca") %>%
  pull(temperature) %>%
  summary()

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


