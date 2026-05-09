
# setup -------------------------------------------------------------------
library(tidyverse)

# get TAD data
tad_mako <- read_csv('~/Desktop/blues-makos/data/clean/TAD/tad_mako.csv')
tad_blue <- read_csv('~/Desktop/blues-makos/data/clean/TAD/tad_blue.csv')

tad_combo <- rbind(tad_mako, tad_blue)

# get series data
combo_series <- read_csv('~/Desktop/blues-makos/data/clean/Series/combo_series.csv')

# headers
combo_header <- read_csv('~/Desktop/blues-makos/data/clean/combo_hdr.csv')

# Table 1 -----------------------------------------------------------------

# calculate tad data coverage for each shark
# leaving this out because we don't use the innate TAD data in this study
# tad_coverage <- 
#   tad_combo %>% 
#   group_by(ptt) %>%
#   summarize(TAD = n()) %>%
#   rename(TagID = ptt)

# calculate series data coverage
series_coverage <- combo_series %>%
  group_by(ptt) %>%
  summarize(Series = unique(Date) %>% length()) %>% 
  rename(TagID = ptt)

meta <- combo_header %>%
  select(TagID = ptt,
         Species = taxonomic_serial_number,
         TagDate = time_coverage_start,
         PopupDate = time_coverage_end,
         TagLat = geospatial_lat_start,
         TagLon = geospatial_lon_start,
         PopupLat = geospatial_lat_end,
         PopupLon = geospatial_lon_end) %>% 
  mutate(DaysAtLiberty = PopupDate - TagDate)

table2 <- 
  meta %>%
  inner_join(series_coverage, by = 'TagID') %>% 
  mutate(ProportionSeries = Series/as.numeric(DaysAtLiberty)) %>% 
  dplyr::select(-c(Series))


# Table 2 -----------------------------------------------------------------

# This table will contain performance metrics for our multinom. logit models

source('scripts/analysis/general_analysis.R')

# median diel depth
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
  group_by(cluster, species, dn) %>% 
  summarize(
    Med.depth = median(depth)) %>% View()

# standard deviation of depth use
high_res %>% 
  # rename the clusters for easy reading
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>%
  # calculate the average % time spent in each depth bin
  group_by(cluster, species) %>% 
  summarize(
    mean.SD.day = mean(d.sd) %>% round(1),
    mean.SD.night = mean(n.sd) %>% round(1)) %>% View()

# percent time in mesopelagic
high_res <- 
  high_res %>% 
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3'))

high_res %>% 
  group_by(cluster, species) %>% count()

PercMeso <- 
  high_res %>% 
  # rename the clusters for easy reading
  mutate(cluster = case_when(
    cluster == 1 ~ 'EPI 2',
    cluster == 2 ~ 'DVM 1',
    cluster == 3 ~ 'EPI 1',
    cluster == 4 ~ 'DVM 2',
    cluster == 5 ~ 'DVM 3')) %>%
  # create columns for % time in each zone
  group_by(kode, species, cluster) %>%
  summarize(
    # daytime depth use
    'day_epi' = sum(c(d.b1, d.b2, d.b3, d.b4)),
    'day_meso' = sum(c(d.b5, d.b6, d.b7, d.b8)),
    'night_epi' = sum(c(n.b1, n.b2, n.b3, n.b4)),
    'night_meso' = sum(c(n.b5, n.b6, n.b7, n.b8))) %>% 
  # calculate the average % time spent in each depth bin
  group_by(cluster, species) %>% 
  summarise(
    # daytime habitat
    # 'day_mean_epi' = mean(day_epi),
    # 'day_sd_epi' = sd(day_epi),
    'day_mean_meso' = mean(day_meso) %>% round(1),
    # 'day_sd_meso' = sd(day_meso),
    # nightly habitat
    # 'night_mean_epi' = mean(night_epi),
    # 'night_sd_epi' = sd(night_epi),
    'night_mean_meso' = mean(night_meso) %>% round(1),
    # 'night_sd_meso' = sd(night_meso),
    # 'overall_mean_meso' = mean(c(day_meso, night_meso)),
    # 'overall_sd_meso' = sd(c(day_meso, night_meso))
    )

# latitude 
complete_series_0.5 %>%
  dplyr::select(cluster, species.x, latitude, longitude) %>% 
  group_by(cluster, species.x) %>% 
  summarize(
    mean.Lat = mean(latitude),
    sd.Lat = sd(latitude))

# distance from shelf
dist_shelf %>% #filter(cluster %in% c('DVM 2', 'DVM 3')) %>% 
  group_by(cluster, species) %>%
  summarize(
    mean.dist.km = mean(distance)/1000,
    sd.dist.km = sd(distance)/1000) %>% View()

