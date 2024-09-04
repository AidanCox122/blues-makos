
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

# Table 2 -----------------------------------------------------------------

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


# Table 3 -----------------------------------------------------------------

# This table will contain performance metrics for our multinom. logit models

source('scripts/analysis/modeling.R')


