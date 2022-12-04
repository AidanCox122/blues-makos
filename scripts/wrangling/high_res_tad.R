
# setups ------------------------------------------------------------------
library(tidyverse)
library(tags2etuff)

fList <- list.files("data/raw/etuff", full.names = T)
blues <- fList[grep('160424', fList)]
makos <- fList[grep('159924', fList)]

combo_series <- read_csv('data/clean/Series/combo_series.csv')

combo_track <- read_csv('data/clean/Tracks/combo_track.csv')

bathy_stmp <- read_rds('data/raw/bathymetry.rds')

# High Resolution TAD -----------------------------------------------------
# this document will create daily time-at-depth summaries based on high-resolution
# series datasets for each of our tagged sharks

# filter out tracking days with less than 60% of their data
high_res <- combo_series %>%
  # group by individual
  group_by(ptt) %>%
  # count number of records from each date (576 is max. possible total)
  count(Date) %>%
  # there are 675 (928) tracking days with more than 75% (60%) of their series data
  filter(n >= 345) %>% 
  left_join(combo_series, by = c("Date", "ptt")) %>% 
  ungroup()

# count the net observations occur in each depth bin for day and night time
high_res <- high_res %>% group_by(ptt, Date, dn) 
daily.count <- high_res %>% count(Date)
high_res <- high_res %>% summarise(
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
  .groups = 'drop'
)

# check that sum of rows matches daily.count values - 
## this should produce a tibble of length zero
high_res %>% dplyr::select(ptt, Date, dn, total) %>% 
  left_join(daily.count,
            by = c('ptt', 'Date', 'dn')) %>% 
  filter(n != total)

# separate day and night time observations for each tracking day
daytime <- high_res %>% filter(dn == "d")
nightime <- high_res %>% filter(dn == "n")

# convert net counts into percentages within each depth bin
daytime <- daytime %>% transmute(
  ptt = ptt,
  Date = Date,
  d.b1 = (b1/total) * 100,
  d.b2 = (b2/total) * 100,
  d.b3 = (b3/total) * 100,
  d.b4 = (b4/total) * 100,
  d.b5 = (b5/total) * 100,
  d.b6 = (b6/total) * 100,
  d.b7 = (b7/total) * 100,
  d.b8 = (b8/total) * 100,
  d.sd = sd
)

nightime <- nightime %>% transmute(
  ptt = ptt,
  Date = Date,
  n.b1 = (b1/total) * 100,
  n.b2 = (b2/total) * 100,
  n.b3 = (b3/total) * 100,
  n.b4 = (b4/total) * 100,
  n.b5 = (b5/total) * 100,
  n.b6 = (b6/total) * 100,
  n.b7 = (b7/total) * 100,
  n.b8 = (b8/total) * 100,
  n.sd = sd
)

high_res <- inner_join(daytime, nightime, by = c("ptt", "Date"))
rm(daytime, nightime)

# now create a species stamp to apply to this object
species_stmp <- data.frame(
  ptt = c(141257, 163098, 163096, 141254, 141256, 141259, 106754, 133016, 133017, 133018, 133021, 141247, 141255, 163097, 141258, 154096, 78680,  78682,  78683, 206771),
  species = c("I.oxyrinchus", "I.oxyrinchus", "I.oxyrinchus", "P.glauca", "P.glauca", "P.glauca", "P.glauca", "P.glauca", "P.glauca", "P.glauca", "P.glauca", "P.glauca", "P.glauca", "P.glauca", "P.glauca", "P.glauca", "I.oxyrinchus", "I.oxyrinchus", "I.oxyrinchus", "I.oxyrinchus")
)
high_res <- left_join(high_res, species_stmp, by = c("ptt"))

# add kode to high_res
high_res <- high_res %>% mutate(
  kode = paste(Date, ptt, sep = "_")
)

# add depth over location of each tracking day to high-resolution
high_res <- high_res %>% 
  left_join(bathy_stmp %>% 
              bind_rows(), by = "kode")


# remove records which occur over <1000m
high_res <- 
  high_res %>% 
  filter(bathy <= -1000)

