## SCRIPT 1: Etuff Setup and Data Formatting

# setup -------------------------------------------------------------------
library(tidyverse)
library(tags2etuff)


# TAD input ------------------------------------------------------------------

# Load the etuff data for blue and mako sharks --
fList <- list.files("data/raw/etuff", full.names = T)
# separate etuffs for blue sharks
blues <- fList[grep('160424', fList)]
# and for makos
makos <- fList[grep('159924', fList)]

# for mako sharks, get on-board calculated time-at-depth summaries
for(i in 1:length(makos)) {
  # read the raw etuff file
  etuff <- read_etuff(makos[i])
  # extract tad data from the etuff 
  tad <- get_tad(etuff)
  # separate date and time data, we only want to know unique dates
  tad <- tad %>% separate(DateTime, into = c("Date", "Time"), sep = "([ ])")
  # pivot data into a tidy format, take the mean of bin values for dates with >1 summary
  tad <- tad %>% pivot_wider(id_cols = Date,
                             names_from = bin,
                             values_from = freq,
                             values_fn = list(freq = mean))
  
  # what is the tag number
  tad$ptt <- etuff$meta$ptt
  # who programmed this tag
  tad$owner <- etuff$meta$person_owner
  # assign output and cleanup
  assign(paste("mtad", i, sep = "_"), tad)
  rm(etuff, tad)
  print(i)
}

mako_raw <- list(
  m1 = mtad_1,
  m2 = mtad_2,
  m3 = mtad_3,
  m4 = mtad_4,
  m5 = mtad_5,
  m6 = mtad_6,
  m7 = mtad_7
)
  
# for blue sharks, get on-board calculated time-at-depth summaries 
for(i in 1:length(blues)) {
  etuff <- read_etuff(blues[i])
  tad <- get_tad(etuff)
  tad <- tad %>% separate(DateTime, into = c("Date", "Time"), sep = "([ ])")
  tad <- tad %>% pivot_wider(id_cols = Date,
                             names_from = bin,
                             values_from = freq,
                             values_fn = list(freq = mean))
  tad$ptt <- etuff$meta$ptt
  tad$owner <- etuff$meta$person_owner
  assign(paste("btad", i, sep = "_"), tad)
  rm(etuff, tad)
  print(i)
}

blue_raw <- 
  list(
    b1 = btad_1,
    b2 = btad_2,
    b3 = btad_3,
    b4 = btad_4,
    b5 = btad_5,
    b6 = btad_6,
    b7 = btad_7,
    b8 = btad_8,
    b9 = btad_9,
    b10 = btad_10,
    b11 = btad_11,
    b12 = btad_12,
    b13 = btad_13
  )


## format bins 

# establish common bins between the different tags

# first make sure that Greg's tags have 14 tad bins and Cam's tags have 12
# makos
mtad_3$TimeAtDepthBin12 <- 0 # Greg's tag
mtad_3$TimeAtDepthBin13 <- 0 # Greg's tag
mtad_3$TimeAtDepthBin14 <- 0 # Greg's tag

mtad_5$TimeAtDepthBin11 <- 0 # Cam's tag
mtad_5$TimeAtDepthBin12 <- 0 # Cam's tag

mtad_6$TimeAtDepthBin12 <- 0 # Cam's tag

mtad_7$TimeAtDepthBin12 <- 0 # Cam's tag

# blues (all Cam's tags)
btad_1$TimeAtDepthBin11 <- 0
btad_1$TimeAtDepthBin12 <- 0

btad_2$TimeAtDepthBin12 <- 0

btad_3$TimeAtDepthBin12 <- 0

btad_7$TimeAtDepthBin12 <- 0

btad_9$TimeAtDepthBin11 <- 0
btad_9$TimeAtDepthBin12 <- 0

btad_11$TimeAtDepthBin12 <- 0

# now combine bins to create the below distribution
bins <- c(0, 10, 50, 100, 200, 300, 400, 500, 2000)
# Camrin's Innate Bins = c(2, 10, 20, 50, 100, 200, 300, 400, 500, 700, 1000, 2000)
# Gregory's Innate Bins = c(5, 10, 25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 2000)

# makos
gregory_mako <- rbind(mtad_1, mtad_2, mtad_3) 
gregory_mako <- 
  gregory_mako %>% 
  transmute(
    Date = Date,
    bin1 = TimeAtDepthBin01 + TimeAtDepthBin02,
    bin2 = TimeAtDepthBin03 + TimeAtDepthBin04,
    bin3 = TimeAtDepthBin05 + TimeAtDepthBin06,
    bin4 = TimeAtDepthBin07 + TimeAtDepthBin08,
    bin5 = TimeAtDepthBin09,
    bin6 = TimeAtDepthBin10,
    bin7 = TimeAtDepthBin11,
    bin8 = TimeAtDepthBin12 + TimeAtDepthBin13 + TimeAtDepthBin14,
    ptt = ptt,
    owner = owner
  )
camrin_mako <- rbind(mtad_4, mtad_5, mtad_6, mtad_7)
camrin_mako <- 
  camrin_mako %>%
  transmute(
    Date = Date,
    bin1 = TimeAtDepthBin01 + TimeAtDepthBin02,
    bin2 = TimeAtDepthBin03 + TimeAtDepthBin04,
    bin3 = TimeAtDepthBin05,
    bin4 = TimeAtDepthBin06,
    bin5 = TimeAtDepthBin07,
    bin6 = TimeAtDepthBin08,
    bin7 = TimeAtDepthBin09,
    bin8 = TimeAtDepthBin10 + TimeAtDepthBin11 + TimeAtDepthBin12,
    ptt = ptt,
    owner = owner
  )

tad_mako <- rbind(gregory_mako, camrin_mako)

rm(gregory_mako, camrin_mako)

# blues 
# all blue sharks from this study were tagged by camrin, combining bins accordingly

tad_blue <- rbind(btad_1, btad_2, btad_3, btad_4, btad_5, btad_6, btad_7,
                   btad_8, btad_9, btad_10, btad_11, btad_12, btad_13)
tad_blue <- 
  tad_blue %>%
  transmute(
    Date = Date,
    bin1 = TimeAtDepthBin01 + TimeAtDepthBin02,
    bin2 = TimeAtDepthBin03 + TimeAtDepthBin04,
    bin3 = TimeAtDepthBin05,
    bin4 = TimeAtDepthBin06,
    bin5 = TimeAtDepthBin07,
    bin6 = TimeAtDepthBin08,
    bin7 = TimeAtDepthBin09,
    bin8 = TimeAtDepthBin10 + TimeAtDepthBin11 + TimeAtDepthBin12,
    ptt = ptt,
    owner = owner
  )

# cleanup
rm(mtad_1, mtad_2, mtad_3, mtad_4, mtad_5, mtad_6, mtad_7, btad_1, btad_2, btad_3, btad_4, btad_5,
   btad_6, btad_7, btad_8, btad_9, btad_10, btad_11, btad_12, btad_13 )

## store these data

# TAD DATA
# mako TAD file
write_csv(tad_mako, 'data/clean/TAD/tad_mako.csv')
# list of individual TAD summaries
write_rds(mako_raw, 'data/clean/TAD/all_mako_tads.rds')

# blue TAD file
write_csv(tad_blue, 'data/clean/TAD/tad_blue.csv')
# list of individual TAD summaries 
write_rds(blue_raw, 'data/clean/TAD/all_blue_tads.rds')

# Track input -------------------------------------------------------------

# extract geolocation data for each individual 

# makos
for(i in 1:length(makos)) {
  # read in the etuff file
  etuff <- read_etuff(makos[i])
  # extract location data from etuff
  track <- get_track(etuff)
  # specify the tag number
  track$ptt <- etuff$meta$ptt
  # separate date and time data
  track$dup <- track$DateTime
  track <- track %>% separate(dup, into = c("Date", "Time"), sep = "([ ])")
  # create a primary key 
  track <- track %>% mutate(
    kode = paste(Date, ptt, sep = "_")
  )
  # assignment and cleanup
  assign(paste("mtrack",i, sep = "_" ), track)
  rm(etuff, track)
  print(i)
}

# combine all location data
mako_tracks <- rbind(mtrack_1, mtrack_2, mtrack_3, mtrack_4,
                     mtrack_5, mtrack_6, mtrack_7)

# add a secondary key for species
mako_tracks$species <- "I.oxyrinchus"

rm(mtrack_1, mtrack_2, mtrack_3, mtrack_4, mtrack_5, mtrack_6, mtrack_7)

# blues
for(i in 1:length(blues)) {
  etuff <- read_etuff(blues[i])
  track <- get_track(etuff)
  track$ptt <- etuff$meta$ptt
  track$dup <- track$DateTime
  track <- track %>% separate(dup, into = c("Date", "Time"), sep = "([ ])")
  track <- track %>% mutate(
    kode = paste(Date, ptt, sep = "_")
  )
  assign(paste("btrack",i, sep = "_" ), track)
  rm(etuff, track)
  print(i)
}

blue_tracks <- rbind(btrack_1, btrack_2, btrack_3, btrack_4, btrack_5, btrack_6,
                     btrack_7, btrack_8, btrack_9, btrack_10, btrack_11,
                     btrack_12, btrack_13)
blue_tracks$species <- "P.glauca"

rm(btrack_1, btrack_2, btrack_3, btrack_4, btrack_5, btrack_6, btrack_7, btrack_8, btrack_9, btrack_10, btrack_11, btrack_12, btrack_13)

# create combo_track object
combo_track <- rbind(blue_tracks, mako_tracks)

shark_tracks <- list(
  mako = mako_tracks,
  blue = blue_tracks
)

## add bathymetry data ----------------------------------------------------------

## load global bathy from file stored on MPG drive
bathy <- raster::raster('data/raw/global_bathy_0.01.nc')

## this is a global grid with pacific-centered coordinates (longitudes 0 to 360)
## raster::rotate converts from 0-360 longitudes to 180 longitudes (atlantic-centered)
bathy <- raster::rotate(bathy)

# create an empty vessel for bathymetry data
bathy_stmp <- 
  vector("list", 
         length = length(unique(combo_track$kode))) %>% 
  setNames(unique(combo_track$kode))

for(x in unique(combo_track$kode)) {
  # subset the data to unique day for unique shark
  subset <- 
    filter(combo_track, kode == x) %>%
    # for sharks with multiple locations per day, just use the first one
    filter(DateTime == min(DateTime))
  
  if(near(subset$latitudeError,0)) {
    # if position error is negligable, do a simple extraction
    bathy_stmp[[x]] <-
      tibble(
        kode = x,
        bathy = raster::extract(bathy, cbind(subset$longitude, subset$latitude))
      )
  } else { 
    # if position error substantial, take mean of depth across possible location
    bathy_stmp[[x]] <- 
      tibble(
        kode = x,
        bathy = subset %>% 
          summarize(xmin = longitude-longitudeError,
                    xmax = longitude+longitudeError,
                    ymin = latitude - latitudeError,
                    ymax = latitude + latitudeError) %>% 
          as.numeric() %>%
          raster::extent() %>% 
          raster::extract(x = bathy, fun=mean, na.rm = TRUE)
      )
  }
}

combo_track <- 
  combo_track %>% 
  # select only one positon for each day to match TAD summary resolution
  group_by(kode) %>% 
  filter(DateTime == min(DateTime)) %>% 
  ungroup() %>% 
  # join tracks w. bathymetry
  left_join(
    bind_rows(bathy_stmp),
    by = 'kode') # %>% mutate(prop.off = 1 - bathy.x/bathy.y) %>% arrange(desc(prop.off))

## store these data

# write the combined track object
write_csv(combo_track, 'data/clean/Tracks/combo_track.csv')

# write a redundant object broken down by species
write_rds(shark_tracks, 'data/clean/Tracks/all_sharks.rds')

write_rds(bathy_stmp, 'data/raw/bathymetry.rds')

# Series input ------------------------------------------------------------------

# makos
for (i in 4:length(makos)) {
  # load the etuff files
  etuff <- read_etuff(makos[i])
  # get the series data
  series <- get_series(etuff)
  # remove incomplete observations
  series <- filter(series, !is.na(depth))
  # add temperature data
  # get temperature data recorded by tag
  temp <- get_pdt(etuff)
  # apply loess smoother to create temperature profile of water column
  interp <- interp_pdt(etuff)
  series <- add_series_temp(series, temp, interp)
  # add daytime or nighttime to each observation
  series$dn <- add_daynight(series, etuff)
  # metadata
  series$species <- "I.oxyrinchus"
  series$ptt <- etuff$meta$ptt
  series$dup <- series$DateTime_local
  series <- series %>% separate(dup, into = c("Date", "Local_Time"), sep = "([ ])")
  series <- series %>% mutate(
    kode = paste(Date, ptt, sep = "_"))
  assign(paste("mseries",i, sep = "_" ), series)
  rm(etuff, series)
  print(i)
}

# mako_7 has 2x the temporal resolution. Let's kick it down a notch
toDelete <- c(seq(0,nrow(mseries_7), 2))
mseries_7.short <- mseries_7[-toDelete,]
rm(toDelete)

mako_series <- rbind(mseries_4, mseries_5, mseries_6, mseries_7.short)

# blues
for (i in 1:length(blues)) {
  etuff <- read_etuff(blues[i])
  series <- get_series(etuff)
  series <- series %>% filter(!is.na(depth))
  # add temperature data
  temp <- get_pdt(etuff)
  interp <- interp_pdt(etuff)
  series <- add_series_temp(series, temp, interp)
  # add day/night
  series$dn <- add_daynight(series, etuff)
  # add metadata
  series$ptt <- etuff$meta$ptt
  series$species <- "P.glauca"
  # separate date from local time
  series$dup <- series$DateTime_local
  series <- series %>% separate(dup, into = c("Date", "Local_Time"), sep = "([ ])")
  series <- series %>% mutate(
    kode = paste(Date, ptt, sep = "_"))
  assign(paste("bseries",i, sep = "_" ), series)
  rm(etuff, series)
  print(i)
}

blue_series <- rbind(bseries_1, bseries_2, bseries_3, bseries_4, bseries_5,
                     bseries_6, bseries_7, bseries_8, bseries_9, bseries_10,
                     bseries_11, bseries_12, bseries_13)

mako_series_raw <- c(mseries_4, mseries_5, mseries_6, mseries_7)
  
shark_series <- 
  list(
  makoSeries = mako_series_raw,
  blueSeries = blue_series
)
combo_series <- rbind(mako_series, blue_series)

# ## add daynight data to combo_series ---------------------------------------
# # another way of adding daynight data
# daynight_stmp <- 
#   vector("list", 
#          length = length(unique(combo_series$ptt))) %>% 
#   setNames(unique(combo_series$ptt))
# 
# for(x in unique(combo_series$ptt)) {
#   file <- fList[grep(paste(x), fList)]
#   etuff <- read_etuff(file)
#   series <- 
#     combo_series %>% 
#     filter(ptt == x)
#   series$dn <- add_daynight(series, etuff)
#   daynight_stmp[[paste(x)]] <-
#     series %>% select(kode, DateTime_local, dn)
# }
# 
# combo_series <- 
#   combo_series %>% 
#   left_join(bind_rows(daynight_stmp),
#             by = c('kode', 'DateTime_local'))

## add bathymetry -----------------------------------------------------------
combo_series <- 
  combo_series %>% 
  left_join(bind_rows(bathy_stmp),
            by = 'kode')

## save these data

# combo series
write_csv(combo_series, 'data/clean/Series/combo_series.csv')

# individual shark series
write_rds(shark_series, 'data/clean/Series/shark_series.rds')


# Header input ------------------------------------------------------------

# makos
for (i in 4:length(makos)) {
  hdr <- get_header(makos[i])
  hdr <- hdr %>% mutate(friendly_name = ifelse(ptt == 206771, 'Muscles', friendly_name))
  hdr <- dplyr::select(hdr, c(instrument_type, length_capture,
                length_type_capture, length_unit_capture, lifestage_capture,
                person_owner, ptt, friendly_name, sex, time_coverage_start,
                time_coverage_end, geospatial_lat_start,
                geospatial_lon_start, geospatial_lat_end,
                geospatial_lon_end, taxonomic_serial_number, waypoints_source))
  assign(paste("mHdr",i, sep = "_" ), hdr)
  rm(hdr)
  print(i)
}

# blues
for (i in 1:length(blues)) {
  hdr <- get_header(blues[i])
  
  hdr <- dplyr::select(hdr, c(instrument_type, length_capture,
                              length_type_capture, length_unit_capture, lifestage_capture,
                              person_owner, ptt, friendly_name, sex, time_coverage_start,
                              time_coverage_end, geospatial_lat_start,
                              geospatial_lon_start, geospatial_lat_end,
                              geospatial_lon_end, taxonomic_serial_number, waypoints_source)) 
  assign(paste("bHdr",i, sep = "_" ), hdr)
  rm(hdr)
  
  # select columns of interest
  
  print(i)}


combo_hdr <- rbind(mHdr_4, mHdr_5, mHdr_6, mHdr_7,
                   bHdr_1, bHdr_2, bHdr_3, bHdr_4,
                   bHdr_5, bHdr_6, bHdr_7, bHdr_8,
                   bHdr_9, bHdr_10, bHdr_11, bHdr_12,
                   bHdr_13)

# store the data
write_csv(combo_hdr, 'data/clean/combo_hdr.csv')


