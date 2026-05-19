### THIS SCRIPT CALCULATES STATS. FOR INDIV. DIVES FROM HIGH-RES. TIMESERIES DATA FOR ALL TAGGED SHARKS
# setup -------------------------------------------------------------------

library(tidyverse)
devtools::load_all('../../Documents/Ecology/Marine/tags2etuff')

cluster_stamp <- 
  read_csv('data/clean/cluster_stamp.csv') %>% 
  mutate(
    cluster_n = as.numeric(cluster),
    cluster = case_when(
      cluster_n == 1 ~ 'EPI 2',
      cluster_n == 2 ~ 'DVM 1',
      cluster_n == 3 ~ 'EPI 1',
      cluster_n == 4 ~ 'DVM 2',
      cluster_n == 5 ~ 'DVM 3')
  )

# create a vector of satellite tag identifiers
ptts <- c(141257, 163096, 163098, '20P1624', #206771,this is the ptt I have been using
          141254, 141256, 141259, 106754, 133016,
          133017, 133018, 133021, 141247, 141255, 141258, 154096, 163097)

# create a list of series files to read in
flist <- list.files('data/clean/Series/all_shark_series', full.names = T)

# read in series data as a list
series_list <- 
  # create a list
  vector(mode = 'list', length = 17) %>%
  # name each object for a ptt
  set_names(ptts)


# get series data ---------------------------------------------------------

for(x in names(series_list)) {
  # find the file associated with x
  file <- flist[grep(x, flist)]
  
  # read the file and store in series_list
  series_list[[x]] <-
    read_csv(file)
  
  # update user
  print(
    glue::glue('Done with {x}')
  )
  
  # cleanup
  rm(x)
}


# get dives ---------------------------------------------------------------

for(x in names(series_list)) {
  
  # extract data from series_list
  series <- series_list[[x]]
  
  # rearrange series based on local time reported by tag x
  series <- series[order(series$DateTime_local),]
  
  # create a diveMove TDR object from series
  std <- diveMove::createTDR(series$DateTime_local, series$depth, #concurrentData = series[,-c(12,4)],
                             file='xxx') # file argument required but doesnt seem to matter in this case
  
  # attempt to identify individual 'dives' within TDR obj: adjust depth of dive threshold
  try.stdc <- try(stdc <- diveMove::calibrateDepth(std, dry.thr=0, wet.cond=rep(TRUE, length.out=nrow(series)),
                                         dive.thr = 200, zoc.method='offset', offset=0, dive.model='smooth.spline'), silent=TRUE)
  
  # if an error occurred try again after offsetting all depths down by 10m
  if (class(try.stdc) == 'try-error'){
    stdc <- diveMove::calibrateDepth(std, dry.thr=0, wet.cond=rep(TRUE, length.out=nrow(series)),
                                     dive.thr = 200, zoc.method='offset', offset=10, dive.model='smooth.spline')
  } else{
    # if try.stdc worked initially, save it as a new object
    stdc <- try.stdc
  }
  
  # cleanup
  rm(try.stdc)
  
  ## visualize dives
  #diveMove::plotTDR(std)
  #plotTDR(stdc, diveNo=30)
  
  # get dive stats 
  stdsumm <- try(diveMove::diveStats(stdc), silent=TRUE)
  
  # match dives to series obs. 
  if (class(stdsumm) != 'try-error'){
    
    # critera req'd such as max dive time is 24 hours (or less?) and min dive time is 30 mins?
    # diveIdx <- which(stdsumm$maxdep >= 200)# & stdsumm$divetim > 1800 & stdsumm$divetim <= (24*3600))
    
    # create an index of dives and a reference index in series
    diveIdx <- 1:nrow(stdsumm)
    series$rownum <- 1:nrow(series)
    
    # create a column to store the final index
    series$mesoBout <- NA
    
    ## find the series record corresponding to start and end of each dive
    for (ii in diveIdx){
      start.idx <- which.min(as.numeric((series$DateTime_local - stdsumm$begdesc[ii])) ^ 2)
      if (!is.finite(start.idx)) start.idx <- 1
      end.idx <- which.min(as.numeric((series$DateTime_local - (stdsumm$begdesc[ii] + stdsumm$divetim[ii]))) ^ 2)
      # if the above step fails, assign end time to last obs. in series
      if (!is.finite(end.idx)) end.idx <- nrow(series)
      
      # if bout starts or ends at deployment edge, discard
      if(start.idx == 1 || end.idx == nrow(series)) next
      
      # if dive starts or ends at NA record (or one adjacent), discard the bout
      if (
        is.na(series$depth[start.idx]) ||
        is.na(series$depth[end.idx]) ||
        is.na(series$depth[start.idx - 1]) ||
        is.na(series$depth[end.idx + 1])) next
      
      # save series index
      series.idx <- c(start.idx:end.idx)
      
      # if bout ends below 200m, and next record is above 200m, modify dive end
      if (series$depth[end.idx] > 200) end.idx <- min(which(series$depth <= 200 & series$rownum > series.idx[length(series.idx)]))      
      
      # add a test
      # if bout ends below 200m, and next record is above 200m, modify dive end
      # if(series$depth[end.idx] > 200 & series$depth[end.idx+1] > 200 & !is.na(series$depth[end.idx+1])) ii <- paste('wtf', ii, sep = '_')
      
      # re-save series index
      series.idx <- c(start.idx:end.idx)
      
      # if only one observation below 200m, discard
      if(sum(series$depth[series.idx] > 200, na.rm = T) <= 1) next
    
      # if ≥50% of observations within a bout are NAs then discard
      if(sum(is.na(series$depth[series.idx]))/length(series.idx) >= 0.5) next
      
      # if any NA runs longer than 15 mins (half median duration), discard bout
      # identify NA runs within the dive
      na.runs <- rle(is.na(series$depth[series.idx]))
      if(any(na.runs$values & na.runs$lengths > 6)) next
      
      # if any missing data periods longer than 15 mins. discard
      if(max(diff(series$DateTime_local[series.idx], units = 'mins')) > 15) next
      
      # if the bout only has three observations and contains any NAs, discard
      if(length(series.idx) <= 3 & anyNA(series$depth[series.idx])) next
      
      # if a bout has already been assigned, move to next (prevent later records from overriding earlier ones)
      # if(sum(!is.na(series$mesoBout[series.idx])) >= 1) next
      
      # if the bout passed all the tests, assign it to series
      series$mesoBout[series.idx] <- ii
    }
  }
    
    # convert date to a date object
    series$date <- as.Date(series$DateTime_local)
    
    # re-store the timeseries with mesobouts added
    series_list[[x]] <- series 
    
    # update the user
    print(
      glue::glue('Done with tag {x}')
    )
    
    # cleanup 
    rm(x, ii, series)
  }


combo_series <- purrr::list_rbind(series_list)


# mesoBout key ------------------------------------------------------------

# pair each mesobout with the date when most observations occur
mesoBout_stamp <-
  combo_series %>%
  filter(!is.na(mesoBout)) %>% 
  # make sure we're counting based on the right date
  mutate(date = lubridate::date(DateTime_local)) %>% 
  group_by(mesoBout, instrument_name) %>%
  count(date) %>% 
  filter(n == max(n)) %>%  # this should leave us with 1 data / bout (n = 3730)
  # there are actually 3734 obs. because 4 mesobouts evenly span two dates
  ungroup() %>% 
  separate(col = instrument_name, into = c('spp', 'year', 'ptt'), sep = '_', remove = F) %>% 
  transmute(mesoBout,
            instrument_name,
            species = if_else(spp == 159924, 'I.oxyrinchus', 'P.glauca'),
            ptt,
            # use my ptt code for 20P1624
            kode = if_else(ptt == '20P1624',
                           paste(as.character(date), '206771', sep = '_'),
                           paste(as.character(date), ptt, sep = '_'))) %>% 
  # add in the cluster codes
  left_join(
    cluster_stamp, # original has 786 unique kodes
    by = c('kode', 'species')) # 1129 NAs, non-NA = 643 unique kodes (18% loss)
# the missing codes are simply not represented in our list of mesoBouts, suggesting that no bouts occurred on those days
  
# bout duration -----------------------------------------------------------

bout_duration <- tibble()

# calculate the duration of each bout for each shark
for(x in names(series_list)) {
  
  # pull the series
  series <- series_list[[x]]
  
  series <-
    series %>%
    # select only obs. within a bout
    filter(!is.na(mesoBout)) %>% 
    # for each mesoBout, calculate the time difference from start to end
    group_by(instrument_name, mesoBout) %>% 
    summarize(
      start = min(DateTime),
      end = max(DateTime),
      max.z = max(depth, na.rm = T),
      .groups = 'drop') %>% 
    # calculate the difftime
    mutate(
      duration.min = lubridate::interval(start = start, end = end) %>% lubridate::time_length(unit = 'minutes'),
      species = if_else(str_detect(string = instrument_name, pattern = '160424'), 'P.glauca', 'I.oxyrinchus'),
      ptt = x)
  
  # store the answer
  bout_duration <- 
    rbind(bout_duration, series)
  
  # update user
  print(
    glue::glue('Done with {x}'))
  
  # cleanup
  rm(x, series)
  
}

# add in clusers
bout_duration <- 
  bout_duration %>% # n = 3730
  left_join(
    mesoBout_stamp,
    by = c('species', 'ptt', 'instrument_name', 'mesoBout')) %>%  # final n = 3734
  mutate(cluster = factor(cluster, levels = c('EPI 1', 'EPI 2', 'DVM 1', 'DVM 2', 'DVM 3 ')))


# visualize distribution
bout_duration %>% 
  filter(!is.na(cluster)) %>%
  ggplot() +
  geom_histogram(aes(x = duration.min), color = 'black', fill = 'grey', binwidth = 5) +
  annotate(geom = 'text', label = 'binwidth = 5 mins.', fontface = 'italic', color = 'grey', x = 750, y = 300, size = 5) +
  coord_cartesian() +
  theme_classic()

# distribution by cluster
duration.freqpoly.cluster <- 
  bout_duration %>% 
  filter(!is.na(cluster)) %>%
  ggplot() + 
  geom_freqpoly(aes(x = log10(duration.min), color = cluster, group = ptt)) +
  geom_vline(data = bout_duration %>% group_by(cluster, species) %>% summarize(median = median(duration.min), .groups = 'drop'), aes(xintercept = log10(median), linetype = species), alpha = 0.5) +
  facet_grid(cluster~.) +
  coord_cartesian() +
  scale_color_manual(values = c("yellow",
                                "yellow3",
                                "#488E9EFF",
                                "#404C8BFF",
                                "#281A2CFF"))
  theme_classic()

# ggsave(
#   duration.freqpoly,
#   file = 'products/figures/revisions/exploration/duration_freqpoly.png',
#   height = 6,
#   width = 8,
#   dpi = 300)

# distribution by species
duration.freqpoly.species <- 
  bout_duration %>% 
  ggplot() + 
  geom_freqpoly(aes(x = log10(duration.min), color = species, group = ptt)) +
  geom_vline(data = bout_duration %>% group_by(species, ptt) %>% summarize(median = median(duration.min), .groups = 'drop'), aes(xintercept = log10(median), color = species), alpha = 0.5) +
  facet_grid(species~.) +
  coord_cartesian() +
  guides(color = 'none') +
  theme_classic()

ggsave(
  duration.freqpoly.species,
  file = 'products/figures/revisions/exploration/duration_freqpoly_species.png',
  height = 6,
  width = 8,
  dpi = 300)


