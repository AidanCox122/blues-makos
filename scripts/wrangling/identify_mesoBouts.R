### THIS SCRIPT CALCULATES STATS. FOR INDIV. DIVES FROM HIGH-RES. TIMESERIES DATA FOR ALL TAGGED SHARKS
# setup -------------------------------------------------------------------

library(tidyverse)
library(glmmTMB)
library(DHARMa)
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

latlon_stamp <-
  read_csv('data/clean/latlon_stamp.csv')

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
      
      # if the gap between the bout start and the previous record is more than 2.5 mins, discard
      if(difftime(series$DateTime_local[start.idx], series$DateTime_local[start.idx-1], units = 'mins') > 2.5) next
      
      # if the gap between the bout end and the next record is more than 2.5 mins, discard
      if(difftime(series$DateTime_local[end.idx+1], series$DateTime_local[end.idx], units = 'mins') > 2.5) next
      
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

# combine series_list into a single dataframe
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
  # add in cluster ID
  left_join(
    mesoBout_stamp,
    by = c('species', 'ptt', 'instrument_name', 'mesoBout')) %>%  # final n = 3734
  filter(!is.na(cluster)) %>% # n = 2605
  # add in lat and lon
  left_join(
    latlon_stamp,
    by = c('kode', 'mesoBout')) %>% # several mesoBouts are associated with multiple location estimates
  # take the average position
  group_by(instrument_name, mesoBout, start, end, max.z, duration.min, species, ptt, kode, cluster, cluster_n) %>% 
  summarize(lat = mean(latitude, na.rm = T), lon = mean(longitude, na.rm = T), .groups = 'drop') %>% 
  separate(kode, into = c('Date', 'ptt'), sep = '_', remove = F) %>% 
  mutate(
    cluster = factor(cluster, levels = c('EPI 1', 'EPI 2', 'DVM 1', 'DVM 2', 'DVM 3'), ordered = T),
    Date = lubridate::ymd(Date),
    month = lubridate::month(Date))


# bout frequency ----------------------------------------------------------

# create a repository
bout_frequency <- 
  tibble()

# extract bout frequencies
for(x in names(series_list)) {
  
  # pull the series
  series <- series_list[[x]]
  
  # detect the interval between records
  interval <- diff(series$DateTime_local, units = 'min') %>% median() %>% as.numeric()
  
  # calculate number of bouts and number of records for each shark-day
  bout_freq_x <- 
    series %>% 
    # standardize date for consistency across code chunks
    mutate(date = lubridate::date(DateTime_local)) %>% 
    # count the number of bouts and number of valid series records in each shark-day
    group_by(instrument_name,
             date) %>% 
    summarize(
      n_bouts = sum(!is.na(unique(mesoBout))),
      n_records = sum(!is.na(depth)),
      n_hours = (n_records * interval)/60,
      .groups = 'drop') %>% 
    # format result
    separate(col = instrument_name, into = c('spp', 'year', 'ptt'), sep = '_', remove = F) %>% 
    transmute(
      instrument_name,
      species = if_else(str_detect(string = instrument_name, pattern = '160424'), 'P.glauca', 'I.oxyrinchus'),
      ptt,
      date,
      kode = if_else(ptt == '20P1624',
                     paste(as.character(date), '206771', sep = '_'),
                     paste(as.character(date), ptt, sep = '_')),
      n_bouts,
      n_records,
      n_hours,
      bout_freq = n_bouts / n_hours)
  
  # save the results
  bout_frequency <- 
    rbind(bout_frequency, bout_freq_x)
  
  # update user
  print(
    glue::glue('Done with {x}')
  )
  
  # cleanup
  rm(x, series, bout_freq_x)
  
}

# add in clusters
bout_frequency <- 
  bout_frequency %>% 
  left_join(
    mesoBout_stamp %>% dplyr::select(-c(mesoBout)) %>% distinct(),
    by = c('instrument_name', 'species', 'ptt', 'kode')) %>% 
  filter(!is.na(cluster)) %>% 
  # add in lat and lon
  left_join(
    latlon_stamp %>% group_by(kode) %>% summarize(lat = mean(latitude, na.rm = T), lon = mean(longitude, na.rm = T)),
    by = c('kode')) %>% # several mesoBouts are associated with multiple location estimates
  mutate(
    month = lubridate::month(date),
    cluster = factor(cluster, levels = c('EPI 1', 'EPI 2', 'DVM 1', 'DVM 2', 'DVM 3'), ordered = T)) 
  

# cold exposure -----------------------------------------------------------

temp_duration <- tibble()

# calculate the duration of each bout for each shark
for(x in names(series_list)) {
  
  # pull the series
  series <- series_list[[x]]
  
  # interval <- diff(series$DateTime_local, units = 'min') %>% median() %>% as.numeric()
  
  temps <- 
    series %>% 
    transmute(
      instrument_name,
      mesoBout) %>% 
    distinct()

  
  for(t in c(14, 16, 18)) {
    # create a column name
    daname <- paste('deg',t,'min', sep = '.')
    
    # get degree minutes
    deg.min <-
      series %>%
      # select only obs. within a bout
      filter(!is.na(mesoBout) & !is.na(t)) %>% 
      # filter(temperature < t) %>% # this was messy because it fudged up the time interval
      mutate(delta.temp = t - temperature,
             # if temperature is above threshold, convert delta temp to 0 (this will zero it in our sum later)
             delta.temp = if_else(
               delta.temp < 0 | is.na(delta.temp),
               0,
               delta.temp)
             ) %>% 
      # for each mesoBout, calculate the time difference from start to end
      group_by(instrument_name, mesoBout) %>%
      mutate(interval = (lead(DateTime_local) - DateTime_local) %>% lubridate::time_length(unit = 'minutes'),
             temp.mins = delta.temp * interval) %>% 
      summarize(
        !!daname := sum(temp.mins, na.rm = T),
        .groups = 'drop') 
    
    # add columns to the repository
    temps <- 
      left_join(
      temps,
      deg.min,
      by = c('instrument_name', 'mesoBout'))
  }
  
  # add some handy columns
  temps <- 
    temps %>% 
    mutate(
      species = if_else(str_detect(string = instrument_name, pattern = '160424'), 'P.glauca', 'I.oxyrinchus'),
      ptt = x) 
    
  # store the answer
  temp_duration <- 
    rbind(temp_duration, temps)
  
  # update user
  print(
    glue::glue('Done with {x}'))
  
  # cleanup
  rm(x, series, deg.min, temps, daname)
  
}

# add in clusters

temp_duration <- 
  temp_duration %>% # n = 3747
  # add in cluster ID
  left_join(
    mesoBout_stamp,
    by = c('species', 'ptt', 'instrument_name', 'mesoBout')) %>%  # final n = 3751
  filter(!is.na(cluster)) %>% # n = 2605
  # add in lat and lon
  left_join(
    latlon_stamp,
    by = c('kode', 'mesoBout')) %>% # several mesoBouts are associated with multiple location estimates
  # take the average position
  group_by(instrument_name, mesoBout, deg.14.min, deg.16.min, deg.18.min, species, ptt, kode, cluster, cluster_n) %>% 
  summarize(lat = mean(latitude, na.rm = T), lon = mean(longitude, na.rm = T), .groups = 'drop') %>% 
  separate(kode, into = c('Date', 'ptt'), sep = '_', remove = F) %>% 
  mutate(
    cluster = factor(cluster, levels = c('EPI 1', 'EPI 2', 'DVM 1', 'DVM 2', 'DVM 3'), ordered = T),
    Date = lubridate::ymd(Date),
    month = lubridate::month(Date))

# pre-/post-dive surfacing ------------------------------------------------

surface_use <- 
  tibble()

baseline_surface_use <- 
  tibble()

for(x in names(series_list)) {
 
  # extract series info
  series <- 
    series_list[[x]] 
  
  # extract start time / end time
  starttime <- series$DateTime_local[1]
  endtime <- series$DateTime_local[nrow(series)]
  
  # extract bouts
  bouts <-
    bout_duration %>%
    # final all bouts from individual x
    filter(str_detect(instrument_name, x)) %>% 
    # consider only bouts one hour after series start and one hour before series end
    mutate(start.diff = (start - starttime) %>% lubridate::time_length(unit = 'minutes') %>% as.numeric(),
           end.diff = (end - endtime) %>% lubridate::time_length(unit = 'minutes') %>% as.numeric()) %>% 
    filter(start.diff >= 60,
           abs(end.diff) >= 60)
  
  # calculate pre-post bout stats for each 
  stats_x <- tibble()
  
  series$pre.post <- NA
    
  for(i in 1:nrow(bouts)) {
    # find observations 
    bout_i <- bouts[i,]
    
    # find observations within an hour of the bout start
    pre.stats <- 
      series %>% 
      mutate(
        difftime = 
          (bout_i$start - DateTime) %>%
          lubridate::time_length(unit = 'minutes') %>%
          as.numeric()) %>% 
      filter(difftime > 0 & difftime <= 60) %>% 
      mutate(prop.10 = depth <= 10,
             prop.50 = depth <= 50) %T>% 
      # save the list of before rows for baseline calc later
      assign('pre.rows', ., envir = .GlobalEnv) %>% 
      # now calculate summary stats.
      summarize(
        pre.prop.10 = mean(prop.10, na.rm = T),
        pre.prop.50 = mean(prop.50, na.rm = T)) 
    
    # find observations within an hour of the bout end
    post.stats <- 
      series %>% 
      mutate(
        difftime = 
          (bout_i$end - DateTime) %>%
          lubridate::time_length(unit = 'minutes') %>%
          as.numeric()) %>% 
      filter(difftime < 0 & difftime >= -60) %>% 
      mutate(prop.10 = depth <= 10,
             prop.50 = depth <= 50) %T>% 
      # save the list of post rows for baseline calc later
      assign('post.rows', ., envir = .GlobalEnv) %>% 
      # now calculate summary stats. 
      summarize(
        post.prop.10 = mean(prop.10, na.rm = T),
        post.prop.50 = mean(prop.50, na.rm = T)) 
    
    # assemble the pre and post stats w. metadata
    bout_i_final <- 
      bout_i %>% 
      dplyr::select(-c(duration.min, max.z)) %>% 
      cbind(pre.stats) %>% 
      cbind(post.stats)
    
    # store the data from each bout
    stats_x <- rbind(stats_x, bout_i_final)
    
    # get rownums for pre- and post-bout surface intervals
    pre.rows <- seq(min(pre.rows$rownum), max(pre.rows$rownum))
    post.rows <- seq(min(post.rows$rownum), max(post.rows$rownum))
    
    # mark what rows were used in this calculation for baseline calc. later
    series$pre.post[pre.rows] <- T
    series$pre.post[post.rows] <- T
  }
  
  # calculate series baseline
  baseline <- 
    series %>% 
    filter(is.na(mesoBout)) %>% 
    filter(is.na(pre.post)) %>% 
    mutate(
      date = lubridate::date(DateTime_local),
      prop.10 = depth <= 10,
      prop.50 = depth <= 50) %>% 
    # calculate baseline surface-use proportion
    group_by(species = platform, instrument_name, date, dn) %>% 
    summarize(
      base.prop.10 = mean(prop.10, na.rm = T),
      base.prop.50 = mean(prop.50, na.rm = T),
      .groups = 'drop') %>% 
    # add information needed for cluster join
    mutate(
      ptt = x,
      species = if_else(species == 'Prionace glauca', "P.glauca", "I.oxyrinchus"))
    
  # store baseline surface-use stats
  baseline_surface_use <- 
    rbind(baseline_surface_use, baseline)

  # store the pre/post-bout surface use stats for each shark
  surface_use <- 
    rbind(surface_use, stats_x)
  
  # update the user
  print(
    glue::glue('Done with {x}')
  )
  
  # cleanup
  rm(series, bouts, stats_x, bout_i, pre.stats, post.stats, bout_i_final, x, i, baseline, pre.rows, post.rows)
}

# add in cluster info to baseline stats
baseline_surface_use <- 
  baseline_surface_use %>% 
  mutate(
    # use my ptt code for 20P1624
    kode = if_else(ptt == '20P1624',
                   paste(as.character(date), '206771', sep = '_'),
                   paste(as.character(date), ptt, sep = '_'))) %>% 
  # add in cluster result for each date
  left_join(
    mesoBout_stamp %>% dplyr::select(-c(mesoBout)) %>% distinct(),
    by = c('species', 'ptt', 'instrument_name', 'kode')) %>% 
  mutate(cluster = factor(cluster, levels = c('EPI 1', 'EPI 2', 'DVM 1', 'DVM 2', 'DVM 3'), ordered = T))

# check the output
surface_use %>% summary()

