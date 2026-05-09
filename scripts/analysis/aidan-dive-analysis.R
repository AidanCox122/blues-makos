library(dplyr)
devtools::load_all('../tags2etuff')
# devtools::load_all('/Users/aidansmacpro/Documents/Ecology/Marine/tags2etuff')


## get meta & filter
meta <- read.table('../nip_drake/RawData/all_tag_meta.csv', sep=',', header=T, blank.lines.skip = F, skip=0, stringsAsFactors = F)
ids <- c("159924_2015_141257", "159924_2017_163096", "159924_2017_163098",
         "159924_2021_20P1624", # which is PTT 206771", ## makos
         "160424_2015_141256", "160424_2015_141254", "160424_2015_141259", 
         "160424_2016_106754", "160424_2016_133016", "160424_2016_133017",
         "160424_2016_133018", "160424_2016_133021", "160424_2016_141247",
         "160424_2016_141255", "160424_2016_141258", "160424_2016_154096",
         "160424_2016_163097") ## blues

which(!(ids %in% meta$instrument_name))
meta_sub <- meta %>% filter(instrument_name %in% ids)

## cam's data directory
existing_dir <- '~/Desktop/data_org/'
## aidan directory for outputs
aidan_dir <- '~/Google Drive/Shared drives/MPG_WHOI/data/etuff_aidan/series_outputs/'

for (i in 1:nrow(meta_sub)){

  if (is.na(meta_sub$tag_recovered[i])){
    fList <- list.files(data_dir)
    meta_sub$tag_recovered[i] <- ifelse(length(grep('-Archive.csv', fList)) > 0, 'yes', 'no')
  }
  
  ## prioritize eTUFF_hdr.txt unless eTUFF.txt is newer
  data_dir <- paste(existing_dir, meta_sub$instrument_name[i], '/', sep='')
  fname1 <- paste(data_dir, meta_sub$instrument_name[i], '_eTUFF.txt', sep='')
  fname2 <- paste(data_dir, meta_sub$instrument_name[i], '_eTUFF_hdr.txt', sep='')
  if (file.exists(fname1) & file.exists(fname2)){
    etuff_file <- ifelse(file.info(fname1)$mtime > file.info(fname2)$mtime,
                         etuff_file <- fname1,
                         etuff_file <- fname2)
  } else{
    etuff_file <- fname1
  }
  
  ## read some initial data
  etuff <- try(read_archival(etuff_file, metaTypes = metaTypes), TRUE)
  if (class(etuff)[1] == 'try-error') etuff <- read_etuff(etuff_file)
  track <- try(get_track(etuff), TRUE)
  if (class(track) == 'try-error') stop('Trouble getting track. Does this eTUFF have location information?')
  #series <- try(get_3d(etuff), TRUE)
  
  #=================================
  ## get PDT
  #=================================
  pdt_try <- try(pdt <- get_pdt(etuff), TRUE)
  if (class(pdt_try) == 'try-error'){
    pdt <- NA
  } else{
    pdt <- get_pdt(etuff)
  }
  
  #=================================
  ## interp PDT
  #=================================
  if (is.na(meta_sub$tag_recovered[i]) & length(grep('Archive', list.files(data_dir))) == 0) meta_sub$tag_recovered[i] <- 'no'
  if (file.exists(paste('./results/', meta_sub$instrument_name[i], '_pdt_interp.RDS', sep='')) & !is.na(pdt)[1]){
    pdt_interp <- readRDS(paste('./results/', meta_sub$instrument_name[i], '_pdt_interp.RDS', sep=''))
  } else if (meta_sub$tag_recovered[i] == 'yes'){
    pdt_interp <- NULL
  } else if (!file.exists(paste('./results/', meta_sub$instrument_name[i], '_pdt_interp.RDS', sep='')) & !is.na(pdt)[1]){
    pdt_interp <- interp_pdt(pdt)
    pdt_interp$id <- meta_sub$instrument_name[i]
    saveRDS(pdt_interp, file=paste('./results/', meta_sub$instrument_name[i], '_pdt_interp.RDS', sep=''))#, sep=',', col.names=F, row.names=F, append=T)
  } else{
    pdt_interp <- NULL
  }
  
  
  #=================================
  ## get series/srss/dives
  #=================================
  
  series <- try(get_series(etuff), TRUE)
  
  if (class(series) != 'try-error' & !is.na(series)[1]){
    if (!any(names(series) %in% c('depth','pressure'))){
      series <- NA
    } else{
      
      if (!is.null(pdt_interp)) series <- add_series_temp(series, pdt, pdt_interp)
      
      ## format series and srss
      series$dn <- add_daynight(series, etuff)
      series$instrument_name <- meta_sub$instrument_name[i]
      srss <- get_srss(etuff)
      srss$instrument_name <- meta_sub$instrument_name[i]
      write.table(srss, file=paste0(aidan_dir, '/', meta_sub$instrument_name[i], '_srss.csv'), sep=',')#, col.names=F, row.names=F, append=T)
      
      ## get dives
      series <- series[order(series$DateTime_local),]
      std <- diveMove::createTDR(series$DateTime_local, series$depth, #concurrentData = series[,-c(12,4)],
                                 file='xxx') # file argument required but doesnt seem to matter in this case
      
      try.stdc <- try(stdc <- calibrateDepth(std, dry.thr=0, wet.cond=rep(TRUE, length.out=nrow(series)),
                                             dive.thr = 50, zoc.method='offset', offset=0, dive.model='smooth.spline'), silent=TRUE)
      
      if (class(try.stdc) == 'try-error'){
        stdc <- diveMove::calibrateDepth(std, dry.thr=0, wet.cond=rep(TRUE, length.out=nrow(series)),
                                         dive.thr = 50, zoc.method='offset', offset=10, dive.model='smooth.spline')
      } else{
        stdc <- try.stdc
      }
      rm(try.stdc); gc()
      #plotTDR(std, interact=FALSE)
      #abline(h=200,col='red')
      #plotTDR(stdc, diveNo=30, interact=F)
      
      # get dive stats and figure out which are useful
      #if (max(stdc@dive.activity[,1] > 0)){
      stdsumm <- try(diveMove::diveStats(stdc), silent=TRUE)
      
      if (class(stdsumm) != 'try-error'){
        
        # critera req'd such as max dive time is 24 hours (or less?) and
        # min dive time is 30 mins?
        diveIdx <- which(stdsumm$maxdep >= 200)# & stdsumm$divetim > 1800 & stdsumm$divetim <= (24*3600))
        
        series$rownum <- 1:nrow(series)
        series$diveNo <- NA
        #tres <- as.numeric(std@time[2] - std@time[1], 'secs')
        ## AN INDEX ONLY FOR DEEP DIVES
        for (ii in diveIdx){
          start.idx <- which.min(as.numeric((series$DateTime_local - stdsumm$begdesc[ii])) ^ 2)
          if (!is.finite(start.idx)) start.idx <- 1
          end.idx <- which.min(as.numeric((series$DateTime_local - (stdsumm$begdesc[ii] + stdsumm$divetim[ii]))) ^ 2)
          if (!is.finite(end.idx)) end.idx <- nrow(series)
          series.idx <- c(start.idx:end.idx)
          
          if (is.na(series$depth[end.idx]) | series$depth[end.idx] > 50) end.idx <- min(which(series$depth <= 50 & series$rownum > series.idx[length(series.idx)]))
          if (!is.finite(end.idx)) next
          
          # get time resolution and set row numbers for indexing
          
          # search for data before and after current dive index that is above min depth threshold (50m)
          # if those data arent within one hour of the current index, we revert back to existing index
          # and discard the dive
          start.idx_new <- max(which(series$depth <= 50 & series$rownum <= series.idx[1]))
          if (!is.finite(start.idx_new)) next
          if (as.numeric(difftime(series$DateTime[start.idx], series$DateTime[start.idx_new], units='mins')) <= 60) start.idx <- start.idx_new
          
          series.idx <- c(start.idx:end.idx)
          series$diveNo[series.idx] <- ii
        }
        
        series$date <- as.Date(series$DateTime_local)
        
        ## how many dives below 200 per day? conservative estimate as it reqs the series data measures fish leaving from and returning to mixed layer (shallower than 50 m)
        stdsumm$date <- as.Date(stdsumm$begdesc)
        dives <- stdsumm %>% group_by(date) %>% summarise(n=n(), dives200 = length(which(maxdep > 200)),
                                                          mean_maxz = round(mean(maxdep[which(maxdep > 200)]), 1),
                                                          median_maxz = round(median(maxdep[which(maxdep > 200)]), 1),
                                                          mean_bottdep = round(mean(bottdep.mean[which(maxdep > 200 & !is.na(botttim))], na.rm=T), 1),
                                                          mean_bottdep.median = round(mean(bottdep.median[which(maxdep > 200 & !is.na(botttim))], na.rm=T), 1),
                                                          mean_bottdep.sd = round(mean(bottdep.sd[which(maxdep > 200 & !is.na(botttim))], na.rm=T), 1)) %>% as.data.frame()
        
        ## here the mean of max depths can be lower than the mean bottom depth because all dives below 200 are included in the max calculation
        ## while only those dives that actually have a bottom time (i.e. dives with multiple depth measurements below 200) are included in the latter
        stdsumm$instrument_name <- meta_sub$instrument_name[i]
        
      } else{
        series$rownum <- 1:nrow(series)
        series$diveNo <- NA
        series$date <- as.Date(series$DateTime_local)
        
      }
      
      #-------------------------
      ## ADD TRACK TO SERIES
      ## get track temporal interval
      names(track)[1] <- 'track_DateTime'
      track_interval <- lubridate::interval(track$track_DateTime[1:(nrow(track)-1)], track$track_DateTime[2:(nrow(track))])
      
      nms <- names(series)
      
      ## identify idx
      idx <- lapply(series$DateTime, FUN=function(x){
        i <- which(x %within% track_interval)
        #print(i)
        if (length(i) == 0) i <- NA
        if (length(i) == 2) i <- i[1]
        i
      }) %>% unlist()
      if (length(idx) != nrow(series)) stop('Series and matching track interval are not of same length.')
      
      ## then combine
      series$idx <- idx
      track$idx <- 1:nrow(track)
      series <- merge(series, track, by='idx', all.x = TRUE)
      series$track_DateTime <- as.POSIXct(series$track_DateTime, origin='1970-01-01', tz='UTC')
      series <- series[,c(nms, 'track_DateTime','latitude','longitude')]
      
      
      #-------------------------
      ## ADD TRACK TO DIVES
      if (class(stdsumm) != 'try-error'){
        ## get track temporal interval
        nms <- names(dives)
        dives$DateTime <- as.POSIXct(dives$date, tz='UTC')
        
        ## identify idx
        idx <- lapply(dives$DateTime, FUN=function(x){
          i <- which(x %within% track_interval)
          #print(i)
          if (length(i) == 0) i <- NA
          if (length(i) == 2) i <- i[1]
          i
        }) %>% unlist()
        if (length(idx) != nrow(dives)) stop('Dives and matching track interval are not of same length.')
        
        ## then combine
        dives$idx <- idx
        track$idx <- 1:nrow(track)
        dives <- merge(dives, track, by='idx', all.x = TRUE)
        dives$track_DateTime <- as.POSIXct(dives$track_DateTime, origin='1970-01-01', tz='UTC')
        dives <- dives[,c(nms, 'track_DateTime','latitude','longitude')]
        dives$instrument_name <-  meta_sub$instrument_name[i]
        dives$platform <-  meta_sub$platform[i]
        
      }
      
      #-------------------------
      ## SUMMARISE TIME > 200 M / DAY
      series_meso200 <- series %>% group_by(date) %>% summarise(n=n(), longitude=mean(longitude), latitude=mean(latitude), meso200 = length(which(depth >= 200)) / n() * 100) 
      series_meso200$type <- 'series'
      
      ## more complete summary of series
      series_summ <- series %>% group_by(date) %>% 
        summarise(n=n(), n_valid = length(which(!is.na(depth))), maxz = max(depth, na.rm=T),
                  longitude=mean(longitude), latitude=mean(latitude), 
                  meso200 = length(which(depth >= 200)) / n() * 100) %>% 
        filter(n_valid != 0)
      series_summ$instrument_name <- meta_sub$instrument_name[i]
      
      write.table(series_summ, file=paste0(aidan_dir, '/', meta_sub$instrument_name[i], '_series_summ.csv'), sep=',', col.names=T, row.names=F, append=F)
      if (class(stdsumm) != 'try-error') write.table(stdsumm, file=paste0(aidan_dir, '/', meta_sub$instrument_name[i], '_stdsumm.csv'), sep=',', col.names=T, row.names=F, append=F)
      if (class(stdsumm) != 'try-error') write.table(dives, file=paste0(aidan_dir, '/', meta_sub$instrument_name[i], '_dives200.csv'), sep=',', col.names=T, row.names=F, append=F)
      
    }
  }
  
  
  #=================================
  ## get mmd
  #=================================
  
  mmd <- try(get_mmd(etuff), TRUE)
  
  if (class(mmd) != 'try-error' & !is.na(mmd)[1]){
    mmd$instrument_name <- meta_sub$instrument_name[i]
    
    ## ADD TRACK TO mmd
    names(track)[1] <- 'track_DateTime'
    track_interval <- lubridate::interval(track$track_DateTime[1:(nrow(track)-1)], track$track_DateTime[2:(nrow(track))])
    
    nms <- names(mmd)
    
    ## identify idx
    idx <- lapply(mmd$DateTime, FUN=function(x){
      i <- which(x %within% track_interval)
      #print(i)
      if (length(i) == 0) i <- NA
      if (length(i) == 2) i <- i[1]
      i
    }) %>% unlist()
    if (length(idx) != nrow(mmd)) stop('Series and matching track interval are not of same length.')
    
    ## then combine
    mmd$idx <- idx
    track$idx <- 1:nrow(track)
    mmd <- merge(mmd, track, by='idx', all.x = TRUE)
    mmd$track_DateTime <- as.POSIXct(mmd$track_DateTime, origin='1970-01-01', tz='UTC')
    mmd <- mmd[,c(nms, 'track_DateTime','latitude','longitude')]
    mmd$date <- as.Date(mmd$DateTime, tz='UTC')
    mmd <- mmd %>% dplyr::select(date, depthMax, longitude, latitude, instrument_name)
    names(mmd)[2] <- 'maxz'
    mmd$type <- 'mmd'
  }
  
  if (class(series)[1] != 'try-error' & length(series) > 1){
    if (class(series_summ)[1] != 'try-error' & length(series_summ) > 1){
      series_mmd <- series_summ %>% dplyr::select(date, maxz, longitude, latitude, instrument_name)
      series_mmd$type <- 'series'
    }
  } else{
    series_mmd <- NA
  }
  
  ## end up with series and/or mmd combined with track data
  if (class(series_mmd)[1] != 'try-error' & !is.na(series_mmd)[1] & class(mmd)[1] != 'try-error'){
    combine <- rbind(series_mmd, mmd)
  } else if ((class(series_mmd)[1] == 'try-error' | is.na(series_mmd)[1]) & class(mmd)[1] != 'try-error'){
    combine <- mmd
  } else if ((class(series_mmd)[1] != 'try-error' & !is.na(series_mmd)[1]) & class(mmd)[1] == 'try-error'){
    combine <- series_mmd
  } else{
    combine <- NA
  }
  
  if (!is.na(combine)[1]){
    combine$instrument_name <- meta_sub$instrument_name[i]
    combine$platform <- meta_sub$platform[i]
    combine <- combine %>% filter(!is.na(longitude))
    
    write.table(combine, file=paste0(aidan_dir, '/', meta_sub$instrument_name[i], '_mmd.csv'), sep=',', col.names = TRUE, row.names = F)

  }
  
  #=================================
  ## write series
  #=================================
  
  if (class(series) != 'try-error' & !is.na(series)[1]){
    series$instrument_name <- meta_sub$instrument_name[i]
    series$platform <-  meta_sub$platform[i]
    write.table(series, file=paste0(aidan_dir, '/', meta_sub$instrument_name[i], '_series.csv'), sep=',', col.names=T, row.names=F, append=F)
  }
  
  rm(etuff); rm(series); rm(mmd); rm(tad); rm(tat); rm(pdt); rm(stdsumm); rm(dives); 
  rm(pdt_interp); rm(series_summ); rm(series_mmd); rm(series_tad); rm(srss); rm(std); 
  rm(stdc); rm(tr); rm(combine); rm(series_meso200); rm(track); rm(series.idx); rm(series.idx_new);
  rm(diveIdx); rm(pdt_try); rm(idx); rm(nms); rm(track_interval); rm(fname1); rm(fname2); gc()
  print(i)
  
}

