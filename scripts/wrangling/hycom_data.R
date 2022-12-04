## GET ENVIRONMENTAL DATA 

# setup -------------------------------------------------------------------
source('scripts/wrangling/high_res_tad.R')

combo_series <- 
  read_csv('data/clean/Series/combo_series.csv')

library(tidyverse)
library(readr)
library(raster)
library(HMMoce)
library(lunar)
library(ncdf4)
devtools::load_all('analyzePSAT')


# GEBCO -------------------------------------------------------------------

# General Bathymetric Chart of the Ocean
## Data from : https://www.gebco.net/data_and_products/gridded_bathymetry_data/

# HYCOM -------------------------------------------------------------------

# create a stamp with coordinates for each tracking day
position.stmp <- combo_track %>% dplyr::select(latitude, longitude, kode)

# remove duplicates in days w. more than one location
position.stmp <- position.stmp[which(!duplicated(position.stmp$kode)),] 

# apply position data to high_res data frame
high_res <- left_join(high_res, position.stmp, by = c("kode"))
# extract environmental files from HYCOM for each unique date
## create an index of new columns where you want the environmental data to go. I happen to know this outputs 15 cols so I cheated...
col_idx <- c((ncol(high_res) + 1):(ncol(high_res) + 16))

t1 <- Sys.time()
for (tt in 1:nrow(high_res)){
  data <- unlist(facet_hycom(xpos = high_res$longitude[tt],
                             ypos = high_res$latitude[tt],
                             tpos = as.Date(high_res$Date[tt]), ## needs to be of class 'Date'
                             xlen = 0.25, ## these "errors" need to be approx this size to allow the calculations to run successfully
                             ylen = 0.25, 
                             varName = c('water_temp', 'water_u', 'water_v','surf_el','salinity')))
  data <- as.data.frame(t(data))
  high_res[tt,col_idx] <- data
  rm(data)
}

t2 <- Sys.time()
t2 - t1


# Chlorophyll -------------------------------------------------------------
source('analyzePSAT/R/vgpm.load.r')
## creating rasters ----
fileList <- list.files("data/raw/chloro/", full.names = TRUE)

tVec <- c(seq.Date(as.Date('2015-10-01'), as.Date('2016-04-01'), by='month'),
          seq.Date(as.Date('2016-09-01'), as.Date('2017-04-01'), by='month'),
          seq.Date(as.Date('2017-11-01'), as.Date('2018-01-01'), by='month'),
          seq.Date(as.Date('2021-10-01'), as.Date('2021-12-01'), by='month'))

for (i in 1:(length(tVec))){
  t1 <- Sys.time()
  doy <- lubridate::yday(tVec[i])
  if (nchar(doy) == 1) doy <- paste('00', doy, sep='')
  if (nchar(doy) == 2) doy <- paste('0', doy, sep='')
  yrday <- paste(lubridate::year(tVec[i]), doy, sep='')
  # load vgpm and rasterize
  vgpm <- vgpm.load(file=paste('data/raw/chloro/vgpm.', yrday, '.all.xyz', sep=''), w.lon = -180, e.lon = 180,
                    n.lat = 90, s.lat = -90, raster = TRUE)
  vgpm <- vgpm / 1000 # convert from mg to g C/ m^2 / d^1
  #t2 <- Sys.time()
  #print(paste(t2 - t1)) # print how long it took
  crs_geo <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  crs(vgpm) <- crs_geo
  t3 <- Sys.time()
  vgpm_proj <- projectRaster(vgpm, crs=crs_geo)
  
  #plot(vgpm_proj)
  
  e <- extent(-80, -35, 0, 50)
  e <- as(e, 'SpatialPolygons')
  crs(e) <- crs_geo
  
  vgpm_new <- mask(vgpm_proj, e)
  vgpm_new <- projectRaster(vgpm_new, crs=crs_geo)
  
  writeRaster(
    vgpm_new,
    file=paste('data/clean/Chlorophyll/vgpm.', yrday, '.all.grd', sep=''),
    overwrite=T)
  rm(vgpm); rm(vgpm_new); rm(vgpm_proj); gc()
  t4 <- Sys.time()
  print(paste(i, 'took', round(t4 - t1), 'secs'))
}

rm(fileList)

## extracting values ----
fileList <- list.files("data/clean/Chlorophyll/", full.names = TRUE)
fileList <- fileList[grep('.gri', fileList)]

data <- tibble()

for(i in 1:(length(tVec))) {
  doy <- lubridate::yday(tVec[i])
  if (nchar(doy) == 1) doy <- paste('00', doy, sep='')
  if (nchar(doy) == 2) doy <- paste('0', doy, sep='')
  yrday <- paste(lubridate::year(tVec[i]), doy, sep='')
  file <- fileList[grep(yrday, fileList)]
  c_brick <- raster::raster(file)
  sub1 <- high_res[high_res$Date >= tVec[i] & high_res$Date < tVec[i+1],]
  extraction_points <- data.frame(
    x = sub1$longitude,
    y = sub1$latitude
  )
  for(tt in 1:nrow(extraction_points)) {
    point_series <- raster::extract(c_brick, extraction_points[tt,], method='simple')
    sub1[tt,42] <- point_series
  }
  
  data <- rbind(data, sub1)
  print(i)
}

high_res <- 
  data %>% 
  filter(!is.na(ptt)) %>% 
  as_tibble() %>% 
  mutate(Date = lubridate::ymd(Date),
         chloro = V42) %>%
  select(-c('V42', water_u, water_v, tpos))

rm(data,
   fileList,
   sub1,
   doy,
   file,
   point_series,
   tt,
   tVec,
   yrday,
   i,
   c_brick,
   extraction_points)

# Lunar illumination ------------------------------------------------------

high_res$lunar <- 
  lunar::lunar.illumination(
    as.Date(high_res$Date))


# HYCOM for all tracking days ---------------------------------------------

positions <- combo_track %>% dplyr::select(Date, ptt, latitude, longitude, kode)

# remove duplicates in days w. more than one location
positions <- 
  positions[which(!duplicated(positions$kode)),] 

# match w. series data
positions <- left_join(combo_series, positions, by = c("kode"))

# extract environmental files from HYCOM for each unique date
## create an index of new columns where you want the environmental data to go. I happen to know this outputs 15 cols so I cheated...
col_idx <- c((ncol(positions) + 1):(ncol(positions) + 16))

t1 <- Sys.time()
for (tt in 1:nrow(positions)){
  data <- unlist(facet_hycom(xpos = high_res$longitude[tt],
                             ypos = high_res$latitude[tt],
                             tpos = as.Date(high_res$Date[tt]), ## needs to be of class 'Date'
                             xlen = 0.25, ## these "errors" need to be approx this size to allow the calculations to run successfully
                             ylen = 0.25, 
                             varName = c('water_temp', 'water_u', 'water_v','surf_el','salinity')))
  data <- as.data.frame(t(data))
  high_res[tt,col_idx] <- data
  rm(data)
}

t2 <- Sys.time()
t2 - t1
