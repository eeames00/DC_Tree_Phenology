## Load in Libraries

library(data.table)
library(tidyverse)

### Load in data folder

in_dir <- "C:/Users/valea/OneDrive/Desktop/extracts_20260121"

files <- list.files(
  path = in_dir,
  pattern = "^crown_ndvi_\\d{8}_\\d{6}\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

## Check how many files, paths, and a quick look at one ex file
length(files)      
files[1:5]         

dt_one <- fread(files[1], na.strings = c("NA", "", "NaN"))
dt_one[1:10]        
names(dt_one)       

## Combine all files in folder

all_dt <- rbindlist(
  lapply(files, \(f) fread(f, na.strings = c("NA", "", "NaN"))),
  use.names = TRUE,
  fill = TRUE
)
## Check # of rows and column names
nrow(all_dt)        
names(all_dt)       


####### For Sourcing Data from Each File if needed #########

all_dt <- rbindlist(
  lapply(files, \(f) {
    dt <- fread(f, na.strings = c("NA", "", "NaN"))
    dt[, source_file := basename(f)]
    dt
  }),
  use.names = TRUE,
  fill = TRUE
)

all_dt[, .N, by = source_file][1:10]   

##############################################################



## Save combined file (two forms)
fwrite(all_dt, file.path(in_dir, "compiled_crown_ndvi.csv"))
saveRDS(all_dt, file.path(in_dir, "compiled_crown_ndvi.rds"))

## Second Set of Data - cleaned/combined time column 


all_dt[, im_datetime := as.POSIXct(
  paste0(im_date, sprintf("%06d", im_time)),
  format = "%Y%m%d%H%M%S",
  tz = "UTC"
)]


all_dt[1:5, .(im_date, im_time, im_datetime, im_date_time)]



## Remove "0" -> these are missing Tree
all_dt <- all_dt[UID != 0]
sum(all_dt$UID == 0)
all_dt[UID == 0, .N]


all_dt[, SpCode := toupper(trimws(SpCode))]


## Save Further Cleaned Files 

fwrite(all_dt, file.path(in_dir, "compiled_crown_ndvi_cleaned.csv"))
saveRDS(all_dt, file.path(in_dir, "compiled_crown_ndvi_cleaned.rds"))



## Check cleaned file with small glimpse

all_dt[, .N, by = im_date][order(im_date)][1:10]

## Check NDVI Range 

all_dt[, .(ndvi_min = min(ndvi, na.rm=TRUE),
           ndvi_max = max(ndvi, na.rm=TRUE),
           ndvi_mean = mean(ndvi, na.rm=TRUE))]

## More detailed data check for Cleaned Dataset

nrow(all_dt)
uniqueN(all_dt$UID)
uniqueN(all_dt$SpCode)
range(all_dt$im_date)

all_dt[, .(
  ndvi_na = sum(is.na(ndvi)),
  r2red_na = sum(is.na(im_r2_red)),
  r2nir_na = sum(is.na(im_r2_nir))
)]

## Sort/key for fast filtering by tree + time 

setkey(all_dt, UID, im_datetime)
all_dt[1:5]


## Daily Per Tree Dataset 
#### one value per UID per datetime (should already be the way the data is but just covers it in case)

tree_ts <- all_dt[, .(
  ndvi = mean(ndvi, na.rm = TRUE)
), by = .(UID, SpCode, im_datetime, im_date, im_doy)]

nrow(tree_ts)
tree_ts[1:5]


## Another quality check 

one_uid <- all_dt$UID[1]

plot(
  all_dt[UID == one_uid]$im_datetime,
  all_dt[UID == one_uid]$ndvi,
  xlab = "Date",
  ylab = "NDVI",
  main = paste("UID", one_uid),
  pch = 16
)

## Save both files (filtered option and also daily tome series)
saveRDS(all_dt,  file.path(in_dir, "compiled_crown_ndvi_filtered.rds"))
saveRDS(tree_ts, file.path(in_dir, "tree_timeseries_daily.rds"))








