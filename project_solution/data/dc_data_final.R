###Full Working Dataset for DC###

library(tidyverse)
library(data.table)
library(readxl)

tree_timeseries_daily_polygon_subset <- readRDS("C:/Users/valea/OneDrive/Desktop/tree_timeseries_daily_polygon_subset.rds")

View(tree_timeseries_daily_polygon_subset)


###dc_trees 

dc_trees <- read_excel("dc_trees.xlsx")
View(dc_trees)

library(dplyr)

# NDVI time-series data
ndvi_dat <- tree_timeseries_daily_polygon_subset

# Clean ArcGIS tree data before join
arc_dc_trees_clean <- dc_trees %>%
  filter(!is.na(UID), UID != 0) %>%
  distinct(UID, .keep_all = TRUE)

# Join ArcGIS attributes onto NDVI observations
dc_data <- ndvi_dat %>%
  left_join(arc_dc_trees_clean, by = "UID")

# If join created duplicate SpCode names, keep the NDVI one
if ("SpCode.x" %in% names(dc_data)) dc_data$SpCode <- dc_data$SpCode.x
if ("SpCode.y" %in% names(dc_data) && !("SpCode" %in% names(dc_data))) dc_data$SpCode <- dc_data$SpCode.y

# Add year
dc_data$year <- as.integer(substr(as.character(dc_data$im_date), 1, 4))

head(dc_data)


#clean

dc_data <- dc_data %>%
  select(-SpCode.x, -SpCode.y)

head(dc_data)
names(dc_data)