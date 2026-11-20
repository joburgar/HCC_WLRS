# Copyright 2019 Province of British Columbia
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at 
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

############################################################
# P2a_build_helibase.R
# Cleaning raw heli data and saving for further processing
############################################################

library(sf)
library(tidyverse)

############################################################
# SETTINGS
############################################################

HELI_DIR <- "raw_data"

Sys.setenv(TZ = "America/Vancouver")

SEASON_MONTHS <- c(12,1,2,3,4)

############################################################
# READ HELI DATA FILES
############################################################

flight_dir <- "./raw_data"

files <- list.files(flight_dir,pattern = "\\.gpkg$",recursive = TRUE,full.names = TRUE)

tibble(file = basename(files),path = files) |>
  arrange(file)

get_operator <- function(file){
  case_when(
    str_detect(basename(file), "Bearpaw") ~ "Bearpaw",
    str_detect(basename(file), "MWHS") ~ "MWHS",
    str_detect(basename(file), "Crescent") ~ "Crescent",
    str_detect(basename(file), "Silvertip") ~ "Silvertip",
    TRUE ~ "Unknown")}

############################################################
# COMBINE ALL HELI POINTS & TRIM DOWN
############################################################

heli <- map_dfr(files,function(f){message("Reading: ", basename(f))
  x <- st_read(f, quiet = TRUE)
  x$source_gpkg <- basename(f)
  x$operator <- get_operator(f)
  x})

names(heli)

heli <- heli %>%
  select(any_of(c("time","ele","alt_m","event","source_file","source_gpkg",
                  "operator","flight_date","tail","geometry")))

############################################################
# STANDARDIZE TIME
############################################################

heli <- heli %>%
  mutate(# raw timestamps are stored in BC local time
    local_time = force_tz(time,"America/Vancouver"),
    # convert to UTC for matching with caribou data later
    t_utc = with_tz(local_time,"UTC"),
    date_local = as.Date(local_time),
    season =year(date_local) +if_else(month(date_local) == 12,1,0))

############################################################
# FILTER TO HELISKI SEASON
############################################################

heli <- heli %>% filter(month(date_local) %in% SEASON_MONTHS)

############################################################
# COORDINATES
############################################################

coords <- st_coordinates(heli)

heli <- heli %>% mutate(x = coords[,1],y = coords[,2])

############################################################
# SAVE INTERMEDIATE FILE
############################################################

saveRDS(heli,"fl_processed/heli_base.rds")
