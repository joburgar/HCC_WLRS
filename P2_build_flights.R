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
# P2_build_flights.R
#
# Translation of build_flights.py (created by Steve Wilson)
#
# Bringing in heli data cleaned using 3.Read_heli_activity.R
############################################################

library(sf)
library(tidyverse)
library(data.table)

############################################################
# SETTINGS
############################################################

HELI_DIR <- "fl_processed"

TIMEZONE <- "America/Vancouver"

SEASON_MONTHS <- c(12,1,2,3,4)

MAX_GAP_MIN <- 10
MAX_JUMP_KM <- 15
AIRBORNE_KMH <- 60

THIN_S <- 120

############################################################
# READ HELI DATA FILES
############################################################

flight_dir <- "./raw_data"

files <- list.files(flight_dir,recursive = TRUE,full.names = TRUE)

tibble(file = basename(files),path = files) |>
  arrange(file)

############################################################
# COMBINE ALL HELI POINTS
############################################################

heli <- map_dfr(gpkg_files[1], # to test for one gpkg
  function(f){
    message("Reading: ", basename(f))
    x <- st_read(f,quiet = TRUE)
    x$source_gpkg <- basename(f)
    x
  }
)

############################################################
# QA & TESTING
############################################################

heli <- st_read("raw_data/Bearpaw_flightlines_GPX_2025.gpkg")

nrow(heli)
names(heli)
summary(heli$time)

class(heli$time)
attr(heli$time, "tzone")
table(is.na(heli$time))
length(unique(heli$track_seg_id))

heli |>st_drop_geometry() |>count(track_seg_id, sort = TRUE) |>head(20)

heli |>  st_drop_geometry() |> mutate( date_local = as.Date(time)) |>
  count(date_local)
range(as.Date(heli$time))

heli_dt <- heli |>st_drop_geometry() |>arrange(time)

heli_dt |>mutate(gap_s = as.numeric(difftime(time,lag(time),units = "secs"))) |>
  summarise(
    min_gap = min(gap_s, na.rm = TRUE),
    median_gap = median(gap_s, na.rm = TRUE),
    p90_gap = quantile(gap_s, 0.9, na.rm = TRUE),
    p99_gap = quantile(gap_s, 0.99, na.rm = TRUE),
    max_gap = max(gap_s, na.rm = TRUE)
  )

heli_dt |>mutate(gap_min = as.numeric(difftime(time,lag(time),units = "mins"))) |>
  filter(gap_min > 10) |>
  select(time, gap_min) |>
  head(20)


heli_dt <- heli |>st_drop_geometry() |>arrange(time) |>
  mutate(gap_min = as.numeric(difftime(time,lag(time),units = "mins")),
    new_segment =is.na(gap_min) |gap_min > 10,
    segment_id =cumsum(new_segment))

heli_dt |> count(segment_id, sort = TRUE)
n_distinct(heli_dt$segment_id)

seg_summary <- heli_dt |> count(segment_id, name = "n_pts")

summary(seg_summary$n_pts)
sum(seg_summary$n_pts < 10)
sum(seg_summary$n_pts < 20)


quantile(seg_summary$n_pts,c(.01,.05,.25,.5,.75,.9,.95,.99))

seg_stats <- heli_dt |>group_by(segment_id) |>
  summarise(n_pts = n(),duration_min =as.numeric(difftime(max(time),min(time),units = "mins")))

summary(seg_stats$duration_min)

seg_stats |> arrange(duration_min) |> head(20)


############################################################
# BUILD FLIGHT SEGMENTS
############################################################

heli_dt <- heli |> st_drop_geometry() |> arrange(time) |>
  mutate(gap_min = as.numeric(difftime(time,lag(time),units = "mins")),
      new_segment = is.na(gap_min) | gap_min > 10,
      segment_id =cumsum(new_segment))

head(heli_dt)
head(heli %>%st_drop_geometry() %>%select(time))

summary(heli$flight_date)
head(heli %>%st_drop_geometry() %>%select(time, flight_date))
range(format(heli$time,"%H:%M:%S"))
heli %>% st_drop_geometry() %>% mutate(hour = hour(time)) %>% count(hour)

############################################################
# STANDARDIZE TIME
############################################################

heli <- heli %>%
  mutate(t_utc = time,date_local =as.Date(with_tz(t_utc,TIMEZONE)),
    season =year(date_local) +if_else(month(date_local) == 12,1,0))

############################################################
# FILTER TO HELISKI SEASON
############################################################

heli <- heli %>%
  filter(
    month(date_local) %in%
      SEASON_MONTHS
  )

############################################################
# COORDINATES
############################################################

coords <- st_coordinates(heli)

heli <- heli %>%
  mutate(
    x = coords[,1],
    y = coords[,2]
  )

############################################################
# CREATE AIRCRAFT-DAY TABLE
############################################################

flying_days <- heli %>%
  st_drop_geometry() %>%
  mutate(
    tail =
      if_else(
        is.na(tail),
        "unknown",
        tail
      )
  ) %>%
  group_by(
    operator,
    tail,
    date_local
  ) %>%
  summarise(
    
    n_points = n(),
    
    t_first_utc =
      min(
        t_utc,
        na.rm = TRUE
      ),
    
    t_last_utc =
      max(
        t_utc,
        na.rm = TRUE
      ),
    
    season =
      first(season),
    
    .groups = "drop"
  )

############################################################
# QA
############################################################

cat(
  "\nAircraft-days:",
  nrow(flying_days),
  "\n"
)

flying_days %>%
  count(
    operator,
    sort = TRUE
  )

summary(
  flying_days$n_points
)

############################################################
# PREP FOR TRACK BUILDING
############################################################

heli_dt <- heli %>%
  st_drop_geometry() %>%
  arrange(
    operator,
    tail,
    date_local,
    t_utc
  ) %>%
  as.data.table()

############################################################
# FLIGHT GROUPS
############################################################

flight_groups <- heli_dt[
  ,
  .(
    n_points = .N,
    t_first = min(t_utc),
    t_last = max(t_utc)
  ),
  by = .(
    operator,
    tail,
    date_local,
    season
  )
]

############################################################
# QA
############################################################

print(head(flight_groups))

summary(
  flight_groups$n_points
)