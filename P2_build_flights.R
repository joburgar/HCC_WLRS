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

HELI_DIR <- "raw_data"

Sys.setenv(TZ = "America/Vancouver")
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

get_operator <- function(file){
  case_when(
    str_detect(basename(file), "Bearpaw") ~ "Bearpaw",
    str_detect(basename(file), "MWHS") ~ "MWHS",
    str_detect(basename(file), "Crescent") ~ "Crescent",
    str_detect(basename(file), "Silvertip") ~ "Silvertip",
    TRUE ~ "Unknown")}

############################################################
# COMBINE ALL HELI POINTS
############################################################

heli <- map_dfr(files[c(3,27)],function(f){message("Reading: ", basename(f))
  x <- st_read(f, quiet = TRUE)
  x$source_gpkg <- basename(f)
  x$operator <- get_operator(f)
  x})

############################################################
# QA & TESTING
############################################################
# nrow(heli)
# names(heli)
# summary(heli$time)
# 
# class(heli$time)
# attr(heli$time, "tzone")
# table(is.na(heli$time))
# length(unique(heli$track_seg_id))
# 
# heli |>st_drop_geometry() |>count(track_seg_id, sort = TRUE) |>head(20)
# 
# heli |>  st_drop_geometry() |> mutate( date_local = as.Date(time)) |>
#   count(date_local)
# range(as.Date(heli$time))
# 
# heli_dt <- heli |>st_drop_geometry() |>arrange(time)
# 
# heli_dt |>mutate(gap_s = as.numeric(difftime(time,lag(time),units = "secs"))) |>
#   summarise(
#     min_gap = min(gap_s, na.rm = TRUE),
#     median_gap = median(gap_s, na.rm = TRUE),
#     p90_gap = quantile(gap_s, 0.9, na.rm = TRUE),
#     p99_gap = quantile(gap_s, 0.99, na.rm = TRUE),
#     max_gap = max(gap_s, na.rm = TRUE)
#   )
# 
# heli_dt |>mutate(gap_min = as.numeric(difftime(time,lag(time),units = "mins"))) |>
#   filter(gap_min > 10) |>
#   select(time, gap_min) |>
#   head(20)
# 
# 
# heli_dt <- heli |>st_drop_geometry() |>arrange(time) |>
#   mutate(gap_min = as.numeric(difftime(time,lag(time),units = "mins")),
#     new_segment =is.na(gap_min) |gap_min > 10,
#     segment_id =cumsum(new_segment))
# 
# heli_dt |> count(segment_id, sort = TRUE)
# n_distinct(heli_dt$segment_id)
# 
# seg_summary <- heli_dt |> count(segment_id, name = "n_pts")
# 
# summary(seg_summary$n_pts)
# sum(seg_summary$n_pts < 10)
# sum(seg_summary$n_pts < 20)
# 
# quantile(seg_summary$n_pts,c(.01,.05,.25,.5,.75,.9,.95,.99))
# 
# seg_stats <- heli_dt |>group_by(segment_id) |>
#   summarise(n_pts = n(),duration_min =as.numeric(difftime(max(time),min(time),units = "mins")))
# 
# summary(seg_stats$duration_min)
# 
# seg_stats |> arrange(duration_min) |> head(20)

############################################################
# STANDARDIZE TIME
############################################################

heli <- heli %>% mutate( # GPKG timestamps are UTC but timezone attribute is lost
  t_utc = force_tz(time, "UTC"),local_time = with_tz(t_utc,TIMEZONE),
  date_local = as.Date(local_time),
  season = year(date_local) + if_else(month(date_local) == 12,1,0))

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
# BUILD FLIGHT SEGMENTS
############################################################

heli_dt <- heli %>% st_drop_geometry() %>%
  mutate(tail = if_else(is.na(tail) | tail == "","unknown",tail)) %>%
  arrange(operator,tail,date_local,t_utc) %>%
  group_by(operator,tail,date_local) %>%
  mutate(gap_min =as.numeric(difftime(t_utc,lag(t_utc),units = "mins")),
         step_m =sqrt((x - lag(x))^2 +(y - lag(y))^2),
         new_file =source_file != lag(source_file),
         break_track =is.na(gap_min) |gap_min > MAX_GAP_MIN |step_m > MAX_JUMP_KM * 1000,
         segment_id =cumsum(break_track)) %>%
  ungroup()


############################################################
# DISTANCE / SPEED QA
############################################################

heli_dt <- heli_dt %>%
  group_by(operator,tail) %>%
  mutate(dt_h = as.numeric(difftime(t_utc,lag(t_utc),units = "hours")),
         dist_km =sqrt((x - lag(x))^2 +(y - lag(y))^2) / 1000,
         speed_kmh =if_else(dt_h > 0,dist_km / dt_h,NA_real_),
         good_jump =is.na(dist_km) |dist_km <= MAX_JUMP_KM,
         airborne =speed_kmh >= AIRBORNE_KMH) %>%
  ungroup()


############################################################
# CREATE AIRCRAFT-DAY TABLE
############################################################
flying_days <- heli_dt %>%
  mutate(tail =if_else(is.na(tail) | tail == "","unknown",tail)) %>%
  group_by(operator,tail,date_local) %>%
  summarise(n_points = n(),t_first_utc =min(t_utc,na.rm = TRUE),
    t_last_utc =max(t_utc,na.rm = TRUE),
    season =first(season),
    .groups = "drop")

# flying_days %>% count(operator, tail, sort = TRUE)
# 
# flying_days %>% mutate(duration_h =as.numeric(difftime(t_last_utc,t_first_utc,units = "hours"))) %>%
#   summarise(
#     min_dur = min(duration_h),
#     median_dur = median(duration_h),
#     p90_dur = quantile(duration_h, 0.9),
#     max_dur = max(duration_h)
#   )

############################################################
# QA
############################################################

cat("\nAircraft-days:",nrow(flying_days),"\n")

flying_days %>%
  count(operator,sort = TRUE)

summary(flying_days$n_points)

############################################################
# FLIGHT GROUPS
############################################################

flight_groups <- heli_dt %>%
  group_by(operator,tail,date_local,season) %>%
  summarise(n_points = n(),t_first = min(t_utc),t_last = max(t_utc),.groups = "drop")

############################################################
# QA
############################################################

print(head(flight_groups))

summary(flight_groups$n_points)

segments <- heli_dt %>%group_by(operator,tail,date_local,season,segment_id) %>%
  summarise(n_points = n(),first_utc = min(t_utc),last_utc  = max(t_utc),
            duration_h =as.numeric(difftime(max(t_utc),min(t_utc),units = "hours")),
            airborne_pts =sum(airborne, na.rm = TRUE),
            .groups = "drop")

summary(segments$duration_h)
summary(segments$n_points)
summary(segments$airborne_pts)
quantile(segments$duration_h,c(.5,.75,.9,.95,.99))

############################################################
# THIN FLIGHT TRACKS
############################################################

heli_thin <- heli_dt %>% group_by(operator,tail,date_local,segment_id) %>%
  mutate(elapsed_s =as.numeric(difftime(t_utc,first(t_utc),units = "secs")),
         thin_bin =floor(elapsed_s / THIN_S)) %>%
  group_by(operator,tail,date_local,segment_id,thin_bin) %>%
  slice(1) %>%
  ungroup()

inherits(heli, "sf")
inherits(heli_dt, "sf")
inherits(heli_thin, "sf")

############################################################
# QA
############################################################

cat("\nOriginal points:",nrow(heli_dt),"\n")
cat("Thinned points:",nrow(heli_thin),"\n")
cat("Reduction:",round(100 * (1 - nrow(heli_thin)/nrow(heli_dt)),1),"%\n")

nrow(heli_dt)
nrow(heli_thin)

summary(heli_thin$speed_kmh)

# heli_thin %>%summarise(
#     median_speed = median(speed_kmh, na.rm = TRUE),
#     p90_speed = quantile(speed_kmh, 0.9, na.rm = TRUE),
#     airborne_prop = mean(airborne, na.rm = TRUE))
# 
# heli_dt %>%summarise(
#     median_speed = median(speed_kmh, na.rm = TRUE),
#     p90_speed = quantile(speed_kmh, 0.9, na.rm = TRUE),
#     airborne_prop = mean(airborne, na.rm = TRUE))

heli_thin <- heli_thin %>%filter(good_jump)

############################################################
# BUILD FLIGHT LINES
############################################################

flight_lines <- heli_thin %>% group_by(operator,tail,date_local,season,segment_id) %>%
  summarise(n_points = n(),first_utc = min(t_utc),last_utc  = max(t_utc),
            duration_h =as.numeric(difftime(max(t_utc),min(t_utc),units = "hours")),
            airborne_pts =sum(airborne,na.rm = TRUE),
            do_union = FALSE,
            .groups = "drop")

############################################################
# KEEP VALID FLIGHT SEGMENTS
############################################################

flight_lines <- flight_lines %>% filter(n_points >= 2)

############################################################
# QA
############################################################

cat("\nFlight segments:",nrow(flight_lines),"\n")

summary(flight_lines$n_points)
summary(flight_lines$duration_h)
summary(flight_lines$airborne_pts)

############################################################
# BUILD FLIGHT LINES
############################################################

heli_thin_sf <- st_as_sf(heli_thin,coords = c("x", "y"),crs = st_crs(heli))

flight_lines <- heli_thin_sf %>% group_by(operator,tail,date_local,season,segment_id) %>%
  summarise(n_points = n(),first_utc = min(t_utc),last_utc = max(t_utc),
            duration_h =as.numeric(difftime(max(t_utc),min(t_utc),units = "hours")),
            airborne_pts =sum(airborne,na.rm = TRUE),
            mean_alt_m =mean(alt_m,na.rm = TRUE),
            n_missing_alt =sum(is.na(alt_m)),
            do_union = FALSE,
            .groups = "drop") %>%
  filter(n_points >= 2)

table(st_geometry_type(flight_lines))
nrow(flight_lines)
summary(flight_lines$n_points)
summary(flight_lines$duration_h)

############################################################
# EARLY SEGMENTS
#
# Steve's script keeps segments occurring before
# the daily anchor cutoff (~07:10 local)
############################################################

EARLY_CUTOFF_HOUR <- 7 + 10/60

early_segments <- flight_lines %>% filter(hour(with_tz(first_utc,TIMEZONE)) < EARLY_CUTOFF_HOUR)

############################################################
# QA
############################################################

nrow(early_segments)

summary(early_segments$duration_h)

############################################################
# FLIGHTS_THIN
############################################################

flights_thin <- flight_lines %>%filter(airborne_pts > 0)

############################################################
# QA
############################################################

nrow(flight_lines)
nrow(flights_thin)
summary(flights_thin$airborne_pts)

nrow(early_segments)
nrow(flights_thin)
table(st_geometry_type(flights_thin))

############################################################
# FLIGHTS_NO_ALT
############################################################

flights_no_alt <- flight_lines %>% filter(n_missing_alt > 0)

summary(flight_lines$mean_alt_m)
summary(flight_lines$n_missing_alt)
nrow(flights_no_alt)

############################################################
# ALTITUDE QA
############################################################

flight_lines <- heli_thin_sf %>%
  group_by(
    operator,
    tail,
    date_local,
    season,
    segment_id
  ) %>%
  summarise(
    n_points = n(),
    
    first_utc = min(t_utc),
    
    last_utc = max(t_utc),
    
    duration_h =
      as.numeric(
        difftime(
          max(t_utc),
          min(t_utc),
          units = "hours"
        )
      ),
    
    airborne_pts =
      sum(
        airborne,
        na.rm = TRUE
      ),
    
    mean_alt_m =
      mean(
        alt_m,
        na.rm = TRUE
      ),
    
    n_missing_alt =
      sum(
        is.na(alt_m)
      ),
    
    do_union = FALSE,
    
    .groups = "drop"
  ) %>%
  filter(
    n_points >= 2
  )

############################################################
# FLIGHTS WITH MISSING ALTITUDE
############################################################

flights_no_alt <- flight_lines %>%
  filter(
    n_missing_alt > 0
  )

############################################################
# EARLY SEGMENTS
############################################################

EARLY_CUTOFF_HOUR <- 7 + 10/60

early_segments <- flight_lines %>%
  filter(
    hour(
      with_tz(
        first_utc,
        TIMEZONE
      )
    ) < EARLY_CUTOFF_HOUR
  )

############################################################
# THINNED FLIGHT SEGMENTS
############################################################

flights_thin <- flight_lines %>%
  filter(
    airborne_pts > 0
  )

############################################################
# QA
############################################################

cat(
  "\nFlight segments:",
  nrow(flight_lines),
  "\n"
)

cat(
  "Airborne segments:",
  nrow(flights_thin),
  "\n"
)

cat(
  "Missing altitude segments:",
  nrow(flights_no_alt),
  "\n"
)

cat(
  "Early segments:",
  nrow(early_segments),
  "\n"
)

############################################################
# PLACEHOLDER FOR NEXT SECTION
############################################################

# TODO
# Revisit Steve's build_tracks() logic
# Compare:
#   nrow(flight_lines)
#   nrow(flights_thin)
#   nrow(early_segments)
#   nrow(flights_no_alt)
# against Python output
#
# Then reproduce:
#   Aircraft-day tracks = 4846
#   Tracks entering study area = 3537
#   Operator-specific thinning
#       time_120s
#       every_52th_vertex
#       every_119th_vertex
#
# before exporting flights.gpkg