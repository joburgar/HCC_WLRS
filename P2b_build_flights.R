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
# P2b_build_flights.R
#
# Translation of build_flights.py (created by Steve Wilson)
#
# Bringing in heli data cleaned using P2a_build_heli_base.R
############################################################

library(sf)
library(tidyverse)
# library(data.table)
library(bcdata)
library(janitor)

############################################################
# SETTINGS
############################################################

MAX_GAP_MIN <- 10

MAX_JUMP_KM <- 15

AIRBORNE_KMH <- 60

THIN_S <- 120

run_id <- paste0(
  "g", MAX_GAP_MIN,
  "_j", MAX_JUMP_KM,
  "_air", AIRBORNE_KMH,
  "_t", THIN_S,
  "_c0710"
)

############################################################
# READ HELI DATA
############################################################

heli <- readRDS("fl_processed/heli_base.rds")

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

glimpse(heli_dt)

############################################################
# THIN FLIGHT TRACKS
############################################################

heli_thin <- heli_dt %>% group_by(operator,tail,date_local,segment_id) %>%
  mutate(elapsed_s =as.numeric(difftime(t_utc,first(t_utc),units = "secs")),
         thin_bin =floor(elapsed_s / THIN_S)) %>%
  group_by(operator,tail,date_local,segment_id,thin_bin) %>%
  slice(1) %>%
  ungroup()

# inherits(heli, "sf")
# inherits(heli_dt, "sf")
# inherits(heli_thin, "sf")

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

heli_thin <- heli_thin %>% filter(good_jump)

heli_thin_sf <- st_as_sf(heli_thin,coords = c("x", "y"),crs = st_crs(heli))

# QA
# inherits(heli_thin_sf, "sf"); table(st_geometry_type(heli_thin_sf))

############################################################
# BUILD FLIGHT LINES
############################################################

flight_lines <- heli_thin_sf %>%
  group_by(operator,tail,date_local,season,segment_id) %>%
  summarise(n_points = n(),first_utc = min(t_utc),last_utc = max(t_utc),
            duration_h =as.numeric(difftime(max(t_utc),min(t_utc),units = "hours")),
    airborne_pts =sum(airborne,na.rm = TRUE),
    mean_alt_m =mean(alt_m,na.rm = TRUE),
    n_missing_alt =sum(is.na(alt_m)),
    do_union = FALSE,
    .groups = "drop") %>%
  st_cast("LINESTRING") %>%
  filter(n_points >= 2)

flight_lines <- flight_lines %>% mutate(length_km =as.numeric(st_length(geometry)) / 1000)

table(st_geometry_type(flight_lines))
nrow(flight_lines)
summary(flight_lines$n_points)
summary(flight_lines$duration_h)

rm(heli_thin_sf) # housekeeping

############################################################
# EARLY SEGMENTS
## Steve's script keeps segments occurring before
# the daily anchor cutoff (~07:10 local)
############################################################

EARLY_CUTOFF_HOUR <- 7 + 10/60

early_segments <- flight_lines %>% 
  mutate(first_hour = hour(with_tz(first_utc, TIMEZONE)) + 
           minute(with_tz(first_utc, TIMEZONE)) / 60 +
           second(with_tz(first_utc, TIMEZONE)) / 3600) %>%
  filter(first_hour < EARLY_CUTOFF_HOUR)

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
# QA
############################################################

cat("\nFlight segments:",nrow(flight_lines),"\n")
cat("Airborne segments:",nrow(flights_thin),"\n")
cat("Missing altitude segments:",nrow(flights_no_alt),"\n")
cat("Early segments:",nrow(early_segments),"\n")

# summary(hour(with_tz(flight_lines$first_local,TIMEZONE)))
# table(hour(with_tz(flight_lines$first_local,TIMEZONE)))

############################################################
# LOAD STUDY AREA
############################################################

herds <- bcdc_get_data("2b217585-f48d-4d9f-b7ba-746909ac35ca") %>% clean_names() %>% filter(ecotype == "Mountain")
study_area <- herds %>%st_buffer(20000)

ggplot() + geom_sf(data = study_area,fill = "grey90",colour = "black") +
  geom_sf(data = flight_lines,aes(colour = operator),linewidth = 0.3,alpha = 0.6) +
  theme_bw()

ggplot() +
  geom_sf(data = study_area,fill = "grey95") +
  geom_sf(data = flight_lines %>%slice_sample(n = 50),aes(colour = factor(segment_id)),linewidth = 0.5,show.legend = FALSE) +
  theme_bw()

############################################################
# STUDY AREA INTERSECTION
############################################################

track_in_study <- st_intersects(flight_lines,study_area,sparse=FALSE)
flight_lines$in_study_area <-apply(track_in_study,1,any)
table(flight_lines$in_study_area)

study_lookup <- flight_lines %>% st_drop_geometry() %>%
  select(operator,tail,date_local,season,segment_id,in_study_area)

############################################################
# TRACK ATTRIBUTES
############################################################

tracks <- heli_dt %>%
  group_by(operator,tail,date_local,season,segment_id) %>%
  summarise(n_points = n(),t_first = min(t_utc),t_last = max(t_utc),
            duration_h = as.numeric(difftime(max(t_utc),min(t_utc),units = "hours")),
            airborne_pts = sum(airborne,na.rm = TRUE),
            median_speed_kmh = if(all(is.na(speed_kmh)))
              {NA_real_} else {median(speed_kmh, na.rm = TRUE)},
            max_speed_kmh = if(all(is.na(speed_kmh)))
              {NA_real_} else {max(speed_kmh, na.rm = TRUE)},
            mean_alt_m = mean(alt_m,na.rm = TRUE),
            source_files = paste(unique(source_file),collapse = "; "),
            .groups = "drop")

track_lookup <- flight_lines %>% st_drop_geometry() %>%
  select(operator,tail,date_local,season,segment_id,length_km,in_study_area)

tracks <- tracks %>%
  left_join(track_lookup,by = c("operator","tail","date_local","season","segment_id"))

tracks <- tracks %>%
  mutate(in_study_area = replace_na(in_study_area,FALSE),
         airborne_track =airborne_pts > 0,
         airborne_study_track = airborne_track &in_study_area)

############################################################
# FIRST AIRBORNE FIX
############################################################

heli_dt <- heli_dt %>%left_join(study_lookup,by = c("operator","tail","date_local","season","segment_id"))

first_airborne <- heli_dt %>% filter(airborne) %>%
  group_by(operator,tail,date_local) %>%
  summarise(first_airborne_utc = min(t_utc),.groups = "drop")

############################################################
# FIRST AIRBORNE FIX IN STUDY AREA
############################################################

first_airborne_study <- heli_dt %>%
  filter(airborne,in_study_area) %>%
  group_by(operator,tail,date_local) %>%
  summarise(first_airborne_study_utc = min(t_utc),.groups = "drop")

nrow(first_airborne)
nrow(first_airborne_study)
head(first_airborne_study)

############################################################
# THINNING METADATA
############################################################
time_resolution <- heli_dt %>%
  group_by(operator,tail,date_local,season,segment_id) %>%
  summarise(time_res_s =median(dt_h * 3600,na.rm = TRUE),.groups = "drop")

tracks <- tracks %>%
  left_join(time_resolution,by = c("operator","tail","date_local","season","segment_id")) %>%
  mutate(thin_method =case_when(time_res_s >= THIN_S ~ "none",TRUE ~ "time_120s"))

############################################################
# TRACK METADATA
############################################################

tracks <- tracks %>%
  left_join(first_airborne,by = c("operator","tail","date_local")) %>%
  left_join(first_airborne_study,by = c("operator","tail","date_local"))

############################################################
# TRACK SUMMARY
############################################################

track_summary <- tracks %>%
  summarise(n_tracks = n(),n_in_study =sum(in_study_area,na.rm = TRUE),
            n_airborne =sum(airborne_track,na.rm = TRUE),
            n_airborne_study =sum(airborne_study_track,na.rm = TRUE),
            median_length_km =median(length_km,na.rm = TRUE),
            median_duration_h =median(duration_h,na.rm = TRUE))

track_summary

tracks %>%count(operator,airborne_study_track)

############################################################
# FINAL TRACK FILTERS
############################################################

tracks_final <- tracks %>% filter(airborne_track,in_study_area)

############################################################
# QA
############################################################

cat("\nTracks entering study area:",sum(tracks$in_study_area),"of",nrow(tracks),"\n")
cat("Airborne tracks in study area:",nrow(tracks_final),"\n")

nrow(tracks_final)
summary(tracks_final$length_km)
summary(tracks_final$duration_h)


############################################################
# FIRST AIRBORNE LOCAL TIME
############################################################

first_airborne_study <- first_airborne_study %>%
  mutate(first_airborne_local =with_tz(first_airborne_study_utc,TIMEZONE),
         first_airborne_hour =hour(first_airborne_local) +minute(first_airborne_local)/60 +second(first_airborne_local)/3600)

############################################################
# QA
############################################################

summary(first_airborne_study$first_airborne_hour)
quantile(first_airborne_study$first_airborne_hour,c(.05,.25,.5,.75,.95))

ggplot(first_airborne_study,aes(first_airborne_hour)) +
  geom_histogram(binwidth = 0.25)

############################################################
# CREATE AIRCRAFT-DAY TABLE
############################################################
flying_days <- heli_dt %>%
  mutate(tail =if_else(is.na(tail) | tail == "","unknown",tail)) %>%
  group_by(operator,tail,date_local) %>%
  summarise(n_points = n(),first_utc =min(t_utc,na.rm = TRUE),
            last_utc =max(t_utc,na.rm = TRUE),
            season =first(season),
            .groups = "drop")

# flying_days %>% count(operator, tail, sort = TRUE)
# 
# flying_days %>% mutate(duration_h =as.numeric(difftime(last_utc,first_utc,units = "hours"))) %>%
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
flying_days %>% count(operator,sort = TRUE)
summary(flying_days$n_points)

############################################################
# TRACK SUMMARY TABLE
############################################################

tracks <- tracks %>%
  select(operator,tail,date_local,season,segment_id,n_points,
         t_first,t_last,first_airborne_utc,first_airborne_study_utc,
         duration_h,length_km,
         airborne_pts,airborne_track,in_study_area,airborne_study_track,
         median_speed_kmh,max_speed_kmh,mean_alt_m,
         thin_method,time_res_s,
         source_files)

glimpse(tracks)

track_summary

output_summary <- tibble(
  object = c("flight_lines","flights_thin","flights_no_alt","early_segments",
             "flying_days","tracks","tracks_final"),
  n = c(nrow(flight_lines),nrow(flights_thin),nrow(flights_no_alt),nrow(early_segments),
        nrow(flying_days),nrow(tracks),nrow(tracks_final)))
output_summary
 
# ###############
# steve_tracks <- flight_lines %>%
#   group_by(operator,tail,date_local,season) %>%
#   summarise(
#     n_segments = n(),
#     n_points = sum(n_points, na.rm = TRUE),
#     length_km = sum(length_km, na.rm = TRUE),
#     duration_h = sum(duration_h, na.rm = TRUE),
#     airborne_pts = sum(airborne_pts, na.rm = TRUE),
#     in_study_area = any(in_study_area),
#     do_union = FALSE,
#     .groups = "drop")
# 
# nrow(steve_tracks) # 3889 vs 4846 from Steve's output

############################################################
# EXPORT OUTPUTS
############################################################

output_dir <- "outputs"
if(!dir.exists(output_dir)){dir.create(output_dir)}

write_csv(tracks,file.path(output_dir,paste0("tracks_", run_id, ".csv")))

gpkg_out <- file.path(output_dir,paste0("flights_",run_id,".gpkg"))
st_write(flight_lines,gpkg_out,layer = "flights",quiet = TRUE)
st_write(flights_thin,gpkg_out,layer = "flights_thin",append = TRUE,quiet = TRUE)
st_write(flights_no_alt,gpkg_out,layer = "flights_no_alt",append = TRUE,quiet = TRUE)
st_write(early_segments,gpkg_out,layer = "early_segments",append = TRUE,quiet = TRUE)
