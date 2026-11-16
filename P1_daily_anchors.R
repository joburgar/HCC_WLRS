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
# P1_daily_anchors.R
#
# Translation of daily_anchors.py (created by Steve Wilson)
#
# Creates one anchor per collar-day:
# Last caribou fix prior to local cutoff time
#
# Joanna Burgar
############################################################

library(sf)
library(tidyverse)
library(data.table)

############################################################
# USER SETTINGS
############################################################

INPUT_FILE <- "//sfp.idir.bcgov/home/perpetual/HCC_flight_data/HCC_data_share_14Sept2026.shp"

OUTPUT_STEM <- "daily_anchors"

COLLAR_FIELD <- "COLLAR_"
TIME_FIELD   <- "datetim"

Sys.setenv(TZ = "America/Vancouver")
TIMEZONE <- "America/Vancouver"

CUTOFF_HOUR <- 7 + 10/60     # 07:10

MAX_AGE_H <- 12

SEASON_MONTHS <- c(12,1,2,3,4) # Dec - April

############################################################
# FUNCTIONS
############################################################

read_fixes <- function(path,
                       collar_field,
                       time_field){
  
  cat("\nReading telemetry...\n")
  
  g <- st_read(path, quiet = TRUE)
  
  if(is.na(st_crs(g))){stop("Input layer has no CRS")}
  
  g <- g %>%
    rename(collar = all_of(collar_field),time_raw = all_of(time_field))
  
  g$collar <- as.integer(round(as.numeric(g$collar)))
  
  # Parse time as UTC
  g$t <- ymd_hms(as.character(g$time_raw),tz = "UTC")
  
  if(any(is.na(g$t))){stop("Failed to parse one or more timestamps")}
  
  g <- g %>% distinct(collar, t, .keep_all = TRUE) %>% arrange(collar, t)
  
  g
}

############################################################

calc_fix_interval <- function(fixes){
  fixes %>%
    st_drop_geometry() %>%
    arrange(collar, t) %>%
    group_by(collar) %>%
    mutate(dt_h = as.numeric(difftime(t,lag(t),units = "hours"))) %>%
    summarise(fix_interval_h =round(median(dt_h, na.rm = TRUE)),.groups = "drop")
}

# sched <- calc_fix_interval(fixes_test)

############################################################

create_candidate_days <- function(fixes){
  
  local_dates <- fixes %>%
    st_drop_geometry() %>%
    mutate(date_local =as.Date(with_tz(t, TIMEZONE))) %>%
    select(collar, date_local)
  
  d1 <- distinct(local_dates)
  
  d2 <- d1 %>%
    mutate(date_local = date_local + 1)
  
  bind_rows(d1, d2) %>%
    distinct() %>%
    filter(month(date_local) %in% SEASON_MONTHS)
}

# days <- create_candidate_days(fixes_test)
# head(days)
# nrow(days)
# range(days$date_local)

############################################################

build_cutoffs <- function(days){
  days %>%
    mutate(
      cutoff_local =ymd_hms(paste(date_local,sprintf("%02d:%02d:00",
              floor(CUTOFF_HOUR),round((CUTOFF_HOUR %% 1) * 60))),tz = TIMEZONE),
      
      cutoff_utc =with_tz(cutoff_local,"UTC")) %>%
    arrange(cutoff_utc)
}

# days <- build_cutoffs(days)
# days %>%
#   select(collar,date_local,cutoff_local,cutoff_utc) %>%
#   head(10)

############################################################

# fx <- fixes_test
# coords <- st_coordinates(fx)
# 
# fx <- fx |>mutate(x = coords[,1],y = coords[,2]) |>
#   st_drop_geometry()
# 
# fx_dt <- as.data.table(
#   fx |>select(collar, t, x, y))
# 
# days_dt <- as.data.table(days)
# 
# setkey(fx_dt, collar, t)
# 
# anc <- fx_dt[days_dt,on = .(collar,t = cutoff_utc),roll = TRUE]
# sum(is.na(anc$t))
# str(anc)

build_anchors <- function(fixes){
  
  sched <- calc_fix_interval(fixes)
  
  days <- create_candidate_days(fixes)
  
  days <- build_cutoffs(days)
  
 
  ##########################################################
  # PREP DATA.TABLES
  ##########################################################
  
  fx <- fixes
  
  coords <- st_coordinates(fx)
  
  fx <- fx %>%
    mutate(x = coords[,1],y = coords[,2]) %>%
    st_drop_geometry()
  
  fx_dt <- as.data.table(fx %>%select(collar,t,x,y))
  
  days_dt <- as.data.table(days)
  
  setkey(fx_dt, collar, t)
  
  ##########################################################
  # LAST FIX BEFORE CUTOFF
  ##########################################################
  
  ##########################################################
  # LAST FIX BEFORE CUTOFF
  ##########################################################
  
  anc <- copy(days_dt)
  
  anc[, anchor_time_utc :=
        as.POSIXct(sapply(seq_len(.N),function(i){
              hits <- fx_dt[collar == anc$collar[i] & t < anc$cutoff_utc[i]]$t
              if(length(hits) == 0) return(NA)
              max(hits)}),
          origin = "1970-01-01",
          tz = "UTC")]
  
  anc[, anchor_age_h :=
        as.numeric(difftime(cutoff_utc,anchor_time_utc,units = "hours"))]
  
  ##########################################################
  # FIRST FIX AFTER CUTOFF
  ##########################################################
  
  anc[, next_fix_utc :=as.POSIXct(sapply(seq_len(.N),function(i){
              hits <- fx_dt[collar == anc$collar[i] &t >= anc$cutoff_utc[i]]$t
              if(length(hits) == 0) return(NA)
              min(hits)}),
          origin = "1970-01-01",
          tz = "UTC")]
  
  anc[, next_fix_gap_h :=as.numeric(difftime(next_fix_utc,anchor_time_utc,units = "hours"))]
  
  ##########################################################
  # ADDITIONAL METRICS
  ##########################################################
  
  anc[, next_fix_gap_h :=as.numeric(difftime(next_fix_utc,anchor_time_utc,units = "hours"))]
  
  anc[, anchor_time_local :=format(with_tz(anchor_time_utc,TIMEZONE),"%Y-%m-%d %H:%M:%S")]
  
  ##########################################################
  # SEASON LABEL
  ##########################################################
  
  anc[, season :=year(cutoff_utc) +ifelse(month(cutoff_utc) == 12,1,0)]
  table(anc$season)
  
  nrow(anc)
  
  summary(anc$anchor_age_h)
  
  summary(anc$next_fix_gap_h)
  
  sum(is.na(anc$anchor_time_utc))
  
  sum(is.na(anc$next_fix_utc))
  
  anc <- merge(anc,sched,by="collar",all.x=TRUE)
  head(anc$fix_interval_h)
  unique(anc$fix_interval_h)
  
  n_all <- nrow(anc)
  
  # FILTERING
  
  anc <- anc %>%
    filter(
      !is.na(anchor_time_utc),
      anchor_age_h <= MAX_AGE_H,
      !is.na(next_fix_utc)
    )
  
  ##########################################################
  # ADD ANCHOR COORDINATES
  ##########################################################
  
  anchor_xy <- fx_dt[, .(
    collar,
    anchor_time_utc = t,
    x,
    y
  )]
  
  anc <- merge(
    anc,
    anchor_xy,
    by = c("collar", "anchor_time_utc"),
    all.x = TRUE
  )
  
  ##########################################################
  # FINAL FORMATTING
  ##########################################################
  
  anc <- anc %>%
    arrange(collar,date_local)
  
  return(
    list(anchors = anc,n_all = n_all))
}
############################################################

report_anchors <- function(anc, n_all){
  
  cat("\n")
  cat("========================================\n")
  cat("ANCHOR SUMMARY\n")
  cat("========================================\n")
  
  cat(
    "\nCollar-days considered:",
    n_all,
    "\n"
  )
  
  cat(
    "Anchors retained:",
    nrow(anc),
    "\n"
  )
  
  cat(
    "Collars represented:",
    length(unique(anc$collar)),
    "\n\n"
  )
  
  print(
    anc %>%
      group_by(fix_interval_h) %>%
      summarise(
        n = n(),
        collars=n_distinct(collar),
        median_age_h=
          median(anchor_age_h),
        p90_age_h=
          quantile(
            anchor_age_h,
            0.90
          ),
        .groups="drop"
      )
  )
  
  cat("\nAnchors per season\n")
  
  print(
    table(
      anc$season
    )
  )
}

############################################################

write_outputs <- function(
    anc,
    fixes_sf,
    output_stem){
  
  csv_out <-paste0(output_stem,".csv")
  
  gpkg_out <-paste0(output_stem,".gpkg")
  
  csv <- anc
  
  csv$cutoff_utc <-format(csv$cutoff_utc,"%Y-%m-%dT%H:%M:%SZ")
  
  csv$anchor_time_utc <-format(csv$anchor_time_utc,"%Y-%m-%dT%H:%M:%SZ")
  
  csv$next_fix_utc <-format(csv$next_fix_utc,"%Y-%m-%dT%H:%M:%SZ")
  
  write_csv(csv,csv_out)
  
  g <- st_as_sf(anc,coords=c("x","y"),crs=st_crs(fixes_sf))
  
  st_write(g,gpkg_out,layer="daily_anchors",delete_layer=TRUE,quiet=TRUE)
  
  cat(
    "\nWrote:\n",
    csv_out,
    "\n",
    gpkg_out,
    "\n"
  )
}

#######################################################
# TESTING WITH ONE COLLAR

# fixes <- read_fixes(INPUT_FILE,COLLAR_FIELD,TIME_FIELD)
# collar_counts <- fixes |>st_drop_geometry() |>count(collar, sort = TRUE)
# head(collar_counts, 20)
# fixes_test <- fixes |>filter(collar == 82515)
# nrow(fixes_test)
# range(fixes_test$t)
# result <- build_anchors(fixes_test)
# head(result$anchors)
# anchors <- result$anchors

############################################################
# MAIN
############################################################

fixes <- read_fixes(INPUT_FILE,COLLAR_FIELD,TIME_FIELD)

result <- build_anchors(fixes)

anchors <- result$anchors

report_anchors(anchors,result$n_all)
# ========================================
#   ANCHOR SUMMARY
# ========================================
#   
#   Collar-days considered: 59607 
# Anchors retained: 58189 
# Collars represented: 200 
# 
# fix_interval_h     n collars median_age_h p90_age_h
# <dbl> <int>   <int>     <dbl>     <dbl>
# 6     4294      16       3.15       3.16
# 8     13764     57       0.156      7.16
# 12    39767     124      3.15       3.16
# 15    347       2        1.16      10.2 
# 24    17        1        2.15       2.16
# 
# Anchors per season
# 
# 2020  2021  2022  2023  2024  2025  2026 
# 729  4338  8586 10668 10935 11736 11197 
 

write_outputs(anchors,fixes,OUTPUT_STEM)

