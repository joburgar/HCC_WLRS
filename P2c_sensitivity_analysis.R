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
# P2c_sensitivity_analysis.R
#
# Sensitivity analysis framework for heli processing
#
# Assumes:
#   P2a_build_heli_base.R has already created
#   fl_processed/heli_base.rds
#
# Runs alternative parameter sets and produces
# comparable summary tables.
############################################################

library(sf)
library(tidyverse)
library(bcdata)
library(janitor)

Sys.setenv(TZ = "America/Vancouver")

############################################################
# LOAD BASE DATA
############################################################

heli <- readRDS("fl_processed/heli_base.rds")

############################################################
# STUDY AREA
############################################################

herds <- bcdc_get_data("2b217585-f48d-4d9f-b7ba-746909ac35ca") %>%
  clean_names() %>% filter(ecotype == "Mountain")

study_area <- herds %>%st_buffer(20000)

############################################################
# PARAMETER GRID
############################################################

parameter_sets <- tribble(~run_id,~MAX_GAP_MIN, ~MAX_JUMP_KM, ~AIRBORNE_KMH, ~THIN_S, ~USE_NEW_FILE_BREAK,
  "base", 10,            15,             60,      120, FALSE,
  "gap5",5,            15,             60,      120, FALSE,
  "gap20",20,            15,             60,      120, FALSE,
  "air50",10,            15,             50,      120, FALSE,
  "air80",10,            15,             80,      120, FALSE,
  "thin60",10,            15,             60,       60, FALSE,
  "thin300",10,            15,             60,      300, FALSE,
  "jump10",10,            10,             60,      120, FALSE,
  "jump20",10,            20,             60,      120, FALSE,
  "filebreak",10,            15,             60,      120, TRUE)

############################################################
# RESULTS OBJECT
############################################################

results <- vector(
  "list",
  nrow(parameter_sets)
)

############################################################
# RUN SENSITIVITY ANALYSES
############################################################

for(i in seq_len(nrow(parameter_sets))){
  
  pars <- parameter_sets[i,]
  
  message(
    "\n==============================",
    "\nRUN: ", pars$run_id,
    "\n==============================\n"
  )
  
  ##########################################################
  # PARAMETERS
  ##########################################################
  
  MAX_GAP_MIN <- pars$MAX_GAP_MIN
  MAX_JUMP_KM <- pars$MAX_JUMP_KM
  AIRBORNE_KMH <- pars$AIRBORNE_KMH
  THIN_S <- pars$THIN_S
  USE_NEW_FILE_BREAK <- pars$USE_NEW_FILE_BREAK
  
  ##########################################################
  # BUILD WORKFLOW
  ##########################################################
  
  heli_dt <- heli %>%
    st_drop_geometry() %>%
    mutate(
      tail = if_else(
        is.na(tail) | tail == "",
        "unknown",
        tail
      )
    ) %>%
    arrange(
      operator,
      tail,
      date_local,
      t_utc
    ) %>%
    group_by(
      operator,
      tail,
      date_local
    ) %>%
    mutate(
      gap_min =
        as.numeric(
          difftime(
            t_utc,
            lag(t_utc),
            units = "mins"
          )
        ),
      
      step_m =
        sqrt(
          (x - lag(x))^2 +
            (y - lag(y))^2
        ),
      
      new_file =
        source_file != lag(source_file),
      
      break_track =
        is.na(gap_min) |
        gap_min > MAX_GAP_MIN |
        step_m > MAX_JUMP_KM * 1000 |
        (
          USE_NEW_FILE_BREAK &
            new_file
        ),
      
      segment_id =
        cumsum(break_track)
    ) %>%
    ungroup()
  
  heli_dt <- heli_dt %>%
    group_by(
      operator,
      tail
    ) %>%
    mutate(
      dt_h =
        as.numeric(
          difftime(
            t_utc,
            lag(t_utc),
            units = "hours"
          )
        ),
      
      dist_km =
        sqrt(
          (x - lag(x))^2 +
            (y - lag(y))^2
        ) / 1000,
      
      speed_kmh =
        if_else(
          dt_h > 0,
          dist_km / dt_h,
          NA_real_
        ),
      
      airborne =
        speed_kmh >= AIRBORNE_KMH
    ) %>%
    ungroup()
  
  heli_thin <- heli_dt %>%
    group_by(
      operator,
      tail,
      date_local,
      segment_id
    ) %>%
    mutate(
      elapsed_s =
        as.numeric(
          difftime(
            t_utc,
            first(t_utc),
            units = "secs"
          )
        ),
      
      thin_bin =
        floor(
          elapsed_s / THIN_S
        )
    ) %>%
    group_by(
      operator,
      tail,
      date_local,
      segment_id,
      thin_bin
    ) %>%
    slice(1) %>%
    ungroup()
  
  heli_thin_sf <- st_as_sf(
    heli_thin,
    coords = c("x","y"),
    crs = st_crs(heli)
  )
  
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
      airborne_pts =
        sum(
          airborne,
          na.rm = TRUE
        ),
      do_union = FALSE,
      .groups = "drop"
    ) %>%
    st_cast("LINESTRING") %>%
    filter(
      n_points >= 2
    ) %>%
    mutate(
      length_km =
        as.numeric(
          st_length(geometry)
        ) / 1000
    )
  
  ##########################################################
  # SAVE ONLY SUMMARY
  ##########################################################
  
  results[[i]] <- tibble(
    run_id =
      pars$run_id,
    
    gap_min =
      MAX_GAP_MIN,
    
    jump_km =
      MAX_JUMP_KM,
    
    airborne_kmh =
      AIRBORNE_KMH,
    
    thin_s =
      THIN_S,
    
    use_new_file_break =
      USE_NEW_FILE_BREAK,
    
    n_segments =
      nrow(
        flight_lines
      ),
    
    median_length_km =
      median(
        flight_lines$length_km,
        na.rm = TRUE
      )
  )
  
  ##########################################################
  # CLEAN UP MEMORY
  ##########################################################
  
  rm(
    heli_dt,
    heli_thin,
    heli_thin_sf,
    flight_lines
  )
  
  gc(verbose = FALSE)
  
}

############################################################
# COMBINE OUTPUTS
############################################################

sensitivity_results <- bind_rows(results)

sensitivity_results

############################################################
# SAVE RESULTS
############################################################

write_csv(
  sensitivity_results,
  "outputs/sensitivity_results.csv"
)