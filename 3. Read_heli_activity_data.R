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

#############################################################
# Helicopter Activity Data (read and collate)
# 5 - October - 2026
# Joanna Burgar
# Helicopter Activity Data
# Rebuild by project, preserve timestamps,
# clip to mountain caribou herds
############################################################

library(sf)
library(tidyverse)
library(bcdata)
library(janitor)

#===========================================================
# DIRECTORIES
#===========================================================

flight_dir <- "./raw_data"

flight_dir_out <- "fl_processed"

#===========================================================
# CARIBOU HERDS
#===========================================================

herds <- bcdc_get_data(
  "2b217585-f48d-4d9f-b7ba-746909ac35ca"
) |>
  clean_names() |>
  filter(ecotype == "Mountain")

herds_small <- herds |>
  select(
    herd_name,
    ecotype
  )

# 20 km buffer
caribou_buff <- st_buffer(
  herds_small,
  20000
)

#===========================================================
# FIND INPUT FILES
#===========================================================

gpkg_files <- list.files(
  flight_dir,
  pattern = "\\.gpkg$",
  full.names = TRUE
)

#===========================================================
# IDENTIFY PROJECT FOR EACH FILE
#===========================================================

project_lookup <- tibble(
  file = gpkg_files,
  fname = tools::file_path_sans_ext(
    basename(gpkg_files)
  )
) |>
  mutate(
    project = str_remove(
      fname,
      "_flightlines_(GPX|kml)_\\d{4}$"
    )
  )

projects <- sort(
  unique(project_lookup$project)
)

print(projects)

#===========================================================
# PROCESS PROJECTS INDIVIDUALLY
#===========================================================

for(p in projects){
  
  cat("\n")
  cat("=====================================\n")
  cat("PROCESSING:", p, "\n")
  cat("=====================================\n")
  
  project_files <- project_lookup |>
    filter(project == p) |>
    pull(file)
  
  project_dat <- purrr::map_dfr(
    project_files,
    function(f){
      
      message("Reading: ", basename(f))
      
      dat <- st_read(
        f,
        quiet = TRUE
      )
      
      fname <- tools::file_path_sans_ext(
        basename(f)
      )
      
      #----------------------------------
      # STANDARDIZE COLUMNS
      #----------------------------------
      
      if(!"time" %in% names(dat))
        dat$time <- as.POSIXct(NA)
      
      if(!"tail" %in% names(dat))
        dat$tail <- NA_character_
      
      if(!"alt_m" %in% names(dat))
        dat$alt_m <- NA_real_
      
      if(!"speed_kn" %in% names(dat))
        dat$speed_kn <- NA_real_
      
      if(!"heading" %in% names(dat))
        dat$heading <- NA_real_
      
      if(!"source_file" %in% names(dat))
        dat$source_file <- basename(f)
      
      if(!"operator" %in% names(dat))
        dat$operator <- NA_character_
      
      if(!"flight_date" %in% names(dat))
        dat$flight_date <- as.Date(NA)
      
      dat |>
        mutate(
          file_name = fname,
          
          source = case_when(
            str_detect(fname, "_GPX_") ~ "GPX",
            str_detect(fname, "_kml_") ~ "KML",
            TRUE ~ NA_character_
          ),
          
          year = str_extract(
            fname,
            "\\d{4}$"
          ),
          
          project = str_remove(
            fname,
            "_flightlines_(GPX|kml)_\\d{4}$"
          )
        ) |>
        
        select(
          any_of(c(
            "project",
            "source",
            "source_file",
            "year",
            "operator",
            "flight_date",
            "time",
            "tail",
            "alt_m",
            "speed_kn",
            "heading"
          ))
        )
    }
  )
  
  cat(
    "\nOriginal points:",
    format(nrow(project_dat), big.mark = ","),
    "\n"
  )
  
  #---------------------------------------------------
  # MATCH CRS
  #---------------------------------------------------
  
  project_dat <- st_transform(
    project_dat,
    st_crs(caribou_buff)
  )
  
  #---------------------------------------------------
  # CLIP TO BUFFERED HERDS
  #---------------------------------------------------
  
  project_dat <- st_filter(
    project_dat,
    caribou_buff,
    .predicate = st_intersects
  )
  
  cat(
    "Clipped points:",
    format(nrow(project_dat), big.mark = ","),
    "\n"
  )
  
  #---------------------------------------------------
  # QA CHECK
  #---------------------------------------------------
  
  cat(
    "Non-missing timestamps:",
    sum(!is.na(project_dat$time)),
    "\n"
  )
  
  print(
    range(
      project_dat$time,
      na.rm = TRUE
    )
  )
  
  #---------------------------------------------------
  # SAVE AS GPKG
  #---------------------------------------------------
  
  out_file <- file.path(
    flight_dir_out,
    paste0(
      p,
      "_flightpoints_caribou_14Oct2026.gpkg"
    )
  )
  
  message(
    "Writing: ",
    basename(out_file)
  )
  
  st_write(
    project_dat,
    out_file,
    delete_dsn = TRUE,
    quiet = TRUE
  )
  
  #---------------------------------------------------
  # CLEAN MEMORY
  #---------------------------------------------------
  
  rm(project_dat)
  
  gc(full = TRUE)
}

cat("\nFinished.\n")
