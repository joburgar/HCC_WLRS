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
############################################################

# ================================
# COMBINE ALL FLIGHTLINE GPKGs
# ================================

library(sf)
library(tidyverse)
library(mapview)

# Folder containing the gpkg files
flight_dir <- "//sfp.idir.bcgov/home/perpetual/HCC_flight_data"

# Find all geopackages
gpkg_files <- list.files(
  flight_dir,
  pattern = "\\.gpkg$",
  full.names = TRUE
)

# Read and combine
flightlines <- map_dfr(gpkg_files, function(f) {
  
  message("Reading: ", basename(f))
  
  dat <- st_read(f, quiet = TRUE)
  
  fname <- tools::file_path_sans_ext(basename(f))
  
  dat %>%
    mutate(
      file_name = fname,
      
      # Source type
      source = case_when(
        str_detect(fname, "_GPX_") ~ "GPX",
        str_detect(fname, "_kml_") ~ "KML",
        TRUE ~ NA_character_
      ),
      
      # Year
      year = str_extract(fname, "\\d{4}$"),
      
      # Project name
      project = str_remove(
        fname,
        "_flightlines_(GPX|kml)_\\d{4}$"
      )
    )
})

# Check results
print(table(flightlines$project))
print(table(flightlines$year))
print(st_geometry_type(flightlines) |> table())


# Save combined geopackage
flightlines_out <- flightlines %>% select(
  project,
  source_file,
  year,
  operator,
  flight_date,
  source_file,
  tail,
  event,
  alt_m,
  speed_kn,
  heading,
  geom)

st_write(
  flightlines_out,
  file.path(flight_dir, "flightlines_5Oct2026.gpkg"),
  layer = "flightlines",
  delete_dsn = TRUE
)


# Optional: save shapefile
st_write(
  flightlines_out,
  file.path(flight_dir, "flightlines_5Oct2026.shp"),
  delete_layer = TRUE
)

# # Quick interactive map
# mapview(
#   flightlines_out,
#   zcol = "project"
# )

################################################################################
## Clipped heli data to caribou herds
library(sf)
library(tidyverse)
library(bcdata)
library(janitor)

# Read flightlines
flightlines <- st_read(
  "//sfp.idir.bcgov/home/perpetual/HCC_flight_data/flightlines_5Oct2026.gpkg",
  quiet = TRUE
)

# Read BC caribou herd polygons
herds <- bcdc_get_data(
  "2b217585-f48d-4d9f-b7ba-746909ac35ca"
) |>
  clean_names() |>
  filter(ecotype == "Mountain")

herds_small <- herds %>%
  select(
    herd_name,
    ecotype
  )

# Match projections
flightlines <- st_transform(flightlines, st_crs(herds_small))

# Buffer herds by 20 km to capture nearby flight tracks
caribou_buff <- st_buffer(herds_small, 20000)

# Keep only points intersecting buffered herd ranges
flightlines_clip <- st_filter(
  flightlines,
  caribou_buff,
  .predicate = st_intersects
)

# Save clipped dataset
st_write(
  flightlines_clip,
  "//sfp.idir.bcgov/home/perpetual/HCC_flight_data/combined_flightlines_caribou.gpkg",
  delete_dsn = TRUE
)

# Quick summary
cat("Original points:", nrow(flightlines), "\n")
cat("Clipped points:", nrow(flightlines_clip), "\n")
