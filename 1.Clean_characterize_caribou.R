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
# Caribou Winter Collar Characterization
# Joanna Burgar
############################################################

library(sf)
library(tidyverse)
library(adehabitatHR)
library(lubridate)
library(bcdata)
library(janitor)

############################################################
# OUTPUT FOLDER
############################################################

output_dir <- "Outputs"
dir.create(output_dir, showWarnings = FALSE)

############################################################
# LOAD GPS DATA
############################################################

gps <- read_csv("Caribou_winter_MCP.csv")

############################################################
# LOAD HERDS
############################################################

herds <- bcdc_get_data(
  "2b217585-f48d-4d9f-b7ba-746909ac35ca"
) |>
  clean_names() |>
  filter(ecotype == "Mountain")

############################################################
# LOAD HELI-SKI TENURES
############################################################

tenures <- bcdc_get_data(
  "3544ad91-0cf2-4926-a08a-bfe42d9a031d"
) |>
  clean_names() |>
  filter(tenure_subpurpose == "HELI SKI")

############################################################
# CLEAN GPS DATA
############################################################

gps_clean <- gps %>%
  filter(
    !is.na(latitude),
    !is.na(longitude)
  ) %>%
  filter(
    dop <= 8
  ) %>%
  mutate(
    datetime = mdy_hms(Date_Time_)
  ) %>%
  filter(!is.na(datetime))

############################################################
# CREATE SF OBJECT
############################################################

pts <- st_as_sf(
  gps_clean,
  coords = c("longitude", "latitude"),
  crs = 4326
)

pts <- st_transform(pts, 3005)

############################################################
# COLLAR SUMMARY
############################################################

collar_summary <- pts %>%
  st_drop_geometry() %>%
  group_by(COLLAR_ID) %>%
  summarise(
    fixes = n(),
    start_date = min(datetime),
    end_date = max(datetime),
    mean_dop = mean(dop, na.rm = TRUE),
    max_dop = max(dop, na.rm = TRUE)
  )

write_csv(
  collar_summary,
  file.path(output_dir,
            "collar_summary.csv")
)

############################################################
# REMOVE POINTS FAR FROM MAIN CLUSTER
############################################################

coords <- st_coordinates(pts)

pts$X <- coords[,1]
pts$Y <- coords[,2]

centers <- pts %>%
  st_drop_geometry() %>%
  group_by(COLLAR_ID) %>%
  summarise(
    med_x = median(X),
    med_y = median(Y)
  )

pts <- left_join(
  pts,
  centers,
  by = "COLLAR_ID"
)

pts$dist_from_center_m <- sqrt(
  (pts$X - pts$med_x)^2 +
    (pts$Y - pts$med_y)^2
)

############################################################
# CHANGE THIS THRESHOLD IF NEEDED
############################################################

distance_threshold_km <- 100

pts$distance_outlier <-
  pts$dist_from_center_m >
  distance_threshold_km * 1000

############################################################
# OUTLIERS
############################################################

outlier_pts <- pts %>%
  filter(distance_outlier)

pts_clean <- pts %>%
  filter(!distance_outlier)

############################################################
# EXPORT POINTS
############################################################

st_write(
  pts_clean,
  file.path(output_dir,
            "pts_clean.gpkg"),
  delete_dsn = TRUE
)

st_write(
  outlier_pts,
  file.path(output_dir,
            "pts_flagged_outliers.gpkg"),
  delete_dsn = TRUE
)

############################################################
# OUTLIER SUMMARY
############################################################

outlier_summary <- pts %>%
  st_drop_geometry() %>%
  group_by(COLLAR_ID) %>%
  summarise(
    total_fixes = n(),
    flagged = sum(distance_outlier),
    percent_flagged =
      round(100 * flagged / total_fixes,1)
  )

write_csv(
  outlier_summary,
  file.path(output_dir,
            "outlier_summary.csv")
)

############################################################
# MCP GENERATION
############################################################

valid_ids <- pts_clean %>%
  st_drop_geometry() %>%
  count(COLLAR_ID) %>%
  filter(n >= 30) %>%
  pull(COLLAR_ID)

mcp_pts <- pts_clean %>%
  filter(COLLAR_ID %in% valid_ids)

sp_pts <- as(mcp_pts, "Spatial")

############################################################
# MCP95
############################################################

mcp95 <- mcp(
  sp_pts["COLLAR_ID"],
  percent = 95
)

mcp95_sf <- st_as_sf(mcp95)

############################################################
# MCP100
############################################################

mcp100 <- mcp(
  sp_pts["COLLAR_ID"],
  percent = 100
)

mcp100_sf <- st_as_sf(mcp100)

############################################################
# MCP AREA
############################################################

mcp95_sf$MCP95_km2 <-
  as.numeric(st_area(mcp95_sf))/1000000

mcp100_sf$MCP100_km2 <-
  as.numeric(st_area(mcp100_sf))/1000000

############################################################
# EXPORT MCPS
############################################################

st_write(
  mcp95_sf,
  file.path(output_dir,
            "mcp95.gpkg"),
  delete_dsn = TRUE
)

st_write(
  mcp100_sf,
  file.path(output_dir,
            "mcp100.gpkg"),
  delete_dsn = TRUE
)

############################################################
# ASSIGN HERD BASED ON POINTS
############################################################

herds_small <- herds %>%
  select(
    herd_name = herd_name,
    ecotype
  )

pts_herd <- st_join(
  pts_clean,
  herds_small,
  join = st_intersects,
  left = TRUE
)

herd_assignment <- pts_herd %>%
  st_drop_geometry() %>%
  count(COLLAR_ID,
        herd_name = herd_name,
        sort = TRUE) %>%
  group_by(COLLAR_ID) %>%
  slice_max(n, n = 1) %>%
  ungroup()

write_csv(
  herd_assignment,
  file.path(output_dir,
            "collar_herd_assignment.csv")
)

############################################################
# HELI-SKI TENURE OVERLAP
############################################################

tenures_small <- tenures %>%
  select(
    tenure_location,
    client_name
  )

pts_tenure <- st_join(
  pts_clean,
  tenures_small,
  join = st_intersects,
  left = TRUE
)

############################################################
# INDIVIDUAL OVERLAP
############################################################

individual_overlap <- pts_tenure %>%
  st_drop_geometry() %>%
  group_by(COLLAR_ID) %>%
  summarise(
    overlap_tenure =
      any(!is.na(tenure_location))
  ) %>%
  left_join(
    herd_assignment,
    by = "COLLAR_ID"
  )

write_csv(
  individual_overlap,
  file.path(output_dir,
            "individual_overlap.csv")
)

############################################################
# HERD SUMMARY
############################################################

herd_summary <- individual_overlap %>%
  group_by(herd_name) %>%
  summarise(
    individuals = n(),
    overlap = sum(overlap_tenure),
    no_overlap = sum(!overlap_tenure)
  )

write_csv(
  herd_summary,
  file.path(output_dir,
            "herd_summary.csv")
)

############################################################
# TENURE SUMMARY
############################################################

tenure_summary <- pts_tenure %>%
  st_drop_geometry() %>%
  filter(!is.na(tenure_location)) %>%
  distinct(COLLAR_ID,
           tenure_location) %>%
  count(tenure_location)

write_csv(
  tenure_summary,
  file.path(output_dir,
            "tenure_summary.csv")
)

############################################################
# EXPORT SPATIAL LAYERS
############################################################

st_write(
  pts_herd,
  file.path(output_dir,
            "pts_clean_herds.gpkg"),
  delete_dsn = TRUE
)

st_write(
  pts_tenure,
  file.path(output_dir,
            "pts_clean_tenures.gpkg"),
  delete_dsn = TRUE
)

############################################################
# SIMPLE QA PLOT
############################################################

plot(st_geometry(mcp100_sf),
     border = "red")

plot(st_geometry(mcp95_sf),
     border = "blue",
     add = TRUE)

plot(st_geometry(pts_clean),
     pch = 16,
     cex = 0.3,
     add = TRUE)

############################################################
# DONE
############################################################
