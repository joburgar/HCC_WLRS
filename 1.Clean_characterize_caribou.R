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
# Joanna Burgar (adapted from Bevan Ernst)
############################################################

# Loads and cleans GPS data
# Creates pts_clean
# Assigns herd attributes to each point
# Creates MCP95s
# Flags MCPs that overlap heli-ski tenures
# Adds the MCP-level tenure flag back to every point
# Exports:
#   pts_clean.shp
# mcp95.shp
# pts_flagged_outliers.shp
# collar_summary.csv
# two histogram PNGs
# Produces a winter collar summary
# Produces a histogram of:
#   active collars by winter
# active collars by winter and tenure overlap

############################################################
# Load libraries
############################################################

library(collar)
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

##################################################
# DOWNLOAD AND PREPARE CARIBOU GPS DATA
##################################################

key_paths <- list.files(
  "KEY_FILES",
  full.names = TRUE,
  pattern = "\\.keyx$",
  recursive = TRUE
)

gps <- purrr::map_dfr(
  key_paths[1:3],
  ~fetch_vectronics(
    key_paths = .x,
    start_date = "2019-10-01T00:00:00",
    which_date = "acquisition"
  ),
  .progress = TRUE
)

gps %>%
  count(idcollar, sort = TRUE)


##################################################
# CLEAN
##################################################
names(gps)
gps <- gps %>%
  filter(
    !is.na(latitude),
    !is.na(longitude)
  ) %>%
  mutate(
    datetime =
      ymd_hms(acquisitiontime, tz = "UTC") -
      hours(8),
    COLLAR_ID = as.numeric(idcollar)
  ) %>%
  select(
    COLLAR_ID,
    datetime,
    latitude,
    longitude,
    dop,
    idmortalitystatus,
    activity,
    temperature,
    mainvoltage,
    backupvoltage
  )


names(gps)

##################################################
# SAVE EXPORTED DATA
##################################################

saveRDS(gps, "gps_download.rds") # save
gps <- readRDS("gps_download.rds") # load

message(
  "Downloaded ",
  nrow(gps),
  " fixes from ",
  n_distinct(gps$COLLAR_ID),
  " collars."
)


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
    !is.na(longitude),
    !is.na(datetime),
    dop <= 8
  )

############################################################
# CREATE SF OBJECT
############################################################

pts <- st_as_sf(
  gps_clean,
  coords = c("longitude", "latitude"),
  crs = 4326
)

pts <- st_transform(pts, 3005)
# dim(pts) # [1] 105838     33


############################################################
# REMOVE POINTS FAR FROM MAIN CLUSTER
############################################################

coords <- st_coordinates(pts)

pts$X <- coords[, 1]
pts$Y <- coords[, 2]

centers <- pts %>%
  st_drop_geometry() %>%
  group_by(COLLAR_ID) %>%
  summarise(
    med_x = median(X),
    med_y = median(Y),
    .groups = "drop"
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
# EXPORT OUTLIER POINTS
############################################################

st_write(
  outlier_pts,
  file.path(output_dir,
            "pts_flagged_outliers.shp"),
  driver = "ESRI Shapefile",
  delete_layer = TRUE
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
      round(100 * flagged / total_fixes, 1),
    .groups = "drop"
  )

write_csv(
  outlier_summary,
  file.path(
    output_dir,
    "outlier_summary.csv"
  )
)

############################################################
# ASSIGN HERD ATTRIBUTES
############################################################

herds_small <- herds %>%
  select(
    herd_name,
    ecotype
  )

pts_clean <- st_join(
  pts_clean,
  herds_small,
  join = st_intersects,
  left = TRUE
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
# MCP-TENURE OVERLAP
############################################################

mcp95_sf <- mcp95_sf %>%
  mutate(
    overlaps_tenure =
      lengths(
        st_intersects(
          .,
          tenures
        )
      ) > 0
  )

############################################################
# ADD TENURE FLAG TO POINTS
############################################################

collar_overlap <- mcp95_sf %>%
  st_drop_geometry() %>%
  transmute(
    COLLAR_ID = as.numeric(id),
    overlaps_tenure
  )

pts_clean <- pts_clean %>%
  left_join(
    collar_overlap,
    by = "COLLAR_ID"
  )

############################################################
# EXPORT SPATIAL LAYERS
############################################################

st_write(
  pts_clean,
  file.path(
    output_dir,
    "pts_clean.shp"),
  driver = "ESRI Shapefile",
  delete_layer = TRUE
)

st_write(
  mcp95_sf,
  file.path(
    output_dir,
    "mcp95.shp"),
  driver = "ESRI Shapefile",
  delete_layer = TRUE
)

st_write(
  mcp100_sf,
  file.path(output_dir,
            "mcp100.shp"),
  driver = "ESRI Shapefile",
  delete_layer = TRUE
)


############################################################
# WINTER SUMMARY
############################################################

winter_pts <- pts_clean %>%
  st_drop_geometry() %>%
  mutate(
    month = month(datetime),
    winter = case_when(
      month %in% c(11, 12) ~ year(datetime) + 1,
      month %in% c(1, 2, 3, 4) ~ year(datetime),
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(winter))

collar_summary <- winter_pts %>%
  group_by(
    COLLAR_ID,
    herd_name,
    overlaps_tenure,
    winter
  ) %>%
  summarise(
    fixes = n(),
    start_date = min(datetime),
    end_date = max(datetime),
    mean_dop = mean(dop, na.rm = TRUE),
    max_dop = max(dop, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  collar_summary,
  file.path(
    output_dir,
    "collar_summary.csv"
  )
)

############################################################
# ACTIVE COLLARS BY WINTER
############################################################

active_collars <- collar_summary %>%
  filter(fixes > 30) %>%
  group_by(winter) %>%
  summarise(
    n_collars = n_distinct(COLLAR_ID),
    .groups = "drop"
  )

p1 <- ggplot(
  active_collars,
  aes(
    x = factor(winter),
    y = n_collars)) +
  geom_col() +
  geom_text(
    aes(label = n_collars),
    vjust = -0.3) +
  labs(
    x = "Winter",
    y = "Number of Active Collars",
    title = "Active Collars by Winter (>30 fixes)") +
  theme_minimal()

ggsave(
  file.path(
    output_dir,
    "active_collars_by_winter.png"),
  p1,
  width = 8,
  height = 5,
  dpi = 300)

############################################################
# ACTIVE COLLARS BY WINTER AND TENURE OVERLAP
############################################################

active_collars_tenure <- collar_summary %>%
  filter(fixes > 30) %>%
  group_by(
    winter,
    overlaps_tenure) %>%
  summarise(
    n_collars = n_distinct(COLLAR_ID),
    .groups = "drop")

p2 <- ggplot(
  active_collars_tenure,
  aes(
    x = factor(winter),
    y = n_collars,
    fill = overlaps_tenure)) +
  geom_col() +
  labs(
    x = "Winter",
    y = "Number of Active Collars",
    fill = "MCP Overlaps Tenure",
    title = "Active Collars by Winter and Tenure Overlap") +
  theme_minimal()

ggsave(
  file.path(
    output_dir,
    "active_collars_by_winter_tenure.png"
  ),
  p2,
  width = 8,
  height = 5,
  dpi = 300
)


herd_overlap_summary <- collar_summary %>%
  mutate(
    overlaps_tenure = replace_na(overlaps_tenure, FALSE)) %>%
  group_by(herd_name, COLLAR_ID) %>%
  summarise(
    overlaps_tenure = any(overlaps_tenure),
    .groups = "drop") %>%
  group_by(herd_name) %>%
  summarise(
    n_collars = n(),
    n_overlap = sum(overlaps_tenure),
    pct_overlap = round(
      100 * n_overlap / n_collars,
      1),
    .groups = "drop")

write_csv(
  herd_overlap_summary,
  file.path(
    output_dir,
    "herd_overlap_summary.csv"))

############################################################
# HERD MONITORING HEATMAP
############################################################

herd_winter <- collar_summary %>%
  filter(
    fixes > 30,
    !is.na(herd_name)
  ) %>%
  group_by(
    herd_name,
    winter
  ) %>%
  summarise(
    n_collars = n_distinct(COLLAR_ID),
    .groups = "drop"
  )

############################################################
# ORDER HERDS BY TOTAL MONITORING EFFORT
############################################################

herd_order <- herd_winter %>%
  group_by(herd_name) %>%
  summarise(
    total_collars = sum(n_collars),
    .groups = "drop"
  ) %>%
  arrange(desc(total_collars)) %>%
  pull(herd_name)

herd_winter$herd_name <- factor(
  herd_winter$herd_name,
  levels = rev(herd_order)
)

############################################################
# NON-TENURE MONITORING HEATMAP
############################################################

herd_winter_no_tenure <- collar_summary %>%
  filter(
    fixes > 30,
    !overlaps_tenure,
    !is.na(herd_name)
  ) %>%
  group_by(
    herd_name,
    winter
  ) %>%
  summarise(
    n_collars = n_distinct(COLLAR_ID),
    .groups = "drop"
  )

herd_winter_no_tenure$herd_name <- factor(
  herd_winter_no_tenure$herd_name,
  levels = levels(herd_winter$herd_name)
)

p_heatmap_no_tenure <- ggplot(
  herd_winter_no_tenure,
  aes(
    x = factor(winter),
    y = herd_name,
    fill = n_collars)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_collars),color = "white",size = 3) +
  scale_fill_gradient(low = "#deebf7",high = "#08519c", name = "Active\ncollars") +
  labs(
    title = "Monitoring Coverage for Herds Without Heli-Ski Tenure Overlap",
    subtitle = "Collared individuals whose MCP does not overlap a heli-ski tenure",
    x = "Winter",
    y = "Herd"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8))

ggsave(
  file.path(
    output_dir,
    "herd_monitoring_heatmap_no_tenure_overlap.png"
  ),
  p_heatmap_no_tenure,
  width = 10,
  height = 7,
  dpi = 300
)
############################################################
# TENURE-OVERLAP MONITORING HEATMAP
############################################################

herd_winter_tenure <- collar_summary %>%
  filter(
    fixes > 30,
    overlaps_tenure,
    !is.na(herd_name)
  ) %>%
  group_by(
    herd_name,
    winter
  ) %>%
  summarise(
    n_collars = n_distinct(COLLAR_ID),
    .groups = "drop"
  )

herd_winter_tenure$herd_name <- factor(
  herd_winter_tenure$herd_name,
  levels = levels(herd_winter$herd_name)
)

p_heatmap_tenure <- ggplot(
  herd_winter_tenure,
  aes(
    x = factor(winter),
    y = herd_name,
    fill = n_collars
  )
) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_collars), color = "white", size = 3)+
  scale_fill_gradient(low = "#fee5d9",high = "#a50f15",name = "Active\ncollars") +
  labs(
    title = "Monitoring Coverage for Herds with Heli-Ski Tenure Overlap",
    subtitle = "Collared individuals whose MCP overlaps at least one heli-ski tenure",
    x = "Winter",
    y = "Herd"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  file.path(
    output_dir,
    "herd_monitoring_heatmap_tenure_overlap.png"
  ),
  p_heatmap_tenure,
  width = 10,
  height = 7,
  dpi = 300
)

############################################################
# QA PLOT
############################################################

plot(
  st_geometry(mcp95_sf),
  border = ifelse(
    mcp95_sf$overlaps_tenure,
    "red",
    "blue"
  )
)

plot(
  st_geometry(pts_clean),
  pch = 16,
  cex = 0.2,
  add = TRUE
)

############################################################
# DONE
############################################################
