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
# 21 - October - 2026
# Joanna Burgar
# Heli & caribou MCP overlap (spatial)
############################################################

library(sf)
library(tidyverse)
library(janitor)
library(bcdata)

#===========================================================
# INPUTS
#===========================================================

target_crs <- 3005
flight_dir <- "fl_processed"

#===========================================================
# CARIBOU
#===========================================================

caribou <- st_read(
  "//sfp.idir.bcgov/home/perpetual/HCC_flight_data/HCC_data_share_14Sept2026.shp",
  quiet = TRUE) %>%
  st_transform(target_crs) %>%
  mutate(datetime = ymd_hms(datetim, tz = "UTC"),year = year(datetime))

caribou_mcp <- caribou %>%
  filter(!is.na(year)) %>%
  group_by(COLLAR_, herd_nm, year) %>%
  summarise(.groups = "drop") %>%
  st_convex_hull() %>%
  mutate(caribou_area_km2 =as.numeric(st_area(.)) / 1e6)

st_write(
  caribou_mcp,
  "Outputs/caribou_mcp_21Oct2026.gpkg",
  delete_dsn = TRUE,
  quiet = TRUE
)

#===========================================================
# TENURES
#===========================================================

tenures <- st_read("Tenures.shp",quiet = TRUE) %>%
  clean_names() %>%
  st_transform(target_crs)

tenuresBCGOV <- bcdc_get_data("3544ad91-0cf2-4926-a08a-bfe42d9a031d") %>%
  clean_names() %>%
  filter(tenure_subpurpose == "HELI SKI") %>%
  st_transform(target_crs)

#--------------------------------------------------------
# LOOKUP TABLE
#--------------------------------------------------------

tenure_lookup <- tenures %>%
  st_drop_geometry() %>%
  transmute(crown_lands_file = cl_file,project = recode(tenure,"CMHGothics" = "Crescent_Spur"))

tenure_lookup <- tenures %>%
  st_drop_geometry() %>%
  transmute(crown_lands_file = cl_file,project = tenure)

# replace CMHGothics record with actual Crescent Spur tenure
tenure_lookup <- tenure_lookup %>%
  filter(project != "CMHGothics") %>%
  bind_rows(tibble(crown_lands_file = "7405593",project = "Crescent_Spur"))

tenure_lookup

op_tenures <- tenuresBCGOV %>%
  filter(crown_lands_file %in% tenure_lookup$crown_lands_file) %>%
  left_join(tenure_lookup,by = "crown_lands_file")

op_tenures %>% count(tenure_location)

# QA

ggplot() +
  geom_sf(data = op_tenures,aes(colour = project),fill = NA) +
  theme_bw()


#===========================================================
# FLIGHT FILES
#===========================================================

gpkg_files <- list.files(flight_dir,pattern = "_flightpoints_caribou_14Oct2026\\.gpkg$",full.names = TRUE)

#===========================================================
# BUILD TENURE-CLIPPED MCPs
#===========================================================

flight_mcp <- map_dfr(
  gpkg_files,
  function(f){
    
    cat("\nProcessing:",basename(f),"\n")
    
    x <- st_read(f,quiet = TRUE) %>%
      st_transform(target_crs)
    
    op <- unique(x$project)
    
    cat("Points before clipping:",format(nrow(x), big.mark = ","),"\n")
    
    tenure_poly <- op_tenures %>%
      filter(project == op) %>%
      summarise()
    
    if(nrow(tenure_poly) == 0){
      
      cat("No matching tenure for",op,"\n")
      
      return(NULL)
      
    }
    
    #----------------------------------
    # CLIP TO TENURE
    #----------------------------------
    
    x <- st_filter(x,tenure_poly,.predicate = st_intersects)
    
    cat("Points after clipping:",format(nrow(x), big.mark = ","),"\n")
    
    if(nrow(x) == 0){
      return(NULL)
    }
    
    x <- x %>%
      mutate(datetime = as.POSIXct(time,tz = "UTC"),year = year(datetime)) %>%
      filter(!is.na(year))
    
    x %>%
      group_by(project, year) %>%
      summarise(.groups = "drop") %>%
      st_convex_hull() %>%
      mutate(flight_area_km2 =as.numeric(st_area(.)) / 1e6)
    
  }
)

# Processing: Bearpaw_flightpoints_caribou_14Oct2026.gpkg 
# Points before clipping: 1,309,353 
# Points after clipping: 1,159,589 
# 
# Processing: Crescent_Spur_flightpoints_caribou_14Oct2026.gpkg 
# Points before clipping: 5,829,852 
# Points after clipping: 4,926,048 
# 
# Processing: MWHS_flightpoints_caribou_14Oct2026.gpkg 
# Points before clipping: 533,567 
# Points after clipping: 515,722 
# 
# Processing: Silvertip_flightpoints_caribou_14Oct2026.gpkg 
# Points before clipping: 161,205 
# Points after clipping: 61,064

#===========================================================
# SAVE MCPs
#===========================================================

st_write(flight_mcp,"Outputs/flight_mcp_tenure_clipped_21Oct2026.gpkg",delete_dsn = TRUE,quiet = TRUE)

#===========================================================
# QA MAP
#===========================================================

ggplot() +
  geom_sf(data = op_tenures,aes(colour = project),fill = NA,linewidth = 0.8) +
  geom_sf(data = flight_mcp,aes(fill = project),alpha = 0.3,colour = NA) +
  theme_bw() +
  labs(title = "Operator tenures and tenure-clipped flight MCPs")

#===========================================================
# MCP OVERLAP
#===========================================================

names(caribou_mcp); names(flight_mcp)

flight_mcp <- st_set_geometry(flight_mcp, "geom")

mcp_overlap <- purrr::map_dfr(
  intersect(unique(caribou_mcp$year),unique(flight_mcp$year)),
  function(y){
    
    cat("Year:", y, "\n")
    
    st_intersection(caribou_mcp %>%filter(year == y),
                    flight_mcp %>%filter(year == y))
    
  }
)


#===========================================================
# SAVE OVERLAP
#===========================================================

st_write(mcp_overlap,"Outputs/DONOTSHARE_mcp_overlap_21Oct2026.gpkg",delete_dsn = TRUE,quiet = TRUE)

write.csv(st_drop_geometry(mcp_overlap),"Outputs/mcp_overlap_summary.csv",row.names = FALSE)

#===========================================================
# OVERLAP MAP
#===========================================================

ggplot() +
  geom_sf(
    data = flight_mcp,
    aes(fill = project),
    alpha = 0.3,
    colour = NA
  ) +
  geom_sf(
    data = caribou_mcp,
    fill = NA,
    colour = "black",
    linewidth = 0.15
  ) +
  facet_wrap(~year) +
  theme_bw() +
  labs(
    title = "Caribou MCPs and tenure-clipped helicopter MCPs"
  )

ggsave(
  "Outputs/DONOTSHARE_mcp_overlap_map.png",
  width = 12,
  height = 8
)

#===========================================================
# SUMMARY TABLE
#===========================================================
mcp_overlap <- mcp_overlap %>%
  mutate(overlap_km2 =as.numeric(st_area(.)) / 1e6,pct_caribou_overlap =100 * overlap_km2 / caribou_area_km2,
         pct_flight_overlap =100 * overlap_km2 / flight_area_km2)

overlap_summary <- mcp_overlap %>%
  st_drop_geometry() %>%
  group_by(project, year) %>%
  summarise(n_overlap = n(),mean_overlap = mean(pct_caribou_overlap,na.rm = TRUE),
    median_overlap = median(pct_caribou_overlap,na.rm = TRUE),
    max_overlap = max(pct_caribou_overlap,na.rm = TRUE),.groups = "drop")

write.csv(overlap_summary,"Outputs/mcp_overlap_by_operator_year.csv",row.names = FALSE)

#===========================================================
# OVERLAP HISTOGRAM
#===========================================================

ggplot(
  st_drop_geometry(mcp_overlap),aes(pct_caribou_overlap)) +
  geom_histogram(bins = 30,colour = "black",fill = "steelblue") +
  theme_bw() +
  labs(x = "% of caribou MCP overlapped",y = "Count",title = "Distribution of caribou MCP overlap")

ggsave("Outputs/mcp_overlap_histogram.png",width = 8,height = 6)

#===========================================================
# OVERLAP HEATMAP
#===========================================================

# For each operator and year, 
# what was the average percentage of an individual caribou's annual MCP that overlapped the operator's annual flight MCP?

ggplot(
  overlap_summary,aes(project,factor(year),fill = mean_overlap)) +
  geom_tile() +
  geom_text(aes(label = round(mean_overlap, 1)),size = 3) +
  scale_fill_viridis_c(name = "Mean % overlap") +
  theme_bw() +
  labs(x = "Operator",y = "Year")

ggsave("Outputs/mcp_overlap_heatmap.png",width = 10,height = 6)

#===========================================================
# TOP 20 OVERLAPS
#===========================================================

mcp_overlap %>%
  st_drop_geometry() %>%
  arrange(
    desc(pct_caribou_overlap)
  ) %>%
  select(
    COLLAR_,
    herd_nm,
    project,
    year,
    overlap_km2,
    pct_caribou_overlap
  ) %>%
  head(20) %>%
  write.csv(
    "Outputs/top20_mcp_overlaps.csv",
    row.names = FALSE
  )

cat("\nFinished.\n")
