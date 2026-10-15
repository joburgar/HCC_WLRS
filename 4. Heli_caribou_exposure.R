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
# Exposure matrix (helicopter and caribou)
# 14 - October - 2026
# Joanna Burgar
############################################################

# ==========================================================
# 1. How many potential helicopter encounters exist
# 2. What is the distance-to-nearest flight distribution
# 3. Does movement change following exposure
# 4. What threshold maximizes signal
# ==========================================================

### LOAD LIBRARIES

library(sf)
library(tidyverse)

### LOAD SF OBJECTS
flight_dir <- "//sfp.idir.bcgov/home/perpetual/HCC_flight_data"
flight_dir_out <- "fl_processed"


# flights <- purrr::map_dfr(list.files(flight_dir_out, pattern = "_flightpoints_caribou_14Oct2026\\.gpkg$", full.names = TRUE),
#                           ~st_read(.x, quiet = TRUE))
# 
# flights %>% count(project)
# range(flights$time, na.rm = TRUE)
# st_write(flights, "combined_flightpoints_caribou_14Oct2026.gpkg", delete_dsn = TRUE)

# if already saved as single objects
flights <- st_read(file.path(flight_dir,"combined_flightpoints_caribou_14Oct2026.gpkg"))
caribou <- st_read("//sfp.idir.bcgov/home/perpetual/HCC_flight_data/HCC_data_share_14Sept2026.shp")
print(caribou)


cat("\n====================\n")
cat("FLIGHT DATA\n")
cat("====================\n")

print(flights)

cat("\nFields:\n")
print(names(flights))

cat("\nCRS:\n")
print(st_crs(flights))

cat("\nGeometry types:\n")
print(table(st_geometry_type(flights)))

#--------------------------------------------------
# REPROJECT TO BC ALBERS
#--------------------------------------------------

target_crs <- 3005

caribou <- st_transform(caribou, target_crs)
flights <- st_transform(flights, target_crs)

#--------------------------------------------------
# SPATIAL EXTENTS
#--------------------------------------------------

cat("\nCARIBOU BBOX\n")
print(st_bbox(caribou))

cat("\nFLIGHTS BBOX\n")
print(st_bbox(flights))

#--------------------------------------------------
# QUICK MAP
#--------------------------------------------------

plot(st_geometry(flights),
     col = "grey",
     main = "Flightlines and Caribou")

plot(st_geometry(caribou),
     pch = 16,
     col = "red",
     add = TRUE)

#--------------------------------------------------
# LOOK FOR DATE/TIME FIELDS
#--------------------------------------------------

cat("\n\nPotential time fields in CARIBOU:\n")
print(names(caribou)[grepl(
  "date|time|datetime|timestamp",
  names(caribou),
  ignore.case = TRUE
)])
# [1] "datetim" # need to check if PST or GMT
print(caribou)

cat("\nPotential time fields in FLIGHTS:\n")
print(names(flights)[grepl(
  "date|time|datetime|timestamp",
  names(flights),
  ignore.case = TRUE
)])
# [1] "flight_date" "time" # time is the better option
print(flights)

#----------------------------------------------------------
# SETTINGS
#----------------------------------------------------------

dist_thresholds <- c(100,500,1000,2000,5000)
time_thresholds <- c(2,6,12,24)

#----------------------------------------------------------
# FUNCTION TO FIND NEAREST HELICOPTER OBSERVATION
#----------------------------------------------------------

class(caribou$datetim)
class(flights$time)

str(caribou$datetim)
str(flights$time)

caribou <- caribou %>% mutate(datetime = ymd_hms(datetim, tz = "UTC"))
caribou <- caribou %>% mutate(year = year(datetime))
flights <- flights %>% mutate(datetime = as.POSIXct(time, tz = "UTC"))
flights <- flights %>% mutate(year = year(datetime))

caribou %>%st_drop_geometry() %>% count(year = year(datetime))
flights %>%st_drop_geometry() %>% count(year = year(datetime))

caribou %>%
  st_drop_geometry() %>%
  ggplot(aes(year(datetime))) +
  geom_bar() +
  labs(
    title = "Caribou fixes by year",
    x = "Year",
    y = "Number of fixes"
  )

flights %>%
  st_drop_geometry() %>%
  ggplot(aes(year(datetime))) +
  geom_bar() +
  labs(
    title = "Flight GPS points by year",
    x = "Year",
    y = "Number of points"
  )


caribou %>%
  st_drop_geometry() %>%
  mutate(
    datetime_bc = with_tz(datetime, "America/Vancouver")
  ) %>%
  count(hour = hour(datetime_bc)) %>%
  ggplot(aes(hour, n)) +
  geom_col() +
  scale_x_continuous(breaks = 0:23) +
  labs(
    title = "Caribou fixes by hour (BC local time)",
    x = "Hour",
    y = "Number of fixes"
  )

flights %>%
  st_drop_geometry() %>%
  mutate(
    datetime_bc = with_tz(datetime, "America/Vancouver")
  ) %>%
  count(hour = hour(datetime_bc)) %>%
  ggplot(aes(hour, n)) +
  geom_col() +
  scale_x_continuous(breaks = 0:23) +
  labs(
    title = "Flight points by hour (BC local time)",
    x = "Hour",
    y = "Number of flight points"
  )

caribou %>%
  st_drop_geometry() %>%
  count(year = year(datetime),
        month = month(datetime, label = TRUE)) %>%
  ggplot(aes(month, n)) +
  geom_col() +
  facet_wrap(~year) +
  labs(
    title = "Caribou fixes by month"
  )

flights %>%
  st_drop_geometry() %>%
  count(year = year(datetime),
        month = month(datetime, label = TRUE)) %>%
  ggplot(aes(month, n)) +
  geom_col() +
  facet_wrap(~year) +
  labs(
    title = "Flight points by month"
  )


nearest_flight <- function(i){
  
  cb <- caribou[i,]
  
  tdiff <- abs(as.numeric(
    difftime(
      flights$datetime,
      cb$datetime,
      units = "hours"
    )
  ))
  
  nearest_time_idx <- which.min(tdiff)
  
  nearest_dist <- as.numeric(
    st_distance(
      cb,
      flights[nearest_time_idx,],
      by_element = TRUE
    )
  )
  
  tibble(
    animal_id = cb$animal_id,
    herd = cb$herd,
    datetime = cb$datetime,
    nearest_hours = tdiff[nearest_time_idx],
    nearest_distance = nearest_dist
  )
}

#----------------------------------------------------------
# CALCULATE NEAREST DISTANCE/TIME FOR EACH CARIBOU FIX
#----------------------------------------------------------

years <- intersect(
  unique(year(caribou$datetime)),
  unique(year(flights$datetime))
)
years <- sort(years)

for(y in years){
  
  cat("\nYear:", y, "\n")
  
  caribou_y <- filter(caribou, year == y)
  
  for(op in unique(flights$project)){
    
    flights_sub <- flights %>%
      filter(
        year == y,
        project == op
      )
    
    cat(
      op,
      "- caribou fixes:",
      nrow(caribou_y),
      "- flight points:",
      nrow(flights_sub),
      "\n"
    )
  }
}

flight_summary_report <- flight_summary %>%
  mutate(
    n_points = format(n_points, big.mark = ","),
    n_days = format(n_days, big.mark = ","),
    start = format(start, "%Y-%m-%d"),
    end = format(end, "%Y-%m-%d")
  )

write.csv(
  flight_summary_report,
  "flight_data_summary_report.csv",
  row.names = FALSE
)


flight_days <- flights %>%
  st_drop_geometry() %>%
  mutate(date = as.Date(datetime)) %>%
  distinct(project, date)


flight_resolution <- flights %>%
  st_drop_geometry() %>%
  filter(!is.na(datetime)) %>%
  arrange(project, source_file, datetime) %>%
  group_by(project, source_file) %>%
  mutate(
    dt_sec = as.numeric(
      difftime(
        datetime,
        lag(datetime),
        units = "secs"
      )
    )
  ) %>%
  ungroup()

flight_resolution_report <- flight_resolution %>%
  group_by(project) %>%
  summarise(
    n = sum(!is.na(dt_sec)),
    min_sec = min(dt_sec, na.rm = TRUE),
    p25_sec = quantile(dt_sec, 0.25, na.rm = TRUE),
    median_sec = median(dt_sec, na.rm = TRUE),
    mean_sec = mean(dt_sec, na.rm = TRUE),
    p75_sec = quantile(dt_sec, 0.75, na.rm = TRUE),
    max_sec = max(dt_sec, na.rm = TRUE))

write.csv(
  flight_resolution_report,
  "flight_resolution_report.csv",
  row.names = FALSE
)

flight_resolution %>%
  filter(dt_sec < 150) %>%   # ignore huge gaps between flights
  ggplot(aes(dt_sec)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~project, scales = "free_y")

flight_resolution %>%
  filter(dt_sec == 0) %>%
  count(project)

# loads of duplicates in MWHS flight data, also very short intervals between points

flights %>%
  st_drop_geometry() %>%
  filter(project == "MWHS") %>%
  arrange(source_file, datetime) %>%
  select(source_file, datetime) %>%
  head(50)

flights_with_time <- flights %>%
  st_drop_geometry() %>%
  group_by(project) %>%
  summarise(
    n_points = n(),
    n_with_time = sum(!is.na(datetime)),
    pct_with_time = round(
      100 * n_with_time / n_points,1))

write.csv(
  flights_with_time,
  "flights_with_time.csv",
  row.names = FALSE
)

flights %>%
  st_drop_geometry() %>%
  mutate(date = as.Date(datetime)) %>%
  count(project, date)

###################################################################

#=================================================
# OPERATOR
#=================================================
op <- "Crescent_Spur"

#=================================================
# THIN TO 60 SECONDS
#=================================================

flights_thin <- flights %>%
  filter(
    project == op,
    !is.na(datetime)
  ) %>%
  arrange(source_file, datetime) %>%
  mutate(
    minute_bin = floor_date(datetime, "1 minute")
  ) %>%
  group_by(source_file, minute_bin) %>%
  slice_head(n = 1) %>%
  ungroup()

flights %>%
  filter(project == op, !is.na(datetime)) %>%
  summarise(n_original = n())

flights_thin %>%
  summarise(n_thin = n())

# silvertip: 161205 to 105091
# Bearpaw: 574933 to 32453
# MWHS: 373795 to 80631
# Crescent_Spur: 5804297 to 97788
#=================================================
# YEARS
#=================================================

yrs <- sort(
  intersect(
    year(caribou$datetime),
    year(flights_thin$datetime)
  )
)

#=================================================
# PROCESS EACH YEAR
#=================================================

all_results <- map_dfr(yrs, function(yr){
  
  cat("Processing", yr, "\n")
  
  caribou_sub <- caribou %>%
    filter(year(datetime) == yr) %>%
    mutate(date = as.Date(datetime))
  
  flights_sub <- flights_thin %>%
    filter(year(datetime) == yr) %>%
    mutate(date = as.Date(datetime))
  
  dates <- intersect(
    unique(caribou_sub$date),
    unique(flights_sub$date)
  )
  
  if(length(dates) == 0) return(NULL)
  
  out <- map_dfr(dates, function(d){
    
    c_day <- filter(caribou_sub, date == d)
    
    f_day <- filter(flights_sub, date == d)
    
    idx <- st_nearest_feature(
      c_day,
      f_day
    )
    
    tibble(
      year = yr,
      nearest_dist_m = as.numeric(
        st_distance(
          c_day,
          f_day[idx, ],
          by_element = TRUE
        )
      )
    )
    
  })
  
  out
  
})

#=================================================
# HISTOGRAMS BY YEAR
#=================================================

ggplot(
  all_results,
  aes(nearest_dist_m)
) +
  geom_histogram(
    bins = 100
  ) +
  facet_wrap(~year, scales = "free_y") +
  theme_bw() +
  labs(
    title = paste(
      op,
      "- Same-day nearest helicopter distance"
    ),
    x = "Nearest helicopter distance (m)",
    y = "Number of caribou fixes"
  )

#=================================================
# ENCOUNTER SUMMARY BY YEAR
#=================================================

thresholds <- c(
  100,
  250,
  500,
  1000,
  2000,
  5000
)

encounters_by_year <- expand.grid(
  year = sort(unique(all_results$year)),
  threshold_m = thresholds
) %>%
  rowwise() %>%
  mutate(
    n_encounters =
      sum(
        all_results$nearest_dist_m <= threshold_m &
          all_results$year == year,
        na.rm = TRUE
      )
  ) %>%
  ungroup()

encounters_by_year

#=================================================
# HEATMAP
#=================================================

ggplot(
  encounters_by_year,
  aes(
    factor(year),
    factor(threshold_m),
    fill = n_encounters
  )
) +
  geom_tile() +
  geom_text(aes(label = scales::comma(n_encounters))) +
  theme_bw() +
  labs(
    title = paste(op, "- Encounter counts"),
    x = "Year",
    y = "Distance threshold (m)"
  )


#=================================================
# THIN ALL FLIGHTS TO 1 POINT PER MINUTE
#=================================================

flights_thin <- flights %>%
  filter(!is.na(datetime)) %>%
  arrange(project, source_file, datetime) %>%
  mutate(
    minute_bin = floor_date(datetime, "1 minute")
  ) %>%
  group_by(project, source_file, minute_bin) %>%
  slice_head(n = 1) %>%
  ungroup()

#=================================================
# OPERATOR-YEAR COMBINATIONS
#=================================================

combos <- flights_thin %>%
  st_drop_geometry() %>%
  mutate(
    year = year(datetime)
  ) %>%
  distinct(
    project,
    year
  ) %>%
  arrange(
    project,
    year
  )

#=================================================
# SAME-DAY DISTANCE/TIME ANALYSIS
#=================================================

all_time_results <- map_dfr(
  seq_len(nrow(combos)),
  function(i){
    
    op <- combos$project[i]
    yr <- combos$year[i]
    
    cat(
      "\nProcessing",
      op,
      yr,
      "\n"
    )
    
    caribou_sub <- caribou %>%
      filter(
        year(datetime) == yr
      ) %>%
      mutate(
        date = as.Date(datetime)
      )
    
    flights_sub <- flights_thin %>%
      filter(
        project == op,
        year(datetime) == yr
      ) %>%
      mutate(
        date = as.Date(datetime)
      )
    
    dates <- intersect(
      unique(caribou_sub$date),
      unique(flights_sub$date)
    )
    
    if(length(dates) == 0){
      return(NULL)
    }
    
    map_dfr(
      dates,
      function(d){
        
        c_day <- caribou_sub %>%
          filter(date == d)
        
        f_day <- flights_sub %>%
          filter(date == d)
        
        idx <- st_nearest_feature(
          c_day,
          f_day
        )
        
        tibble(
          project = op,
          year = yr,
          date = d,
          
          nearest_dist_m = as.numeric(
            st_distance(
              c_day,
              f_day[idx, ],
              by_element = TRUE
            )
          ),
          
          time_diff_hr = abs(
            as.numeric(
              difftime(
                c_day$datetime,
                f_day$datetime[idx],
                units = "hours"
              )
            )
          )
        )
      }
    )
  }
)

#=================================================
# TIME BINS
#=================================================

all_time_results <- all_time_results %>%
  mutate(
    time_bin = cut(
      time_diff_hr,
      breaks = c(0, 1, 2, 6, 12, 24, Inf),
      labels = c(
        "<1 h",
        "1-2 h",
        "2-6 h",
        "6-12 h",
        "12-24 h",
        ">24 h"
      ),
      include.lowest = TRUE
    )
  )

#=================================================
# TIME BIN SUMMARY TABLE
#=================================================

time_bin_summary <- all_time_results %>%
  count(
    project,
    year,
    time_bin
  )

print(time_bin_summary)

#=================================================
# SAVE SUMMARY TABLE
#=================================================

write.csv(
  time_bin_summary,
  "time_bin_summary_all_operators_years.csv",
  row.names = FALSE
)

#=================================================
# DISTANCE BINS TOO (MAY BE USEFUL)
#=================================================

all_time_results <- all_time_results %>%
  mutate(
    distance_bin = cut(
      nearest_dist_m,
      breaks = c(0,100,250,500,1000,2000,5000,Inf),
      labels = c(
        "<100m",
        "100-250m",
        "250-500m",
        "500m-1km",
        "1-2km",
        "2-5km",
        ">5km"
      ),
      include.lowest = TRUE
    )
  )

distance_bin_summary <- all_time_results %>%
  count(
    project,
    year,
    distance_bin
  )

write.csv(
  distance_bin_summary,
  "distance_bin_summary_all_operators_years.csv",
  row.names = FALSE
)

#=================================================
# HISTOGRAM OF TIME DIFFERENCES
#=================================================

ggplot(
  all_time_results,
  aes(time_diff_hr)
) +
  geom_histogram(bins = 40) +
  facet_grid(
    project ~ year,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    title = "Time difference to nearest same-day helicopter point",
    x = "Time difference (hours)",
    y = "Number of caribou fixes"
  )

#=================================================
# HISTOGRAM OF DISTANCES
#=================================================

ggplot(
  all_time_results,
  aes(nearest_dist_m)
) +
  geom_histogram(bins = 100) +
  coord_cartesian(xlim = c(0, 20000)) +
  facet_grid(
    project ~ year,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    title = "Distance to nearest same-day helicopter point",
    x = "Nearest distance (m)",
    y = "Number of caribou fixes"
  )

#=================================================
# SAVE MASTER RESULTS
#=================================================

saveRDS(
  all_time_results,
  "all_same_day_encounters.rds"
)

#########################################################

caribou_same_day %>%
  filter(nearest_dist_m <= 5000) %>%
  st_drop_geometry() %>%
  summarise(
    n_caribou = n_distinct(COLLAR_)
  )

caribou_same_day %>%
  filter(nearest_dist_m <= 5000) %>%
  st_drop_geometry() %>%
  count(COLLAR_, sort = TRUE)


caribou_same_day %>%
  filter(nearest_dist_m <= 5000) %>%
  mutate(
    distance_bin = cut(
      nearest_dist_m,
      breaks = c(0,100,250,500,1000,2000,5000),
      labels = c(
        "<100m",
        "100-250m",
        "250-500m",
        "500m-1km",
        "1-2km",
        "2-5km"
      ),
      include.lowest = TRUE
    )
  ) %>%
  st_drop_geometry() %>%
  group_by(distance_bin) %>%
  summarise(
    n_encounters = n(),
    n_caribou = n_distinct(COLLAR_)
  )


caribou_same_day %>%
  filter(nearest_dist_m <= 5000) %>%
  st_drop_geometry() %>%
  count(COLLAR_, herd_nm, sort = TRUE) %>%
  mutate(
    pct = round(100 * n / sum(n), 1)
  )

caribou_same_day %>%
  filter(nearest_dist_m <= 5000) %>%
  st_drop_geometry() %>%
  summarise(
    n_encounters = n(),
    n_caribou = n_distinct(COLLAR_),
    n_herds = n_distinct(herd_nm)
  )

for(d in c(2000, 5000, 10000, 20000)) {
  
  tmp <- caribou_same_day %>%
    filter(nearest_dist_m <= d) %>%
    st_drop_geometry()
  
  cat(
    "\nDistance:", d, "m",
    "\nEncounters:", nrow(tmp),
    "\nCaribou:", n_distinct(tmp$COLLAR_),
    "\nHerds:", n_distinct(tmp$herd_nm),
    "\n"
  )
  
}