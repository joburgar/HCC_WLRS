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
# Exploratory Analysis - Obj 1 (HMM)
# 17 - July - 2026
# Joanna Burgar
############################################################

# Phase 1. Explore whether HMMs are feasible (using a subset of data)
# Question: Can Hart Ranges telemetry data be partitioned into biologically meaningful movement states using only GPS locations?
#   Do step lengths suggest more than one movement mode?
#   Are collar fix schedules sufficiently regular for HMM analysis?
#   Are turning angles informative?
#   Is there enough data for a 2-state HMM?
  
############################################################
# Load libraries
############################################################

library(tidyverse)
library(lubridate)
library(bcdata)
library(janitor)
library(sf)

############################################################
# HMM Feasability Assessment
############################################################

### Load GPS (telemetry) data

telemetry <- st_read("Outputs/pts_clean.shp")
names(telemetry)
# [1] "COLLAR_"  "datetim"  "dop"      "idmrtlt"  "activty"  "temprtr"  "manvltg"  "bckpvlt"  "X"       
# [10] "Y"        "med_x"    "med_y"    "dst_f__"  "dstnc_t"  "herd_nm"  "ecotype"  "ovrlps_"  "geometry"

collar <- st_drop_geometry(telemetry) %>% clean_names()


collar <- collar %>% 
  mutate(year = year(datetim))
collar <- collar %>% 
  mutate(winter_year =if_else(month(datetim) >= 11,year(datetim) + 1,year(datetim)))

names(collar)
### 1. Create movement metrics for Hart Ranges collared caribou

movement <- collar %>%
  filter(herd_nm == "Hart Ranges") %>%
  arrange(collar, datetim) %>%
  group_by(collar) %>%
  mutate(dt_hours =as.numeric(difftime(datetim,lag(datetim),units = "hours")),
         dx = x - lag(x),dy = y - lag(y),step_length =sqrt(dx^2 + dy^2),bearing =atan2(dy, dx)) %>%
  mutate(turning_angle =bearing - lag(bearing)) %>%
  ungroup()

### 2. Wrap turning angles

movement <- movement %>%
  mutate(turning_angle = atan2(sin(turning_angle),cos(turning_angle)))

### 3. Assess sample size

movement %>% 
  summarise(individuals = n_distinct(collar),locations = n(),steps = sum(!is.na(step_length)))

# individuals locations  steps
#          84    153393 153309

### 4. Evaluate fix intervals

ggplot(movement,aes(dt_hours)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(x = "Hours between fixes",y = "Count")


movement %>% 
  summarise(mean_fix = mean(dt_hours, na.rm = TRUE),
                       median_fix = median(dt_hours, na.rm = TRUE),
                       sd_fix = sd(dt_hours, na.rm = TRUE))

# mean_fix median_fix sd_fix
#     9.75          8   72.4

### 5. Explore step lengths

ggplot(movement,aes(step_length)) +
  geom_histogram(bins = 60) +
  scale_x_log10() +
  theme_bw()

ggplot(movement,aes(step_length,fill = ovrlps)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +
  theme_bw()

### 6. Explore turning angles

ggplot(movement,aes(turning_angle)) +
  geom_histogram(bins = 50) +
  theme_bw()

ggplot(movement,aes(turning_angle,fill = ovrlps)) +
  geom_density(alpha = 0.4) +
  theme_bw()

### 7. Examine annual variation
# full winter
hart_full <- collar %>%
  filter(herd_nm == "Hart Ranges", month(datetim) %in% c(11,12,1,2,3,4))
# core winter
hart_core <- collar %>%filter(herd_nm == "Hart Ranges",month(datetim) %in% c(12,1,2,3))


hart_full %>%
  group_by(winter_year) %>%
  summarise(collars =n_distinct(collar))

names(hart_full)
hart_full %>%
  group_by(winter_year) %>%
  summarise(fixes = n())
