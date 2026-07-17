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
library(moveHMM)


############################################################
# HMM Feasability Assessment
############################################################

### Load GPS (telemetry) data

telemetry <- st_read("Outputs/pts_clean.shp")
names(telemetry)
# [1] "COLLAR_"  "datetim"  "dop"      "idmrtlt"  "activty"  "temprtr"  "manvltg"  "bckpvlt"  "X"       
# [10] "Y"        "med_x"    "med_y"    "dst_f__"  "dstnc_t"  "herd_nm"  "ecotype"  "ovrlps_"  "geometry"

collar <- st_drop_geometry(telemetry) %>% clean_names()
rm(telemetry)

collar <- collar %>% 
  mutate(year = year(datetim))
collar <- collar %>% 
  mutate(winter_year =if_else(month(datetim) >= 11,year(datetim) + 1,year(datetim)))

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

hart_core %>%
  group_by(winter_year) %>%
  summarise(collars =n_distinct(collar))

hart_full %>% group_by(winter_year,collar) %>%
  summarise(fixes = n(),.groups = "drop")

hart_full %>%
  group_by(winter_year,collar) %>%
  summarise(fixes = n(),.groups = "drop") %>%
  ggplot(aes(factor(winter_year),fixes)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Winter",y = "Fixes per collar")

movement %>%
  group_by(winter_year) %>%
  summarise(median_dt =median(dt_hours,na.rm = TRUE),
    mean_dt =mean(dt_hours,na.rm = TRUE),
    sd_dt =sd(dt_hours,na.rm = TRUE))

ggplot(
  movement,aes(step_length)) +
  geom_histogram(bins = 100) +
  coord_cartesian(xlim = c(0, 5000)) +
  theme_bw()

quantile(movement$step_length,probs = c(0.1,0.25,0.5,0.75,0.9,0.95),na.rm = TRUE)

movement %>%summarise(collars = n_distinct(collar),locations = n(),steps = sum(!is.na(step_length)))
summary(movement$dt_hours)

movement %>%
  mutate(year = lubridate::year(datetim)) %>%
  group_by(year) %>%
  summarise(median_dt = median(dt_hours, na.rm = TRUE),mean_dt = mean(dt_hours, na.rm = TRUE),n = n())

movement %>%
  group_by(collar) %>%
  summarise(median_dt = median(dt_hours, na.rm = TRUE),n = n()) %>%
  arrange(median_dt)

# For HMM feasibility pilot, will use Hart Ranges data from the full winter (Nov-Apr) for 2023-2026
# 1. consistent fix rate,
# 2. largest sample sizes,
# 3. most collars,
# 4. most overlap/non-overlap animals,
# 5. retains November and April movements,
# 6. still directly relevant to heli-ski operations.
# 
# Then we can compare against the core winter dataset to see how much information is lost before deciding which seasonal definition to carry forward


############################################################
# HMM Pilot
############################################################
# Objective:
# Determine whether Hart Ranges telemetry data
# can be classified into biologically meaningful
# movement states using Hidden Markov Models.
#
# Clean data
# Create movement metrics
# Standardize to 8-hour fixes
# Examine movement characteristics
# Fit a 2-state HMM
# Decode movement states
# Compare tenure-overlap vs non-overlap animals
# Fit a 3-state HMM for comparison

hmm_data <- movement %>%
  filter(year %in% c(2022, 2023, 2024, 2025),
         dt_hours >= 6,
         dt_hours <= 10)

summary(hmm_data$dt_hours)

hmm_data %>%
  summarise(collars = n_distinct(collar),steps = n())

# ###############################################
# # 1. PREPARE DATA
# ###############################################
# 
# collar <- telemetry %>%
#   st_drop_geometry() %>%
#   clean_names() %>%
#   rename(
#     animal_id = collar_,
#     datetime = datetim,
#     overlap = ovrlps_
#   )
# 
###############################################
# 2. CREATE MOVEMENT METRICS
###############################################

movement <- collar %>%
  filter(herd_nm == "Hart Ranges") %>%
  arrange(animal_id, datetime) %>%
  group_by(animal_id) %>%
  
  mutate(
    
    year = year(datetime),
    
    dt_hours =
      as.numeric(
        difftime(
          datetime,
          lag(datetime),
          units = "hours"
        )
      ),
    
    dx = x - lag(x),
    dy = y - lag(y),
    
    step_length =
      sqrt(dx^2 + dy^2),
    
    bearing =
      atan2(dy, dx)
    
  ) %>%
  
  mutate(
    turning_angle =
      bearing - lag(bearing)
  ) %>%
  
  ungroup()

###############################################
# WRAP TURNING ANGLES
###############################################

movement <- movement %>%
  mutate(
    turning_angle =
      atan2(
        sin(turning_angle),
        cos(turning_angle)
      )
  )

###############################################
# 3. EXAMINE FIX INTERVALS
###############################################

summary(movement$dt_hours)

movement %>%
  group_by(year) %>%
  summarise(
    median_dt =
      median(dt_hours,
             na.rm = TRUE),
    mean_dt =
      mean(dt_hours,
           na.rm = TRUE),
    n = n()
  )

ggplot(
  movement,
  aes(dt_hours)
) +
  geom_histogram(
    bins = 100
  ) +
  theme_bw()

###############################################
# 4. STANDARDIZE TO 8-HOUR FIXES
###############################################

hmm_data <- movement %>%
  filter(
    year %in% c(2022, 2023, 2024, 2025),
    dt_hours >= 6,
    dt_hours <= 10,
    !is.na(step_length),
    step_length > 0
  )

###############################################
# SAMPLE SIZE CHECK
###############################################

hmm_data %>%
  summarise(
    collars =
      n_distinct(animal_id),
    
    locations = n(),
    
    median_dt =
      median(dt_hours)
  )

###############################################
# 5. STEP LENGTH EXPLORATION
###############################################

summary(hmm_data$step_length)

quantile(
  hmm_data$step_length,
  probs = c(
    0.25,
    0.50,
    0.75,
    0.90
  ),
  na.rm = TRUE
)

ggplot(
  hmm_data,
  aes(step_length)
) +
  geom_histogram(
    bins = 100
  ) +
  coord_cartesian(
    xlim = c(0,5000)
  ) +
  theme_bw()

ggplot(
  hmm_data,
  aes(
    step_length,
    fill = overlap
  )
) +
  geom_density(
    alpha = 0.4
  ) +
  scale_x_log10() +
  theme_bw()

###############################################
# 6. PREPARE HMM DATA
###############################################

mod_data <- data.frame(
  ID = hmm_data$animal_id,
  step = hmm_data$step_length,
  angle = 0
)

###############################################
# STARTING VALUES
#
# Adjust if necessary after
# inspecting quantiles.
###############################################

stepPar0 <- c(
  200,
  1200,
  100,
  500
)

###############################################
# 7. FIT 2-STATE MODEL
###############################################

mod2 <- fitHMM(
  data = mod_data,
  
  nbStates = 2,
  
  dist = list(
    step = "gamma",
    angle = "wrpcauchy"
  ),
  
  Par0 = list(
    step = stepPar0,
    angle = c(0.1, 0.1)
  )
)

###############################################
# INSPECT RESULTS
###############################################

print(mod2)

plot(mod2)

###############################################
# DECODE STATES
###############################################

hmm_data$state <- viterbi(mod2)

table(hmm_data$state)

###############################################
# VISUALIZE STATES
###############################################

ggplot(
  hmm_data,
  aes(
    factor(state),
    step_length
  )
) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(
    x = "Decoded State",
    y = "Step Length (m)"
  )

###############################################
# STATE OCCUPANCY
###############################################

occupancy <- hmm_data %>%
  count(
    overlap,
    state
  ) %>%
  
  group_by(overlap) %>%
  
  mutate(
    proportion =
      n / sum(n)
  )

occupancy

###############################################
# PLOT OCCUPANCY
###############################################

ggplot(
  occupancy,
  aes(
    factor(state),
    proportion,
    fill = overlap
  )
) +
  geom_col(
    position = "dodge"
  ) +
  theme_bw()

###############################################
# 8. FIT 3-STATE MODEL
###############################################

stepPar0_3 <- c(
  150,
  800,
  2500,
  100,
  300,
  800
)

mod3 <- fitHMM(
  data = mod_data,
  
  nbStates = 3,
  
  dist = list(
    step = "gamma",
    angle = "wrpcauchy"
  ),
  
  Par0 = list(
    step = stepPar0_3,
    angle = c(
      0.1,
      0.1,
      0.1
    )
  )
)

###############################################
# MODEL COMPARISON
###############################################

AIC(mod2, mod3)

###############################################
# DECODE 3-STATE MODEL
###############################################

hmm_data$state3 <- viterbi(mod3)

table(hmm_data$state3)

###############################################
# VISUALIZE 3 STATES
###############################################

ggplot(
  hmm_data,
  aes(
    factor(state3),
    step_length
  )
) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw()

###############################################
# END
###############################################
