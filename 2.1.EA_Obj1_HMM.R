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

movement %>% summarise(collars = n_distinct(collar),locations = n(),steps = sum(!is.na(step_length)))
summary(movement$dt_hours)

saveRDS(movement, "movement.rds") # save

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

###############################################
# 1. PREPARE DATA
###############################################

movement <- readRDS("movement.rds") # load

hmm_data <- movement %>%
  filter(year %in% c(2022, 2023, 2024, 2025),
         dt_hours >= 6,
         dt_hours <= 10)

summary(hmm_data$dt_hours)

hmm_data %>%
  summarise(collars = n_distinct(collar),steps = n())

hmm_data %>%
  summarise(collars =n_distinct(collar),locations = n(),median_dt =median(dt_hours))

fixes_per_collar <- hmm_data %>%
  count(collar, name = "n_fixes") %>%
  arrange(n_fixes)

summary(fixes_per_collar$n_fixes)

ggplot(fixes_per_collar,
       aes(n_fixes)) +
  geom_histogram(binwidth = 25) +
  theme_bw()

fixes_per_collar %>%
  summarise(collars = n(),median = median(n_fixes),min = min(n_fixes),max = max(n_fixes))

good_collars <- fixes_per_collar %>%
  filter(n_fixes >= 50) %>%
  pull(collar)

hmm_sub <- hmm_data %>%
  filter(collar %in% good_collars)

hmm_track <- hmm_sub %>%
  mutate(ID = as.factor(collar),
    datetim = as.POSIXct(datetim)) %>%
  select(ID, datetim, x, y, ovrlps) %>%
  arrange(ID, datetim)

hmm_track %>% count(ID) # 18 collars with n fixes ranging from 55 to 2098

hmm_prep <- prepData(trackData = hmm_track,type = "UTM",coordNames = c("x", "y"))

class(hmm_prep)
names(hmm_prep)

###############################################
# 7. FIT 2-STATE MODEL
# Adjust starting values if necessary after inspecting quantiles
###############################################

mod2 <- fitHMM(
  data = hmm_prep,
  nbStates = 2,
  stepPar0 = c(200,1200,100,500),
  angleDist = "none"
)

###############################################
# INSPECT RESULTS
###############################################

print(mod2)

# State 1
# Mean step ≈ 511 m
# Persistent (96.3% chance of remaining in state)
# This looks like:localized use / foraging / residency
# State 2
# Mean step ≈ 2215 m
# Persistent (89.8% chance of remaining in state)
# This looks like: travelling or broader-scale movement

plot(mod2)

###############################################
# DECODE STATES
###############################################

hmm_track$state <- viterbi(mod2)

table(hmm_track$state)
prop.table(table(hmm_track$state))

table(hmm_track$state, hmm_track$ovrlps)
prop.table(table(hmm_track$state, hmm_track$ovrlps),margin = 2)

###############################################
# VISUALIZE STATES
###############################################
nrow(hmm_track); nrow(hmm_prep)
hmm_track$step <- hmm_prep$step

ggplot(hmm_track,
  aes(factor(state),step)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Decoded State",y = "Step Length (m)")

###############################################
# STATE OCCUPANCY
###############################################

occupancy <- hmm_track %>%
  count(ovrlps,state) %>%
  group_by(ovrlps) %>%
  mutate(proportion =n / sum(n))

occupancy

###############################################
# PLOT OCCUPANCY
###############################################

ggplot(occupancy,aes(factor(state),proportion, fill = ovrlps)) +
  geom_col(position = "dodge") +
  theme_bw()

summary(mod2)
print(mod2)
AIC(mod2)

###############################################
# 8. FIT 3-STATE MODEL
###############################################


mod3 <- fitHMM(
  data = hmm_prep,
  nbStates = 3,
  stepPar0 = c(150,800,2500,100,300,800),
  angleDist = "none")

print(mod3)

AIC(mod2)
AIC(mod3)

# mod3_1 <- fitHMM(
#   data = hmm_prep,
#   nbStates = 3,
#   stepPar0 = c(150,800,2500,100,300,800),
#   angleDist = "none"
# )
# 
# mod3_2 <- fitHMM(
#   data = hmm_prep,
#   nbStates = 3,
#   stepPar0 = c(100,500,2000,100,200,1000),
#   angleDist = "none"
# )
# 
# mod3_3 <- fitHMM(
#   data = hmm_prep,
#   nbStates = 3,
#   stepPar0 = c(300,1000,3000,200,500,1500),
#   angleDist = "none"
# )
# 
# AIC(mod3_1)
# AIC(mod3_2)
# AIC(mod3_3)
# 
# Multiple starting values were examined to assess sensitivity of model fitting. 
# Three independent runs converged on the same maximum likelihood solution.

###############################################
# MODEL COMPARISON
###############################################

AIC(mod2, mod3)

###############################################
# DECODE 3-STATE MODEL
###############################################

hmm_track$state3 <- viterbi(mod3)

table(hmm_track$state3)

hmm_track$state3 <- viterbi(mod3)

prop.table(table(hmm_track$state3))
prop.table(table(hmm_track$state3,hmm_track$ovrlps),margin = 2)

###############################################
# VISUALIZE 3 STATES
###############################################

ggplot(hmm_track,
  aes(factor(state3),step)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_bw()


# Occupancy table
occupancy3 <- hmm_track %>%
  count(ovrlps, state3) %>%
  group_by(ovrlps) %>%
  mutate(proportion = n / sum(n))

occupancy3

occupancy3 <- occupancy3 %>%
  mutate(state_label = factor(state3,levels = c(1, 2, 3),labels = c(
        "Localized\n(263 m)",
        "Moderate\n(661 m)",
        "Travelling\n(2373 m)")))

p_occ <- ggplot(
  occupancy3,aes(x = state_label,y = proportion,fill = ovrlps)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent,limits = c(0, 0.8)) +
  scale_fill_manual(values = c("FALSE" = "steelblue","TRUE" = "darkorange"),labels = c("Non-overlap","Overlap")) +
  labs(x = "Decoded Movement State",y = "Proportion of Locations",fill = "",title = "Occupancy of Three HMM Movement States") +
  theme_bw(base_size = 12)


ggsave("HMM_3state_occupancy.png",plot = p_occ,width = 7,height = 5,dpi = 300)

# Non-overlap animals are dominated by the moderate movement state.
# Overlap animals are dominated by the travelling state.
# Travelling-state occupancy increases from 17% to 65%.

###############################################
# END
###############################################

prop.table(table(hmm_track$state3,hmm_track$ovrlps),margin = 2)

library(lme4)
m3_high <- glmer(
  I(state3 == 3) ~ ovrlps + (1 | ID),
  data = hmm_track,
  family = binomial)

summary(m3_high)

exp(fixef(m3_high))
exp(confint(m3_high,parm = "beta_",method = "Wald"))

# ovrlpsTRUE = 2.345
exp(2.345) # 10.4; the overlap effect is substantial

# Conclusions
# A 2-state HMM successfully classified movements into low- and high-mobility states.
# A 3-state HMM was strongly supported over the 2-state model (ΔAIC ≈ 420).
# The three states corresponded to localized (263 m), moderate (661 m), and travelling (2373 m) movement behaviour.
# Overlap animals occupied the travelling state much more frequently (65%) than non-overlap animals (17%).
# A mixed-effects model indicated substantially greater odds of occupying the travelling state for overlap animals (OR ≈ 32, p = 0.005).
# These results suggest HMM-derived state metrics are likely suitable for subsequent disturbance-response analyses

