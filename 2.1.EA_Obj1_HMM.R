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

caribou <- collar %>%
  arrange(collar, datetim) %>%
  group_by(collar) %>%
  mutate(dt_hours =as.numeric(difftime(datetim,lag(datetim),units = "hours")),
         dx = x - lag(x),dy = y - lag(y),step_length =sqrt(dx^2 + dy^2),bearing =atan2(dy, dx)) %>%
  mutate(turning_angle =bearing - lag(bearing)) %>%
  ungroup()

saveRDS(caribou, "caribou.rds") # save

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


##################################################################################
##################################################################################
### Continue validation now that HMM seems useful with Hart Ranges herd
##################################################################################
##################################################################################
# STEP 1: Continued HMM validation
#
# Goal:
# Determine whether the 3-state winter movement structure is robust across:
#   1) Herds
#   2) Leave-one-herd-out analyses
#   3) Winter years
#   4) Individuals
#
##################################################################################
################################################################################

################################################################################
# USER INPUTS
################################################################################

focal_herds <- c(
  "Hart Ranges",
  "Columbia North",
  "Central Selkirks",
  "Wells Gray North",
  "North Cariboo",
  "Itcha-Ilgachuz"
)

################################################################################
# PREP DATA
################################################################################

caribou <- readRDS("caribou.rds") # load
names(caribou)
# [1] "collar"        "datetim"       "dop"           "idmrtlt"       "activty"       "temprtr"      
# [7] "manvltg"       "bckpvlt"       "x"             "y"             "med_x"         "med_y"        
# [13] "dst_f"         "dstnc_t"       "herd_nm"       "ecotype"       "ovrlps"        "year"         
# [19] "winter_year"   "dt_hours"      "dx"            "dy"            "step_length"   "bearing"      
# [25] "turning_angle"

caribou_hmm <- caribou %>%
  filter(herd_nm %in% focal_herds) %>%
  mutate(ID = as.factor(collar),
         datetim = as.POSIXct(datetim)) %>%
  mutate(month  = month(datetim),year   = year(datetim),
         # define biological winter
    winter = case_when(month >= 11 ~ year + 1,TRUE ~ year))%>%
  filter(dt_hours >= 6, dt_hours <= 12) %>%
  filter(month(datetim) %in% c(11,12,1,2,3,4)) %>%  # clarify winter months
  filter(winter %in% c(2021,2022,2023,2024,2025,2026)) # clarify winters for analysis

summary(caribou_hmm$dt_hours)
summary(caribou_hmm$winter)
caribou_hmm %>% count(month)

caribou_hmm %>% group_by(herd_nm) %>%
  summarise(collars = n_distinct(collar),steps = n()) 

caribou_hmm %>% group_by(herd_nm, winter) %>%
  summarise(collars =n_distinct(collar),locations = n(),median_dt =median(dt_hours))%>% print(n=40)

fixes_per_collar <- caribou_hmm %>% group_by(herd_nm) %>%
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

caribou_track <- caribou_hmm %>%
  filter(collar %in% good_collars) %>%
  mutate(ID = as.factor(collar),
         datetim = as.POSIXct(datetim)) %>%
  select(herd_nm, ID, datetim, x, y, ovrlps) %>%
  arrange(ID, datetim)

################################################################################
# STEP 1: VALIDATE 3-STATE WINTER HMM
#
# Goal:
# Determine whether the 3-state movement-state structure is robust across:
#
# 1. Herd-specific models
# 2. Leave-one-herd-out (LOHO) models
# 3. Winter-year models
# 4. Individual-level models
#
################################################################################
# EXTRACT MODEL SUMMARY
################################################################################

extract_summary <- function(mod,
                            group_name,
                            group_type){
  if(is.null(mod)){
    return(NULL)
  }
  state_means <- mod$mle$stepPar[1,]
  state_order <- order(state_means)
  state_means <- state_means[state_order]
  tibble(
    group_type = group_type,
    group      = group_name,
    state      = 1:3,
    mean_step  = state_means,
    AIC        = AIC(mod)
  )
}

################################################################################
# EXTRACT STATE OCCUPANCY
################################################################################

extract_occupancy <- function(mod,
                              group_name,
                              group_type){
  if(is.null(mod)){
    return(NULL)
  }
  v <- viterbi(mod)
  occ <- prop.table(table(v))
  tibble(
    group_type = group_type,
    group      = group_name,
    state      = names(occ),
    occupancy  = as.numeric(occ)
  )
}

################################################################################
# 1. HERD-SPECIFIC MODELS
################################################################################

herd_models <- list()

for(h in focal_herds){
  cat("Running:", h, "\n")
  # h = focal_herds[1]
  
  hmm_track <- caribou_track %>%
    filter(herd_nm == h) %>%
    select(ID, datetim, x, y, ovrlps)
  
  dat_h <- prepData(trackData = hmm_track,type = "UTM",coordNames = c("x", "y"))
  dat_h <- dat_h %>% filter(step > 0)
  
  herd_models[[h]] <- fitHMM(
    data = dat_h,
    nbStates = 3,
    # stepPar0 = c(300,1000,3000,200,500,1500), # stepPar3
    # stepPar0 = c(100,500,2000,100,200,1000), # stepPar2
    stepPar0 = c(150,800,2500,100,300,800), # stepPar1
    angleDist = "none")
}

herd_results <- bind_rows(
  
  lapply(
    names(herd_models),
    function(h){
      
      extract_summary(
        herd_models[[h]],
        h,
        "herd")}))


herd_models_stepPar3 <- herd_models
herd_results_stepPar3 <- herd_results

herd_models_stepPar2 <- herd_models
herd_results_stepPar2 <- herd_results

herd_models_stepPar1 <- herd_models
herd_results_stepPar1 <- herd_results

### exporting / saving so can start from here without re-running models ###
# saveRDS(herd_models_stepPar1, "Outputs/herd_models_stepPar1.rds")
# saveRDS(herd_models_stepPar2, "Outputs/herd_models_stepPar2.rds")
# saveRDS(herd_models_stepPar3, "Outputs/herd_models_stepPar3.rds")
# 
# saveRDS(herd_results_stepPar1, "Outputs/herd_results_stepPar1.rds")
# saveRDS(herd_results_stepPar2, "Outputs/herd_results_stepPar2.rds")
# saveRDS(herd_results_stepPar3, "Outputs/herd_results_stepPar3.rds")

herd_models_stepPar1 <- readRDS("Outputs/herd_models_stepPar1.rds")
herd_models_stepPar2 <- readRDS("Outputs/herd_models_stepPar2.rds")
herd_models_stepPar3 <- readRDS("Outputs/herd_models_stepPar3.rds")

herd_results_stepPar1 <- readRDS("Outputs/herd_results_stepPar1.rds")
herd_results_stepPar2 <- readRDS("Outputs/herd_results_stepPar2.rds")
herd_results_stepPar3 <- readRDS("Outputs/herd_results_stepPar3.rds")

compare_starts <- bind_rows(
  herd_results_stepPar1 %>% mutate(start = "par1"),
  herd_results_stepPar2 %>% mutate(start = "par2"),
  herd_results_stepPar3 %>% mutate(start = "par3")
)

compare_starts_wide <- compare_starts %>%
  select(group, start, state, mean_step, AIC) %>%
  tidyr::pivot_wider(
    names_from = state,
    values_from = mean_step,
    names_prefix = "State_"
  ) %>%
  arrange(group, start)

compare_starts_wide # all look nearly identicial, except Columbia North

# going with starting levels for stepPar1
# stepPar0 = c(150,800,2500,100,300,800), # stepPar1

compare_starts_wide %>%
  group_by(group) %>%
  summarise(
    mean_state1 = mean(State_1),
    mean_state2 = mean(State_2),
    mean_state3 = mean(State_3),
    cv_state1 = sd(State_1)/mean(State_1),
    cv_state2 = sd(State_2)/mean(State_2),
    cv_state3 = sd(State_3)/mean(State_3)
  )

compare_starts_wide %>%
  filter(group == "Columbia North")

caribou %>%
  filter(herd_nm %in% focal_herds) %>%
  group_by(herd_nm) %>%
  summarise(
    mean_step = mean(step_length, na.rm = TRUE),
    median_step = median(step_length, na.rm = TRUE),
    p95 = quantile(step_length, 0.95, na.rm = TRUE)
  )

write.csv(compare_starts_wide,"outputs/compare_starts_wide.csv",row.names = FALSE)

# Passing the starting-value test: Yes.
# Evidence of instability: Only Columbia North's par3 fit, and the huge AIC difference makes it easy to reject.
# Biggest issue revealed so far: Not starting values. It's that Central Selkirks and especially Itcha-Ilgachuz appear to have much larger movement scales than the other herds.

################################################################################
# 2. LEAVE-ONE-HERD-OUT MODELS
################################################################################

loho_models <- list()

for(h in focal_herds){
  
  cat("Leaving out:", h, "\n")
  
  hmm_track <- caribou_track %>%
    filter(herd_nm != h) %>%
    select(herd_nm, ID, datetim, x, y, ovrlps)
  
  dat_h <- prepData(trackData = hmm_track,type = "UTM",coordNames = c("x", "y"))
  dat_h <- dat_h %>% filter(step > 0)
  
  loho_models[[h]] <- fitHMM(
    data = dat_h,
    nbStates = 3,
    stepPar0 = c(150,800,2500,100,300,800), # stepPar1 (from initial tests)
    angleDist = "none")
  
}

loho_results <- bind_rows(
  lapply(
    names(loho_models),
    function(h){
      
      extract_summary(
        loho_models[[h]],
        paste0("Exclude_", h),
        "LOHO")})
)


# saveRDS(loho_results, "Outputs/loho_results.rds")
loho_results <- readRDS("Outputs/loho_results.rds")


loho_wide <- loho_results %>%
  select(group, state, mean_step, AIC) %>%
  pivot_wider(
    names_from = state,
    values_from = mean_step,
    names_prefix = "State_"
  ) %>%
  rename(AIC = AIC)

loho_wide

loho_wide %>%
  summarise(State1_min = min(State_1),State1_max = max(State_1),
    State2_min = min(State_2),State2_max = max(State_2),
    State3_min = min(State_3),
    State3_max = max(State_3))

loho_plot <- ggplot(loho_results,aes(group, mean_step)) +
  geom_point(size = 3) +
  facet_wrap(~state, scales = "free_y") +
  coord_flip() +
  theme_bw()

ggsave("loho_plot.png",plot = loho_plot,width = 10,height = 5,dpi = 300)

# Leave-one-herd-out analyses recovered the same three-state structure regardless of which herd was excluded.
# Estimated step lengths for the localized and intermediate states were relatively stable among analyses (~15-179 m and ~497-701 m, respectively)
# The travelling state showed greater variability (~1968-3281 m), largely due to the influence of the Itcha-Ilgachuz herd, which exhibited substantially larger movement scales than other focal herds.
# Exclusion of any single herd did not alter recovery of the three-state solution, suggesting that HMM state structure is broadly robust across study areas.

################################################################################
# 3. WINTER-YEAR MODELS
################################################################################

winter_models <- list()

winters <- sort(unique(hmm_data$winter_year))

for(w in winters){
  
  cat("Winter:", w, "\n")
  
  dat_w <- hmm_data %>%
    filter(winter_year == w)
  
  winter_models[[as.character(w)]] <- fit_3state_hmm(
    dat_w,
    stepPar0 = stepPar0
  )
  
}

winter_results <- bind_rows(
  
  lapply(
    names(winter_models),
    function(w){
      
      extract_summary(
        winter_models[[w]],
        w,
        "winter"
      )
      
    }
  )
  
)

################################################################################
# 4. INDIVIDUAL-LEVEL HMMS
################################################################################

good_ids <- hmm_data %>%
  count(ID) %>%
  filter(n >= min_fixes_individual) %>%
  pull(ID)

individual_models <- list()

for(id in good_ids){
  
  cat("Individual:", as.character(id), "\n")
  
  dat_i <- hmm_data %>%
    filter(ID == id)
  
  individual_models[[as.character(id)]] <-
    fit_3state_hmm(
      dat_i,
      stepPar0 = stepPar0
    )
  
}

individual_results <- bind_rows(
  
  lapply(
    names(individual_models),
    function(id){
      
      extract_summary(
        individual_models[[id]],
        id,
        "individual"
      )
      
    }
  )
  
)

################################################################################
# COMBINE RESULTS
################################################################################

robustness_summary <- bind_rows(
  herd_results,
  loho_results,
  winter_results,
  individual_results
)

################################################################################
# STATE STABILITY
################################################################################

state_stability <- robustness_summary %>%
  group_by(group_type, state) %>%
  summarise(
    mean_step = mean(mean_step, na.rm = TRUE),
    sd_step   = sd(mean_step, na.rm = TRUE),
    cv_step   = sd_step / mean_step,
    .groups = "drop"
  )

################################################################################
# OCCUPANCY SUMMARIES
################################################################################

herd_occupancy <- bind_rows(
  
  lapply(
    names(herd_models),
    function(h){
      
      extract_occupancy(
        herd_models[[h]],
        h,
        "herd"
      )
      
    }
  )
  
)

################################################################################
# REVIEW OUTPUTS
################################################################################

print(state_stability)

head(robustness_summary)

head(herd_occupancy)

################################################################################
# SAVE OUTPUTS
################################################################################

dir.create("outputs", showWarnings = FALSE)

write.csv(
  robustness_summary,
  "outputs/HMM_robustness_summary.csv",
  row.names = FALSE
)

write.csv(
  state_stability,
  "outputs/HMM_state_stability.csv",
  row.names = FALSE
)

write.csv(
  herd_occupancy,
  "outputs/HMM_state_occupancy.csv",
  row.names = FALSE
)

################################################################################
# OPTIONAL VISUAL CHECK
################################################################################

library(ggplot2)

ggplot(
  robustness_summary,
  aes(state,
      mean_step,
      colour = group_type)
) +
  geom_point(size = 2) +
  geom_line(aes(group = group)) +
  theme_bw()

################################################################################