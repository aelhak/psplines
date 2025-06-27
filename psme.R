
#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(tidyverse); library(psme); library(growthstandards)

#--------------------#
#### LOAD DATASET ####
#--------------------#

load("gusto_ps_dat.RData")

#-----------------------------------#
#### FIT P-SPLINE MODELS IN PSME ####
#-----------------------------------#

psme_ht_M <- psme(
  cht ~ s(sqrt_age, bs = 'ps', k = 10, m = c(2, 2)) + 
    s(sqrt_age, id, bs = "fs", xt = "ps", k = 10,  m = c(2, 1)), 
  data = dat_M)

psme_wt_M <- psme(
  cwt ~ s(sqrt_age, bs = 'ps', k = 10, m = c(2, 2)) + 
    s(sqrt_age, id, bs = "fs", xt = "ps", k = 10,  m = c(2, 1)), 
  data = dat_M)

psme_bmi_M <- psme(
  cbmi ~ s(sqrt_age, bs = 'ps', k = 10, m = c(2, 2)) + 
    s(sqrt_age, id, bs = "fs", xt = "ps", k = 10,  m = c(2, 1)), 
  data = dat_M)

psme_ht_F <- psme(
  cht ~ s(sqrt_age, bs = 'ps', k = 10, m = c(2, 2)) + 
    s(sqrt_age, id, bs = "fs", xt = "ps", k = 10,  m = c(2, 1)), 
  data = dat_F)

psme_wt_F <- psme(
  cwt ~ s(sqrt_age, bs = 'ps', k = 10, m = c(2, 2)) + 
    s(sqrt_age, id, bs = "fs", xt = "ps", k = 10,  m = c(2, 1)), 
  data = dat_F)

psme_bmi_F <- psme(
  cbmi ~ s(sqrt_age, bs = 'ps', k = 10, m = c(2, 2)) + 
    s(sqrt_age, id, bs = "fs", xt = "ps", k = 10,  m = c(2, 1)), 
  data = dat_F)

save.image("gusto_psme_mods.RData")

rm(list = ls())

#----------------------------------#
#### GET PSME MODEL PREDICTIONS ####
#----------------------------------#

load("gusto_psme_mods.RData")
source("src/psme_helper_functions.R")

# MALES

age.M <- with(dat_M, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))

psme_pred_M <- as.data.frame(xTraj(psme_ht_M, age.M)$sub) %>% mutate(age = age.M^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cht")

psme_pred_M <- as.data.frame(xTraj(psme_wt_M, age.M)$sub) %>% mutate(age = age.M^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cwt") %>% 
  full_join(psme_pred_M)

psme_pred_M <- as.data.frame(xTraj(psme_bmi_M, age.M)$sub) %>% mutate(age = age.M^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cbmi") %>% 
  full_join(psme_pred_M)

psme_pred_M$id <- as.integer(sub("^\\D+", "", psme_pred_M$id))

# FEMALES

age.F <- with(dat_F, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))

psme_pred_F <- as.data.frame(xTraj(psme_ht_F, age.F)$sub) %>% mutate(age = age.F^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cht")

psme_pred_F <- as.data.frame(xTraj(psme_wt_F, age.F)$sub) %>% mutate(age = age.F^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cwt") %>% 
  full_join(psme_pred_F)

psme_pred_F <- as.data.frame(xTraj(psme_bmi_F, age.F)$sub) %>% mutate(age = age.F^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cbmi") %>% 
  full_join(psme_pred_F)

psme_pred_F$id <- as.integer(sub("^\\D+", "", psme_pred_F$id))

rm(dat_M, dat_F, age.M, age.F, psme_ht_M, psme_ht_F, psme_wt_M, psme_wt_F, 
   psme_bmi_M, psme_bmi_F, xTraj, AUC, peak)

save.image("gusto_psme_preds.RData")

rm(list = ls())

#---------------------------#
#### GET GROWTH FEATURES ####
#---------------------------#

load("gusto_ps_dat.RData")
load("gusto_psme_preds.RData")
source("src/psme_helper_functions.R")

#### INFANT PEAK BMI - BOYS ####

peak_M <- psme_pred_M %>% select(id, age, cbmi) %>% filter(age <= 1.5)
peak_M <- split(peak_M, peak_M$id)

est_M <- map_df(peak_M, ~ sitar::getPeak(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = unique(psme_pred_M$id)) %>% relocate(id) %>% rename(
    PeakBMIAge = x, PeakBMI = y) %>% mutate(sex = "Boys")

rm(peak_M)
hist(est_M$PeakBMIAge)
summary(est_M$PeakBMIAge) # NA = 19

#### CHILD REBOUND BMI - BOYS ####

rebound_M <- psme_pred_M %>% select(id, age, cbmi) %>% filter(age > 2.5)
rebound_M <- split(rebound_M, rebound_M$id)

est_M <- map_df(rebound_M, ~ sitar::getTrough(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = unique(est_M$id)) %>% relocate(id) %>% rename(
    ReboundBMIAge = x, ReboundBMI = y) %>% full_join(est_M)

rm(rebound_M)
hist(est_M$ReboundBMIAge)
summary(est_M$ReboundBMIAge) # NA = 64

#### INFANT PEAK HEIGHT VELOCITY - BOYS ####

est_M <- psme_pred_M %>% select(id, age, cht) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PHV = max(dY)) %>% full_join(est_M)

summary(est_M$PHV)
hist(est_M$PHV)

#### INFANT PEAK WEIGHT VELOCITY - BOYS ####

est_M <- psme_pred_M %>% select(id, age, cwt) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PWV = max(dY)) %>% full_join(est_M)

summary(est_M$PWV)
hist(est_M$PWV)

#### HAD, WAD - BOYS ####

load("gusto_psme_mods.RData")

age.M <- c(
  sqrt(0.0833334), sqrt(0.5), sqrt(1), sqrt(1.5), sqrt(2), sqrt(2.5),
  sqrt(3), sqrt(3.5), sqrt(4), sqrt(4.5), sqrt(5))

psme_pred_M <- as.data.frame(xTraj(psme_ht_M, age.M)$sub) %>% mutate(age = age.M^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cht")

psme_pred_M <- as.data.frame(xTraj(psme_wt_M, age.M)$sub) %>% mutate(age = age.M^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cwt") %>% 
  full_join(psme_pred_M)

psme_pred_M$id <- as.integer(sub("^\\D+", "", psme_pred_M$id))

psme_pred_M <- psme_pred_M %>% 
  mutate(age = ifelse(
    age == '0.0833334', 1, ifelse(
      age == '0.5', 2, ifelse(
        age == '1', 3, ifelse(
          age == '1.5', 4, ifelse(
            age == '2', 5, ifelse(
              age == '2.5', 6, ifelse(
                age == '3', 7, ifelse(
                  age == '3.5', 8, ifelse(
                    age == '4', 9, ifelse(
                      age == '4.5', 10, 11
                    )))))))))))

table(psme_pred_M$age)

psme_pred_M <- psme_pred_M %>% mutate(
  HAD_1m = ifelse(age == 1, cht - 54.7, NA),
  HAD_6m = ifelse(age == 2, cht - 67.6, NA),
  HAD_12m = ifelse(age == 3, cht - 75.7, NA),
  HAD_18m = ifelse(age == 4, cht - 82.3, NA),
  HAD_24m = ifelse(age == 5, cht - 87.8, NA),
  HAD_30m = ifelse(age == 6, cht - 91.9, NA),
  HAD_36m = ifelse(age == 7, cht - 96.1, NA),
  HAD_42m = ifelse(age == 8, cht - 99.9, NA),
  HAD_48m = ifelse(age == 9, cht - 103.3, NA),
  HAD_54m = ifelse(age == 10, cht - 106.7, NA),
  HAD_60m = ifelse(age == 11, cht - 110, NA),
  WAD_1m = ifelse(age == 1, cwt - 4.5, NA),
  WAD_6m = ifelse(age == 2, cwt - 7.9, NA),
  WAD_12m = ifelse(age == 3, cwt - 9.6, NA),
  WAD_18m = ifelse(age == 4, cwt - 10.9, NA),
  WAD_24m = ifelse(age == 5, cwt - 12.2, NA),
  WAD_30m = ifelse(age == 6, cwt - 13.3, NA),
  WAD_36m = ifelse(age == 7, cwt - 14.3, NA),
  WAD_42m = ifelse(age == 8, cwt - 15.3, NA),
  WAD_48m = ifelse(age == 9, cwt - 16.3, NA),
  WAD_54m = ifelse(age == 10, cwt - 17.3, NA),
  WAD_60m = ifelse(age == 11, cwt - 18.3, NA),
) %>% select(-cwt, -cht, -age)

est_M <- psme_pred_M %>% select(id, HAD_1m, WAD_1m) %>% drop_na() %>% full_join(est_M)
est_M <- psme_pred_M %>% select(id, HAD_6m, WAD_6m) %>% drop_na() %>% full_join(est_M)
est_M <- psme_pred_M %>% select(id, HAD_12m, WAD_12m) %>% drop_na() %>% full_join(est_M)
est_M <- psme_pred_M %>% select(id, HAD_18m, WAD_18m) %>% drop_na() %>% full_join(est_M)
est_M <- psme_pred_M %>% select(id, HAD_24m, WAD_24m) %>% drop_na() %>% full_join(est_M)
est_M <- psme_pred_M %>% select(id, HAD_30m, WAD_30m) %>% drop_na() %>% full_join(est_M)
est_M <- psme_pred_M %>% select(id, HAD_36m, WAD_36m) %>% drop_na() %>% full_join(est_M)
est_M <- psme_pred_M %>% select(id, HAD_42m, WAD_42m) %>% drop_na() %>% full_join(est_M)
est_M <- psme_pred_M %>% select(id, HAD_48m, WAD_48m) %>% drop_na() %>% full_join(est_M)
est_M <- psme_pred_M %>% select(id, HAD_54m, WAD_54m) %>% drop_na() %>% full_join(est_M)
est_M <- psme_pred_M %>% select(id, HAD_60m, WAD_60m) %>% drop_na() %>% full_join(est_M)

#### INFANT PEAK BMI - GIRLS ####

peak_F <- psme_pred_F %>% select(id, age, cbmi) %>% filter(age <= 1.5)
peak_F <- split(peak_F, peak_F$id)

est_F <- map_df(peak_F, ~ sitar::getPeak(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = unique(psme_pred_F$id)) %>% relocate(id) %>% rename(
    PeakBMIAge = x, PeakBMI = y) %>% mutate(sex = "Girls")

rm(peak_F)
hist(est_F$PeakBMIAge)
summary(est_F$PeakBMIAge) # NA=13

#### CHILD REBOUND BMI - GIRLS ####

rebound_F <- psme_pred_F %>% select(id, age, cbmi) %>% filter(age >= 3)
rebound_F <- split(rebound_F, rebound_F$id)

est_F <- map_df(rebound_F, ~ sitar::getTrough(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = unique(est_F$id)) %>% relocate(id) %>% rename(
    ReboundBMIAge = x, ReboundBMI = y) %>% full_join(est_F) 

rm(rebound_F)
hist(est_F$ReboundBMIAge)
summary(est_F$ReboundBMIAge) # NA=80

#### INFANT PEAK HEIGHT VELOCITY - GIRLS ####

est_F <- psme_pred_F %>% select(id, age, cht) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PHV = max(dY)) %>% full_join(est_F) 

summary(est_F$PHV)
hist(est_F$PHV)

#### INFANT PEAK WEIGHT VELOCITY - GIRLS ####

est_F <- psme_pred_F %>% select(id, age, cwt) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PWV = max(dY)) %>% full_join(est_F)

summary(est_F$PWV)
hist(est_F$PWV)

#### HAD, WAD - GIRLS ####

age.F <- c(
  sqrt(0.0833334), sqrt(0.5), sqrt(1), sqrt(1.5), sqrt(2), sqrt(2.5),
  sqrt(3), sqrt(3.5), sqrt(4), sqrt(4.5), sqrt(5))

psme_pred_F <- as.data.frame(xTraj(psme_ht_F, age.F)$sub) %>% mutate(age = age.F^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cht")

psme_pred_F <- as.data.frame(xTraj(psme_wt_F, age.F)$sub) %>% mutate(age = age.F^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cwt") %>% 
  full_join(psme_pred_F)

psme_pred_F$id <- as.integer(sub("^\\D+", "", psme_pred_F$id))

psme_pred_F <- psme_pred_F %>% 
  mutate(age = ifelse(
    age == '0.0833334', 1, ifelse(
      age == '0.5', 2, ifelse(
        age == '1', 3, ifelse(
          age == '1.5', 4, ifelse(
            age == '2', 5, ifelse(
              age == '2.5', 6, ifelse(
                age == '3', 7, ifelse(
                  age == '3.5', 8, ifelse(
                    age == '4', 9, ifelse(
                      age == '4.5', 10, 11
                    )))))))))))

table(psme_pred_F$age)

psme_pred_F <- psme_pred_F %>% mutate(
  HAD_1m = ifelse(age == 1, cht - 53.7, NA),
  HAD_6m = ifelse(age == 2, cht - 65.7, NA),
  HAD_12m = ifelse(age == 3, cht - 74.0, NA),
  HAD_18m = ifelse(age == 4, cht - 80.7, NA),
  HAD_24m = ifelse(age == 5, cht - 86.4, NA),
  HAD_30m = ifelse(age == 6, cht - 90.7, NA),
  HAD_36m = ifelse(age == 7, cht - 95.1, NA),
  HAD_42m = ifelse(age == 8, cht - 99.0, NA),
  HAD_48m = ifelse(age == 9, cht - 102.7, NA),
  HAD_54m = ifelse(age == 10, cht - 106.2, NA),
  HAD_60m = ifelse(age == 11, cht - 109.4, NA),
  WAD_1m = ifelse(age == 1, cwt - 4.2, NA),
  WAD_6m = ifelse(age == 2, cwt - 7.3, NA),
  WAD_12m = ifelse(age == 3, cwt - 8.9, NA),
  WAD_18m = ifelse(age == 4, cwt - 10.2, NA),
  WAD_24m = ifelse(age == 5, cwt - 11.5, NA),
  WAD_30m = ifelse(age == 6, cwt - 12.7, NA),
  WAD_36m = ifelse(age == 7, cwt - 13.9, NA),
  WAD_42m = ifelse(age == 8, cwt - 15, NA),
  WAD_48m = ifelse(age == 9, cwt - 16.1, NA),
  WAD_54m = ifelse(age == 10, cwt - 17.2, NA),
  WAD_60m = ifelse(age == 11, cwt - 18.2, NA),
) %>% select(-cwt, -cht, -age)

est_F <- psme_pred_F %>% select(id, HAD_1m, WAD_1m) %>% drop_na() %>% full_join(est_F)
est_F <- psme_pred_F %>% select(id, HAD_6m, WAD_6m) %>% drop_na() %>% full_join(est_F)
est_F <- psme_pred_F %>% select(id, HAD_12m, WAD_12m) %>% drop_na() %>% full_join(est_F)
est_F <- psme_pred_F %>% select(id, HAD_18m, WAD_18m) %>% drop_na() %>% full_join(est_F)
est_F <- psme_pred_F %>% select(id, HAD_24m, WAD_24m) %>% drop_na() %>% full_join(est_F)
est_F <- psme_pred_F %>% select(id, HAD_30m, WAD_30m) %>% drop_na() %>% full_join(est_F)
est_F <- psme_pred_F %>% select(id, HAD_36m, WAD_36m) %>% drop_na() %>% full_join(est_F)
est_F <- psme_pred_F %>% select(id, HAD_42m, WAD_42m) %>% drop_na() %>% full_join(est_F)
est_F <- psme_pred_F %>% select(id, HAD_48m, WAD_48m) %>% drop_na() %>% full_join(est_F)
est_F <- psme_pred_F %>% select(id, HAD_54m, WAD_54m) %>% drop_na() %>% full_join(est_F)
est_F <- psme_pred_F %>% select(id, HAD_60m, WAD_60m) %>% drop_na() %>% full_join(est_F)

#-----------------------#
#### COMBINE & MERGE ####
#-----------------------#

est_M <- dat_M %>% mutate(id = as.integer(id)) %>% select(
  -cht, -cwt, -cbmi, -age, -sqrt_age, -sex) %>% distinct(id, .keep_all = TRUE) %>% 
  full_join(est_M) %>% select(-id)

est_F <- dat_F %>% mutate(id = as.integer(id)) %>% select(
  -cht, -cwt, -cbmi, -age, -sqrt_age, -sex) %>% distinct(id, .keep_all = TRUE) %>% 
  full_join(est_F) %>% select(-id)

psme_est_MF <- est_M %>% bind_rows(est_F) %>% mutate(id = row_number())

rm(est_M, est_F, dat_M, dat_F, psme_pred_M, psme_pred_F, peak_M, peak_F,
   rebound_M, rebound_F, AUC, peak, xTraj, age.F, age.M, psme_bmi_F, psme_bmi_M,
   psme_ht_F, psme_ht_M, psme_wt_F, psme_wt_M)

save.image("psme_features.RData")
