#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(tidyverse)

#--------------------#
#### LOAD DATASET ####
#--------------------#

load("gusto_ps_dat.RData")

#----------------------#
#### LOAD FUNCTIONS ####
#----------------------#

source("src/double_penalty_helper_functions.R")

age.M <-  with(dat_M, seq.int(min(sqrt_age), max(sqrt_age), length.out = 100))
age.F <- with(dat_F, seq.int(min(sqrt_age), max(sqrt_age), length.out = 100))

#------------------------------------#
#### FIT P-SPLINE SNH (DP) MODELS ####
#------------------------------------#

nseg <- 10      # Number of segments
bdeg <- 3       # polynomial degree
pord.1_dp <- 1  # Penalty order SNH 1
pord.2_dp <- 2  # Penalty order SNH 2

smooth.pop <- smooth.ind <- list(
  nseg = nseg, 
  bdeg = bdeg, 
  pord1 = pord.1_dp,
  pord2 = pord.2_dp
)

SNH_ht_M <- try(SubjectSpecificCurves.DP(
  x = dat_M$sqrt_age, y = dat_M$cht, case = dat_M$id, 
  xgrid = age.M, smooth.pop = smooth.pop, smooth.ind = smooth.ind))  

SNH_wt_M <- try(SubjectSpecificCurves.DP(
  x = dat_M$sqrt_age, y = dat_M$cwt, case = dat_M$id, 
  xgrid = age.M, smooth.pop = smooth.pop, smooth.ind = smooth.ind))  

SNH_bmi_M <- try(SubjectSpecificCurves.DP(
  x = dat_M$sqrt_age, y = dat_M$cbmi, case = dat_M$id, 
  xgrid = age.M, smooth.pop = smooth.pop, smooth.ind = smooth.ind))  

SNH_ht_F <- try(SubjectSpecificCurves.DP(
  x = dat_F$sqrt_age, y = dat_F$cht, case = dat_F$id, 
  xgrid = age.F, smooth.pop = smooth.pop, smooth.ind = smooth.ind))  

SNH_wt_F <- try(SubjectSpecificCurves.DP(
  x = dat_F$sqrt_age, y = dat_F$cwt, case = dat_F$id, 
  xgrid = age.F, smooth.pop = smooth.pop, smooth.ind = smooth.ind))  

SNH_bmi_F <- try(SubjectSpecificCurves.DP(
  x = dat_F$sqrt_age, y = dat_F$cbmi, case = dat_F$id, 
  xgrid = age.F, smooth.pop = smooth.pop, smooth.ind = smooth.ind))  

save.image("gusto_psal_mods.RData")

#-----------------------------------#
#### FIT P-SPLINE MODELS IN PSME ####
#-----------------------------------#

library(psme)

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

save.image("gusto_psal_mods.RData")

rm(list = ls())

#-----------------------------#
#### GET MODEL PREDICTIONS ####
#-----------------------------#

load("gusto_psal_mods.RData")

# SNH

snh_pred_M <- as.data.frame(SNH_ht_M$fitted.ind.grid) %>% 
  mutate(age = xgrid^2, cht = fitted, id = case) %>% 
  select(id, age, cht)

snh_pred_M <- as.data.frame(SNH_wt_M$fitted.ind.grid) %>% 
  mutate(age = xgrid^2, cwt = fitted, id = case) %>% 
  select(id, age, cwt) %>% full_join(snh_pred_M)

snh_pred_M <- as.data.frame(SNH_bmi_M$fitted.ind.grid) %>% 
  mutate(age = xgrid^2, cbmi = fitted, id = case) %>% 
  select(id, age, cbmi) %>% full_join(snh_pred_M) %>% 
  mutate(sex = "Boys", method = "SNH")

snh_pred_F <- as.data.frame(SNH_ht_F$fitted.ind.grid) %>% 
  mutate(age = xgrid^2, cht = fitted, id = case) %>% 
  select(id, age, cht)

snh_pred_F <- as.data.frame(SNH_wt_F$fitted.ind.grid) %>% 
  mutate(age = xgrid^2, cwt = fitted, id = case) %>% 
  select(id, age, cwt) %>% full_join(snh_pred_F)

snh_pred_F <- as.data.frame(SNH_bmi_F$fitted.ind.grid) %>% 
  mutate(age = xgrid^2, cbmi = fitted, id = case) %>% 
  select(id, age, cbmi) %>% full_join(snh_pred_F) %>% 
  mutate(sex = "Girls", method = "SNH")

# STANDARD (PSME)

# utility function 1

xTraj <- function (model, new.x) {
  pop <- model$pcoef[[1]] + psme:::EvalSmooth(model$smooth[[1]], new.x)
  sub <- pop + psme:::EvalSmooth(model$smooth[[2]], new.x)
  list(pop = pop, sub = sub)
}

psme_pred_M <- as.data.frame(xTraj(psme_ht_M, age.M)$sub) %>% mutate(age = age.M^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cht")

psme_pred_M <- as.data.frame(xTraj(psme_wt_M, age.M)$sub) %>% mutate(age = age.M^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cwt") %>% 
  full_join(psme_pred_M)

psme_pred_M <- as.data.frame(xTraj(psme_bmi_M, age.M)$sub) %>% mutate(age = age.M^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cbmi") %>% 
  full_join(psme_pred_M) %>% mutate(sex = "Boys", method = "Standard")

psme_pred_M$id <- as.integer(sub("^\\D+", "", psme_pred_M$id))

psme_pred_F <- as.data.frame(xTraj(psme_ht_F, age.F)$sub) %>% mutate(age = age.F^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cht")

psme_pred_F <- as.data.frame(xTraj(psme_wt_F, age.F)$sub) %>% mutate(age = age.F^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cwt") %>% 
  full_join(psme_pred_F)

psme_pred_F <- as.data.frame(xTraj(psme_bmi_F, age.F)$sub) %>% mutate(age = age.F^2) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cbmi") %>% 
  full_join(psme_pred_F) %>% mutate(sex = "Girls", method = "Standard")

psme_pred_F$id <- as.integer(sub("^\\D+", "", psme_pred_F$id))

rm(SNH_bmi_M, SNH_bmi_F, SNH_ht_M, SNH_ht_F, SNH_wt_M, SNH_wt_F, dat_M, age.M, 
   age.F, psme_ht_M, psme_ht_F, psme_wt_M, psme_wt_F, dat_F, psme_bmi_M, psme_bmi_F, 
   xTraj, AUC, peak, smooth.ind, smooth.pop, bdeg, nseg, pord_dc, pord.1_dp, 
   pord.2_dp, Lambda.matrices, predict.bbase, pspline.comp, SubjectSpecificCurves.DP, msgEnv)

pred_M <- snh_pred_M %>% mutate(id = as.integer(id)) %>% full_join(psme_pred_M)
pred_F <- snh_pred_F %>% mutate(id = as.integer(id)) %>% full_join(psme_pred_F)

rm(dch_pred_M, dch_pred_F, snh_pred_M, snh_pred_F, psme_pred_M, psme_pred_F, msgEnv)

save.image("psal_preds.RData")

#---------------------------#
#### GET GROWTH FEATURES ####
#---------------------------#

load("psal_preds.RData")

#### INFANT PEAK BMI - BOYS ####

peak_M <- pred_M %>% filter(method == "Standard") %>% select(id, age, cbmi) %>% filter(age <= 1.5)
peak_M <- split(peak_M, peak_M$id)

est_M <- map_df(peak_M, ~ sitar::getPeak(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = unique(pred_M$id)) %>% relocate(id) %>% rename(
    PeakBMIAge_psme = x, PeakBMI_psme = y) %>% mutate(sex = "Boys")

peak_M <- pred_M %>% filter(method == "SNH") %>% select(id, age, cbmi) %>% filter(age <= 1.5)
peak_M <- split(peak_M, peak_M$id)

est_M <- map_df(peak_M, ~ sitar::getPeak(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = unique(pred_M$id)) %>% relocate(id) %>% rename(
    PeakBMIAge_dp = x, PeakBMI_dp = y) %>% full_join(est_M)

rm(peak_M)

#### CHILD REBOUND BMI - BOYS ####

rebound_M <- pred_M %>% filter(method == "Standard") %>% select(id, age, cbmi) %>% filter(age > 2.5)
rebound_M <- split(rebound_M, rebound_M$id)

est_M <- map_df(rebound_M, ~ sitar::getTrough(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = est_M$id) %>% relocate(id) %>% rename(
    ReboundBMIAge_psme = x, ReboundBMI_psme = y) %>% full_join(est_M)

rebound_M <- pred_M %>% filter(method == "SNH") %>% select(id, age, cbmi) %>% filter(age > 2.5)
rebound_M <- split(rebound_M, rebound_M$id)

est_M <- map_df(rebound_M, ~ sitar::getTrough(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = est_M$id) %>% relocate(id) %>% rename(
    ReboundBMIAge_dp = x, ReboundBMI_dp = y) %>% full_join(est_M)

rm(rebound_M)

#### INFANT PEAK HEIGHT VELOCITY - BOYS ####

est_M <- pred_M %>% filter(method == "Standard") %>% select(id, age, cht) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PHV_psme = max(dY)) %>% full_join(est_M)

est_M <- pred_M %>% filter(method == "SNH") %>% select(id, age, cht) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PHV_dp = max(dY)) %>% full_join(est_M)

#### INFANT PEAK WEIGHT VELOCITY - BOYS ####

est_M <- pred_M %>% filter(method == "Standard") %>% select(id, age, cwt) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PWV_psme = max(dY)) %>% full_join(est_M)

est_M <- pred_M %>% filter(method == "SNH") %>% select(id, age, cwt) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PWV_dp = max(dY)) %>% full_join(est_M)

rm(pred_M)

#### INFANT PEAK BMI - GIRLS ####

peak_F <- pred_F %>% filter(method == "Standard") %>% select(id, age, cbmi) %>% filter(age <= 1.5)
peak_F <- split(peak_F, peak_F$id)

est_F <- map_df(peak_F, ~ sitar::getPeak(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = unique(pred_F$id)) %>% relocate(id) %>% rename(
    PeakBMIAge_psme = x, PeakBMI_psme = y) %>% mutate(sex = "Girls")

peak_F <- pred_F %>% filter(method == "SNH") %>% select(id, age, cbmi) %>% filter(age <= 1.5)
peak_F <- split(peak_F, peak_F$id)

est_F <- map_df(peak_F, ~ sitar::getPeak(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = unique(pred_F$id)) %>% relocate(id) %>% rename(
    PeakBMIAge_dp = x, PeakBMI_dp = y) %>% full_join(est_F)

rm(peak_F)

#### CHILD REBOUND BMI - GIRLS ####

rebound_F <- pred_F %>% filter(method == "Standard") %>% select(id, age, cbmi) %>% filter(age > 2.5)
rebound_F <- split(rebound_F, rebound_F$id)

est_F <- map_df(rebound_F, ~ sitar::getTrough(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = est_F$id) %>% relocate(id) %>% rename(
    ReboundBMIAge_psme = x, ReboundBMI_psme = y) %>% full_join(est_F)

rebound_F <- pred_F %>% filter(method == "SNH") %>% select(id, age, cbmi) %>% filter(age > 2.5)
rebound_F <- split(rebound_F, rebound_F$id)

est_F <- map_df(rebound_F, ~ sitar::getTrough(x = .x$age, y = .x$cbmi)) %>% 
  mutate(id = est_F$id) %>% relocate(id) %>% rename(
    ReboundBMIAge_dp = x, ReboundBMI_dp = y) %>% full_join(est_F)

rm(rebound_F)

#### INFANT PEAK HEIGHT VELOCITY - GIRLS ####

est_F <- pred_F %>% filter(method == "Standard") %>% select(id, age, cht) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PHV_psme = max(dY)) %>% full_join(est_F)

est_F <- pred_F %>% filter(method == "SNH") %>% select(id, age, cht) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PHV_dp = max(dY)) %>% full_join(est_F)

#### INFANT PEAK WEIGHT VELOCITY - GIRLS ####

est_F <- pred_F %>% filter(method == "Standard") %>% select(id, age, cwt) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PWV_psme = max(dY)) %>% full_join(est_F)

est_F <- pred_F %>% filter(method == "SNH") %>% select(id, age, cwt) %>% filter(age < 1) %>% 
  group_by(id) %>% reframe(dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2))) %>% 
  ungroup() %>% group_by(id) %>% summarise(PWV_dp = max(dY)) %>% full_join(est_F)

rm(pred_F)

## SAVE

est_M <- est_M %>% select(-id)
est_MF <- est_F %>% select(-id) %>% bind_rows(est_M)
rm(est_M, est_F)

save.image("psal_features.RData")
