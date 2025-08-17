
#------------------------#
#### UTILITY FUNCTION ####
#------------------------#

xTraj <- function (model, new.x) {
  pop <- model$pcoef[[1]] + psme:::EvalSmooth(model$smooth[[1]], new.x)
  sub <- pop + psme:::EvalSmooth(model$smooth[[2]], new.x)
  list(pop = pop, sub = sub)
}

#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(psme)
library(sitar)
library(tidyverse)
library(growthstandards)

#-----------------------------------#
#### FIT P-SPLINE MODELS IN PSME ####
#-----------------------------------#

load("gusto_ps_dat.RData")

fit_psme <- function(response, data) {
  psme(
    as.formula(paste0(
      response,
      " ~ s(sqrt_age, bs = 'ps', k = 10, m = c(2, 2)) + ",
      "s(sqrt_age, id, bs = 'fs', xt = 'ps', k = 10, m = c(2, 1))"
    )),
    data = data
  )
}

vars <- c(ht = "cht", wt = "cwt", bmi = "cbmi")
datasets <- list(M = dat_M, F = dat_F)

for (sex in names(datasets)) {
  for (vname in names(vars)) {
    obj_name <- paste0("psme_", vname, "_", sex)
    assign(obj_name,
           fit_psme(vars[[vname]], datasets[[sex]]),
           envir = .GlobalEnv)
  }
}

rm(fit_psme, vars, datasets)
save.image("gusto_psme_mods.RData")
rm(list = ls())

#----------------------------------#
#### GET PSME MODEL PREDICTIONS ####
#----------------------------------#

load("gusto_psme_mods.RData")

make_preds <- function(models, age_seq) {
  vars <- c("cht", "cwt", "cbmi")
  out <- NULL
  for (i in seq_along(models)) {
    tmp <- as.data.frame(xTraj(models[[i]], age_seq)$sub) %>%
      mutate(age = age_seq^2) %>% pivot_longer(
        cols = starts_with("V"), names_to = "id", values_to = vars[i])
    out <- if (is.null(out)) tmp else full_join(tmp, out)
  }
  out$id <- as.integer(sub("^\\D+", "", out$id))
  out
}

age.M <- with(dat_M, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))
psme_pred_M <- make_preds(list(psme_ht_M, psme_wt_M, psme_bmi_M), age.M)

age.F <- with(dat_F, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))
psme_pred_F <- make_preds(list(psme_ht_F, psme_wt_F, psme_bmi_F), age.F)

rm(dat_M, dat_F, age.M, age.F, psme_ht_M, psme_ht_F, psme_wt_M, psme_wt_F, 
   psme_bmi_M, psme_bmi_F, xTraj, make_preds)

save.image("gusto_psme_preds.RData")
rm(list = ls())

#------------------------------------------------------#
#### GET GROWTH FEATURES: PV, PEAK BMI, REBOUND BMI ####
#------------------------------------------------------#

load("gusto_ps_dat.RData")
load("gusto_psme_preds.RData")

get_peak <- function(df, age_limit) {
  tmp <- df %>% filter(age <= age_limit) %>% split(.$id)
  map_df(tmp, ~ sitar::getPeak(x = .x$age, y = .x$cbmi))
}

get_trough <- function(df, age_min) {
  tmp <- df %>% filter(age >= age_min) %>% split(.$id)
  map_df(tmp, ~ sitar::getTrough(x = .x$age, y = .x$cbmi))
}

#### BOYS ####

est_M <- get_peak(psme_pred_M, 1.5) %>%
  mutate(id = unique(psme_pred_M$id), sex = "Boys") %>%
  rename(PeakBMIAge = x, PeakBMI = y)

est_M <- get_trough(psme_pred_M, 2.5) %>%
  mutate(id = unique(psme_pred_M$id)) %>%
  rename(ReboundBMIAge = x, ReboundBMI = y) %>% 
  full_join(est_M)

est_M <- psme_pred_M %>% select(id, age, cht, cwt) %>% 
  filter(age < 1) %>% group_by(id) %>% summarise(
    PHV = max(diff(cht) / diff(age)), PWV = max(diff(cwt) / diff(age))) %>% 
  full_join(est_M)

#### GIRLS ####

est_F <- get_peak(psme_pred_F, 1.5) %>%
  mutate(id = unique(psme_pred_F$id), sex = "Girls") %>%
  rename(PeakBMIAge = x, PeakBMI = y)

est_F <- get_trough(psme_pred_F, 2.5) %>%
  mutate(id = unique(psme_pred_F$id)) %>%
  rename(ReboundBMIAge = x, ReboundBMI = y) %>% 
  full_join(est_F)

est_F <- psme_pred_F %>% select(id, age, cht, cwt) %>% 
  filter(age < 1) %>% group_by(id) %>% summarise(
    PHV = max(diff(cht) / diff(age)), PWV = max(diff(cwt) / diff(age))) %>% 
  full_join(est_F)

rm(get_peak, get_trough, psme_pred_M, psme_pred_F)

#-----------------------------------------------------#
#### GET HEIGHT & WEIGHT DIFFERENCES VS WHO MEDIAN ####
#-----------------------------------------------------#

load("gusto_psme_mods.RData")

age_vec <- sqrt(c(0.0833334, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5))

calc_had_wad <- function(model_ht, model_wt, ages_sqrt, had_ref, wad_ref) {
  pred_ht <- as.data.frame(xTraj(model_ht, ages_sqrt)$sub) %>%
    mutate(age = ages_sqrt^2) %>%
    pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cht")
  pred_wt <- as.data.frame(xTraj(model_wt, ages_sqrt)$sub) %>%
    mutate(age = ages_sqrt^2) %>%
    pivot_longer(cols = starts_with("V"), names_to = "id", values_to = "cwt")
  pred <- full_join(pred_ht, pred_wt)
  pred$id <- as.integer(sub("^\\D+", "", pred$id))
  pred$age <- as.integer(factor(pred$age, levels = sort(unique(pred$age))))
  for (i in seq_along(had_ref)) {
    pred[[paste0("HAD_", names(had_ref)[i])]] <- ifelse(pred$age == i, pred$cht - had_ref[i], NA)
    pred[[paste0("WAD_", names(wad_ref)[i])]] <- ifelse(pred$age == i, pred$cwt - wad_ref[i], NA)
  }
  pred %>% select(-cwt, -cht, -age)
}

#### BOYS ####

had_ref_M <- c(`1m`=54.7, `6m`=67.6, `12m`=75.7, `18m`=82.3, `24m`=87.8, `30m`=91.9,
               `36m`=96.1, `42m`=99.9, `48m`=103.3, `54m`=106.7, `60m`=110)
wad_ref_M <- c(`1m`=4.5, `6m`=7.9, `12m`=9.6, `18m`=10.9, `24m`=12.2, `30m`=13.3,
               `36m`=14.3, `42m`=15.3, `48m`=16.3, `54m`=17.3, `60m`=18.3)

whad_M <- calc_had_wad(psme_ht_M, psme_wt_M, age_vec, had_ref_M, wad_ref_M)

for (col in grep("^HAD_|^WAD_", names(whad_M), value = TRUE)) {
  est_M <- full_join(est_M, whad_M %>% select(id, !!col) %>% drop_na())
}

rm(had_ref_M, wad_ref_M, whad_M, psme_ht_M, psme_wt_M, psme_bmi_M)

#### GIRLS ####

had_ref_F <- c(`1m`=53.7, `6m`=65.7, `12m`=74.0, `18m`=80.7, `24m`=86.4, `30m`=90.7,
               `36m`=95.1, `42m`=99.0, `48m`=102.7, `54m`=106.2, `60m`=109.4)
wad_ref_F <- c(`1m`=4.2, `6m`=7.3, `12m`=8.9, `18m`=10.2, `24m`=11.5, `30m`=12.7,
               `36m`=13.9, `42m`=15.0, `48m`=16.1, `54m`=17.2, `60m`=18.2)

whad_F <- calc_had_wad(psme_ht_F, psme_wt_F, age_vec, had_ref_F, wad_ref_F)

for (col in grep("^HAD_|^WAD_", names(whad_F), value = TRUE)) {
  est_F <- full_join(est_F, whad_F %>% select(id, !!col) %>% drop_na())
}

rm(had_ref_F, wad_ref_F, whad_F, psme_ht_F, psme_wt_F, psme_bmi_F,
   age_vec, col, calc_had_wad, xTraj)

#-----------------------#
#### COMBINE & MERGE ####
#-----------------------#

est_M <- dat_M %>% mutate(id = as.integer(id)) %>% select(
  -cht, -cwt, -cbmi, -age, -sqrt_age, -sex) %>% distinct(
    id, .keep_all = TRUE) %>% full_join(est_M) %>% select(-id)

est_F <- dat_F %>% mutate(id = as.integer(id)) %>% select(
  -cht, -cwt, -cbmi, -age, -sqrt_age, -sex) %>% distinct(
    id, .keep_all = TRUE) %>% full_join(est_F) %>% select(-id)

psme_est_MF <- est_M %>% bind_rows(est_F) %>% mutate(id = row_number())

rm(est_M, est_F, dat_M, dat_F)

save.image("psme_features.RData")

summary(psme_est_MF)
psme_est_MF %>% filter(is.na(PeakBMI) | is.na(ReboundBMI)) %>% summarise(N = n_distinct(id)) # N=153
psme_est_MF %>% filter(is.na(PeakBMI) & !is.na(ReboundBMI)) %>% summarise(N = n_distinct(id)) # N=23
psme_est_MF %>% filter(!is.na(PeakBMI) & is.na(ReboundBMI)) %>% summarise(N = n_distinct(id)) # N=121
psme_est_MF %>% filter(is.na(PeakBMI) & is.na(ReboundBMI)) %>% summarise(N = n_distinct(id)) # N=9
