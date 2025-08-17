
#-----------------------------#
#### PSME UTILITY FUNCTION ####
#-----------------------------#

xTraj <- function (model, new.x) {
  pop <- model$pcoef[[1]] + psme:::EvalSmooth(model$smooth[[1]], new.x)
  sub <- pop + psme:::EvalSmooth(model$smooth[[2]], new.x)
  list(pop = pop, sub = sub)
}

#---------------------------#
#### DP HELPER FUNCTIONS ####
#---------------------------#

source("src/6_1_double_penalty_helper_functions.R")

#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(tidyverse)
library(psme)

#--------------------#
#### LOAD DATASET ####
#--------------------#

load("gusto_ps_dat.RData")

age.M <-  with(dat_M, seq.int(min(sqrt_age), max(sqrt_age), length.out = 100))
age.F <- with(dat_F, seq.int(min(sqrt_age), max(sqrt_age), length.out = 100))

#--------------------------------#
#### FIT P-SPLINE (DP) MODELS ####
#--------------------------------#

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

datasets  <- list(M = dat_M, F = dat_F)
age_grids <- list(M = age.M, F = age.F)
vars      <- c(ht = "cht", wt = "cwt", bmi = "cbmi")

dp_models <- function(data, yvar, age_grid) {
  try(
    SubjectSpecificCurves.DP(
      x = data$sqrt_age,
      y = data[[yvar]],
      case = data$id,
      xgrid = age_grid,
      smooth.pop = smooth.pop,
      smooth.ind = smooth.ind
    )
  )
}

for (sex in names(datasets)) {
  for (vname in names(vars)) {
    obj_name <- paste0("DP_", vname, "_", sex)
    assign(obj_name,
           dp_models(datasets[[sex]], vars[[vname]], age_grids[[sex]]),
           envir = .GlobalEnv)
  }
}

rm(msgEnv, smooth.ind, smooth.pop, bdeg, nseg, obj_name, pord.1_dp, pord.2_dp,
   sex, vname, dp_models, Lambda.matrices, predict.bbase, pspline.comp,
   SubjectSpecificCurves.DP)

save.image("psdp_mods.RData")

#-----------------------------------#
#### FIT P-SPLINE MODELS IN PSME ####
#-----------------------------------#

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

save.image("psdp_mods.RData")

#-----------------------------#
#### GET MODEL PREDICTIONS ####
#-----------------------------#

load("psdp_mods.RData")

# DP models

sex <- list(
  M = list(suffix = "M", label = "Boys"),
  F = list(suffix = "F", label = "Girls"))

dp_pred_MF <- lapply(sex, function(params) {
  
  measurements <- c("ht", "wt", "bmi")
  
  df_list <- map(measurements, ~{
    object_name <- paste0("DP_", .x, "_", params$suffix)
    new_col_name <- paste0("c", .x)
    
    as.data.frame(get(object_name)$fitted.ind.grid) %>%
      transmute(id = case, age = xgrid^2, !!new_col_name := fitted)
  })
  
  reduce(df_list, full_join, by = c("id", "age")) %>%
    mutate(sex = params$label, method = "DP")
})

dp_pred_M <- dp_pred_MF$M
dp_pred_F <- dp_pred_MF$F

dp_pred_M$id <- as.integer(dp_pred_M$id)
dp_pred_F$id <- as.integer(dp_pred_F$id)

rm(sex, dp_pred_MF, DP_bmi_M, DP_bmi_F, DP_ht_M, DP_ht_F, DP_wt_M, DP_wt_F)

# PSME models

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

psme_pred_M <- make_preds(list(psme_ht_M, psme_wt_M, psme_bmi_M), age.M)
psme_pred_F <- make_preds(list(psme_ht_F, psme_wt_F, psme_bmi_F), age.F)

rm(dat_M, age.M, age.F, make_preds, psme_ht_M, psme_ht_F, psme_wt_M, 
   psme_wt_F, dat_F, psme_bmi_M, psme_bmi_F, xTraj, age_grids)

save.image("psdp_preds.RData")

#---------------------------#
#### GET GROWTH FEATURES ####
#---------------------------#

load("psdp_preds.RData")

get_peak <- function(df, age_limit) {
  tmp <- df %>% filter(age <= age_limit) %>% split(.$id)
  map_df(tmp, ~ sitar::getPeak(x = .x$age, y = .x$cbmi))
}

get_trough <- function(df, age_min) {
  tmp <- df %>% filter(age >= age_min) %>% split(.$id)
  map_df(tmp, ~ sitar::getTrough(x = .x$age, y = .x$cbmi))
}

# DP

dp_est_M <- get_peak(dp_pred_M, 1.5) %>% mutate(
  id = unique(dp_pred_M$id), sex = "Boys", method = "double penalty") %>%
  rename(PeakBMIAge = x, PeakBMI = y)

dp_est_M <- get_trough(dp_pred_M, 2.5) %>% mutate(
  id = unique(dp_pred_M$id)) %>% rename(ReboundBMIAge = x, ReboundBMI = y) %>% 
  full_join(dp_est_M)

dp_est_M <- dp_pred_M %>% select(id, age, cht, cwt) %>% filter(age < 1) %>% group_by(id) %>% 
  summarise(PHV = max(diff(cht) / diff(age)), PWV = max(diff(cwt) / diff(age))) %>% 
  full_join(dp_est_M)

# PSME

psme_est_M <- get_peak(psme_pred_M, 1.5) %>% mutate(
  id = unique(psme_pred_M$id), sex = "Boys", method = "psme penalty") %>%
  rename(PeakBMIAge = x, PeakBMI = y)

psme_est_M <- get_trough(psme_pred_M, 2.5) %>% mutate(
  id = unique(psme_pred_M$id)) %>% rename(ReboundBMIAge = x, ReboundBMI = y) %>% 
  full_join(psme_est_M)

psme_est_M <- psme_pred_M %>% select(id, age, cht, cwt) %>% filter(age < 1) %>% group_by(id) %>% 
  summarise(PHV = max(diff(cht) / diff(age)), PWV = max(diff(cwt) / diff(age))) %>% 
  full_join(psme_est_M)

# JOIN

psme_dp_est_M <- full_join(dp_est_M, psme_est_M)

rm(get_peak, get_trough, dp_est_M, psme_est_M, dp_pred_M, psme_pred_M, dp_pred_F, psme_pred_F)

save.image("psdp_features.RData")

#--------------------------------#
#### PLOT PSME VS. DP FEATRES ####
#--------------------------------#

load("psdp_features.RData")

print(psme_dp_est_M %>% group_by(method) %>% summarise(
  across(c(PHV, PWV, PeakBMI, ReboundBMI, PeakBMIAge, ReboundBMIAge), list(
    mean = ~ mean(.x, na.rm = T), median = ~ median(.x, na.rm = T),
    min = ~ min(.x, na.rm = T), max = ~ max(.x, na.rm = T)), 
    .names = "{.col}_{.fn}")) %>% pivot_longer(
      -method, names_to = c("variable", "stat"), names_sep = "_") %>%
    pivot_wider(names_from = method, values_from = value), n= 24)

psme_dp_est_M <- psme_dp_est_M %>%
  pivot_wider(id_cols = id, names_from = method, values_from = c(
    PHV, PWV, ReboundBMIAge, ReboundBMI, PeakBMIAge, PeakBMI))

graphics.off()
grDevices::cairo_pdf("res/fig_psdp_cor.pdf", height = 7.5, width = 6.5)
par(mgp=c(2,1,0), mfrow = c(3, 2), mar = c(4, 4.6, 1.6, 1), oma = c(0, 0, 0, 0))

with(psme_dp_est_M, plot(`PHV_psme penalty`, `PHV_double penalty`, xlab = "psme penalty", ylab = "double penalty", xlim=c(60, 97), ylim=c(60, 97)))
title(paste("a. Peak height velocity (cm/yr), r =", round(with(psme_dp_est_M, cor(`PHV_psme penalty`, `PHV_double penalty`)), 2)), adj = 0, line = 0.6)
abline(a = 0, b = 1, col = "red")

with(psme_dp_est_M, plot(`PWV_psme penalty`, `PWV_double penalty`, xlab = "psme penalty", ylab = "double penalty", xlim=c(9, 36), ylim=c(9, 36)))
title(paste("b. Peak weight velocity (kg/yr), r =", round(with(psme_dp_est_M, cor(`PWV_psme penalty`, `PWV_double penalty`)), 2)), adj = 0, line = 0.6)
abline(a = 0, b = 1, col = "red")

with(psme_dp_est_M, plot(`PeakBMI_psme penalty`, `PeakBMI_double penalty`, xlab = "psme penalty", ylab = "double penalty", xlim=c(12, 24), ylim=c(12, 24)))
title(paste("c. Peak BMI (kg/m²), r =", round(with(psme_dp_est_M, cor(`PeakBMI_psme penalty`, `PeakBMI_double penalty`, use='complete.obs')), 3)), adj = 0, line = 0.6)
abline(a = 0, b = 1, col = "red")

with(psme_dp_est_M, plot(`ReboundBMI_psme penalty`, `ReboundBMI_double penalty`, xlab = "psme penalty", ylab = "double penalty", xlim=c(10, 26), ylim=c(10, 26)))
title(paste("d. Rebound BMI (kg/m²), r =", round(with(psme_dp_est_M, cor(`ReboundBMI_psme penalty`, `ReboundBMI_double penalty`, use='complete.obs')), 3)), adj = 0, line = 0.6)
abline(a = 0, b = 1, col = "red")

with(psme_dp_est_M, plot(`PeakBMIAge_psme penalty`*12, `PeakBMIAge_double penalty`*12, xlab = "psme penalty", ylab = "double penalty", xlim=c(2, 18), ylim=c(2, 18)))
title(paste("e. Age at peak BMI (month), r =", round(with(psme_dp_est_M, cor(`PeakBMIAge_psme penalty`*12, `PeakBMIAge_double penalty`*12, use='complete.obs')), 2)), adj = 0, line = 0.6)
abline(a = 0, b = 1, col = "red")

with(psme_dp_est_M, plot(`ReboundBMIAge_psme penalty`, `ReboundBMIAge_double penalty`, xlab = "psme penalty", ylab = "double penalty", xlim=c(2, 10), ylim=c(2, 10)))
title(paste("f. Age at rebound BMI (year), r =", round(with(psme_dp_est_M, cor(`ReboundBMIAge_psme penalty`, `ReboundBMIAge_double penalty`, use='complete.obs')), 2)), adj = 0, line = 0.6)
abline(a = 0, b = 1, col = "red")

dev.off()

