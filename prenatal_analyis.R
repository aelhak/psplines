
#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(tidyverse)
library(scales)
library(broom)

#-----------------------------#
#### BWP_CAT RESULTS EUCCN ####
#-----------------------------#

load("psme_features.RData")

psme_est_MF <- psme_est_MF %>% mutate_at(
  c("m_age", "m_wt", "m_ht", "GA", "BWT", "BLT"), .funs = list(z = ~ scale(.))) %>% 
  as.data.frame()

exp <- c("m_age_z", "m_wt_z", "m_ht_z", "GA_z", "BWT_z", "BLT_z")
out <- c("PHV", "PWV", "PeakBMI", "PeakBMIAge * 12", "ReboundBMI", "ReboundBMIAge")

(psme_assoc_res <- expand.grid(out, exp) %>% group_by(Var1) %>% 
  rowwise() %>% summarise(frm = paste0(Var1, " ~ sex + ", Var2)) %>% 
  group_by(model_id = row_number(), frm) %>% do(cbind(tidy(
    lm(.$frm, data = psme_est_MF)))) %>% mutate(
      lci = estimate - (1.96 * std.error), 
      uci = estimate + (1.96 * std.error)) %>% 
  filter(term != "sexGirls" & term != "(Intercept)") %>% 
  select(-std.error, -std.error, -statistic) %>% as.data.frame())

psme_assoc_res$term <- dplyr::recode(
  psme_assoc_res$term, 
  "m_bmi_z" = "i. Maternal BMI", 
  "f_bmi_z" = "ii. Paternal BMI", 
  "m_age_z" = "iii. Maternal age", 
  "f_age_z" = "iv. Paternal age", 
  "GA_z" = "v. Gestational age",
  "BW_z" = "vi. Birth weight"
)

psme_assoc_res$outs <- rep(c(
  "a. Peak height velocity (cm/yr)",
  "b. Peak weight velocity (kg/yr)",
  "c. Height growth: 0-2y (cm)",
  "d. Weight growth: 0-2y (kg)",
  "e. Height growth: 2-10y (cm)",
  "f. Weight growth: 2-10y (kg)",
  "g. Peak BMI (kg/m²)",
  "h. Age at peak BMI (months)",
  "i. Rebound BMI (kg/m²)",
  "j. Age at rebound BMI (years)"), 6)


rm(psme_est_MF, exp, out)

save.image("psme_assoc_res.RData")

rm(list = ls())
