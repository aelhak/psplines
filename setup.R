
#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(tidyverse); library(data.table); library(growthcleanr)

#-------------------------------------------------------#
#### LOAD ANALYSIS DATASET & CLEAN WITH GROWTHCLEANR ####
#-------------------------------------------------------#

gusto_dat <- read.table(
  "dat/gusto_birthlength.csv", header = TRUE, sep = ",", na.strings=c("","NA"))

gusto_dat <- read.table(
  "dat/gusto_long.csv", header = TRUE, sep = ",", na.strings=c("","NA")) %>%
  full_join(gusto_dat) %>% filter(age >= 0.038 & age <= 10) %>% group_by(id) %>% 
  mutate(id = cur_group_id()) %>% ungroup() %>% group_by(id) %>% 
  distinct(age, .keep_all = T) %>% ungroup() %>% mutate(
    subjid = id, agedays = age * 365, sex = ifelse(sex == 1, 0, 1)) %>% 
  select(-cbmi, -age) %>% filter(!is.na(agedays) & !is.na(cht) & !is.na(cwt)) %>% 
  pivot_longer(cols = c(cht, cwt), names_to = c("param"), values_to = "measurement") %>% 
  mutate(param = ifelse(param == "cwt", "WEIGHTKG", ifelse(
    param == "cht" & agedays <= 730, "LENGTHCM", "HEIGHTCM")))

summary(gusto_dat$id) # 1034

gusto_dat <- as.data.table(gusto_dat)
setkey(gusto_dat, subjid, param, agedays)

cleaned_dat <- gusto_dat[, gcr_result := cleangrowth(subjid, param, agedays, sex, measurement)]
gusto_dat <- cleaned_dat[gcr_result == "Include"]

gusto_dat <- gusto_dat %>% mutate(
  param = ifelse(param == "WEIGHTKG", "cwt", "cht")) %>% as.data.frame() %>% 
  pivot_wider(names_from = param, values_from = measurement) %>% mutate(
    age = agedays / 365, sqrt_age = sqrt(age), cbmi = cwt / ((cht/100)^2)) %>% 
  select(-gcr_result, -agedays, -subjid, -visit) %>% filter(
    !is.na(age) & !is.na(cbmi) & !is.na(cht) & !is.na(cwt))

#--------------------------------#
#### PREPARE ANALYSIS DATASET ####
#--------------------------------#

dat_M <- gusto_dat %>% filter(sex == 0)
dat_M <- subset(dat_M, id %in% with(rle(dat_M$id), values[lengths > 1]))
dat_M <- dat_M %>% group_by(id) %>% mutate(id = cur_group_id()) %>% ungroup()

dat_M %>% group_by(id) %>% count() %>% ungroup() %>% 
  summarise(n_distinct(id), median(n), min(n), max(n), IQR(n))

dat_M$id <- as.factor(dat_M$id)

dat_F <- gusto_dat %>% filter(sex == 1)
dat_F <- subset(dat_F, id %in% with(rle(dat_F$id), values[lengths > 1]))
dat_F <- dat_F %>% group_by(id) %>% mutate(id = cur_group_id()) %>% ungroup()

dat_F %>% group_by(id) %>% count() %>% ungroup() %>% 
  summarise(n_distinct(id), median(n), min(n), max(n), IQR(n))

dat_F$id <- as.factor(dat_F$id)

rm(cleaned_dat, gusto_dat)

#------------------------------#
#### SAVE ANALYSIS DATASETS ####
#------------------------------#

save.image("gusto_ps_dat.RData")
