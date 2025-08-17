
#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(tidyverse)
library(data.table)
library(growthcleanr)

#-------------------------------------------------------#
#### LOAD ANALYSIS DATASET & CLEAN WITH GROWTHCLEANR ####
#-------------------------------------------------------#

gusto_dat <- read.table(
  "dat/gusto_birthlength.csv", header = TRUE, 
  sep = ",", na.strings = c("","NA")
  )

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

#---------------------------------#
#### PREPARE ANALYSIS DATASETS ####
#---------------------------------#

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

#-----------------------------------#
#### PLOT OBSERVED GROWTH VALUES ####
#-----------------------------------#

load("gusto_ps_dat.RData")

dat_MF <- dat_M %>% select(cwt, cht, cbmi, age, sex)
dat_MF <- dat_F %>% select(cwt, cht, cbmi, age, sex) %>% bind_rows(dat_MF)

rm(dat_M, dat_F)

dat_MF$sex <- dplyr::recode(
  dat_MF$sex, '0' = "boys", '1' = "girls")

dat_MF <- dat_MF %>% pivot_longer(
  cols = c(cwt, cht, cbmi), 
  names_to = c("outs"), 
  values_to = "count")

dat_MF$outs <- dplyr::recode(
  dat_MF$outs, 
  cwt = "Weight (kg)",
  cht = "Height (cm)",
  cbmi = "BMI (kg/m²)"
)

dat_MF <- dat_MF %>% mutate(
  outsex = ifelse(
    sex == "boys" & outs == "Height (cm)", "a. Height (cm) - boys", ifelse(
      sex == "girls" & outs == "Height (cm)", "b. Height (cm) - girls", ifelse(
        sex == "boys" & outs == "Weight (kg)", "c. Weight (kg) - boys", ifelse(
          sex == "girls" & outs == "Weight (kg)", "d. Weight (kg) - girls", ifelse(
            sex == "boys" & outs == "BMI (kg/m²)", "e. BMI (kg/m²) - boys", "f. BMI (kg/m²) - girls"
          ))))))

graphics.off()
grDevices::cairo_pdf("res/fig_observed_values.pdf", height = 6, width = 5)

ggplot(data = dat_MF, aes(x = age, y = count)) + theme_classic() + 
  geom_point(size = 0.01)  + facet_wrap(
    outsex ~ ., strip.position="top", nrow = 3, scales = "free") + theme(
      legend.title = element_blank(), panel.grid.minor.x = element_blank(), 
      legend.position = "none", strip.text = element_text(hjust = 0, face="bold"),
      strip.background = element_blank()) + ylab('observed values') + xlab("age - years") +
  scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10))

dev.off()

rm(list = ls())
