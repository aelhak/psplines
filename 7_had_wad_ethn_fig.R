
#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(ggeffects)
library(tidyverse) 
library(splines)
library(scales)
library(lme4)

#------------------------------#
#### LOAD + PREPARE DATASET ####
#------------------------------#

load("psme_features.RData")

psme_est_MF$m_ethn <- dplyr::recode(
  psme_est_MF$m_ethn, 
  "1" = "Chinese", 
  "2" = "Malay", 
  "3" = "Indian"
)

psme_est_MF$m_ethn <- as.factor(psme_est_MF$m_ethn)

table(psme_est_MF$m_ethn)

ethn_HAD_dat <- psme_est_MF %>% select(
  id, sex, m_ethn, HAD_1m, HAD_6m, HAD_12m, HAD_18m, HAD_24m, 
  HAD_30m, HAD_36m, HAD_42m, HAD_48m, HAD_54m, HAD_60m) %>% pivot_longer(
    cols = starts_with("HAD"), names_to = "age", values_to = "had") %>% 
  mutate(t = ifelse(age == "HAD_1m", 1, ifelse(
    age == "HAD_6m", 2, ifelse(
      age == "HAD_12m", 3, ifelse(
        age == "HAD_18m", 4, ifelse(
          age == "HAD_24m", 5, ifelse(
            age == "HAD_30m", 6, ifelse(
              age == "HAD_36m", 7, ifelse(
                age == "HAD_42m", 8, ifelse(
                  age == "HAD_48m", 9, ifelse(
                    age == "HAD_54m", 10, 11
                  )))))))))))

ethn_WAD_dat <- psme_est_MF %>% select(
  id, sex, m_ethn, WAD_1m, WAD_6m, WAD_12m, WAD_18m, WAD_24m, 
  WAD_30m, WAD_36m, WAD_42m, WAD_48m, WAD_54m, WAD_60m) %>% pivot_longer(
    cols = starts_with("WAD"), names_to = "age", values_to = "wad") %>% 
  mutate(t = ifelse(age == "WAD_1m", 1, ifelse(
    age == "WAD_6m", 2, ifelse(
      age == "WAD_12m", 3, ifelse(
        age == "WAD_18m", 4, ifelse(
          age == "WAD_24m", 5, ifelse(
            age == "WAD_30m", 6, ifelse(
              age == "WAD_36m", 7, ifelse(
                age == "WAD_42m", 8, ifelse(
                  age == "WAD_48m", 9, ifelse(
                    age == "WAD_54m", 10, 11
                  )))))))))))

rm(psme_est_MF)

#-----------------#
#### ME MODELS ####
#-----------------#

#### HAD

HAD_mod <- lmer(
  had ~ (ns(t, df = 3) * m_ethn) + sex + ( 1 | id ), data = ethn_HAD_dat
)

HAD_mod_P <- ggemmeans(HAD_mod, terms = c('m_ethn [all]', 't [all]'))
HAD_mod_P$out <- "a. Height difference, cm"

#### WAD

WAD_mod <- lmer(
  wad ~ (ns(t, df = 3) * m_ethn) + sex + ( 1 | id ), data = ethn_WAD_dat
  )

WAD_mod_P <- ggemmeans(WAD_mod, terms = c('m_ethn [all]', 't [all]'))
WAD_mod_P$out <- "b. Weight difference, kg"

#------------#
#### PLOT ####
#------------#

mods_P <- bind_rows(HAD_mod_P, WAD_mod_P)

mods_P$group <- dplyr::recode(
  mods_P$group, 
  "1" = "1", 
  "2" = "6", 
  "3" = "12", 
  "4" = "18", 
  "5" = "24", 
  "6" = "30", 
  "7" = "36", 
  "8" = "42", 
  "9" = "48",
  "10" = "54",
  "11" = "60"
  
)

mods_P$group <- factor(
  mods_P$group, levels=c(
    "1",
    "6",
    "12",
    "18",
    "24",
    "30",
    "36",
    "42",
    "48",
    "54",
    "60"
  ))

graphics.off()
grDevices::cairo_pdf("res/fig_WHAD_ME.pdf", height = 3.8, width = 7.5)

ggplot(mods_P, aes(x = group, y = predicted, group = x, col = x, fill = x)) + theme_classic() + 
  geom_line(linewidth = 1.5) + geom_ribbon(aes(ymin = conf.low, ymax = conf.high), colour = NA, alpha = 0.15) +
  scale_color_brewer(palette = "Set2") +  ylab("Mean difference vs. WHO standards ") + 
  xlab("age, months") + facet_wrap(. ~ out, strip.position = "top", ncol = 2, scales = "free") + theme(
    legend.position = "bottom", legend.title = element_blank(), strip.background = element_blank(),
    strip.text = element_text(face="bold", hjust = 0), axis.ticks.y = element_blank()) + guides(
      override.aes = list(size = 1))

dev.off()
