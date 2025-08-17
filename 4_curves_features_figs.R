
# utility function 1

xTraj <- function (model, new.x) {
  pop <- model$pcoef[[1]] + psme:::EvalSmooth(model$smooth[[1]], new.x)
  sub <- pop + psme:::EvalSmooth(model$smooth[[2]], new.x)
  list(pop = pop, sub = sub)
}

#---------------------#
#### LOAD PACKAGES ####
#---------------------#

library(RColorBrewer)
library(tidyverse)
library(patchwork)
library(corrplot)
library(cowplot)
library(stringr)
library(ggh4x)
library(ggpubr)
library(scales)
library(sitar)
library(broom)
library(psme)

#-------------------------------#
#### PLOT PSME GROWTH CURVES ####
#-------------------------------#

load("gusto_psme_preds.RData")
load("gusto_psme_mods.RData")

age.M <- with(dat_M, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))
age.F <- with(dat_F, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))

# F3: DISTANCE

graphics.off()
grDevices::cairo_pdf("res/fig_psme_dist_curves.pdf", height = 6.5, width = 5.5)
par(mfrow = c(3, 2), mar = c(4, 4.6, 1.6, 1), oma = c(0, 0, 0, 0))

mplot(x = age, y = cht, id = id, data = psme_pred_M, col = id, las = 1, ylim = c(35, 165), xlab = 'age', ylab = 'cm')
lines(age.M^2, xTraj(psme_ht_M, age.M)$pop, type = "l", col = 'black', lwd = 3.5)
title('a. Height - boys', adj = 0, line = 0.6)

mplot(x = age, y = cht, id = id, data = psme_pred_F, col = id, las = 1, ylim = c(35, 165), xlab = 'age', ylab = 'cm')
lines(age.F^2, xTraj(psme_ht_F, age.F)$pop, type = "l", col = 'black', lwd = 3.5)
title('b. Height - girls', adj = 0, line = 0.6)

mplot(x = age, y = cwt, id = id, data = psme_pred_M, col = id, las = 1, ylim = c(0, 80), xlab = 'age', ylab = 'kg')
lines(age.M^2, xTraj(psme_wt_M, age.M)$pop, type = "l", col = 'black', lwd = 3.5)
title('c. Weight - boys', adj = 0, line = 0.6)

mplot(x = age, y = cwt, id = id, data = psme_pred_F, col = id, las = 1, ylim = c(0, 80), xlab = 'age', ylab = 'kg')
lines(age.F^2, xTraj(psme_wt_F, age.F)$pop, type = "l", col = 'black', lwd = 3.5)
title('d. Weight - girls', adj = 0, line = 0.6)

mplot(x = age, y = cbmi, id = id, data = psme_pred_M, col = id, las = 1, ylim = c(8, 35), xlab = 'age', ylab = 'kg/m²')
lines(age.M^2, xTraj(psme_bmi_M, age.M)$pop, type = "l", col = 'black', lwd = 3.5)
title('e. BMI - boys', adj = 0, line = 0.6)

mplot(x = age, y = cbmi, id = id, data = psme_pred_F, col = id, las = 1, ylim = c(8, 35), xlab = 'age', ylab = 'kg/m²')
lines(age.F^2, xTraj(psme_bmi_F, age.F)$pop, type = "l", col = 'black', lwd = 3.5)
title('f. BMI - girls', adj = 0, line = 0.6)

dev.off()

#---------------------------------#
#### PLOT PSME VELOCITY CURVES ####
#---------------------------------#

psme_htV_M <- psme_pred_M %>% select(id, age, cht) %>% group_by(id) %>% reframe(
  dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2))) %>% ungroup()

psme_htV_Mm <- data.frame(cht = xTraj(psme_ht_M, age.M)$pop) %>% mutate(age = age.M^2) %>% 
  reframe(dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2)))

psme_wtV_M <- psme_pred_M %>% select(id, age, cwt) %>% group_by(id) %>% reframe(
  dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2))) %>% ungroup()

psme_wtV_Mm <- data.frame(cwt = xTraj(psme_wt_M, age.M)$pop) %>% mutate(age = age.M^2) %>% 
  reframe(dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2)))

psme_htV_F <- psme_pred_F %>% select(id, age, cht) %>% group_by(id) %>% reframe(
  dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2))) %>% ungroup()

psme_htV_Fm <- data.frame(cht = xTraj(psme_ht_F, age.M)$pop) %>% mutate(age = age.M^2) %>% 
  reframe(dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2)))

psme_wtV_F <- psme_pred_F %>% select(id, age, cwt) %>% group_by(id) %>% reframe(
  dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2))) %>% ungroup()

psme_wtV_Fm <- data.frame(cwt = xTraj(psme_wt_F, age.M)$pop) %>% mutate(age = age.M^2) %>% 
  reframe(dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2)))

graphics.off()
grDevices::cairo_pdf("res/fig_psme_vel_curves.pdf", height = 6, width = 6.5)
par(mfrow = c(2, 2), mar = c(4, 4.6, 1.6, 1), oma = c(0, 0, 0, 0))

mplot(x = dX, y = dY, id = id, data = psme_htV_M, col = id, las = 1, ylim=c(0, 100), xlab = 'age', ylab = 'cm/year')
title('a. Height velocity - boys', adj = 0, line = 0.6)

mplot(x = dX, y = dY, id = id, data = psme_htV_F, col = id, las = 1, ylim=c(0, 100), xlab = 'age', ylab = 'cm/year')
title('b. Height velocity - girls', adj = 0, line = 0.6)

mplot(x = dX, y = dY, id = id, data = psme_wtV_M, col = id, las = 1, ylim=c(-5, 38), xlab = 'age', ylab = 'kg/year')
title('c. Weight velocity - boys', adj = 0, line = 0.6)

mplot(x = dX, y = dY, id = id, data = psme_wtV_F, col = id, las = 1, ylim=c(-5, 38), xlab = 'age', ylab = 'kg/year')
title('d. Weight velocity - girls', adj = 0, line = 0.6)

dev.off()

rm(list = ls())

#-----------------------------------------#
#### PLOT PREDICTED VS OBSERVED VALUES ####
#-----------------------------------------#

load("gusto_ps_dat.RData")
load("gusto_psme_preds.RData")

# MALES

#dat_M %>% group_by(id) %>% filter(n() == 5) %>% ungroup %>% sample_n(1)
# N=5: 285
# N=10: 401
# N=15: 396

psme_pred_M <- psme_pred_M %>% filter(id == 285 | id == 401 | id == 396) %>% 
  rename(pred_ht = cht, pred_wt = cwt, pred_bmi = cbmi)

psme_pred_M <- dat_M %>% filter(
  id == "285" | id == "401" | id == "396") %>% 
  select(id, age, cht, cwt, cbmi) %>% mutate(id = as.integer(id)) %>% 
  full_join(psme_pred_M) %>% mutate(Nobs = ifelse(
    id == 285, "Nj: 5", ifelse(
      id == 401, "Nj: 10", "Nj: 15")
    ), sex = "Boys")

# FEMALES

#dat_F %>% group_by(id) %>% filter(n() == 10) %>% ungroup %>% sample_n(1)
# N=5: 44
# N=10: 477
# N=15: 421

psme_pred_F <- psme_pred_F %>% filter(id == 44 | id == 477 | id == 421) %>% 
  rename(pred_ht = cht, pred_wt = cwt, pred_bmi = cbmi)

psme_pred_F <- dat_F %>% filter(
  id == "44" | id == "477" | id == "421") %>% 
  select(id, age, cht, cwt, cbmi) %>% mutate(id = as.integer(id)) %>% 
  full_join(psme_pred_F) %>% mutate(Nobs = ifelse(
    id == 44, "Nj: 5", ifelse(
      id == 477, "Nj: 10", "Nj: 15")
      ), sex = "Girls")

psme_pred_MF <- bind_rows(psme_pred_M, psme_pred_F)

psme_pred_MF$Nobs <- factor(
  psme_pred_MF$Nobs, levels=c(
    "Nj: 5",
    "Nj: 10",
    "Nj: 15"
  ))

(ht_plot_MF <- ggplot(
  psme_pred_MF, aes(x = age, col = sex)) + theme_classic() + 
    geom_line(aes(y = pred_ht)) + geom_point(aes(y = cht), size = 0.6) +
    facet_wrap(. ~ Nobs, strip.position = "top", ncol = 4) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10)) +
    scale_color_brewer(palette = "Set1") +
    ggtitle("a. Height") + ylab("cm") + theme(
      legend.position = "bottom", strip.background = element_blank(),
      legend.title = element_blank(), strip.text = element_text(face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=3))))

(wt_plot_MF <- ggplot(
  psme_pred_MF, aes(x = age, col = sex)) + theme_classic() + 
    geom_line(aes(y = pred_wt)) + geom_point(aes(y = cwt), size = 0.6) +
    facet_wrap(. ~ Nobs, strip.position = "top", ncol = 4) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10)) +
    scale_color_brewer(palette = "Set1") +
    ggtitle("b. Weight") + ylab("kg") + theme(
      legend.position = "bottom", strip.background = element_blank(),
      legend.title = element_blank(), strip.text = element_text(face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=3))))

(bmi_plot_MF <- ggplot(
  psme_pred_MF, aes(x = age, col = sex)) + theme_classic() + 
    geom_line(aes(y = pred_bmi)) + geom_point(aes(y = cbmi), size = 0.6) +
    facet_wrap(. ~ Nobs, strip.position = "top", ncol = 4) +
    scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10)) +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    scale_color_brewer(palette = "Set1") +
    ggtitle("c. BMI") + ylab("kg/m²") + theme(
      legend.position = "bottom", strip.background = element_blank(),
      legend.title = element_blank(), strip.text = element_text(face="bold")) +
    guides(colour = guide_legend(override.aes = list(size=3))))

graphics.off()
grDevices::cairo_pdf("res/fig_pred_obs_MF.pdf", height = 5.8, width = 4.5)
(ht_plot_MF / wt_plot_MF / bmi_plot_MF) + plot_layout(
  guides = "collect", axis_titles = "collect", axes = "collect") & theme(
    legend.position = 'bottom')
dev.off()

rm(list = ls())

#---------------------------------#
#### BOXPLOTS GROWTH FEATURES  ####
#--------====---------------------#

load("psme_features.RData")

psme_est_MF <- psme_est_MF %>% select(
  sex, PWV, PHV, PeakBMI, PeakBMIAge, ReboundBMI, ReboundBMIAge)

psme_est_MF <- psme_est_MF %>% 
  mutate(PeakBMIAge = PeakBMIAge * 12) %>% pivot_longer(cols = c(
    PWV, PHV, PeakBMI, PeakBMIAge, ReboundBMI, ReboundBMIAge), 
    names_to = c("outs"), values_to = "count")

psme_est_MF$outs <- dplyr::recode(
  psme_est_MF$outs, 
  PHV = "a. Peak height velocity (cm/yr)",
  PWV = "b. Peak weight velocity (kg/yr)",
  PeakBMI = "c. Peak BMI (kg/m²)",
  ReboundBMI = "d. Rebound BMI (kg/m²)",
  PeakBMIAge = "e. Age at peak BMI (months)",
  ReboundBMIAge = "f. Age at rebound BMI (years)"
)

graphics.off()
grDevices::cairo_pdf("res/fig_boxplot_features_raw.pdf", height = 5.5, width = 5.5)

ggplot(psme_est_MF) + theme_classic() + 
  geom_boxplot(aes(x = count, col = sex), outlier.size = 0.3) + 
  scale_color_brewer(palette = "Set1") + 
  facet_wrap(outs ~ ., scales = "free", ncol = 2) + coord_flip() + theme(
    legend.title = element_blank(), legend.position = "right", 
    axis.title = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), strip.background = element_blank(),
    strip.text = element_text(face="bold", hjust = 0)) + guides(
      colour = guide_legend(override.aes = list(
        size = 5)))

dev.off()

psme_est_MF %>% group_by(outs, sex) %>% summarise_all(list(
  median = ~round(median(.,na.rm = T), 1),
  iqr = ~round(IQR(.,na.rm = T), 1),
  mean = ~round(mean(.,na.rm = T), 1),
  sd = ~round(sd(.,na.rm = T), 1))) %>% 
  write.csv("res/features.csv")

rm(list = ls())

#----------------------------------------#
#### GROWTH FEATURES CORRELATION PLOT ####
#----------------------------------------#

load("psme_features.RData")

psme_est_MF <- psme_est_MF %>% select(
  sex, PWV, PHV, PeakBMI, PeakBMIAge, ReboundBMI, ReboundBMIAge) %>% drop_na()

psme_est_MF <- psme_est_MF %>% rename(
  "Peak weight velocity" = PWV,
  "Peak height velocity" = PHV,
  "Peak BMI" = PeakBMI,
  "Rebound BMI" = ReboundBMI,
  "Age at peak BMI" = PeakBMIAge,
  "Age at rebound BMI" = ReboundBMIAge
)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

graphics.off()
grDevices::cairo_pdf("res/fig_corplot_features.pdf", height = 7.5, width = 7.5)
par(mfrow = c(2,1), mar = c(4, 3.6, 1.6, 0.5), oma = c(0, 4, 0, 0))

corrplot(cor(psme_est_MF %>% filter(sex == "Boys") %>% select(-sex)), 
         method = "color", col=col(200),  
         type = "upper", order = "FPC",
         tl.col = "black", tl.srt = 45,
         diag=FALSE,
         addCoef.col = "black"
)

mtext("a. Boys", side = 3, line = 0.5, font = 2, adj = 0, cex = 1.5)

corrplot(cor(psme_est_MF %>% filter(sex == "Girls") %>% select(-sex)), 
         method = "color", col=col(200),  
         type = "upper", order = "FPC",
         tl.col = "black", tl.srt = 45,
         diag=FALSE,
         addCoef.col = "black"
)

mtext("b. Girls", side = 3, line = 0.5, font = 2, adj = 0, cex = 1.5)

dev.off()

rm(list = ls())