
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
library(ggpubr)
library(scales)
library(broom)
library(psme)

#---------------------#
#### F1: SCHEMATIC ####
#---------------------#

load("gusto_psme_mods.RData")

age.M <- with(dat_M, seq.int(min(sqrt_age), sqrt(1), length.out = 150))

# HEIGHT

ht_pred_M <- data.frame(cht = xTraj(psme_ht_M, age.M)$pop) %>% mutate(age = age.M^2)
htv_pred_M <- ht_pred_M %>% reframe(dY = diff(cht)/diff(age), dX = rowMeans(embed(age, 2)))

summary(ht_pred_M$cht) # 51.50   58.09   64.35   64.15   70.29   76.05
summary(htv_pred_M$dY/12) # 1.188   1.512   2.112   2.642   3.343   7.035

graphics.off()
grDevices::cairo_pdf("res/F1.schematic_ht.pdf", height = 3, width = 4)
par(mar = c(3, 3, 2, 3) + 0.2)

ticks.1 <- seq(35,75,10)
ticks.2 <- seq(1,8,2)

plot(ht_pred_M$age*12, ht_pred_M$cht, lwd = 1.5, type = "l", yaxt="n", xlab="", ylab="", col="blue")
axis(2, at=ticks.1, col.ticks="blue", col.axis="blue")
par(new=T)
plot(htv_pred_M$dX*12, htv_pred_M$dY/12, lwd = 1.5, type = "l", xaxt="n", yaxt="n", xlab="", ylab="", col="red", lty=2)
axis(4, at=ticks.2, col.ticks="red", col.axis="red")
mtext("age - months", side=1, line=2.1, col = "black")
mtext("length - cm", side=2, line=2.1, col = "blue")
mtext("velocity - cm per month", side=4, line=2.1, col = "red")
title('a. Infant height (length)', adj = 0, line = 0.6)

dev.off()

# WEIGHT

wt_pred_M <- data.frame(cwt = xTraj(psme_wt_M, age.M)$pop) %>% mutate(age = age.M^2)
wtv_pred_M <- wt_pred_M %>% reframe(dY = diff(cwt)/diff(age), dX = rowMeans(embed(age, 2)))

summary(wt_pred_M$cwt) # 3.468   5.428   7.156   6.929   8.509   9.652
summary(wtv_pred_M$dY/12) # 0.2343  0.3131  0.5356  0.7097  0.9761  2.0961

graphics.off()
grDevices::cairo_pdf("res/F1.schematic_wt.pdf", height = 3, width = 4)
par(mar = c(3, 3, 2, 3) + 0.2)

ticks.1 <- seq(3, 10, 1)
ticks.2 <- seq(0.2, 2.1, 0.3)

plot(wt_pred_M$age*12, wt_pred_M$cwt, lwd = 1.5, type = "l", yaxt="n", xlab="", ylab="", col="blue")
axis(2, at=ticks.1, col.ticks="blue", col.axis="blue")
par(new=T)
plot(wtv_pred_M$dX*12, wtv_pred_M$dY/12, lwd = 1.5, type = "l", xaxt="n", yaxt="n", xlab="", ylab="", col="red", lty=2)
axis(4, at=ticks.2, col.ticks="red", col.axis="red")
mtext("age - months", side=1, line=2.1, col = "black")
mtext("weight - kg", side=2, line=2.1, col = "blue")
mtext("velocity - kg per month", side=4, line=2.1, col = "red")
title('b. Infant weight', adj = 0, line = 0.6)

dev.off()

# BMI

age.M <- with(dat_M, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))
bmi_pred_M <- data.frame(cbmi = xTraj(psme_bmi_M, age.M)$pop) %>% mutate(age = age.M^2)

graphics.off()
grDevices::cairo_pdf("res/F1.schematic_bmi.pdf", height = 3, width = 4)
par(mfrow = c(1, 1), mar = c(4, 4.6, 1.6, 1), oma = c(0, 0, 0, 0))

plot(
  bmi_pred_M$age, bmi_pred_M$cbmi, lwd = 1.5, type = "l", 
  xlab = '', ylab = '', xaxt="n", yaxt="n", ylim = c(13, 18.9))
axis(1, at=seq(0, 10, 1))
axis(2, at=seq(13, 19, 2))
mtext("age - years", side=1, line=2.1)
mtext("BMI - kg/m²", side=2, line=2.1)
title('c. Infant to childhood BMI', adj = 0, line = 0.6)

dev.off()

#---------------------------------------#
#### F2: PLOT OBSERVED GROWTH VALUES ####
#---------------------------------------#

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
grDevices::cairo_pdf("res/F2.observed_values.pdf", height = 6, width = 5)

ggplot(data = dat_MF, aes(x = age, y = count)) + theme_classic() + 
  geom_point(size = 0.01)  + facet_wrap(
  outsex ~ ., strip.position="top", nrow = 3, scales = "free") + theme(
    legend.title = element_blank(), panel.grid.minor.x = element_blank(), 
    legend.position = "none", strip.text = element_text(face="bold"),
        strip.background = element_blank()) + ylab('observed values') + xlab("age - years") +
  scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10))

dev.off()

rm(list = ls())

#-------------------------------------------------#
#### F3/F4: PLOT PSME DISTANCE/VELOCITY CURVES ####
#-------------------------------------------------#

load("gusto_psme_preds.RData")
load("gusto_psme_mods.RData")

age.M <- with(dat_M, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))
age.F <- with(dat_F, seq.int(min(sqrt_age), max(sqrt_age), length.out = 550))

# F3: DISTANCE

graphics.off()
grDevices::cairo_pdf("res/F3. psme_dist_curves.pdf", height = 6.5, width = 5.5)
par(mfrow = c(3, 2), mar = c(4, 4.6, 1.6, 1), oma = c(0, 0, 0, 0))

sitar::mplot(x = age, y = cht, id = id, data = psme_pred_M, col = id, las = 1, ylim = c(35, 165), xlab = 'age', ylab = 'cm')
lines(age.M^2, xTraj(psme_ht_M, age.M)$pop, type = "l", col = 'black', lwd = 3.5)
title('a. Height - boys', adj = 0, line = 0.6)

sitar::mplot(x = age, y = cht, id = id, data = psme_pred_F, col = id, las = 1, ylim = c(35, 165), xlab = 'age', ylab = 'cm')
lines(age.F^2, xTraj(psme_ht_F, age.F)$pop, type = "l", col = 'black', lwd = 3.5)
title('b. Height - girls', adj = 0, line = 0.6)

sitar::mplot(x = age, y = cwt, id = id, data = psme_pred_M, col = id, las = 1, ylim = c(0, 80), xlab = 'age', ylab = 'kg')
lines(age.M^2, xTraj(psme_wt_M, age.M)$pop, type = "l", col = 'black', lwd = 3.5)
title('c. Weight - boys', adj = 0, line = 0.6)

sitar::mplot(x = age, y = cwt, id = id, data = psme_pred_F, col = id, las = 1, ylim = c(0, 80), xlab = 'age', ylab = 'kg')
lines(age.F^2, xTraj(psme_wt_F, age.F)$pop, type = "l", col = 'black', lwd = 3.5)
title('d. Weight - girls', adj = 0, line = 0.6)

sitar::mplot(x = age, y = cbmi, id = id, data = psme_pred_M, col = id, las = 1, ylim = c(8, 35), xlab = 'age', ylab = 'kg/m²')
lines(age.M^2, xTraj(psme_bmi_M, age.M)$pop, type = "l", col = 'black', lwd = 3.5)
title('e. BMI - boys', adj = 0, line = 0.6)

sitar::mplot(x = age, y = cbmi, id = id, data = psme_pred_F, col = id, las = 1, ylim = c(8, 35), xlab = 'age', ylab = 'kg/m²')
lines(age.F^2, xTraj(psme_bmi_F, age.F)$pop, type = "l", col = 'black', lwd = 3.5)
title('f. BMI - girls', adj = 0, line = 0.6)

dev.off()

# F4: VELOCITY

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
grDevices::cairo_pdf("res/F5. psme_vel_curves.pdf", height = 6, width = 6.5)
par(mfrow = c(2, 2), mar = c(4, 4.6, 1.6, 1), oma = c(0, 0, 0, 0))

sitar::mplot(x = dX, y = dY, id = id, data = psme_htV_M, col = id, las = 1, ylim=c(0, 100), xlab = 'age', ylab = 'cm/year')
title('a. Height velocity - boys', adj = 0, line = 0.6)

sitar::mplot(x = dX, y = dY, id = id, data = psme_htV_F, col = id, las = 1, ylim=c(0, 100), xlab = 'age', ylab = 'cm/year')
title('b. Height velocity - girls', adj = 0, line = 0.6)

sitar::mplot(x = dX, y = dY, id = id, data = psme_wtV_M, col = id, las = 1, ylim=c(-5, 38), xlab = 'age', ylab = 'kg/year')
title('c. Weight velocity - boys', adj = 0, line = 0.6)

sitar::mplot(x = dX, y = dY, id = id, data = psme_wtV_F, col = id, las = 1, ylim=c(-5, 38), xlab = 'age', ylab = 'kg/year')
title('d. Weight velocity - girls', adj = 0, line = 0.6)

dev.off()

rm(list = ls())

#---------------------------------------------#
#### F5: PLOT PREDICTED VS OBSERVED VALUES ####
#---------------------------------------------#

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
grDevices::cairo_pdf("res/F4.pred_obs_MF.pdf", height = 5.8, width = 4.5)
(ht_plot_MF / wt_plot_MF / bmi_plot_MF) + plot_layout(
  guides = "collect", axis_titles = "collect", axes = "collect") & theme(legend.position = 'bottom')
dev.off()

rm(list = ls())

#-------------------------------------#
#### F6: BOXPLOTS GROWTH FEATURES  ####
#--------====-------------------------#

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

#

graphics.off()
grDevices::cairo_pdf("res/F6.boxplot_features_raw.pdf", height = 5.5, width = 5.5)

ggplot(psme_est_MF) + theme_classic() + 
  geom_boxplot(aes(x = count, col = sex), outlier.size = 0.3) + 
  scale_color_brewer(palette = "Set1") + 
  facet_wrap(outs ~ ., scales = "free", ncol = 2) + coord_flip() + theme(
    legend.title = element_blank(), legend.position = "right", 
    axis.title = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), strip.background = element_blank(),
    strip.text = element_text(face="bold")) + guides(
      colour = guide_legend(override.aes = list(
        size = 5)))

dev.off()

rm(list = ls())

#----------------------------------------------#
#### F7: GROWTH FEATURES CORRELATION MATRIX ####
#----------------------------------------------#

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
grDevices::cairo_pdf("res/corplot_features.pdf", height = 7.5, width = 7.5)
par(mfrow = c(2,1), mar = c(4, 3.6, 1.6, 0.5), oma = c(0, 4, 0, 0))

#

corrplot(cor(psme_est_MF %>% filter(sex == "Boys") %>% select(-sex)), 
         method = "color", col=col(200),  
         type = "upper", order = "FPC",
         tl.col = "black", tl.srt = 45,
         diag=FALSE,
         addCoef.col = "black"
)

mtext("a. Boys", side = 3, line = 0.5, font = 2, adj = 0, cex = 1.5)

#

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

#-----------------------------------------#
#### F8: PLOT PSME ASSOCIATION RESULTS ####
#-----------------------------------------#

load("psme_features.RData")

psme_est_MF <- psme_est_MF %>% mutate_at(
  c("m_ht", "m_wt", "m_age", "GA", "BLT", "BWT"), .funs = list(z = ~ scale(.))) %>% 
  as.data.frame()

exp <- c("m_ht_z", "m_wt_z", "m_age_z", "GA_z", "BLT_z", "BWT_z")
out <- c("PHV", "PWV", "PeakBMI", "ReboundBMI", "PeakBMIAge * 12", "ReboundBMIAge")

(psme_assoc_res <- expand.grid(out, exp) %>% group_by(Var1) %>% 
    rowwise() %>% summarise(frm = paste0(Var1, " ~ sex + ", Var2)) %>% 
    group_by(model_id = row_number(), frm) %>% do(cbind(tidy(
      lm(.$frm, data = psme_est_MF)))) %>% mutate(
        lci = estimate - (1.96 * std.error), 
        uci = estimate + (1.96 * std.error)) %>% 
    filter(term != "sexGirls" & term != "(Intercept)") %>% 
    select(-std.error, -std.error, -statistic) %>% 
    as.data.frame())

psme_assoc_res$term <- dplyr::recode(
  psme_assoc_res$term, 
  "m_age_z" = "Maternal age", 
  "m_wt_z" = "Maternal weight", 
  "m_ht_z" = "Maternal height", 
  "GA_z" = "Gestational age", 
  "BWT_z" = "Birth weight",
  "BLT_z" = "Birth length"
)

psme_assoc_res$term <- factor(
  psme_assoc_res$term, levels=c(
    "Maternal age",
    "Maternal weight",
    "Maternal height",
    "Gestational age",
    "Birth weight",
    "Birth length"    
  ))

psme_assoc_res$outs <- rep(c(
  "a. Peak height velocity (cm/yr)",
  "b. Peak weight velocity (kg/yr)",
  "c. Peak BMI (kg/m²)",
  "d. Rebound BMI (kg/m²)",
  "h. Age at peak BMI (months)",
  "i. Age at rebound BMI (years)"), 6)

rm(psme_est_MF, exp, out)

graphics.off()
grDevices::cairo_pdf("res/F8.assoc_res.pdf", height = 6, width = 4)

ggplot(
  data = psme_assoc_res, aes(x = term, y = estimate, ymin = lci, ymax = uci)) + 
  theme_classic() + geom_pointrange(aes(col = term)) + coord_flip() +
  scale_color_brewer(palette = "Dark2") + 
  facet_wrap(outs ~ ., scales = "free", ncol = 2) +
  geom_hline(yintercept = 0, linewidth = 0.2, col = "red")  + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), legend.title = element_blank(), 
        legend.position  = "bottom",
        strip.background = element_blank(),, strip.text = element_text(face="bold"),
        strip.text.y.left = element_text(angle = 0), plot.title = element_text(
          face = 'bold', hjust = 0),axis.line.y = element_blank(), 
        axis.title.x = element_text(size = 9)) + guides(
          colour = guide_legend(override.aes = list(size = 0.5))) + 
  scale_x_discrete(limits = rev) + ylab("Mean difference") +
  scale_y_continuous(labels = label_number(accuracy = 0.1))

dev.off()

rm(list = ls())

#-----------------------------------#
#### F9: PLOT PSAL CORRELATIONS ####
#-----------------------------------#

load("psal_features.RData")

summary(est_MF[est_MF$sex == "Boys",]$PHV_psme) # 67 95 
summary(est_MF[est_MF$sex == "Boys",]$PHV_dp) # 61 96 

summary(est_MF[est_MF$sex == "Boys",]$PWV_psme) # 14 35
summary(est_MF[est_MF$sex == "Boys",]$PWV_dp) # 10 34 

summary(est_MF[est_MF$sex == "Boys",]$PeakBMI_psme) # 13 23
summary(est_MF[est_MF$sex == "Boys",]$PeakBMI_dp) # 13 23 

summary(est_MF[est_MF$sex == "Boys",]$ReboundBMI_psme) # 12 25
summary(est_MF[est_MF$sex == "Boys",]$ReboundBMI_dp) # 12 25

summary(est_MF[est_MF$sex == "Boys",]$PeakBMIAge_psme*12) # 3 15
summary(est_MF[est_MF$sex == "Boys",]$PeakBMIAge_dp*12) # 3 17

summary(est_MF[est_MF$sex == "Boys",]$ReboundBMIAge_psme) # 2 9
summary(est_MF[est_MF$sex == "Boys",]$ReboundBMIAge_dp) # 2 9

graphics.off()
grDevices::cairo_pdf("res/F9. psal_cor.pdf", height = 7, width = 6)
par(mgp=c(2,1,0), mfrow = c(3, 2), mar = c(4, 4.6, 1.6, 1), oma = c(0, 0, 0, 0))

with(est_MF[est_MF$sex == "Boys",], plot(PHV_psme, PHV_dp, xlab = "standard method", ylab = "double penalty", xlim=c(60, 97), ylim=c(60, 97)))
title(paste("Peak height velocity (cm/yr), r =", round(with(est_MF[est_MF$sex == "Boys",], cor(PHV_psme, PHV_dp)), 2)), adj = 0, line = 0.6)

with(est_MF[est_MF$sex == "Boys",], plot(PWV_psme, PWV_dp, xlab = "standard method", ylab = "double penalty", xlim=c(9, 36), ylim=c(9, 36)))
title(paste("Peak weight velocity (kg/yr), r =", round(with(est_MF[est_MF$sex == "Boys",], cor(PWV_psme, PWV_dp)), 2)), adj = 0, line = 0.6)

with(est_MF[est_MF$sex == "Boys",], plot(PeakBMI_psme, PeakBMI_dp, xlab = "standard method", ylab = "double penalty", xlim=c(12, 24), ylim=c(12, 24)))
title(paste("Peak BMI (kg/m²), r =", round(with(est_MF[est_MF$sex == "Boys",], cor(PeakBMI_psme, PeakBMI_dp, use='complete.obs')), 3)), adj = 0, line = 0.6)

with(est_MF[est_MF$sex == "Boys",], plot(ReboundBMI_psme, ReboundBMI_dp, xlab = "standard method", ylab = "double penalty", xlim=c(10, 26), ylim=c(10, 26)))
title(paste("Rebound BMI (kg/m²), r =", round(with(est_MF[est_MF$sex == "Boys",], cor(ReboundBMI_psme, ReboundBMI_dp, use='complete.obs')), 3)), adj = 0, line = 0.6)

with(est_MF[est_MF$sex == "Boys",], plot(PeakBMIAge_psme*12, PeakBMIAge_dp*12, xlab = "standard method", ylab = "double penalty", xlim=c(2, 18), ylim=c(2, 18)))
title(paste("Age at peak BMI (mo), r =", round(with(est_MF[est_MF$sex == "Boys",], cor(PeakBMIAge_psme*12, PeakBMIAge_dp*12, use='complete.obs')), 2)), adj = 0, line = 0.6)

with(est_MF[est_MF$sex == "Boys",], plot(ReboundBMIAge_psme, ReboundBMIAge_dp, xlab = "standard method", ylab = "double penalty", xlim=c(2, 10), ylim=c(2, 10)))
title(paste("Age at rebound BMI (yr), r =", round(with(est_MF[est_MF$sex == "Boys",], cor(ReboundBMIAge_psme, ReboundBMIAge_dp, use='complete.obs')), 2)), adj = 0, line = 0.6)

dev.off()
