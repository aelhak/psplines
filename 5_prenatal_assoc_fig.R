#---------------------# 
#### LOAD PACKAGES #### 
#---------------------# 

library(RColorBrewer) 
library(tidyverse) 
library(patchwork) 
library(scales) 
library(broom) 

#-----------------------------------#
#### ASSOCIATION ANALYSIS & PLOT #### 
#-----------------------------------# 

load("psme_features.RData") 

psme_est_MF <- psme_est_MF %>% mutate_at( c("m_ht", "m_wt", "m_age", "GA", "BLT", "BWT"), .funs = list(z = ~ scale(.))) %>% as.data.frame() 

exp <- c("m_ht_z", "m_wt_z", "m_age_z", "GA_z", "BLT_z", "BWT_z") 
out <- c("PHV", "PWV", "PeakBMI", "ReboundBMI", "PeakBMIAge * 12", "ReboundBMIAge * 12") 

(psme_assoc_res <- expand.grid(out, exp) %>% group_by(Var1) %>% rowwise() %>% 
    summarise(frm = paste0(Var1, " ~ sex + ", Var2)) %>% group_by(model_id = row_number(), frm) %>% 
    do(cbind(tidy( lm(.$frm, data = psme_est_MF)))) %>% mutate( 
      lci = estimate - (1.96 * std.error), uci = estimate + (1.96 * std.error)) %>% 
    filter(term != "sexGirls" & term != "(Intercept)") %>% select(-std.error, -std.error, -statistic) %>% 
    as.data.frame()) 

psme_assoc_res$term <- dplyr::recode(
  psme_assoc_res$term, "m_age_z" = "Maternal age", "m_wt_z" = "Maternal weight", 
  "m_ht_z" = "Maternal height", "GA_z" = "Gestational age", "BWT_z" = "Birth weight", 
  "BLT_z" = "Birth length" ) 

psme_assoc_res$term <- factor(
  psme_assoc_res$term, levels=c(
    "Maternal age", 
    "Maternal weight", 
    "Maternal height", 
    "Gestational age", 
    "Birth weight", 
    "Birth length" )) 

psme_assoc_res$outs <- rep(c(
  "a. Peak height velocity", 
  "b. Peak weight velocity", 
  "c. Peak BMI", 
  "d. Rebound BMI", 
  "e. Age at peak BMI", 
  "f. Age at rebound BMI"), 6) 

rm(psme_est_MF, exp, out) 

# PLOT

assoc_plot <- function(data, outcome_filter, y_lab, plot_title) {
  
  data %>%
    filter(outs == outcome_filter) %>% ggplot(
      aes(x = term, y = estimate, ymin = lci, ymax = uci)) +
    geom_pointrange(aes(color = term)) + geom_hline(
      yintercept = 0, linewidth = 0.2, color = "red") +
    coord_flip() + scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    scale_x_discrete(limits = rev) + labs(
      title = plot_title, y = y_lab, x = NULL) +
    theme_classic() + theme(
      axis.ticks.y = element_blank(),
      legend.position = "none")
}

(p1 <- assoc_plot(psme_assoc_res, "a. Peak height velocity", "cm/yr", "a. Peak height velocity"))
(p2 <- assoc_plot(psme_assoc_res, "b. Peak weight velocity", "kg/yr", "b. Peak weight velocity"))
(p3 <- assoc_plot(psme_assoc_res, "c. Peak BMI", "kg/m²", "c. Peak BMI"))
(p4 <- assoc_plot(psme_assoc_res, "d. Rebound BMI", "kg/m²", "d. Rebound BMI"))
(p5 <- assoc_plot(psme_assoc_res, "e. Age at peak BMI", "months", "e. Age at peak BMI"))
(p6 <- assoc_plot(psme_assoc_res, "f. Age at rebound BMI", "months", "f. Age at rebound BMI"))

graphics.off() 
grDevices::cairo_pdf("res/fig.prenatal_assoc.pdf", height = 3, width = 13.5) 
(p1 | p2 | p3 | p4 | p5 | p6) + plot_layout( 
  guides = "collect", axis_titles = "collect", axes = "collect") & 
  theme( legend.position = 'none') 
dev.off()