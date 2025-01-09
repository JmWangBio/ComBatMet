
## Author: Junmin Wang
## Date: January 8th, 2025
## This script takes simulated p-values as input and calculates the TPR and FPR of the lrt applied to the simulated data of different sample sizes.
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr and ggplot2 R packages installed.
## Make sure to change the path of input to where you have saved it (line 13).

## load libraries
library(dplyr)
library(ggplot2)

## read data
sim_pval_df_all_total <- readRDS("path/to/sim_biseq_lrt_pval_data.rds")

## calculate tpr and fpr
rate.total <- sim_pval_df_all_total %>% 
  group_by(sim, batch_effect, disp_batch_effect, num.rep) %>%
  summarize(tpr = sum(pval < 0.05 & TP, na.rm = TRUE) / 
              sum(TP, na.rm = TRUE),
            fpr = sum(pval < 0.05 & (!TP), na.rm = TRUE) / 
              sum(!TP, na.rm = TRUE))

## calculate median tpr and fpr
rate.averaged <- rate.total %>% 
  group_by(batch_effect, disp_batch_effect, num.rep) %>%
  summarize(tpr_median = median(tpr, na.rm = TRUE),
            fpr_median = median(fpr, na.rm = TRUE)) %>%
  ungroup()

## plot median tpr and fpr for each combination of mean and dispersion batch effects
rate.averaged.reformatted <- rate.averaged %>%
  mutate(
    batch_effect = case_when(batch_effect == '0' ~ 'No mean\nbatch effect',
                             batch_effect == '2' ~ 'Mean batch\npercent change 2%',
                             batch_effect == '5' ~ 'Mean batch\npercent change 5%',
                             batch_effect == '10' ~ 'Mean batch\npercent change 10%'),
    disp_batch_effect = case_when(disp_batch_effect == '1' ~ 'No precision\nbatch effect',
                                  disp_batch_effect == '2' ~ 'Precision batch\nfold change 2',
                                  disp_batch_effect == '5' ~ 'Precision batch\nfold change 5',
                                  disp_batch_effect == '10' ~ 'Precision batch\nfold change 10'),
    num.rep = paste(num.rep, "samples"),
    batch_effect = factor(batch_effect, 
                          levels = c("No mean\nbatch effect",
                                     "Mean batch\npercent change 2%",
                                     "Mean batch\npercent change 5%",
                                     "Mean batch\npercent change 10%")),
    disp_batch_effect = factor(disp_batch_effect,
                        levels = c("No precision\nbatch effect",
                                   "Precision batch\nfold change 2",
                                   "Precision batch\nfold change 5",
                                   "Precision batch\nfold change 10")),
    num.rep = factor(num.rep, levels = c("20 samples", 
                                         "100 samples"))
  )

## make plot
p1 <- ggplot(data = rate.averaged.reformatted,
             aes(x = fpr_median, y = tpr_median, color = num.rep,
                 shape = num.rep)) +
  geom_point(size = 4) +
  facet_grid(disp_batch_effect ~ batch_effect) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "lightgray") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.spacing = unit(1, "lines"),
        legend.position = "bottom") +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
            colour = "black", fill = NA, inherit.aes = FALSE) +
  labs(x = "False Positive Rate (FPR)",
       y = "True Positive Rate (TPR)") +
  scale_x_continuous(limits = c(0, 0.17)) +
  scale_y_continuous(limits = c(0.4, 1.05))
