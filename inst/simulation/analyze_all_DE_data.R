
## Author: Junmin Wang
## Date: July 18th, 2024
## This script takes simulated p-values as input and calculated the TPR and FPR of the statistical methods applied to the simulation data.
## IMPORTANT NOTE: To run this script successfully, you need to have the dplyr and ggplot2 R packages installed.
## Make sure to change the path of input to where you have saved it (line 13).

## load libraries
library(dplyr)
library(ggplot2)

## read data
sim_pval_df_all_total <- readRDS("path/to/sim_all_DE_pval_data.rds")

## calculate tpr and fpr
rate.total <- sim_pval_df_all_total %>% 
  group_by(sim, batch_effect, disp_batch_effect, method) %>%
  summarize(tpr = sum(pval < 0.05 & TP, na.rm = TRUE) / 
              sum(TP, na.rm = TRUE),
            fpr = sum(pval < 0.05 & (!TP), na.rm = TRUE) / 
              sum(!TP, na.rm = TRUE))

## calculate median tpr and fpr
rate.averaged <- rate.total %>% 
  group_by(batch_effect, disp_batch_effect, method) %>%
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
    disp_batch_effect = case_when(disp_batch_effect == '1' ~ 'No dispersion\nbatch effect',
                                  disp_batch_effect == '2' ~ 'Dispersion batch\nfold change 2',
                                  disp_batch_effect == '5' ~ 'Dispersion batch\nfold change 5',
                                  disp_batch_effect == '10' ~ 'Dispersion batch\nfold change 10'),
    method = case_when(method == 'SVA + t-test' ~ 'SVA',
                       method == "covariate in t-test" ~ 'Including batch as a covariate in DE',
                       method == 'ComBat-Met + t-test' ~ 'ComBat-met',
                       method == 'M-value ComBat + t-test' ~ 'M-value ComBat',
                       method == 'raw + t-test' ~ 'No adjustment',
                       method == 'ComBat-biseq + lrt' ~ 'ComBat-biseq'),
    batch_effect = factor(batch_effect, 
                          levels = c("No mean\nbatch effect",
                                     "Mean batch\npercent change 2%",
                                     "Mean batch\npercent change 5%",
                                     "Mean batch\npercent change 10%")),
    disp_batch_effect = factor(disp_batch_effect,
                        levels = c("No dispersion\nbatch effect",
                                   "Dispersion batch\nfold change 2",
                                   "Dispersion batch\nfold change 5",
                                   "Dispersion batch\nfold change 10")),
    method = factor(method,
                    levels = c("No adjustment",
                               "Including batch as a covariate in DE",
                               "SVA",
                               "ComBat-met",
                               "M-value ComBat",
                               "ComBat-biseq"))
  )

## make plot
p1 <- ggplot(data = rate.averaged.reformatted,
             aes(x = fpr_median, y = tpr_median, color = method,
                 shape = method)) +
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
  scale_x_continuous(limits = c(0, 0.17))
