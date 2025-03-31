
## Author: Junmin Wang
## Date: March 30th, 2025

## Load libraries
library(dplyr)
library(ggplot2)

## Read data
accuracies_comparison <- read.delim("/path/to/accuracies_comparison.csv", sep = ",")

## Calculate difference pre vs post
accuracies_comparison <- accuracies_comparison %>%
  mutate(diff = Post.Adjustment.Accuracy - Pre.Adjustment.Accuracy)

## Paired t-test
t.test(accuracies_comparison$Pre.Adjustment.Accuracy, 
       accuracies_comparison$Post.Adjustment.Accuracy,
       paired = TRUE, alternative = "two.sided")

## Box plot plus dots
ggplot(accuracies_comparison,
       aes(x = "", y = diff)) +
  geom_boxplot(color = "black") +
  geom_jitter(width = 0.1, color = "lightblue") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Change in accuracy") +
  scale_y_continuous(limits = c(-0.28, 0.62)) +
  annotate(geom = "text", x = 1, y = 0.6,
           label = "Mean accuracy without adjustment: 0.70\nMean accuracy post ComBat-met: 0.81", size = 3) +
  annotate(geom = "text", x = 1, y = 0.51,
           label = "p = 3e-05", size = 3)
