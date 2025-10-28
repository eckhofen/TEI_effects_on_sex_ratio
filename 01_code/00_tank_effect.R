# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: February 2025
# Description: Exploring tank effect on sex ratios and growth variables

# Settings ----------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(rstatix)

# Loading data ------------------------------------------------------------

# Main data
load("00_data/data_main.RData")

# Loading color palettes
load("00_data/color_palettes.RData")


# Tank effect -------------------------------------------------------------

# Preparing data for tank comparison for each group
data_tank <-
  data_main %>%
  group_by(treatment, tank_pregrow) %>%
  mutate(n = n(),
         n_male = sum(sex == "M"),
         mr = mean(n_male / n)) %>% ungroup() %>%
  group_by(treatment) %>%
  mutate(sd = sd(mr),
         se = sd / sqrt(2),
         mr_mean = mean(mr)) %>%
  select(treatment, tank_pregrow, n, n_male, mr, mr_mean, sd, se) %>%
  distinct() %>%
  arrange(treatment) %>%
  ungroup()

# Testing significance with Fisher's exact
sex_counts <-
  data_main %>%
  count(treatment, tank_pregrow, sex) %>%
  pivot_wider(names_from = sex, values_from = n, values_fill = 0)

# Looping through treatment to run Fisher's exact test for each group
fisher_results <- 
  sex_counts %>%
  group_split(treatment) %>%
  map_df(~{
    if (nrow(.x) == 2) {
      tbl <- matrix(c(.x$M, .x$F), nrow = 2, byrow = TRUE)
      ft <- fisher.test(tbl)
      tibble(group = unique(.x$treatment),
             tank_1 = .x$tank_pregrow[1],
             tank_2 = .x$tank_pregrow[2],
             n_1 = sum(.x$F[1], .x$M[1]), 
             n_2 = sum(.x$F[2], .x$M[2]), 
             odds_ratio = ft$estimate,
             p = ft$p.value)
    } else {
      NULL  
    }
  })
fisher_results

# All comparisons not significant

# Fig. S1 --------------------------------------------------------------------

p_tank_sex_by_group <- 
  ggplot(data_tank, aes(x = treatment, y = mr, fill = tank_pregrow)) +
  geom_col(position = "dodge") +
  labs(y = "F2 male proportion", 
       fill = "Tank", 
       x = "Temperature history (F0, F1, F2)") +
  scale_fill_manual(values = color_tank) +
  theme_pubr() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggsave("02_plots/00_tank_sex_overview_groups.pdf", p_tank_sex_by_group, width = 6, height = 5)
ggsave("02_plots/00_tank_sex_overview_groups.png", p_tank_sex_by_group, width = 6, height = 5)
