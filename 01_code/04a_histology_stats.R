# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: October 2025
# Description: Histological analyses

# Settings ----------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(ggpubr)
library(rstatix)
library(tidyverse)
library(gt)

# Loading data ------------------------------------------------------------

# Loading Histology data
# Compiled RData from two .csv files with females and males
load("00_data/data_histology.RData")

# Loading main data for comparison
load("00_data/data_main.RData")

# Loading color palettes
load("00_data/color_palettes.RData")


# Summary -----------------------------------------------------------------

# Females
summary_hist_f <- data_histology_f %>% 
  group_by(sample, treatment, sample_id) %>%
  summarise(across(c(og, pg), sum), .groups = "drop") %>% 
  arrange(treatment, sample) %>% 
  rename(
    Og = og, 
    Pg = pg
    ) %>% 
  left_join(data_main %>% select(-treatment), by = "sample_id") %>% 
  ungroup()

# Defining new cell types of females
cell_females <- c("Og", "Pg")

# Long format
summary_hist_f_long <- summary_hist_f %>% 
  pivot_longer(cols = all_of(cell_females), names_to = "cell_type", values_to = "count") %>% 
  mutate(cell_type = factor(cell_type, levels = cell_females)) %>% 
  mutate(treatment = factor(treatment, levels = unique(treatment))) %>% 
  mutate(sample = factor(sample, levels = unique(sample))) %>% 
  ungroup()

# Reordering the samples accordingly
summary_hist_f_long$sample <- 
  factor(summary_hist_f_long$sample,levels = unique(summary_hist_f_long$sample))

# Comparing size of histology samples with main ones

# Histology samples
stat_hist_f <- summary_hist_f %>%
  select(temp_f2, sample_id, length_standard_mm, weight_g, Og, Pg) %>%
  group_by(temp_f2) %>% 
  summarise(
    across(where(is.numeric), 
           list(stat = ~ paste(round(mean(.x), 2), "±", round(sd(.x), 2)))), 
    n = n()
    )

# All samples
data_main %>% 
  filter(sex == "F") %>% 
  group_by(temp_f2) %>% 
  select(length_standard_mm, weight_g) %>% 
  summarise(
    n = n(), 
    across(c(length_standard_mm, weight_g), 
           list(stat = ~ paste(round(mean(.x), 2), "±", round(sd(.x), 2))))
    )

# Pretty table
tab_summary_females <- 
  stat_hist_f %>% 
  gt() %>% 
  cols_label(
    temp_f2 = "Temp F2",
    Og_stat = "Og (n)",
    Pg_stat = "Pg (n)",
    weight_g_stat = "Weight (g)",
    length_standard_mm_stat = "Length (mm)"
    ) %>% 
  opt_table_font(font = "Times New Roman", size = 10) %>% 
  tab_style(style = cell_text(weight = "bold"), 
            locations = cells_column_labels()) %>% 
  tab_options() %>% 
  tab_style(style = cell_borders(sides = "bottom", 
                                 color = "black", 
                                 weight = px(1)), 
            locations = cells_column_labels())

gtsave(tab_summary_females, "03_tables/04_summary_females.html")
gtsave(tab_summary_females, "03_tables/04_summary_females.docx")

# Males
summary_hist_m <- data_histology_m %>% 
  group_by(sample, treatment, sample_id) %>%
  summarise(across(c(aund, adiff, b, spc, spt), sum), .groups = "drop") %>% 
  arrange(treatment, sample) %>% 
  rename(
    Aund = aund,
    Adiff = adiff,
    B = b, 
    Spc = spc, 
    Spt = spt
    ) %>% 
  left_join(data_main %>% select(-treatment), by = "sample_id") %>% 
  ungroup()

# Defining cells for males
cell_males <- c("Aund", "Adiff", "B", "Spc", "Spt")

# Long format
summary_hist_m_long <- summary_hist_m %>% 
  pivot_longer(cols = all_of(cell_males), names_to = "cell_type", values_to = "count") %>% 
  mutate(cell_type = factor(cell_type, levels = cell_males)) %>% 
  mutate(treatment = factor(treatment, levels = unique(treatment))) %>% 
  mutate(sample = factor(sample, levels = unique(sample))) %>% 
  ungroup()

# Reordering the samples accordingly
summary_hist_m_long$sample <- factor(summary_hist_m_long$sample, 
                                     levels = unique(summary_hist_m_long$sample))

# Comparing size of histology samples with main ones

# Histology samples
stat_hist_m <- summary_hist_m %>%
  select(temp_f2, sample_id, length_standard_mm, weight_g, Aund, Adiff, B, Spc, Spt) %>%
  group_by(temp_f2) %>% 
  summarise(
    across(where(is.numeric), 
           list(stat = ~ paste(round(mean(.x), 2), "±", round(sd(.x), 2)))), 
    n = n())

# All samples
data_main %>% 
  filter(sex == "M") %>% 
  group_by(temp_f2) %>% 
  select(length_standard_mm, weight_g) %>% 
  summarise(
    n = n(), 
    across(c(length_standard_mm, weight_g), 
           list(stat = ~ paste(round(mean(.x), 2), "±", round(sd(.x), 2))))
    )

# Pretty table
tab_summary_males <- 
  stat_hist_m %>% 
  gt() %>% 
  cols_label(
    temp_f2 = "Temp F2",
    Aund_stat = "Aund (n)",
    Adiff_stat = "Adiff (n)",
    B_stat = "B (n)",
    Spc_stat = "Spc (n)",
    Spt_stat = "Spt (n)",
    weight_g_stat = "Weight (g)",
    length_standard_mm_stat = "Length (mm)"
    ) %>% 
  opt_table_font(font = "Times New Roman", size = 10) %>% 
  tab_style(style = cell_text(weight = "bold"), 
            locations = cells_column_labels()) %>% 
  tab_options() %>% 
  tab_style(style = cell_borders(sides = "bottom", 
                                 color = "black", 
                                 weight = px(1)), 
            locations = cells_column_labels())

gtsave(tab_summary_males, "03_tables/04_summary_males.html")
gtsave(tab_summary_males, "03_tables/04_summary_males.docx")



# Summary statistics ------------------------------------------------------

# Females
sum_stat_females <- summary_hist_f %>% 
  summarise(length_mean = mean(length_standard_mm),
            length_sd = sd(length_standard_mm),
            weight_mean = mean(weight_g),
            weight_sd = sd(weight_g), 
            n = n()) %>% 
  mutate(sample_set = "Histology Females", .before = 1)

# Males
sum_stat_males <- summary_hist_m %>% 
  summarise(length_mean = mean(length_standard_mm),
            length_sd = sd(length_standard_mm),
            weight_mean = mean(weight_g),
            weight_sd = sd(weight_g), 
            n = n()) %>% 
  mutate(sample_set = "Histology Males", .before = 1)

# Summary of all fish 
# Females 
sum_stat_females_all <- data_main %>% 
  filter(sex == "F") %>% 
  summarise(length_mean = mean(length_standard_mm),
            length_sd = sd(length_standard_mm),
            weight_mean = mean(weight_g),
            weight_sd = sd(weight_g), 
            n = n()) %>% 
  mutate(sample_set = "All females", .before = 1)

# Males 
sum_stat_males_all <- data_main %>% 
  filter(sex == "M") %>% 
  summarise(length_mean = mean(length_standard_mm),
            length_sd = sd(length_standard_mm),
            weight_mean = mean(weight_g),
            weight_sd = sd(weight_g), 
            n = n()) %>% 
  mutate(sample_set = "All Males", .before = 1)

sum_stat <- rbind(sum_stat_males, sum_stat_males_all, sum_stat_females, sum_stat_females_all)
write_csv(sum_stat, file = "00_data/01_results/04_hist_summary_stats.csv")


# Comparing histology samples with all samples ----------------------------

# whole sample body size
growth_whole_sample <- 
  data_main %>% 
  select(sex, length_standard_mm) %>% 
  mutate(group = "whole_sample")

# males
growth_selected_m <- 
  data_histology_m %>%
  select(sex, length_standard_mm, sample_id, sample) %>% 
  unique() %>% 
  mutate(group = "selected")

# males
growth_selected_f <- 
  data_histology_f %>%
  select(sex, length_standard_mm, sample_id, sample) %>% 
  unique() %>% 
  mutate(group = "selected")

# combined
growth_selected <- 
  rbind(growth_selected_f, growth_selected_m)

# Testing if histology samples length have different distributions than all samples
# Weight was not considered due to larger variability 

# Females
ks.test(growth_whole_sample %>% filter(sex == "F") %>% pull(length_standard_mm), growth_selected_f$length_standard_mm)

# Males
ks.test(growth_whole_sample %>% filter(sex == "M") %>% pull(length_standard_mm), growth_selected_m$length_standard_mm)


# Saving data for plotting ------------------------------------------------

save(summary_hist_m_long, summary_hist_f_long, file = "00_data/00_processed/data_histology.RData")
