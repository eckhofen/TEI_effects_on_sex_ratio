# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: October 2025
# Description: exploratory analysis of TEI chemical composition data 

# Settings ----------------------------------------------------------------

library(tidyverse)
library(gt)
library(outliers)

# Loading data ------------------------------------------------------------

# Loading main data
load("00_data/data_cc.RData")

# Loading color palettes
load("00_data/color_palettes.RData")

# Data preparation -----------------------------------------------------------------

# Creating function to remove outliers based on Dixon's Q test
remove_outlier <- function(x, p_value_threshold = 0.05) {
  # Checking for same value
  if (sum(x[1] == x) == length(x)) {
    return(x)
  }
  else {
    dixon_results <- dixon.test(x)
    if (dixon_results$p.value < p_value_threshold) {
      mean_val <- mean(x, na.rm = TRUE)
      outlier_index <- which.max(abs(x - mean_val))
      x[outlier_index] <- NA
    }
    return(x)
  }
}

# Visually checking for outliers in technical replicates
data_cc %>% group_by(fish) %>% 
  ggplot(aes(x = fish, y = ca_mass_per, label = replicate)) +
  geom_text()

# Removing replicates
data_cc_filtered <- 
  data_cc %>%
  group_by(fish, sex_binomial) %>%
  mutate(
    across(where(is.numeric), ~ remove_outlier(.))
    ) %>%
  ungroup()

table(is.na(data_cc_filtered))

# Taking mean from filtered replicates
data_cc_filtered <- data_cc %>% 
  group_by(fish, sex, treatment) %>% 
  summarise(across(where(is.numeric), ~ mean(x = ., na.rm =TRUE))) %>% 
  ungroup()

# Checking for outliers in technical replicates
data_cc_filtered %>% 
  ggplot(aes(x = fish, y = na_mass_per, label = fish)) +
  geom_text()

# Separating into mass and atomic (without the sum column)
data_cc_mass <- data_cc_filtered %>% select(fish, sex, treatment, contains("mass_"))
meta_data_cc <- data_cc[1:5] %>% group_by(fish) %>% select(-replicate) %>% unique()

# Summary -----------------------------------------------------------------

# Mass normalized
cc_summary_mass <- data_cc_mass %>%
  group_by(sex, treatment) %>%
  summarize(n = n(), across(ends_with("_mass_per"), list(mean = mean, sd = sd), .names = "{.col}_{.fn}")) %>%
  ungroup()

cc_summary_mass

# Saving as csv
write.csv(cc_summary_mass, "00_data/01_results/03_summary_mass.csv", row.names = FALSE)

# Tidy table
cc_summary_mass_tab <- 
  data_cc_filtered %>% 
  # select(sex, treatment) %>% 
  # left_join(data_cc_mass) %>% 
  group_by(sex, treatment) %>% 
  summarise(
    across(
      where(is.numeric) & !matches("^fish$"),
      list(comb = ~ paste0(
        round(mean(.x, na.rm = TRUE), 3),
        " ± ",
        round(sd(.x, na.rm = TRUE), 3), "%")), 
      .names = "{.col}")) %>% 
  select(-sex_binomial)

tab_cc_summary_mass <- 
  cc_summary_mass_tab %>%
  select(where(is.character)) %>% 
  mutate(N = 5, 
         Replicates = 3, 
         .after = treatment, 
         sex = factor(sex, levels = c("M", "F"))) %>% 
  ungroup() %>% 
  arrange(sex) %>% 
  gt() %>% 
  cols_label(sex = "Sex", 
             treatment = "Temp F2", 
             c_mass_per = "C", 
             n_mass_per = "N", 
             o_mass_per = "O", 
             na_mass_per = "Na", 
             p_mass_per = "P", 
             s_mass_per = "S", 
             ca_mass_per = "Ca") %>% 
  opt_table_font(font = "Times New Roman", size = 10) %>% 
  tab_style(style = cell_text(weight = "bold"), 
            locations = cells_column_labels()) %>% 
  tab_options() %>% 
  tab_style(style = cell_borders(sides = "bottom", 
                                 color = "black", 
                                 weight = px(1)), 
            locations = cells_column_labels())

gtsave(tab_cc_summary_mass, "00_data/01_results/03_cc_mass.html")

# Saving files for plots
save(data_cc_mass, meta_data_cc, file = "00_data/00_processed/data_chemistry.RData")
