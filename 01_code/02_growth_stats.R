# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: October 2025
# Description: Analysis of growth data 

# Settings ----------------------------------------------------------------

# Libraries
library(tidyverse)
library(rstatix)
library(gt)
library(gtsummary)
library(gtExtras)

# Loading data ------------------------------------------------------------

# All data as
load("00_data/data_main.RData")

# Loading color palettes
load("00_data/color_palettes.RData")

# General statistics ------------------------------------------------------

# Summary
growth_summary <- data_main %>%
  group_by(treatment, temp_f2, sex) %>%
  summarize(n = n(),
            BW_mean = mean(weight_g, na.rm = TRUE),
            BW_sd = sd(weight_g, na.rm = TRUE),
            SL_mean = mean(length_standard_mm, na.rm = TRUE),
            SL_sd = sd(length_standard_mm, na.rm = TRUE),
            .groups = "drop")

write.csv(growth_summary, "00_data/01_results/02_growth_summary.csv", row.names = FALSE)

# Even more detailed summary
growth_summary_det <- data_main %>%
  group_by(treatment, temp_f2, temp_f1, temp_f0, line_f0, sex) %>%
  summarize(n = n(),
            BW_mean = mean(weight_g, na.rm = TRUE),
            BW_sd = sd(weight_g, na.rm = TRUE),
            SL_mean = mean(length_standard_mm, na.rm = TRUE),
            SL_sd = sd(length_standard_mm, na.rm = TRUE),
            .groups = "drop")

write.csv(growth_summary_det, "00_data/01_results/02_growth_summary_det.csv", row.names = FALSE)

# Summary based on F0 line
growth_summary_by_f0 <- 
  data_main %>% 
  group_by(line_f0, sex) %>% 
  summarize(n = n(), 
            bw_mean = round(mean(weight_g), 1), 
            bw_sd = round(sd(weight_g), 1),
            weight_w_sd = paste(bw_mean, "±", bw_sd), 
            sl_mean = round(mean(length_standard_mm), 1),
            sl_sd = round(sd(length_standard_mm), 1), 
            length_w_sd = paste(sl_mean, "±", sl_sd), 
            .groups = "drop")

# Summary more detailed
growth_summary_by_f0_det <- 
  data_main %>% 
  group_by(line_f0, sex, temp_f2, temp_f1) %>% 
  summarize(n = n(), 
            bw_mean = round(mean(weight_g), 1), 
            bw_sd = round(sd(weight_g), 1),
            weight_w_sd = paste(bw_mean, "±", bw_sd), 
            sl_mean = round(mean(length_standard_mm), 1),
            sl_sd = round(sd(length_standard_mm), 1), 
            length_w_sd = paste(sl_mean, "±", sl_sd),
            bw_values_for_plot = list(weight_g),
            sl_values_for_plot = list(length_standard_mm),
            .groups = "drop")

# Creating readable table
# Setting up summary table function 
table_styles <- 
  function(table) {
    table %>% 
      opt_horizontal_padding(scale = 3) %>% 
      opt_vertical_padding(scale = .1) %>% 
      opt_table_font(font = "Times New Roman", 
                     size = 10, 
                     color = "black") %>% 
      tab_style(cell_text(weight = "bold", align = "center"), 
                locations = cells_row_groups()) %>%
      tab_style(cell_fill(color = "grey90"), 
                locations = cells_row_groups()) %>%
      tab_style(cell_text(weight = "bold"), 
                locations = cells_column_labels()) %>% 
      opt_table_lines(extent = "none")
    
  }

# General summary
tab_growth_summary <- 
  growth_summary_by_f0 %>% 
  select(line_f0, sex, n, weight_w_sd, length_w_sd) %>% 
  arrange(sex) %>% 
  gt(groupname_col = "line_f0") %>%
  cols_label(sex = "Sex", 
             n = "N",
             weight_w_sd = "Average weight ± SD",
             length_w_sd = "Average length ± SD") %>% 
  table_styles

# Saving table
gtsave(data = tab_growth_summary, filename = "00_data/01_results/02_growth_summary_by_f0.html")

# More detailed summary
tab_growth_summary_det <- 
  growth_summary_by_f0_det %>% 
  select(line_f0, sex, n, temp_f1, temp_f2, weight_w_sd, length_w_sd) %>% 
  arrange(sex, temp_f1) %>% 
  gt(groupname_col = "line_f0") %>%
  cols_label(sex = "Sex", 
             n = "N", 
             temp_f1 = "Temp F1", 
             temp_f2 = "Temp F2",
             weight_w_sd = "Average weight ± SD",
             length_w_sd = "Average length ± SD") %>% 
  table_styles

gtsave(data = tab_growth_summary_det, filename = "00_data/01_results/02_growth_summary_by_f0_det.html")

# Simplifying data for comparison and anova
# Creating pivot longer table 
data_growth <- data_main %>% 
  select(sex, line_f0, temp_f2, temp_f1, temp_f0, weight_g, length_standard_mm) %>% 
  pivot_longer(cols = c(weight_g, length_standard_mm), 
               names_to = "variable_growth", 
               values_to = "values")

# Comparisons -------------------------------------------------------------
# 1) Males vs females
# General stats
growth_M_F <- 
  data_growth %>% 
  group_by(sex, variable_growth) %>%
  summarise(mean = mean(values), 
            sd = sd(values))

data_growth %>% 
  group_by(variable_growth) %>% 
  shapiro_test(values)
# Normality not met (p < 0.05), but sample size that high that T test is preferred

# Comparison
test_M_F <- 
  data_growth %>% 
  group_by(variable_growth) %>%
  t_test(values ~ sex) %>% 
  add_significance()

# comp_M_F <- 
comp_M_F <- left_join(growth_M_F, test_M_F) 

comp_M_F %>% 
  arrange(variable_growth) %>%
  select(-c(".y.", group1, group2, p.signif, n1, n2)) %>% 
  gt(groupname_col = "variable_growth") %>% 
  fmt_scientific(columns = c(statistic, p)) %>% 
  fmt_number(columns = c(mean, sd), 
             decimals = 1) 

# 2) C vs T
# General stats
growth_C_T <- 
  data_growth %>% 
  group_by(sex, variable_growth, temp_f2) %>%
  summarise(mean = mean(values), 
            sd = sd(values)) %>% 
  arrange(sex)

# T_test
test_C_T <- 
  data_growth %>% 
  group_by(sex, variable_growth) %>%
  t_test(values ~ temp_f2) %>% 
  add_significance() %>% 
  select(-".y.") %>% 
  arrange(sex)

# Combining general stats with test
comp_C_T <- left_join(growth_C_T, test_C_T) 
comp_C_T %>% 
  select(-c(group1, group2, n1, n2, p.signif)) %>% 
  gt(groupname_col = "variable_growth") %>% 
  fmt_scientific(columns = c(statistic)) %>% 
  fmt_auto(columns = c(p)) %>% 
  fmt_number(columns = c(mean, sd), 
             decimals = 1) %>% 
  tab_style(style = cell_text(weight = "bold"), 
            locations = cells_body(columns = p, 
                                   rows = p < 0.05))

# 3) Comparing F0 specific size difference to temp f2 
# General stats
growth_f0_C_T <- 
  data_growth %>% 
  group_by(line_f0, sex, variable_growth, temp_f2) %>%
  summarise(mean = mean(values), 
            sd = sd(values)) %>% 
  group_by(line_f0, sex, variable_growth) %>%
  mutate(diff = round(diff(mean), 1)) %>% 
  ungroup() %>% 
  arrange(sex)

# T_test Assuming that normality is also not met, but sample size large enough to make meaningful comparisons
test_f0_C_T <- 
  data_growth %>% 
  group_by(line_f0, sex, variable_growth) %>%
  t_test(values ~ temp_f2) %>% 
  add_significance() %>% 
  select(-".y.") %>% 
  arrange(sex)

# Combining general stats with test
comp_f0_C_T <- left_join(growth_f0_C_T, test_f0_C_T) 
comp_f0_C_T %>% 
  select(-c(group1, group2, n1, n2, p.signif)) %>% 
  gt(groupname_col = "variable_growth") %>% 
  fmt_scientific(columns = c(statistic)) %>% 
  fmt_auto(columns = c(p)) %>% 
  fmt_number(columns = c(mean, sd, df), 
             decimals = 1) %>% 
  tab_style(style = cell_text(weight = "bold"), 
            locations = cells_body(columns = p, 
                                   rows = p < 0.05)) %>% 
  table_styles()

# ANOVA ------------------------------------------------------------------
# Including all categorical variables with interaction effects
# No interaction terms were found between temp and temp_f1, hence was excluded
anova_general <- data_growth %>% 
  group_by(variable_growth) %>% 
  anova_test(values ~ sex * temp_f2 * line_f0 + temp_f1, type = "3") %>% 
  as_tibble() %>% 
  # Make Effect rows more readable
  mutate(Effect = case_when(Effect == "sex" ~ "Sex",
                            Effect == "temp_f2" ~ "Temp F2",
                            Effect == "temp_f1" ~ "Temp F1",
                            Effect == "sex:temp_f2" ~ "Sex:Temp F2",
                            Effect == "line_f0" ~ "F0-line",
                            Effect == "temp_f2:line_f0" ~ "Temp F2:F0-line",
                            Effect == "sex:temp_f2:line_f0" ~ "Sex:Temp F2:F0-line", 
                            Effect == "sex:line_f0" ~ "Sex:F0-line"))

# Here a type III ANOVA is performed. rstatix::anova_test is a wrapper for the car::Anova() function
# No interaction terms were found between temp and temp_f1, hence was excluded

# ANOVA
anova_bw_by_sex <- data_growth %>% 
  group_by(line_f0, variable_growth) %>%
  anova_test(values ~ sex * temp_f2 + temp_f1, type = "3") %>% 
  as_tibble() %>% 
  # Make Effect rows more readable
  mutate(Effect = case_when(Effect == "sex" ~ "Sex",
                            Effect == "temp_f2" ~ "Temp F2",
                            Effect == "temp_f1" ~ "Temp F1",
                            Effect == "sex:temp_f2" ~ "Sex:Temp F2"))

# Anova results
# Weight
anova_weight <- 
  anova_bw_by_sex %>% 
  filter(variable_growth == "weight_g") %>% 
  select(-variable_growth)

# Length
anova_length <- 
  anova_bw_by_sex %>% 
  filter(variable_growth == "length_standard_mm") %>% 
  select(-variable_growth)


# ANOVA tables ------------------------------------------------------------

## Creating table style ----------------------------------------------------
anova_to_table <- function(table) {
  table %>% 
    select(-`p<.05`) %>%
    gt(groupname_col = "line_f0") %>% 
    cols_label(Effect = "Effect",
               DFn = html("df<sub>num</sub>"),
               DFd = html("df<sub>den</sub>"),
               F = "F-value",
               p = "p-value",
               ges = html("η<sup>2</sup><sub>g</sub>")) %>% 
    fmt_scientific(columns = c(ges)) %>% 
    fmt_auto(columns = c(p, F)) %>% 
    tab_style(style = cell_text(weight = "bold", align = "center"),
              locations = cells_row_groups()) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_body(columns = p,
                                     rows = p < 0.05)) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels()) %>%
    tab_style(style = cell_fill(color = "grey90"), 
              locations = cells_row_groups()) %>% 
    opt_table_font(font = "Times New Roman", 
                   size = 10, 
                   color = "black") %>% 
    opt_table_lines("none") %>%
    opt_vertical_padding(scale = .1) %>% 
    opt_horizontal_padding(scale = 3)
}


## Creating tables ---------------------------------------------------------
# General ANOVA
# Weight
tab_gen_anova_weigth <- 
  anova_general %>% 
  filter(variable_growth == "weight_g") %>% 
  select(-variable_growth) %>% 
  anova_to_table()

# Saving table
gtsave(data = tab_gen_anova_weigth, filename = "00_data/01_results/02_tab_gen_anova_weight.html")

# Length
tab_gen_anova_length <- 
  anova_general %>% 
  filter(variable_growth == "length_standard_mm") %>% 
  select(-variable_growth) %>% 
  anova_to_table()

# Saving table
gtsave(data = tab_gen_anova_length, filename = "00_data/01_results/02_tab_gen_anova_length.html")

# Grandsire specific ANOVA
# Weight
tab_anova_weight <- 
  anova_weight %>% anova_to_table()

# Saving table
gtsave(data = tab_anova_weight, filename = "00_data/01_results/02_tab_anova_weight.html")

# Length
tab_anova_length <- 
  anova_length %>% anova_to_table()

# Saving table
gtsave(data = tab_anova_length, filename = "00_data/01_results/02_tab_anova_length.html")

### ANOVA per sex and F0 line
# using type II ANOVA due to missing S2 tank (unbalanced data)
anova_specific <- data_growth %>% 
  group_by(variable_growth, line_f0, sex) %>% 
  anova_test(values ~ temp_f2 * temp_f1, type = "2") %>% 
  as_tibble() %>% 
  # Make Effect rows more readable
  mutate(
    Effect = case_when(Effect == "sex" ~ "Sex",
                       Effect == "temp_f2" ~ "Temp F2",
                       Effect == "temp_f1" ~ "Temp F1",
                       Effect == "temp_f2:temp_f1" ~ "Temp F2:Temp F1"), 
    variable_growth = ifelse(variable_growth == "weight_g", "Body weight (g)", "Standard length (mm)")
  )

t_anova_specific <- 
  anova_specific %>% 
  select(-`p<.05`) %>%
  gt(groupname_col = "line_f0") %>% 
  cols_label(Effect = "Effect",
             sex = "Sex",
             variable_growth = "Variable",
             DFn = html("df<sub>num</sub>"),
             DFd = html("df<sub>den</sub>"),
             F = "F-value",
             p = "p-value",
             ges = html("η<sup>2</sup><sub>g</sub>")) %>% 
  fmt_scientific(columns = c(ges)) %>% 
  fmt_auto(columns = c(p, F)) %>% 
  tab_style(style = cell_text(weight = "bold", align = "center"),
            locations = cells_row_groups()) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = p,
                                   rows = p < 0.05)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) %>%
  tab_style(style = cell_fill(color = "grey90"), 
            locations = cells_row_groups()) %>% 
  opt_table_font(font = "Times New Roman", 
                 size = 10, 
                 color = "black") %>% 
  opt_table_lines("none") %>%
  opt_vertical_padding(scale = .1) %>% 
  opt_horizontal_padding(scale = 3) 

gtsave(data = t_anova_specific, filename = "03_tables/t_anova_specific.html")
t_anova_specific %>% gtsave("03_tables/t_anova_specific.docx")
