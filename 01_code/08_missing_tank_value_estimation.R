# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: October 2025
# Description: Missing tank value imputation

# Settings ----------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(gt)
library(gtsummary)
library(gtExtras)

# Loading data ------------------------------------------------------------

# Loading main data
load("00_data/data_main.RData")

# Loading color palettes
load("00_data/color_palettes.RData")

# Summarizing key statistics ----------------------------------------------
# Simple overview
sex_ratio_summary <- data_main %>% 
  group_by(treatment, temp_f2) %>% 
  summarise(n = n(), n_female = sum(sex == "F"), ratio_female = sum(sex == "F")/n,
            n_male = sum(sex == "M"), ratio_male = sum(sex == "M")/n, .groups = "drop") %>% 
  arrange(treatment)

# Sex ratio more detailed
sex_ratio_summary_det <- data_main %>% 
  group_by(line_f0, treatment, temp_f0, temp_f1, temp_f2) %>%
  summarise(n = n(), n_female = sum(sex == "F"), ratio_female = sum(sex == "F")/n, 
            n_male = sum(sex == "M"), ratio_male = sum(sex == "M")/n, .groups = "drop") %>%
  arrange(treatment)

# Checking that male proportions are normally distributed
sex_ratio_summary_det %>% 
  select(ratio_male) %>% 
  shapiro_test(ratio_male)

shapiro.test(sex_ratio_summary_det %>% pull(ratio_male))

## Adding estimation value (missing tank S2-CTT)
# Missing tank with equal sex ratio
tank_missing <- sex_ratio_summary_det[6,] %>% 
  mutate(treatment = "CTT", 
         temp_f2 = "T", 
         n_female = 54, 
         n_male = 46, 
         ratio_female = 0.54, 
         ratio_male = 0.46)

sex_ratio_summary_det <- rbind(sex_ratio_summary_det,tank_missing) %>% 
  mutate(code_2 = paste0(temp_f1, temp_f2), 
         treatment = paste0(temp_f0, temp_f1, temp_f2))

# Sex ratio comparison for known data
sex_ratio_diff <- sex_ratio_summary_det %>% 
  group_by(temp_f2, temp_f1, line_f0) %>% 
  filter(code_2 %in% c("TC", "TT")) %>% 
  group_by(line_f0) %>% 
  mutate(ratio_diff = diff(ratio_male)) %>% 
  select(line_f0, ratio_diff) %>% 
  unique() %>% 
  filter(line_f0 != "S2") %>% 
  ungroup() %>% 
  mutate(mean = mean(ratio_diff), 
         sd = sd(ratio_diff))

# Male ratio difference (masculinisation) from TC to TT: 0.244, 0.225, 0.306

sex_ratio_mean <-  sex_ratio_diff$mean %>% unique
sex_ratio_sd <-  sex_ratio_diff$sd %>% unique

# Plotting sex ratio per sire with different standard deviations for estimated value
p_sex_ratio_per_sire <- 
  ggplot(sex_ratio_summary_det, aes(x = treatment, y = ratio_male, fill = treatment)) +
  stat_summary(fun = mean,
               geom = "col") +
  geom_errorbar(data = . %>% filter(treatment == "CTT", line_f0 == "S2"), 
                aes(x = treatment, 
                    ymin = 0.2 + sex_ratio_mean - sex_ratio_sd*3,
                    ymax = 0.2 + sex_ratio_mean + sex_ratio_sd*3), 
                width = .5, 
                colour = "grey70") +
  geom_errorbar(data = . %>% filter(treatment == "CTT", line_f0 == "S2"), 
                aes(x = treatment, 
                    ymin = 0.2 + sex_ratio_mean - sex_ratio_sd*2,
                    ymax = 0.2 + sex_ratio_mean + sex_ratio_sd*2), 
                width = .35, 
                colour = "grey40") +
  geom_errorbar(data = . %>% filter(treatment == "CTT", line_f0 == "S2"), 
                aes(x = treatment, 
                    ymin = 0.2 + sex_ratio_mean - sex_ratio_sd,
                    ymax = 0.2 + sex_ratio_mean + sex_ratio_sd), 
                width = .25) +
  facet_grid(~ line_f0, 
             scales = "free_x",
             space = "free") +
  theme_pubr() +
  labs(x = "Temperature history (F0, F1, F2)", 
       y = "F2 male proportion", 
       fill = "Temperature") +
  scale_fill_manual(values = c("CCC" = "royalblue", "CCT" = "#ad52b5", "CTC" = "#da4c72", "CTT" = "orangered", 
                               "TCC" = "royalblue", "TCT" = "#ad52b5", "TTC" = "#da4c72", "TTT" = "orangered")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = "none", 
        strip.background = element_rect(linewidth = NA), 
        strip.text = element_text(face = "bold"), 
        axis.text.x = element_text(size = 9))

ggsave("02_plots/08_sex_ratio_per_sire.pdf", p_sex_ratio_per_sire, width = 6, height = 5)
ggsave("02_plots/08_sex_ratio_per_sire.png", p_sex_ratio_per_sire, width = 6, height = 5)

# Plotting sex ratio effect plots with uncertainty
# F2 effect plot
# Getting mean value for T
sex_ratio_mean_f2 <- sex_ratio_summary_det %>% 
  filter(line_f0 == "S2", temp_f2 == "T") %>% 
  summarise(ratio_mean = mean(ratio_male)) %>% as.numeric()

p_developmental_effect <- ggplot(sex_ratio_summary_det, 
       aes(x = temp_f2, y = ratio_male, group = line_f0, color = line_f0)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, show.legend = FALSE) +
  geom_errorbar(data = . %>% filter(temp_f2 == "T", line_f0 == "S2"), 
                aes(x = temp_f2, 
                    ymin = sex_ratio_mean_f2 - sex_ratio_sd*1.5,
                    ymax = sex_ratio_mean_f2 + sex_ratio_sd*1.5), 
                width = .15, 
                colour = "grey70") +
  geom_errorbar(data = . %>% filter(temp_f2 == "T", line_f0 == "S2"), 
                aes(x = temp_f2, 
                    ymin = sex_ratio_mean_f2 - sex_ratio_sd,
                    ymax = sex_ratio_mean_f2 + sex_ratio_sd), 
                width = .125, 
                colour = "grey40") +
  geom_errorbar(data = . %>% filter(temp_f2 == "T", line_f0 == "S2"), 
                aes(x = temp_f2, 
                    ymin = sex_ratio_mean_f2 - sex_ratio_sd/2,
                    ymax = sex_ratio_mean_f2 + sex_ratio_sd/2), 
                width = .1, 
                colour = "black") +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(x = "Temperature F2", 
       y = "F2 Male proportion",
       color = "F0 line", 
       title = "Developmental effect") +
  scale_color_manual(values = color_grand) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_pubr() +
  theme(legend.position = "right",
        legend.direction = "vertical", 
        panel.grid.major.y = element_line(), 
        plot.title = element_text(face = "bold", hjust = .5))

ggsave("02_plots/08_developmental_effect.pdf", p_developmental_effect, width = 6, height = 5)
ggsave("02_plots/08_developmental_effect.png", p_developmental_effect, width = 6, height = 5)

# F1 effect plot
# Getting mean value for T
sex_ratio_mean_f1 <- sex_ratio_summary_det %>% 
  filter(line_f0 == "S2", temp_f1 == "T") %>% 
  summarise(ratio_mean = mean(ratio_male)) %>% as.numeric()

p_parental_effect <- ggplot(sex_ratio_summary_det,
       aes(x = temp_f1, y = ratio_male, group = line_f0, color = line_f0)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, show.legend = FALSE) +
  geom_errorbar(data = . %>% filter(temp_f1 == "T", line_f0 == "S2"), 
                aes(x = temp_f1, 
                    ymin = sex_ratio_mean_f1 - sex_ratio_sd*1.5,
                    ymax = sex_ratio_mean_f1 + sex_ratio_sd*1.5), 
                width = .15, 
                colour = "grey70") +
  geom_errorbar(data = . %>% filter(temp_f1 == "T", line_f0 == "S2"), 
                aes(x = temp_f1, 
                    ymin = sex_ratio_mean_f1 - sex_ratio_sd,
                    ymax = sex_ratio_mean_f1 + sex_ratio_sd), 
                width = .125, 
                colour = "grey40") +
  geom_errorbar(data = . %>% filter(temp_f1 == "T", line_f0 == "S2"), 
                aes(x = temp_f1, 
                    ymin = sex_ratio_mean_f1 - sex_ratio_sd/2,
                    ymax = sex_ratio_mean_f1 + sex_ratio_sd/2), 
                width = .1, 
                colour = "black") +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(x = "Temperature F1", 
       y = "F2 Male proportion", 
       color = "F0 line",
       title = "Parental effect") +
  scale_color_manual(values = color_grand) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_pubr() +
  theme(legend.position = "right",
        legend.direction = "vertical", 
        panel.grid.major.y = element_line(), 
        plot.title = element_text(face = "bold", hjust = .5))

ggsave("02_plots/08_parental_effect.pdf", p_parental_effect, width = 6, height = 5)
ggsave("02_plots/08_parental_effect.png", p_parental_effect, width = 6, height = 5)

# Sex proportion plot
sex_prop_f1_alt <- sex_ratio_summary_det %>%
  group_by(line_f0, temp_f1, temp_f2) %>%
  summarise(sex_prop = mean(ratio_male, na.rm = TRUE),
            n = sum(n),
            n_male = sum(n_male),
            n_female = sum(n_female),
            .groups = "drop") %>%
  pivot_wider(names_from = temp_f1,
              values_from = c(sex_prop, n, n_male, n_female),
              names_glue = "{.value}_{temp_f1}") %>%
  mutate(sex_prop_diff = sex_prop_T - sex_prop_C)

# Adding type column for "real" values and simulated ones
sex_prop_f1_alt$type <- "observed"  
sex_prop_f1_alt[4, "type"] <- "estimated"

# Converting type column into ordered factor
sex_prop_f1_alt$type <- factor(sex_prop_f1_alt$type, levels = c("observed", "estimated"))

# Variables
an_x_pos_w_est <- 0.6
an_y_pos_w_est <- 1.9
an_y_increment <- 0.3
an_text_size <- 3

fill_w_est <- "green4"
label_alpha <- 0.3
text_color <- "black"

# Plot
p_effect_plot <-
  ggplot(sex_prop_f1_alt,
         aes(
           x = qnorm(sex_prop_C),
           y = qnorm(sex_prop_T) - qnorm(sex_prop_C),
           color = line_f0
         )) +
  geom_abline(
    slope = 0,
    intercept = 0,
    linewidth = .5,
    linetype = 3,
    color = "black",
    alpha = .5
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = .5,
    linetype = 3,
    color = "black",
    alpha = .5
  ) +
  # For data with estimate
  geom_smooth(
    data = sex_prop_f1_alt,
    color = fill_w_est,
    method = "lm",
    fill = fill_w_est,
    linewidth = 1,
    alpha = .1,
    show.legend = FALSE
  ) +
  stat_regline_equation(
    aes(label = after_stat(eq.label)),
    formula = y ~ x,
    label.x = an_x_pos_w_est,
    label.y = an_y_pos_w_est - 2 * an_y_increment,
    size = an_text_size,
    color = "black"
  ) +
  stat_cor(
    aes(label = after_stat(r.label)),
    method = "pearson",
    label.x = an_x_pos_w_est,
    label.y = an_y_pos_w_est - 3 * an_y_increment,
    color = "black",
    size = an_text_size
  ) +
  stat_cor(
    aes(label = after_stat(p.label)),
    method = "pearson",
    label.x = an_x_pos_w_est,
    label.y = an_y_pos_w_est - 4 * an_y_increment,
    color = "black",
    digits = 3,
    size = an_text_size
  ) +
  geom_errorbar(data = . %>% filter(temp_f2 == "T", line_f0 == "S2"),
                 aes(x = qnorm(sex_prop_C),
                     ymin = qnorm(sex_prop_T - sex_ratio_sd*3) - qnorm(sex_prop_C),
                     ymax = qnorm(sex_prop_T + sex_ratio_sd*3) - qnorm(sex_prop_C)), 
                 width = .05, 
                 colour = "grey70") +
  geom_errorbar(data = . %>% filter(temp_f2 == "T", line_f0 == "S2"),
                 aes(x = qnorm(sex_prop_C),
                     ymin = qnorm(sex_prop_T - sex_ratio_sd*2) - qnorm(sex_prop_C),
                     ymax = qnorm(sex_prop_T + sex_ratio_sd*2) - qnorm(sex_prop_C)), 
                width = .025, 
                 colour = "grey40") +
  geom_errorbar(data = . %>% filter(temp_f2 == "T", line_f0 == "S2"),
                 aes(x = qnorm(sex_prop_C),
                     ymin = qnorm(sex_prop_T - sex_ratio_sd) - qnorm(sex_prop_C),
                     ymax = qnorm(sex_prop_T + sex_ratio_sd) - qnorm(sex_prop_C)), 
                width = .015, 
                 colour = "black") +
  geom_point(
    data = sex_prop_f1_alt,
    aes(shape = type),
    show.legend = TRUE,
    size = 2
  ) +
  theme_pubr() +
  labs(
    x = "Probit F2 male proportion (F1 XCX)",
    y = "Probit F2 male proportion delta (F1 XTX - F1 XCX)",
    color = "F0 line",
    shape = "Type",
  ) +
  scale_color_manual(values = color_grand) +
  scale_shape_manual(values = c("observed" = 16, "estimated" = 4)) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    plot.title = element_text(face = "bold", hjust = .5)
  )

p_effect_plot

ggsave("02_plots/08_effect_plot_f1_with_errorbar.pdf", p_effect_plot, width = 6, height = 5)
ggsave("02_plots/08_effect_plot_f1_with_errorbar.png", p_effect_plot, width = 6, height = 5)

## Creating df with range of possible male proportions
# New df
df_sex_prop <- sex_prop_f1_alt %>% 
  select(line_f0, temp_f2, sex_prop_C, sex_prop_T, type)

# Range of possible values up to 3 times sd
male_prop_values <- seq(from = tank_missing$ratio_male - 3*sex_ratio_sd, 
                        to = tank_missing$ratio_male + 3*sex_ratio_sd,
                        by = sex_ratio_sd)

group_names <- c("-3*sd", "-2*sd", "-1*sd", "mean", "+1*sd", "+2*sd", "+3*sd") 

# Creating new values in long format
df_sex_prop_it <- 
  tibble(group_id = group_names,
         new_sex_prop_T = male_prop_values) %>%
  rowwise() %>%
  mutate(modified_data = list({
      df_temp <- df_sex_prop
      # Adding possible estimate values
      df_temp$sex_prop_T[4] <- new_sex_prop_T
      df_temp
    })
  ) %>%
  ungroup() %>% 
  unnest(modified_data) %>%
  rename(group_name = group_id) %>% 
  select(-new_sex_prop_T) %>% 
  mutate(
    # adding normal difference 
    sex_prop_diff = sex_prop_T - sex_prop_C, 
    # probit transforming for testing correlation
    probit_sr_C = qnorm(sex_prop_C), 
    probit_sr_T = qnorm(sex_prop_T), 
    probit_sr_diff = qnorm(sex_prop_T) - qnorm(sex_prop_C)
         )

# Testing correlation for all estimated values
cor_test_estimates <- 
  df_sex_prop_it %>% 
  group_by(group_name) %>% 
  cor_test(probit_sr_C, probit_sr_diff) %>% 
  select(-c(var1, var2, conf.low, conf.high, method)) %>% 
  arrange(match(group_name, group_names)) %>% 
  mutate(est_male_prop = male_prop_values, .after = group_name)

## Combined plots
# Variables
an_x_pos_w_est <- 0.3
an_y_pos_w_est <- 2
an_y_increment <- 0.2
an_text_size <- 3

# Plotting all possible plots 
p_effect_plot_all_sd <-
  ggplot(df_sex_prop_it,
         aes(
           x = probit_sr_C,
           y = probit_sr_diff,
           color = line_f0
         )) +
  geom_abline(
    slope = 0,
    intercept = 0,
    linewidth = .5,
    linetype = 3,
    color = "black",
    alpha = .5
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = .5,
    linetype = 3,
    color = "black",
    alpha = .5
  ) +
  # For data with estimate
  geom_smooth(
    data = df_sex_prop_it,
    color = fill_w_est,
    method = "lm",
    fill = fill_w_est,
    linewidth = 1,
    alpha = .1,
    show.legend = FALSE
  ) +
  stat_regline_equation(
    aes(label = after_stat(eq.label)),
    formula = y ~ x,
    label.x = an_x_pos_w_est,
    label.y = an_y_pos_w_est - 2 * an_y_increment,
    size = an_text_size,
    color = "black"
  ) +
  stat_cor(
    aes(label = after_stat(r.label)),
    method = "pearson",
    label.x = an_x_pos_w_est,
    label.y = an_y_pos_w_est - 3 * an_y_increment,
    color = "black",
    size = an_text_size
  ) +
  stat_cor(
    aes(label = after_stat(p.label)),
    method = "pearson",
    label.x = an_x_pos_w_est,
    label.y = an_y_pos_w_est - 4 * an_y_increment,
    color = "black",
    digits = 3,
    size = an_text_size
  ) +
  facet_wrap(~ group_name) +
  geom_errorbar(data = df_sex_prop_it %>% filter(temp_f2 == "T", line_f0 == "S2"),
                aes(x = probit_sr_C,
                    ymin = qnorm(tank_missing$ratio_male - sex_ratio_sd*3) - probit_sr_C,
                    ymax = qnorm(tank_missing$ratio_male + sex_ratio_sd*3) - probit_sr_C), 
                width = .05, 
                colour = "grey70") +
  geom_errorbar(data = df_sex_prop_it %>% filter(temp_f2 == "T", line_f0 == "S2"),
                aes(x = probit_sr_C,
                    ymin = qnorm(tank_missing$ratio_male - sex_ratio_sd*2) - probit_sr_C,
                    ymax = qnorm(tank_missing$ratio_male + sex_ratio_sd*2) - probit_sr_C), 
                width = .025, 
                colour = "grey40") +
  geom_errorbar(data = df_sex_prop_it %>% filter(temp_f2 == "T", line_f0 == "S2"),
                aes(x = probit_sr_C,
                    ymin = qnorm(tank_missing$ratio_male - sex_ratio_sd) - probit_sr_C,
                    ymax = qnorm(tank_missing$ratio_male + sex_ratio_sd) - probit_sr_C), 
                width = .015, 
                colour = "black") +
  geom_point(
    aes(shape = type),
    show.legend = TRUE,
    size = 2
  ) +
  theme_pubr() +
  labs(
    x = "Probit F2 male proportion (XCX)",
    y = "Probit F2 male proportion delta (XTX - XCX)",
    color = "F0 line",
    shape = "Type",
  ) +
  scale_color_manual(values = color_grand) +
  scale_shape_manual(values = c("observed" = 16, "estimated" = 4)) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    plot.title = element_text(face = "bold", hjust = .5)
  )

p_effect_plot_all_sd

ggsave("02_plots/08_effect_plot_f1_all_sd.pdf", p_effect_plot_all_sd, width = 10, height = 10)
ggsave("02_plots/08_effect_plot_f1_all_sd.png", p_effect_plot_all_sd, width = 10, height = 10)

p_summary_plot <- 
  p_sex_ratio_per_sire + p_developmental_effect + p_parental_effect + p_effect_plot +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", 
                                size = 18), 
        plot.title = element_blank(), 
        legend.direction = "vertical")

ggsave("02_plots/08_summary_plot.pdf", p_summary_plot, width = 12, height = 10)
ggsave("02_plots/08_summary_plot.png", p_summary_plot, width = 12, height = 10)

p_effects <- 
  p_developmental_effect + p_parental_effect +
  plot_layout(guides = "collect", axes = "collect_y") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", 
                                size = 18), 
        plot.title = element_blank(), 
        legend.direction = "vertical")

ggsave("02_plots/08_effects.pdf", p_effects, width = 8, height = 5)
ggsave("02_plots/08_effects.png", p_effects, width = 8, height = 5)



# Tables ------------------------------------------------------------------

# Creating tables for correlation tests
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

# Correlation tests for estimates
tab_cor_test <- 
  cor_test_estimates %>% 
  gt() %>% 
  cols_label(group_name = "Estimate", 
             est_male_prop = "Estimate value", 
             cor = "R", 
             statistic = "Statistic", 
             p = "p-value") %>% 
  fmt_number(columns = c(est_male_prop, statistic), 
             decimals = 3) %>% 
  table_styles

# Saving table
gtsave(data = tab_cor_test, filename = "00_data/01_results/08_cor_test_estimate.html")
gtsave(data = tab_cor_test, filename = "03_tables/08_cor_test_estimate.docx")
