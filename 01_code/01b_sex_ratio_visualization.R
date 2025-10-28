# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: October 2025
# Description: Sex ratio analysis  

# Settings ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(ggpubr)
library(legendry)
library(car)
library(lmtest)

# Loading data ------------------------------------------------------------

# Main data
load("00_data/data_main.RData")

# Loading color palettes
load("00_data/color_palettes.RData")

# Data preparation -----------------------------------------------------------------

# Simple overview
sex_ratio_summary <- data_main %>% 
  group_by(treatment, temp_f2) %>%
  summarise(
    n = n(),
    n_female = sum(sex == "F"),
    ratio_female = sum(sex == "F") / n,
    n_male = sum(sex == "M"),
    ratio_male = sum(sex == "M") / n,
    .groups = "drop"
  ) %>%
  arrange(treatment)
  
write.csv(sex_ratio_summary, "00_data/01_results/01_sex_ratio_summary.csv", row.names = FALSE)

# Sex ratio more detailed
sex_ratio_summary_det <- data_main %>%
  group_by(line_f0, treatment, temp_f0, temp_f1, temp_f2) %>%
  summarise(
    n = n(),
    n_female = sum(sex == "F"),
    ratio_female = sum(sex == "F") / n,
    n_male = sum(sex == "M"),
    ratio_male = sum(sex == "M") / n,
    .groups = "drop"
  ) %>%
  arrange(treatment)

write.csv(sex_ratio_summary_det, "00_data/01_results/01_sex_ratio_summary_det.csv", row.names = FALSE)

## Adding estimation value (missing tank S2-TT)
# See explanation: "08_missing_tank_value_estimation.R"
# Missing tank
tank_missing <- sex_ratio_summary_det[6,] %>% 
  mutate(temp_f2 = "T", 
         n_female = 54, 
         n_male = 46, 
         ratio_female = 0.54, 
         ratio_male = 0.46)

sex_ratio_summary_det_w_est <- rbind(sex_ratio_summary_det, tank_missing)
  
# Adding SE and mean ratio (with estimated value)
sex_ratio_summary_det_2 <- 
  sex_ratio_summary_det %>% 
  ungroup() %>% 
  group_by(treatment, temp_f2, temp_f1, temp_f0) %>% 
  mutate(mean_ratio_male = mean(n_male/n), sd_ratio_male = sd(ratio_male, na.rm = TRUE)) %>% 
  replace_na(list(sd_ratio_male = 0)) %>% 
  mutate(se = sd_ratio_male / sqrt(2))

# Adding SE and mean ratio and group by line_f0 (with estimated value)
sex_ratio_summary_det_3 <- 
  sex_ratio_summary_det %>% 
  ungroup() %>% 
  group_by(line_f0) %>% 
  mutate(
    mean_ratio_male = mean(n_male/n), 
    sd_ratio_male = sd(ratio_male, na.rm = TRUE)
    ) %>% 
  replace_na(list(sd_ratio_male = 0)) %>% 
  mutate(se = sd_ratio_male / sqrt(2))

sex_ratio_summary_det_3_w_est <- 
  sex_ratio_summary_det_w_est %>% 
  ungroup() %>% 
  group_by(line_f0) %>% 
  mutate(
    mean_ratio_male = mean(n_male/n), 
    sd_ratio_male = sd(ratio_male, na.rm = TRUE)
  ) %>% 
  replace_na(list(sd_ratio_male = 0)) %>% 
  mutate(se = sd_ratio_male / sqrt(2))

# Plots -------------------------------------------------------------------

# Fig. 1a -----------------------------------------------------------------

p_sex_col <- 
  ggplot(
    sex_ratio_summary_det_2, 
    aes(x = interaction(temp_f2, temp_f1, temp_f0), 
        y = ratio_male, 
        fill = temp_f2)
  ) +
  stat_summary(fun.data = mean_se, 
               geom = "errorbar",
               width = 0.2) +
  stat_summary(fun = mean, 
               geom = "col") +
  geom_abline(slope = 0, 
              intercept = .5, 
              linewidth = .5, 
              linetype = 3, 
              color = "black") +
  labs(y = "F2 male proportion", 
       x = "Temperature history",
       fill = "F2 temp") +
  scale_fill_manual(values = color_temp) +
  theme_pubr() +
  scale_x_discrete(guide = guide_axis_nested (title = "Temperature history")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(x = guide_axis_nested(pad_discrete = .1))
p_sex_col

ggsave("02_plots/01_sex_overview_error.pdf", p_sex_col, width = 6, height = 5)


# Fig. 1b -----------------------------------------------------------------

# Data preparation
sex_ratio_summary_overview <- 
  sex_ratio_summary_det %>% 
  group_by(temp_f2) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(ratio_male = n_male/n, 
         ratio_female = n_female/n) %>% 
  pivot_longer(cols = c(n_male, n_female), names_to = "n_sex") %>% 
  mutate(sex = ifelse(n_sex == "n_male", "M", "F"), 
         pos_text = ifelse(sex == "F", ratio_male/2, ratio_male/4), 
         sex_full = ifelse(sex == "F", "Females", "Males")) %>% 
  mutate(sex = factor(sex, levels = c("F", "M")))


p_sex_summary_temp_f2 <- 
  ggplot(sex_ratio_summary_overview, aes(x = temp_f2, y = value, fill = sex)) +
  geom_col(position = "fill") + 
  geom_abline(slope = 0, intercept = .5, linewidth = .5, linetype = 3, color = "black") +
  geom_text(aes(label = paste0("n = ", value), fontface = 1),
            color = "black", position = position_fill(vjust = 0.5), check_overlap = TRUE) +
  geom_text(aes(label = paste0("Total = ", n), y = 1.00, fontface = 2),
            color = "black",
            position = "identity", 
            vjust = -.4,
            check_overlap = TRUE) +
  labs(y = "F2 sex proportion", 
       fill = "Sex", 
       x = "F2 temperature") +
  scale_fill_manual(values = color_sex) +
  theme_pubr() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)), breaks = c(0,.2,.4,.6,.8,1))
p_sex_summary_temp_f2

ggsave("02_plots/01_sex_summary_temp_f2.pdf", p_sex_summary_temp_f2, width = 4.5, height = 4.5)
ggsave("02_plots/01_sex_summary_temp_f2.png", p_sex_summary_temp_f2, width = 4.5, height = 4.5)

# Fig. 1c -----------------------------------------------------------------

# Point range plot for different F0 lines
p_sex_pointrange <- 
  ggplot(sex_ratio_summary_det_3, aes(x = line_f0, color = line_f0)) +
  geom_point(aes(y = mean_ratio_male), size = 3, show.legend = FALSE) +
  geom_errorbar(
    aes(y = mean_ratio_male, 
        ymax = mean_ratio_male+se, 
        ymin = mean_ratio_male-se), 
    width = .2,
    show.legend = FALSE
  ) +
  labs(y = "F2 male proportion", 
       x = "F0 line") +
  ylim(c(0.2, 0.83)) +
  theme_pubr()

ggsave("02_plots/01_sex_pointrange.pdf", p_sex_pointrange, width = 4.5, height = 4.5)
ggsave("02_plots/01_sex_pointrange_bigger.pdf", p_sex_pointrange, width = 4, height = 4)


# Fig. 1 ------------------------------------------------------------------

# Final version manuscript
plot_spacer() + p_sex_summary_temp_f2 + p_sex_col + free(p_sex_pointrange, type = "panel", side = "b") +
  plot_layout(guides = "collect", axis_titles = "collect_x", widths = c(3,2)) &
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 22), legend.direction = "vertical")

ggsave("02_plots/fig_1.pdf", width = 9, height = 7)

# Fig. 2a -----------------------------------------------------------------

# Adding codes to data
sex_ratio_summary_code <- 
  sex_ratio_summary_det_w_est %>%
  mutate(code_2 = paste0(temp_f1, temp_f2),
         treatment = paste0(temp_f0, temp_f1, temp_f2)) %>% 
  filter(line_f0 != "S2" | treatment != "CTT")

p_sex_by_sire <- 
  ggplot(sex_ratio_summary_code, aes(x = treatment, y = ratio_male, fill = code_2)) +
  stat_summary(fun = mean, geom = "col", show.legend = FALSE) +
  facet_grid(~ line_f0, scales = "free_x", space = "free") +
  theme_pubr() +
  labs(x = "Temperature regime (F0, F1, F2)", 
       y = "F2 male proportion", 
       fill = "Temperature") +
  scale_fill_manual(values = c("CC" = "royalblue", "CT" = "#ad52b5", "TC" = "#da4c72", "TT" = "orangered")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none", 
        strip.background = element_rect(linewidth = NA), 
        axis.text.x = element_text(size = 8))

ggsave("02_plots/01_sex_ratio_F0_sep_wo_est.pdf", p_sex_by_sire, width = 8, height = 6)


# Fig. 2b -----------------------------------------------------------------

# F2
p_sex_effect_F2 <- 
  ggplot(sex_ratio_summary_code, aes(x = temp_f2, y = ratio_male, group = line_f0, color = line_f0)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(x = "F2 temperature", 
       y = "F2 male proportion",
       color = "F0 line", 
       title = "Developmental effect") +
  scale_color_manual(values = color_grand) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_pubr() +
  theme(legend.position = "right",
        legend.direction = "vertical", 
        panel.grid.major.y = element_line(), 
        plot.title = element_text(face = "bold", hjust = .5))
p_sex_effect_F2

# F1
p_sex_effect_F1 <- 
  ggplot(sex_ratio_summary_code, aes(x = temp_f1, y = ratio_male, group = line_f0, color = line_f0)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  labs(x = "F1 temperature", 
       y = "F2 male proportion", 
       color = "F0 line",
       title = "Parental effect") +
  scale_color_manual(values = color_grand) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_pubr() +
  theme(legend.position = "right",
        legend.direction = "vertical", 
        panel.grid.major.y = element_line(), 
        plot.title = element_text(face = "bold", hjust = .5))
p_sex_effect_F1

# Summary plot
p_sex_effect_both <- p_sex_effect_F2 + p_sex_effect_F1 +
  plot_layout(axes = "collect_y", guides = "collect") & theme(legend.position =  "right")

ggsave("02_plots/01_sex_effect_both.pdf", p_sex_effect_both, width = 6, height = 5)
ggsave("02_plots/01_sex_effect_f2.pdf", p_sex_effect_F2, width = 4, height = 4.5)
ggsave("02_plots/01_sex_effect_f1.pdf", p_sex_effect_F1, width = 4, height = 4.5)

# Fig. 2c -----------------------------------------------------------------

# Sex proportions F1
data_probit_f1 <- sex_ratio_summary_det_3_w_est %>%
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

# Testing assumption for linear model 
model <- lm(qnorm(sex_prop_C)~(qnorm(sex_prop_T)-qnorm(sex_prop_C)), data_probit_f1)

shapiro.test(residuals(model))
bptest(model)
durbinWatsonTest(model)

# Adding type column for "real" values and simulated ones
data_probit_f1$type <- "observed"  
data_probit_f1[4, "type"] <- "estimated"

# Converting type column into ordered factor
data_probit_f1$type <- factor(data_probit_f1$type, levels = c("observed", "estimated"))

# Variables
an_x_pos_w_est <- 0.6
an_y_pos_w_est <- 1.9
an_x_pos <- -0.65
an_y_pos <- -0.5
an_y_increment <- 0.3
an_text_size <- 3

fill_w_est <- "green4"
fill_wo_est <- "orange"
label_alpha <- 0.3
text_color <- "black"

# Combined plot
p_probit_sex_prop_both <-
  ggplot(data_probit_f1,
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
    data = data_probit_f1,
    color = fill_w_est,
    method = "lm",
    fill = fill_w_est,
    linewidth = 1,
    alpha = .1,
    show.legend = FALSE
  ) +
  annotate(
    "label",
    label = "With estimate",
    x = an_x_pos_w_est,
    y = an_y_pos_w_est - 1 * an_y_increment,
    hjust = 0,
    fontface = "bold",
    size = an_text_size,
    fill = fill_w_est,
    alpha = label_alpha,
    color = text_color,
    label.size = NA
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
  # For data without estimate
  geom_smooth(
    data = data_probit_f1 %>% filter(type != "estimated"),
    aes(
      x = qnorm(sex_prop_C),
      y = qnorm(sex_prop_T) - qnorm(sex_prop_C)
    ),
    color = fill_wo_est,
    linetype = "dashed",
    method = "lm",
    fill = fill_wo_est,
    alpha = .1,
    show.legend = FALSE
  ) +
  annotate(
    "label",
    label = "Without estimate",
    x = an_x_pos,
    y = an_y_pos,
    hjust = 0,
    fontface = "bold",
    size = an_text_size,
    fill = fill_wo_est,
    alpha = label_alpha,
    color = text_color,
    label.size = NA
  ) +
  stat_regline_equation(
    data = data_probit_f1 %>% filter(type != "estimated"),
    aes(label = after_stat(eq.label)),
    formula = y ~ x,
    label.x = an_x_pos,
    label.y = an_y_pos - 1 * an_y_increment,
    color = "black",
    size = an_text_size
  ) +
  stat_cor(
    data = data_probit_f1 %>% filter(type != "estimated"),
    aes(label = after_stat(r.label)),
    method = "pearson",
    label.x = an_x_pos,
    label.y = an_y_pos - 2 * an_y_increment,
    color = "black",
    size = an_text_size
  ) +
  stat_cor(
    data = data_probit_f1 %>% filter(type != "estimated"),
    aes(label = after_stat(p.label)),
    method = "pearson",
    label.x = an_x_pos,
    label.y = an_y_pos - 3 * an_y_increment,
    color = "black",
    digits = 3,
    size = an_text_size
  ) +
  # For both
  geom_point(
    data = data_probit_f1,
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

p_probit_sex_prop_both

# Fig. 2 ------------------------------------------------------------------

# Summary plot of sex by sire and effect by sire
p_sex_effect_combined <-  p_sex_by_sire + p_sex_effect_both + p_probit_sex_prop_both +
  plot_layout(axis_titles = "collect_y", ncol = 1, guides = "collect") &
  theme(legend.position = "right",
        legend.direction = "vertical", 
        plot.title = element_blank())

# Legend is bugged, needs to be fixed after
ggsave("02_plots/fig_2.pdf", p_sex_effect_combined, width = 7, height = 10)
