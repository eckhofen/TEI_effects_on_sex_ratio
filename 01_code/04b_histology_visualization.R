# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: October 2025
# Description: Histological analyses

# Settings ----------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(MASS)
library(broom)

# Loading data ------------------------------------------------------------

# Loading in sample data
load("00_data/00_processed/data_histology.RData")

# Loading color palettes
load("00_data/color_palettes.RData")

# Plots -------------------------------------------------------------------

# Custom function to avoid lower error bars
mean_sd_upper <- function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  data.frame(ymin = 0, ymax = mean_val + sd_val)
}

# Creating ggplot style object 
ggstyle <- function(color,
                    x = "Temperature F2",
                    y = "Cell count",
                    fill = "Temperature F2",
                    x_axis_angle = 90) {
  list(
    facet_wrap(
      ~ cell_type,
      nrow = 1,
      strip.position = "top",
      scales = "free"
    ),
    scale_fill_manual(values = color),
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))),
    labs(x = x, y = y, fill = fill),
    theme_pubr(),
    theme(
      strip.background = element_rect(linewidth = NA),
      axis.text.x = element_text(angle = x_axis_angle, vjust = .5)
    ),
    guides(fill = guide_legend(nrow = 1))
  )
}

# Fig. 3 ------------------------------------------------------------------

## Males F2 temp comparison with negative binomial GLM due to zero counts
# Fitting separate models for each cell_type
labels_hist_m_temp_f2 <- 
  summary_hist_m_long %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm.nb(count ~ temp_f2, data = .x)),
    results = map(model, tidy)
  ) %>%
  unnest(results) %>%
  filter(term == "temp_f2T") %>%
  mutate(
    label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    group1 = "C",
    group2 = "T",
    p.adj = p.value
  ) %>%
  dplyr::select(cell_type, group1, group2, p.adj, label) %>%
  ungroup() %>% 
  mutate(y.position = c(250, 75, 250, 150, 15)) %>% 
  filter(cell_type != "Spt")

# Manual labels for super and subscript
costum_labels_m <- list("Aund" = expression("A"["und"]), 
                        "Adiff" = expression("A"["diff"]), 
                        "B" = "B",
                        "Spc" = "Spc",
                        "Spt" = "Spt")

# Plot
p_hist_m <-
  ggplot(summary_hist_m_long, aes(x = temp_f2, y = count, fill = temp_f2)) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.25) +
  stat_summary(fun = mean, geom = "col") +
  facet_wrap(~ cell_type,
             scales = "free_y",
             labeller = labeller(cell_type = as_labeller(costum_labels_m, label_parsed))) +
  ggstyle(color = color_temp,
          fill = "Temperature F2",
          x_axis_angle = 0) +
  stat_pvalue_manual(
    data = labels_hist_m_temp_f2,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    inherit.aes = FALSE
  )

ggsave("02_plots/04_hist_males.pdf", p_hist_m, width = 6, height = 5)

## Males F2 temp comparison with negative binomial GLM
# Fitting separate models for each cell_type
labels_hist_f_temp_f2 <- 
  summary_hist_f_long %>%
  group_by(cell_type) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm.nb(count ~ temp_f2, data = .x)),
    results = map(model, tidy),
    max_y = map_dbl(data, ~ max(.x$count, na.rm = TRUE))
  ) %>%
  unnest(results) %>%
  filter(term == "temp_f2T") %>%
  mutate(
    signif_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    group1 = "C",
    group2 = "T",
    p.adj = p.value,
    y.position = max_y * 1.05
  ) %>%
  dplyr::select(cell_type, group1, group2, p.adj, label = signif_label, y.position) %>%
  ungroup()

# Plot
p_hist_f <-
  ggplot(summary_hist_f_long, aes(x = temp_f2, y = count, fill = temp_f2)) +
  stat_summary(fun.data = mean_sd_upper,
               geom = "errorbar",
               width = 0.25) +
  stat_summary(fun = mean, geom = "col") +
  facet_wrap( ~ cell_type, scales = "free_y") +
  ggstyle(color = color_temp,
          fill = "Temperature F2",
          x_axis_angle = 0) +
  stat_pvalue_manual(
    data = labels_hist_f_temp_f2,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    inherit.aes = FALSE
  )

ggsave("02_plots/04_hist_females.pdf", p_hist_f, width = 6, height = 5)

# Final figure
# Counts per temp f2 
p_hist_m + plot_spacer() + p_hist_f + 
  plot_layout(widths = c(3,.5,2), guides = "collect") &
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18), 
        legend.position = "none", legend.direction = "vertical", 
        strip.text = element_text(face = "bold"))  

ggsave("02_plots/fig_3.pdf", width = 8, height = 3.5)
ggsave("02_plots/fig_3.png", width = 8, height = 3.5)
