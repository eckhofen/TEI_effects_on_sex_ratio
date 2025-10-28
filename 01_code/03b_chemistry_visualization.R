# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: October 2025
# Description: exploratory analysis of TEI chemical composition data 

# Settings ----------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(patchwork)
library(compositions)
library(outliers)

# Loading data ------------------------------------------------------------

# Loading data
load("00_data/00_processed/data_chemistry.RData")

# Loading color palettes
load("00_data/color_palettes.RData")

# Fig. 4a, c, e (PCAs) ---------------------------------------------------------------------

# Creating function to transform, process and plot the data into a PCA
PCA_comp_fun <- function(data, comp_vars, pc1 = 1, pc2 = 2, meta_data = NULL, plot = TRUE, 
                         title = "Compositional PCA", color = "sex", color_legend = "F2 temperature", 
                         color_pal = color_sex, zero_replacement = 0.01) {
  # Select only the compositional variables
  comp_data <- data %>% select(all_of(comp_vars))
  
  # Handling zeros
  if (any(comp_data == 0)) {
    comp_data[comp_data == 0] <- zero_replacement
  }
  
  # Apply the Centered Log-Ratio transformation
  clr_data <- clr(acomp(comp_data))

  pc <- prcomp(clr_data, center = TRUE, scale. = TRUE)
  
  explained_variance <- pc$sdev^2 / sum(pc$sdev^2) * 100
  
  # Passing as label
  x_label <- paste0("PC", pc1, " (", round(explained_variance[pc1], 2), "%)")
  y_label <- paste0("PC", pc2, " (", round(explained_variance[pc2], 2), "%)")
  
  # Adding meta data
  pc_scores <- cbind(meta_data, pc$x)
  
  pca_results <- list()
  
  if(plot == TRUE) {
    # Defining PCs for plotting
    x <- paste0("PC", pc1)
    y <- paste0("PC", pc2)
    
    # Plotting
    pca_plot <- ggplot(pc_scores, aes(x = .data[[x]], y = .data[[y]], fill = .data[[color]])) +
      geom_point(shape = 21, size = 3, stroke = .5, color = "white", alpha = 0.8) +
      labs(
        x = x_label,
        y = y_label, 
        fill = color_legend,
        title = title
      ) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
      scale_fill_manual(values = color_pal) +
      scale_x_continuous(expand = expansion(c(.1,.1))) +
      scale_y_continuous(expand = expansion(c(.1,.1))) +
      theme_pubr(border = TRUE) +
      theme(plot.title = element_text(face = "bold", hjust = .5), 
            legend.position = "bottom", 
            panel.background = element_rect(linewidth = 1.5, color = "black")) 
     
    pca_results <- list(pca_object = pc, data = pc_scores, plot = pca_plot)
    
    return(pca_results) 
    
  } else {
    pca_results <- list(pca_object = pc, data = pc_scores)
    return(pca_results)
  }
}

## PCA for all
# Mass
cc_pca_mass <-
  PCA_comp_fun(
    data = data_cc_mass,
    comp_vars = data_cc_mass %>% select(where(is.numeric)) %>% colnames(),
    meta_data = meta_data_cc,
    title = "PCA Mass",
    color_legend = "Sex"
  )

ggsave(filename = "02_plots/03_cc_PCA_mass.pdf",
       plot = cc_pca_mass$plot,
       width = 6,
       height = 6)

# Males
data_cc_mass_m <- data_cc_mass %>% filter(sex == "M")
meta_data_cc_m <- meta_data_cc %>% filter(sex == "M")

# Mass
cc_pca_mass_m <-
  PCA_comp_fun(
    data = data_cc_mass_m,
    comp_vars = data_cc_mass_m %>% select(where(is.numeric)) %>% colnames(),
    meta_data = meta_data_cc_m,
    title = "PCA mass (males)",
    color_legend = "Temperature", 
    color = "treatment", 
    color_pal = color_temp
  )

ggsave(filename = "02_plots/03_cc_PCA_mass_m.pdf",
       plot = cc_pca_mass_m$plot,
       width = 6,
       height = 6)

# Females

data_cc_mass_f <- data_cc_mass %>% filter(sex == "F")
meta_data_cc_f <- meta_data_cc %>% filter(sex == "F")

# Mass
cc_pca_mass_f <-
  PCA_comp_fun(
    data = data_cc_mass_f,
    comp_vars = data_cc_mass_f %>% select(where(is.numeric)) %>% colnames(),
    meta_data = meta_data_cc_f,
    title = "PCA mass (females)",
    color_legend = "Temperature", 
    color = "treatment", 
    color_pal = color_temp
  )

ggsave(filename = "02_plots/03_cc_PCA_mass_f.pdf",
       plot = cc_pca_mass_f$plot,
       width = 6,
       height = 6)

# Fig. 4b, d, f (boxplots) ---------------------------------------------------------------------

# Pivot longer
data_cc_mass_long <- 
  pivot_longer(
    data_cc_mass, 
    cols = where(is.numeric),
    names_to = "element", 
    values_to = "percent"
  ) %>% 
  mutate(
    element = str_to_title(str_replace(element, "_mass_per", "")),
    # Handling zero-reads
    percent = ifelse(percent == 0, 0.01, percent)
  )

# Testing normality
data_cc_mass_long %>% 
  group_by(element) %>% 
  shapiro_test(percent)

# Normality not given with exceptions
# Proceeding with Wilcoxon test

## Sex
cc_bp_sex <- ggplot(data_cc_mass_long, aes(x = sex, y = percent, fill = sex)) +
  geom_boxplot(show.legend = FALSE) +
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif",
    comparisons = list(c("F", "M"))) +
  labs(x = "Sex", 
       y = "Element mass (%)", 
       fill = "Sex", 
       title = "Element mass") +
  facet_wrap(~ element, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = color_sex) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  theme_classic2() +
  theme(plot.title = element_text(face = "bold", hjust = .5),
        legend.position = "bottom", legend.direction = "horizontal",
        strip.background = element_rect(linewidth = NA, fill = "grey95"), 
        strip.text = element_text(face = "bold"))

ggsave(filename = "02_plots/03_cc_bp_sex.pdf", plot = cc_bp_sex, width = 9, height = 5)


## Males
cc_bp_mass_m <- ggplot(data_cc_mass_long %>% filter(sex == "M"), aes(x = treatment, y = percent, fill = treatment)) +
  geom_boxplot(show.legend = FALSE) +
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif",
    comparisons = list(c("C", "T"))) +
  labs(x = "Temperature", 
       y = "Element mass (%)", 
       fill = "Temperature", 
       title = "Element mass (Males)") +
  facet_wrap(~ element, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = color_temp) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  theme_classic2() +
  theme(plot.title = element_text(face = "bold", hjust = .5),         legend.position = "bottom", legend.direction = "horizontal",         strip.background = element_rect(linewidth = NA, fill = "grey95"),          strip.text = element_text(face = "bold"))

ggsave(filename = "02_plots/03_cc_bp_mass_m.pdf", plot = cc_bp_mass_m, width = 9, height = 5)

## Females
cc_bp_mass_f <- ggplot(data_cc_mass_long %>% filter(sex == "F"), aes(x = treatment, y = percent, fill = treatment)) +
  geom_boxplot(show.legend = FALSE) +
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif",
    comparisons = list(c("C", "T"))) +
  labs(x = "Temperature", 
       y = "Element mass (%)", 
       fill = "Temperature", 
       title = "Element mass (Females)") +
  facet_wrap(~ element, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = color_temp) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  theme_classic2() +
  theme(plot.title = element_text(face = "bold", hjust = .5),legend.position = "bottom", legend.direction = "horizontal",         strip.background = element_rect(linewidth = NA, fill = "grey95"),          strip.text = element_text(face = "bold"))

ggsave(filename = "02_plots/03_cc_bp_mass_f.pdf", plot = cc_bp_mass_f, width = 9, height = 5)

## Final plot
free(cc_pca_mass$plot) + cc_bp_sex + free(cc_pca_mass_m$plot) + cc_bp_mass_m + free(cc_pca_mass_f$plot) + cc_bp_mass_f +
  plot_layout(nrow = 3, widths = c(1,3), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 28), plot.title = element_blank(), 
        legend.position = "right", legend.direction = "vertical", strip.text = element_text(size = 12), 
        axis.title = element_text(size = 14))

ggsave("02_plots/fig_4.pdf", width = 13, height = 9)
ggsave("02_plots/fig_4.png", width = 13, height = 9)
  