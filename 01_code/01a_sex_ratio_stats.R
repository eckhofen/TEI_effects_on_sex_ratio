# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: October 2025
# Description: Sex ratio analysis  

# Settings ----------------------------------------------------------------

library(tidyverse)
library(gt)
library(rstatix)
library(flextable)

# Loading data_main ------------------------------------------------------------

# Main data
load("00_data/data_main.RData")

# Loading color palettes
load("00_data/color_palettes.RData")

# Summary -----------------------------------------------------------------

# Simple overview
sex_ratio_summary <- data_main %>% 
  group_by(treatment, temp_f2) %>% 
  summarise(n = n(), n_female = sum(sex == "F"), ratio_female = sum(sex == "F")/n,
            n_male = sum(sex == "M"), ratio_male = sum(sex == "M")/n, .groups = "drop") %>% 
  arrange(treatment)
  
write.csv(sex_ratio_summary, "00_data/01_results/01_sex_ratio_summary.csv", row.names = FALSE)

# Sex ratio more detailed
sex_ratio_summary_det <- data_main %>% 
  group_by(line_f0, treatment, temp_f0, temp_f1, temp_f2) %>%
  summarise(n = n(), n_female = sum(sex == "F"), ratio_female = sum(sex == "F")/n, 
            n_male = sum(sex == "M"), ratio_male = sum(sex == "M")/n, .groups = "drop") %>%
  arrange(treatment)

write.csv(sex_ratio_summary_det, "00_data/01_results/01_sex_ratio_summary_det.csv", row.names = FALSE)

## Adding estimation value (missing tank S2-CTT)
# See explanation: "08_missing_tank_value_estimation.R"

# Missing tank
tank_missing <- sex_ratio_summary_det[6,] %>% 
  mutate(treatment = "G4", 
         temp_f2 = "T", 
         n_female = 54, 
         n_male = 46, 
         ratio_female = 0.54, 
         ratio_male = 0.46)

sex_ratio_summary_det_w_est <- rbind(sex_ratio_summary_det,tank_missing)

# Getting sex ratio for all groups 
sex_ratio <- data_main %>% 
  group_by(temp_f0, temp_f1, temp_f2, line_f0) %>% 
  summarise(male_ratio = 1-sum(sex_binary)/length(sex_binary))

# Checking average sex ratio change for F2 T
sex_ratio_diff_F2 <- data_main %>% 
  group_by(temp_f2) %>% mutate(male_ratio = 1-sum(sex_binary)/length(sex_binary), n = n()) %>% 
  distinct(temp_f2, male_ratio, n) %>% 
  arrange(temp_f2)

# Taking all C vs T directly (treating them as two groups only) 
diff(sex_ratio_diff_F2$male_ratio)
# 0.2192214

sex_ratio %>% 
  group_by(temp_f2) %>% 
  summarise(male_ratio_mean = mean(male_ratio)) %>% 
  mutate(male_ratio_diff = male_ratio_mean - lag(male_ratio_mean))
# 0.218

# Checking average male proportion for each replicate (F0 lines)
sex_ratio_diff_F0 <- 
  data_main %>% 
  group_by(temp_f0, line_f0) %>% 
  mutate(male_ratio = 1-sum(sex_binary)/length(sex_binary), n = n()) %>% 
  distinct(temp_f0, line_f0, male_ratio) %>% 
  arrange(temp_f0)

# Comparison of groups with Fisher's exact test ---------------------------

## F2
# Creating contingency table for the f2 
cont_table_temp_f2 <-
  sex_ratio_summary %>%
  group_by(temp_f2) %>%
  summarise(
    total_females = sum(n_female),
    total_males = sum(n_male),
    .groups = 'drop'
  ) %>%
  select(total_females, total_males) %>% # Select the count columns
  as.matrix()

rownames(cont_table_temp_f2) <- c("C", "T")

test_fisher_TE_f2 <- 
  cont_table_temp_f2 %>% 
  fisher.test()

## F1
# Creating contingency table for the f1 
cont_table_temp_f1 <- 
  sex_ratio_summary_det %>%
  group_by(temp_f1) %>%
  summarise(total_females = sum(n_female),
            total_males = sum(n_male),
            .groups = 'drop')  %>% 
  select(total_females, total_males) %>% # Select the count columns
  as.matrix()

rownames(cont_table_temp_f1) <- c("C", "T")

test_fisher_TE_f1 <- 
  cont_table_temp_f1 %>% 
  fisher.test()

## F0
# Creating contingency table for the f0 
cont_table_temp_f0 <- 
  sex_ratio_summary_det %>%
  group_by(temp_f0) %>%
  summarise(total_females = sum(n_female),
            total_males = sum(n_male),
            .groups = 'drop')  %>% 
  select(total_females, total_males) %>% # Select the count columns
  as.matrix()

rownames(cont_table_temp_f0) <- c("C", "T")

test_fisher_TE_f0 <-
  cont_table_temp_f0 %>% 
  fisher.test()

test_fisher_TE <- 
  list(test_fisher_TE_f2, test_fisher_TE_f1, test_fisher_TE_f0)

## Specific comparisons (X = either C or T; TE = temperature effect)
# Parental effect (temperature for F1)

# TE F1 for all F0 lines
f0_lines <- c("S1", "S2", "S3", "S4")

for(f0_line in f0_lines) {
  print(f0_line)
  x <- 
    sex_ratio_summary_det %>%
    filter(line_f0 == f0_line) %>% 
    group_by(temp_f1) %>%
    summarise(total_females = sum(n_female),
              total_males = sum(n_male),
              .groups = 'drop')  %>% 
    print()
  x %>% 
    select(total_females, total_males) %>%
    as.matrix() %>% 
    fisher.test() %>% 
    print()
} 

# TE F1 S3 (CC + CT vs TC + TT)
sex_ratio_summary_det %>%
  filter(line_f0 == "S3") %>% 
  group_by(temp_f1) %>%
  summarise(total_females = sum(n_female),
            total_males = sum(n_male),
            .groups = 'drop')  %>% 
  select(total_females, total_males) %>%
  as.matrix() %>% 
  fisher.test()

# TE F1 S2 (CC vs TC)
sex_ratio_summary_det %>%
  filter(line_f0 == "S2", temp_f2 == "C") %>% 
  group_by(temp_f1) %>%
  summarise(total_females = sum(n_female),
            total_males = sum(n_male),
            .groups = 'drop')  %>% 
  select(total_females, total_males) %>%
  as.matrix() %>% 
  fisher.test()

# XCC vs XCT
sex_ratio_summary_det %>%
  filter(temp_f1 == "C") %>% 
  group_by(temp_f2) %>%
  summarise(total_females = sum(n_female),
            total_males = sum(n_male),
            .groups = 'drop')  %>% 
  select(total_females, total_males) %>%
  as.matrix() %>% 
  fisher.test()

## Comparing SR 
sex_ratio_summary_det %>% 
  group_by(temp_f2) %>% 
  summarise(ratio_male = mean(ratio_male)) %>% 
  ungroup() %>% 
  mutate(diff_sr = diff(ratio_male))

# Inferential statistics --------------------------------------------------

# Probit regression
data_probit <- 
  data_main %>% 
  select(line_f0, sex, sex_binary, temp_f2, temp_f1, temp_f0)

# Finding best fit based on AIC
glm(sex_binary ~ temp_f2 * temp_f1 + temp_f0, 
    data = data_probit, 
    family = binomial(link = "probit"))$aic

glm(sex_binary ~ temp_f2 * temp_f1 * temp_f0, 
    data = data_probit, 
    family = binomial(link = "probit"))$aic

glm(sex_binary ~ temp_f2 * temp_f1 + line_f0, 
    data = data_probit, 
    family = binomial(link = "probit"))$aic

glm(sex_binary ~ temp_f2 * temp_f1 * line_f0, 
    data = data_probit, 
    family = binomial(link = "probit"))$aic
# sex_binary ~ temp_f2 * temp_f1 * line_f0 has lowest AIC (1823)

# Fitting the model
probit_model <- glm(sex_binary ~ temp_f2 * temp_f1 * line_f0, 
                    data = data_probit, 
                    family = binomial(link = "probit"))

summary(probit_model)

probit_model_tidy <- 
  probit_model %>% 
  tidy()

# Saving results
write_csv(x = probit_model_tidy, 
          file = "00_data/01_results/02_sex_ratio_probit_model.csv")

# Clean table for word
t_probit_model <- probit_model_tidy %>%
  select(term, estimate, std.error, statistic, p.value) %>%
  flextable() %>%
  set_header_labels(
    term = "Predictor",
    estimate = "Estimate",
    std.error = "Std. Error",
    statistic = "z-value",
    p.value = "p-value"
  ) %>%
  colformat_double(
    j = c("estimate", "std.error", "statistic"),
    digits = 3 
  ) %>%
  colformat_double(
    j = "p.value",
    digits = 3,
    na_str = "NA*"
  ) %>%
  set_caption(caption = "Table 1: Results of Probit Regression on Sex Proportion.") %>%
  theme_booktabs() %>%
  add_footer_lines(
    values = paste0("Samples n = ", probit_model$df.null, 
                    "; Degrees of freedom = ", probit_model$df.residual, 
                    "; AIC = ", round(probit_model$aic, 2), 
                    ". *S2 line has a missing group (CTT)."
  )) %>%
  autofit() %>%
  align(j = c("estimate", "std.error", "statistic", "p.value"), align = "center", part = "all")

t_probit_model

# To save the table to a Word document:
save_as_docx(t_probit_model, path = "03_tables/01_probit_model.docx")


# Probit regression with estimate -----------------------------------------

# Creating estimated data_main 
data_est <- tibble(
  line_f0 = "S2",
  sex = c(rep("M", tank_missing$n_male), 
          rep("F", tank_missing$n_female)),
  sex_binary = ifelse(sex == "M", 0, 1),
  temp_f2 = "T",
  temp_f1 = "T",
  temp_f0 = "C"
)
# Adding estimate data
data_probit_w_est <- 
  data_probit %>% 
  bind_rows(data_est)

## Verifying
# Creating ratios from estimated data
sr_data_probit_w_est <- 
  data_probit_w_est %>% 
  group_by(line_f0, temp_f2, temp_f1, temp_f0) %>% 
  summarise(ratio_male = sum(sex == "M")/n()) %>% 
  ungroup() %>% 
  arrange(line_f0, temp_f2, temp_f1, temp_f0)

# Creating ratios from original estimated data
sr_data_orig <- 
  sex_ratio_summary_det_w_est %>% 
  arrange(line_f0, temp_f2, temp_f1, temp_f0)

# Check if match
sr_data_probit_w_est$ratio_male == sr_data_orig$ratio_male

## Fitting the model
# Checking fit 
glm(sex_binary ~ temp_f2 * temp_f1 + temp_f0, 
    data = data_probit_w_est, 
    family = binomial(link = "probit"))$aic

glm(sex_binary ~ temp_f2 * temp_f1 * temp_f0, 
    data = data_probit_w_est, 
    family = binomial(link = "probit"))$aic

glm(sex_binary ~ temp_f2 * temp_f1 + line_f0, 
    data = data_probit_w_est, 
    family = binomial(link = "probit"))$aic

glm(sex_binary ~ temp_f2 * temp_f1 * line_f0, 
    data = data_probit_w_est, 
    family = binomial(link = "probit"))$aic

# Model with interaction term still ahs lowest AIC (1963)

# Fitting actual model
probit_model_w_est <- glm(sex_binary ~ temp_f2 * temp_f1 * line_f0, 
                          data = data_probit_w_est, 
                          family = binomial(link = "probit"))

summary(probit_model_w_est)

probit_model_w_est_tidy <- 
  probit_model_w_est %>% 
  tidy()

# Saving results
write_csv(x = probit_model_w_est_tidy, 
          file = "00_data/01_results/01_sex_ratio_probit_model_w_est.csv")

# Clean table for word
t_probit_model_w_est <- probit_model_w_est_tidy %>%
  select(term, estimate, std.error, statistic, p.value) %>%
  flextable() %>%
  set_header_labels(
    term = "Predictor",
    estimate = "Estimate",
    std.error = "Std. Error",
    statistic = "z-value",
    p.value = "p-value"
  ) %>%
  colformat_double(
    j = c("estimate", "std.error", "statistic"),
    digits = 3 
  ) %>%
  colformat_double(
    j = "p.value",
    digits = 3,
    na_str = "NA"
  ) %>%
  set_caption(caption = "Table 1: Results of Probit Regression on Sex Proportion with estimate.") %>%
  theme_booktabs() %>%
  add_footer_lines(
    values = paste0("Samples n = ", probit_model_w_est$df.null, 
                    "; Degrees of freedom = ", probit_model_w_est$df.residual, 
                    "; AIC = ", round(probit_model_w_est$aic, 2), 
                    "."
    )) %>%
  autofit() %>%
  align(j = c("estimate", "std.error", "statistic", "p.value"), align = "center", part = "all")

t_probit_model_w_est

# To save the table to a Word document:
save_as_docx(t_probit_model_w_est, path = "03_tables/01_probit_model_w_est.docx")