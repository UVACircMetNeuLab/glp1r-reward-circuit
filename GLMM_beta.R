# GLMM Beta Regression
# Isabelle Sajonia

# Load libraries
library(dplyr)
library(tidyr)
library(glmmTMB)
library(emmeans)
library(ggplot2)

setwd("C:/Users/Owner/Desktop/ML")

# -------------------------------
# 1. Preprocess behavioral data (time spent performing grouped behaviors)

# Load data
behavior_data <- read.csv('behavior_time_v3.csv')

# Clean names for columns 3 through 8
colnames(behavior_data)[3:8] <- c("total_drink", "total_food", "total_groom", "total_move", "total_move_shelter", "total_rest")

# Remove IDs excluded 
behavior_data <- behavior_data %>%
  filter(!id %in% c(24, 11, 93))

# Normalize numeric columns (excluding ID) by total time to get proportion (0-1)
behavior_data <- behavior_data %>%
  rowwise() %>%
  dplyr::mutate(
    total_time = sum(c_across(c(total_drink, total_food, total_groom, total_move, total_move_shelter, total_rest)), na.rm = TRUE)
  ) %>%
  dplyr::mutate(across(
    c(total_drink, total_food, total_groom, total_move, total_move_shelter, total_rest),
    ~ .x / total_time
  )) %>%
  ungroup()

# Combine move out of shelter and move in shelter columns and remove 
# total_move_shelter
behavior_data <- behavior_data %>%
  dplyr::mutate(total_move = total_move + total_move_shelter) %>%
  dplyr::select(-total_move_shelter)

# Filter for groups of interest
behavior_data <- behavior_data %>%
  filter(group %in% c("LiCl", "saline_licl","dan_sd", "veh_sd","lira_sd", 
                      "saline_lira","orfo", "saline_orfo" ,"fed"))

# Check n per group
behavior_data %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n_unique_ids = n_distinct(id)) %>%
  print()

# Group and ID must be factors
behavior_data$group <- factor(behavior_data$group)
behavior_data$id <- factor(behavior_data$id)

# -------------------------------
# 2. Process behavioral proportion data 

data_long <- behavior_data %>%
  pivot_longer(cols = 3:7, names_to = "measurement", values_to = "value") %>%
  filter(value >= 0 & value <= 1) %>%
  mutate(value = ifelse(value == 0, 1e-6, value))  # Avoid zeros for beta reg.

# -------------------------------
# 3. Density plots used to guide link function selection based on distribution shape
ggplot(data_long, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~measurement, scales = "free_y") +
  xlim(0, 1) +
  theme_minimal() +
  labs(title = "Density of behavior values by group", x = "Proportion", y = "Density")


# -------------------------------
# 4. Define best link functions for each measurement

link_map <- list(
  "total_rest"         = "logit",
  "total_move"         = "logit",
  "total_food"         = "cloglog",
  "total_groom"        = "cloglog",
  "total_drink"        = "cloglog"
  #,"total_move_shelter" = "cloglog"
)

# -------------------------------
# 5. GLMM for saline vs LiCl (this is an example, you can run for different
# group comparisons)

licl_data <- data_long %>%
  filter(group %in% c("LiCl", "saline_licl"))

measurements <- unique(licl_data$measurement)

for (m in measurements) {
  cat("\nBeta Regression for:", m, "\n")
  
  filtered <- licl_data %>%
    filter(measurement == m & value > 0 & value < 1) %>%
    droplevels()
  
  link_fn <- link_map[[m]]
  
  model <- glmmTMB(
    value ~ group + (1 | id), # Adds a random intercept for mouse (paired)
    family = beta_family(link = link_fn),
    data = filtered
  )
  
  print(summary(model))
  print(summary(model)$coefficients$cond[, "Pr(>|z|)"])
  
}

# -------------------------------
# 6. GLMM across all treatment groups (optional)

multi_group_data <- data_long %>%
  filter(group %in% c("fed", "dan_sd", "LiCl", "lira_sd", "orfo")) %>%
  mutate(group = factor(group, levels = c("fed", "dan_sd", "LiCl", "lira_sd", "orfo")))

measurements <- unique(multi_group_data$measurement)

for (m in measurements) {
  cat("\nBeta Regression for:", m, "\n")
  
  filtered <- multi_group_data %>%
    filter(measurement == m & value > 0 & value < 1) %>%
    droplevels()
  
  if (nrow(filtered) < 5) {
    cat("Not enough data for", m, "\n")
    next
  }
  
  link_fn <- link_map[[m]]
  
  model <- glmmTMB(
    value ~ group,
    family = beta_family(link = link_fn),
    data = filtered
  )
  
  print(summary(model))
  print(summary(model)$coefficients$cond[, "Pr(>|z|)"])
  
  # Compare to "fed"
  emm <- emmeans(model, ~ group)
  fed_contrasts <- contrast(emm, method = "trt.vs.ctrl", ref = "fed", 
                            adjust = "holm")
  print(fed_contrasts)
}
