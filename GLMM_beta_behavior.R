# GLMM Beta Regression - Behavioral Time Data

# Fits beta regression models (via glmmTMB) to proportion of time spent in
# each behavior. Beta regression for bounded (0,1) proportional data. 
# Link functions were selected per-behavior based on distribution shape.

# Two comparison types:
#   - Paired (vehicle vs. drug): random intercept for mouse ID
#   - Multi-group (all drug groups vs. fed): fixed effects only, emmeans
#     contrasts vs. fed as reference

library(dplyr)
library(tidyr)
library(glmmTMB)
library(emmeans)
library(ggplot2)

setwd("C:/Users/Owner/Desktop/ML")

# 1. Load and process behavior time data
behavior_data <- read.csv("behavior_time_dark_cycle.csv")

colnames(behavior_data)[3:8] <- c(
  "total_drink", "total_food", "total_groom",
  "total_move", "total_move_shelter", "total_rest"
)

behavior_data <- behavior_data %>%
  filter(!id %in% c(24, 11, 93)) %>%
  # Normalize to proportion of total observed time (keep as 0-1 for beta reg.)
  rowwise() %>%
  mutate(total_time = sum(c_across(c(
    total_drink, total_food, total_groom,
    total_move, total_move_shelter, total_rest
  )), na.rm = TRUE)) %>%
  mutate(across(
    c(total_drink, total_food, total_groom, total_move, total_move_shelter, total_rest),
    ~ .x / total_time
  )) %>%
  ungroup() %>%
  # Combine cage floor and shelter movement
  dplyr::mutate(total_move = total_move + total_move_shelter) %>%
  dplyr::select(-total_move_shelter, -total_time) %>%
  filter(group %in% c(
    "LiCl", "saline_licl",
    "dan_sd", "veh_sd",
    "lira_sd", "saline_lira",
    "orfo", "saline_orfo",
    "fed"
  )) %>%
  mutate(
    group = factor(group),
    id    = factor(id)
  )

# Sanity check
behavior_data %>%
  group_by(group) %>%
  summarise(n_mice = n_distinct(id)) %>%
  print()

sum(behavior_data[, c("total_drink", "total_food", "total_groom", "total_move", "total_rest")] == 0, na.rm = TRUE)

# Reshape + add a small epsilon to replace exact zeros (though none found here)
data_long <- behavior_data %>%
  pivot_longer(cols = c(total_drink, total_food, total_groom, total_move, total_rest),
               names_to = "measurement", values_to = "value") %>%
  filter(value >= 0 & value <= 1) %>%
  mutate(value = ifelse(value == 0, 1e-6, value))  

# 2. Distribution plot (guide link function selection)
ggplot(data_long, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~measurement, scales = "free_y") +
  xlim(0, 1) +
  theme_minimal() +
  labs(title = "Behavior proportion distributions by group",
       x = "Proportion of time", y = "Density")

# 3. Link function map

# Logit: symmetric distributions centered away from bounds
# Cloglog: right-skewed distributions with mass near 0 (rare behaviors)

link_map <- list(
  total_rest  = "logit",
  total_move  = "logit",
  total_food  = "cloglog",
  total_groom = "cloglog",
  total_drink = "cloglog"
)


# 4. PAIRED COMPARISONS - vehicle vs. drug (random intercept for mouse ID)

# Edit the group filter below to run for different vehicle/drug pairs.
# Current example: saline_licl vs. LiCl

paired_groups <- c("LiCl", "saline_licl")  # <- change this for other pairs

paired_data <- data_long %>%
  filter(group %in% paired_groups) %>%
  droplevels()

for (m in names(link_map)) {
  cat("\n--- Paired GLMM:", m, "(", paste(paired_groups, collapse = " vs. "), ") ---\n")
  
  filtered <- paired_data %>%
    filter(measurement == m, value > 0, value < 1) %>%
    droplevels()
  
  model <- glmmTMB(
    value ~ group + (1 | id),  # random intercept = paired design
    family = beta_family(link = link_map[[m]]),
    data = filtered
  )
  
  print(summary(model)$coefficients$cond)
}

# 5. Multi-group comparison - all drug groups vs. fed (reference)
multi_data <- data_long %>%
  filter(group %in% c("fed", "dan_sd", "LiCl", "lira_sd", "orfo")) %>%
  mutate(group = factor(group, levels = c("fed", "dan_sd", "LiCl", "lira_sd", "orfo")))

for (m in names(link_map)) {
  cat("\n--- Multi-group GLMM:", m, "---\n")
  
  filtered <- multi_data %>%
    filter(measurement == m, value > 0, value < 1) %>%
    droplevels()
  
  if (nrow(filtered) < 5) {
    cat("Skipping", m, "- insufficient data\n")
    next
  }
  
  model <- glmmTMB(
    value ~ group,
    family = beta_family(link = link_map[[m]]),
    data = filtered
  )
  
  print(summary(model)$coefficients$cond)
  
  # Pairwise contrasts vs. fed (Holm correction for multiple comparisons)
  emm <- emmeans(model, ~ group)
  print(contrast(emm, method = "trt.vs.ctrl", ref = "fed", adjust = "holm"))
}
