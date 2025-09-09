# PCA and MANOVA Analysis Pipeline
# Isabelle Sajonia
# Description: Behavioral data analysis using PCA, MANOVA, and permutation testing

# Load libraries
library(MANOVA.RM)
library(tidyverse)
library(plotly)
library(scales)

# Set working directory
setwd("C:/Users/Owner/Desktop/ML")

# -------------------------------
# 1. Load and combine datasets

# Time spent performing grouped behavior 
behavior_data <- read.csv("behavior_time_v3.csv") %>%
  filter(!id %in% c(24, 11, 93),
         group %in% c("dan_sd", "lira_sd", "orfo", "LiCl", "fed")) %>%
  mutate(total_move = move.explore.sec + move.explore...shelter.sec) %>%
  select(-move.explore.sec, -move.explore...shelter.sec)

# Grouped behavior (normalized) transitions 
transitions <- read_csv("normalized_transitions_v4.csv") %>%
  filter(!id %in% c(24, 11, 93),
         group %in% c("dan_sd", "lira_sd", "orfo", "LiCl", "fed"))

# Combine reciprocal transitions
reciprocal_pairs <- list(
  c("groom_drink", "drink_groom"),
  c("food motivated_move/explore", "move/explore_food motivated"),
  c("groom_food motivated", "food motivated_groom"),
  c("groom_move/explore", "move/explore_groom"),
  c("move/explore_drink", "drink_move/explore"),
  c("move/explore_shelter", "shelter_move/explore")
)

for (pair in reciprocal_pairs) {
  transitions[[paste0(pair[1], "_combined")]] <- rowSums(transitions[, pair], na.rm = TRUE)
}

transitions_combined <- transitions %>%
  select(id, group, ends_with("_combined")) %>%
  select(-matches("groom_drink|drink_groom|groom_food motivated|food motivated_groom"))

# Grouped behavior bouts
bouts <- read_csv("bouts_v3.csv") %>%
  filter(!id %in% c(24, 11, 93),
         group %in% c("dan_sd", "lira_sd", "orfo", "LiCl", "fed")) %>%
  select(-starts_with("max")) # variable, just use average and # bouts

# -------------------------------
# 2. Merge and process dataset

# Make sure id and group columns are factors in all datasets
colnames(behavior_data)[1] <- "id"
behavior_data <- behavior_data %>% mutate(across(c(id, group), as.factor))
transitions_combined <- transitions_combined %>% mutate(across(c(id, group), as.factor))
bouts <- bouts %>% mutate(across(c(id, group), as.factor))

# Merge 
combined_data <- behavior_data %>%
  inner_join(transitions_combined, by = c("id", "group")) %>%
  inner_join(bouts, by = c("id", "group"))

# Clean feature names
names(combined_data) <- make.names(names(combined_data))

# Remove food motivated behavior features
feature_cols <- setdiff(names(combined_data), c("id", "group"))
filtered_cols <- feature_cols[!grepl("food\\.motivated", feature_cols)]

# Final dataset
selected_features <- combined_data[, c("id", filtered_cols, "group")]

# -------------------------------
# 3. PCA

# Scale features
combined_scaled <- selected_features %>%
  mutate(across(all_of(filtered_cols), scale))

# Run PCA
pca_res <- prcomp(combined_scaled[, filtered_cols], center = FALSE, scale. = FALSE)

# Variance explained
explained_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
cumulative_var <- cumsum(explained_var)
n_pc_90 <- which(cumulative_var >= 0.90)[1]

# Scree plot
pca_df <- data.frame(PC = 1:length(explained_var), CumulativeVariance = cumulative_var)

ggplot(pca_df, aes(x = PC, y = CumulativeVariance)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(color = "blue", size = 2) +
  geom_vline(xintercept = n_pc_90, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
  annotate("text", x = n_pc_90 + 0.5, y = 0.85, label = paste(n_pc_90, "PCs"), color = "red") +
  scale_y_continuous(labels = percent_format()) +
  labs(title = "Cumulative Variance Explained by PCA", x = "Principal Component", y = "Cumulative Variance") +
  theme_minimal()

# -------------------------------
# 4. MANOVA (parametric bootstrap)
pc_scores <- as.data.frame(pca_res$x[, 1:3]) %>% # PC 1-3
  mutate(group = combined_scaled$group,
         id = combined_scaled$id)

# Group order
pc_scores$group <- factor(pc_scores$group, levels = c("orfo", "fed","dan_sd", 
                                                      "lira_sd", "LiCl"))

formula_pc <- as.formula(paste("cbind(", paste(names(pc_scores)[1:3], 
                                               collapse = ", "), ") ~ group"))

set.seed(120) # for reproducibility 
fit_pc <- MANOVA.wide(
  formula = formula_pc,
  data = pc_scores,
  resampling = "paramBS",
  iter = 5000, # 5000 iterations
  alpha = 0.05
)

summary(fit_pc)

# Pairwise comparisons (orfo vs all other groups to minimize comparisons given
# moderate amount of groups and low sample size)
pairwise_pc <- simCI(fit_pc, contrast = "pairwise", type = "Dunnett", 
                     ref = "orfo")

# -------------------------------
# 5. PCA Visualization

# Recode group names (remove '_sd')
pc_scores <- pc_scores %>%
  mutate(group = recode(group, "dan_sd" = "dan", "lira_sd" = "lira"))

# Assign colors for PCA
group_colors <- c("dan" = "#ED1C24", "lira" = "#26A7DF", "LiCl" = "#6ebe45", 
                  "fed" = "orange", "orfo" = "#F06CA8")

group_levels <- names(group_colors)

# Compute group means
pc_summary <- pc_scores %>%
  filter(group %in% group_levels) %>%
  group_by(group) %>%
  summarise(across(starts_with("PC"), mean), .groups = "drop") %>%
  mutate(group = factor(group, levels = group_levels))

# Prepare data for plotting
mean_points <- pc_summary %>%
  rename_with(~ gsub("_mean", "", .x)) %>%
  mutate(point_type = "mean", id = NA)

individual_points <- pc_scores %>%
  select(id, group, PC1, PC2, PC3) %>%
  mutate(point_type = "individual")

all_points <- bind_rows(individual_points, mean_points)

# 3D PCA plot
plot_ly(all_points,
        x = ~PC1, y = ~PC2, z = ~PC3,
        color = ~group, colors = group_colors,
        symbol = ~point_type, symbols = c("circle", "circle"),
        size = ~ifelse(point_type == "mean", 12, 5),
        opacity = ~ifelse(point_type == "mean", 1, 0.4),
        text = ~paste("ID:", id, "<br>Group:", group)) %>%
  layout(scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3")))

# 2D PCA plot 
ggplot() +
  geom_point(data = individual_points, aes(x = PC1, y = PC2, color = group), shape = 16, size = 2.5, alpha = 0.5) +
  geom_point(data = mean_points, aes(x = PC1, y = PC2, fill = group), shape = 21, size = 5, color = "black", stroke = 1.2) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(title = "PCA: PC1 vs PC2", x = "PC1", y = "PC2", color = "Group", fill = "Group") +
  theme_minimal()

# -------------------------------
# 7. PCA Loadings Heatmap
# -------------------------------

# Extract loadings for PC1-PC3
loadings <- as.data.frame(pca_res$rotation[, 1:3])
loadings$behavior <- rownames(loadings)

# Reshape for plotting
loadings_long <- loadings %>%
  pivot_longer(cols = starts_with("PC"),
               names_to = "Principal_Component",
               values_to = "Loading") %>%
  mutate(Principal_Component = factor(Principal_Component, levels = paste0("PC", 1:3)))

# Set max absolute loading for color scale
max_abs <- 0.4

# Plot heatmap
ggplot(loadings_long, aes(x = Principal_Component, y = behavior, fill = Loading)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    name = "Loading",
    low = "red", mid = "white", high = "blue",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    oob = squish
  ) +
  theme_minimal() +
  labs(
    title = "PCA Loadings for Behaviors (PC1-PC3)",
    x = "Principal Component",
    y = "Behavior"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10)
  )

# -------------------------------
# 8. Permutation Test on PC1 (you can manually modify this code and run again 
# for PC2, PC3)

set.seed(120)  # reproducibility

# Extract PC1 scores and group
pc_scores_perm <- data.frame(PC1 = pca_res$x[, 1],
                             group = combined_scaled$group)

# Observed F-statistic
observed_stat <- summary(aov(PC1 ~ group, data = pc_scores_perm))[[1]]["group", "F value"]

# Permutation test
n_perm <- 1000
perm_stats <- replicate(n_perm, {
  perm_group <- sample(pc_scores_perm$group)
  summary(aov(PC1 ~ perm_group, data = pc_scores_perm))[[1]]["perm_group", "F value"]
})

# p-value
# p_value <- mean(perm_stats >= observed_stat)
# Adjusted p-value to avoid reporting 0; reflects minimum detectable 
# significance from 1,000 permutations
p_value <- (sum(perm_stats >= observed_stat) + 1) / (n_perm + 1)

# Plot permutation distribution
hist(perm_stats, breaks = 30, col = "gray", main = "Permutation Test (PC1)",
     xlab = "F-statistic", xlim = range(c(perm_stats, observed_stat)))
abline(v = observed_stat, col = "red", lwd = 2)
legend("topright", legend = paste("Observed F =", round(observed_stat, 2),
                                  "\np-value =", round(p_value, 4)),
       bty = "n")
