# PCA plot and loadings plots 

# Load libraries
library(tidyverse)
library(plotly)
library(scales)

# Set working directory
setwd("C:/Users/Owner/Desktop/ML")

# -------------------------------
# Load and combine datasets

# Time spent performing grouped behavior 
behavior_data <- read.csv("behavior_time.csv") %>% 
  filter(group %in% c("dan_sd", "lira_sd", "orfo", "LiCl", "fed")) %>%
  mutate(total_move = move.explore.sec + move.explore...shelter.sec) %>%
  select(-move.explore.sec, -move.explore...shelter.sec)

# Grouped behavior transitions 
transitions <- read_csv("transitions.csv") %>%
  filter(group %in% c("dan_sd", "lira_sd", "orfo", "LiCl", "fed"))

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

transitions_combined <- transitions %>% # rare transitions
  select(id, group, ends_with("_combined")) %>%
  select(-matches("groom_drink|drink_groom|groom_food motivated|food motivated_groom"))

# Grouped behavior bouts
bouts <- read_csv("bouts.csv") %>%
  filter(!id %in% c(24, 11, 93),
         group %in% c("dan_sd", "lira_sd", "orfo", "LiCl", "fed")) %>%
  select(-starts_with("max")) # variable, just use average and # bouts

# -------------------------------
# Merge and process dataset

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
# PCA

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
# PCA Visualization

# Recode group names (remove '_sd')
pc_scores <- pc_scores %>%
  mutate(group = recode(group, "dan_sd" = "dan", "lira_sd" = "lira"))

# Assign colors for PCA
group_colors <- c("dan" = "#ED1C24", "lira" = "#26A7DF", "LiCl" = "#6ebe45", 
                  "fed" = "orange", "orfo" = "#F06CA8")

group_levels <- names(group_colors)

# Group means
pc_summary <- pc_scores %>%
  filter(group %in% group_levels) %>%
  group_by(group) %>%
  summarise(across(starts_with("PC"), mean), .groups = "drop") %>%
  mutate(group = factor(group, levels = group_levels))

# Data for plotting
mean_points <- pc_summary %>%
  rename_with(~ gsub("_mean", "", .x)) %>%
  mutate(point_type = "mean", id = NA)

individual_points <- pc_scores %>%
  dplyr::select(id, group, PC1, PC2, PC3) %>%
  mutate(point_type = "individual")

all_points <- bind_rows(individual_points, mean_points)

# 2D plot 
ggplot() +
  geom_point(data = individual_points, aes(x = PC1, y = PC2, color = group), shape = 16, size = 2.5, alpha = 0.5) +
  geom_point(data = mean_points, aes(x = PC1, y = PC2, fill = group), shape = 21, size = 5, color = "black", stroke = 1.2) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(title = "PCA: PC1 vs PC2", x = "PC1", y = "PC2", color = "Group", fill = "Group") +
  theme_minimal()

# -------------------------------
# Add loadings
loadings <- as.data.frame(pca_res$rotation[, 1:2]) 

# Convert to matrix with feature names
loadings_mat <- as.matrix(loadings[, c("PC1", "PC2")])
rownames(loadings_mat) <- rownames(loadings)  

# Function to get top pos and top neg per PC
get_top_pos_neg <- function(pc_vector, n = 1) { 
  if (is.null(names(pc_vector))) stop("pc_vector must have feature names")
  
  top_pos <- names(sort(pc_vector, decreasing = TRUE))[1:n]
  top_neg <- names(sort(pc_vector, decreasing = FALSE))[1:n]
  
  unique(c(top_pos, top_neg))  # return feature names
}

pc1_vec <- loadings_mat[, "PC1"]
pc2_vec <- loadings_mat[, "PC2"]

# Top features for PC1 and PC2
top_pc1 <- get_top_pos_neg(pc1_vec)
top_pc2 <- get_top_pos_neg(pc2_vec)

# Combine 
top_features <- unique(c(top_pc1, top_pc2))

top_loadings <- loadings_mat[top_features, c("PC1", "PC2")]

# Rescale loadings to fit plot
mult <- min(
  (max(individual_points$PC1) - min(individual_points$PC1)) / 
    (max(top_loadings[, "PC1"]) - min(top_loadings[, "PC1"])),
  (max(individual_points$PC2) - min(individual_points$PC2)) / 
    (max(top_loadings[, "PC2"]) - min(top_loadings[, "PC2"]))
) * 0.4

top_loadings_scaled <- top_loadings * mult
top_loadings_scaled <- as.data.frame(top_loadings_scaled)
top_loadings_scaled$feature <- rownames(top_loadings_scaled)

# Plot
ggplot() +
  geom_point(data = individual_points, aes(x = PC1, y = PC2, color = group), 
             shape = 16, size = 2.5, alpha = 0.5) +
  geom_point(data = mean_points, aes(x = PC1, y = PC2, fill = group), 
             shape = 21, size = 5, color = "black", stroke = 1.2) +
  geom_segment(data = top_loadings_scaled, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black") +
  geom_text(data = top_loadings_scaled, 
            aes(x = PC1, y = PC2, label = feature),
            vjust = 1.5, color = "black") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(title = "PCA: PC1 vs PC2 with Top Loadings", 
       x = "PC1", y = "PC2", color = "Group", fill = "Group") +
  theme_minimal()

# -------------------------------
# PCA loadings heatmap

# Loadings for PC1-PC3
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