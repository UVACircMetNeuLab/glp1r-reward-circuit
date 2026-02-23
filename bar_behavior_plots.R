# Behavior Analysis - Paper Figures

# Figure overview:
#   - Fig. 2 (main): Normalized bar plots pairing each vehicle with its
#     drug (saline_lira/lira_sd, veh_sd/dan_sd, saline_orfo/orfo,
#     saline_licl/LiCl) for select behaviors (and more behaviors in extended). 
#     Each drug group is scaled to its own vehicle mean for visualization. 

#   - Extended figures cont.: Raw bar plots for drug-only groups + fed control
   
#     NOTE: IDs 1789 and 1836 are excluded ONLY from wheel cycle plots
#           because they lack wheel data but are included in other behaviors.

library(patchwork)
library(tidyverse)
library(ggplot2)
library(writexl)

setwd("C:/Users/Owner/Desktop/ML")

# Load and filter sensor data
home_cage <- read_csv("sensor_data_dark_cycle.csv") %>%
  filter(Genotype == "Homo") %>% # should all be homo in this data anyway
  filter(Treatment %in% c(
    "Saline_orfo", "Orfo",
    "Vehicle_SD",  "Dan_SD",
    "Saline_lira", "Lira",
    "Saline_licl", "LiCl", 
    "Fed"
  ))

# Per-drug inclusion lists (animals with complete/valid data)
keep_ids <- list(
  orfo      = c(96, 1710, 1788, 1789, 1790, 1836, 1837, 2784, 2782, 2785),
  dan       = c(8, 69, 253, 287, 712, 718, 720, 786, 1006),
  lira      = c(4, 5, 7, 8, 13, 1207, 1208, 1613, 1615),
  licl      = c(1790, 27, 20, 26, 23, 1710),
  fed       = c(1, 2, 3, 9, 10, 12, 13, 23, 25, 26) 

)

home_cage <- home_cage %>%
  filter(
    !(Treatment %in% c("Saline_orfo", "Orfo"))       | ID %in% keep_ids$orfo,
    !(Treatment %in% c("Dan_SD", "Vehicle_SD")) | ID %in% keep_ids$dan,
    !(Treatment %in% c("Lira", "Saline_lira")) | ID %in% keep_ids$lira,
    !(Treatment %in% c("LiCl", "Saline_licl")) | ID %in% keep_ids$licl,
    !(Treatment %in% c("Fed")) | ID %in% keep_ids$fed
    
  )

# Make sure sample sizes are correct!
home_cage %>%
  dplyr::group_by(Treatment) %>%
  dplyr::summarise(n_unique_ids = n_distinct(ID)) %>%
  print()

# Aggregate sensor data (first 2 hours)
#detach("package:plyr", unload = TRUE)
twohours <- home_cage %>%
  filter(Time_interval %in% c("1hr", "2hr"))

# Aggregate
sensor_twohours <- twohours %>%
  group_by(ID, Treatment) %>%
  summarise(
    Lickometer      = sum(Lickometer,      na.rm = TRUE),
    Food_nose_pokes = sum(Food_nose_pokes, na.rm = TRUE),
    Wheel_cycles    = sum(Wheel_cycles,    na.rm = TRUE),
    .groups = "drop"
  )

# Check we didn't drop columns we need
names(sensor_twohours)

# Load and clean time spent performing behaviors
behavior_data <- read.csv("behavior_time_dark_cycle.csv") %>%
  rename(
    total_drink        = 3,
    total_food         = 4,
    total_groom        = 5,
    total_move         = 6,
    total_move_shelter = 7,
    total_rest         = 8
  ) %>%
  filter(!id %in% c(24, 11, 93)) %>% 
  filter(group %in% c(
    "LiCl",       "dan_sd",   "lira_sd",  "orfo",
    "saline_licl", "veh_sd",  "saline_orfo", "saline_lira", "fed"
  ))

# Normalize each row to proportion of total observed time, then convert to %
behavior_data <- behavior_data %>%
  rowwise() %>%
  mutate(total_time = sum(c_across(c(
    total_drink, total_food, total_groom,
    total_move, total_move_shelter, total_rest
  )), na.rm = TRUE)) %>%
  mutate(across(
    c(total_drink, total_food, total_groom, total_move, total_move_shelter, total_rest),
    ~ (.x / total_time) * 100
  )) %>%
  ungroup() %>%
  # Combine cage floor and shelter movement into one column
  dplyr::mutate(total_move = total_move + total_move_shelter) %>%
  dplyr::select(-total_move_shelter, -total_time)

# Sanity check
behavior_data %>%
  group_by(group) %>%
  summarise(n_mice = n_distinct(id)) %>%
  print()

# Combine behavior + sensor data
names(sensor_twohours)[names(sensor_twohours) == "ID"] <- "id"
names(sensor_twohours)[names(sensor_twohours) == "Treatment"] <- "group"

# Make names match 
sensor_twohours$group <- dplyr::recode(sensor_twohours$group,
                                       "Vehicle_SD"  = "veh_sd",
                                       "Dan_SD"      = "dan_sd",
                                       "Lira"        = "lira_sd",
                                       "Saline_lira" = "saline_lira",
                                       "LiCl"        = "LiCl",
                                       "Saline_licl" = "saline_licl",
                                       "Orfo"        = "orfo",
                                       "Saline_orfo" = "saline_orfo",
                                       "Fed"         = "fed" 
)

behavior_data$id <- as.character(behavior_data$id)
sensor_twohours$id <- as.character(sensor_twohours$id)

combined_behavior <- left_join(behavior_data, sensor_twohours, 
                               by = c("id", "group"))

# -----------------------------------------------------------------------------
# Shared plot aesthetics
group_colors <- c(
  "veh_sd"      = "black",
  "dan_sd"      = "red",
  "saline_orfo" = "grey30",
  "orfo"        = "hotpink",
  "saline_lira" = "#00008B",
  "lira_sd"     = "#1E90FF",
  "saline_licl" = "grey70",
  "LiCl"        = "#66CD00",
  "fed"         = "#FFA500"
)

# Shared plot theme
theme_pub <- function() {
  theme_minimal(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      legend.position = "none"
    )
}

# Helper: bar + SEM + jitter for a single behavior column
make_bar_plot <- function(df, behavior_col, y_label = behavior_col,
                          group_order = NULL) {
  if (!is.null(group_order)) {
    df$group <- factor(df$group, levels = group_order)
  }
  ggplot(df, aes(x = group, y = .data[[behavior_col]],
                 fill = group, color = group)) +
    stat_summary(fun = mean, geom = "bar",
                 alpha = 0.5, width = 0.6, linewidth = 1) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 width = 0.2, linewidth = 0.8) +
    geom_jitter(width = 0.15, size = 2, show.legend = FALSE) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    theme_pub() +
    ylab(y_label)
}

# ------------------------------------------------------------------------------
# Figure 2: Normalized (vehicle-relative) plots
#    Behaviors: Food_nose_pokes, food motivated, drink, move, rest/groom in
#    the shelter.

#    Modify norm_cols for plotting other behaviors, such as those in
#    extended fig 5a-5e.

# Scale each drug + its vehicle to the vehicle mean
group_pairs <- list(
  saline_lira  = "lira_sd",
  veh_sd       = "dan_sd",
  saline_orfo  = "orfo",
  saline_licl  = "LiCl"
)

# Columns to normalize 
norm_cols <- c("Food_nose_pokes", "total_food", "total_drink",
               "total_move", "total_rest")

scaled_list <- lapply(names(group_pairs), function(ctrl) {
  drug <- group_pairs[[ctrl]]
  ctrl_df <- combined_behavior %>% filter(group == ctrl)
  drug_df <- combined_behavior %>% filter(group == drug)
  if (nrow(ctrl_df) == 0 || nrow(drug_df) == 0) return(NULL)
  
  # Only normalize columns that exist and have non-zero means
  cols_present <- intersect(norm_cols, names(combined_behavior))
  ctrl_means   <- colMeans(ctrl_df[, cols_present], na.rm = TRUE)
  ctrl_means[ctrl_means == 0] <- NA  
  
  ctrl_df[, cols_present] <- sweep(ctrl_df[, cols_present], 2, ctrl_means, "/")
  drug_df[, cols_present] <- sweep(drug_df[, cols_present], 2, ctrl_means, "/")
  bind_rows(ctrl_df, drug_df)
})

scaled_df <- bind_rows(scaled_list)

# Save scaled data for sharing / archiving
write_xlsx(scaled_df, "fig2_behavior_scaled.xlsx")

# Make main plots  
drug_order <- c(
  "veh_sd",      "dan_sd",
  "saline_orfo", "orfo",
  "saline_lira", "lira_sd",
  "saline_licl", "LiCl"
)

main_behaviors <- c("Food_nose_pokes", "total_food", "total_drink",
                          "total_move", "total_rest")

main_plots <- lapply(main_behaviors, function(b) {
  # For wheel cycles: exclude mice missing wheel data
  df <- if (b == "Wheel_cycles") {
    scaled_df %>% filter(!id %in% c("1789", "1836")) 
  } else {
    scaled_df
  }
  make_bar_plot(df, b, y_label = b, group_order = drug_order)
})

combined_main_plot <- wrap_plots(main_plots, nrow = 1)
print(combined_main_plot)

ggsave("fig_main_normalized.png", plot = combined_main_plot, # eps/svg for adobe
       device = "png", width = 24, height = 4, units = "in", dpi = 300)

#-----------------------------------------------------------------------------
# Extended Figures 5g-n: Raw data, drug groups only (+ fed control)
# Groups shown: fed, LiCl, lira_sd, dan_sd, orfo

raw_order  <- c("fed", "LiCl", "lira_sd", "dan_sd", "orfo")

raw_df <- combined_behavior %>%
  filter(group %in% raw_order)

raw_behaviors <- c("total_rest","Food_nose_pokes", 
                   "Lickometer", "Wheel_cycles") # modify for other behaviors

extend_fig_plots <- lapply(raw_behaviors, function(b) {
  # For wheel cycles: exclude mice missing wheel data
  df <- if (b == "Wheel_cycles") {
    raw_df %>% filter(!id %in% c("1789", "1836"))
  } else {
    raw_df
  }
  make_bar_plot(df, b, y_label = b, group_order = raw_order)
})

combined_extend_fig <- wrap_plots(extend_fig_plots, nrow = 1)
print(combined_extend_fig)

ggsave("fig_raw_bar_extended.png", plot = combined_extend_fig, #eps/svg for adobe
       device = "png", width = 24, height = 4, units = "in", dpi = 300)


# -----------------------------------------------------------------------------
# Transitions, Bouts, & Distance Traveled (extended fig 5f, 6)

# Load & clean distance (paired vehicle + drug, raw/unscaled)
distance_data <- read.csv("distance_new.csv") %>%
  rename(
    id       = 1,
    group    = 2,
    distance = 3
  ) 
# fix naming
distance_data$group <- dplyr::recode(distance_data$group,
                                       "veh"  = "veh_sd",
                                       "dan"      = "dan_sd",
                                       "lira"        = "lira_sd",
                                       "licl"        = "LiCl"
                                     )

distance_data <- distance_data %>%  # no fed for this plot
  filter(group %in% c(  "veh_sd",      "dan_sd",
                        "saline_orfo", "orfo",
                        "saline_lira", "lira_sd",
                        "saline_licl", "LiCl"))

# Plot distance - paired order (vehicle next to drug)
distance_plot <- make_bar_plot(distance_data, "distance",
                               y_label = "Distance travelled (cm)",
                               group_order = drug_order)  
print(distance_plot)

ggsave("fig_distance.png", plot = distance_plot,
       device = "png", width = 6, height = 4, units = "in", dpi = 300)

# Load & clean transitions (drug groups + fed)
transitions <- read_csv("norm_transitions_dark_cycle.csv") %>%
  filter(group %in% c("dan_sd", "lira_sd", "orfo", "LiCl", "fed"))

# Combine reciprocal transition pairs (e.g. groom->drink + drink->groom)
reciprocal_pairs <- list(
  c("groom_drink",                 "drink_groom"),
  c("food motivated_move/explore", "move/explore_food motivated"),
  c("groom_food motivated",        "food motivated_groom"),
  c("groom_move/explore",          "move/explore_groom"),
  c("move/explore_drink",          "drink_move/explore"),
  c("move/explore_shelter",        "shelter_move/explore")
)

for (pair in reciprocal_pairs) {
  transitions[[paste0(pair[1], "_combined")]] <- rowSums(transitions[, pair], na.rm = TRUE)
}

# Keep only combined columns
transitions_combined <- transitions %>%
  dplyr::select(id, group, ends_with("_combined")) 

# Load & clean bouts 
bouts <- read_csv("bouts_dark_cycle.csv") %>%
  filter(
    !id %in% c(24, 11, 93),
    group %in% c("dan_sd", "lira_sd", "orfo", "LiCl", "fed")
  ) 

# Plot transitions
transition_cols <- setdiff(names(transitions_combined), c("id", "group"))

transition_plots <- lapply(transition_cols, function(b) {
  make_bar_plot(transitions_combined, b, y_label = b, group_order = raw_order)
})

combined_transitions <- wrap_plots(transition_plots, nrow = 1)
print(combined_transitions)
ggsave("fig_transitions.png", plot = combined_transitions,
       device = "png", width = 24, height = 4, units = "in", dpi = 300)

# Plot bouts
bout_cols <- setdiff(names(bouts), c("id", "group"))

bout_plots <- lapply(bout_cols, function(b) {
  make_bar_plot(bouts, b, y_label = b, group_order = raw_order)
})

combined_bouts <- wrap_plots(bout_plots, nrow = 1)
print(combined_bouts)
ggsave("fig_bouts.png", plot = combined_bouts,
       device = "png", width = 24, height = 4, units = "in", dpi = 300)

# -----------------------------------------------------------------------------
# Stats note:

# Sensor data (Lickometer, Food_nose_pokes, Wheel_cycles), distance
# travelled: t-test or one-way ANOVA depending on comparison (see paper for 
# details).

# Behavioral time data (% time move, rest, etc.): GLMM with beta
# regression to account for bounded proportional data (code available, see 
# paper for details).

# All statistical tests were run on raw, un-normalized data.
# Normalization to vehicle mean is for visualization purposes only.
