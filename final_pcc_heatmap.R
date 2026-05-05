### Species Heatmap ###

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tibble)

# 1) Pick the top 12 species by number of observations
top12_species <- dc_tree_pheno_climate %>%
  count(SpCode, sort = TRUE) %>%
  slice_head(n = 12) %>%
  pull(SpCode)

# 2) Define the heatmap rows

row_spec <- tribble(
  ~phase,      ~row_label,         ~phenology_var,   ~climate_var,
  "Greenup",   "spring_tavg",      "sos_early",      "spring_tavg",
  "Greenup",   "spring_prcp",      "sos_early",      "spring_prcp",
  "Greenup",   "winter_tavg",      "sos_early",      "winter_tavg",
  "Greenup",   "winter_snow",      "sos_early",      "winter_snow",
  
  "Maturity",  "spring_tavg",      "summer_shoulder","spring_tavg",
  "Maturity",  "spring_prcp",      "summer_shoulder","spring_prcp",
  "Maturity",  "winter_tavg",      "summer_shoulder","winter_tavg",
  "Maturity",  "winter_snow",      "summer_shoulder","winter_snow",
  
  "Senescence","spring_tavg",      "aut_shoulder",   "spring_tavg",
  "Senescence","spring_prcp",      "aut_shoulder",   "spring_prcp",
  "Senescence","winter_tavg",      "aut_shoulder",   "winter_tavg",
  "Senescence","winter_snow",      "aut_shoulder",   "winter_snow",
  
  "Dormancy",  "spring_tavg",      "eos_late",       "spring_tavg",
  "Dormancy",  "spring_prcp",      "eos_late",       "spring_prcp",
  "Dormancy",  "winter_tavg",      "eos_late",       "winter_tavg",
  "Dormancy",  "winter_snow",      "eos_late",       "winter_snow"
)

# 3) Correlation helper
safe_cor <- function(df, x, y) {
  ok <- complete.cases(df[[x]], df[[y]])
  if (sum(ok) < 3) return(NA_real_)
  cor(df[[x]][ok], df[[y]][ok], method = "pearson")
}

# 4) Compute correlations for each species + overall mean
plot_species <- c(top12_species, "Mean")

corr_heatmap <- crossing(
  row_spec,
  SpCode_plot = plot_species
) %>%
  rowwise() %>%
  mutate(
    correlation = {
      dat <- if (SpCode_plot == "Mean") {
        dc_tree_pheno_climate %>%
          filter(SpCode %in% top12_species)
      } else {
        dc_tree_pheno_climate %>%
          filter(SpCode == SpCode_plot)
      }
      safe_cor(dat, phenology_var, climate_var)
    }
  ) %>%
  ungroup()

# 5) Order axes
corr_heatmap <- corr_heatmap %>%
  mutate(
    SpCode_plot = factor(SpCode_plot, levels = plot_species),
    row_label = factor(row_label, levels = rev(unique(row_spec$row_label))),
    phase = factor(phase, levels = c("Greenup", "Maturity", "Senescence", "Dormancy"))
  )

# 6) Plot
ggplot(corr_heatmap, aes(x = SpCode_plot, y = row_label, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
  facet_grid(rows = vars(phase), scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "red3",
    mid = "white",
    high = "darkgreen",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90"
  ) +
  labs(
    title = "Pearson correlations between phenology and climate by species",
    x = NULL,
    y = NULL,
    fill = "r"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y.right = element_text(angle = 0, face = "bold"),
    strip.background = element_blank()
  )