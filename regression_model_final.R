library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)


model_data <- dc_pheno_climate_nlcd %>%
  mutate(
    SpCode = as.factor(SpCode),
    nlcd_class = as.factor(nlcd_class)
  ) %>%
  drop_na(sos_early, mean_tavg_calc, total_prcp, DBH, SpCode, nlcd_class)

nrow(model_data)

model_data %>%
  distinct(SpCode) %>%
  nrow()


###First Model 

model1 <- lm(
  sos_early ~ mean_tavg_calc + total_prcp + DBH + SpCode + nlcd_class,
  data = model_data
)

summary(model1)


###Second Model (more detailed)

model2 <- lm(
  sos_early ~ mean_tavg_calc * SpCode + total_prcp * SpCode + 
    DBH + nlcd_class,
  data = model_data
)

summary(model2)

### Both Summaries 

summary(model1)
summary(model2)

### Overall Temp Effect

ggplot(model_data, aes(x = mean_tavg_calc, y = sos_early)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Relationship Between Mean Temperature and Start of Season",
    x = "Mean Temperature",
    y = "Start of Season (SOS Early)"
  ) +
  theme_minimal()


### Species Specific 


ggplot(model_data, aes(x = mean_tavg_calc, y = sos_early, color = SpCode)) +
  geom_point(alpha = 0.08, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Species-Specific Temperature Sensitivity of Start of Season",
    x = "Mean Temperature",
    y = "Start of Season (SOS Early)",
    color = "Species"
  ) +
  theme_minimal()


###Cleaner Version 

ggplot(model_data, aes(x = mean_tavg_calc, y = sos_early)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ SpCode) +
  labs(
    title = "Temperature Effects on Start of Season by Species",
    x = "Mean Temperature",
    y = "Start of Season (SOS Early)"
  ) +
  theme_minimal()


### Precipitation Plot

ggplot(model_data, aes(x = total_prcp, y = sos_early)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Relationship Between Precipitation and Start of Season",
    x = "Total Precipitation",
    y = "Start of Season (SOS Early)"
  ) +
  theme_minimal()

### Specie Specific  (temp on sos)

ggplot(model_data, aes(x = mean_tavg_calc, y = sos_early)) +
  geom_point(alpha = 0.12, size = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ SpCode) +
  labs(
    title = "Species-Specific Relationship Between Temperature and Start of Season",
    x = "Mean Temperature",
    y = "Start of Season (SOS Early)"
  ) +
  theme_minimal()

ggsave("species_temperature_sos_plot.png", width = 12, height = 10, dpi = 300)



### Cleaner Visual

library(broom)

# Extract coefficients from model2
coefs <- tidy(model2)

# Get temperature interaction terms
temp_slopes <- coefs %>%
  filter(grepl("mean_tavg_calc:SpCode", term)) %>%
  mutate(
    SpCode = gsub("mean_tavg_calc:SpCode", "", term),
    slope = estimate
  )

# Add baseline slope
baseline <- coefs %>%
  filter(term == "mean_tavg_calc") %>%
  pull(estimate)

temp_slopes <- temp_slopes %>%
  mutate(total_slope = slope + baseline)

# Plot
ggplot(temp_slopes, aes(x = reorder(SpCode, total_slope), y = total_slope)) +
  geom_point() +
  coord_flip() +
  labs(
    title = "Species Differences in Temperature Sensitivity",
    x = "Species",
    y = "Effect of Temperature on SOS (days per °C)"
  ) +
  theme_minimal()

###tem


coefs <- tidy(model2)

# baseline slope
baseline <- coefs %>%
  filter(term == "mean_tavg_calc") %>%
  pull(estimate)

# species-specific slopes
temp_slopes <- coefs %>%
  filter(grepl("mean_tavg_calc:SpCode", term)) %>%
  mutate(
    SpCode = gsub("mean_tavg_calc:SpCode", "", term),
    slope = estimate
  ) %>%
  mutate(
    total_slope = slope + baseline
  )

ggplot(temp_slopes, aes(x = reorder(SpCode, total_slope), y = total_slope)) +
  geom_point(size = 3) +
  geom_hline(yintercept = baseline, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Species Differences in Temperature Sensitivity",
    subtitle = "Effect of temperature on start of season (SOS)",
    x = "Species",
    y = "Change in SOS (days per °C)"
  ) +
  theme_minimal(base_size = 13)


### Top 5 Speices


top_species <- temp_slopes %>%
  arrange(total_slope) %>%
  slice(1:5)

bottom_species <- temp_slopes %>%
  arrange(desc(total_slope)) %>%
  slice(1:5)

top_species
bottom_species




temp_slopes <- temp_slopes %>%
  mutate(
    group = case_when(
      SpCode %in% top_species$SpCode ~ "Most Sensitive",
      SpCode %in% bottom_species$SpCode ~ "Least Sensitive",
      TRUE ~ "Other"
    )
  )

ggplot(temp_slopes, aes(x = reorder(SpCode, total_slope), y = total_slope, color = group)) +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    title = "Variation in Temperature Sensitivity Across Species",
    x = "Species",
    y = "Change in SOS (days per °C)",
    color = ""
  ) +
  theme_minimal(base_size = 13)





library(ggrepel)

ggplot(temp_slopes, aes(x = reorder(SpCode, total_slope), y = total_slope, color = group)) +
  geom_point(size = 3) +
  geom_hline(yintercept = baseline, linetype = "dashed") +
  geom_text_repel(
    data = temp_slopes %>% filter(group != "Other"),
    aes(label = SpCode),
    size = 4
  ) +
  coord_flip() +
  labs(
    title = "Species Differences in Temperature Sensitivity",
    subtitle = "Highlighted species show strongest and weakest responses",
    x = "Species",
    y = "Change in SOS (days per °C)",
    color = ""
  ) +
  theme_minimal(base_size = 13)