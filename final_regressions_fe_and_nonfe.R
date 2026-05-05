library(dplyr)

model_data <- dc_pheno_climate_nlcd_imp %>%
  filter(
    !is.na(sos_early),
    !is.na(summer_shoulder),
    !is.na(aut_shoulder),
    !is.na(eos_late),
    !is.na(spring_tavg),
    !is.na(spring_prcp),
    !is.na(DBH),
    !is.na(SpCode),
    !is.na(impervious_90m),
    !is.na(nlcd_landcover_class)
  ) %>%
  mutate(
    SpCode = as.factor(SpCode),
    run_year = as.factor(run_year),
    nlcd_landcover_class = as.factor(nlcd_landcover_class)
  )

nrow(model_data)


write.csv(model_data, "model_data_final_with_impervious.csv", row.names = FALSE)
saveRDS(model_data, "model_data_final_with_impervious.rds")

#Baseline

model_sos <- lm(
  sos_early ~ spring_tavg + spring_prcp + impervious_90m + DBH + SpCode,
  data = model_data
)

summary(model_sos)



#Fuller Baseline w/ NLCD

model_sos_nlcd <- lm(
  sos_early ~ spring_tavg + spring_prcp + impervious_90m + DBH +
    SpCode + nlcd_landcover_class,
  data = model_data
)

summary(model_sos_nlcd)



#Fixed Effects Model (to check for robustness and controlling for differences across years)

model_sos_yearFE <- lm(
  sos_early ~ spring_tavg + spring_prcp + impervious_90m + DBH +
    SpCode + nlcd_landcover_class + run_year,
  data = model_data
)

summary(model_sos_yearFE)





# Summer shoulder model
model_summer <- lm(
  summer_shoulder ~ spring_tavg + spring_prcp + impervious_90m + DBH + SpCode,
  data = model_data
)

summary(model_summer)


# Autumn shoulder model
model_autumn <- lm(
  aut_shoulder ~ spring_tavg + spring_prcp + impervious_90m + DBH + SpCode,
  data = model_data
)

summary(model_autumn)


# EOS model
model_eos <- lm(
  eos_late ~ spring_tavg + spring_prcp + impervious_90m + DBH + SpCode,
  data = model_data
)

summary(model_eos)


### Robustness Check w/ FE for all models


model_summer_yearFE <- lm(
  summer_shoulder ~ spring_tavg + spring_prcp + impervious_90m + DBH +
    SpCode + nlcd_landcover_class + run_year,
  data = model_data
)

model_autumn_yearFE <- lm(
  aut_shoulder ~ spring_tavg + spring_prcp + impervious_90m + DBH +
    SpCode + nlcd_landcover_class + run_year,
  data = model_data
)

model_eos_yearFE <- lm(
  eos_late ~ spring_tavg + spring_prcp + impervious_90m + DBH +
    SpCode + nlcd_landcover_class + run_year,
  data = model_data
)

summary(model_summer_yearFE)
summary(model_autumn_yearFE)
summary(model_eos_yearFE)





library(modelsummary)


models <- list(
  "SOS" = model_sos,
  "Summer" = model_summer,
  "Autumn" = model_autumn,
  "EOS" = model_eos,
  "SOS (FE)" = model_sos_yearFE,
  "Summer (FE)" = model_summer_yearFE,
  "Autumn (FE)" = model_autumn_yearFE,
  "EOS (FE)" = model_eos_yearFE
)

modelsummary(
  models,
  stars = TRUE,
  statistic = "std.error",
  coef_omit = "run_year",   # 🔥 removes ALL year FE rows (including NA ones)
  gof_omit = "AIC|BIC|Log.Lik",
  output = "regression_table_clean.docx"
)




modelsummary(
  models,
  stars = TRUE,
  statistic = "std.error",
  coef_omit = "run_year",
  gof_omit = "AIC|BIC|Log.Lik",
  add_rows = data.frame(
    term = c("Species FE", "Year FE"),
    "SOS" = c("Yes", "No"),
    "Summer" = c("Yes", "No"),
    "Autumn" = c("Yes", "No"),
    "EOS" = c("Yes", "No"),
    "SOS (FE)" = c("Yes", "Yes"),
    "Summer (FE)" = c("Yes", "Yes"),
    "Autumn (FE)" = c("Yes", "Yes"),
    "EOS (FE)" = c("Yes", "Yes")
  ),
  output = "regression_table_clean.docx"
)


modelsummary(
  models,
  stars = TRUE,
  statistic = "std.error",
  coef_omit = "run_year",
  output = "regression_table.xlsx"
)