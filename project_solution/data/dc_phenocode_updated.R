# Purpose: Fit logistics to each crown, calculate SOS/EOS/phenology metrics
# DC version
# Edited to:
#   - run one target year with winter buffers
#   - process UIDs in batches
#   - remove bootstrap loop
#   - remove graph saving
#   - keep real spline + double-logistic fitting
#   - preserve SpCode in failed rows when possible
#   - recompute final spline after logistic filtering

library(data.table)
library(doSNOW)
library(psych)
library(matrixStats)
library(minpack.lm)
library(quantreg)
library(sf)
library(dplyr)
library(caTools)

# 1 Setup ------------------------------------------------------------

# Output directory
output_dir <- "C:/Users/valea/OneDrive/Documents/dc_pheno_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Output name prefix (NO .csv here)
user_output_filename <- "dc_pheno_metrics_2022"

# Year to analyze
user_year <- 2023

# Buffered date window
start_date <- 20221201
end_date   <- 20240331

# Winter shoulder year target
winter_shoulders_years <- user_year

# Parallel settings
NumberOfCluster <- 6

# Sensor / image settings
user_sensor <- "pss4"
image_type <- "surface_reflectance"

# Quality / smoothing options
user_r_red <- 0.80
user_mad_red <- 400
quant_filt <- TRUE
df <- 10

# Optional shapefile join
do_join <- "no"
shapefile_filename <- NA
shp_merge_filename <- NA

# Optional fit_stats workflow
use_fit_stats <- FALSE

# Batch controls
batch_size   <- 5000
batch_number <-2

# 1.2 Helper functions -----------------------------------------------

apply_winter_shoulders <- function(tree_dat, years_vec) {
  if (length(years_vec) == 0) {
    return(tree_dat)
  }
  
  years_vec <- unique(as.integer(years_vec))
  
  for (y in years_vec) {
    if (is.na(y)) next
    
    prev_idx <- which(tree_dat$im_date < as.integer(paste0(y, "0000")))
    next_idx <- which(tree_dat$im_date > as.integer(paste0(y + 1, "0000")))
    
    if (length(prev_idx) > 0) {
      tree_dat$im_doy[prev_idx] <- tree_dat$im_doy[prev_idx] - 365
    }
    
    if (length(next_idx) > 0) {
      tree_dat$im_doy[next_idx] <- tree_dat$im_doy[next_idx] + 365
    }
  }
  
  tree_dat
}

calc_spline_metrics <- function(sm_spline) {
  maxNDVI <- max(sm_spline$y, na.rm = TRUE)
  max_green_doy_idx <- which.min(abs(sm_spline$y - maxNDVI))
  max_green_doy <- sm_spline$x[max_green_doy_idx]
  
  df_greenup_idx <- which(sm_spline$x < max_green_doy)
  minNDVI <- min(sm_spline$y[df_greenup_idx], na.rm = TRUE)
  
  greenup_15 <- (maxNDVI - minNDVI) * 0.15 + minNDVI
  greenup_50 <- (maxNDVI - minNDVI) * 0.50 + minNDVI
  greenup_90 <- (maxNDVI - minNDVI) * 0.90 + minNDVI
  
  minNDVI_DOY_idx <- which(sm_spline$y == minNDVI)
  minNDVI_DOY <- sm_spline$x[minNDVI_DOY_idx]
  df_greenup_idx <- which(sm_spline$x < max_green_doy & sm_spline$x > minNDVI_DOY)
  
  g15 <- approx(x = sm_spline$y[df_greenup_idx], y = sm_spline$x[df_greenup_idx], xout = greenup_15)$y
  g50 <- approx(x = sm_spline$y[df_greenup_idx], y = sm_spline$x[df_greenup_idx], xout = greenup_50)$y
  g90 <- approx(x = sm_spline$y[df_greenup_idx], y = sm_spline$x[df_greenup_idx], xout = greenup_90)$y
  
  df_browndown_idx <- which(sm_spline$x >= max_green_doy)
  minNDVI_bd <- min(sm_spline$y[df_browndown_idx], na.rm = TRUE)
  
  brown_10 <- (maxNDVI - minNDVI_bd) * 0.90 + minNDVI_bd
  brown_50 <- (maxNDVI - minNDVI_bd) * 0.50 + minNDVI_bd
  brown_85 <- (maxNDVI - minNDVI_bd) * 0.15 + minNDVI_bd
  
  b10 <- approx(x = sm_spline$y[df_browndown_idx], y = sm_spline$x[df_browndown_idx], xout = brown_10)$y
  b50 <- approx(x = sm_spline$y[df_browndown_idx], y = sm_spline$x[df_browndown_idx], xout = brown_50)$y
  b85 <- approx(x = sm_spline$y[df_browndown_idx], y = sm_spline$x[df_browndown_idx], xout = brown_85)$y
  
  list(
    maxNDVI = maxNDVI,
    minNDVI_pregreen = minNDVI,
    g15 = g15,
    g50 = g50,
    g90 = g90,
    b10 = b10,
    b50 = b50,
    b85 = b85
  )
}

v <- function(m, t) {
  m[1] + (m[2] - m[7] * t) *
    ((1 / (1 + exp((m[3] - t) / m[4]))) - (1 / (1 + exp((m[5] - t) / m[6]))))
}

# 2 Data prep ---------------------------------------------------------

dc_dat <- as.data.table(dc_data)

required_cols <- c("UID", "SpCode", "ndvi", "im_date", "im_doy", "year")
missing_cols <- setdiff(required_cols, names(dc_dat))

if (length(missing_cols) > 0) {
  stop(paste("dc_dat is missing required columns:", paste(missing_cols, collapse = ", ")))
}

dc_dat <- dc_dat[
  !is.na(UID) &
    !is.na(ndvi) &
    im_date >= start_date &
    im_date <= end_date
]

print(head(dc_dat))
print(dim(dc_dat))

if (use_fit_stats) {
  if (!exists("fit_stats_filename")) {
    stop("use_fit_stats is TRUE but fit_stats_filename is not defined.")
  }
  if (!exists("fit_stats_helper_path")) {
    stop("use_fit_stats is TRUE but fit_stats_helper_path is not defined.")
  }
  
  source(fit_stats_helper_path)
  fit_stats_summary <- load_fit_stats_summary(fit_stats_filename)
  dc_dat <- merge(dc_dat, fit_stats_summary, by.x = "im_datetime", by.y = "folder")
}

# 3 Parallel loop -----------------------------------------------------

uid_unique <- unique(dc_dat$UID)
uid_unique <- uid_unique[!is.na(uid_unique)]
uid_unique <- sort(uid_unique)

iterations <- length(uid_unique)

batch_starts <- seq(1, iterations, by = batch_size)
batch_ends   <- pmin(batch_starts + batch_size - 1, iterations)

if (batch_number > length(batch_starts)) {
  stop(paste0(
    "batch_number = ", batch_number,
    " is too large. Max batch_number is ", length(batch_starts)
  ))
}

batch_start_idx <- batch_starts[batch_number]
batch_end_idx   <- batch_ends[batch_number]

uid_batch <- uid_unique[batch_start_idx:batch_end_idx]
iter_subset <- seq_along(uid_batch)

cat("Total unique UIDs in filtered run window:", iterations, "\n")
cat("Running batch", batch_number, "of", length(batch_starts), "\n")
cat("Batch index range:", batch_start_idx, "to", batch_end_idx, "\n")
cat("UIDs in this batch:", length(uid_batch), "\n")

cl <- makeCluster(NumberOfCluster, outfile = "NUL")
registerDoSNOW(cl)

start_time <- Sys.time()

results_list <- foreach(
  i = iter_subset,
  .packages = c("data.table", "psych", "matrixStats",
                "minpack.lm", "quantreg", "splines", "dplyr", "caTools"),
  .export = c("apply_winter_shoulders", "calc_spline_metrics", "v")
) %dopar% {
  
  result_tc <- tryCatch(
    expr = {
      
      user_uid <- uid_batch[i]
      
      idx_uid_starter <- which(
        dc_dat$UID == user_uid &
          dc_dat$im_date >= start_date &
          dc_dat$im_date <= end_date
      )
      
      tree_dat <- dc_dat[idx_uid_starter, ]
      tree_dat <- tree_dat[!is.na(tree_dat$ndvi), ]
      
      if (nrow(tree_dat) < 8) {
        stop("Too few observations after initial filtering")
      }
      
      # Preserve species early in case the fit fails later
      this_spcode <- if ("SpCode" %in% names(tree_dat) && nrow(tree_dat) > 0) {
        unique(tree_dat$SpCode)[1]
      } else {
        NA_character_
      }
      
      # Scale NDVI
      tree_dat$ndvi <- as.integer(tree_dat$ndvi * 10000)
      
      # Drop duplicate dates, keeping highest absolute NDVI
      tree_dat <- tree_dat[order(tree_dat$im_date, -abs(tree_dat$ndvi)), ]
      dup_bin <- duplicated(tree_dat$im_date)
      tree_dat <- tree_dat[!dup_bin, ]
      
      if (nrow(tree_dat) < 8) {
        stop("Too few unique observation dates after deduplication")
      }
      
      # Apply winter shoulders
      tree_dat <- apply_winter_shoulders(tree_dat, winter_shoulders_years)
      
      # Sort for fitting
      tree_dat_sort <- tree_dat[order(tree_dat$im_doy), ]
      y0 <- tree_dat_sort$ndvi
      x0 <- tree_dat_sort$im_doy
      year0 <- tree_dat_sort$year
      
      image_dates_initial <- length(x0)
      
      # Optional running quantile smoothing
      if (quant_filt) {
        y0 <- runquantile(y0, k = 3, probs = 0.75, endrule = "quantile")
      }
      
      # Initial spline + spline-based outlier filtering
      x0_dith <- dither(x0)
      sm_spline <- smooth.spline(x0_dith, y0, cv = FALSE, df = df)
      fit_resids <- sm_spline$y - y0
      
      outliers <- boxplot(fit_resids, plot = FALSE)$out
      inliers <- setdiff(fit_resids, outliers)
      inliers_idx <- which(fit_resids %in% inliers)
      
      x0 <- x0[inliers_idx]
      y0 <- y0[inliers_idx]
      year0 <- year0[inliers_idx]
      
      if (length(x0) < 8) {
        stop("Too few observations after spline outlier filtering")
      }
      
      # Refit spline after spline-based filtering
      x0_dith <- dither(x0)
      sm_spline <- smooth.spline(x0_dith, y0, cv = FALSE, df = df)
      
      spl_metrics_op <- calc_spline_metrics(sm_spline)
      maxNDVI <- spl_metrics_op$maxNDVI
      minNDVI_pregreen <- spl_metrics_op$minNDVI_pregreen
      g15_spl_op <- spl_metrics_op$g15
      g50_spl_op <- spl_metrics_op$g50
      g90_spl_op <- spl_metrics_op$g90
      b10_spl_op <- spl_metrics_op$b10
      b50_spl_op <- spl_metrics_op$b50
      b85_spl_op <- spl_metrics_op$b85
      
      # Initial parameter guesses
      m1_idx <- which(x0 < 90 | x0 > 330)
      if (length(m1_idx) == 0) {
        m1_idx <- seq_along(x0)
      }
      
      m1 <- as.numeric(quantile(y0[m1_idx], 0.50, na.rm = TRUE))
      m2 <- as.numeric(quantile(y0, 0.99, na.rm = TRUE) - m1)
      m3 <- 110
      m4 <- 10
      m5 <- 310
      m6 <- 10
      m7 <- 4
      
      m <- c(m1, m2, m3, m4, m5, m6, m7)
      
      # First logistic fit
      mnew1 <- nlsLM(
        y0 ~ a + (b - g * x0) * ((1 / (1 + exp((c - x0) / d))) - (1 / (1 + exp((e - x0) / f)))),
        start = list(a = m[1], b = m[2], c = m[3], d = m[4], e = m[5], f = m[6], g = m[7])
      )
      mnew1_params <- mnew1$m$getPars()
      
      # Weights from residuals
      w0 <- max(abs(residuals(mnew1)))^2 - abs(residuals(mnew1))^2
      
      # Second weighted fit
      mnew2 <- nlsLM(
        y0 ~ a + (b - g * x0) * ((1 / (1 + exp((c - x0) / d))) - (1 / (1 + exp((e - x0) / f)))),
        start = list(
          a = mnew1_params[1], b = mnew1_params[2], c = mnew1_params[3], d = mnew1_params[4],
          e = mnew1_params[5], f = mnew1_params[6], g = mnew1_params[7]
        ),
        weights = w0
      )
      mnew2_params <- mnew2$m$getPars()
      
      # Logistic residual filtering
      mnew2_resids <- abs(residuals(mnew2))
      outliers <- boxplot(mnew2_resids, plot = FALSE)$out
      inliers <- setdiff(mnew2_resids, outliers)
      inliers_idx <- which(mnew2_resids %in% inliers)
      
      x0 <- x0[inliers_idx]
      y0 <- y0[inliers_idx]
      w0 <- w0[inliers_idx]
      
      if (length(x0) < 8) {
        stop("Too few observations after logistic outlier filtering")
      }
      
      # Final weighted fit
      mnew3 <- nlsLM(
        y0 ~ a + (b - g * x0) * ((1 / (1 + exp((c - x0) / d))) - (1 / (1 + exp((e - x0) / f)))),
        start = list(
          a = mnew2_params[1], b = mnew2_params[2], c = mnew2_params[3], d = mnew2_params[4],
          e = mnew2_params[5], f = mnew2_params[6], g = mnew2_params[7]
        ),
        weights = w0
      )
      mnew3_params <- mnew3$m$getPars()
      
      image_dates_final <- length(x0)
      
      # Recompute final spline on final retained observations
      x0_dith_final <- dither(x0)
      sm_spline_final <- smooth.spline(x0_dith_final, y0, cv = FALSE, df = df)
      fit_resids_final <- sm_spline_final$y - y0
      
      spl_metrics_final <- calc_spline_metrics(sm_spline_final)
      maxNDVI_final <- spl_metrics_final$maxNDVI
      minNDVI_pregreen_final <- spl_metrics_final$minNDVI_pregreen
      g15_spl_final <- spl_metrics_final$g15
      g50_spl_final <- spl_metrics_final$g50
      g90_spl_final <- spl_metrics_final$g90
      b10_spl_final <- spl_metrics_final$b10
      b50_spl_final <- spl_metrics_final$b50
      b85_spl_final <- spl_metrics_final$b85
      
      # Curvature-based DL dates
      x_interp <- seq(min(x0), max(x0), 0.1)
      ndvi_interp <- v(mnew3_params, x_interp)
      ndvi_curve <- diff(diff(ndvi_interp))
      
      early_season_idx <- which(x_interp < 180)
      late_season_idx  <- which(x_interp > 180)
      
      sos_concavity_idx <- which(ndvi_curve == max(ndvi_curve[early_season_idx], na.rm = TRUE))
      sos_convexity_idx <- which(ndvi_curve == min(ndvi_curve[early_season_idx], na.rm = TRUE))
      eos_concavity_idx <- which(ndvi_curve == max(ndvi_curve[late_season_idx], na.rm = TRUE))
      eos_convexity_idx <- which(ndvi_curve == min(ndvi_curve[late_season_idx], na.rm = TRUE))
      
      sos_early <- x_interp[sos_concavity_idx][1]
      summer_shoulder <- x_interp[sos_convexity_idx][1]
      aut_shoulder <- x_interp[eos_convexity_idx][1]
      eos_late <- x_interp[eos_concavity_idx][1]
      
      # Full-year fitted DL curve
      x_365 <- seq(1, 365)
      mnew3_fitted_pts <- data.frame(
        doy_365 = x_365,
        ndvi = v(mnew3_params, x_365)
      )
      
      # DL green-up metrics
      green_up_pts <- mnew3_fitted_pts %>% filter(doy_365 < 180)
      minNDVI_green <- min(green_up_pts$ndvi, na.rm = TRUE)
      maxNDVI_green <- max(green_up_pts$ndvi, na.rm = TRUE)
      
      greenup_15 <- (maxNDVI_green - minNDVI_green) * 0.15 + minNDVI_green
      greenup_50 <- (maxNDVI_green - minNDVI_green) * 0.50 + minNDVI_green
      greenup_90 <- (maxNDVI_green - minNDVI_green) * 0.90 + minNDVI_green
      
      g15_dl <- approx(x = green_up_pts$ndvi, y = green_up_pts$doy_365, xout = greenup_15)$y
      g50_dl <- approx(x = green_up_pts$ndvi, y = green_up_pts$doy_365, xout = greenup_50)$y
      g90_dl <- approx(x = green_up_pts$ndvi, y = green_up_pts$doy_365, xout = greenup_90)$y
      
      # DL brown-down metrics
      brown_down_pts <- mnew3_fitted_pts %>% filter(doy_365 >= 180)
      minNDVI_brown <- min(brown_down_pts$ndvi, na.rm = TRUE)
      maxNDVI_brown <- max(brown_down_pts$ndvi, na.rm = TRUE)
      
      brown_10 <- (maxNDVI_brown - minNDVI_brown) * 0.90 + minNDVI_brown
      brown_50 <- (maxNDVI_brown - minNDVI_brown) * 0.50 + minNDVI_brown
      brown_85 <- (maxNDVI_brown - minNDVI_brown) * 0.15 + minNDVI_brown
      
      b10_dl <- approx(x = brown_down_pts$ndvi, y = brown_down_pts$doy_365, xout = brown_10)$y
      b50_dl <- approx(x = brown_down_pts$ndvi, y = brown_down_pts$doy_365, xout = brown_50)$y
      b85_dl <- approx(x = brown_down_pts$ndvi, y = brown_down_pts$doy_365, xout = brown_85)$y
      
      # AUC
      auc_dl <- sum(mnew3_fitted_pts$ndvi, na.rm = TRUE) / nrow(mnew3_fitted_pts)
      auc_spl <- sum(sm_spline_final$y, na.rm = TRUE) / length(sm_spline_final$y)
      
      # Fit error metrics
      dl_fit <- v(mnew3_params, x0)
      dl_fit_resids <- dl_fit - y0
      rmse_dl <- sqrt(mean(dl_fit_resids^2, na.rm = TRUE))
      mad_dl <- median(abs(dl_fit_resids), na.rm = TRUE)
      
      rmse_spl <- sqrt(mean(fit_resids_final^2, na.rm = TRUE))
      mad_spl <- median(abs(fit_resids_final), na.rm = TRUE)
      
      data.frame(
        UID = user_uid,
        SpCode = this_spcode,
        
        m1 = mnew3_params[1],
        m2 = mnew3_params[2],
        m3 = mnew3_params[3],
        m4 = mnew3_params[4],
        m5 = mnew3_params[5],
        m6 = mnew3_params[6],
        m7 = mnew3_params[7],
        
        sos_early = sos_early,
        summer_shoulder = summer_shoulder,
        aut_shoulder = aut_shoulder,
        eos_late = eos_late,
        
        minNDVI_pregreen = minNDVI_pregreen_final,
        maxNDVI = maxNDVI_final,
        
        g15_spl = g15_spl_final,
        g50_spl = g50_spl_final,
        g90_spl = g90_spl_final,
        b10_spl = b10_spl_final,
        b50_spl = b50_spl_final,
        b85_spl = b85_spl_final,
        
        g15_dl = g15_dl,
        g50_dl = g50_dl,
        g90_dl = g90_dl,
        b10_dl = b10_dl,
        b50_dl = b50_dl,
        b85_dl = b85_dl,
        
        auc_dl = auc_dl,
        auc_spl = auc_spl,
        rmse_dl = rmse_dl,
        mad_dl = mad_dl,
        rmse_spl = rmse_spl,
        mad_spl = mad_spl,
        
        m1_op = mnew3_params[1],
        m2_op = mnew3_params[2],
        m3_op = mnew3_params[3],
        m4_op = mnew3_params[4],
        m5_op = mnew3_params[5],
        m6_op = mnew3_params[6],
        m7_op = mnew3_params[7],
        
        g15_spl_op = g15_spl_op,
        g50_spl_op = g50_spl_op,
        g90_spl_op = g90_spl_op,
        b10_spl_op = b10_spl_op,
        b50_spl_op = b50_spl_op,
        b85_spl_op = b85_spl_op,
        
        imCountInit = image_dates_initial,
        imCountFinal = image_dates_final,
        error_message = NA_character_
      )
    },
    
    error = function(e) {
      data.frame(
        UID = uid_batch[i],
        SpCode = if (exists("this_spcode")) this_spcode else NA_character_,
        
        m1 = NA_real_,
        m2 = NA_real_,
        m3 = NA_real_,
        m4 = NA_real_,
        m5 = NA_real_,
        m6 = NA_real_,
        m7 = NA_real_,
        
        sos_early = NA_real_,
        summer_shoulder = NA_real_,
        aut_shoulder = NA_real_,
        eos_late = NA_real_,
        
        minNDVI_pregreen = NA_real_,
        maxNDVI = NA_real_,
        
        g15_spl = NA_real_,
        g50_spl = NA_real_,
        g90_spl = NA_real_,
        b10_spl = NA_real_,
        b50_spl = NA_real_,
        b85_spl = NA_real_,
        
        g15_dl = NA_real_,
        g50_dl = NA_real_,
        g90_dl = NA_real_,
        b10_dl = NA_real_,
        b50_dl = NA_real_,
        b85_dl = NA_real_,
        
        auc_dl = NA_real_,
        auc_spl = NA_real_,
        rmse_dl = NA_real_,
        mad_dl = NA_real_,
        rmse_spl = NA_real_,
        mad_spl = NA_real_,
        
        m1_op = NA_real_,
        m2_op = NA_real_,
        m3_op = NA_real_,
        m4_op = NA_real_,
        m5_op = NA_real_,
        m6_op = NA_real_,
        m7_op = NA_real_,
        
        g15_spl_op = NA_real_,
        g50_spl_op = NA_real_,
        g90_spl_op = NA_real_,
        b10_spl_op = NA_real_,
        b50_spl_op = NA_real_,
        b85_spl_op = NA_real_,
        
        imCountInit = NA_integer_,
        imCountFinal = NA_integer_,
        error_message = as.character(e$message)
      )
    }
  )
  
  return(result_tc)
}

stopCluster(cl)

# 4 Output ------------------------------------------------------------

df_results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

# Optional backfill in case any failed rows still lost SpCode
uid_to_species <- unique(dc_dat[, .(UID, SpCode)])
df_results <- merge(df_results, uid_to_species, by = "UID", all.x = TRUE, suffixes = c("", "_orig"))
df_results[, SpCode := fifelse(is.na(SpCode), SpCode_orig, SpCode)]
df_results[, SpCode_orig := NULL]

batch_tag <- paste0("batch_", sprintf("%03d", batch_number))
output_csv <- file.path(output_dir, paste0(user_output_filename, "_", batch_tag, ".csv"))

fwrite(df_results, output_csv)

cat("Saved batch results to:\n", output_csv, "\n")
cat("Runtime:\n")
print(Sys.time() - start_time)

# Optional combine step ----------------------------------------------

all_batch_files <- list.files(
  path = output_dir,
  pattern = paste0("^", user_output_filename, "_batch_.*\\.csv$"),
  full.names = TRUE
)

all_results <- rbindlist(lapply(all_batch_files, fread), use.names = TRUE, fill = TRUE)

# Backfill species again for combined file, just to be safe
all_results <- merge(all_results, uid_to_species, by = "UID", all.x = TRUE, suffixes = c("", "_orig"))
all_results[, SpCode := fifelse(is.na(SpCode), SpCode_orig, SpCode)]
all_results[, SpCode_orig := NULL]

combined_csv <- file.path(output_dir, paste0(user_output_filename, "_ALL_BATCHES.csv"))
fwrite(all_results, combined_csv)

cat("Saved combined results to:\n", combined_csv, "\n")

# 5 Optional shapefile join ------------------------------------------

if (do_join == "yes") {
  shp <- st_read(shapefile_filename)
  shp_merg <- merge(shp, all_results, by = "UID")
  st_write(shp_merg, shp_merge_filename, delete_dsn = TRUE)
}
