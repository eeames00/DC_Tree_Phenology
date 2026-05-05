library(data.table)

# Install if needed
if (!requireNamespace("callr", quietly = TRUE)) install.packages("callr")
if (!requireNamespace("phenofit", quietly = TRUE)) install.packages("phenofit")

out_dir  <- "C:/Users/valea/OneDrive/Desktop/pheno_sensitivity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

qc_csv  <- file.path(out_dir, "qc_results_sample_all_years.csv")
met_csv <- file.path(out_dir, "qc_metrics_beck_sample_all_years.csv")

# start fresh
if (file.exists(qc_csv))  file.remove(qc_csv)
if (file.exists(met_csv)) file.remove(met_csv)

# dt_sample must exist in your session already
stopifnot(exists("dt_sample"))
setDT(dt_sample)
stopifnot(all(c("year","UID","t","y") %in% names(dt_sample)))

# minimal cleaning
dt_sample <- dt_sample[!is.na(t) & is.finite(y)]
dt_sample[, y := as.numeric(y)]
dt_sample[, t := as.Date(t)]
setkey(dt_sample, year, UID, t)

years <- sort(unique(dt_sample$year))

# ---------- helper: run ONE UID-year in a separate R process ----------
run_one_uid_year <- function(d_sub, uid, yy, min_n = 25L, min_amp = 0.08) {
  
  # This function runs in a *child* R session via callr::r()
  # so it must load its own packages.
  callr::r(function(d_sub, uid, yy, min_n, min_amp) {
    library(data.table)
    library(phenofit)
    
    setDT(d_sub)
    
    # collapse dup dates
    dd <- d_sub[!is.na(t) & is.finite(y), .(y = mean(y)), by = t][order(t)]
    if (nrow(dd) < min_n) {
      return(list(ok=FALSE, reason="too_few_points", out=NULL,
                  qc=list(n_points=nrow(dd), amp=NA_real_)))
    }
    
    amp <- diff(range(dd$y, na.rm = TRUE))
    if (!is.finite(amp) || amp < min_amp) {
      return(list(ok=FALSE, reason="low_amplitude", out=NULL,
                  qc=list(n_points=nrow(dd), amp=amp)))
    }
    
    # phenofit
    INPUT <- phenofit::check_input(dd$t, dd$y)
    brks  <- phenofit::season(INPUT)
    fits  <- phenofit::curvefits(INPUT, brks,
                                 options = list(methods="Beck", wFUN="wTSM", verbose=FALSE))
    PHE   <- phenofit::get_pheno(fits)
    
    out <- as.data.table(PHE$doy$Beck)
    out[, `:=`(UID=uid, year=yy, method="Beck")]
    
    list(ok=TRUE, reason=NA_character_, out=out,
         qc=list(n_points=nrow(dd), amp=amp))
  },
  args = list(d_sub=d_sub, uid=uid, yy=yy, min_n=min_n, min_amp=min_amp),
  # If a child process crashes, callr will throw an error here,
  # which we handle in the parent loop.
  show = FALSE)
}

# ---------- main loop ----------
wrote_qc  <- FALSE
wrote_met <- FALSE

for (yy in years) {
  message("YEAR ", yy)
  
  uids <- unique(dt_sample[year == yy, UID])
  uids <- sort(uids)
  
  for (i in seq_along(uids)) {
    uid <- uids[i]
    d_sub <- dt_sample[year == yy & UID == uid, .(t, y)]
    
    fit <- tryCatch(
      run_one_uid_year(d_sub, uid, yy),
      error = function(e) {
        # includes child crash info
        list(ok=FALSE,
             reason=paste0("R_ABORT_or_child_error: ", conditionMessage(e)),
             out=NULL,
             qc=list(n_points=nrow(d_sub), amp=NA_real_))
      }
    )
    
    # QC row
    qc_row <- data.table(
      UID      = uid,
      year     = yy,
      n_points = fit$qc$n_points,
      amp      = fit$qc$amp,
      success  = isTRUE(fit$ok),
      reason   = if (isTRUE(fit$ok)) NA_character_ else fit$reason
    )
    fwrite(qc_row, qc_csv, append = wrote_qc, col.names = !wrote_qc)
    wrote_qc <- TRUE
    
    # Metrics row(s)
    if (isTRUE(fit$ok) && !is.null(fit$out) && nrow(fit$out) > 0) {
      fwrite(fit$out, met_csv, append = wrote_met, col.names = !wrote_met)
      wrote_met <- TRUE
    }
    
    if (i %% 25 == 0) { gc(); message("  done ", i, " / ", length(uids)) }
  }
}

# ---------- summaries ----------
qc <- fread(qc_csv)

qc_by_year <- qc[, .(
  trees = .N,
  success_n = sum(success),
  success_rate = mean(success)
), by = year][order(year)]
print(qc_by_year)

fail_reasons <- qc[success == FALSE, .N, by = .(year, reason)][order(year, -N)][1:40]
print(fail_reasons)

if (file.exists(met_csv)) {
  met <- fread(met_csv)
  print(met[, .N, by = year][order(year)])
}