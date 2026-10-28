# Analysis functions (see PREREG.md for the protocol each implements).

suppressPackageStartupMessages({
  library(EnvStats)
  library(NADA2)
})

by_series <- function(x) split(x, paste(x$point, x$analyte, sep = "|"), drop = TRUE)

# A. Occurrence: one value per point (its max detect, or "<RL" if never detected).
occurrence <- function(x, probs = c(0.90, 0.95, 0.99)) {
  do.call(rbind, lapply(split(x, x$analyte), function(a) {
    pm <- tapply(ifelse(a$censored, 0, a$value), a$point, max)  # 0 = never detected
    q <- quantile(pm, probs, type = 1, names = FALSE)
    data.frame(analyte = a$analyte[1], points = length(pm),
               points_detected = sum(pm > 0), pct_points_detected = 100 * mean(pm > 0),
               pct_results_detected = 100 * mean(!a$censored),
               censored_pct = 100 * mean(a$censored),
               p90 = q[1], p95 = q[2], p99 = q[3])  # 0 means "<RL"
  }))
}

# B. Trailing-365-day running annual average at each sample date.
raa <- function(date, value, censored, nd_value, min_n = 1) {
  v <- ifelse(censored, nd_value, value)
  sapply(seq_along(date), function(i) {
    w <- date > date[i] - 365 & date <= date[i]
    if (sum(w) < min_n) NA else mean(v[w])
  })
}

raa_screen <- function(x, mcl) {
  do.call(rbind, lapply(by_series(x[x$analyte %in% names(mcl), ]), function(s) {
    s <- s[order(s$date), ]
    lim <- mcl[[s$analyte[1]]]
    mx <- function(...) suppressWarnings(max(raa(s$date, s$value, s$censored, ...), na.rm = TRUE))
    full <- raa(s$date, s$value, s$censored, 0, min_n = 4)
    data.frame(point = s$point[1], HexID = s$HexID[1], analyte = s$analyte[1], mcl = lim,
               n = nrow(s),
               any_result_over = any(!s$censored & s$value > lim),
               raa_nd0 = mx(0), raa_ndRL = mx(RL),
               raa4_nd0 = mx(0, min_n = 4), raa4_ndRL = mx(RL, min_n = 4),
               # exploratory (not pre-registered): latest full-year RAA = current status
               raa4_latest_nd0 = if (any(!is.na(full))) full[max(which(!is.na(full)))] else NA,
               latest_date = max(s$date))
  }))
}
RL <- 2  # every nondetect in this dataset is "<2"

# C. Akritas-Theil-Sen trend (censored Kendall), slope in ppt/yr.
trend_tests <- function(x, min_dates = 8, min_span_days = 730, min_detects = 3) {
  out <- lapply(by_series(x), function(s) {
    if (nrow(s) < min_dates || diff(range(s$date)) < min_span_days || sum(!s$censored) < min_detects) return(NULL)
    yrs <- as.numeric(s$date) / 365.25
    r <- tryCatch(ATSmini(s$value, s$censored, yrs), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    data.frame(point = s$point[1], analyte = s$analyte[1], n = nrow(s), detects = sum(!s$censored),
               slope_ppt_yr = r$slope, tau = r$tau, p = r$pval)
  })
  out <- do.call(rbind, out)
  out$q <- p.adjust(out$p, "BH")
  out
}

# D. Intrawell detection monitoring with verification resampling.
#    Returns status: "SSI" (verified), "pending" (unverified on last sample), "none".
ssi_check <- function(value, censored, n_base = 8) {
  n <- length(value)
  after <- seq_len(n)[-seq_len(n_base)]
  if (all(censored[seq_len(n_base)])) {
    rule <- "double quantification"
    hit <- !censored[after]
  } else {
    rule <- "NP prediction limit (baseline max), 1-of-2"
    pl <- max(value[seq_len(n_base)][!censored[seq_len(n_base)]])
    hit <- !censored[after] & value[after] > pl
  }
  verified <- which(hit[-length(hit)] & hit[-1])[1]   # exceedance confirmed by the next result
  status <- if (!is.na(verified)) "SSI"
    else if (hit[length(hit)]) "pending"               # last result exceeds, no resample yet
    else "none"                                        # no exceedance, or every one disconfirmed
  data.frame(rule = rule, status = status, ssi_index = n_base + verified)
}

ssi_screen <- function(x, n_base = 8) {
  do.call(rbind, lapply(by_series(x), function(s) {
    if (nrow(s) <= n_base) return(NULL)
    s <- s[order(s$date), ]
    cbind(data.frame(point = s$point[1], analyte = s$analyte[1], n = nrow(s)),
          ssi_check(s$value, s$censored, n_base))
  }))
}

# E. Simulation: mean estimators under single-RL left censoring.
sim_estimators <- function(n, cens_rate, reps = 2000, sdlog = 1) {
  rl <- exp(qnorm(cens_rate) * sdlog)          # meanlog 0 -> RL at the cens_rate quantile
  truth <- exp(sdlog^2 / 2)
  est <- t(replicate(reps, {
    y <- rlnorm(n, 0, sdlog)
    cen <- y < rl
    obs <- ifelse(cen, rl, y)
    km <- if (sum(!cen) >= 2) tryCatch(enparCensored(obs, cen)$parameters[["mean"]], error = function(e) NA) else NA
    ros <- if (sum(!cen) >= 2) tryCatch(elnormAltCensored(obs, cen, method = "rROS")$parameters[["mean"]], error = function(e) NA) else NA
    c(sub0 = mean(ifelse(cen, 0, y)), subhalf = mean(ifelse(cen, rl / 2, y)),
      subRL = mean(obs), KM = km, ROS = ros)
  }))
  # Score every estimator on the same replicates: KM/ROS need >= 2 detects,
  # and keeping only their successful runs would select high-valued samples.
  ok <- stats::complete.cases(est)
  do.call(rbind, lapply(colnames(est), function(m) {
    e <- est[ok, m]
    data.frame(n = n, cens_rate = cens_rate, sdlog = sdlog, estimator = m,
               failed_pct = 100 * mean(!ok),
               bias_pct = 100 * (mean(e) - truth) / truth,
               rmse_pct = 100 * sqrt(mean((e - truth)^2)) / truth)
  }))
}
