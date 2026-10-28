# Cleaning: parse censored results, drop duplicate rows, reshape long.
# Every rule below is counted in the QC log so the denominator is auditable.

ANALYTES <- c("HFPODA", "PFBS", "PFHxA", "PFHxS", "PFNA", "PFOA", "PFOS")

# Michigan MCLs (R 325.10604g, effective 2020-08-03), ppt.
MI_MCL <- c(HFPODA = 370, PFBS = 420, PFHxA = 400000, PFHxS = 51, PFNA = 6, PFOA = 8, PFOS = 16)
# Federal MCLs (April 2024 rule) for PFOA/PFOS, kept in EPA's 2025 proposal.
US_MCL <- c(PFOA = 4, PFOS = 4)

# "<2" -> value 2, censored TRUE; "5.1" -> 5.1, FALSE; "" / NA -> NA.
parse_result <- function(x) {
  x <- trimws(as.character(x))
  cens <- grepl("^<", x)
  val <- suppressWarnings(as.numeric(sub("^<\\s*", "", x)))
  val[is.na(x) | x == ""] <- NA
  data.frame(value = val, censored = cens & !is.na(val))
}

clean_results <- function(d) {
  log <- data.frame(step = "raw rows", rows = nrow(d))
  step <- function(name) log <<- rbind(log, data.frame(step = name, rows = nrow(d)))

  d$SampleDate <- as.Date(d$SampleDate)
  d$LocName <- trimws(d$LocName)
  d$LocName[d$LocName == ""] <- NA
  res_key <- do.call(paste, c(d[c("WSSN", "SampleDate")], d[ANALYTES], sep = "|"))

  # 1. Same sample loaded twice (EGLE lab rows mirrored as "Lansing Lab" rows
  #    with no sampling point): drop the unlabeled copy when a labeled row has
  #    identical results for the same system and date.
  labeled <- res_key[!is.na(d$LocName)]
  d <- d[!(is.na(d$LocName) & res_key %in% labeled), ]
  step("drop unlabeled mirror rows")

  # 2. Exact duplicates of point + date + all results.
  d <- d[!duplicated(do.call(paste, d[c("WSSN", "LocName", "SampleDate", ANALYTES)])), ]
  step("drop exact duplicates")

  # 3. Rows still without a sampling point can't enter a per-point series.
  d <- d[!is.na(d$LocName), ]
  step("drop rows with no sampling point")

  d$point <- paste(d$WSSN, d$LocName, sep = ":")
  long <- do.call(rbind, lapply(ANALYTES, function(a) {
    p <- parse_result(d[[a]])
    data.frame(WSSN = d$WSSN, point = d$point, system = d$SystemName, HexID = d$HexID,
               date = d$SampleDate, lab = d$LabNameCode, analyte = a, p)
  }))
  long <- long[!is.na(long$value), ]

  # 3b. Detects < 0.1 ppt sit below any plausible EPA 537.1 detection limit;
  #     they look like ug/L entered as ng/L (x1000 gives 1-1.6 ppt). Excluded,
  #     not converted: the unit is a guess.
  bad <- !long$censored & long$value < 0.1
  log <- rbind(log, data.frame(step = "results < 0.1 ppt excluded (suspected unit error)", rows = sum(bad)))
  long <- long[!bad, ]

  # 4. Remaining same point/date/analyte repeats (field duplicates or
  #    confirmation samples) are averaged, as R 325.10717d(12) does for
  #    confirmation samples. If any repeat is a detect, the average is over
  #    detects only (a "<2" is not a measured 2); all-censored stays "<RL".
  key <- paste(long$point, long$date, long$analyte)
  rep_keys <- unique(key[duplicated(key)])
  mixed <- 0
  out <- long[!duplicated(key), ]
  idx <- match(rep_keys, key[!duplicated(key)])
  for (i in seq_along(rep_keys)) {
    g <- long[key == rep_keys[i], ]
    det <- !g$censored
    if (any(det) && !all(det)) mixed <- mixed + 1
    out$value[idx[i]] <- if (any(det)) mean(g$value[det]) else min(g$value)
    out$censored[idx[i]] <- !any(det)
  }
  log <- rbind(log, data.frame(step = c("point-date-analyte groups collapsed", "  of which mixed detect/nondetect"),
                               rows = c(length(rep_keys), mixed)))
  list(data = out[order(out$point, out$analyte, out$date), ], log = log)
}
