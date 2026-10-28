# Full pipeline: Rscript R/run.R  (from repo root; needs data/raw from R/ingest.R)
source("R/clean.R")
source("R/stats.R")
suppressPackageStartupMessages({ library(ggplot2); library(sf) })
dir.create("outputs", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
set.seed(20261028)

r <- clean_results(read.csv("data/raw/pws_results.csv", stringsAsFactors = FALSE))
x <- r$data
saveRDS(x, "data/clean.rds")
write.csv(r$log, "outputs/qc_log.csv", row.names = FALSE)

# A
occ <- occurrence(x)
write.csv(occ, "outputs/A_occurrence.csv", row.names = FALSE)

# B
mi <- raa_screen(x, MI_MCL)
us <- raa_screen(x, US_MCL)
summ_b <- function(s, label) do.call(rbind, lapply(split(s, s$analyte), function(a) data.frame(
  standard = label, analyte = a$analyte[1], mcl = a$mcl[1], points = nrow(a),
  any_result_over = sum(a$any_result_over),
  raa_over_nd0 = sum(a$raa_nd0 > a$mcl), raa_over_ndRL = sum(a$raa_ndRL > a$mcl),
  raa4_over_nd0 = sum(a$raa4_nd0 > a$mcl, na.rm = TRUE), raa4_over_ndRL = sum(a$raa4_ndRL > a$mcl, na.rm = TRUE),
  points_with_full_year = sum(is.finite(a$raa4_nd0)),
  latest_raa4_over_nd0 = sum(a$raa4_latest_nd0 > a$mcl, na.rm = TRUE))))
b <- rbind(summ_b(mi, "Michigan"), summ_b(us, "Federal 2024"))
b$flips_nd_convention <- b$raa_over_ndRL - b$raa_over_nd0
write.csv(b, "outputs/B_raa_screen.csv", row.names = FALSE)
write.csv(mi, "outputs/B_raa_points_michigan.csv", row.names = FALSE)
write.csv(us, "outputs/B_raa_points_federal.csv", row.names = FALSE)

# C
tr <- trend_tests(x)
write.csv(tr, "outputs/C_trends.csv", row.names = FALSE)
sig <- tr[tr$q < 0.10, ]
up <- sum(sig$slope_ppt_yr > 0); down <- sum(sig$slope_ppt_yr < 0)
sign_p <- if (up + down > 0) binom.test(down, up + down)$p.value else NA
c_summary <- data.frame(tests = nrow(tr), fdr_significant = nrow(sig), increasing = up, decreasing = down,
                        sign_test_p = sign_p, points = length(unique(tr$point)))
write.csv(c_summary, "outputs/C_trend_summary.csv", row.names = FALSE)
write.csv(as.data.frame(table(analyte = sig$analyte, direction = ifelse(sig$slope_ppt_yr > 0, "up", "down"))),
          "outputs/C_trend_by_analyte.csv", row.names = FALSE)

# D
ss <- ssi_screen(x)
write.csv(ss, "outputs/D_ssi.csv", row.names = FALSE)
write.csv(as.data.frame(table(analyte = ss$analyte, rule = ss$rule, status = ss$status)),
          "outputs/D_ssi_summary.csv", row.names = FALSE)

# E
sim_grid <- function(sdlog) do.call(rbind, lapply(c(8, 20, 50), function(n)
  do.call(rbind, lapply(c(0.5, 0.7, 0.85), function(p) sim_estimators(n, p, sdlog = sdlog)))))
sim <- sim_grid(1)                                    # pre-registered
write.csv(sim, "outputs/E_simulation.csv", row.names = FALSE)
sim_x <- rbind(sim_grid(0.5), sim_grid(2))            # exploratory: is sdlog = 1 a lucky case?
write.csv(sim_x, "outputs/E_simulation_sdlog_exploratory.csv", row.names = FALSE)

# F: maps by EGLE hexbin
hex <- st_read("data/raw/hexbins.geojson", quiet = TRUE)
pt <- unique(x[, c("point", "HexID")])
pt$detected <- pt$point %in% x$point[!x$censored]
pt$over_mi <- pt$point %in% mi$point[mi$raa4_nd0 > mi$mcl]
pt$over_us <- pt$point %in% us$point[us$raa4_nd0 > us$mcl]
hx <- aggregate(cbind(points = 1, detected, over_mi, over_us) ~ HexID, data = pt, FUN = sum)
hx$pct_detected <- 100 * hx$detected / hx$points
hex <- merge(hex, hx, by = "HexID")
write.csv(st_drop_geometry(hex), "outputs/F_hexbins.csv", row.names = FALSE)
theme_set(theme_minimal(base_size = 11))
ggsave("figures/map_detection.png", width = 6, height = 7, dpi = 160,
  ggplot(hex) + geom_sf(aes(fill = pct_detected), colour = NA) +
    scale_fill_viridis_c("% points with\nany detect", option = "magma", direction = -1) +
    labs(title = "Michigan public water supplies: PFAS detection",
         subtitle = "Share of sampling points with >= 1 detect (any of 7 PFAS), 2020-2025",
         caption = "Data: EGLE PWS PFAS open data; EGLE hexbins"))
ggsave("figures/map_exceed.png", width = 6, height = 7, dpi = 160,
  ggplot(hex) + geom_sf(fill = "grey90", colour = NA) +
    geom_sf(data = transform(hex[hex$over_us > 0, ], status = ifelse(over_mi > 0, "Over a Michigan MCL", "Over federal 4 ppt only")),
            aes(fill = status), colour = NA) +
    scale_fill_manual(NULL, values = c("Over a Michigan MCL" = "#b2182b", "Over federal 4 ppt only" = "#f4a582")) +
    theme(legend.position = "bottom") +
    labs(title = "Screening full-year RAA over a PFAS MCL, 2020-2025",
         subtitle = "Hexbins with >= 1 sampling point over (nondetects = 0, >= 4 samples in window)",
         caption = "Screening estimate, not an EGLE compliance determination"))

# Figures: trend example + simulation
series_plot <- function(pt, an, title, subtitle) {
  ex <- x[x$point == pt & x$analyte == an, ]
  ggplot(ex, aes(date, value, shape = censored)) + geom_point(size = 2) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("detect", "<2 (nondetect)"), name = NULL) +
    geom_hline(yintercept = MI_MCL[[an]], linetype = 2) +
    geom_hline(yintercept = US_MCL[[an]], linetype = 3) +
    expand_limits(y = 0) +
    labs(x = NULL, y = paste(an, "(ng/L)"), title = title, subtitle = subtitle)
}
pf <- sig[sig$analyte %in% c("PFOA", "PFOS"), ]
for (dir in c("up", "down")) {
  cand <- pf[if (dir == "up") pf$slope_ppt_yr > 0 else pf$slope_ppt_yr < 0, ]
  best <- cand[order(cand$q), ][1, ]
  ggsave(sprintf("figures/trend_%s.png", dir), width = 7, height = 3.4, dpi = 160,
    series_plot(best$point, best$analyte,
      sprintf("Most significant %s %s trend (one sampling point)", if (dir == "up") "increasing" else "decreasing", best$analyte),
      sprintf("ATS slope %+.2f ppt/yr, BH q = %.2g. Dashed: Michigan MCL; dotted: federal 4 ppt",
              best$slope_ppt_yr, best$q)))
}
ex_ssi <- ss[ss$status == "SSI" & ss$analyte %in% c("PFOA", "PFOS") & ss$rule != "double quantification", ]
ex_ssi$peak <- mapply(function(p, a) max(x$value[x$point == p & x$analyte == a]), ex_ssi$point, ex_ssi$analyte)
ex_ssi <- ex_ssi[order(-ex_ssi$peak), ][1, ]
ggsave("figures/ssi_example.png", width = 7, height = 3.4, dpi = 160,
  series_plot(ex_ssi$point, ex_ssi$analyte,
    sprintf("Verified SSI: %s, intrawell prediction limit", ex_ssi$analyte),
    sprintf("Grey line ends the 8-date baseline; result %d exceeds its max, confirmed by result %d",
            ex_ssi$ssi_index, ex_ssi$ssi_index + 1)) +
    geom_vline(xintercept = sort(x$date[x$point == ex_ssi$point & x$analyte == ex_ssi$analyte])[8] + 0.5, colour = "grey50"))
ggsave("figures/simulation_rmse.png", width = 7.5, height = 5.5, dpi = 160,
  ggplot(transform(rbind(sim, sim_x), n = factor(paste("n =", n), paste("n =", c(8, 20, 50)))), aes(factor(100 * cens_rate), rmse_pct, colour = estimator, group = estimator)) +
    geom_line() + geom_point() + facet_grid(paste("sdlog =", sdlog) ~ n, scales = "free_y") +
    labs(x = "% censored", y = "RMSE of mean (% of truth)", title = "Estimating a mean under a single reporting limit",
         subtitle = "Lognormal, 2,000 common replicates per cell. sdlog = 1 pre-registered; 0.5, 2 exploratory"))

cat("done\n")
print(occ); print(b); print(c_summary); print(table(ss$rule, ss$status))
