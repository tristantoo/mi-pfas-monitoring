source("../../R/clean.R")
source("../../R/stats.R")

test_that("censored results parse", {
  p <- parse_result(c("<2", "5.1", "", NA, " <2 "))
  expect_equal(p$value, c(2, 5.1, NA, NA, 2))
  expect_equal(p$censored, c(TRUE, FALSE, FALSE, FALSE, TRUE))
})

toy <- function() {
  base <- data.frame(WSSN = 1, SystemName = "A", LocName = "EP01", SysSampleCode = "s1",
                     SampleDate = "2022-01-01", LabNameCode = "EGLE", LabNumber = "1", HexID = 9,
                     HFPODA = "<2", PFBS = "<2", PFHxA = "<2", PFHxS = "<2", PFNA = "<2", PFOA = "5", PFOS = "<2")
  mirror <- transform(base, LocName = NA, SysSampleCode = "s1b", LabNameCode = "Lansing Lab")
  rep1 <- transform(base, SampleDate = "2022-04-01", SysSampleCode = "s2", PFOA = "4")
  rep2 <- transform(rep1, SysSampleCode = "s3", PFOA = "<2")   # same point/date, mixed
  rbind(base, mirror, rep1, rep2)
}

test_that("cleaning drops mirrors and averages repeats over detects only", {
  r <- clean_results(toy())
  x <- r$data[r$data$analyte == "PFOA", ]
  expect_equal(nrow(x), 2)                     # mirror dropped, repeat collapsed
  expect_equal(x$value, c(5, 4))               # "<2" not averaged in as 2
  expect_false(any(x$censored))
  expect_equal(r$log$rows[r$log$step == "  of which mixed detect/nondetect"], 1)
})

test_that("running annual average uses a trailing 365-day window", {
  d <- as.Date(c("2022-01-01", "2022-06-01", "2022-12-01", "2023-06-01"))
  v <- c(10, 2, 6, 4); cen <- c(FALSE, TRUE, FALSE, FALSE)
  expect_equal(raa(d, v, cen, nd_value = 0), c(10, 5, 16 / 3, 5))  # Jun-2023 window: Dec + Jun only
  expect_equal(raa(d, v, cen, nd_value = 2)[2], 6)
  expect_true(is.na(raa(d, v, cen, 0, min_n = 4)[3]))
})

test_that("verification resampling", {
  nd8 <- rep(TRUE, 8)
  v <- c(rep(2, 8), 3, 3.5)
  expect_equal(ssi_check(v, c(nd8, FALSE, FALSE))$status, "SSI")        # two consecutive detects
  expect_equal(ssi_check(v, c(nd8, FALSE, TRUE))$status, "none")         # disconfirmed
  expect_equal(ssi_check(v[1:9], c(nd8, FALSE))$status, "pending")       # no resample yet
  b <- c(3, 4, 5, 2, 2, 6, 2, 2); bc <- c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE)
  s <- ssi_check(c(b, 7, 8), c(bc, FALSE, FALSE))
  expect_equal(s$rule, "NP prediction limit (baseline max), 1-of-2")
  expect_equal(s$status, "SSI")
  expect_equal(ssi_check(c(b, 5.5, 7), c(bc, FALSE, FALSE))$status, "pending")  # 5.5 < PL 6; 7 unverified
  expect_equal(ssi_check(c(b, 7, 5.5), c(bc, FALSE, FALSE))$status, "none")     # 7 disconfirmed by 5.5
})

test_that("ATS recovers a known increasing trend with censoring", {
  set.seed(1)
  yrs <- seq(2020, 2025, by = 0.25)
  y <- 2 + 1.5 * (yrs - 2020) + rnorm(length(yrs), 0, 0.5)
  cen <- y < 2.5; y[cen] <- 2.5
  r <- ATSmini(y, cen, yrs)
  expect_gt(r$slope, 1); expect_lt(r$slope, 2); expect_lt(r$pval, 0.01)
})
