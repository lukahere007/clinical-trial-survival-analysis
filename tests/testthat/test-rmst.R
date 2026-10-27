test_that("RMST matches a hand-computed example without censoring", {
  # Events at 1, 2, 3 in 3 patients: S = 1 on [0,1), 2/3 on [1,2), 1/3 on [2,3), 0 after.
  r <- rmst_km(time = c(1, 2, 3), status = c(1, 1, 1), tau = 3)
  expect_equal(r$rmst, 1 + 2/3 + 1/3)
  # Truncating at 1.5: area = 1 + 0.5 * 2/3.
  expect_equal(rmst_km(c(1, 2, 3), c(1, 1, 1), tau = 1.5)$rmst, 1 + 1/3)
})

test_that("RMST matches a hand-computed example with censoring", {
  # Times 1 (event), 2 (censored), 3 (event), 4 (censored).
  # S: 1 on [0,1); 3/4 on [1,3); 3/4 * 1/2 = 3/8 from 3 on.
  r <- rmst_km(c(1, 2, 3, 4), c(1, 0, 1, 0), tau = 4)
  expect_equal(r$rmst, 1 + 2 * 3/4 + 1 * 3/8)
})

test_that("RMST equals tau when nobody has an event", {
  expect_equal(rmst_km(c(2, 5, 7), c(0, 0, 0), tau = 4)$rmst, 4)
  expect_equal(rmst_km(c(2, 5, 7), c(0, 0, 0), tau = 4)$se, 0)
})

test_that("RMST and SE agree with survival::survfit(rmean = tau)", {
  d <- survival::lung
  d$status01 <- as.integer(d$status == 2)
  tau <- 500
  ours <- rmst_km(d$time, d$status01, tau)
  ref <- summary(survival::survfit(survival::Surv(time, status01) ~ 1, data = d),
                 rmean = tau)$table
  expect_equal(ours$rmst, unname(ref["rmean"]), tolerance = 1e-8)
  expect_equal(ours$se, unname(ref["se(rmean)"]), tolerance = 1e-6)
})

test_that("RMST difference is comparator minus reference", {
  time <- c(1, 2, 3, 2, 4, 6)
  status <- c(1, 1, 1, 1, 1, 1)
  grp <- rep(c("A", "B"), each = 3)
  out <- rmst_diff(time, status, grp, tau = 3, ref = "A", comparator = "B")
  expect_equal(out$diff, out$rmst_comp - out$rmst_ref)
  expect_gt(out$diff, 0)
  expect_lt(out$lower, out$diff)
  expect_gt(out$upper, out$diff)
})
