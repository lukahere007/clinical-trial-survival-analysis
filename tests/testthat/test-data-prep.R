dat <- make_analysis_data(survival::colon)

test_that("one row per randomised patient", {
  expect_equal(nrow(dat), length(unique(survival::colon$id)))
  expect_equal(nrow(dat), 929L)
  expect_false(anyDuplicated(dat$id) > 0)
  expect_equal(as.vector(table(dat$arm)), c(315L, 310L, 304L))
})

test_that("long format has one row per patient per endpoint", {
  long <- endpoints_long(dat)
  expect_equal(nrow(long), 3L * nrow(dat))
  expect_true(all(table(long$id, long$endpoint) == 1L))
})

test_that("endpoint definitions are internally consistent", {
  # RFS cannot end after death / last follow-up.
  expect_true(all(dat$rfs_time <= dat$os_time + 1e-12))
  # Every death and every recurrence is an RFS event.
  expect_true(all(dat$rfs_event[dat$os_event == 1] == 1))
  expect_true(all(dat$rfs_event[dat$ttr_event == 1] == 1))
  # RFS events = recurrences + deaths without recurrence.
  expect_equal(sum(dat$rfs_event), sum(dat$ttr_event) +
                 sum(dat$cr_status == "Death without recurrence"))
  # Arm and outcomes are complete.
  expect_false(anyNA(dat[, c("arm", "os_time", "os_event", "rfs_time", "rfs_event")]))
})

test_that("missing-indicator variable keeps every patient", {
  expect_false(anyNA(dat$differ_mi))
  expect_equal(sum(dat$differ_mi == "Missing"), sum(is.na(dat$differ)))
})

test_that("bad input is rejected", {
  expect_error(make_analysis_data(survival::colon[, -1]), "Missing columns")
  expect_error(make_analysis_data(survival::colon[-1, ]), "both a recurrence and a death")
})
