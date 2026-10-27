# Data preparation for the colon adjuvant chemotherapy trial (survival::colon).
#
# survival::colon has two rows per patient: etype = 1 (recurrence) and
# etype = 2 (death). make_analysis_data() turns this into one row per
# patient with one time/event pair per endpoint, in years.

DAYS_PER_YEAR <- 365.25

#' One row per patient, with OS, RFS and TTR endpoints and labelled covariates.
#'
#' @param colon_raw data frame in the layout of survival::colon.
#' @return data frame with one row per randomised patient.
make_analysis_data <- function(colon_raw) {
  required <- c("id", "rx", "etype", "time", "status", "sex", "age", "obstruct",
                "perfor", "adhere", "nodes", "node4", "differ", "extent", "surg")
  missing_cols <- setdiff(required, names(colon_raw))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }

  rec <- colon_raw[colon_raw$etype == 1, ]
  dth <- colon_raw[colon_raw$etype == 2, ]
  if (anyDuplicated(rec$id) || anyDuplicated(dth$id)) {
    stop("Expected at most one recurrence and one death record per patient.")
  }
  if (!setequal(rec$id, dth$id)) {
    stop("Every patient needs both a recurrence and a death record.")
  }
  rec <- rec[order(rec$id), ]
  dth <- dth[order(dth$id), ]

  os_time  <- dth$time
  os_event <- as.integer(dth$status == 1)
  rec_event <- as.integer(rec$status == 1)

  # Recurrence-free survival: first of recurrence or death.
  rfs_event <- as.integer(rec_event == 1 | os_event == 1)
  rfs_time  <- ifelse(rec_event == 1, rec$time, os_time)

  # Competing-risks status for time to recurrence.
  cr_status <- ifelse(rec_event == 1, "Recurrence",
                      ifelse(os_event == 1, "Death without recurrence", "Censored"))

  out <- data.frame(
    id        = dth$id,
    arm       = factor(dth$rx, levels = c("Obs", "Lev", "Lev+5FU")),
    os_time   = os_time / DAYS_PER_YEAR,
    os_event  = os_event,
    rfs_time  = rfs_time / DAYS_PER_YEAR,
    rfs_event = rfs_event,
    # TTR: recurrence is the event; death without recurrence censors (cause-specific).
    ttr_time  = rfs_time / DAYS_PER_YEAR,
    ttr_event = rec_event,
    cr_status = factor(cr_status,
                       levels = c("Censored", "Recurrence", "Death without recurrence")),
    age       = dth$age,
    sex       = factor(dth$sex, levels = c(0, 1), labels = c("Female", "Male")),
    nodes     = dth$nodes,
    node4     = factor(dth$node4, levels = c(0, 1), labels = c("1-4", ">4")),
    extent    = factor(cut(dth$extent, breaks = c(0, 2, 3, 4),
                           labels = c("Submucosa/muscle", "Serosa", "Contiguous structures"))),
    differ    = factor(dth$differ, levels = 1:3, labels = c("Well", "Moderate", "Poor")),
    obstruct  = factor(dth$obstruct, levels = c(0, 1), labels = c("No", "Yes")),
    perfor    = factor(dth$perfor, levels = c(0, 1), labels = c("No", "Yes")),
    adhere    = factor(dth$adhere, levels = c(0, 1), labels = c("No", "Yes")),
    surg      = factor(dth$surg, levels = c(0, 1), labels = c("Short", "Long")),
    stringsAsFactors = FALSE
  )

  # Sensitivity analysis S2: explicit missing-indicator level.
  out$differ_mi <- addNA(out$differ)
  levels(out$differ_mi)[is.na(levels(out$differ_mi))] <- "Missing"

  rownames(out) <- NULL
  out
}

#' Stack the per-patient endpoints into long format: one row per patient per endpoint.
endpoints_long <- function(dat) {
  eps <- c(OS = "os", RFS = "rfs", TTR = "ttr")
  do.call(rbind, lapply(names(eps), function(nm) {
    p <- eps[[nm]]
    data.frame(id = dat$id, arm = dat$arm, endpoint = nm,
               time = dat[[paste0(p, "_time")]],
               event = dat[[paste0(p, "_event")]])
  }))
}
