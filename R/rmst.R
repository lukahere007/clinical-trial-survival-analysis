# Restricted mean survival time (RMST) from the Kaplan-Meier estimator.
#
# RMST(tau) = integral_0^tau S(t) dt, the mean event-free time over the
# first tau time units. It needs no proportional-hazards assumption.

#' Kaplan-Meier based RMST with a Greenwood-type standard error.
#'
#' Variance: sum over event times t_j <= tau of A_j^2 * d_j / (n_j (n_j - d_j)),
#' where A_j is the area under the KM curve between t_j and tau.
#'
#' @param time follow-up times.
#' @param status 1 = event, 0 = censored.
#' @param tau truncation time.
#' @return list with rmst, se, tau.
rmst_km <- function(time, status, tau) {
  stopifnot(length(time) == length(status), all(status %in% c(0, 1)),
            is.numeric(tau), length(tau) == 1, tau > 0)
  if (any(time < 0)) stop("Negative times are not allowed.")

  ev_times <- sort(unique(time[status == 1 & time <= tau]))
  n_risk <- vapply(ev_times, function(t) sum(time >= t), numeric(1))
  n_ev   <- vapply(ev_times, function(t) sum(time == t & status == 1), numeric(1))
  surv   <- cumprod(1 - n_ev / n_risk)

  # Step function: S = 1 on [0, t_1), surv[j] on [t_j, t_{j+1}), up to tau.
  knots   <- c(0, ev_times, tau)
  heights <- c(1, surv)
  widths  <- diff(knots)
  rmst <- sum(heights * widths)

  # Area from each event time to tau.
  piece <- heights[-1] * widths[-1]
  area_after <- rev(cumsum(rev(piece)))
  denom <- n_risk * (n_risk - n_ev)
  terms <- ifelse(denom > 0, area_after^2 * n_ev / denom, 0)
  se <- sqrt(sum(terms))

  list(rmst = rmst, se = se, tau = tau)
}

#' RMST difference between two arms (arm2 - arm1) with 95% CI and z-test.
rmst_diff <- function(time, status, group, tau, ref, comparator, conf = 0.95) {
  g1 <- group == ref
  g2 <- group == comparator
  if (!any(g1) || !any(g2)) stop("Both groups must be present.")
  r1 <- rmst_km(time[g1], status[g1], tau)
  r2 <- rmst_km(time[g2], status[g2], tau)
  diff <- r2$rmst - r1$rmst
  se <- sqrt(r1$se^2 + r2$se^2)
  z <- qnorm(1 - (1 - conf) / 2)
  data.frame(
    tau = tau, comparison = paste(comparator, "vs", ref),
    rmst_ref = r1$rmst, rmst_comp = r2$rmst,
    diff = diff, lower = diff - z * se, upper = diff + z * se,
    p = 2 * pnorm(-abs(diff / se))
  )
}
