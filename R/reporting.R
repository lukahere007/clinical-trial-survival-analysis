# Tables and figures used by analysis/colon_trial_report.Rmd.

ARM_COLOURS <- c("Obs" = "#2a78d6", "Lev" = "#eb6834", "Lev+5FU" = "#1baf7a")
ARM_LINETYPES <- c("Obs" = "solid", "Lev" = "dashed", "Lev+5FU" = "solid")

theme_report <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e1e0d9", linewidth = 0.3),
      axis.line = ggplot2::element_line(colour = "#c3c2b7", linewidth = 0.4),
      legend.position = "top",
      legend.justification = "left",
      plot.title.position = "plot",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
}

fmt_p_eq <- function(p) {
  ifelse(p < 0.001, "p < 0.001", paste("p =", formatC(p, format = "f", digits = 3)))
}

fmt_ci <- function(est, lo, hi, digits = 2) {
  f <- function(x) formatC(x, format = "f", digits = digits)
  paste0(f(est), " (", f(lo), " to ", f(hi), ")")
}

#' Table 1: baseline characteristics by arm.
baseline_table <- function(dat, arm = "arm") {
  arms <- levels(dat[[arm]])
  groups <- c(split(dat, dat[[arm]]), list(Total = dat))
  header <- c(arms, "Total")
  ns <- vapply(groups, nrow, integer(1))

  med_iqr <- function(x) {
    q <- stats::quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE)
    sprintf("%.0f (%.0f-%.0f)", q[1], q[2], q[3])
  }
  n_pct <- function(x, lvl) {
    n <- sum(x == lvl, na.rm = TRUE)
    sprintf("%d (%.1f%%)", n, 100 * n / sum(!is.na(x)))
  }
  rows <- list()
  add_row <- function(label, values) {
    rows[[length(rows) + 1]] <<- c(Characteristic = label, stats::setNames(values, header))
  }
  add_row("N", as.character(ns))
  add_row("Age, years, median (IQR)", vapply(groups, function(g) med_iqr(g$age), ""))
  add_row("Positive nodes, median (IQR)", vapply(groups, function(g) med_iqr(g$nodes), ""))
  cats <- c(sex = "Sex", node4 = "Positive nodes", extent = "Extent of local spread",
            differ = "Differentiation", obstruct = "Obstruction", perfor = "Perforation",
            adhere = "Adherence to nearby organs", surg = "Time from surgery to registration")
  for (v in names(cats)) {
    lv <- levels(dat[[v]])
    if (length(lv) == 2 && all(lv == c("No", "Yes"))) {
      add_row(paste0(cats[[v]], ", yes"), vapply(groups, function(g) n_pct(g[[v]], "Yes"), ""))
    } else {
      add_row(cats[[v]], rep("", length(header)))
      for (l in lv) {
        add_row(paste0("&nbsp;&nbsp;", l), vapply(groups, function(g) n_pct(g[[v]], l), ""))
      }
    }
  }
  for (v in c("nodes", "differ")) {
    add_row(paste0("Missing: ", if (v == "nodes") "node count" else "differentiation"),
            vapply(groups, function(g) as.character(sum(is.na(g[[v]]))), ""))
  }
  as.data.frame(do.call(rbind, rows), check.names = FALSE)
}

#' Kaplan-Meier (or 1 - cumulative incidence) plot with a numbers-at-risk table.
km_plot <- function(fit, title, ylab, times = 0:9, xmax = 9) {
  td <- broom::tidy(fit)
  td$arm <- factor(sub("^arm=", "", td$strata), levels = names(ARM_COLOURS))
  # Add time 0 so every curve starts at 1.
  start <- data.frame(time = 0, estimate = 1, conf.high = 1, conf.low = 1,
                      n.censor = 0, arm = factor(levels(td$arm), levels = levels(td$arm)))
  td <- rbind(start, td[, names(start)])
  td <- td[order(td$arm, td$time), ]

  p <- ggplot2::ggplot(td, ggplot2::aes(time, estimate, colour = arm, fill = arm)) +
    ggplot2::geom_step(ggplot2::aes(linetype = arm), linewidth = 0.8) +
    ggplot2::geom_point(data = td[td$n.censor > 0, ], shape = 3, size = 1.2,
                        show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = ARM_COLOURS, name = NULL) +
    ggplot2::scale_fill_manual(values = ARM_COLOURS, name = NULL) +
    ggplot2::scale_linetype_manual(values = ARM_LINETYPES, name = NULL) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                                limits = c(0, 1), expand = c(0, 0.01)) +
    ggplot2::scale_x_continuous(breaks = times, limits = c(0, xmax)) +
    ggplot2::labs(title = title, x = NULL, y = ylab) +
    theme_report()

  s <- summary(fit, times = times, extend = TRUE)
  risk <- data.frame(time = s$time, n = s$n.risk,
                     arm = factor(sub("^arm=", "", s$strata), levels = names(ARM_COLOURS)))
  rt <- ggplot2::ggplot(risk, ggplot2::aes(time, factor(arm, levels = rev(levels(arm))), label = n)) +
    ggplot2::geom_text(size = 3.2, colour = "#333333") +
    ggplot2::scale_x_continuous(breaks = times, limits = c(0, xmax)) +
    ggplot2::labs(x = "Years since randomisation", y = NULL, subtitle = "Number at risk") +
    theme_report(10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(colour = "#333333"))
  patchwork::wrap_plots(p, rt, ncol = 1, heights = c(4, 1))
}

#' Forest plot of hazard ratios. `hr` needs columns label, group, estimate, conf.low, conf.high.
forest_plot <- function(hr, title, xlab = "Hazard ratio (95% CI, log scale)") {
  hr$label <- factor(hr$label, levels = rev(unique(hr$label)))
  hr$txt <- fmt_ci(hr$estimate, hr$conf.low, hr$conf.high)
  ggplot2::ggplot(hr, ggplot2::aes(estimate, label)) +
    ggplot2::geom_vline(xintercept = 1, colour = "#888888", linewidth = 0.4) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = conf.low, xmax = conf.high),
                            height = 0.2, colour = "#2a78d6", linewidth = 0.6) +
    ggplot2::geom_point(size = 2.6, colour = "#2a78d6") +
    ggplot2::geom_text(ggplot2::aes(x = Inf, label = txt), hjust = 1.05, size = 3.1,
                       colour = "#333333") +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.05, 0.6))) +
    ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::labs(title = title, x = xlab, y = NULL) +
    theme_report() +
    ggplot2::theme(strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 90, face = "bold"))
}
