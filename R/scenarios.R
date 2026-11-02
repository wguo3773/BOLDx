#' Run a panel of dose-toxicity scenarios
#'
#' Results include scenario-level estimates and mean/SD summaries across random
#' scenarios within each method and MTD category. Use explicit MTD labels for
#' all-too-toxic and below-range cases.
#'
#' @param scenarios Numeric matrix with one nondecreasing risk vector per row.
#' @param mtd Optional integer vector giving each scenario's optimal dose.
#' @param method,config,nsim,seed As in simulate_bold(). Seeds increment by
#'   scenario index; methods within a scenario share the same seed.
#' @return List with summary and trials data frames, each containing scenario,
#'   and aggregate containing means and between-scenario standard deviations.
#' @export
#' @examples
#' x <- rbind(c(.1,.25,.45,.6),c(.05,.1,.25,.45))
#' simulate_scenarios(x, mtd=c(2,3), method="prediction", nsim=5)$summary
simulate_scenarios <- function(scenarios, mtd = NULL,
                               method = c("prediction", "tite", "combined"),
                               config = bold_config(), nsim = 100, seed = 1) {
  if (!is.matrix(scenarios) || !is.numeric(scenarios) || nrow(scenarios) < 1)
    stop("scenarios must be a nonempty numeric matrix")
  n <- nrow(scenarios)
  if (!is.null(mtd) && length(mtd) != n) stop("mtd must match the number of scenarios")
  if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed) || seed < 0 ||
      seed != floor(seed) || seed+n-1 > .Machine$integer.max) stop("Invalid seed")
  out <- lapply(seq_len(n), function(i) {
    z <- simulate_bold(scenarios[i,], method, config, nsim, seed+i-1,
      if (is.null(mtd)) NULL else mtd[i])
    z$summary$scenario <- i
    z$summary$mtd <- z$mtd
    z$trials$scenario <- i
    z$trials$mtd <- z$mtd
    z
  })
  sm <- do.call(rbind,lapply(out,`[[`,"summary"))
  list(summary=sm, trials=do.call(rbind,lapply(out,`[[`,"trials")),
       aggregate=summarise_scenarios(sm))
}

#' Summarise operating characteristics across random scenarios
#'
#' Compute the equally weighted mean and sample standard deviation of
#' scenario-level estimates, each obtained from its simulated trial replicates.
#'
#' @section Calculation:
#' For J scenario-level estimates Q_j, the reported values are
#' \deqn{\overline{Q}=J^{-1}\sum_{j=1}^{J} Q_j,\qquad
#' SD(Q)=\sqrt{\frac{\sum_{j=1}^{J}(Q_j-\overline{Q})^2}{J-1}}.}
#'
#' @details The SD measures between-scenario dispersion, not a standard error.
#' For 50 scenarios with 100 replicates each, the SD is calculated from the
#' 50 scenario-level estimates, not from the 5,000 individual trials. It includes
#' both scenario heterogeneity and finite-replicate simulation noise. Group by
#' MTD category, upper toxicity gap, and design settings where appropriate.
#' A group with only one scenario has undefined SD, reported as NA.
#' Above-MTD percentages are first patient-weighted within each scenario, then
#' equally weighted across scenarios. Missing duration for complete-follow-up
#' methods is retained as NA.
#'
#' @param x A data frame of scenario-level summaries, or the result returned
#'   by simulate_scenarios().
#' @param group_by Columns defining the reporting groups; must include method.
#'   Additional columns, such as delta, can be attached to x before summarising.
#' @return A data frame of point-estimate means and matching *_sd columns,
#'   with n_scenarios and the minimum/maximum replicates per scenario.
#' @export
#' @examples
#' panel <- simulate_scenarios(rbind(c(.1,.25,.45), c(.15,.25,.4)),
#'   mtd=c(2,2), method="prediction", nsim=5)
#' panel$aggregate
summarise_scenarios <- function(x, group_by=c("method","mtd")) {
  if (is.list(x) && !is.data.frame(x) && !is.null(x$summary)) x <- x$summary
  metrics <- c("accuracy_pct","mean_n","actual_dlts","above_mtd_pct","duration_months")
  if (!is.data.frame(x) || !nrow(x)) stop("x must contain scenario-level summaries")
  if (!is.character(group_by) || !length(group_by) || anyNA(group_by) ||
      !"method" %in% group_by || anyDuplicated(group_by))
    stop("group_by must contain unique column names including method")
  needed <- unique(c(group_by,"scenario","nsim",metrics))
  if (!all(needed %in% names(x))) stop("Missing required summary columns")
  if (anyNA(x[c(group_by,"scenario")])) stop("Grouping and scenario identifiers cannot be missing")
  if (anyDuplicated(x[unique(c(group_by,"scenario"))]))
    stop("Each scenario must appear once per reporting group")
  if (!is.numeric(x$nsim) || any(!is.finite(x$nsim) | x$nsim < 2 | x$nsim != floor(x$nsim)))
    stop("nsim must record at least two replicates per scenario")
  for (nm in metrics) {
    if (!is.numeric(x[[nm]]) || any(is.infinite(x[[nm]])) ||
        (nm != "duration_months" && anyNA(x[[nm]]))) stop("Invalid metric: ",nm)
  }
  groups <- interaction(lapply(x[group_by], function(v) factor(v, levels=unique(v))),
                        drop=TRUE, lex.order=TRUE)
  out <- lapply(split(x,groups), function(d) {
    z <- d[1,group_by,drop=FALSE]
    z$n_scenarios <- nrow(d)
    z$replicates_min <- min(d$nsim)
    z$replicates_max <- max(d$nsim)
    for (nm in metrics) {
      if (anyNA(d[[nm]]) && !all(is.na(d[[nm]])))
        stop("Partially missing metric within a group: ",nm)
      z[[nm]] <- if (all(is.na(d[[nm]]))) NA_real_ else mean(d[[nm]])
      sd_name <- switch(nm, accuracy_pct="accuracy_sd", mean_n="mean_n_sd",
        actual_dlts="actual_dlts_sd", above_mtd_pct="above_mtd_sd",
        duration_months="duration_sd")
      z[[sd_name]] <- if (all(is.na(d[[nm]]))) NA_real_ else stats::sd(d[[nm]])
    }
    z
  })
  out <- do.call(rbind,out)
  rownames(out) <- NULL
  out
}
