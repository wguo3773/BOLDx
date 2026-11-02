#' Simulate BOLD extensions and matched comparators
#'
#' Evaluate dose-finding accuracy, enrollment, actual toxicity, allocation above
#' the MTD, and calendar duration for BOLD extensions and their comparators.
#'
#' @details Fixed cohorts of three are used. Intervention failure is an observed
#' DLT, not an extra latent branch. In combined mode, intervention follow-up
#' resolves before the next dose decision. Mixture weights use sequential
#' child-specific plug-in scores; the full calculation and assumptions are
#' documented in the package's online Methods document.
#'
#' @param truth Nondecreasing vector of 2 to 10 untreated DLT probabilities,
#'   strictly between zero and one.
#' @param p.true Alias for truth, for a user-specified dose-toxicity scenario.
#'   Supply one name, or identical values under both names.
#' @param method One or more of prediction, tite, combined, bold, boin, tite_boin.
#' @param config A configuration from bold_config().
#' @param nsim Number of replicate trials per method, default 1000.
#' @param seed Nonnegative integer seed. The caller's random state is restored.
#' @param mtd True optimal dose index for evaluation. NULL chooses the dose
#'   closest to target (lower index breaks ties). Set 0 for all-too-toxic;
#'   K+1 means all doses are below target and the highest dose is optimal.
#'   This value affects summaries only, never dose decisions.
#' @return A bold_simulation list containing summary, trials, config, truth,
#'   mtd, and seed. Summary contains scenario-level point estimates. Optional
#'   Monte Carlo precision diagnostics are retained separately in mcse.
#' @section Uncertainty:
#'   For one fixed scenario, mcse reports Monte Carlo standard errors across
#'   the nsim trials, not between-scenario SD. Accuracy SE is 100 times the
#'   sample SD of trial-level correct-selection indicators divided by
#'   sqrt(nsim). For mean N, actual DLTs and duration, SE is the trial-level
#'   sample SD divided by sqrt(nsim). Above-MTD SE uses the influence values
#'   of the patient-weighted ratio, not an unweighted mean of trial ratios.
#'   Raw trial records are retained in trials for other summaries.
#' @export
#' @examples
#' z <- simulate_bold(c(.10, .25, .45, .60),
#'   method = c("prediction", "bold", "boin"), nsim = 10, seed = 1)
#' z$summary
simulate_bold <- function(truth = NULL, method = c("prediction", "tite", "combined"),
                          config = bold_config(), nsim = 1000L, seed = 1L,
                          mtd = NULL, p.true = NULL) {
  if (!is.null(p.true)) {
    if (!is.null(truth) && !isTRUE(all.equal(truth,p.true)))
      stop("truth and p.true must agree when both are supplied")
    truth <- p.true
  }
  if (!inherits(config, "bold_config")) stop("Use bold_config() for config")
  config <- do.call(bold_config, unclass(config))
  if (!is.numeric(truth) || length(truth) < 2 || length(truth) > 10 ||
      any(!is.finite(truth) | truth <= 0 | truth >= 1) || any(diff(truth) < 0))
    stop("truth must contain 2-10 nondecreasing probabilities strictly in (0,1)")
  k <- length(truth)
  if (config$start_dose > k) stop("start_dose exceeds the number of doses")
  expand <- function(x, name) {
    if (!length(x) %in% c(1L, k)) stop(name, " must have length 1 or number of doses")
    rep(x, length.out = k)
  }
  pm <- expand(if (is.null(config$prior_mean)) config$target else config$prior_mean, "prior_mean")
  pp <- expand(config$prior_pess, "prior_pess")
  for (nm in c("nsim", "seed")) {
    x <- get(nm)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 ||
        x != floor(x) || x > .Machine$integer.max) stop(nm, " must be a nonnegative integer")
  }
  if (nsim < 2) stop("nsim must be at least 2 for Monte Carlo standard errors")
  modes <- c(prediction = 5L, tite = 1L, combined = 6L,
             bold = 4L, boin = 3L, tite_boin = 0L)
  if (!is.character(method) || !length(method) || anyNA(method) ||
      any(!method %in% names(modes))) stop("Unknown method")
  method <- unique(method)
  if (is.null(mtd)) mtd <- which.min(abs(truth - config$target))
  if (!is.numeric(mtd) || length(mtd) != 1 || !is.finite(mtd) ||
      mtd != floor(mtd) || mtd < 0 || mtd > k+1) stop("mtd must be in 0:(K+1)")
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  trials <- lapply(method, function(m) {
    z <- as.data.frame(calendar_combined(truth, phi = config$target,
      nsim = as.integer(nsim), mode = unname(modes[m]), seed = as.integer(seed),
      nmax = config$nmax, window = config$window, rate = config$accrual_rate,
      gamma = config$overdose_cutoff, tau = config$tau,
      sensitivity = config$sensitivity, specificity = config$specificity,
      timing = match(config$timing, c("weibull", "uniform", "late")) - 1L,
      early_n = as.integer(min(config$nmax + 1, config$early_n)),
      min_two = config$min_two, eligibility = config$eligibility, efficacy = config$rho,
      prior_means = pm, prior_sizes = pp, start = config$start_dose,
      pending_limit = config$pending_limit, arrival = as.integer(config$arrival == "exponential")))
    z$above_mtd <- if (mtd == 0) z$n else if (mtd >= k) 0 else
      rowSums(z[paste0("n", seq.int(mtd+1, k))])
    z$correct <- z$selected == min(mtd, k)
    z$method <- m
    z$replicate <- seq_len(nsim)
    z
  })
  trials <- do.call(rbind, trials)
  sm <- do.call(rbind, lapply(method, function(m) {
    z <- trials[trials$method == m, ]
    se <- function(x) stats::sd(x)/sqrt(nsim)
    over <- sum(z$above_mtd)/sum(z$n)
    timed <- m %in% c("tite", "combined", "tite_boin")
    data.frame(method = m, nsim = nsim,
      accuracy_pct = 100*mean(z$correct), accuracy_se = 100*se(z$correct),
      mean_n = mean(z$n), mean_n_se = se(z$n),
      actual_dlts = mean(z$dlt), actual_dlts_se = se(z$dlt),
      above_mtd_pct = 100*over,
      above_mtd_se = 100*se(z$above_mtd-over*z$n)/mean(z$n),
      duration_months = if (timed) mean(z$duration) else NA_real_,
      duration_se = if (timed) se(z$duration) else NA_real_)
  }))
  rownames(trials) <- rownames(sm) <- NULL
  diagnostic <- grepl("_se$", names(sm))
  mcse <- sm[c("method", "nsim", names(sm)[diagnostic])]
  sm <- sm[!diagnostic]
  structure(list(summary = sm, mcse = mcse, trials = trials, config = config,
    truth = truth, mtd = mtd, seed = seed), class = "bold_simulation")
}

#' @export
print.bold_simulation <- function(x, ...) {
  print(x$summary, row.names = FALSE, ...)
  invisible(x)
}
