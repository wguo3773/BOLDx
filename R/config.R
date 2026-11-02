#' Configure a BOLD extension simulation
#'
#' The prior at dose j is Beta(prior_mean[j] * prior_pess[j],
#' (1 - prior_mean[j]) * prior_pess[j]). True untreated risks are separate
#' simulation inputs and are never supplied to dose-assignment calculations.
#'
#' @param target Target untreated DLT probability, strictly between 0 and 5/7.
#' @param nmax Patient cap, a positive multiple of three.
#' @param tau BOLD decision-tail target in (0,1); not the DLT target.
#' @param overdose_cutoff Posterior overdose elimination cutoff in (0,1).
#' @param early_n Stop when this many patients have received the current dose
#'   and the next recommendation is unchanged. Use Inf to disable.
#' @param sensitivity,specificity Prediction sensitivity and specificity for
#'   latent untreated DLT status, each in [0,1].
#' @param eligibility Probability a patient is eligible for prediction.
#' @param rho Probability an intervention prevents an otherwise-occurring DLT.
#' @param window DLT assessment window, in months.
#' @param accrual_rate Mean number of arrivals per month.
#' @param arrival Arrival-gap distribution, uniform or exponential.
#' @param timing Generative conditional DLT-time distribution. The working
#'   follow-up model remains uniform for all choices.
#' @param prior_mean,prior_pess Dose-specific vectors, or scalars recycled to
#'   all doses. NULL prior_mean uses target; PESS defaults to 3.
#' @param start_dose Starting dose index, validated against each scenario.
#' @param pending_limit Maximum pending fraction at the current dose before
#'   suspension; defaults to 0.5.
#' @param min_two Require two completed outcomes before escalation.
#' @return A validated configuration list of class bold_config.
#' @export
#' @examples
#' cfg <- bold_config(nmax = 21, rho = 0.8)
#' cfg$target
bold_config <- function(target = 0.25, nmax = 21L, tau = 0.50,
                        overdose_cutoff = 0.90, early_n = 12L,
                        sensitivity = 0.90, specificity = 0.98,
                        eligibility = 1, rho = 0.80, window = 3,
                        accrual_rate = 2, arrival = c("uniform", "exponential"),
                        timing = c("weibull", "uniform", "late"),
                        prior_mean = NULL, prior_pess = 3,
                        start_dose = 1L, pending_limit = 0.5, min_two = FALSE) {
  prob <- function(x, name, open = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        if (open) x <= 0 || x >= 1 else x < 0 || x > 1)
      stop(name, " must be a probability", call. = FALSE)
  }
  for (nm in c("target", "tau", "overdose_cutoff")) prob(get(nm), nm, TRUE)
  if (target >= 5/7) stop("target must be below 5/7 for BOIN boundaries")
  for (nm in c("sensitivity", "specificity", "eligibility", "rho", "pending_limit"))
    prob(get(nm), nm)
  for (nm in c("nmax", "start_dose")) {
    x <- get(nm)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 ||
        x != floor(x)) stop(nm, " must be a positive integer")
  }
  if (nmax %% 3 != 0 || nmax > 1000) stop("nmax must be a multiple of 3, at most 1000")
  if (!is.numeric(early_n) || length(early_n) != 1L || is.na(early_n) ||
      early_n < 1 || (is.finite(early_n) && early_n != floor(early_n)))
    stop("early_n must be a positive integer or Inf")
  for (nm in c("window", "accrual_rate")) {
    x <- get(nm)
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0)
      stop(nm, " must be positive and finite")
  }
  if (!is.logical(min_two) || length(min_two) != 1L || is.na(min_two))
    stop("min_two must be TRUE or FALSE")
  if (!is.null(prior_mean) && (!is.numeric(prior_mean) || !length(prior_mean) ||
      any(!is.finite(prior_mean) | prior_mean <= 0 | prior_mean >= 1)))
    stop("prior_mean values must be strictly between 0 and 1")
  if (!is.numeric(prior_pess) || !length(prior_pess) ||
      any(!is.finite(prior_pess) | prior_pess <= 0)) stop("prior_pess must be positive")
  structure(list(target = target, nmax = as.integer(nmax), tau = tau,
    overdose_cutoff = overdose_cutoff, early_n = early_n,
    sensitivity = sensitivity, specificity = specificity, eligibility = eligibility,
    rho = rho, window = window, accrual_rate = accrual_rate,
    arrival = match.arg(arrival), timing = match.arg(timing),
    prior_mean = prior_mean, prior_pess = prior_pess,
    start_dose = as.integer(start_dose), pending_limit = pending_limit,
    min_two = min_two), class = "bold_config")
}
