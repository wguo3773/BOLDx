#' Run a panel of dose-toxicity scenarios
#'
#' Results are returned separately for each scenario, without averaging across
#' MTD categories. Use explicit MTD labels for all-too-toxic and below-range cases.
#'
#' @param scenarios Numeric matrix with one nondecreasing risk vector per row.
#' @param mtd Optional integer vector giving each scenario's optimal dose.
#' @param method,config,nsim,seed As in simulate_bold(). Seeds increment by
#'   scenario index; methods within a scenario share the same seed.
#' @return List with summary and trials data frames, each containing scenario.
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
  list(summary=do.call(rbind,lapply(out,`[[`,"summary")),
       trials=do.call(rbind,lapply(out,`[[`,"trials")))
}
