#' Generate ordered toxicity scenarios from adjacent-dose gaps
#'
#' Generate a reproducible panel for one MTD category and gap combination.
#' This is a defined uniform-tail generator, not a claim to sample all possible
#' dose-toxicity curves or to reproduce previously saved manuscript scenarios.
#' @param ndoses Number of dose levels, 2 to 10; default 4.
#' @param target Target untreated DLT probability; default 0.25.
#' @param mtd MTD category. For 1:ndoses, fix that dose at target. Zero means
#'   all too toxic; ndoses+1 means all doses below target. Default 2.
#' @param lower_gap,upper_gap Positive adjacent-dose risk gaps, default 0.10
#'   and 0.15. For an interior MTD j, set p[j-1]=target-lower_gap and
#'   p[j+1]=target+upper_gap when those neighbors exist.
#' @param risk_bounds Lower and upper limits for generated risks, default
#'   c(0.01,0.70). Must lie strictly inside (0,1) and contain target.
#' @param nscenarios Number of scenarios for this category and gap setting,
#'   default 50.
#' @param seed Nonnegative integer seed; the caller's random state is restored.
#' @return Numeric matrix of true risk vectors, one scenario per row.
#' @details Remaining lower doses are sorted uniform draws between the lower
#'   risk bound and the fixed lower neighbor; remaining upper doses are sorted
#'   uniform draws between the fixed upper neighbor and upper bound. For
#'   mtd=0, the first risk is target+upper_gap and remaining risks are uniform
#'   draws above it. For mtd=ndoses+1, the last risk is target-lower_gap and
#'   remaining risks are uniform draws below it. Gaps requiring a neighbor
#'   must leave a nonempty interval inside risk_bounds. Some boundary settings
#'   have no free doses and produce identical curves; inspect the matrix.
#' @export
#' @examples
#' generate_scenarios(nscenarios=3, mtd=2, upper_gap=.15, seed=42)
generate_scenarios <- function(ndoses=4L, target=.25, mtd=2L,
                               lower_gap=.10, upper_gap=.15,
                               risk_bounds=c(.01,.70), nscenarios=50L, seed=1L) {
  int <- function(x,lo,hi) is.numeric(x) && length(x)==1L && is.finite(x) &&
    x==floor(x) && x>=lo && x<=hi
  if (!int(ndoses,2,10) || !int(mtd,0,ndoses+1) ||
      !int(nscenarios,1,.Machine$integer.max) || !int(seed,0,.Machine$integer.max))
    stop("Invalid ndoses, mtd, nscenarios or seed")
  if (!is.numeric(target) || length(target)!=1L || !is.finite(target) ||
      target<=0 || target>=1) stop("target must be in (0,1)")
  if (!is.numeric(risk_bounds) || length(risk_bounds)!=2L ||
      any(!is.finite(risk_bounds)) || risk_bounds[1]<=0 || risk_bounds[2]>=1 ||
      risk_bounds[1]>=target || risk_bounds[2]<=target)
    stop("risk_bounds must lie in (0,1) and contain target")
  for (g in list(lower_gap,upper_gap))
    if (!is.numeric(g) || length(g)!=1L || !is.finite(g) || g<=0)
      stop("Gaps must be positive finite scalars")
  lo <- target-lower_gap; hi <- target+upper_gap
  if ((mtd>1 && lo<=risk_bounds[1]) ||
      (mtd<ndoses && hi>=risk_bounds[2]))
    stop("Required adjacent risks must be strictly inside risk_bounds")
  had <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if (had) old <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if (had) assign(".Random.seed",old,envir=.GlobalEnv)
    else if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(as.integer(seed))
  out <- t(vapply(seq_len(nscenarios),function(i) {
    if (mtd==0) return(c(hi,sort(stats::runif(ndoses-1,hi,risk_bounds[2]))))
    if (mtd==ndoses+1)
      return(c(sort(stats::runif(ndoses-1,risk_bounds[1],lo)),lo))
    lower <- if (mtd>1) c(sort(stats::runif(mtd-2,risk_bounds[1],lo)),lo) else numeric()
    upper <- if (mtd<ndoses) c(hi,sort(stats::runif(ndoses-mtd-1,hi,risk_bounds[2]))) else numeric()
    c(lower,target,upper)
  },numeric(ndoses)))
  colnames(out) <- paste0("p",seq_len(ndoses))
  out
}

#' Simulate random scenario panels defined by toxicity gaps
#'
#' Generate and simulate nscenarios scenarios for each MTD and gap combination.
#' Defaults are 50 scenarios and 1000 trials per scenario per method.
#' @param ndoses,mtd,lower_gap,upper_gap,risk_bounds,nscenarios Arguments to
#'   generate_scenarios(). Here mtd and either gap may be vectors; all distinct
#'   combinations are evaluated separately.
#' @param method,config,nsim,seed Arguments to simulate_bold(). The target
#'   probability is taken from config. Scenario seeds increment globally;
#'   a separate deterministic generation seed is used for each setting.
#' @return A list with scenarios (matrix), settings (scenario labels and gaps),
#'   summary (per-scenario estimates), mcse (within-scenario Monte Carlo SE),
#'   trials (individual trials), and aggregate (mean and SD across scenarios
#'   grouped separately by method, MTD category and gap combination).
#' @details Increasing nscenarios explores scenario heterogeneity; increasing
#'   nsim improves within-scenario Monte Carlo precision. Aggregate SD is not
#'   the SE of a mean and includes finite-trial simulation noise. Patient caps,
#'   predictor performance, rho and timing are set through config. To vary
#'   configurations, run this function separately and label results accordingly.
#' @export
#' @examples
#' z <- simulate_random_scenarios(mtd=2, upper_gap=.15,
#'   nscenarios=2, nsim=5, method="prediction", seed=42)
#' z$settings
#' z$aggregate
simulate_random_scenarios <- function(ndoses=4L, mtd=2L,
                                      lower_gap=.10, upper_gap=.15,
                                      risk_bounds=c(.01,.70), nscenarios=50L,
                                      method=c("bold","prediction","tite","combined"),
                                      config=bold_config(), nsim=1000L, seed=1L) {
  if (!inherits(config,"bold_config")) stop("Use bold_config() for config")
  config <- do.call(bold_config,unclass(config))
  for (x in list(mtd,lower_gap,upper_gap))
    if (!is.numeric(x) || !length(x) || any(!is.finite(x))) stop("Invalid category/gap vectors")
  grid <- expand.grid(mtd=unique(mtd),lower_gap=unique(lower_gap),
                      upper_gap=unique(upper_gap),KEEP.OUT.ATTRS=FALSE)
  if (!is.numeric(seed) || length(seed)!=1L || !is.finite(seed) || seed<0 ||
      seed!=floor(seed) || !is.numeric(nscenarios) || length(nscenarios)!=1L ||
      !is.finite(nscenarios) || nscenarios<1 || nscenarios!=floor(nscenarios) ||
      seed+nrow(grid)*nscenarios-1>.Machine$integer.max) stop("Invalid seed/scenario count")
  defs <- lapply(seq_len(nrow(grid)),function(i)
    generate_scenarios(ndoses,config$target,grid$mtd[i],grid$lower_gap[i],
      grid$upper_gap[i],risk_bounds,nscenarios,seed+i-1))
  scenarios <- do.call(rbind,defs)
  settings <- grid[rep(seq_len(nrow(grid)),each=nscenarios),,drop=FALSE]
  rownames(settings) <- NULL
  settings$scenario <- seq_len(nrow(settings))
  z <- simulate_scenarios(scenarios,mtd=settings$mtd,method=method,
                          config=config,nsim=nsim,seed=seed)
  attach <- function(d) {
    d$lower_gap <- settings$lower_gap[d$scenario]
    d$upper_gap <- settings$upper_gap[d$scenario]
    d
  }
  z$summary <- attach(z$summary); z$mcse <- attach(z$mcse)
  z$trials <- attach(z$trials)
  z$aggregate <- summarise_scenarios(z$summary,c("method","mtd","lower_gap","upper_gap"))
  z$scenarios <- scenarios; z$settings <- settings
  z
}
