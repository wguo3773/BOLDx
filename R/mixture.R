#' Inspect sequential weighted Beta components
#'
#' Compute mixture weights, Beta component parameters, the mixture mean,
#' and the overdose probability for prediction and partial-follow-up histories.
#'
#' @details Components are weighted by sequential child-specific plug-in scores.
#' Patient order is therefore part of the calculation. Known outcomes should
#' first be incorporated into a and b. An observed
#' intervention failure is a known DLT: increment a rather than branching.
#'
#' @param a,b Positive Beta parameters after resolved, nonlatent outcomes.
#' @param q Remaining DLT-window fraction for each latent patient in [0,1].
#'   Use 1 for an intervention patient without an observed DLT.
#' @param labels Prediction labels: 1 positive, 0 negative, -1 uninformative.
#' @param sensitivity,specificity Prediction operating characteristics.
#' @param target Tail threshold for computing overdose probability.
#' @return List containing components, mean, and overdose_probability.
#' @export
#' @examples
#' beta_mixture(.75, 2.25, q = c(1, .5), labels = c(1L, -1L))
beta_mixture <- function(a, b, q = numeric(), labels = rep(-1L, length(q)),
                         sensitivity = .9, specificity = .98, target = .25) {
  for (nm in c("a", "b")) {
    z <- get(nm)
    if (!is.numeric(z) || length(z) != 1 || !is.finite(z) || z <= 0)
      stop(nm, " must be positive")
  }
  if (!is.numeric(q) || any(!is.finite(q) | q < 0 | q > 1)) stop("q must be in [0,1]")
  if (length(labels) != length(q) || anyNA(labels) || any(!labels %in% c(-1,0,1)))
    stop("labels must be -1, 0, or 1 and match q")
  for (nm in c("sensitivity", "specificity", "target")) {
    z <- get(nm)
    if (!is.numeric(z) || length(z) != 1 || !is.finite(z) || z < 0 || z > 1)
      stop(nm, " must be in [0,1]")
  }
  p <- posterior_check(a, b, q, as.integer(labels), sensitivity, specificity, target)
  k <- seq.int(0, length(q))
  list(components = data.frame(weight = as.numeric(p$weights),
    a = a+k, b = b+length(q)-k), mean = p$mean, overdose_probability = p$cp)
}
