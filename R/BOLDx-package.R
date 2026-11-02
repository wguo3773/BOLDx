#' BOLDx: prediction and delayed-outcome dose-finding simulations
#'
#' Simulation-based operating-characteristic evaluation of BOLD extensions
#' with prediction-guided intervention and delayed toxicity outcomes.
#'
#' @references Wang, G.-M. and Tatsuoka, C. (2026). Bayesian Ordered Lattice
#' Design for Phase I Clinical Trials. Statistics in Medicine.
#' \doi{10.1002/sim.70456}.
#' @seealso \code{\link{bold_config}}, \code{\link{simulate_bold}},
#'   \code{\link{simulate_scenarios}}, \code{\link{beta_mixture}}
#'   \code{\link{recommend_dose}}, \code{\link{plan_intervention}}
#' @useDynLib BOLDx, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
