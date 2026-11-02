#' Plan the next dose and prediction-guided intervention together
#'
#' Combine an observed-history dose recommendation with intervention flags for
#' incoming patients. Incoming predictions never alter the historical posterior.
#' @param history,ndoses,current_dose,config,excluded Arguments passed to
#'   recommend_dose().
#' @param method One of bold, prediction, tite or combined.
#' @param prediction External binary predictions for the incoming cohort, or
#'   NULL when predictions have not yet been obtained. A supplied cohort must
#'   contain three predictions; NA denotes an unavailable prediction.
#' @param eligible Logical scalar or vector for the incoming cohort.
#' @return A list with dose_recommendation, intervention_plan and ready.
#'   intervention_plan contains the candidate dose, intervene flags and status.
#'   No intervention is scheduled when dose allocation is waiting or stopped.
#'   ready is TRUE only when a dose is assigned and all required incoming
#'   prediction decisions are available. Actual interventions must subsequently
#'   be recorded in history, not inferred from these proposed flags.
#' @details Original BOLD and TITE-BOLD do not use preventive intervention.
#'   Prediction and combined modes flag eligible positive predictions. If
#'   predictions depend on dose, obtain them at the returned candidate dose
#'   before treating patients. This function neither trains a classifier nor
#'   supplies a final MTD declaration. Protocol-specific statistical review and
#'   validation are required before clinical application.
#' @export
#' @examples
#' cfg <- bold_config(gamma=.90, tau=.50)
#' h <- data.frame(dose=c(1,1,1), dlt=c(0,0,0), followup=c(3,3,3))
#' p <- plan_bold(h, 4, 1, method="prediction", config=cfg,
#'                prediction=c(1L,0L,0L))
#' p$dose_recommendation
#' p$intervention_plan
plan_bold <- function(history=NULL, ndoses, current_dose=NULL,
                      method=c("bold","prediction","tite","combined"),
                      config=bold_config(), excluded=integer(),
                      prediction=NULL, eligible=TRUE) {
  method <- match.arg(method)
  d <- recommend_dose(history, ndoses, current_dose, method, config, excluded)
  uses_prediction <- method %in% c("prediction","combined")
  if (is.null(prediction)) prediction <- rep(NA_integer_,3L)
  if (length(prediction)!=3L) stop("Supply three incoming predictions, or NULL")
  p <- plan_intervention(prediction, eligible)
  assigned <- d$action %in% c("start","escalate","deescalate","stay")
  p$dose <- rep(d$next_dose,3L)
  if (!uses_prediction) {
    p$intervene <- FALSE
    p$status <- "no_intervention_in_design"
  }
  if (!assigned) {
    p$intervene <- FALSE
    p$status <- paste0("dose_",d$action)
  }
  list(dose_recommendation=d, intervention_plan=p,
       ready=assigned && !anyNA(p$intervene))
}
