#' Evaluate the next cohort dose from entered patient histories
#'
#' Apply the same dose-allocation rules used by the simulator to observed
#' patient-level data. No true dose-toxicity probabilities are required.
#'
#' @param history Data frame in enrollment order with dose, dlt and followup
#'   columns. dlt=1 means a DLT has been observed; dlt=0 means none observed
#'   to date. NA is permitted only before follow-up reaches the assessment
#'   window. followup is elapsed follow-up in months. Optional prediction
#'   is 0, 1 or NA; eligible and intervened are logical columns. If omitted,
#'   eligible is inferred from prediction availability and intervened is FALSE.
#'   Use NULL for no enrolled patients.
#' @param ndoses Number of ordered dose levels, from 2 to 10.
#' @param current_dose Dose of the most recently enrolled cohort; required for
#'   nonempty history. For empty history, defaults to the configured start dose.
#' @param method Design to evaluate: bold, prediction, tite, combined, boin,
#'   or tite_boin. Default is Original BOLD.
#' @param config Design settings from bold_config().
#' @param excluded Previously eliminated dose indices. Pass the previous
#'   result's excluded_doses to preserve irreversible elimination.
#' @return A bold_recommendation list with action (start, escalate, stay,
#'   deescalate, wait or stop), next_dose (NA when waiting/stopped), reason,
#'   excluded_doses, dose_summary and patient_status. It is a next-cohort
#'   allocation result, not a final MTD selection.
#' @details Evaluation is supported at cohort boundaries of three patients.
#'   Original BOLD and prediction-assisted complete-follow-up BOLD wait until
#'   observed follow-up resolves. TITE methods use partial follow-up. Combined
#'   mode additionally waits for intervention patients to resolve.
#'
#'   Without intervention, a DLT or completed non-DLT follow-up replaces any
#'   pending prediction contribution. After intervention, an observed DLT is
#'   a known failure; no observed DLT leaves untreated status latent even after
#'   full follow-up. An intervention must have a positive prediction and an
#'   eligible patient. Actual intervention decisions must be explicitly recorded.
#'
#'   sensitivity and specificity weight available predictions. rho, eligibility
#'   probability, timing and accrual settings generate simulation histories;
#'   they do not replace observed patient records in this function. No unknown
#'   untreated outcomes should be entered. This interface implements the package's
#'   research protocol; application in a clinical trial requires protocol-specific
#'   statistical review and validation.
#' @export
#' @examples
#' h <- data.frame(dose=c(1,1,1), dlt=c(0,0,0), followup=c(3,3,3))
#' recommend_dose(h, ndoses=4, current_dose=1, method="bold")
#'
#' hp <- transform(h, prediction=c(1L,0L,0L), eligible=TRUE,
#'                 intervened=c(TRUE,FALSE,FALSE))
#' recommend_dose(hp, ndoses=4, current_dose=1, method="prediction")
recommend_dose <- function(history=NULL, ndoses, current_dose=NULL,
                           method=c("bold","prediction","tite","combined","boin","tite_boin"),
                           config=bold_config(), excluded=integer()) {
  method <- match.arg(method)
  if (!inherits(config,"bold_config")) stop("Use bold_config() for config")
  config <- do.call(bold_config,unclass(config))
  integer_one <- function(x,name,lo,hi) {
    if (!is.numeric(x) || length(x)!=1L || !is.finite(x) || x!=floor(x) || x<lo || x>hi)
      stop(name," must be an integer in ",lo,":",hi)
  }
  integer_one(ndoses,"ndoses",2,10)
  ndoses <- as.integer(ndoses)
  if (!is.numeric(excluded) || any(!is.finite(excluded) | excluded!=floor(excluded) |
                                  excluded<1 | excluded>ndoses)) stop("Invalid excluded doses")
  blocked <- rep(FALSE,ndoses)
  if (length(excluded)) blocked[seq.int(min(excluded),ndoses)] <- TRUE
  expand <- function(x,name) {
    if (!length(x) %in% c(1L,ndoses)) stop(name," must have length 1 or ndoses")
    rep(x,length.out=ndoses)
  }
  pm <- expand(if (is.null(config$prior_mean)) config$target else config$prior_mean,"prior_mean")
  ps <- expand(config$prior_pess,"prior_pess")
  if (config$start_dose>ndoses) stop("Configured start_dose exceeds ndoses")
  if (is.null(history)) history <- data.frame(dose=integer(),dlt=integer(),followup=numeric())
  if (!is.data.frame(history) || !all(c("dose","dlt","followup") %in% names(history)))
    stop("history must have dose, dlt and followup columns")
  if (anyDuplicated(names(history))) stop("history column names must be unique")
  h <- history
  N <- nrow(h)
  if (is.null(current_dose)) {
    if (N) stop("Supply current_dose for a nonempty history")
    current_dose <- config$start_dose
  }
  integer_one(current_dose,"current_dose",1,ndoses)
  current_dose <- as.integer(current_dose)
  if (!is.numeric(h$dose) || any(!is.finite(h$dose) | h$dose!=floor(h$dose) |
                                h$dose<1 | h$dose>ndoses)) stop("Invalid history dose indices")
  if ((!is.numeric(h$dlt) && !is.logical(h$dlt)) || any(!is.na(h$dlt) & !h$dlt %in% c(0,1)))
    stop("dlt must be 0, 1 or NA")
  if (!is.numeric(h$followup) || any(!is.finite(h$followup) | h$followup<0))
    stop("followup must be nonnegative and finite")
  if (any(is.na(h$dlt) & h$followup>=config$window))
    stop("Record dlt as 0 or 1 when the assessment window is complete")
  if (N>config$nmax) stop("History exceeds the configured patient cap")
  if (N %% 3 != 0) stop("Only between-cohort histories (multiples of three) are supported")
  if (N && any(utils::tail(h$dose,3)!=current_dose))
    stop("current_dose must match the most recently enrolled cohort of three")
  if (!"prediction" %in% names(h)) h$prediction <- rep(NA_integer_,N)
  if ((!is.numeric(h$prediction) && !is.logical(h$prediction)) ||
      any(!is.na(h$prediction) & !h$prediction %in% c(0,1)))
    stop("prediction must be 0, 1 or NA")
  if (!"eligible" %in% names(h)) h$eligible <- !is.na(h$prediction)
  if (!"intervened" %in% names(h)) h$intervened <- rep(FALSE,N)
  for (nm in c("eligible","intervened"))
    if (!is.logical(h[[nm]]) || anyNA(h[[nm]])) stop(nm," must be TRUE or FALSE")
  intervention_mode <- method %in% c("prediction","combined")
  if (any(h$intervened) && !intervention_mode)
    stop("Intervention histories require method prediction or combined")
  if (any(h$intervened & (!h$eligible | is.na(h$prediction) | h$prediction!=1)))
    stop("Every intervened patient requires eligibility and a positive prediction")
  dlt <- !is.na(h$dlt) & h$dlt==1
  resolved <- dlt | h$followup>=config$window
  latent <- h$intervened & !dlt
  known <- resolved & !latent
  n <- tabulate(as.integer(h$dose),nbins=ndoses)
  x <- tabulate(as.integer(h$dose[dlt]),nbins=ndoses)
  complete <- tabulate(as.integer(h$dose[known]),nbins=ndoses)
  pending <- tabulate(as.integer(h$dose[!resolved]),nbins=ndoses)
  ds <- data.frame(dose=seq_len(ndoses),n=n,observed_dlts=x,pending=pending,
    latent_intervention=tabulate(as.integer(h$dose[latent]),nbins=ndoses),
    mean_untreated_risk=NA_real_,overdose_probability=NA_real_,excluded=blocked)
  status <- data.frame(patient=seq_len(N),dose=h$dose,
    status=ifelse(latent,"latent_intervention",ifelse(dlt,"observed_dlt",
      ifelse(resolved,"observed_no_dlt","pending"))))
  finish <- function(action,next_dose=NA_integer_,reason) {
    structure(list(action=action,next_dose=as.integer(next_dose),reason=reason,
      method=method,current_dose=current_dose,excluded_doses=which(blocked),
      dose_summary=ds,patient_status=status,config=config),class="bold_recommendation")
  }
  if (blocked[1]) return(finish("stop",reason="Lowest dose has been eliminated."))
  if (N>=config$nmax) return(finish("stop",reason="Patient cap reached; no further cohort dose is assigned."))
  if (!N) {
    if (blocked[current_dose]) return(finish("stop",reason="The configured starting dose is eliminated."))
    return(finish("start",current_dose,reason="No patients enrolled; use the configured starting dose."))
  }
  if (method %in% c("bold","boin","prediction") && any(!resolved))
    return(finish("wait",reason="Complete-follow-up mode: await all pending outcomes."))
  if (method=="combined" && any(h$intervened & !resolved))
    return(finish("wait",reason="Combined protocol: await intervention-patient outcome resolution."))
  q <- labels <- vector("list",ndoses)
  stft <- numeric(ndoses)
  for (j in seq_len(ndoses)) {
    rows <- which(h$dose==j & !known)
    q[[j]] <- ifelse(latent[rows],1,pmax(0,1-h$followup[rows]/config$window))
    labels[[j]] <- as.integer(ifelse(latent[rows],1,
      ifelse(method=="combined" & h$eligible[rows] & !is.na(h$prediction[rows]),
             h$prediction[rows],-1)))
    p <- which(h$dose==j & !resolved & !h$intervened)
    stft[j] <- sum(pmin(1,h$followup[p]/config$window))
  }
  mode <- c(bold=4L,prediction=5L,tite=1L,combined=6L,boin=3L,tite_boin=0L)[method]
  z <- observed_history_decision(as.integer(n),as.integer(x),as.integer(complete),
    as.integer(pending),stft,q,labels,pm*ps,(1-pm)*ps,blocked,current_dose,
    config$target,config$overdose_cutoff,config$tau,unname(mode),
    config$sensitivity,config$specificity,config$min_two,config$pending_limit)
  blocked <- as.logical(z$excluded)
  ds$excluded <- blocked
  ds$mean_untreated_risk <- as.numeric(z$mean)
  ds$overdose_probability <- as.numeric(z$overdose_probability)
  if (z$safety) return(finish("stop",reason="Overdose elimination excludes the lowest dose."))
  if (z$suspend) return(finish("wait",reason="Pending-outcome or escalation-readiness rule requires more follow-up."))
  if (n[current_dose]>=config$early_n && z$candidate==current_dose && !blocked[current_dose])
    return(finish("stop",reason="The configured early-stopping criterion is met."))
  action <- if (z$candidate>current_dose) "escalate" else if (z$candidate<current_dose) "deescalate" else "stay"
  finish(action,z$candidate,reason=paste0("Next cohort: dose ",z$candidate," (",action,")."))
}

#' @export
print.bold_recommendation <- function(x, ...) {
  cat("Method:",x$method,"| Action:",x$action)
  if (!is.na(x$next_dose)) cat(" | Next dose:",x$next_dose)
  cat("\n",x$reason,"\n",sep="")
  invisible(x)
}

#' Flag intervention eligibility from external predictions
#'
#' Apply the modeled positive-prediction intervention policy. This function
#' does not fit a classifier or perform an intervention.
#' @param prediction Vector of binary predictions (0/1); NA means unavailable.
#' @param eligible Logical vector, or scalar recycled to all predictions.
#' @return A data frame containing prediction, eligible, intervene (logical,
#'   NA if an eligible patient's prediction is unavailable), and status.
#' @export
#' @examples
#' plan_intervention(c(1L,0L,NA_integer_),eligible=c(TRUE,TRUE,FALSE))
plan_intervention <- function(prediction,eligible=TRUE) {
  if ((!is.numeric(prediction) && !is.logical(prediction)) ||
      any(!is.na(prediction) & !prediction %in% c(0,1))) stop("prediction must be 0, 1 or NA")
  if (!is.logical(eligible) || anyNA(eligible) ||
      !length(eligible) %in% c(1L,length(prediction))) stop("eligible must be logical and match prediction")
  eligible <- rep(eligible,length.out=length(prediction))
  intervene <- ifelse(!eligible,FALSE,ifelse(is.na(prediction),NA,prediction==1))
  status <- ifelse(!eligible,"not_eligible",ifelse(is.na(prediction),"await_prediction",
                     ifelse(prediction==1,"positive","negative")))
  data.frame(prediction=as.integer(prediction),eligible=eligible,intervene=intervene,status=status)
}
