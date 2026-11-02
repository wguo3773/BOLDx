# Methods and implementation contract

## Untreated risk, prediction, and intervention

At dose j, p_j = Pr(Z=1 | dose=j, no preventive intervention) is the
data-generating probability of a DLT within the assessment window. It is not
the post-intervention risk or a fitted posterior mean. Draw latent Z from
Bernoulli(p_j). Draw eligibility E from Bernoulli(e), independently. For eligible
patients, Pr(L=1 | Z=1)=sensitivity and Pr(L=0 | Z=0)=specificity. Intervene only
when E=1 and L=1. A true DLT is prevented with probability rho. Therefore

    Pr(actual DLT | dose=j) = p_j * (1 - e * sensitivity * rho).

Specificity affects the number of false-positive interventions and mixture
updates, but not this marginal risk under the assumption of no intervention harm.
Intervention success is randomized conditional on a true-positive prediction.
Latent outcomes and true dose probabilities are never passed to decision rules.

An intervention failure produces an observed DLT and a regular Beta update.
Intervention patients without an observed DLT have unknown untreated status and
enter two latent-outcome branches. The simulation knows whether prevention
occurred, but the dose-assignment rule does not.

## Mixture rule: reproduce the implemented child-specific scores

After known outcomes, let a and b be the Beta parameters. For latent patient i,
parent component k is Beta(a+k,b+i-1-k), with weight w_(i-1,k). Denote its
parameters by A and B. For child z, its mean is theta_z=(A+z)/(A+B+1).
With remaining-window fraction q and prediction likelihoods l1 and l0:

    score_1 = l1*q*theta_1 / (l1*q*theta_1 + l0*(1-theta_1))
    score_0 = l0*(1-theta_0) / (l1*q*theta_0 + l0*(1-theta_0))
    child raw weight = parent weight * score_z

Normalize across both children of all current components; merge components
with identical Beta parameters. Positive prediction: (l1,l0)=(s,1-c);
negative prediction: (1-s,c); uninformative: (1,1). For untreated pending
follow-up, q=1-t/window under a uniform working model. For intervention patients
without observed DLT, q=1 even after the window. The final mixture has at most
m+1 components for m latent patients at a dose.

**This is the historical sequential child-specific plug-in rule. It is not
exact Bayesian marginal-evidence weighting and can depend on processing order.**
In particular, child-based splitting without informative evidence need not
preserve the parent distribution. The local document "Parent versus child
posterior branch weights" discusses a different, parent-evidence construction;
that construction has not silently replaced the October 31 simulation engine.
The figure and beta_mixture() document the implemented rule, not that alternative.
Predictor characteristics and rho are fixed simulation inputs; their uncertainty
is not estimated. A fully specified joint intervention/observation likelihood
would be needed to claim exact Bayesian inference under this observation process.

## Allocation, safety, and final selection

BOLD uses the mixture upper-tail probability Pr(p_j > target), isotonic
regression, and the closest tail probability to tau among the current dose and
adjacent levels. It starts at start_dose and escalates by at most one level.
At least three treated patients are required for overdose elimination. A dose
and higher levels are excluded when its overdose probability exceeds the
cutoff; exclusion of dose 1 stops the trial. Early stopping applies when the
current dose has early_n patients and is recommended again. Final selection
uses isotonic mixture means nearest the DLT target among treated, nonexcluded
levels adjacent to the current dose, following the saved engine. Cohort size
is fixed at three. This implementation should not be presumed identical to
every published BOLD or BOIN implementation.

BOIN uses interval boundaries based on 0.6*target and 1.4*target, a Beta(1,1)
safety prior, and isotonic final selection. TITE-BOIN uses the saved
pending-follow-up estimator. BOLD priors do not alter the BOIN comparator.

## Calendar assumptions

Uniform arrival gaps are U(0,2/rate); exponential gaps have mean 1/rate. A
cohort's patients arrive sequentially at its assigned dose. Untreated DLT times
can be uniform, calibrated Weibull, or deliberately late. The Weibull has the
specified unconditional DLT probability at window and half of within-window
DLTs by window/2. The working timing model is always uniform. No-event
follow-up resolves at window; an observed DLT resolves at its event time.

Complete-follow-up BOLD and BOIN wait for current outcomes before assigning a
new cohort dose. TITE modes ordinarily suspend when more than pending_limit of
patients at the current dose remain pending; optional min_two requires two
resolved outcomes before escalation. Safety de-escalation may proceed without
waiting. **Combined mode also waits for all intervened patients to resolve
before the next dose decision**, implementing the intervention-failure protocol.
Patients are not lost during suspension; a new arrival gap follows the next
enrollment. Duration extends through the last patient's outcome resolution.

## Outcomes

Accuracy is 100 times the fraction of replicate trials selecting the true MTD.
An explicit MTD=0 defines success as recommending no dose. MTD=K+1 defines
success as selecting K. Otherwise the default MTD is the dose closest to target.
Mean N and actual DLTs are replicate averages. Above MTD % is patient-weighted:
100*sum(patients above MTD)/sum(enrolled patients). All patients count as above
MTD in an all-too-toxic scenario. Duration is reported for TITE modes; it is
retained in raw complete-follow-up trials for optional calendar comparisons.
Monte Carlo SEs accompany all summaries. Scenario summaries are not silently
averaged over MTD categories; the analyst must specify the desired estimand.

## Provenance and limits

The native engine derives from Intervention failure 20261031/engine.cpp.
The default prior and calendar settings preserve that engine. Additions expose
dose-specific priors, starting dose, pending threshold, and exponential accrual;
zero-probability children now receive zero weight instead of dividing by zero.
In combined mode, ineligible patients provide no prediction label to pending
updates; this fixes the historical engine's use of generated labels for such
patients when eligibility was below one. Default eligibility=1 is unchanged.
The original October 31 engine and default settings are regression-tested.
This is a research release, not a validated clinical decision-support device.
