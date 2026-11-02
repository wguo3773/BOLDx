# Research background and supplied study summary

The following background preserves the numerical results supplied by the study
author for the 80% intervention-success setting. These are **reported research
results**, not claims that the newly packaged example has independently rerun
and verified every scenario. The source abstract specifies 100 trials per method
and scenario; the saved October 31 run script defaults to 1,000. Reconcile that
replicate count and the aggregation/percentage definitions before using these
numbers as a formal package benchmark or manuscript result.

## Background

In phase I oncology trials, preventing severe toxicity and accommodating delayed
outcomes may improve patient safety and dose finding. Bayesian Ordered Lattice
Designs (BOLD) have been shown to improve identification of maximum tolerable
dose levels (MTDs) over popular models. Here, we extend BOLD to incorporate
prediction-guided intervention, delayed toxicity outcomes, or both, using
weighted Beta-mixture models to represent unobserved untreated toxicity.

## Simulation setting

Simulations included four doses, target DLT probability 0.25, BOLD decision
target tau 0.50, and patient caps of 21 or 30. For each MTD category, 50 random
scenarios were evaluated using upper toxicity gaps of 0.15 and 0.20. All patients
were eligible for prediction, with sensitivity/specificity 0.70/0.90 or 0.90/0.98.
The main imperfect-intervention setting used 80% effectiveness (20% failure);
100% effectiveness was an additional setting, not another 20%-failure setting.
Calendar simulations used uniform arrival gaps averaging two patients monthly,
three-month follow-up, Weibull DLT times, and a uniform working timing model.

For scenario-panel summaries, report mean (SD) across the 50 random scenarios,
with each scenario's operating characteristics estimated from 100 simulated
trials. SD describes dispersion across scenario-level estimates, not a
standard error or variability pooled across individual trials.

## Reported results at 80% effectiveness

Compared with BOIN, BOLD-prediction accuracy at sensitivity/specificity 0.90/0.98
increased by 4.0-6.0% and 3.2-5.3% overall, and 10.6-13.0% and 9.5-12.1%
excluding all-too-toxic scenarios, along with 60.6-61.0% and 62.3-62.8% overall
DLT reductions at caps 21 and 30, respectively. At sensitivity/specificity
0.70/0.90, accuracy increased by 6.8-8.8% and 6.2-8.2% compared with BOIN at caps
21 and 30, respectively, after excluding all-too-toxic scenarios, while mean
actual DLTs per trial decreased by 43.1-43.7% at cap 21 and 45.4-46.2% at cap 30.

For TITE-BOLD, accuracy increased by 3.5-4.6% at cap 21 and 4.1-5.5% at cap 30
across all categories versus TITE-BOIN. Excluding all-too-toxic scenarios, the
improvements were 6.6-8.0% and 7.2-8.9%, respectively. Compared with Original
BOLD, prediction-assisted BOLD reduced actual DLTs by 69.37-69.44% at
sensitivity/specificity 0.90/0.98 with cap 21, and by 54.0-54.1% at 0.70/0.90,
with little impact on accuracy.

TITE-BOLD with overdose control tau=0.45 achieved reported accuracy improvements
of 4.4% and 5.6% compared with TITE-BOIN, with similar actual DLTs and mean trial
duration at caps 21 and 30. It also had similar accuracy and actual DLTs to
Original BOLD while shortening mean trial duration by 30.6% and 33.5% at caps 21
and 30. Compared with TITE-BOIN, combined prediction and TITE-BOLD reduced DLTs
by 45.0-65.9% while having higher accuracy.

The supplied percentage changes are retained as reported. They should not be
interpreted as percentage-point differences or significance tests without
checking the original comparison calculations.

## Interpretation

Under imperfect intervention effectiveness, prediction-assisted BOLD reduced
actual DLT burden in the simulated settings while outperforming BOIN and
preserving Original BOLD accuracy. TITE-BOLD shortened duration relative to
complete-follow-up BOLD and improved MTD identification relative to TITE-BOIN.
Combining prediction with TITE-BOLD offered additional reductions in actual DLTs.
These findings motivate further methodological and real-world evaluation;
clinical effectiveness has not been established by these simulations.
