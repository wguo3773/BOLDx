# BOLDx

### Prediction-assisted intervention and time-to-event extensions of BOLD

**An R package for simulation-based phase I dose-finding research.**

BOLDx supports three extensions of the Bayesian Ordered Lattice Design:
**prediction-guided preventive intervention**, **delayed toxicity outcomes
(TITE-BOLD)**, and **their combination**. Original BOLD, BOIN, and TITE-BOIN
implementations are included for matched comparisons.
The package also evaluates next-cohort dose recommendations from entered patient
histories, using the same allocation rules as its simulations.

[Installation](#installation) | [Quick start](#quick-start) |
[Patient-history workflow](#patient-history-workflow) |
[Inputs](#design-and-prediction-inputs) | [Methods](docs/methods.md) |
[Reference manual](docs/BOLDx-reference-manual.pdf) |
[Research background](docs/research-background.md) |
[Related Shiny app](https://wguo3.shinyapps.io/bold_boin_phase1_app/)

![Prediction-assisted BOLD: persistent untreated-outcome branches](docs/prediction-branching-schematic.png)

![TITE-BOLD with prediction: observed-outcome replacement and persistent intervention counterfactuals](docs/tite-prediction-branching-schematic.png)

For untreated patients, observed outcomes replace pending prediction/follow-up
contributions. For intervened patients, an observed failure is a known DLT;
no observed DLT leaves the untreated counterfactual unresolved. This distinction
applies in both prediction-assisted modes.

## Background

In phase I oncology trials, preventing severe toxicity and accommodating delayed
outcomes may improve patient safety and dose finding. BOLDx extends BOLD to
prediction-guided intervention, delayed outcomes, or both, representing
unobserved untreated toxicity with weighted Beta components. Prediction
sensitivity/specificity and intervention effectiveness are separate inputs:
an intervention can fail even after a correct positive prediction.

The accompanying study examined four doses, a target DLT probability of 0.25,
patient caps of 21 and 30, upper toxicity gaps of 0.15 and 0.20, and prediction
sensitivity/specificity of 0.70/0.90 and 0.90/0.98. Its imperfect-intervention
scenario used 80% effectiveness. The supplied research summary reports lower
actual DLT burden with prediction-assisted intervention and shorter trial
duration with delayed-outcome BOLD relative to complete-follow-up BOLD.
The [full supplied background and numerical results](docs/research-background.md)
are retained separately from the package's reproducible examples.

## Installation

R >= 4.1 and a C++11 compiler are required. Install Rtools on Windows or Xcode
Command Line Tools on macOS if a compiler is not already available.

```r
install.packages("remotes")
remotes::install_github("wguo3773/BOLDx", build_vignettes = FALSE)
library(BOLDx)
```

Alternatively, download the [R source tarball](docs/BOLDx_0.2.0.tar.gz) and run:

```r
install.packages("Rcpp")
install.packages("BOLDx_0.2.0.tar.gz", repos = NULL, type = "source")
```

## Patient-history workflow

In addition to simulating trials, BOLDx can evaluate the next cohort's dose from
entered patient histories. `recommend_dose()` returns **start, escalate, stay,
deescalate, wait, or stop**, together with the next dose, its reason, and
dose-level summaries. It uses the same allocation routine as the simulator;
no assumed true toxicity curve is needed.

### Original BOLD: up, down, or stay

```r
cfg <- bold_config(target=.25, nmax=21, tau=.50)
h <- data.frame(
  dose     = c(1, 1, 1),
  dlt      = c(0, 0, 0),
  followup = c(3, 3, 3)
)
d <- recommend_dose(h, ndoses=4, current_dose=1, method="bold", config=cfg)
print(d)
# Method: bold | Action: escalate | Next dose: 2
# Next cohort: dose 2 (escalate).
d$next_dose
d$dose_summary
```

Enter the next cohort's actual dose and outcomes, then call the function again.
Preserve eliminated doses with `excluded = d$excluded_doses`. `current_dose`
must be the dose of the most recently enrolled cohort, not a proposed future
dose. Evaluation is supported between cohorts of three patients.

### Prediction-guided intervention

Supply predictions from an external classifier; sensitivity and specificity
describe that classifier's performance and do not predict an individual
patient by themselves.

```r
plan_intervention(c(1L, 0L, 0L), eligible=TRUE)
#   prediction eligible intervene   status
# 1          1     TRUE      TRUE positive
# 2          0     TRUE     FALSE negative
# 3          0     TRUE     FALSE negative

# Record interventions actually performed, not just the proposed flags.
hp <- transform(h,
  prediction = c(1L, 0L, 0L), eligible=TRUE,
  intervened = c(TRUE, FALSE, FALSE)
)
recommend_dose(hp, 4, 1, method="prediction", config=cfg)
```

An observed intervention failure is entered as `dlt=1`. If no DLT is observed
after intervention, enter the observed record; **do not invent the untreated
outcome or mark the intervention as proven successful**. The function retains
that patient's latent untreated-outcome branches automatically.

### TITE-BOLD: update partial follow-up

```r
ht <- h
ht$followup <- c(3, 3, 1)  # one patient is still within the 3-month window
recommend_dose(ht, 4, 1, method="tite", config=cfg)
# Method: tite | Action: escalate | Next dose: 2
```

Update `followup` as time elapses and set `dlt=1` if a DLT is observed.
`dlt=0` means no DLT observed **to date**, so it does not imply complete
follow-up when `followup < cfg$window`. The function returns `wait` when the
pending-outcome rule requires it.

### Combined prediction and TITE-BOLD

```r
hc <- transform(h, prediction=c(0L, 1L, 0L), eligible=TRUE,
                intervened=c(FALSE, TRUE, FALSE))
hc$followup <- c(3, 1, 1)
recommend_dose(hc, 4, 1, method="combined", config=cfg)
# Method: combined | Action: wait
# Combined protocol: await intervention-patient outcome resolution.

# Later, this example's intervened patient completes follow-up without DLT.
hc$followup[2] <- 3
recommend_dose(hc, 4, 1, method="combined", config=cfg)
# Method: combined | Action: stay | Next dose: 1
```

Combined mode waits for intervention-patient resolution before another dose
decision, while using partial follow-up for other patients. Observed untreated
outcomes replace their pending prediction contribution; an intervened patient
without an observed DLT retains the untreated-outcome mixture.

The recommendations above are outputs for these example histories, not rules
to escalate whenever three patients have no observed DLT. The result depends
on the complete entered history and configured settings. `wait` and `stop`
have `next_dose=NA`. This function evaluates next-cohort allocation, not final
MTD selection. `rho` and the simulation's eligibility probability are used to
generate simulated histories; entered histories use the actual eligibility,
predictions, interventions, outcomes, and follow-up supplied by the user.

## Quick start

```r
library(BOLDx)

cfg <- bold_config(
  target = 0.25, nmax = 21, tau = 0.50,
  overdose_cutoff = 0.90, early_n = 12,
  sensitivity = 0.90, specificity = 0.98,
  eligibility = 1.00, rho = 0.80,
  window = 3, accrual_rate = 2
)

result <- simulate_bold(
  truth = c(0.10, 0.25, 0.45, 0.60),
  method = c("prediction", "tite", "combined", "bold", "boin", "tite_boin"),
  config = cfg, nsim = 100, seed = 20261031,
  mtd = 2
)

result$summary  # operating characteristics for this scenario
result$trials   # replicate-level outcomes and dose-specific patient/DLT counts
write.csv(result$summary, "BOLDx-summary.csv", row.names = FALSE)
```

### Single-scenario output

The quick-start code above produces the following summary for 100 simulated
trials per method (seed `20261031`, true MTD = dose 2). Values are rounded for
display. Between-scenario SD requires a panel of scenarios, as shown below.

| Method | Accuracy (%) | Mean N | Actual DLTs per trial | Above MTD (%) | Mean duration (months) |
|:--|--:|--:|--:|--:|--:|
| BOLD prediction | 59.0 | 19.65 | 1.56 | 36.18 | - |
| TITE-BOLD | 58.0 | 18.84 | 4.84 | 24.84 | 16.83 |
| TITE-BOLD + prediction | 56.0 | 19.59 | 1.50 | 32.16 | 21.62 |
| Original BOLD | 58.0 | 19.77 | 5.39 | 26.86 | - |
| BOIN | 55.0 | 19.05 | 4.26 | 17.17 | - |
| TITE-BOIN | 48.0 | 18.54 | 4.11 | 16.02 | 15.13 |

The dash denotes duration not displayed in the complete-follow-up summaries.
These are results for one example scenario, not averages over the study's
scenario panel. Download the [full output CSV](docs/example-summary.csv).

For example, inspect prediction-assisted BOLD alone:

```r
subset(result$summary, method == "prediction",
       select = c(accuracy_pct, mean_n, actual_dlts, above_mtd_pct))
#   accuracy_pct mean_n actual_dlts above_mtd_pct
# 1           59  19.65        1.56      36.18321
```

### Across 50 random scenarios: mean (SD)

Each of these 50 saved random scenarios uses 100 simulated trials per method.
This example uses MTD = dose 2, upper toxicity gap = 0.15, patient cap = 21,
sensitivity/specificity = 0.90/0.98, and intervention effectiveness = 0.80.

```r
sc <- read.csv(system.file("extdata", "example-panel.csv", package = "BOLDx"))
panel <- simulate_scenarios(
  scenarios = as.matrix(sc[paste0("p", 1:4)]), mtd = sc$mtd,
  method = c("prediction", "tite", "combined", "bold", "boin", "tite_boin"),
  config = cfg, nsim = 100, seed = 20261031
)
panel$aggregate  # means and SDs across the 50 scenario-level estimates
```

| Method | Accuracy (%) | Mean N | Actual DLTs per trial | Above MTD (%) | Mean duration (months) |
|:--|--:|--:|--:|--:|--:|
| BOLD prediction | 47.26 (6.11) | 19.34 (0.29) | 1.59 (0.11) | 38.55 (2.97) | - |
| TITE-BOLD | 45.64 (4.75) | 18.33 (0.43) | 4.73 (0.16) | 25.88 (2.51) | 16.33 (0.39) |
| TITE-BOLD + prediction | 46.72 (5.94) | 19.15 (0.28) | 1.50 (0.10) | 33.88 (2.82) | 20.77 (0.40) |
| Original BOLD | 49.14 (4.78) | 19.36 (0.28) | 5.17 (0.12) | 28.87 (2.37) | - |
| BOIN | 40.56 (4.95) | 17.62 (0.50) | 4.11 (0.13) | 17.19 (2.06) | - |
| TITE-BOIN | 39.54 (3.40) | 17.73 (0.37) | 4.09 (0.13) | 16.66 (2.49) | 14.44 (0.36) |

Values are **mean (SD) across scenarios**. For each outcome, SD is computed
from 50 scenario-level estimates, not from the 5,000 individual simulated
trials. It describes between-scenario dispersion, including simulation noise
from the finite number of trials within each scenario.

[Scenario-level estimates](docs/example-panel-estimates.csv) |
[Mean/SD summary](docs/example-panel-mean-sd.csv) |
[The 50 scenario definitions](inst/extdata/example-panel.csv).

`nsim = 100` is a runnable example, not a recommendation for precise design
comparison. Increase it to reduce Monte Carlo error. The caller's RNG state is
restored after simulation. Methods use the same starting seed within a scenario;
they need not consume identical random draws after protocol differences arise.

### Three extension modes

| `method` | Prediction/intervention | Delayed-outcome decisions |
|:--|:--:|:--:|
| `"prediction"` | Yes | No; complete follow-up |
| `"tite"` | No | Yes |
| `"combined"` | Yes | Yes, with intervention-resolution guard |
| `"bold"` | No | No; Original BOLD comparator |
| `"boin"` | No | No; BOIN comparator |
| `"tite_boin"` | No | Yes; TITE-BOIN comparator |

### Requested outputs

| Summary column | Definition |
|:--|:--|
| `accuracy_pct` | Trials selecting the true MTD (%) |
| `mean_n` | Mean enrolled sample size |
| `actual_dlts` | Mean observed, post-intervention DLTs per trial |
| `above_mtd_pct` | Patients treated above the true MTD / all enrolled patients (%) |
| `duration_months` | Mean trial duration, for TITE and combined modes and TITE-BOIN |

Panel summaries include an SD across scenario-level estimates for each outcome.
The SD is undefined for a group containing only one scenario. Complete-follow-up
summary duration is `NA`, not zero; its calendar duration is available in
`result$trials$duration` for optional comparisons. Raw output also includes
prevented DLTs, intervention failures, latent DLTs, eligibility, interventions,
and allocation counts. `mtd = 0` includes all-too-toxic scenarios; `mtd = K+1`
represents all doses below target. These evaluation labels never enter decisions.

## Design and prediction inputs

```r
cfg <- bold_config(
  target = .25, nmax = 30, tau = .45,
  overdose_cutoff = .90, early_n = Inf,
  sensitivity = .70, specificity = .90,
  eligibility = 1, rho = .80,
  prior_mean = c(.15, .20, .25, .35),
  prior_pess = c(3, 3, 3, 3),
  start_dose = 1, window = 3, accrual_rate = 2,
  arrival = "uniform", timing = "weibull",
  pending_limit = .50, min_two = FALSE
)
```

- `truth`: 2-10 ordered untreated dose-risk probabilities; not prior means.
- `eligibility`: fraction eligible for prediction (e.g., nonmissing predictors).
- `rho`: intervention success probability; `.80` means 20% failure, `1` means
  perfect effectiveness. It is not sensitivity or specificity.
- Prior at dose j: `Beta(prior_mean[j] * prior_pess[j],
  (1 - prior_mean[j]) * prior_pess[j])`. Scalars apply to every dose.
- `tau`: BOLD's decision-tail target, separate from the DLT target.
- Cohorts contain three patients; `nmax` must be a multiple of three.
- `arrival`: `"uniform"` or `"exponential"`; `timing`: `"weibull"`, `"uniform"`,
  or the `"late"` stress test. Follow-up uses a uniform working timing model.

An observed intervention failure is a regular DLT update, not an additional
branch. Only intervention patients without observed DLT branch over their
unknown untreated status. Combined mode waits for intervened patients to resolve
before the next dose decision. False-positive interventions have no simulated
toxicity harm or other cost in this release.

### Multiple scenarios

```r
scenarios <- rbind(
  c(.10, .25, .45, .60),
  c(.05, .10, .25, .45),
  c(.40, .50, .60, .70)
)
panel <- simulate_scenarios(
  scenarios, mtd = c(2, 3, 0), config = cfg,
  method = c("prediction", "bold", "boin"), nsim = 100, seed = 1
)
panel$summary
panel$aggregate
```

`panel$summary` retains each scenario; `panel$aggregate` reports equally weighted
means and SDs separately by method and MTD category, including all-too-toxic
cases. To keep additional settings separate, add their labels to the summary
and use `summarise_scenarios(panel$summary, group_by = c("method", "mtd", "delta"))`.

### Inspect the mixture

```r
beta_mixture(
  a = .75, b = 2.25,
  q = c(1, .5), labels = c(1L, -1L),
  sensitivity = .90, specificity = .98, target = .25
)
```

Download the branching schematics as vector PDFs:
[prediction-assisted BOLD](docs/prediction-branching-schematic.pdf) and
[TITE-BOLD with prediction](docs/tite-prediction-branching-schematic.pdf).
The [original shared-update PDF](docs/beta-mixture-formulas.pdf) and
[SVG](docs/beta-mixture-formulas.svg) remain available. All three preserve the
implemented formulas, including global normalization and merged branches.

## Documentation and reproducibility

- [R reference manual (PDF)](docs/BOLDx-reference-manual.pdf), generated by
  `R CMD Rd2pdf` from the included `man/*.Rd` help pages.
- [Methods, assumptions, outcome definitions, and provenance](docs/methods.md).
- [Full research background supplied for this release](docs/research-background.md).
- [Runnable quick-start script](inst/examples/quick-start.R).
- Tests cover all modes, accounting identities, random-state restoration,
  limiting cases, scenario boundary labels, priors, and input validation.
- [Regression check](docs/regression-results.txt): 100 seeded trials in each
  of the six modes exactly match the October 31 engine at the default settings.
- [Single-scenario output](docs/example-summary.csv) and
  [50-scenario mean/SD output](docs/example-panel-mean-sd.csv) contain results
  generated by the included example code.

```r
install.packages(c("roxygen2", "testthat", "devtools"))
devtools::test()
```

From the repository parent directory:

```sh
R CMD build BOLDx
R CMD check BOLDx_0.2.0.tar.gz
```

To regenerate the reference manual with a LaTeX installation available:

```sh
R CMD Rd2pdf --no-preview --force --output=BOLDx-reference-manual.pdf BOLDx
```

For CRAN submission, use the [R source package](docs/BOLDx_0.2.0.tar.gz).
It includes the code, Rd documentation, examples, and tests. The PDF is the
readable reference generated from those documentation sources.

Maintainer: Wanru Guo, [wanru.guo@som.umaryland.edu](mailto:wanru.guo@som.umaryland.edu).

## Citation and acknowledgments

Wang G-M and Tatsuoka C (2026). *Bayesian Ordered Lattice Design for Phase I
Clinical Trials*. Statistics in Medicine. [doi:10.1002/sim.70456](https://doi.org/10.1002/sim.70456).

This package extends the research simulation engine maintained by Wanru Guo.
Please also cite the BOLDx software version used. Publication details for the
prediction/TITE extension manuscript will be added when available.

## License

GPL (>= 3). Research software provided without warranty.
