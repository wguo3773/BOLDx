# BOLDx 0.2.0

* Added `recommend_dose()` for patient-history next-cohort recommendations.
* Added `plan_intervention()` for external prediction labels and eligibility.
* Shared native allocation rules keep the history interface and simulator aligned.

# BOLDx 0.1.0

* Initial research release based on the October 31, 2026 intervention-failure engine.
* Prediction, TITE, combined, Original BOLD, BOIN, and TITE-BOIN simulation modes.
* Configurable dose-specific Beta priors, starting dose, arrival distribution,
  eligibility, prediction accuracy, and intervention effectiveness.
* Auditable trial-level output and Monte Carlo operating characteristics.
* Panel summaries report means and SDs across random scenarios; individual
  trial replicates are not used as the units for between-scenario SD.
# BOLDx 0.4.0

- Organized Part 1 (dose/intervention planning) and Part 2 (simulation) with
  input/output tables and reproducible numerical examples.
- Added p.true as an alias for truth and set simulation defaults to 1000 trials.
- Added generate_scenarios() and simulate_random_scenarios() with documented
  uniform-tail generation, adjacent-dose gaps and separate reporting groups.
- Retained within-scenario MCSE in panels, separate from between-scenario SD.

# BOLDx 0.3.0

- Added plan_bold() to return dose recommendations and incoming-patient
  intervention flags together for all four BOLD designs.
- Exposed gamma (default 0.90) as an alias for the configurable posterior
  overdose cutoff; tau remains configurable with default 0.50.
- Preserved simulation and dose-allocation calculations.

