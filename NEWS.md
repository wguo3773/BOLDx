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
