## Test Environment

* macOS Sequoia 15.5, Apple arm64
* R 4.5.1
* Apple clang 17.0.0
* TeX Live 2026 (PDF manual compilation enabled)

## R CMD Check Results

`R CMD check --as-cran BOLDx_0.3.0.tar.gz`

0 errors | 0 warnings | 2 notes

* This is a new submission.
* Local HTML validation was skipped because a recent HTML Tidy and the optional
  V8 math-rendering tool were unavailable. Rd checks and PDF manual compilation
  passed. These are local validation-tool notes, not disabled package checks.

All 149 test expectations passed. The six native simulation modes reproduce the
saved October 31 engine at its original default settings. Panel SD is calculated
across scenario-level estimates, with 100 trial replicates within each scenario
in the displayed 50-scenario example. The patient-history interface uses the
same native allocation routine as the simulator.

## Submission Metadata

Maintainer: Wanru Guo <wanru.guo@som.umaryland.edu>.

Windows/Linux and current-R-release checks remain to be performed before final
submission. The archive and manual have not been submitted to CRAN.
