# Default native simulation must reproduce the October 31 research engine.
# Run with: Rscript tools/regression.R /path/to/original/engine.cpp
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Supply the original October 31 engine.cpp")
old <- new.env(parent=globalenv())
Rcpp::sourceCpp(args[1],env=old)
pkgload::load_all(".",quiet=TRUE)
results <- character()
for (mode in c(0L,1L,3L,4L,5L,6L)) {
  args <- list(truth=c(.1,.25,.45,.6),phi=.25,nsim=100L,mode=mode,
    seed=20261031L,nmax=21L,window=3,rate=2,gamma=.9,tau=.5,
    sensitivity=.9,specificity=.98,early_n=12L,min_two=FALSE,
    eligibility=1,efficacy=.8)
  before <- do.call(old$calendar_combined,args)
  after <- do.call(BOLDx:::calendar_combined,args)
  stopifnot(identical(before,after))
  results <- c(results,paste("Mode",mode,": 100 seeded trials identical to October 31 engine"))
}
writeLines(results,"docs/regression-results.txt")
cat(paste(results,collapse="\n"),"\n")
