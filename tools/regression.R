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
curves <- list(c(.4,.5,.6,.7),c(.05,.1,.15,.2),c(.1,.25,.4,.6),c(.01,.1,.2,.25))
checked <- 0L
for (mode in c(0L,1L,3L,4L,5L,6L)) for (cap in c(21L,30L))
  for (tau in c(.45,.50)) for (i in seq_along(curves)) {
    a <- list(truth=curves[[i]],phi=.25,nsim=30L,mode=mode,seed=20261101L+i,
      nmax=cap,window=3,rate=2,gamma=.9,tau=tau,sensitivity=.7,specificity=.9,
      early_n=12L,min_two=i%%2L==0L,eligibility=1,efficacy=if(cap==21L) .8 else 1)
    stopifnot(identical(do.call(old$calendar_combined,a),do.call(BOLDx:::calendar_combined,a)))
    checked <- checked+1L
  }
results <- c(results,paste(checked,"additional settings x 30 trials: identical to October 31 engine"),
  "Additional settings cover four toxicity curves, caps 21/30, tau .45/.50,",
  "all six modes, with/without the two-completed guard, and rho .8/1.")
writeLines(results,"docs/regression-results.txt")
cat(paste(results,collapse="\n"),"\n")
