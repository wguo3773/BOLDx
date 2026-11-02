test_that("all modes conserve actual and latent events", {
  z <- simulate_bold(c(.1,.25,.45,.6), nsim = 30, seed = 42,
    method = c("prediction","tite","combined","bold","boin","tite_boin"))
  x <- z$trials
  expect_equal(rowSums(x[paste0("n",1:4)]), x$n)
  expect_equal(rowSums(x[paste0("y",1:4)]), x$dlt)
  expect_equal(x$dlt+x$prevented, x$latent_dlt)
  expect_true(all(x$n <= 21 & x$n %% 3 == 0))
  expect_true(all(x$intervention_failures <= x$dlt))
  expect_equal(x$mixture_patients, x$interventions-x$intervention_failures)
  expect_true(all(x$duration >= x$last_enrollment))
  expect_true(all(is.na(z$summary$duration_months[c(1,4,5)])))
})

test_that("failure and eligibility limiting cases reduce to original BOLD", {
  cfg <- bold_config(sensitivity=1, specificity=1, rho=0)
  z <- simulate_bold(c(.1,.25,.45,.6), c("prediction","bold"), cfg, nsim=40)
  a <- subset(z$trials, method=="prediction")
  b <- subset(z$trials, method=="bold")
  for (nm in c("selected","n","dlt","duration")) expect_equal(a[[nm]],b[[nm]])
  z <- simulate_bold(c(.1,.25,.45,.6), c("prediction","bold"),
    bold_config(eligibility=0), nsim=40)
  a <- subset(z$trials, method=="prediction")
  b <- subset(z$trials, method=="bold")
  expect_equal(a$selected,b$selected)
  expect_equal(a$dlt,b$dlt)
  z <- simulate_bold(c(.1,.25,.45,.6), c("combined","tite"),
    bold_config(eligibility=0), nsim=40)
  a <- subset(z$trials, method=="combined")
  b <- subset(z$trials, method=="tite")
  for (nm in c("selected","n","dlt","duration")) expect_equal(a[[nm]],b[[nm]])
})

test_that("perfect prevention yields no observed events", {
  z <- simulate_bold(c(.1,.25,.45,.6), "prediction",
    bold_config(sensitivity=1,specificity=1,rho=1), nsim=40)
  expect_true(all(z$trials$dlt==0))
  expect_true(all(z$trials$intervention_failures==0))
})

test_that("summaries use the supplied evaluation MTD", {
  z <- simulate_bold(c(.4,.5,.6,.7), "prediction", nsim=10, mtd=0)
  expect_equal(z$summary$above_mtd_pct,100)
  expect_equal(z$summary$accuracy_pct,100*mean(z$trials$selected==0))
  z <- simulate_bold(c(.05,.1,.15,.2), "tite", nsim=10, mtd=5)
  expect_equal(z$summary$above_mtd_pct,0)
  expect_equal(z$summary$accuracy_pct,100*mean(z$trials$selected==4))
})

test_that("random state is restored and runs reproduce", {
  set.seed(9); saved <- .Random.seed
  a <- simulate_bold(c(.1,.25,.45),nsim=4,seed=8)
  expect_identical(.Random.seed,saved)
  b <- simulate_bold(c(.1,.25,.45),nsim=4,seed=8)
  expect_identical(a,b)
})

test_that("all configurable timing and prior settings run", {
  z <- simulate_bold(c(.1,.25,.45), config=bold_config(
    prior_mean=c(.15,.25,.35),prior_pess=c(2,3,4),start_dose=2,
    arrival="exponential",timing="late",pending_limit=.4,min_two=TRUE),nsim=5)
  expect_equal(nrow(z$summary),3)
  expect_true(all(is.finite(z$trials$duration)))
})

test_that("invalid configurations are rejected before native calls", {
  expect_error(bold_config(nmax=22))
  expect_error(bold_config(target=.8))
  expect_error(bold_config(rho=NA_real_))
  expect_error(simulate_bold(c(.3,.2)))
  expect_error(simulate_bold(c(.1,.25),nsim=1))
  expect_error(simulate_bold(c(.1,.25),method="unknown"))
  expect_error(simulate_bold(c(.1,.25),config=bold_config(start_dose=3)))
  expect_error(simulate_bold(c(.1,.25),config=bold_config(prior_pess=1:3)))
})

test_that("scenario panels retain MTD labels", {
  z <- simulate_scenarios(rbind(c(.1,.25,.45), c(.4,.5,.6)),
    mtd=c(2,0),method="prediction",nsim=5)
  expect_equal(z$summary$mtd,c(2,0))
  expect_equal(z$summary$above_mtd_pct[2],100)
})

test_that("mixtures normalize and match a hand-calculated child update", {
  z <- beta_mixture(1,3,q=1,labels=1,sensitivity=.9,specificity=.98)
  t1 <- 2/5; t0 <- 1/5
  s1 <- .9*t1/(.9*t1+.02*(1-t1))
  s0 <- .02*(1-t0)/(.9*t0+.02*(1-t0))
  expect_equal(z$components$weight,c(s0,s1)/(s0+s1))
  expect_equal(sum(z$components$weight),1)
  expect_equal(beta_mixture(1,3)$mean,.25)
  expect_equal(beta_mixture(1,3,q=1,labels=1,sensitivity=1,specificity=1)$components$weight,c(0,1))
  expect_error(beta_mixture(1,3,q=2))
})
