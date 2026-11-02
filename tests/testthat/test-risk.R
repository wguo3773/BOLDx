test_that("fixed-dose actual risk follows the intervention success model", {
  cfg <- bold_config(nmax=3, sensitivity=.9, specificity=.98, rho=.8)
  z <- simulate_bold(c(.25,.25),"prediction",cfg,nsim=3000,seed=502)
  observed <- sum(z$trials$dlt)/sum(z$trials$n)
  theoretical <- .25*(1-.9*.8)
  expect_lt(abs(observed-theoretical),.02)
})
