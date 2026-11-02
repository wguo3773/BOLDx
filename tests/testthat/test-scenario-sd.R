test_that("scenario SD is computed from scenario estimates, not pooled trials", {
  x <- data.frame(method="prediction", mtd=2,scenario=1:2,nsim=100,
    accuracy_pct=c(50,60),mean_n=c(18,20),actual_dlts=c(1,3),
    above_mtd_pct=c(20,40),duration_months=NA_real_)
  z <- summarise_scenarios(x)
  expect_equal(z$n_scenarios,2L)
  expect_equal(z$accuracy_pct,55)
  expect_equal(z$accuracy_sd,stats::sd(c(50,60)))
  expect_equal(z$mean_n_sd,stats::sd(c(18,20)))
  expect_equal(z$actual_dlts_sd,stats::sd(c(1,3)))
  expect_equal(z$above_mtd_sd,stats::sd(c(20,40)))
  expect_true(is.na(z$duration_sd))
  expect_error(summarise_scenarios(rbind(x,x)))
  expect_true(is.na(summarise_scenarios(x[1,])$accuracy_sd))
})

test_that("scenario panel returns mean and SD without SE columns", {
  z <- simulate_scenarios(rbind(c(.1,.25,.45),c(.15,.25,.4)),mtd=c(2,2),
    method=c("prediction","tite"),nsim=5)
  expect_equal(nrow(z$aggregate),2L)
  expect_true(all(z$aggregate$n_scenarios==2L))
  expect_false(any(grepl("_se$",names(z$summary))))
  expect_false(any(grepl("_se$",names(z$aggregate))))
  expect_equal(summarise_scenarios(z),z$aggregate)
})
