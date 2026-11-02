test_that("gamma is configurable without breaking the previous name", {
  expect_equal(bold_config()$overdose_cutoff,.90)
  expect_equal(bold_config(gamma=.80)$overdose_cutoff,.80)
  expect_equal(bold_config(tau=.45)$tau,.45)
  expect_equal(bold_config(overdose_cutoff=.80)$overdose_cutoff,.80)
  expect_error(bold_config(gamma=.80,overdose_cutoff=.90),"must agree")
  expect_error(bold_config(gamma=1),"probability")
  cfg <- bold_config(gamma=.80)
  expect_equal(do.call(bold_config,unclass(cfg)),cfg)
})

test_that("all four designs plan doses and only intervention designs intervene", {
  h <- data.frame(dose=rep(1,3),dlt=rep(0,3),followup=rep(3,3))
  for (m in c("bold","prediction","tite","combined")) {
    p <- plan_bold(h,4,1,m,prediction=c(1L,0L,0L))
    expect_equal(p$dose_recommendation$next_dose,2L)
    expect_equal(p$intervention_plan$dose,rep(2L,3))
    expect_true(p$ready)
    expect_equal(p$intervention_plan$intervene,
                 c(m %in% c("prediction","combined"),FALSE,FALSE))
  }
  expect_true(plan_bold(ndoses=4,method="bold")$ready)
  expect_false(plan_bold(ndoses=4,method="prediction")$ready)
  expect_true(plan_bold(ndoses=4,method="prediction",eligible=FALSE)$ready)
  expect_error(plan_bold(ndoses=4,prediction=1),"three")
})

test_that("waiting or stopping never schedules an intervention", {
  h <- data.frame(dose=rep(1,3),dlt=rep(0,3),followup=rep(1,3))
  for (m in c("bold","prediction","tite","combined")) {
    p <- plan_bold(h,4,1,m,prediction=rep(1L,3))
    expect_false(p$ready)
    expect_true(all(is.na(p$intervention_plan$dose)))
    expect_false(any(p$intervention_plan$intervene))
    p <- plan_bold(ndoses=4,method=m,excluded=1,prediction=rep(1L,3))
    expect_equal(p$dose_recommendation$action,"stop")
    expect_false(p$ready)
    expect_false(any(p$intervention_plan$intervene))
  }
})
