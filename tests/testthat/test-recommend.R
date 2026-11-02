test_that("Original BOLD supports starting, escalation, de-escalation and stopping", {
  cfg <- bold_config(early_n=Inf)
  expect_equal(recommend_dose(ndoses=4,config=cfg)$next_dose,1L)
  h <- data.frame(dose=rep(1,3),dlt=rep(0,3),followup=rep(3,3))
  d <- recommend_dose(h,4,1,config=cfg)
  expect_equal(d$action,"escalate")
  expect_equal(d$next_dose,2L)
  h$dose <- 2;h$dlt <- 1
  d <- recommend_dose(h,4,2,config=cfg)
  expect_equal(d$action,"deescalate")
  expect_equal(d$next_dose,1L)
  expect_true(all(2:4 %in% d$excluded_doses))
  h$dose <- 1
  d <- recommend_dose(h,4,1,config=cfg)
  expect_equal(d$action,"stop")
  expect_true(is.na(d$next_dose))
})

test_that("observed failures are known events and prevention branches persist", {
  h <- data.frame(dose=rep(1,3),dlt=rep(0,3),followup=rep(3,3),
    prediction=rep(1L,3),eligible=TRUE,intervened=TRUE)
  d <- recommend_dose(h,4,1,"prediction")
  p <- beta_mixture(.75,2.25,q=rep(1,3),labels=rep(1L,3))
  expect_equal(d$dose_summary$mean_untreated_risk[1],p$mean)
  expect_equal(d$dose_summary$overdose_probability[1],p$overdose_probability)
  expect_equal(d$dose_summary$latent_intervention[1],3L)
  h$dlt <- c(1,0,0);h$intervened <- c(TRUE,FALSE,FALSE)
  d <- recommend_dose(h,4,1,"prediction")
  expect_equal(d$dose_summary$mean_untreated_risk[1],1.75/6)
  expect_equal(d$dose_summary$latent_intervention[1],0L)
  expect_equal(d$patient_status$status[1],"observed_dlt")
})

test_that("TITE uses partial follow-up and observed outcomes replace pending labels", {
  h <- data.frame(dose=rep(1,3),dlt=rep(0,3),followup=c(3,3,1),
    prediction=c(0L,0L,0L),eligible=TRUE,intervened=FALSE)
  t <- recommend_dose(h,4,1,"tite")
  expect_false(t$action=="wait")
  b <- recommend_dose(h,4,1,"bold")
  expect_equal(b$action,"wait")
  h$followup <- 1
  expect_equal(recommend_dose(h,4,1,"tite")$action,"wait")
  h$followup <- 3
  a <- recommend_dose(h,4,1,"combined")
  h$prediction <- 1L
  b <- recommend_dose(h,4,1,"combined")
  expect_equal(a$dose_summary,b$dose_summary)
  expect_equal(a$dose_summary$mean_untreated_risk[1],.75/6)
})

test_that("combined mode waits for intervened patients even below pending threshold", {
  h <- data.frame(dose=rep(1,3),dlt=rep(0,3),followup=c(1,3,3),
    prediction=c(1L,0L,0L),eligible=TRUE,intervened=c(TRUE,FALSE,FALSE))
  expect_equal(recommend_dose(h,4,1,"combined")$action,"wait")
  h$followup <- 3
  d <- recommend_dose(h,4,1,"combined")
  expect_false(d$action=="wait")
  expect_equal(d$dose_summary$latent_intervention[1],1L)
})

test_that("elimination, cohort boundaries, cap and invalid history are explicit", {
  h <- data.frame(dose=rep(1,3),dlt=rep(0,3),followup=rep(3,3))
  d <- recommend_dose(h,4,1,excluded=2L)
  expect_equal(d$action,"stay")
  expect_equal(d$excluded_doses,2:4)
  expect_equal(recommend_dose(h,4,1,config=bold_config(nmax=3))$action,"stop")
  expect_equal(recommend_dose(h,4,1,config=bold_config(early_n=3),excluded=2)$action,"stop")
  expect_error(recommend_dose(h,4))
  expect_error(recommend_dose(h[1:2,],4,1))
  expect_error(recommend_dose(h,4,2))
  h$dlt[1] <- NA
  expect_error(recommend_dose(h,4,1))
  h$dlt[1] <- 0;h$intervened <- TRUE
  expect_error(recommend_dose(h,4,1,"tite"))
  expect_error(recommend_dose(h,4,1,"prediction"))
})

test_that("intervention plans are based on supplied predictions, not hidden outcomes", {
  p <- plan_intervention(c(1L,0L,NA_integer_,NA_integer_),c(TRUE,TRUE,FALSE,TRUE))
  expect_equal(p$intervene,c(TRUE,FALSE,FALSE,NA))
  expect_equal(p$status,c("positive","negative","not_eligible","await_prediction"))
  expect_error(plan_intervention(2))
})

test_that("BOIN comparator recommendations use the matched allocation rules", {
  h <- data.frame(dose=rep(1,3),dlt=rep(0,3),followup=rep(3,3))
  for (m in c("boin","tite_boin")) {
    expect_equal(recommend_dose(h,4,1,m)$next_dose,2L)
    hh <- h;hh$dlt <- 1
    expect_equal(recommend_dose(hh,4,1,m)$action,"stop")
  }
  set.seed(42);saved <- .Random.seed
  recommend_dose(h,4,1)
  expect_identical(.Random.seed,saved)
})
