test_that("p.true is a validated alias and defaults use 1000 trials", {
  a <- simulate_bold(p.true=c(.1,.25,.45),method="bold",nsim=5)
  b <- simulate_bold(c(.1,.25,.45),method="bold",nsim=5)
  expect_identical(a,b)
  expect_error(simulate_bold(c(.1,.2),p.true=c(.1,.3)),"must agree")
  expect_equal(formals(simulate_bold)$nsim,1000L)
  expect_equal(formals(simulate_scenarios)$nsim,1000)
  expect_equal(a$mcse$actual_dlts_se,sd(a$trials$dlt)/sqrt(5))
  expect_equal(a$mcse$accuracy_se,100*sd(a$trials$correct)/sqrt(5))
})

test_that("gap-based curves are ordered, reproducible and respect category", {
  set.seed(19); old <- .Random.seed
  a <- generate_scenarios(mtd=2,seed=42)
  expect_identical(.Random.seed,old)
  expect_equal(dim(a),c(50L,4L))
  expect_true(all(a[,1]==.15 & a[,2]==.25 & a[,3]==.40))
  expect_true(all(a[,4]>.40 & a[,4]<.70))
  expect_identical(a,generate_scenarios(mtd=2,seed=42))
  for (j in 0:5) {
    x <- generate_scenarios(mtd=j,nscenarios=3)
    expect_true(all(apply(x,1,function(p) all(diff(p)>=0))))
    expect_true(all(x>0 & x<1))
    if (j==0) expect_true(all(x>.25))
    else if (j==5) expect_true(all(x<.25))
    else expect_true(all(x[,j]==.25))
  }
  expect_equal(dim(generate_scenarios(nscenarios=1)),c(1L,4L))
  expect_error(generate_scenarios(lower_gap=.3),"adjacent")
  expect_error(generate_scenarios(upper_gap=.6),"adjacent")
  expect_error(generate_scenarios(mtd=6),"Invalid")
  expect_error(generate_scenarios(risk_bounds=c(.3,.7)),"contain target")
  expect_error(generate_scenarios(seed=-1),"Invalid")
})

test_that("random panels keep settings and two uncertainty levels separate", {
  z <- simulate_random_scenarios(mtd=c(0,2),upper_gap=c(.15,.20),
    method=c("bold","prediction"),nscenarios=3,nsim=5,seed=42)
  expect_equal(nrow(z$scenarios),12L)
  expect_equal(nrow(z$summary),24L)
  expect_equal(nrow(z$mcse),24L)
  expect_equal(nrow(z$trials),120L)
  expect_equal(nrow(z$aggregate),8L)
  expect_true(all(z$aggregate$n_scenarios==3))
  expect_true(all(z$aggregate$replicates_min==5))
  expect_true(all(c("upper_gap","lower_gap","scenario") %in% names(z$mcse)))
  d <- subset(z$summary,method=="bold" & mtd==2 & upper_gap==.15)
  a <- subset(z$aggregate,method=="bold" & mtd==2 & upper_gap==.15)
  expect_equal(a$accuracy_pct,mean(d$accuracy_pct))
  expect_equal(a$accuracy_sd,sd(d$accuracy_pct))
  expect_identical(z,simulate_random_scenarios(mtd=c(0,2),upper_gap=c(.15,.20),
    method=c("bold","prediction"),nscenarios=3,nsim=5,seed=42))
  expect_error(simulate_random_scenarios(seed=.1),"Invalid")
})
