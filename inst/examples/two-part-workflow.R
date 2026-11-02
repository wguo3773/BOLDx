library(BOLDx)

# PART 1: observed-history dose and intervention recommendations.
cfg <- bold_config(target=.25, nmax=21, gamma=.90, tau=.50)
h <- data.frame(dose=c(1,1,1), dlt=c(0,0,0), followup=c(3,3,3))
p <- plan_bold(h,4,1,method="prediction",config=cfg,
               prediction=c(1L,0L,0L),eligible=TRUE)
print(p$dose_recommendation)
print(p$intervention_plan)
print(p$ready)

# PART 2A: specified risk curve, 1000 trials per method.
methods <- c("bold","prediction","tite","combined")
single <- simulate_bold(p.true=c(.05,.15,.30,.45,.60),
  config=bold_config(target=.30),method=methods,nsim=1000,seed=42,mtd=3)
print(single$summary)
print(single$mcse)

# PART 2B: 50 random scenarios, 1000 trials per scenario per method.
panel <- simulate_random_scenarios(ndoses=4,mtd=2,lower_gap=.10,upper_gap=.15,
  risk_bounds=c(.01,.70),nscenarios=50,nsim=1000,
  config=cfg,method=methods,seed=42)
print(head(panel$scenarios))
print(panel$aggregate)

# PART 2C: a matrix of user-specified scenarios.
curves <- rbind(c(.10,.25,.45,.60),c(.05,.10,.25,.45),c(.40,.50,.60,.70))
custom <- simulate_scenarios(curves,mtd=c(2,3,0),config=cfg,
                             method=methods,nsim=1000,seed=42)
print(custom$summary)
print(custom$aggregate)
