library(BOLDx)
truth <- c(.10,.25,.45,.60)
cfg <- bold_config(nmax=21, target=.25, tau=.50,
  sensitivity=.90, specificity=.98, eligibility=1, rho=.80)
result <- simulate_bold(truth,
  method=c("prediction","tite","combined","bold","boin","tite_boin"),
  config=cfg, nsim=100, seed=20261031)
print(result)
write.csv(result$summary, "BOLDx-summary.csv", row.names=FALSE)
write.csv(result$trials, "BOLDx-trials.csv", row.names=FALSE)

sc <- read.csv(system.file("extdata","example-panel.csv",package="BOLDx"))
panel <- simulate_scenarios(as.matrix(sc[paste0("p",1:4)]),mtd=sc$mtd,
  method=c("prediction","tite","combined","bold","boin","tite_boin"),
  config=cfg,nsim=100,seed=20261031)
print(panel$aggregate)
write.csv(panel$aggregate,"BOLDx-panel-mean-sd.csv",row.names=FALSE)
