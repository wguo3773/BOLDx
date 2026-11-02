pkgload::load_all(".",quiet=TRUE)
methods <- c("bold","prediction","tite","combined")
single <- simulate_bold(p.true=c(.05,.15,.30,.45,.60),method=methods,
  config=bold_config(target=.30),nsim=1000,seed=42,mtd=3)
random <- simulate_random_scenarios(mtd=2,upper_gap=.15,method=methods,
  config=bold_config(target=.25),nscenarios=50,nsim=1000,seed=42)
write.csv(single$summary,"docs/single-1000-summary.csv",row.names=FALSE)
write.csv(single$mcse,"docs/single-1000-mcse.csv",row.names=FALSE)
write.csv(random$summary,"docs/random-50x1000-estimates.csv",row.names=FALSE)
write.csv(random$mcse,"docs/random-50x1000-mcse.csv",row.names=FALSE)
write.csv(random$aggregate,"docs/random-50x1000-mean-sd.csv",row.names=FALSE)
write.csv(cbind(random$settings,random$scenarios),
          "docs/random-50x1000-scenarios.csv",row.names=FALSE)
saveRDS(list(single=single,random=random),"/tmp/BOLDx-workflow-examples.rds")
cat("SINGLE SUMMARY\n"); print(single$summary,digits=5)
cat("SINGLE MCSE\n"); print(single$mcse,digits=5)
cat("PANEL MEAN SD\n"); print(random$aggregate,digits=5)
