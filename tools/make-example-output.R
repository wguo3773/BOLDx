pkgload::load_all(".",quiet=TRUE)
methods <- c("prediction","tite","combined","bold","boin","tite_boin")
cfg <- bold_config(nmax=21,target=.25,tau=.50,sensitivity=.90,specificity=.98,
                   eligibility=1,rho=.80)
single <- simulate_bold(c(.10,.25,.45,.60),methods,cfg,nsim=100,seed=20261031,mtd=2)
write.csv(single$summary,"docs/example-summary.csv",row.names=FALSE)

sc <- read.csv("inst/extdata/example-panel.csv")
stopifnot(nrow(sc)==50,all(sc$mtd==2),all(sc$delta==.15))
panel <- simulate_scenarios(as.matrix(sc[paste0("p",1:4)]),mtd=sc$mtd,
                           method=methods,config=cfg,nsim=100,seed=20261031)
panel$summary$delta <- sc$delta[panel$summary$scenario]
aggregate <- summarise_scenarios(panel$summary,c("method","mtd","delta"))
write.csv(panel$summary,"docs/example-panel-estimates.csv",row.names=FALSE)
write.csv(aggregate,"docs/example-panel-mean-sd.csv",row.names=FALSE)
print(aggregate)
