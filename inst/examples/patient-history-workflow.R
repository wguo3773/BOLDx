library(BOLDx)
cfg <- bold_config(target=.25,nmax=21,tau=.50)

# Between-cohort evaluation: three patients with complete follow-up at dose 1.
h <- data.frame(dose=c(1,1,1),dlt=c(0,0,0),followup=c(3,3,3))
original <- recommend_dose(h,ndoses=4,current_dose=1,method="bold",config=cfg)
print(original)
original$dose_summary

# External classifier labels; record the interventions actually performed.
plan_intervention(c(1L,0L,0L),eligible=TRUE)
hp <- transform(h,prediction=c(1L,0L,0L),eligible=TRUE,
                intervened=c(TRUE,FALSE,FALSE))
print(recommend_dose(hp,4,1,"prediction",cfg))

# TITE uses observed follow-up durations and observed DLTs to date.
ht <- h
ht$followup <- c(3,3,1)
print(recommend_dose(ht,4,1,"tite",cfg))

# Combined mode first waits for an intervened patient to resolve.
hc <- transform(h,prediction=c(0L,1L,0L),eligible=TRUE,
                intervened=c(FALSE,TRUE,FALSE))
hc$followup <- c(3,1,1)
print(recommend_dose(hc,4,1,"combined",cfg))
hc$followup[2] <- 3
combined <- recommend_dose(hc,4,1,"combined",cfg)
print(combined)
combined$patient_status

# With each later evaluation, append/update the observed history, specify the
# latest enrolled cohort's dose, and carry forward eliminated doses:
# recommend_dose(updated_history,4,current_dose=latest_dose,
#                method="combined",config=cfg,excluded=combined$excluded_doses)
