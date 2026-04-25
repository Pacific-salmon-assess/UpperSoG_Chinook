# Extract ERs

ERM_QuinsamCampbell <- readRDS("CM/QuinsamCampbell_04.21.26.rds")
report <- salmonMSE:::get_report(ERM_QuinsamCampbell)

ERM_Adam <- readRDS("CM/Adam_04.24.26.rds")
report_Adam <- salmonMSE:::get_report(ERM_Adam)

ERM_AdamPhillips <- readRDS("CM/AdamPhillips_CM_04.25.26.rds")
report_AdamPhillips <- salmonMSE:::get_report(ERM_AdamPhillips)

ERM_Salmon <- readRDS("CM/Salmon_04.24.26.rds")
report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

ERM_Woss <- readRDS("CM/Woss_04.24.26.rds")
report_Woss <- salmonMSE:::get_report(ERM_Woss)

ERM_Phillips <- readRDS("CM/Phillips_CM_04.25.26.rds")
report_Phillips <- salmonMSE:::get_report(ERM_Phillips)

# ER by age class - array by MCMC simulation x year x age
FPT <- sapply(report_Phillips, function(i) {
  outer(i$FPT, i$vulPT)
}, simplify = "array") %>%
  aperm(c(3, 1, 2))
UPT <- 1 - exp(-FPT)
dim(UPT)
# Summarize by age, over draws and years
apply(UPT, 3, quantile, probs = c(0.025, 0.5, 0.975), na.rm=T)


# Alternative using code extracted from CM_ER() in CM_fun.R in salmonMSE
# accounts for AEQ
# https://github.com/Blue-Matter/salmonMSE/blob/85e10c26277f3799755edadbeae3f147f67063b5/R/CMfun.R#L1038

# However, I think this code makes the assumption that:
# total returns = PT_catch + T_catch + escapement , ie. the denominator
# For Upper SoG terminal ER is set to zero in the code, but s_enroute accounts
# for both terminal ER and return migration.
# I'm not sure if s_enroute is considered in 'escapement' time-series in the
# denominator, i.e., if it is escapement before or after enroute mortality?


brood <- FALSE
type <- "PT"
year1 <-  1
ci <-  TRUE
r <-  1

# Make quantile figure time series
salmonMSE:::CM_ER(report_AdamPhillips, brood, type, year1, ci, r, at_age = FALSE)

# Get individual values by MCMC simulation x year
ER_list <- lapply(report_AdamPhillips, function(i) {

  esc <- apply(i$escyear, 1:2, sum)
  morts_PT <- apply(i$cyearPT, 1:2, sum)
  morts_T <- apply(i$cyearT, 1:2, sum)

  AEQ_PT <- salmonMSE:::calc_AEQ(i)[, , r] # Always by release year
  AEQ_T <- array(1, dim(esc))

  if (brood) {
    esc <- CY2BY(esc)
    morts_PT <- CY2BY(morts_PT)
    morts_T <- CY2BY(morts_T)

    denom <- rowSums(morts_PT * AEQ_PT + morts_T * AEQ_T + esc)

  } else {
    AEQ_PT2 <- array(NA_real_, dim(esc)) # Re-index to align with calendar year
    nt <- nrow(AEQ_PT2)
    na <- ncol(AEQ_PT2)
    for (t in 1:nt) {
      for (a in 1:na) if (t-a+1>0) AEQ_PT2[t, a] <- AEQ_PT[t-a+1, a]
    }
    denom <- rowSums(morts_PT * AEQ_PT2 + morts_T * AEQ_T + esc)
  }

  if (type == "PT") {
    num <- rowSums(morts_PT * AEQ_PT)
  } else if (type == "T") {
    num <- rowSums(morts_T * AEQ_T)
  } else {
    num <- rowSums(morts_PT * AEQ_PT + morts_T * AEQ_T)
  }

  list(ER = num/denom)
})

# Pull elements out of ER_list and put in a data frame
df <- as.data.frame(do.call(
  rbind,
  lapply(ER_list, function(x) x[["ER"]])#[35:40])# could look at just the last 6 years
))


# Plot all simulations
as.matrix(df) |> t() |> matplot(typ = 'l')

# Take the mean/quantiles over all years and posterior draws
quantile(as.matrix(df), probs = c(0.025, 0.5, 0.975), na.rm=T)
quantile(as.matrix(df), probs = 0.5, na.rm=T)
#0.47
