# Extract ERs

ERM_QuinsamCampbell <- readRDS("CM/QuinsamCampbell_05.01.26.rds")
report_QC <- salmonMSE:::get_report(ERM_QuinsamCampbell)

ERM_Adam <- readRDS("CM/Adam_04.24.26.rds")
report_Adam <- salmonMSE:::get_report(ERM_Adam)

ERM_AdamPhillips <- readRDS("CM/AdamPhillips_CM_05.06.26.rds")
report_AdamPhillips <- salmonMSE:::get_report(ERM_AdamPhillips)

ERM_SalmonPhillips <- readRDS("CM/SalmonPhillips_CM_05.06.26.rds")
report_SalmonPhillips <- salmonMSE:::get_report(ERM_SalmonPhillips)

ERM_WossPhillips <- readRDS("CM/WossPhillips_CM_05.06.26.rds")
report_WossPhillips <- salmonMSE:::get_report(ERM_WossPhillips)

ERM_Salmon <- readRDS("CM/Salmon_04.24.26.rds")
report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

ERM_Woss <- readRDS("CM/Woss_04.24.26.rds")
report_Woss <- salmonMSE:::get_report(ERM_Woss)

ERM_Phillips <- readRDS("CM/Phillips_CM_05.06.26.rds")
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
year1 <-  2010 # get this for each pop: min(full_matrix$RELEASE_YEAR) or min(full_table$RELEASE_YEAR)
ci <-  TRUE
r <-  1

# Make quantile figure time series
salmonMSE:::CM_ER(report_SalmonPhillips, brood, type, year1, ci, r, at_age = FALSE)

x_QC <- CM_ERv2(report_Adam, brood, type, year1, ci, r, at_age = FALSE)
x_Phillips <- CM_ERv2(report_AdamPhillips, brood, type, year1, ci, r, at_age = FALSE)

ts_QC <- sapply(x_QC, getElement, "ER") %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 + year1 - 1) %>%
  reshape2::dcast(list("Year", "Var1"), value.var = "value")
ts_QC$CWT <- "Q/C"

ts_Phillips <- sapply(x_QC, getElement, "ER") %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 + year1 - 1) %>%
  reshape2::dcast(list("Year", "Var1"), value.var = "value")
ts_Phillips$CWT <- "Phillips"

ts <- bind_rows(ts_QC, ts_Phillips)

# THIS DOESNOT WORK YET
g <- ggplot(ts, aes(.data$Year, .data$`50%`, fill=factor(CWT))) +
  geom_line() +
  labs(x = xlab, y = "AEQ pre-terminal ER") +
  expand_limits(y = 0) +
  geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), fill = alpha("grey", 0.5))
g

# Get individual values by MCMC simulation x year
ER_list <- lapply(report_SalmonPhillips, function(i) {

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




#_________________________________________________\

CM_ERv2 <- function(report, brood = TRUE, type = c("PT", "T", "all"), year1 = 1, ci = TRUE, at_age = TRUE, r = 1) {

  type <- match.arg(type)

  if (at_age) {

    ER_list <- lapply(report, function(i) {
      name <- ifelse(type == "PT", "survPT", "survT")
      CYER <- 1 - i[[name]]
      if (brood) {
        ER <- CY2BY(CYER)
      } else {
        ER <- CYER
      }
      return(list(ER = ER))
    })

    g <- salmonMSE:::.CM_statevarage(
      ER_list,
      year1 = if (brood) year1 - 1 else year1,
      ci,
      "ER",
      ylab = switch(
        type,
        "PT" = "Preterminal exploitation rate",
        "T" = "Terminal exploitation rate"
      ),
      xlab = ifelse(brood, "Brood Year", "Calendar Year"),
      scales = "fixed"
    )

  } else {

    ER_list <- lapply(report, function(i) {

      esc <- apply(i$escyear, 1:2, sum)
      morts_PT <- apply(i$cyearPT, 1:2, sum)
      morts_T <- apply(i$cyearT, 1:2, sum)

      AEQ_PT <- salmonMSE:::calc_AEQ(i)[, , r] # Always by release year
      AEQ_T <- array(1, dim(esc))

      if (brood) {
        esc <- salmonMSE:::CY2BY(esc)
        morts_PT <- salmonMSE:::CY2BY(morts_PT)
        morts_T <- salmonMSE:::CY2BY(morts_T)

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
    return(ER_list)
    # g <- .CM_ts(
    #   ER_list,
    #   year1 = if (brood) year1 - 1 else year1,
    #   ci,
    #   var = "ER",
    #   xlab = ifelse(brood, "Brood Year", "Calendar Year"),
    #   ylab = switch(
    #     type,
    #     "PT" = "Preterminal exploitation rate",
    #     "T" = "Terminal exploitation rate",
    #     "all" = "Total exploitation rate"
    #   )
    # )
  }

  #g
}
