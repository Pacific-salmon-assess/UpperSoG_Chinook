# Extract ERs

ERM_QuinsamCampbell <- readRDS("CM/QuinsamCampbell_06.15.26.rds")
report_QC <- salmonMSE:::get_report(ERM_QuinsamCampbell)

ERM_Adam <- readRDS("CM/Adam_06.15.26.rds")#readRDS("CM/Adam_06.01.26.JSpt.rds")
report_Adam <- salmonMSE:::get_report(ERM_Adam)

ERM_AdamPhillips <- readRDS("CM/AdamPhillips_CM_05.09.26.rds")
report_AdamPhillips <- salmonMSE:::get_report(ERM_AdamPhillips)

ERM_SalmonPhillips <- readRDS("CM/SalmonPhillips_CM_05.09.26.rds")
report_SalmonPhillips <- salmonMSE:::get_report(ERM_SalmonPhillips)

ERM_WossPhillips <- readRDS("CM/WossPhillips_CM_05.09.26.rds")
report_WossPhillips <- salmonMSE:::get_report(ERM_WossPhillips)

ERM_Salmon <- readRDS("CM/Salmon_06.15.26.x2esc.rds") #file:///C:/github/UpperSoG_Chinook/CM/Salmon_06.10.JSesc.NGTSt.html
report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

ERM_Woss <- readRDS("CM/Woss_06.15.26.rds")
report_Woss <- salmonMSE:::get_report(ERM_Woss)

ERM_Phillips <- readRDS("CM/Phillips_CM_05.08.26.rds")
report_Phillips <- salmonMSE:::get_report(ERM_Phillips)

# ER by age class - array by MCMC simulation x year x age
FPT <- sapply(report_Salmon, function(i) { # Change to Adam, Salmon, Woss
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


samp <- ERM_Salmon#readRDS(paste("CM/",pop,"_06.01.26.JSpt.rds", sep=""))#paste("ERM_", pop, "Phillips", sep="")
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
brood <- FALSE
type <- "PT"# "T" #
year1 <-  2002#2002#2001#1984#2010 # get this for each pop: min(full_matrix$RELEASE_YEAR) or min(full_table$RELEASE_YEAR)
ci <-  TRUE
r <-  1
report <- report_Salmon



# Repeat to fill in AEQs for incomplete brood years
year <- year1 + seq(1, d$Ldyr) - 1
year_borrow <- seq(max(year) - 9, max(year) - 5)
g <- CM_ER(report, brood, type, year1, ci, r, at_age = FALSE, index_AEQ = match(year_borrow, year))
g

# Get the quantiles of what we just plotted
if (packageVersion("ggplot2") >= "4.0") {
  g@data
} else {
  g$data
}
# Make quantile figure time series
CM_ER(report, brood, type="PT", year1, ci, r, at_age = FALSE, index_AEQ = match(year_borrow, year))
CM_ER(report, brood, type="T", year1, ci, r, at_age = FALSE, index_AEQ = match(year_borrow, year))

# Get matrices of individual AEQ ER values by year x MCMC simulation
# See ?.CM_ER
df.PT <- salmonMSE:::.CM_ER(report, type = "PT", r, brood, index_AEQ = match(year_borrow, year))
df.T <- salmonMSE:::.CM_ER(report, type = "T", r, brood, index_AEQ = match(year_borrow, year))
df.all <- salmonMSE:::.CM_ER(report, type = "all", r, brood, index_AEQ = match(year_borrow, year))

# # Plot all simulations
# as.matrix(df) |> t() |> matplot(typ = 'l')

# Take the mean/quantiles over all years and posterior draws
quantile(df.PT, probs = c(0.025, 0.5, 0.975), na.rm=T)
quantile(df.T, probs = c(0.025, 0.5, 0.975), na.rm=T)
quantile(df.all, probs = c(0.025, 0.5, 0.975), na.rm=T)

# Get time-series of ERs
ER.pt <- apply(df.PT, 1, function(row) quantile(row, probs = c(0.025, 0.5, 0.975), na.rm=T))
ER.t <- apply(df.T, 1, function(row) quantile(row, probs = c(0.025, 0.5, 0.975), na.rm=T))
ER.all <- apply(df.all, 1, function(row) quantile(row, probs = c(0.025, 0.5, 0.975), na.rm=T))

# ER <- as.data.frame(t(1-(1-ER.pt)*(1-ER.t)))
ER <- as.data.frame(t(ER.all))
ER <- ER %>% rename("lower"="2.5%", "median"="50%" ,"upper"="97.5%") %>%
  mutate (label="ER", Year = seq(from = year1, length.out = length(ER[,1])))



# Compare to UMSY (see file 98-ExtractMSY.R)

#----------------------------------------------------------------------
# Plot comparing AEQ ERs for Adams River CM based on either Q/C or Phillips ERs
# Get ER_list for each version of CM that includes AEQ ERs (see function below)

 # x_QC <- CM_ERv2(report_Adam, brood, type, 2002, ci, r, at_age = FALSE)
 # x_Phillips <- CM_ERv2(report_AdamPhillips, brood, type, 2010, ci, r, at_age = FALSE)
# x_QC <- CM_ERv2(report_Salmon, brood, type, 2002, ci, r, at_age = FALSE)
# x_Phillips <- CM_ERv2(report_SalmonPhillips, brood, type, 2010, ci, r, at_age = FALSE)
x_QC <- CM_ERv2(report_Woss, brood, type, 2002, ci, r, at_age = FALSE)
x_Phillips <- CM_ERv2(report_WossPhillips, brood, type, 2010, ci, r, at_age = FALSE)

ts_QC <- sapply(x_QC, getElement, "ER") %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 + 2002 - 1) %>% #Time-series for Q/C starts in 2002
  reshape2::dcast(list("Year", "Var1"), value.var = "value")
ts_QC$CWT <- "Quinsam/Campbell"

ts_Phillips <- sapply(x_Phillips, getElement, "ER") %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 + 2010 - 1) %>%#Time-series for Q/C starts in 2010
  reshape2::dcast(list("Year", "Var1"), value.var = "value")
ts_Phillips$CWT <- "Phillips"

ts <- bind_rows(ts_QC, ts_Phillips)


g <- ggplot(ts, aes(.data$Year, .data$`50%`, colour = factor(.data$CWT), fill=factor(.data$CWT))) +
  geom_line(aes(colour=factor(.data$CWT))) +
  labs(x = xlab, y = "AEQ pre-terminal ER") +
  expand_limits(y = 0) +
  geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), colour=NA, alpha=0.2) +
  theme(legend.position = "bottom", legend.title = element_blank())
g
r <- cor(ts_QC$'50%'[9:22],ts_Phillips$'50%'[1:14], use="pairwise.complete.obs")
r^2


#_________________________________________________
# Function to get ER_list as input to time-series plots of AEQ ERs

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


#----------------------------------------------------------------
# Extract maturity

report <- report_SalmonPhillips#report_WossPhillips#report_QC
year1 <- 2010#1981
samp <- ERM_SalmonPhillips
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
r <- d$r_matt
brood <- TRUE
rs_names <- "smolt0+"


  n_r <- d$n_r
  if (missing(rs_names)) rs_names <- seq(1, n_r)


    bmatt <- data.frame(Age = 1:d$Nages, value = d$bmatt) %>%
      dplyr::filter(Age > 1)

      matt <- sapply(report, function(i) salmonMSE:::CY2BY(i[["matt"]][, , r]), simplify = 'array')
    matt_q <- apply(matt, 1:2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
      reshape2::melt() %>%
      mutate(Year = Var2 + year1 - 1) %>%
      rename(Age = Var3) %>%
      reshape2::dcast(Age + Var2 + Year ~ Var1) %>%
      filter(!is.na(`50%`))

      matt_q$Year <- matt_q$Year - 1 # Release year offset by 1

    g <- matt_q %>%
      dplyr::filter(Age > 1) %>%
      ggplot(aes(Year, .data$`50%`, fill = factor(.data$Age), colour = factor(.data$Age))) +
      geom_line() +
      geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
      geom_hline(data = bmatt, linetype = 2, aes(yintercept = .data$value, colour = factor(.data$Age))) +
      labs(x = ifelse(brood, "Brood Year", "Return Year"), y = "Proportion mature", colour = "Age", fill = "Age")
  g


  # For Q/C take >2010 values (approximately last 10 years)
  matt_q <- rename(matt_q, median='50%')#|> filter(Year>2010)
  out <- matt_q |> group_by(Age) |> summarise(mean_ppn = mean(median, na.rm=T))
  matt_q
  out
  age3 <- out$mean_ppn[3]-out$mean_ppn[2]
  age4 <- out$mean_ppn[4]-out$mean_ppn[3]
  age5 <- out$mean_ppn[5]-out$mean_ppn[4]
  age3
  age4
  age5

#----------------------------------------------------------------------\
# plot prod
   samp <- ERM_Salmon

  d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
  prod <- salmonMSE:::.CM_prod(report_Salmon, d)

  prod_q <- apply(prod, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1"))

  g <- prod_q %>%
    ggplot(aes(Year, .data$`50%`)) +
    geom_line() +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
    labs(x = "Calendar Year", y = "Productivity") +
    geom_hline(yintercept=1)
  g
