# calc_eq_PT calculates preterminal UMSY with fixed terminal exploitation rate
# updated from salmonMSE:::.calc_eq()

calc_eq_PT <- function(.F, Mjuv, fec, p_female, gamma = 1, U_T = 0.1, vulPT, vulT, p_mature, s_enroute = 1, n_g = 1,
                       p_LHG = rep(1/n_g, n_g),
                       SRRpars, maximize = c("MSY", "MER"), AEQ = TRUE, opt = TRUE, aggregate_age = TRUE) {

  maximize <- match.arg(maximize)

  if (!is.matrix(Mjuv)) Mjuv <- matrix(Mjuv, length(Mjuv), n_g)
  if (!is.matrix(p_mature)) p_mature <- matrix(p_mature, length(p_mature), n_g)
  if (!is.matrix(fec)) fec <- matrix(fec, length(fec), n_g)

  FPT <- vulPT * .F[1]
  surv_juv <- sapply(1:n_g, function(g) p_LHG[g] * salmonMSE:::calc_survival(Mjuv[, g] + FPT[-length(FPT)], p_mature[, g]))
  surv_return <- surv_juv * p_mature

  if (U_T > 0) {
    opt_fn <- function(.F, U_T = 0.1) {
      catch <- surv_return * (1 - exp(-vulT * .F))
      sum(catch)/sum(surv_return) - U_T
    }
    opt_UT <- uniroot(opt_fn, c(1e-8, 10))
    Fopt <- opt_UT$root
  } else {
    Fopt <- 0
  }

  FT <- vulT * Fopt
  surv_esc <- surv_return * exp(-FT)
  surv_spawn <- surv_esc * s_enroute

  EPR <- surv_spawn * gamma * p_female * fec # Egg per smolt at F

  Smolt <- calc_smolt(
    sum(EPR),
    kappa = SRRpars["kappa"], capacity = SRRpars["capacity"], Smax = SRRpars["Smax"], phi = SRRpars["phi"],
    tau = SRRpars["tau"], fitness_loss = 1, SRrel = as.character(SRRpars["SRrel"]), per_recruit = TRUE
  ) |>
    as.numeric()

  Njuv <- Smolt * surv_juv
  Return <- Smolt * surv_return
  KPT <- Njuv * (1 - exp(-FPT))
  KT <- Return * (1 - exp(-FT))

  AEQ_PT <- array(NA_real_, dim(KPT))
  AEQ_PT[nrow(AEQ_PT), ] <- 1
  for (a in seq(nrow(AEQ_PT) - 1, 1)) {
    AEQ_PT[a, ] <- AEQ_PT[a+1, ] * (1 - p_mature[a, ]) * exp(-Mjuv[a, ] - FPT[a]) + p_mature[a, ]
  }

  KPT_AEQ <- KPT * AEQ_PT

  Escapement <- Smolt * surv_esc
  Spawners <- Smolt * surv_spawn

  if (opt) {

    if (maximize == "MSY") {

      if (AEQ) {
        return(sum(KPT_AEQ, KT))
      } else {
        return(sum(KPT, KT))
      }

    } else if (maximize == "MER") {

      if (AEQ) {
        return(sum(Return, KPT_AEQ) - sum(Spawners))
      } else {
        return(sum(Return, KPT) - sum(Spawners))
      }

    }

  } else {

    output <- list(
      KPT = KPT,
      KPT_AEQ = KPT_AEQ,
      KT = KT,
      FPT = max(FPT),
      FT = max(FT),
      UPT_EQ = sum(KPT_AEQ)/sum(KPT_AEQ, Return),
      UT = sum(KT)/sum(Return),
      Smolt = Smolt,
      Njuv = Njuv,
      Return = Return,
      Escapement = Escapement,
      Spawners = Spawners,
      Egg = Smolt * EPR
    )

    if (aggregate_age) output <- lapply(output, sum)

    return(output)
  }
}

compare_ref <- function(a = 3, # Units of adults/spawner
                        Smax = 1000, # Units of spawners
                        U_T = 0.1, # Fixed terminal ER
                        maxage = 5,
                        p_female = 1,
                        vulPT = rep(0, maxage),
                        vulT = c(0, 0, 0.2, 0.5, 1),
                        M = c(1, 0.8, 0.6, 0.4),
                        fec = c(0, 1500, 3000, 3200, 3500),
                        p_mature = c(0, 0.1, 0.2, 0.3, 1),
                        maximize = c("MSY", "MER")) {

  require(salmonMSE)

  maximize <- match.arg(maximize)

  # Calculate eggs per smolt
  phi <- salmonMSE:::calc_phi(
    Mjuv = M,
    p_mature = p_mature,
    p_female = p_female,
    fec = fec,
    s_enroute = 1
  )

  # To convert Smax from units of spawners to eggs:
  # 1. Calculate Srep from Ricker SRR, set R = S --> S0 = log(a)/b
  # 2. Calculate spawners per juvenile (tau)
  # 3. Calculate smolt per egg (phi)
  # 4. Egg production at replacement (Erep) = unfished spawners * juvenile per spawner * egg per smolt
  # 5. Emax = log(a) / unfished eggs
  tau <- salmonMSE:::calc_phi(
    Mjuv = M,
    p_mature = p_mature,
    p_female = p_female,
    fec = fec,
    s_enroute = 1,
    output = "spawner"
  )

  Srep <- Smax * log(a)
  Erep <- Srep / tau * phi
  beta <- log(a)/Erep
  Emax <- 1/beta

  SRRpars <- data.frame(
    kappa = a,
    Smax = Smax,
    Emax = Emax,
    tau = tau,
    phi = phi,
    SRrel = "Ricker"
  )

  # salmonMSE ref with fixed terminal ER
  MSY_opt <- optimize(
    calc_eq_PT, interval = c(1e-8, 10),
    Mjuv = M,
    fec = fec,
    p_female = p_female,
    U_T = U_T,
    vulPT = vulPT,
    vulT = vulT,
    p_mature = p_mature,
    s_enroute = 1,
    SRRpars = SRRpars,
    maximize = maximize,
    AEQ = TRUE,
    maximum = TRUE
  )

  ref <- calc_eq_PT(
    .F = MSY_opt$maximum,
    Mjuv = M,
    fec = fec,
    p_female = p_female,
    vulPT = vulPT,
    vulT = vulT,
    p_mature = p_mature,
    s_enroute = 1,
    SRRpars = SRRpars,
    maximize = maximize,
    AEQ = TRUE,
    opt = FALSE
  ) |>
    unlist()

  names(ref) <- paste0(names(ref), "_MSY")
  if (ref["KPT_AEQ_MSY"] == 0) ref["UPT_AEQ_MSY"] <- 0
  if (ref["KT_MSY"] == 0) ref["UT_MSY"] <- 0

  #ref["Sgen"] <- calc_Sgen(
  #  M, fec, p_female, rel_F = c(0, 1), vulPT, vulT, p_mature,
  #  s_enroute = 1, SRRpars = SRRpars, SMSY = ref["Spawners_MSY"]
  #)
  ref["Catch"] <- sum(ref["KPT_AEQ_MSY"], ref["KT_MSY"])
  ref["Catch/Return"] <- sum(ref["KPT_AEQ_MSY"], ref["KT_MSY"])/(sum(ref["KPT_AEQ_MSY"], ref["Return_MSY"]))

  ref_salmonMSE <- structure(
    ref[c("UPT_EQ_MSY", "UT_MSY", "Catch", "Catch/Return", "Spawners_MSY", "Sgen")],
    names = c("UMSY (preterminal)", "UMSY (terminal)", "Catch", "Catch/Return", "SMSY", "Sgen")
  ) |>
    round(4)

  # Ricker SRR calcs
  ref_Ricker <- local({
    umsy <- calc_Umsy_Ricker(log(a))
    smsy <- calc_Smsy_Ricker(log(a), 1/Smax)
    sgen <- calc_Sgen_Ricker(log(a), 1/Smax)
    structure(c(umsy, smsy, sgen), names = c("UMSY", "SMSY", "Sgen"))
  }) |>
    round(4)

  output <- list(
    salmonMSE = ref_salmonMSE,
    Ricker = ref_Ricker
  )
  return(output)
}


compare_ref(
  a = 3,
  Smax = 1000,
  maxage = 5,
  p_female = 1,
  vulPT = c(0, 0.1, 0.2, 0.4, 1),
  vulT = c(0, 0.1, 0.2, 0.4, 1),
  M = c(1, 0.3, 0.2, 0.1),
  fec = rep(1, 5),
  p_mature = c(0, 0.1, 0.2, 0.3, 1)
)


#############################################################################
# Apply to Salmon River
# Input Conditioning Model results
maxage <- 5
nsim <- 500
ERM <- readRDS("CM/Salmon_10.22.26.base.rds")
report <- salmonMSE:::get_report(ERM)
d <- salmonMSE:::get_CMdata(ERM@.MISC$CMfit)

# Take maturity average from the 5 most recent brood years (2020-2024)
matt <- sapply(report, getElement, "matt", simplify = "array") %>%
  aperm(c(1, 2, 4, 3))
matt_avg <- apply(matt[seq((d$Ldyr - 4), d$Ldyr), , ,1], c(2,3), mean)
p_mature <- apply (matt_avg,1,median)


### Natural mortality - NOS ----
# Marine survival from CTC 23-06 p. 9
M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1))
# Age 1 value of Marine survival
df <- sapply(report, getElement, "mo", simplify = "array")[,1,] # year x age x sim, only age 1
Mjuv <- median(df)

# Age 2-5
Mjuv[2] <- M_CTC[2]
Mjuv[3] <- M_CTC[3]
Mjuv[4] <- M_CTC[4]


## Eggs/total spawner: Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.
aveFecFemale <- 5871 # From Brendan Zoehner, SEP, assuming for age 4
fec_QC <- c(0, 0, 800, 2000, 2500) # fecundity by age for Quinsam
ppn4 <- fec_QC/fec_QC[4] # proportion of fecundity rleative to age 4
fec_Salmon <- round(aveFecFemale * 0.5 * ppn4) # fecundity per total spawner by age for Salmon River

# This is assumed to already account for ppn female as age-3 fecundity is so low
# Accounts for pre-spawn mortality (5%)
fec <- fec_Salmon * 0.95


phi <- sapply(report, getElement, "epro")
tau <- sapply(report, getElement, "spro")

SRbeta <- median(sapply(report, getElement, "beta"))
Emax <- 1/SRbeta
Smax <- median(Emax * tau/phi)
kappa <- median(as.numeric(exp(pars$log_cr))) #a


vulPT <- sapply(report, getElement, "vulPT") %>% apply(1, median)
vulT <- sapply(report, getElement, "vulT") %>% apply(1, median)


# Terminal ER, median overall years and MC trials
year <- year1 + seq(1, d$Ldyr) - 1
year_borrow <- seq(max(year) - 9, max(year) - 5)


UT <- salmonMSE:::.CM_ER(report, brood = FALSE, type = "T",  r = 1,
                         index_AEQ = match(year_borrow, year)) %>% median(na.rm=T)


compare_ref(
  a = kappa,
  Smax = Smax,
  maxage = 5,
  p_female = 1,
  vulPT = vulPT,
  vulT = vulT,
  M = Mjuv,
  fec = fec,
  p_mature = p_mature,
  U_T = UT
)
