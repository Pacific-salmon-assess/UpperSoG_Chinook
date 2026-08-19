
# Make OM

# Libraries
library(tidyverse)
library(salmonMSE)

# OM specificiations
maxage <- 5
nsim <- 500
proyears <- 30
n_g <- 1
year1 <- 2002 # from conditioning model
msurv <- "low"#"base" #"low" #"high"
# for estimated value or 0.001 or 0.1 for lower and upper sens. anal.
# also used for file name of SOM
som_name <- "Salmon low msurv"#"Salmon base" #"low msurv" # "high msurv" # for internal naming

# Load exploitation rate model - Salmon River
# ERM <- readRDS("CM/Salmon_08.06.26.prior.rds")
ERM <- readRDS("CM/Salmon_08.18.26.prior.rds")

report <- salmonMSE:::get_report(ERM)
set.seed(24)
sim_samp <- sample(seq(1, length(report)), nsim)
d <- salmonMSE:::get_CMdata(ERM@.MISC$CMfit)
pars <- rstan::extract(ERM)

# Take maturity average from the 6 most recent brood years (2020-2024)
# Consider changing to earlier/longer period if recent ppns are less reliable
matt <- sapply(report, getElement, "matt", simplify = "array") %>%
  aperm(c(1, 2, 4, 3))
matt_avg <- apply(matt[seq((d$Ldyr - 4), d$Ldyr), , ,1], c(2,3), mean)
# apply (matt_avg,1,median)


### Natural mortality - NOS ----
Mjuv_NOS <- array(0, c(nsim, maxage-1, proyears, n_g))

# Marine survival from CTC 23-06 p. 9
M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1))

# Age 1 value of Marine survival
# Draw random values
df <- sapply(report, getElement, "mo", simplify = "array")[,1,] # year x age x sim, only age 1
#exp(-range(apply(df, 1, median, na.rm=T)))
set.seed(23)
sim_samp_long <- sample(seq(1, length(report))*d$Ldyr, nsim)
m1_sample <- df[sim_samp_long]
exp(-quantile(df[sim_samp_long], prob = c(0.25, 0.5, 0.75)))

# Alternative approach: making a parametric distribution from all years
if(msurv %in% c("low", "high")){
  df <- sapply(report, getElement, "mo", simplify = "array") %>%
    apply(1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    reshape2::melt() %>%
    mutate(Year = Var2) %>%
    mutate(Age = paste("Age", Var3)) %>%
    reshape2::dcast(Year + Age ~ Var1, value.var = "value") %>%
    filter(Age=='Age 1')

  # Transform percentage ER to logistic distribution before sampling

  if(msurv == "low") s1_mu <- mean(qlogis(0.002))
  if(msurv == "high") s1_mu <- mean(qlogis(0.01))
  s1_sd <- sd(qlogis(exp(-df$'50%')))

  set.seed(234)
  s1_sim  <- matrix(rnorm(nsim * proyears, s1_mu, s1_sd), nsim, proyears)
  # Back transform ERs to % and then to instantaneous mortality
  s1_sim_trans <- plogis(s1_sim)
  m1_sample <- -log(s1_sim_trans)
}

Mjuv_NOS[, 1, , 1] <- m1_sample


# Age 2-5
Mjuv_NOS[, 2, , ] <- M_CTC[2]
Mjuv_NOS[, 3, , ] <- M_CTC[3]
Mjuv_NOS[, 4, , ] <- M_CTC[4]

### Maturity ----
# Assume average maturity
p_mature <- array(0, c(nsim, maxage, proyears))
p_mature[] <- matt_avg[, sim_samp] %>% t() %>%
  array(c(nsim, maxage,  proyears))

# matt.x <- apply(p_mature,2,mean)
# matt_ppn <- NA
# matt_ppn[1] <- matt.x[1]
# matt_ppn[2] <- matt.x[2]-matt.x[1]
# matt_ppn[3] <- matt.x[3]-matt.x[2]
# matt_ppn[4] <- matt.x[4]-matt.x[3]
# matt_ppn[5] <- matt.x[5]-matt.x[4]

# Eggs/total spawner: Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.
aveFecFemale <- 5871 # From Brendan Zoehner, SEP, assuming for age 4
fec_QC <- c(0, 0, 800, 2000, 2500) # fecundity by age for Quinsam
ppn4 <- fec_QC/fec_QC[4] # proportion of fecundity rleative to age 4
fec_Salmon <- round(aveFecFemale * 0.5 * ppn4) # fecundity per total spawner by age for Salmon River

# This is assumed to already account for ppn female as age-3 fecundity is so low
# Accounts for pre-spawn mortality (5%)
fec <- fec_Salmon * 0.95



phi <- sapply(report[sim_samp], getElement, "epro")
tau <- sapply(report[sim_samp], getElement, "spro")

SRbeta <- sapply(report[sim_samp], getElement, "beta")
Emax <- 1/SRbeta
Smax <- Emax * tau/phi
kappa <- as.numeric(exp(pars$log_cr[sim_samp]))

Bio <- new(
  "Bio",
  maxage = maxage,
  n_g = n_g,
  p_LHG = 1,
  p_mature = p_mature,# Includes random variability
  # Freshwater survival and pre-spawn M included in habitat object
  SRrel = "Ricker",
  kappa = kappa,# Includes random variability
  Smax = Smax,# Includes random variability
  phi = phi,
  tau = tau,
  Mjuv_NOS = Mjuv_NOS, # Includes random variability in 1st year M
  p_female = 1, # p_female specified in fecundity
  fec = fec,
  s_enroute = (1 - 0.1)# M.  Clarke pers. comm. (2025), same as CM
)

### Harvest, fishery vulnerability ----
# To estimate AEQ ERs, borrow information on population parameters
# for the last generation
year <- year1 + seq(1, d$Ldyr) - 1
year_borrow <- seq(max(year) - 9, max(year) - 5)

# AEQ ERs (matrix of dimensions Mcdraws x years)
UPT <- salmonMSE:::.CM_ER(report, brood = FALSE, type = "PT",  r = 1,
                          index_AEQ = match(year_borrow, year)) %>% t()
UT <- salmonMSE:::.CM_ER(report, brood = FALSE, type = "T",  r = 1,
                         index_AEQ = match(year_borrow, year)) %>% t()

# Take median AEQ UPT and UT in the last 5 years, and take the average of those values to project forward
UPTrange <-  UPT %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975), na.rm=T) %>%
  t()
UPTmed <- UPTrange[seq((d$Ldyr - 4),d$Ldyr), "50%"] %>%
  mean()
UTrange <-  UT %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975), na.rm=T) %>%
  t()
UTmed <- UTrange[seq((d$Ldyr - 4),d$Ldyr), "50%"] %>%
  mean()

u_preterminal <- matrix(UPTmed, nrow=nsim, ncol=proyears) #No random variability included
u_terminal <- matrix(UTmed, nrow=nsim, ncol=proyears) #No random variability included

vulPT <- sapply(report[sim_samp], getElement, "vulPT")
vulT <- sapply(report[sim_samp], getElement, "vulT")

# Mean over MC samples
# sapply(report, getElement, "vulPT") %>% apply(1, mean)
# sapply(report, getElement, "vulT") %>% apply(1, mean)


Harvest <- new(
  "Harvest",
  type_PT = "u",
  type_T = "u",
  u_preterminal = u_preterminal,
  u_terminal = u_terminal,
  MSF_PT = FALSE,
  MSF_T = FALSE,
  release_mort = c(0.1,0.1), # mark-selective fishery not included for this stock
  vulPT = t(vulPT),
  vulT = t(vulT)
)


### Hatchery ----
# Number of release groups
n_r <- 1

# Natural mortality HOS
Mjuv_HOS <- array(0, c(nsim, maxage-1, proyears, n_r))

# Age 1 survival
# Assume same as NOS- including internannual variability
Mjuv_HOS[, 1, , 1] <- Mjuv_NOS[, 1, , 1]

# Age 2 survival
Mjuv_HOS[, 2, , ] <- M_CTC[2]
Mjuv_HOS[, 3, , ] <- M_CTC[3]
Mjuv_HOS[, 4, , ] <- M_CTC[4]

# Total hatchery releases from Quinsam and Campbell release sites, all release types
rel.x <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-10-08-SalmonRiverChinookReleases-AllYears.xlsx"),
  sheet = "Actual Release"
)

# Suggestion from Brendan Zoehner, SEP and Matt Clarke, StAD to include these
# July 2026
rel <- rel.x %>%
  filter(str_starts(RELEASE_SITE_NAME, "Salmon R/JNST") |
           str_starts(RELEASE_SITE_NAME, "Salmon R Up/JNST") |
           str_starts(RELEASE_SITE_NAME, "Pallans Cr")) %>%
  filter(RELEASE_STAGE_NAME %in% c("Smolt 0+", "Seapen 0+")) %>%
  summarise(n_rel = sum(TotalRelease), .by = c(RELEASE_YEAR)) %>%
  arrange(RELEASE_YEAR)

# Average hatchery releases from Q/C over last 5 years
hatch_rel <- rel %>%
  filter(RELEASE_YEAR %in% c(max(RELEASE_YEAR) - seq(4, 0))) %>%
  pull(n_rel) %>%
  mean()

## Strays removed
# stray <- c(0, 0, 7, 28,  17) #Annual strays,"45 (Age-2:0, Age 3: 7, Age 4: 28, Age 5: 17)",Hatchery,Weil et al. 2025; Proportions mean age-at-maturity Quinsam/Campbell 2018-2023,,,,,,,,,
# OM documentation says:
# @slot stray Matrix `[np, np]` where `np = length(Bio)` and row `p` indicates
# the re-assignment of hatchery fish to each population when they mature (at the
# recruitment life stage). For example,
# `SOM@stray <- matrix(c(0.75, 0.25, 0.25, 0.75), 2, 2)` indicates that 75 %
# of mature fish return to their natal river and 25 % stray in both populations.
# By default, an identity matrix is used (no straying).


# No Removals for other purposes (p_remove), e.g., No CWT sampling

h2 <- EnvStats::rnormTrunc(nsim, 0.25, 0.15, min = 0, max = 0.5)
Hatchery <- new(
  "Hatchery",
  n_r = n_r,
  n_yearling = hatch_rel, # Quinsam traditionals
  n_subyearling = 0,
  s_prespawn = 1, # Fecundity already adjusted for pre-spawn survival from M. Clarke life-cycle table, both hatchery and natural origin (M. Clarke DFO Science pers. comm.)
  s_egg_smolt = 0.75411, # B. Zoehner pers. comm. 16.21% mortality release-smolt x 20% mortality release-smolt: (1-0.1621)*(1 - 0.1) = 0.75411 surv
  s_egg_subyearling = 1,
  Mjuv_HOS = Mjuv_HOS,
  p_mature_HOS = p_mature, # Assume same maturity at NOS
  # stray_external = matrix(c(rep(0, 5), stray), maxage, 2),
  gamma = 0.8,  # HSRG standard, Sarita AHA inputs
  m = 1, # for purposes of premove_HOS for CWT sampling with is selective for marks (but premoveHOS = 0 here). Marking is also selective in f_brood - brood take rule is selective for marks
  pmax_esc = 0.33, # Ignored if using Brood rule in projection.R file #SEP guideline = 0.33. but removals are up to to 63% of total returns to Quinsam (1-esc$p_spawn),
  pmax_NOB = 1.0, # Ignored if using Brood rule in projection.R file #SEP guideline 0.5, suggested by Lian
  #f_brood = f_brood,  # Function defined in script 4
  ptarget_NOB = NA_real_,  # TBD
  phatchery = NA_real_, #proportion of escapement that actually spawn from input data to CM #CHECK # Stand-in for ESSR fishery with HOS exploitation rates of 0.5, 0.75, or 1
  hatchery_MSF = FALSE,
  premove_HOS = 0, #Fish removed for non-brood purposes, e.g. CWT. Alternative option in a function described in script 4.
  premove_NOS = 0,
  fec_brood = fec, #Assume the same age-specific fecundity in brood as in escapement
  fitness_type = c("Ford", "none"),
  theta = c(100, 80),
  rel_loss = rep(0, 3),
  zbar_start = c(90, 80),
  fitness_variance = 100,
  phenotype_variance = 10,
  heritability = h2,
  fitness_floor = 0.01
)

### Historical object ----

# Calculate Initial Njuv
# Sample abundance from final year with complete brood life cycle
nyears_cm <- d$Ldyr - 5

Njuv <- sapply(report[sim_samp], getElement, "N", simplify = "array") # year x age x origin x sim

# Matrix sims x ages x 1 (last dimension for number of life-history or hatchery release groups)
Njuv_NOS <- array(Njuv[nyears_cm, , 1, ], c(maxage, nsim, 1)) %>% aperm(c(2,1,3))
Njuv_HOS <- array(Njuv[nyears_cm, , 2, ], c(maxage, nsim, 1)) %>% aperm(c(2,1,3))



Historical <- new(
  "Historical",
  InitNjuv_NOS = Njuv_NOS,
  InitNjuv_HOS = Njuv_HOS
)



### Habitat ----

# megg is realized egg-smolt mortality rate, and exp(-megg) is egg-smolt survival
megg <- sapply(report[sim_samp], getElement, 'megg')
fry_sdev <- exp(-megg)

# # set.seed(234)
# # If I have fry_surv_mu and fry_surv_sd from another source (as in below)
# fry_surv_sim  <- matrix(rnorm(nsim * proyears, fry_surv_mu, fry_surv_sd), nsim, proyears)
# fry_surv_sim_trans <- plogis(fry_surv_sim)
#
# df <- exp(-megg) %>%
#   apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
#   reshape2::melt() %>%
#   mutate(Year = Var2 ) %>%
#   reshape2::dcast(list("Year", "Var1"), value.var = "value")
# df$'50%'[df$'50%' > 1] <- 0.99
# fpe <- df$'50%'
#
# fry_surv_mu <- df %>%
#   summarise(m = mean(qlogis(fpe)), sd = sd(qlogis(fpe))) %>%
#   pull(m) %>%
#   round(2)
#
# fry_surv_sd <- df %>%
#   summarise(m = mean(qlogis(fpe)), sd = sd(qlogis(fpe))) %>%
#   pull(sd) %>%
#   round(2)


Habitat <- new(
  "Habitat",
  use_habitat = FALSE, # Habitat optoin not used
  prespawn_rel = "BH",
  prespawn_prod = 1,# 5% pre-spawn mortality already accounted for in fecundity
  prespawn_capacity = Inf,
  fry_rel = "BH",
  fry_prod = 1,# To confirm with Quang that fry_surv_sim_trans will include average fw surv
  fry_capacity = Inf,
  fry_sdev = fry_sdev
)


SOM <- new("SOM",
           Name = som_name,
           nsim = nsim,
           proyears = proyears,
           seed = 1,
           Bio = Bio,
           Habitat = Habitat,
           Hatchery = Hatchery,
           Harvest = Harvest,
           Historical = Historical)
saveRDS(SOM, paste0("SOM/SOM_Salmon_", msurv, ".rds"))

