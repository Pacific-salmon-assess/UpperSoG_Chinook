
# Make OM
library(tidyverse)
library(salmonMSE)

maxage <- 5
nsim <- 10 #100
# nyears <- 2
proyears <- 50
n_g <- 1

# Load exploitation rate model - Quinsam/Campbell
ERM_QuinsamCampbell <- readRDS("CM/QuinsamCampbell_02.03.26.rds")
report_QC <- salmonMSE:::get_report(ERM_QuinsamCampbell)

# Take maturity average from the 6 most recent brood years (2018-2024)
# Consider changing to earlier/longer period if recent ppns are less reliable
matt <- sapply(report_QC, getElement, "matt", simplify = "array") %>%
  aperm(c(1, 2, 4, 3))
matt_avg <-   apply(matt[seq(14, 19), , ,1], c(2,3), mean)
# apply (matt_avg,1,mean)

# Sarita code for refernce:
# matt_dev <- readRDS("CM/Sarita_maturity.rds")
# matt_avg <- sapply(matt_dev, function(x) {
#   apply(x$matt[seq(9, 14), , ], 2:3, mean)
# }, simplify = "array")
##apply(matt_avg, 1:2, mean)

set.seed(24)
sim_samp <- sample(seq(1, length(report_QC)), nsim)


### Natural mortality - NOS ----
Mjuv_NOS <- array(0, c(nsim, maxage-1, proyears, n_g))

# Survival from CTC 23-06 p. 9 for ages 2+
M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1))

# Age 1 value of Marine survival
# Option1 for natural populations:
# from life-cycle spreadsheet for natural-origin fish in Salmon River
# (M. Clarke pers .comm. 4 Dec 2025)
# this estimate is > than CTC estimate

# surv1 <- 1 - 0.991
# Mjuv_NOS[, 1, , 1] <- -log(surv1)
Mjuv_NOS[, 1, , 1] <- -log(0.009)#


# Option2 for QC:
# Draw random values
#
# df <- sapply(report_QC, getElement, "mo", simplify = "array") %>%
#   apply(1:2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
#   reshape2::melt() %>%
#   mutate(Year = Var2) %>%
#   mutate(Age = paste("Age", Var3)) %>%
#   reshape2::dcast(Year + Age ~ Var1, value.var = "value") %>%
#   filter(Age=='Age 1')
#
# s1_mu <- mean(qlogis(exp(-df$'50')))
# s1_sd <- sd(qlogis(exp(-df$'50%')))
#
# set.seed(234)
# s1_sim  <- matrix(rnorm(nsim * proyears, s1_mu, s1_sd), nsim, proyears)
# s1_sim_trans <- plogis(s1_sim)
# m1_sim <- -log(s1_sim_trans)
#
# Mjuv_NOS[, 1, , 1] <- m1_sim


# Age 2-5
Mjuv_NOS[, 2, , ] <- M_CTC[2]
Mjuv_NOS[, 3, , ] <- M_CTC[3]
Mjuv_NOS[, 4, , ] <- M_CTC[4]

### Maturity ----
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

fec_QC <- c(0, 0, 800, 2000, 2500) # Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.
p_female <- c(0, 0.01, 0.1, 0.55, 0.8) # Brown et al. in press, WCVI CK
# "The average percent of WCVI Chinook spawners that are female at each age is
# <1% at age two (called ‘jills’—i.e. female jacks), 10% female at age 3, 55%
# female at age 4, and 80% female at ages 5–7. These averages are based mostly
# on Robertson Creek Hatchery broodstock and Stamp River deadpitch sampling,
# but appear to be indicative of most WCVI Chinook populations." (p.26)

fec <- fec_QC * p_female

# Is Walters and Korman (2024) fecundity actually eggs/total spawner, so no need to multiply by p_female?


epro <- mean(sapply(report_QC, getElement, "epro"))
# Should epro be revised to be internally consistent with pars in the OM?

Bio <- new(
  "Bio",
  maxage = maxage,
  n_g = n_g,
  p_LHG = 1,
  p_mature = p_mature,
  # Freshwater survival and pre-spawn M included in habitat object
  # SRrel = "BH",
  # capacity = Inf,
  # kappa = (1 - 0.865) * epro, #0.865 = egg-smolt mortality from life-cycle table. Could use CM instead of life-cycle table estimate
  phi = epro,
  Mjuv_NOS = Mjuv_NOS,
  p_female = 1, # p_female specified in fecundity
  fec = fec,
  s_enroute = 0.855# covers return migration mortality and terminal fisheries (not captured with CWTs)
)

### Harvest, fishery vulnerability ----
# Extract FPT, only one life history group:
FPT <- sapply(report_QC, getElement, "FPT") %>% aperm(c(2,1))
UPT <- 1 - exp(-FPT)
# mean(UPT[14:19,])#Mean over most recent 6 years

vulPT <- sapply(report_QC[sim_samp], getElement, "vulPT")
vulT <- sapply(report_QC[sim_samp], getElement, "vulT")

# sapply(report_QC, getElement, "vulPT") %>% apply(1, mean)
# sapply(report_QC, getElement, "vulT") %>% apply(1, mean)


Harvest <- new(
  "Harvest",
  type_PT = "u",
  type_T = "u",
  u_preterminal = mean(UPT),
  u_terminal = 0,
  MSF_PT = FALSE,
  MSF_T = FALSE,
  release_mort = c(0.1,0.1),
  vulPT = t(vulPT),
  vulT = t(vulT)
)


### Hatchery ----
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

p_mature_RS <- array(0, c(nsim, maxage, proyears, n_r)) # Traditionals mature earlier
p_mature_RS[] <- array(matt_avg[, sim_samp], c(maxage, n_r, nsim,  proyears)) %>%
  aperm(c(3, 1, 4, 2))


# Total hatchery releases from Quinsam and Campbell release sites, all release types
rel_Quinsam.x <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-07-23-Quinsam_Chinook_Releases_1970-2024.xlsx"),
  sheet = "Actual Release"
)

# Suggestion from Brendan Zoehner, SEP, to include all Quinsam sites and Cold
# Creek (for Quinsam River) and Campbell sites, Elk R Chanel sites and Second Is
# (for Campbell River). Discovery Pass and Orange Pt are seapen releases, likely
# impacting both Quinsam and Campbell (8 Oct 2025)
rel_Quinsam <- rel_Quinsam.x %>%
  filter(str_starts(RELEASE_SITE_NAME, "Quinsam") |
           str_starts(RELEASE_SITE_NAME, "Cold") |
           str_starts(RELEASE_SITE_NAME, "Campbell") |
           str_starts(RELEASE_SITE_NAME, "Elk") |
           str_starts(RELEASE_SITE_NAME, "Second") |
           str_starts(RELEASE_SITE_NAME, "Discovery") |
           str_starts(RELEASE_SITE_NAME, "Orange")) %>%
  summarise(n_rel = sum(TotalRelease), .by = c(RELEASE_YEAR)) %>%
  arrange(RELEASE_YEAR)

# Average hatchery releases from Q/C over last 5 years
hatch_rel <- rel_Quinsam %>%
  filter(RELEASE_YEAR %in% c(max(RELEASE_YEAR) - seq(6, 1))) %>%
  pull(n_rel) %>%
  mean()

#Strays removed
stray <- c(0, 0, 7, 28,  17) #Annual strays,"45 (Age-2:0, Age 3: 7, Age 4: 28, Age 5: 17)",Hatchery,Weil et al. 2025; Proportions mean age-at-maturity Quinsam/Campbell 2018-2023,,,,,,,,,
# OM documenation says:
# @slot stray Matrix `[np, np]` where `np = length(Bio)` and row `p` indicates
# the re-assignment of hatchery fish to each population when they mature (at the
# recruitment life stage). For example,
# `SOM@stray <- matrix(c(0.75, 0.25, 0.25, 0.75), 2, 2)` indicates that 75 %
# of mature fish return to their natal river and 25 % stray in both populations.
# By default, an identity matrix is used (no straying).


# Removals for other purposes (p_remove), e.g., CWT sampling
# Non-brood removals are primarily HO fish. Assume 100% HOS
# Non-brood removals are also selective for older male fish, but this is not accounted for
rem_Quinsam.x <- readxl::read_excel(
  file.path("data", "Quinsam", "2026-02-11-QuinsamRiverChinook-RiverReturns-1981-2023.xlsx"),
  sheet = "River Returns"
)
names(rem_Quinsam.x)[names(rem_Quinsam.x)=="Sex/Maturity"] <- "SexMaturity"

# Sum removals for non-brood purposes (Brendan Zoehner pers. comm. 12 Feb 26)
rem_Quinsam <- rem_Quinsam.x %>%
  filter(SexMaturity == "Male"| SexMaturity == "Female") %>%
  mutate(OtherRemovals = rowSums(across(c("04Excess To Spawning Req (ESSR)",
                                          "05Sold",
                                          "06Holding Mortalities",
                                          "07Used For Other Purposes/Other Mortalities")),
                                 na.rm=TRUE)) %>%
  rename("year" = "Recovery Year") %>%
  group_by(year) %>% summarise(OtherRemovals_Sum= sum(OtherRemovals, na.rm=T))

# What ppn are non-brood removals to natural spawners?
# (assuming taken after brood removals and enroute M, so ppn is the ppn of removals to natural spawners not total escapaement)
# First get total escapement (used in CM), and 'natural spawners' which account for all removals
esc_all <- readxl::read_excel(
  file.path("data", "Quinsam", "fsar-sog-cn-cq-nuseds.xlsx"),
  sheet = "Data") %>%
  filter(StAD_Use == 1) %>%
  rename(WaterbodyName = "Waterbody Name") %>%
  filter(str_starts(WaterbodyName, pop.cap[1]) |
           str_starts(WaterbodyName, pop.cap[2])) %>%
  rename(year = "Analysis Year") %>%
  rename(esc.x = "Max Estimate") %>%
  mutate(nat_spawners =
           ifelse(is.na(`Natural Spawners Adult`), `Natural Spawners Total`, `Natural Spawners Adult`)) %>%
  summarise(
    escapement = sum(esc.x),
    nat_spawners = sum(nat_spawners),
    .by = c(year)
  ) %>%
  select (year, escapement, nat_spawners)
# Then combine with removals
esc_rem <- left_join(rem_Quinsam, esc_all) %>%
  mutate(nat_spawners_wOR = OtherRemovals_Sum + nat_spawners) %>%
  mutate(ppnOR = OtherRemovals_Sum/nat_spawners_wOR)
ppnRemoveHOS <- mean(esc_rem$ppnOR[(length(esc_rem$ppnOR-6)):length(esc_rem$ppnOR)])
# This is too high, as ppn should be ppn of HO spawners, not total spawners
# Get HO spawners from CM and repeat?


h2 <- EnvStats::rnormTrunc(nsim, 0.25, 0.15, min = 0, max = 0.5)
Hatchery <- new(
  "Hatchery",
  n_r = n_r,
  n_yearling = hatch_rel, # Quinsam traditionals
  n_subyearling = 0,
  s_prespawn = 0.95,  # M. Clarke DFO Science pers. comm.
  s_egg_smolt = 0.68, #Assumed 20% M shortly after release to match CM and life-cycle table (0.8). Could lower to account for incubation mortality = 0.68
  s_egg_subyearling = 1,
  Mjuv_HOS = Mjuv_HOS,
  p_mature_HOS = p_mature_RS,
  # stray_external = matrix(c(rep(0, 5), stray), maxage, 2),
  gamma = 0.8,  # HSRG standard, Sarita AHA inputs
  m = 1, # for purposes of premove_HOS. A different assumption (m=0) is made for f_brood - brood take rule
  pmax_esc = 0.7, #SEP guideline = 0.33. but removals are up to to 63% of total returns to Quinsam (1-esc$p_spawn),
  pmax_NOB = 1.0, #SEP guideline 0.5, suggested by Lian #Brood rule in projection.R file
  #f_brood = f_brood,  # Function defined in script 4
  ptarget_NOB = 0,  # TBD
  phatchery = NA_real_, #proportion of escapement that actually spawn from input data to CM #CHECK # Stand-in for ESSR fishery with HOS exploitation rates of 0.5, 0.75, or 1
  hatchery_MSF = FALSE,
  premove_HOS = 0,#0.05#ppnRemoveHOS,
  premove_NOS = 0,
  fec_brood = fec,
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


# Calculate HistNjuv_NOS (which return year?)
# Sample abundance from final year with complete brood life cycle
nyears_cm <- dim(sapply(report_QC, getElement, 'egg'))[1] - 5

NOS <- sapply(report_QC[sim_samp], function(x) x$syear[nyears_cm, , 1])
HOS <- sapply(report_QC[sim_samp], function(x) x$syear[nyears_cm, , 2])

NOS_g <- sapply(1:Bio@n_g, function(g) {
  1.0 * NOS
}, simplify = "array") %>%
   aperm(c(2, 1, 3))

HOS_r <- sapply(1:Hatchery@n_r, function(r) {
  1.0 * HOS
}, simplify = "array") %>%
  aperm(c(2, 1, 3))

# Assign 50% of N juveniles to HOS and 50% fo NOS
Njuv <- sapply(report_QC[sim_samp], getElement, "N", simplify = "array") # year x age x origin x sim
Njuv_NOS <- sapply(1:Bio@n_g, function(g) {
  N <- 0.5 * array(Njuv[nyears_cm, , 1, ], c(maxage, nsim))
  return(N)
}, simplify = 'array') %>%
  aperm(c(2, 1, 3))

Njuv_HOS <- sapply(1:Hatchery@n_r, function(r) {
  N <- 0.5 * array(Njuv[nyears_cm, , 2, ], c(maxage, nsim))
  return(N)
}, simplify = 'array') %>%
  aperm(c(2, 1, 3))


Historical <- new(
  "Historical",
  #InitNOS = NOS_g, # Remove for salmonMSE 1.0.0
  #InitHOS = HOS_r,
  InitNjuv_NOS = Njuv_NOS,
  InitNjuv_HOS = Njuv_HOS
)



### Habitat ----


# fry_surv <- read.csv("data/Sarita/fry_surv.csv")
# fry_surv_year <- read.csv("data/Sarita/fry_surv_year.csv")
# instead of sarita fry, take time-seies of median fw surv (exp(-megg) from CM, 'CMoutputs.R')

# megg is realized egg-smolt mortality rate, and exp(-megg) is egg-smolt survival
megg <- sapply(report_QC, getElement, 'megg')
df <- exp(-megg) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 ) %>%
  reshape2::dcast(list("Year", "Var1"), value.var = "value")
df$'50%'[df$'50%' > 1] <- 0.99
# fpe <- df$'50%'[(length(df$'50%')-20):length(df$'50%')]
fpe <- df$'50%'

fry_surv_mu <- df %>%
  summarise(m = mean(qlogis(fpe)), sd = sd(qlogis(fpe))) %>%
  pull(m) %>%
  round(2)

fry_surv_sd <- df %>%
  summarise(m = mean(qlogis(fpe)), sd = sd(qlogis(fpe))) %>%
  pull(sd) %>%
  round(2)

set.seed(234)
# Remove interannual variabitly in fw survival for now
# fry_surv_sim  <- matrix(rnorm(nsim * proyears, fry_surv_mu, fry_surv_sd), nsim, proyears)
fry_surv_sim  <- matrix(rnorm(nsim * proyears, fry_surv_mu, 0), nsim, proyears)
fry_surv_sim_trans <- plogis(fry_surv_sim)

Habitat <- new(
  "Habitat",
  use_habitat = TRUE,
  prespawn_rel = "BH",
  prespawn_prod = 0.95,# Adding 5% pre-spawn mortality
  prespawn_capacity = Inf,
  fry_rel = "BH",
  fry_prod = 1,# To confirm with Quang that fry_surv_sim_trans will include average fw surv
  fry_capacity = Inf,
  fry_sdev = fry_surv_sim_trans
)


SOM <- new("SOM",
           Name = "QC base",
           nsim = nsim,
           proyears = proyears,
           seed = 1,
           Bio = Bio,
           Habitat = Habitat,
           Hatchery = Hatchery,
           Harvest = Harvest,
           Historical = Historical)
saveRDS(SOM, "SOM/SOM_base.rds")

