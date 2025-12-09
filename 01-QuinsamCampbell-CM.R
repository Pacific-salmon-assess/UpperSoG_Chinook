
library(tidyverse)
library(readxl)
library(salmonMSE)

# Quinsam CWT recovery
rec <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-02-17-QuinsamChinook_Analyses_2005-2024.xlsx"),
  sheet = "Estimated"
) %>%
  mutate(is_catch = TotCatch > 0, is_esc = Escape > 0)


# CWT by release strategy, removing fed fry (only traditionals: seapen 0+ and smolt 0+)
cwt_rs <- rec %>%
  filter(RELEASE_STAGE_NAME %in% c("Seapen 0+", "Smolt 0+")) %>%
  # mutate(RS = "Seapen/Smolt 0+") %>%
  summarise(
    n_catch = sum(TotCatch),
    n_esc = sum(Escape),
    .by = c(Age, BROOD_YEAR)#, RS)
  )


# Quinsam - CWT releases
rel <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-02-17-QuinsamChinook_Analyses_2005-2024.xlsx"),
  sheet = "Releases"
)

rel_rs <- rel %>%
  filter(RELEASE_STAGE_NAME %in% c("Smolt 0+", "Seapen 0+")) %>%
  # mutate(RS = "Seapen/Smolt 0+") %>%
  summarise(n_CWT = sum(TaggedNum) - sum(ShedTagNum), .by = c(BROOD_YEAR))#, RS))

# Set up matrices
full_table <- expand.grid(
  BROOD_YEAR = seq(min(cwt_rs$BROOD_YEAR), 2023), # 2005 - 2023
  Age = seq(1, 5)#6)
  # RS = c( "Seapen/Smolt 0+")
) %>%
  left_join(cwt_rs)

cwt_catch <- reshape2::acast(full_table, list("BROOD_YEAR", "Age"), value.var = "n_catch", fill = 0)
cwt_esc <- reshape2::acast(full_table, list("BROOD_YEAR", "Age"), value.var = "n_esc", fill = 0)

cwt_rel <- left_join(
  full_table %>% filter(Age == 1) %>% select(BROOD_YEAR),# RS),
  rel_rs
) %>%
  reshape2::acast(list("BROOD_YEAR"), fill = 0)

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
  summarise(n_rel = sum(TotalRelease), .by = c(BROOD_YEAR)) %>%
  arrange(BROOD_YEAR)

# We start the model at brood year with CWT recoveries: min(cwt_rs$BROOD_YEAR) = 2005
# However, hatchery releases begin prior to that, so we want to have a spool up with a hatchery population
# going into 2005
hatch_init <- rel_Quinsam %>%
  filter(BROOD_YEAR %in% c(min(cwt_rs$BROOD_YEAR) - seq(5, 1))) %>%
  pull(n_rel) %>%
  mean()

# Crop to years with CWT data PLUS an extra year (to confirm)
full_year <- data.frame(BROOD_YEAR= seq(min(cwt_rs$BROOD_YEAR),
                                        max(cwt_rs$BROOD_YEAR) + 1))
rel_Quinsam <- left_join(full_year, rel_Quinsam, by = "BROOD_YEAR")
rel_Quinsam$n_rel[is.na(rel_Quinsam$n_rel)] <- 0

# Escapement time-series
pop <- c("Quinsam", "Campbell")


pop.cap <- str_to_upper(pop)
esc_all <- readxl::read_excel(
  file.path("data", "Quinsam", "fsar-sog-cn-cq-nuseds.xlsx"),
  sheet = "Data") %>%
  filter(StAD_Use == 1) %>%
  rename(WaterbodyName = "Waterbody Name") %>%
  filter(str_starts(WaterbodyName, pop.cap[1]) |
         str_starts(WaterbodyName, pop.cap[2])) %>%
  rename(year = "Analysis Year") %>%
  rename(esc.x = "Max Estimate") %>%
  summarise(
    escapement = sum(esc.x),
    nat_spawners = sum(`Natural Spawners Adult`),
    .by = c(year)
  ) %>%
  select (year, escapement, nat_spawners)

filter(esc_all, year %in% seq(2000, 2004, 1)) %>%
  pull(nat_spawners) %>%
  mean()
plot(nat_spawners ~ year, esc_all, typ = "o")

esc <- esc_all %>%
  right_join(
    full_table %>% filter(Age == 1) %>% select(BROOD_YEAR),
    by = c("year" = "BROOD_YEAR")
  ) %>%
  mutate(p_spawn = nat_spawners/escapement)

# pHOS data (use Quinsam as it's 5-10x larger than Campbell)
# 2024 value is much different compared to previous years!
pHOS_df_all <- readxl::read_excel(
  file.path("data", "2025-10-09-UpperSOGChinook-PNI.xlsx"),
  sheet = "Quinsam"
) %>%
  select(`Brood Year`, `pHOS`) %>%
  mutate(pHOS = as.numeric(pHOS))

#filter(pHOS_df, `Brood Year` %in% seq(2000, 2004, 1)) %>%
#  pull(pHOS) %>%
#  mean()
#plot(pHOS ~ `Brood Year`, pHOS_df, typ = "o")

pHOS_df <- pHOS_df_all %>%
  right_join(full_year, by = c("Brood Year" = "BROOD_YEAR"))

# Data object for model
Ldyr <- dim(cwt_esc)[1]
Nages <- 5#6

#mat <- c(0, 0.1, 0.4, 0.95, 1) # from WCVI = c(0, 0.1, 0.4, 0.95, 1)
mat <- c(0, 0.01, 0.1, 0.6, 1) # Need to tune this vector for initial abundance
vulPT <- c(0, 0.075, 0.9, 0.9, 1) #  from WCVI = c(0, 0.075, 0.9, 0.9, 1)
vulT <- rep(0, Nages)

M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1)) # CTC 23-06 p.9; CWT Exploitation Rate analyses
M_CTC[1] <- 5 # Need to tune this value for initial abundance

fec_Quinsam <- c(0, 0, 800, 2000, 2500) # Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.
# Eggs/total spawner (not female spawner)

# CWT expansion factor, the reciprocal is the likelihood weighting factor
# cwtExp of 0.1 implies weighting scalar of 10, cwtExp = 10 implies weighting scalar of 0.1
# Since cwt_esc and cwt_catch are fully expanded values, adjust the inputs matrices by the reciprocal of cwtExp
cwtExp <- 1

#pHOS_df$pHOS[1:5] <- NA_real_ # Potentially remove first 5 years of pHOS observations to allow model spool-up

d <- list(
  Nages = Nages,
  Ldyr = Ldyr,
  lht = 1,
  n_r = 1,
  cwtrelease = as.vector(cwt_rel),
  cwtesc = array(round(cwt_esc/cwtExp), c(Ldyr, Nages, 1)),
  cwtcatPT = array(round(cwt_catch/cwtExp), c(Ldyr, Nages, 1)),
  cwtcatT = NULL,
  bvulPT = vulPT,
  bvulT = vulT,
  RelRegFPT = rep(1, Ldyr),
  RelRegFT = rep(1, Ldyr),
  mobase = M_CTC,
  bmatt = mat,
  hatchsurv = 0.5, #Walters and Korman (2024); 1 used for WCVI Chinook
  pHOS_init = 0.75,
  spawn_init = 8400,
  gamma = 0.8,
  ssum = 1, # ppn female. Fecundity is eggs/total spawner, so this is set to 1.
  fec = fec_Quinsam,
  obsescape = esc$escapement,
  propwildspawn = esc$p_spawn, # This is the proportion of the natural spawners/return to river
  hatchrelease = rel_Quinsam$n_rel, #rep(0, Ldyr + 1),
  #obs_pHOS = pHOS_df$pHOS[1:Ldyr], # By brood year!
  #pHOS_sd = 0.5,
  finitPT = 0.8, # Walters and Korman (2024)
  finitT = 0,
  cwtExp = cwtExp
)

# Fix these parameters
map <- list()

# Fix maturity
#map$sd_matt <- factor(rep(NA, Nages-2)) # Not estimating year-specific maturity
#map$logit_matt <- factor(rep(NA, Ldyr * (Nages - 2)))

# Fix additional age-1 M
#map$moadd <- factor(NA)

# Fix age-1 density-independent M deviates
#map$wto <- factor(rep(NA, Ldyr))
#map$wto_sd <- factor(NA)

# Fix density dependent egg-smolt M deviates
#map$wt <- factor(rep(NA, Ldyr))
#map$wt_sd <- factor(NA)

# Fix observation error of Sarita escapement (needed, otherwise model can't separate process from obs error)
map$lnE_sd <- factor(NA)

start <- list(log_so = log(3 * max(d$obsescape)))

#### Fit with estimated productivity parameter (log_cr)
fit <- fit_CM(d, start = start, map = map, do_fit = TRUE, lower = list(moadd = -Inf))
samp <- sample_CM(fit, chains = 4, cores = 4, iter = 10000, thin = 5, seed = 1,
                  control=list(adapt_delta = 0.999,
                               stepsize = 0.01,
                               max_treedepth = 20))
saveRDS(samp, file = "CM/QuinsamCampbell_12.05.25.rds")

samp <- readRDS(file = "CM/QuinsamCampbell_12.05.25.rds")


#### Fit with log productivity = 1 ----
fit_p1 <- local({
  map$log_cr <- factor(NA) # Fixes the parameter
  start$log_cr <- 1        # Gives the starting value
  fit_CM(d, start = start, map = map, do_fit = TRUE, lower = list(moadd = -Inf))
})
samp_p1 <- sample_CM(fit_p1, chains = 4, cores = 4, iter = 10000, thin = 5, seed = 1,
                     control=list(adapt_delta = 0.999,
                                  stepsize = 0.01,
                                  max_treedepth = 20))
saveRDS(samp_p1, file = "CM/QuinsamCampbell_12.05.25_p1.rds")

#### Fit with higher older M for age 1+
fit2 <- local({
  d$covariate <- matrix(1, d$Ldyr, 1)
  fit_CM(d, start = start, map = map, do_fit = TRUE, lower = list(moadd = -Inf))
})
samp2 <- sample_CM(fit2, chains = 4, cores = 4, iter = 10000, thin = 5, seed = 1,
                   control=list(adapt_delta = 0.999,
                                stepsize = 0.01,
                                max_treedepth = 20))
saveRDS(samp2, file = "CM/QuinsamCampbell_12.05.25_age2M.rds")

if (FALSE) { # Diagnostic figures do not run when sourcing file

  #samp <- samp2

  # Check fits quickly
  report <- salmonMSE::get_report(samp)
  d <- salmonMSE::get_CMdata(samp@.MISC$CMfit)
  CM_fit_esc(report, d)
  CM_fit_pHOS(report, d) # Figure only when fitted to pHOS observations, otherwise nothing

  # Compare brood year pHOS when not fitted
  salmonMSE:::.CM_ts(report, year1 = pHOS_df$`Brood Year`[1], var = "pHOScensus_brood", ci = TRUE, ylab = "pHOScensus") +
    geom_point(data = pHOS_df, aes(x = `Brood Year`, y = pHOS)) +
    geom_line(data = pHOS_df, aes(x = `Brood Year`, y = pHOS), linetype = 3) +
    labs(x = "Brood Year")

  salmonMSE:::.CM_ts(report, year1 = pHOS_df$`Brood Year`[1], var = "pHOScensus", ci = TRUE, ylab = "pHOS")

  CM_F(report)
  CM_surv2(report) # Survival to age 2

  CM_maturity(report, d, brood = FALSE)
  CM_M(report)
  CM_SRR(report, year1 = pHOS_df$`Brood Year`[1])

  # Quickly check convergence
  CM_trace(samp)
  CM_pairs(samp, c("log_so", "log_cr"))

  # Launch full Stan app, this function frees up console whilst using the Shiny app
  launch_shinystan2 <- function(fit) {
    require(future)
    plan(multisession)
    future(
      shinystan::launch_shinystan(fit)
    )
  }
  launch_shinystan2(samp)
  #shinystan::launch_shinystan(samp)

  rs_names <- c("Smolt 0+")
  salmonMSE::report_CM(
    samp,
    rs_names = rs_names, name = "Quinsam/Campbell", year = unique(full_table$BROOD_YEAR),
    dir = "CM", filename = "QuinsamCampbell_12.05"
  )

  SMSY <- salmonMSE:::.CM_SMSY(report, d)
  Srep <- salmonMSE:::.CM_Srep(report, d)
  Sgen <- salmonMSE:::.CM_Sgen(report, d)

  range(Sgen/SMSY, na.rm = TRUE)
  range(SMSY/Srep, na.rm = TRUE)
  #If you remove the dot from the function name CM_SMSY instead of .CM_SMSY you will get the plotting function for the posterior

}
