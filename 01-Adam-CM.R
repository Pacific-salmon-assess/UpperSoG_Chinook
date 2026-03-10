
library(tidyverse)
library(readxl)
library(salmonMSE)

# Quinsam CWT recovery

rec <- local({
  rec1 <- readxl::read_excel(
    file.path("data", "Quinsam", "2025-02-17-QuinsamChinook_Analyses_2005-2024.xlsx"),
    sheet = "Estimated"
  ) %>%
    select(RELEASE_STAGE_NAME, TotCatch, Escape, Age, BROOD_YEAR, RELEASE_YEAR, RecovYear)

  # rec2 <- readxl::read_excel(
  #   file.path("data", "Quinsam", "2025-12-18-QuinsamCN-CatchbyStat-1980-2005.xlsx"),
  #   sheet = "Estimated"
  # ) %>%
  #   select(RELEASE_STAGE_NAME, TotCatch, Escape, Age, BROOD_YEAR, RELEASE_YEAR, RecovYear)

  rbind(rec1)#, rec2)
}) %>%
  mutate(is_catch = TotCatch > 0, is_esc = Escape > 0)

# CWT by release strategy, removing fed fry (only traditionals: seapen 0+ and smolt 0+)
cwt_rs <- rec %>%
  filter(RELEASE_STAGE_NAME %in% c("Seapen 0+", "Smolt 0+")) %>%
  # mutate(RS = "Seapen/Smolt 0+") %>%
  summarise(
    n_catch = sum(TotCatch),
    n_esc = sum(Escape),
    .by = c(Age, BROOD_YEAR)#, RS)
  ) %>%
  rename("RELEASE_YEAR"="BROOD_YEAR")


# Quinsam - CWT releases
rel <- local({

  rel1 <- readxl::read_excel(
    file.path("data", "Quinsam", "2025-02-17-QuinsamChinook_Analyses_2005-2024.xlsx"),
    sheet = "Releases"
  )

  # rel2 <- readxl::read_excel(
  #   file.path("data", "Quinsam", "2025-12-18-QuinsamCN-Releases-1980-2005.xlsx"),
  #   sheet = "Actual Release"
  # )

  rbind(rel1)#, rel2)
})


rel_rs <- rel %>%
  filter(RELEASE_STAGE_NAME %in% c("Smolt 0+", "Seapen 0+")) %>%
  # mutate(RS = "Seapen/Smolt 0+") %>%
  summarise(n_CWT = sum(TaggedNum) - sum(ShedTagNum), .by = c(RELEASE_YEAR)) %>%
  arrange(RELEASE_YEAR)

# Set up matrices
full_table <- expand.grid(
  RELEASE_YEAR = seq(min(cwt_rs$RELEASE_YEAR), 2023), # 2005 - 2023
  Age = seq(1, 5)#6)
  # RS = c( "Seapen/Smolt 0+")
) %>%
  left_join(cwt_rs)

cwt_catch <- reshape2::acast(full_table, list("RELEASE_YEAR", "Age"), value.var = "n_catch", fill = 0)
cwt_esc <- reshape2::acast(full_table, list("RELEASE_YEAR", "Age"), value.var = "n_esc", fill = 0)

cwt_rel <- left_join(
  full_table %>% filter(Age == 1) %>% select(RELEASE_YEAR),
  rel_rs
) %>%
  reshape2::acast(list("RELEASE_YEAR"), fill = 0)


# No hatchery releases from Adam River

# Escapement time-series
pop <- "Adam" #Campbell, Adam, Nimpkish, Salmon

if(pop %in% c("Adam", "Nimpkish", "Salmon")) {
  esc <- readxl::read_excel(
    file.path("data", "SOG_N_Escapement-Salmon_Adam_Nimpkish.xlsx"),
    sheet = "Data") %>%
    filter(str_starts(Description, pop)) %>%
    rename(year = "Analysis Year") %>%
    rename(escapement="Max Estimate") %>%
    select (year, escapement) %>%
    right_join(
      full_table %>% filter(Age == 1) %>% select(RELEASE_YEAR),
      by = c("year" = "RELEASE_YEAR")
    )
}


# Data object for model
Ldyr <- dim(cwt_esc)[1]
Nages <- 5#6

mat <- c(0, 0.01, 0.05, 0.2, 1) # Need to tune this vector for initial abundance #from WCVI = c(0, 0.1, 0.4, 0.95, 1)
vulPT <- c(0, 0.075, 0.9, 0.9, 1) #  from WCVI = c(0, 0.075, 0.9, 0.9, 1)
vulT <- rep(0, Nages)

M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1)) # CTC 23-06 p.9; CWT Exploitation Rate analyses
M_CTC[1] <- 4 # Need to tune this value for initial abundance

fec_Quinsam <- c(0, 0, 800, 2000, 2500) # Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.
# Eggs/female spawner?
p_female <- c(0, 0.01, 0.1, 0.55, 0.8) # Brown et al. in press, WCVI CK
# "The average percent of WCVI Chinook spawners that are female at each age is
# <1% at age two (called ‘jills’—i.e. female jacks), 10% female at age 3, 55%
# female at age 4, and 80% female at ages 5–7. These averages are based mostly
# on Robertson Creek Hatchery broodstock and Stamp River deadpitch sampling,
# but appear to be indicative of most WCVI Chinook populations." (p.26)
fec_Quinsam <- fec_Quinsam * p_female

d <- list(
  Nages = Nages,
  Ldyr = Ldyr,
  lht = 1,
  n_r = 1,
  s_enroute = 0.855, # From  M. Clarke life-cycle table = 10% return migration M followed by 20% terminal ER
  cwtrelease = as.vector(cwt_rel),
  cwtesc = array(round(cwt_esc), c(Ldyr, Nages, 1)),
  cwtcatPT = array(round(cwt_catch), c(Ldyr, Nages, 1)),
  cwtcatT = NULL,
  bvulPT = vulPT,
  bvulT = vulT,
  RelRegFPT = rep(1, Ldyr),
  RelRegFT = rep(1, Ldyr),
  mobase = M_CTC,
  bmatt = mat,
  hatchsurv = 0.8, #From M. Clarke life-cycle table. Walters and Korman (2024) used 0.5; 1 used for WCVI Chinook
  pHOS_init = 0,
  spawn_init = 67,
  gamma = 0.8,
  ssum = 1, # ppn female. Fecundity is eggs/total spawner, so this is set to 1.
  fec = fec_Quinsam*0.95,
  obsescape = esc$escapement,
  propwildspawn = rep(1, Ldyr),
  hatchrelease = rep(0, Ldyr + 1),
  finitPT = 0.8, # Walters and Korman (2024)
  finitT = 0,  # Walters and Korman (2024)
  cwtExp = 1 # Sarita used 1 #Walters and Korman (2024) used 0.1
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

start <- list(log_so = log(1.5 * max(d$obsescape)))

# Fit with sampling rate = 1
fit <- fit_CM(d, start = start, map = map, do_fit = TRUE)
samp <- sample_CM(fit, chains = 4, cores = 4, iter = 10000, thin = 5, seed=1,
                  control=list(adapt_delta = 0.999,
                               stepsize = 0.01,
                               max_treedepth = 20))
saveRDS(samp, file = "CM/Adam_03.07.26.rds")

samp <- readRDS(file = "CM/Adam_03.07.26.rds")
report <- salmonMSE:::get_report(samp)
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
#shinystan::launch_shinystan(samp)

rs_names <- c("Smolt 0+")
salmonMSE::report_CM(
  samp,
  rs_names = rs_names, name = "Adam", year = unique(full_table$RELEASE_YEAR),
  dir = "CM", filename = "Adam_03.07"
)

SMSY <- salmonMSE:::.CM_SMSY(report, d)
Srep <- salmonMSE:::.CM_Srep(report, d)
Sgen <- salmonMSE:::.CM_Sgen(report, d)
# alpha <- sapply(report, getElement, "alpha"), but need epro to get UMSY (alpha*epro gives full life cycle prod)
calc_Umsy_Ricker(log(2.369))
