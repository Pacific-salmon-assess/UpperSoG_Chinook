
library(tidyverse)
library(readxl)

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
  # Add CWT catches from PFMA 13 and 15 to escapement
  mutate(TotEscSum = rowSums(select(., starts_with("3"),
                                    starts_with('1PFMA 13'),
                                    starts_with('1PFMA 15')), na.rm=T)) %>%
  # Remove CWT catchs from PFMA 13 and 15 to catches
  select(!starts_with('1PFMA 13')) %>%
  select(!starts_with('1PFMA 15')) %>%
  mutate(TotCatchSum = rowSums(select(., starts_with("1"), starts_with("2")), na.rm=T)) %>%
  summarise(
    n_esc = sum(TotEscSum),
    n_catch = sum(TotCatchSum),
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

# Total hatchery releases from Nimpkish/Woss release sites, all release types
rel_total.x <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-07-23-NimpkishWoss_Chinook_Releases_1970-2024.xlsx"),
  sheet = "Actual Release"
)


rel_total <- rel_total.x %>%
  filter(str_starts(RELEASE_SITE_NAME, "Woss R") |
           str_starts(RELEASE_SITE_NAME, "Woss Lk") |
           str_starts(RELEASE_SITE_NAME, "Nimpkish R") |
           str_starts(RELEASE_SITE_NAME, "Nimpkish R Up")) %>%
  summarise(n_rel = sum(TotalRelease), .by = c(BROOD_YEAR)) %>%
  arrange(BROOD_YEAR)

# Crop to years with CWT data PLUS an extra year (to confirm)
full_year <- data.frame(BROOD_YEAR= seq(min(cwt_rs$BROOD_YEAR),
                                        max(cwt_rs$BROOD_YEAR) + 1))
rel_total <- left_join(full_year, rel_total, by = "BROOD_YEAR")
rel_total$n_rel[is.na(rel_total$n_rel)] <- 0

# Escapement time-series
pop <- "Nimpkish" #Campbell, Adam, Nimpkish, Salmon

if(pop %in% c("Adam", "Nimpkish", "Salmon")) {
  esc <- readxl::read_excel(
    file.path("data", "SOG_N_Escapement-Salmon_Adam_Nimpkish.xlsx"),
    sheet = "Data") %>%
    filter(str_starts(Description, pop)) %>%
    rename(year = "Analysis Year") %>%
    rename(escapement="Max Estimate") %>%
    select (year, escapement) %>%
    right_join(
      full_table %>% filter(Age == 1) %>% select(BROOD_YEAR),
      by = c("year" = "BROOD_YEAR")
    ) %>%
    arrange(year)
}


# Data object for model
Ldyr <- dim(cwt_esc)[1]
Nages <- 5#6

mat <- c(0, 0.1, 0.4, 0.95, 1) # from WCVI = c(0, 0.1, 0.4, 0.95, 1)
vulPT <- c(0, 0.075, 0.9, 0.9, 1) #  from WCVI = c(0, 0.075, 0.9, 0.9, 1)
vulT <- rep(0, Nages)

M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1)) # CTC 23-06 p.9; CWT Exploitation Rate analyses

fec_Quinsam <- c(0, 0, 800, 2000, 2500) # Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.
# Eggs/total spawner (not female spawner)

d <- list(
  Nages = Nages,
  Ldyr = Ldyr,
  lht = 1,
  n_r = 1,
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
  pHOS_init = 0.13,
  spawn_init = 500/0.4, #Expanded to account for Woss:Nimpkish ratio, M. Clarke pers. comm. 21 Nov 2025
  gamma = 0.8,
  ssum = 1, # ppn female. Fecundity is eggs/total spawner, so this is set to 1.
  fec = fec_Quinsam,
  obsescape = esc$escapement/0.4,#Expanded to account for Woss:Nimpkish ratio, M. Clarke pers. comm. 21 Nov 2025
  propwildspawn = rep(1, Ldyr),
  hatchrelease = rel_total$n_rel, #rep(0, Ldyr + 1),
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

start <- list(log_so = log(3 * max(d$obsescape, na.rm=TRUE)))

# Fit with sampling rate = 1
fit <- fit_CM(d, start = start, map = map, do_fit = TRUE)
samp <- sample_CM(fit, chains = 4, cores = 4, iter = 10000, thin = 5,
                  control=list(adapt_delta = 0.999, stepsize = 0.01,
                               max_treedepth = 20))
saveRDS(samp, file = "CM/Woss_12.14.25.rds")

samp <- readRDS(file = "CM/Woss_12.14.25.rds")
report <- salmonMSE:::get_report(samp)
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
#shinystan::launch_shinystan(samp)

rs_names <- c("Smolt 0+")
salmonMSE::report_CM(
  samp,
  rs_names = rs_names, name = "Woss", year = unique(full_table$BROOD_YEAR),
  dir = "CM", filename = "Woss_12.14"
)
