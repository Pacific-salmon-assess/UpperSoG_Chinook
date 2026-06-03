
library(tidyverse)
library(readxl)
library(salmonMSE)


#### Data ----
# Escapement time-series
pop <- "Salmon" #Campbell, Adam, Nimpkish, Salmon


esc_all <- readxl::read_excel(
  file.path("data", "SOG_N_Escapement-Salmon_Adam_Nimpkish.xlsx"),
  sheet = "Data") %>%
  filter(str_starts(Description, pop)) %>%
  rename(year = "Analysis Year") %>%
  # rename(escapement="Max Estimate") %>%
  # select (year, escapement)
  rename(esc.x = "Max Estimate") %>%
  mutate(nat_spawners =
           ifelse(is.na(`Natural Adult Spawners`), `Total Natural Spawners`, `Natural Adult Spawners`)) %>%
  summarise(
    escapement = sum(esc.x),
    nat_spawners = sum(nat_spawners),
    .by = c(year)
  ) %>%
  select (year, escapement, nat_spawners)

esc_2024 <-  readxl::read_excel(
  file.path("data", "PereboomA_20260409_091124 - Salmon and Adam 2024.xlsx"),
  sheet = "Data") %>%
  filter(str_starts(Description, pop)) %>%
  filter(Species == "Chinook") %>%
  rename(year = "Analysis Year") %>%
  # rename(escapement="Max Estimate") %>%
  # select (year, escapement)
  rename(esc.x = "Max Estimate") %>%
  mutate(nat_spawners =
           ifelse(is.na(`Natural Spawners Adult`), `Natural Spawners Total`, `Natural Spawners Adult`)) %>%
  summarise(
    escapement = sum(esc.x),
    nat_spawners = sum(nat_spawners),
    .by = c(year)
  ) %>%
  select (year, escapement, nat_spawners)

esc_2025 <- data.frame(year= 2025, escapement = (2532 + 63), nat_spawners=2532) # From Andrew Pereboom 9 April 2026
esc_all <- rbind(esc_all, esc_2024, esc_2025)


# Quinsam - CWT releases 1975-2025 (note, recoveries up to 2025).
# Remove releases prior to beginning of escapement (2002)

rel <-  readxl::read_excel(
  file.path("data", "Quinsam", "2025-07-23-Quinsam_Chinook_Releases_1970-2024.xlsx"),
  sheet = "Actual Release"
) %>%
  filter(RELEASE_YEAR >= min(esc_all$year))

# Quinsam - CWT releases of smolt0+ and seapen0+ only (traditionals).
cwt_rel <- rel %>%
  filter(RELEASE_STAGE_NAME %in% c("Smolt 0+", "Seapen 0+")) %>%
  summarise(n_CWT = sum(TaggedNum) - sum(ShedTagNum),
            .by = c(RELEASE_YEAR)) %>%
  arrange(RELEASE_YEAR)


# lastCWTreleaseYr <- max(cwt_rel %>% filter(n_CWT>0) %>% pull(RELEASE_YEAR))


# Releases aligned by BY
g <- ggplot(cwt_rel, aes(RELEASE_YEAR - 1, n_CWT)) +
  geom_point() +
  geom_line() +
  labs(x = "Brood Year", y = "Quinsam River CWT releases")



# Get tag codes for CWT releases, including smolt0+, seapen0+, excl. Fed Fry

cwt_rel_tags <- rel %>%
  filter(RELEASE_STAGE_NAME %in% c("Smolt 0+", "Seapen 0+")) %>%
  # filter(RELEASE_YEAR >= min(full_year) & RELEASE_YEAR <= max(full_year)) %>%
  summarise(n_CWT = sum(TaggedNum) - sum(ShedTagNum), .by = c(RELEASE_YEAR, MRP_TAGCODE)) %>%
  filter(n_CWT > 0) %>%
  arrange(RELEASE_YEAR) %>%
  rename(tag_code = "MRP_TAGCODE") %>%
  select(tag_code, RELEASE_YEAR)

cwt_rel_tags$tag_code <- as.numeric(cwt_rel_tags$tag_code)

# CWT recoveries
cwt_dat <- readr::read_csv("data/Quinsam/QUI_camp_recovery_wfisheries.csv") %>%
  rename(Age = "age")

# Only recoveries from release type smolt 0+ and seapen 0+ included by matching
# tag code in the recoveries with tagcodes from releases for smolts0+ and
# seapen 0+ (~300 tag_codes from releases removed)

# Used RELEASE_YEAR from release file (cwt_rel_tags) instead of
# RELEASE_YEAR = runyear - age + 1 in recoveries (cwt_dat), as there are some
# mistakes in aging

cwt_dat_subset <- inner_join(cwt_rel_tags, cwt_dat, by=c("tag_code"))



# Set up matrices
# We start the model at brood year max (min(cwt_dat_subset$RELEASE_YEAR), min(esc_all$year))
# Full matrix of ages (1-5) and years 2002 (earliest escapement) - 2025


full_matrix <- expand.grid(
  RELEASE_YEAR = (max(min(cwt_dat_subset$RELEASE_YEAR), min(esc_all$year)) ): 2025,
  Age = seq(1, 5)#6)
  # RS = c( "Seapen/Smolt 0+")
)
full_year <- data.frame(RELEASE_YEAR =
                          (max(min(cwt_dat_subset$RELEASE_YEAR), min(esc_all$year)) ):2025)


# Escapement CWT
# Two terminal fisheries are considered as 'escapement' for this population,
# as these fish are not vulnerable to these fisheries, so added here:
# "TGS FS" North Georgia Strait Freshwater Sport (6 tags)
# "TGEO ST TERM S" North Georgia Strait Terminal Sport (562 tags)

cwt_esc <- cwt_dat_subset %>%
  filter( (fishery_type == "escapement" &
             Coarse_description %in% c("Escapement", "Subsistence")) |
            fishery_era_name %in% c("TGS FS", "TGEO ST TERM S")) %>%
  summarise(n = sum(adjusted_estimated_number), .by = c(RELEASE_YEAR, Age)) %>%
  right_join(full_matrix, by = c("RELEASE_YEAR", "Age")) %>%
  reshape2::acast(list("RELEASE_YEAR", "Age"), value.var = "n", fill = 0)


# Preterminal CWT
# Three terminal fisheries are considered pre-terminal here:
# "TAK TERM T" Alaska Terminal Troll (15 tags)
# "TWCVI TERM N" Southwest WCVI Terminal Net (1 tag)
# "TNORTH FS" North Freshwater Sport (1 tag)

cwt_pt <- cwt_dat_subset %>%
  filter(fishery_type == "pre-terminal" |
           (fishery_type == "terminal" &
              fishery_era_name %in% c("TAK TERM T",
                                      "TWCVI TERM N",
                                      "TNORTH FS"))) %>%
  summarise(n = sum(adjusted_estimated_number), .by = c(RELEASE_YEAR, Age)) %>%
  right_join(full_matrix, by = c("RELEASE_YEAR", "Age")) %>%
  reshape2::acast(list("RELEASE_YEAR", "Age"), value.var = "n", fill = 0)


# Terminal CWTs
# Three terminal fisheries are considered pre-terminal here, so removed:
# "TAK TERM T" Alaska Terminal Troll (15 tags)
# "TWCVI TERM N" Southwest WCVI Terminal Net (1 tag)
# "TNORTH FS" North Freshwater Sport (1 tag)
# Two other pre-terminal fisheries are considered escapement, so removed:
# "TGS FS" North Georgia Strait Freshwater Sport (6 tags)
# "TGEO ST TERM S" North Georgia Strait Terminal Sport (562 tags)
cwt_t <- cwt_dat_subset %>%
  filter( fishery_type == "terminal" &
            !fishery_era_name %in% c("TAK TERM T",
                                     "TWCVI TERM N",
                                     "TNORTH FS",
                                     "TGS FS",
                                     "TGEO ST TERM S")) %>%
  summarise(n = sum(adjusted_estimated_number), .by = c(RELEASE_YEAR, Age)) %>%
  right_join(full_matrix, by = c("RELEASE_YEAR", "Age")) %>%
  reshape2::acast(list("RELEASE_YEAR", "Age"), value.var = "n", fill = 0)

# Note terminal fishery is:
# "TJNST TERM S" Johnstone Strait Terminal Sport (246 tags)

# Total hatchery releases from Salmon
rel_total.x <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-10-08-SalmonRiverChinookReleases-AllYears.xlsx"),
  sheet = "Actual Release"
)

rel_total <- rel_total.x %>%
  filter(str_starts(RELEASE_SITE_NAME, "Salmon R/JNST") |
           str_starts(RELEASE_SITE_NAME, "Salmon R Up/JNST")) %>%
  summarise(n_rel = sum(TotalRelease), .by = c(RELEASE_YEAR)) %>%
  arrange(RELEASE_YEAR)

rel_total <- left_join(full_year, rel_total, by = "RELEASE_YEAR")
rel_total$n_rel[is.na(rel_total$n_rel)] <- 0



esc <- esc_all %>%
  right_join(
    full_matrix %>% filter(Age == 1) %>% select(RELEASE_YEAR),
    by = c("year" = "RELEASE_YEAR")
  ) %>%
  arrange(year) %>%
  mutate(p_spawn = nat_spawners/escapement)
esc$p_spawn[is.na(esc$p_spawn)] <- na.omit(esc$p_spawn)[1]


# Data object for model
Ldyr <- nrow(cwt_esc)
Nages <- 5

mat <- c(0, 0.1, 0.4, 0.95, 1)#c(0, 0.01, 0.05, 0.2, 1) # Need to tune this vector for initial abundance #from WCVI = c(0, 0.1, 0.4, 0.95, 1)
vulPT <- c(0, 0.075, 0.9, 0.9, 1)
vulT <- vulPT

M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1)) # CTC 23-06 p.9; CWT Exploitation Rate analyses
M_CTC[1] <- 4 # Need to tune this value for initial abundance

# Fecundity eggs/adult spawner
fec_Salmon <- c(0, 0, 1174, 2936, 3669) # B. Zoehner, DFO, pers. comm. eggs/female = 5871, ppnal reductions by age from W&K (2024) and 50% female for age 4
fec_Quinsam <- c(0, 0, 800, 2000, 2500) # Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.


# Model assumption of catch expansion factor
# Use alternative values to change data weighting of CWT (re-adjust numbers accordingly)
cwtExp <- 1

d <- list(
  Nages = Nages,
  Ldyr = Ldyr,
  lht = 1,
  n_r = 1,
  s_enroute = 0.9,
  cwtrelease = as.vector(cwt_rel$n_CWT),
  cwtesc = array(round(cwt_esc/cwtExp), c(Ldyr, Nages, 1)),
  cwtcatPT = array(round(cwt_pt/cwtExp), c(Ldyr, Nages, 1)),
  cwtcatT = array(round(cwt_t/cwtExp), c(Ldyr, Nages, 1)),
  bvulPT = vulPT,
  bvulT = vulT,
  RelRegFPT = rep(1, Ldyr),
  RelRegFT = rep(1, Ldyr),
  bmatt = mat,
  mobase = M_CTC,
  hatchsurv =  0.8379,#B. Zoehner, DFO pers. comm., estimated 2017-2024
  gamma = 0.8,
  ssum = 1,
  fec = fec_Salmon*0.95,
  obsescape = esc$escapement,
  propwildspawn = round(esc$p_spawn, 2),
  pHOS_init = 0.18, # From SEP average over available time-series 2007-2022
  hatchrelease = c(rel_total$n_rel, rel_total$n_rel[length(rel_total$n_rel)]),
  finitPT = 0.4,
  finitT = 0.1,#,0.8,
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

start <- list(log_so = log(1 * max(d$obsescape, na.rm = TRUE)))

#### Fit with estimated productivity parameter (log_cr)
fit <- fit_CM(d, start = start, map = map, do_fit = TRUE, lower = list(moadd = -Inf))
samp <- sample_CM(fit, chains = 4, cores = 4, iter = 10000, thin = 5, seed = 1,
                  control=list(adapt_delta = 0.999,
                               stepsize = 0.01,
                               max_treedepth = 20))
saveRDS(samp, file = "CM/Salmon_06.03.26.rds")

samp <- readRDS(file = "CM/Salmon_06.03.26.rds")

year <- unique(full_matrix$RELEASE_YEAR)
rs_names <- c("Smolt 0+")
salmonMSE::report_CM(
  samp,
  rs_names = rs_names, name = "Salmon", year = year,
  dir = "CM", filename = "Salmon_06.03"
)

