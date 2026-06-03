library(tidyverse)
library(readxl)
library(salmonMSE)


#### Data ----
# Escapement time-series 1984-2025
pop <- c("Quinsam", "Campbell")

pop.cap <- str_to_upper(pop)
esc_all <- readxl::read_excel(
  # file.path("data", "Quinsam", "fsar-sog-cn-cq-nuseds.xlsx"),
  file.path("data", "Quinsam", "fsar-sog-cn-cq-nuseds_2025.xlsx"),
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





## pHOS data (use Quinsam as it's 5-10x larger than Campbell)
## 2024 value is much different compared to previous years!
pHOS_df_all <- readxl::read_excel(
  file.path("data", "2025-10-09-UpperSOGChinook-PNI.xlsx"),
  sheet = "Quinsam"
) %>%
  mutate(Release_Year = `Brood Year` + 1) %>%
  select(Release_Year, `pHOS`) %>%
  mutate(pHOS = as.numeric(pHOS))


# pHOS_df <- pHOS_df_all %>%
#   right_join(full_year, by = c("Release_Year" = "RELEASE_YEAR")) %>%
#   filter(Release_Year <= 2025)


# Quinsam - CWT releases 1975-2025 (note, recoveries up to 2025).
# Remove releases prior to 1984, earliest escapaement year

rel <-  readxl::read_excel(
  file.path("data", "Quinsam", "2025-07-23-Quinsam_Chinook_Releases_1970-2024.xlsx"),
  sheet = "Actual Release"
) %>%
  filter(RELEASE_YEAR >= min(esc_all$year))

# Quinsam - CWT releases of smolt0+ and seapen0+ only (traditionals).
# Removing fed fry
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
# Full matrix of ages (1-5) and years 1984 (earliest escapement) - 2025

full_matrix <- expand.grid(
  RELEASE_YEAR = (max(min(cwt_dat_subset$RELEASE_YEAR), min(esc_all$year)) ): 2025, #max(cwt_dat_subset$RELEASE_YEAR),# 1980 - 2025
  Age = seq(1, 5)#6)
  # RS = c( "Seapen/Smolt 0+")
)
full_year <- data.frame(RELEASE_YEAR =
                          (max(min(cwt_dat_subset$RELEASE_YEAR), min(esc_all$year)) ):2025)

# Shortened full_matrix to 1984 (min esc year)-2025, to cover only time-period with good escapement (unlike Walters and Korman 2024)


cwt_esc <- cwt_dat_subset %>%
  filter(fishery_type == "escapement",
         Coarse_description %in% c("Escapement", "Subsistence")) %>%
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
# Three terminal fisheries are considered pre-terminal here:
# "TAK TERM T" Alaska Terminal Troll (15 tags)
# "TWCVI TERM N" Southwest WCVI Terminal Net (1 tag)
# "TNORTH FS" North Freshwater Sport (1 tag)
cwt_t <- cwt_dat_subset %>%
  filter(fishery_type == "terminal" &
           !fishery_era_name %in% c("TAK TERM T", "TWCVI TERM N", "TNORTH FS")) %>%
  summarise(n = sum(adjusted_estimated_number), .by = c(RELEASE_YEAR, Age)) %>%
  right_join(full_matrix, by = c("RELEASE_YEAR", "Age")) %>%
  reshape2::acast(list("RELEASE_YEAR", "Age"), value.var = "n", fill = 0)

# Note terminal fisheries are:
# "TGS FS" North Georgia Strait Freshwater Sport (6 tags)
# "TGEO ST TERM S" North Georgia Strait Terminal Sport (562 tags)
# "TJNST TERM S" Johnstone Strait Terminal Sport (246 tags)

esc <- esc_all %>%
  right_join(
    full_matrix %>% filter(Age == 1) %>% select(RELEASE_YEAR),
    by = c("year" = "RELEASE_YEAR")
  ) %>%
  arrange(year) %>%
  mutate(p_spawn = nat_spawners/escapement)
esc$p_spawn[is.na(esc$p_spawn)] <- na.omit(esc$p_spawn)[1]


# Total hatchery releases from Quinsam and Campbell release sites, all release types

rel_Quinsam.x <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-07-23-Quinsam_Chinook_Releases_1970-2024.xlsx"),
  # File includes 2025
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

rel_Quinsam <- left_join(full_year, rel_Quinsam, by = "RELEASE_YEAR")
rel_Quinsam$n_rel[is.na(rel_Quinsam$n_rel)] <- 0

# Data object for model
Ldyr <- dim(cwt_esc)[1]
Nages <- 5#6

#mat <- c(0, 0.1, 0.4, 0.95, 1) # from WCVI = c(0, 0.1, 0.4, 0.95, 1)
mat <- c(0, 0.01, 0.05, 0.2, 1) # Need to tune this vector for initial abundance
vulPT <- c(0, 0.075, 0.9, 0.9, 1) #  from WCVI = c(0, 0.075, 0.9, 0.9, 1)
vulT <- vulPT#rep(0, Nages)

M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1)) # CTC 23-06 p.9; CWT Exploitation Rate analyses
M_CTC[1] <- 4 # Need to tune this value for initial abundance

covariate1 <- readxl::read_excel(
  file.path("data", "Quinsam", "covariate1.xlsx"),
  sheet = "covariate1"
)
covariate1 <- as.matrix(covariate1)

fec_Q <- c(0, 0, 800, 2000, 2500) # Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.
# Eggs/total spawner (not female spawner)
# p_female <- c(0, 0.01, 0.1, 0.55, 0.8) # Brown et al. in press, WCVI CK
# "The average percent of WCVI Chinook spawners that are female at each age is
# <1% at age two (called ‘jills’—i.e. female jacks), 10% female at age 3, 55%
# female at age 4, and 80% female at ages 5–7. These averages are based mostly
# on Robertson Creek Hatchery broodstock and Stamp River deadpitch sampling,
# but appear to be indicative of most WCVI Chinook populations." (p.26)

fec_Quinsam <- fec_Q #* p_female
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
  s_enroute = 0.9, # From  M. Clarke life-cycle table = 10% return migration M followed by 15% terminal ER #0.765
  cwtrelease = as.vector(cwt_rel$n_CWT),#Aligned by RY
  cwtesc = array(round(cwt_esc/cwtExp), c(Ldyr, Nages, 1)), #Algined by RY
  cwtcatPT = array(round(cwt_pt/cwtExp), c(Ldyr, Nages, 1)), #Algined by RY
  cwtcatT =  array(round(cwt_t/cwtExp), c(Ldyr, Nages, 1)), #Algined by RY
  bvulPT = vulPT,
  bvulT = vulT,
  RelRegFPT = rep(1, Ldyr),
  RelRegFT = rep(1, Ldyr),
  mobase = M_CTC,
  # covariate1 = covariate1,
  bmatt = mat,
  hatchsurv = 0.9, #Increased from 0.8 for Salmon River (M. Clarke life-cycle table) due to shorter river migration
  pHOS_init = 0.75,
  spawn_init = 4000,
  gamma = 0.8,
  ssum = 1, # ppn female. Fecundity is eggs/total spawner, so this is set to 1.
  fec = fec_Quinsam *0.95, # adjusting for pre-spawn survival from M. Clarke life-cycle table
  obsescape = esc$escapement, #Algined by BY
  propwildspawn = esc$p_spawn, # This is the proportion of the natural spawners/return to river
  hatchrelease = c(rel_Quinsam$n_rel, rel_Quinsam$n_rel[length(rel_Quinsam$n_rel)]), #rep(0, Ldyr + 1),#Algined by RY
  # obs_pHOS = pHOS_df$pHOS[1:Ldyr], # By brood year!
  # pHOS_sd = 0.5,
  finitPT = 0.8, # Walters and Korman (2024)
  finitT = 0.2,
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

start <- list(log_so = log(3 * max(d$obsescape, na.rm = TRUE)))

#### Fit with estimated productivity parameter (log_cr)
fit <- fit_CM(d, start = start, map = map, do_fit = TRUE, lower = list(moadd = -Inf))
samp <- sample_CM(fit, chains = 4, cores = 4, iter = 10000, thin = 5, seed = 1,
                  control=list(adapt_delta = 0.999,
                               stepsize = 0.01,
                               max_treedepth = 20))
saveRDS(samp, file = "CM/QuinsamCampbell_05.29.26.rds")

samp <- readRDS(file = "CM/QuinsamCampbell_05.29.26.rds")

year <- unique(full_matrix$RELEASE_YEAR)
rs_names <- c("Smolt 0+")
salmonMSE::report_CM(
  samp,
  rs_names = rs_names, name = "Quinsam/Campbell", year = year,
  dir = "CM", filename = "QuinsamCampbell_05.29"
)

if (FALSE) { # Diagnostic figures do not run when sourcing file

  #samp <- samp2

  # Check fits quickly
  report <- salmonMSE::get_report(samp)
  d <- salmonMSE::get_CMdata(samp@.MISC$CMfit)
  CM_fit_esc(report, d)
  CM_fit_pHOS(report, d) # Figure only when fitted to pHOS observations, otherwise nothing

  # Compare brood year pHOS when not fitted
  year <- unique(full_table$RELEASE_YEAR)
  salmonMSE:::.CM_ts(report, year1 = min(year), var = "pHOScensus_brood", ci = TRUE, ylab = "pHOScensus") +
    geom_point(data = pHOS_df, aes(x = `Release_Year`, y = pHOS)) +
    geom_line(data = pHOS_df, aes(x = `Release_Year`, y = pHOS), linetype = 3) +
    labs(x = "Release Year")

  salmonMSE:::.CM_ts(report, year1 = min(year), var = "pHOScensus", ci = TRUE, ylab = "pHOS")

  CM_F(report)
  CM_surv2(report) # Survival to age 2

  CM_maturity(report, d, brood = FALSE)
  CM_M(report)
  CM_SRR(report, year1 = min(year))

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
    rs_names = rs_names, name = "Quinsam/Campbell", year = year,
    dir = "CM", filename = "QuinsamCampbell_05.29"
  )

  SMSY <- salmonMSE:::.CM_SMSY(report, d)
  Srep <- salmonMSE:::.CM_Srep(report, d)
  Sgen <- salmonMSE:::.CM_Sgen(report, d)

  range(Sgen/SMSY, na.rm = TRUE)
  range(SMSY/Srep, na.rm = TRUE)
  #If you remove the dot from the function name CM_SMSY instead of .CM_SMSY you will get the plotting function for the posterior

}

