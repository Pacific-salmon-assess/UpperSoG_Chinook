
library(tidyverse)
library(readxl)
library(salmonMSE)



#### Data ----
# Escapement time-series
pop <- "Adam" #Campbell, Adam, Nimpkish, Salmon


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

esc_2025 <- data.frame(year= 2025, escapement = 256, nat_spawners=256) # From Andrew Pereboom 9 April 2026
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
  RELEASE_YEAR = (max(min(cwt_dat_subset$RELEASE_YEAR), min(esc_all$year)) ): 2024,
  Age = seq(1, 5)#6)
  # RS = c( "Seapen/Smolt 0+")
)
full_year <- data.frame(RELEASE_YEAR =
                          (max(min(cwt_dat_subset$RELEASE_YEAR), min(esc_all$year)) ):2024)


# Escapement CWT
# Two terminal fisheries are considered as 'escapement' for this population,
# as these fish are not vulnerable to these fisheries, so added here:
# "TGS FS" North Georgia Strait Freshwater Sport (6 tags)
# "TGEO ST TERM S" North Georgia Strait Terminal Sport (562 tags)
# "TJNST TERM S" Johnstone Strait Terminal Sport (246 tags)

cwt_esc <- cwt_dat_subset %>%
  filter( (fishery_type == "escapement" &
         Coarse_description %in% c("Escapement", "Subsistence")) |
           fishery_era_name %in% c("TGS FS",
                                   "TGEO ST TERM S")) %>%
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
                                      "TNORTH FS",
                                      "TJNST TERM S"))) %>%
  summarise(n = sum(adjusted_estimated_number), .by = c(RELEASE_YEAR, Age)) %>%
  right_join(full_matrix, by = c("RELEASE_YEAR", "Age")) %>%
  reshape2::acast(list("RELEASE_YEAR", "Age"), value.var = "n", fill = 0)


# Terminal CWTs-
# Three terminal fisheries are considered pre-terminal here, so removed:
# "TAK TERM T" Alaska Terminal Troll (15 tags)
# "TWCVI TERM N" Southwest WCVI Terminal Net (1 tag)
# "TNORTH FS" North Freshwater Sport (1 tag)
# Two other pre-terminal fisheries are considered escapement, so removed:
# "TGS FS" North Georgia Strait Freshwater Sport (6 tags)
# "TGEO ST TERM S" North Georgia Strait Terminal Sport (562 tags)

# cwt_t <- cwt_dat_subset %>%
#   filter( fishery_type == "terminal" &
#            !fishery_era_name %in% c("TAK TERM T",
#                                     "TWCVI TERM N",
#                                     "TNORTH FS",
#                                     "TGS FS",
#                                     "TGEO ST TERM S")) %>%
#   summarise(n = sum(adjusted_estimated_number), .by = c(RELEASE_YEAR, Age)) %>%
#   right_join(full_matrix, by = c("RELEASE_YEAR", "Age")) %>%
#   reshape2::acast(list("RELEASE_YEAR", "Age"), value.var = "n", fill = 0)

# Note terminal fishery is:
# "TJNST TERM S" Johnstone Strait Terminal Sport (246 tags)

esc <- esc_all %>%
  right_join(
    full_matrix %>% filter(Age == 1) %>% select(RELEASE_YEAR),
    by = c("year" = "RELEASE_YEAR")
  ) %>%
  arrange(year) %>%
  mutate(p_spawn = nat_spawners/escapement)
esc$p_spawn[is.na(esc$p_spawn)] <- na.omit(esc$p_spawn)[1]

cwt_rel <- left_join(full_year, cwt_rel,by = "RELEASE_YEAR")

# Data object for model
Ldyr <- nrow(cwt_esc)
Nages <- 5

### Get initial maturity, by tuning in the model
# ERM_tuned <-  readRDS("CM/Adam_08.06.26.prior.rds")
# ERM_tuned <-  readRDS("CM/Adam_09.09.26.rds")
# report_tuned <- salmonMSE:::get_report(ERM_tuned)
# matt <- sapply(report_tuned, function(i)
#   salmonMSE:::CY2BY(i[["matt"]][, , 1]), simplify = 'array') %>%
#   apply(2, quantile,
#         probs =  0.5, na.rm = TRUE) #median over years and MC trials
# matt = c(0, 0.00599, 0.123, 0.690, 1) # First tuning
# matt = c(0, 0.00552, 0.121, 0.676, 1) # Second tuning
mat <- c(0, 0.00552, 0.121, 0.676, 1)# Tune this vector for initial abundance #from WCVI = c(0, 0.1, 0.4, 0.95, 1)

vulPT <- c(0, 0.075, 0.9, 0.9, 1)
vulT <- vulPT#rep(0, Nages)#

M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1)) # CTC 23-06 p.9; CWT Exploitation Rate analyses
### Get initial M in year 1, by tuning in the model
# ERM_tuned <-  readRDS("CM/Adam_08.06.26.prior.rds")
# ERM_tuned <-  readRDS("CM/Adam_09.09.26.rds")
# report_tuned <- salmonMSE:::get_report(ERM_tuned)
# mo <- sapply(report_tuned, function(x) x$mo[, 1]) %>%
# quantile(probs = 0.5)
# mo = 5.560  # First tuning
# mo = 5.552  # Second tuning
M_CTC[1] <- 5.552 #4 # Tune this value for initial abundance

fec_Quinsam <- c(0, 0, 800, 2000, 2500) # Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.


# Model assumption of catch expansion factor
# Use alternative values to change data weighting of CWT (re-adjust numbers accordingly)
cwtExp <- 1

# # Srep prior
# data_Srep_prior <- as.data.frame( read.csv(
#   ("data/UpperSoGChinook_out_posteriorpredictive_NEWWArev.csv"))
# )
# Srep_prior <- data_Srep_prior %>% filter(Stock=="Adam/Eve") %>% pull(SREP_median)
# logSrep_prior_sd <- data_Srep_prior %>% filter(Stock=="Adam/Eve") %>%
#   mutate(sigma=(log(SREP_upr95)-log(SREP_median))/2) %>%
#   pull(sigma)

data_Smax_prior <- as.data.frame( read.csv(
  ("data/UpperSoGChinook_out_posteriorpredictive_NEWWArev.csv"))
)
med_Smax_prior <- data_Smax_prior %>% filter(Stock == "Adam/Eve") %>% pull(SMAX_median)
logSmax_prior_sd <- data_Smax_prior %>% filter(Stock == "Adam/Eve") %>%
  mutate(sigma=(log(SMAX_upr95)-log(SMAX_median))/2) %>%
  pull(sigma)

# Productivity for Cowichan Chinook (= mean prod for Fraser river, 3 stocks)
# (Greenberg et al. in prep, Table S4)
mean_logalpha <- 0.87
logalpha_sig <- 0.23

logalpha <- rnorm(2000, mean_logalpha, logalpha_sig)
logSmax <- rnorm(2000, log(med_Smax_prior), logSmax_prior_sd)

Srep_prior <-  logalpha * exp(logSmax) # Srep = log(alpha)/beta, beta = 1/Smax
logSrep_prior_sd <- sd(log(Srep_prior))

d <- list(
  Nages = Nages,
  Ldyr = Ldyr,
  lht = 1,
  n_r = 1,
  s_enroute = 0.9,
  cwtrelease = as.vector(cwt_rel$n_CWT),
  cwtesc = array(round(cwt_esc/cwtExp), c(Ldyr, Nages, 1)),
  cwtcatPT = array(round(cwt_pt/cwtExp), c(Ldyr, Nages, 1)),
  cwtcatT = NULL, #array(round(cwt_t/cwtExp), c(Ldyr, Nages, 1)), #NULL,
  bvulPT = vulPT,
  bvulT = vulT,
  RelRegFPT = rep(1, Ldyr),
  RelRegFT = rep(1, Ldyr),
  bmatt = mat,
  mobase = M_CTC,
  #hatchsurv = 0.8,#From M. Clarke life-cycle table. Walters and Korman (2024) used 0.5; 1 used for WCVI Chinook
  gamma = 0.8,
  ssum = 1,
  fec = fec_Quinsam*0.95,
  obsescape = esc$escapement,
  propwildspawn = round(esc$p_spawn, 2),
  hatchrelease =  rep(0, Ldyr + 1),
  finitPT = 0.4,
  finitT = 0,#0.1,#, #0,#,0.8,
  cwtExp = cwtExp,
  so_mu =  NULL, #mean(log(Srep_prior)),#log(3 * max(esc$escapement, na.rm = TRUE)), #prior on S0, reduce from default 3x to 1.5x
  so_sd = NULL, #round(logSrep_prior_sd, 2)# #SD of prior on S0, reduce from default 0.5 to 0.2. Change to uncertainty in logSmax from IWAM
  smax_mu = log(med_Smax_prior), # To be updated by Tor
  smax_sd = logSmax_prior_sd # To be updated by Tor

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

start <- list(log_so = log(2 * max(d$obsescape, na.rm = TRUE)))

#### Fit with estimated productivity parameter (log_cr)
fit <- fit_CM(d, start = start, map = map, do_fit = TRUE, lower = list(moadd = -Inf))
samp <- sample_CM(fit, chains = 4, cores = 4, iter = 10000, thin = 5, seed = 1,
                  control=list(adapt_delta = 0.999,
                               stepsize = 0.01,
                               max_treedepth = 20))
saveRDS(samp, file = "CM/Adam_09.19.26.rds")

samp <- readRDS(file = "CM/Adam_09.19.26.rds")

year <- unique(full_matrix$RELEASE_YEAR)
rs_names <- c("Smolt 0+")
salmonMSE::report_CM(
  samp,
  rs_names = rs_names, name = "Adam", year = year,
  dir = "CM", filename = "Adam_09.19"
)


