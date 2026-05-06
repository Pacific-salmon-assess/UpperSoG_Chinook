library(tidyverse)
library(reshape2)
library(salmonMSE)


#### Data ----
# Phillips escapement 2009-2025.
# Data prior to 2009 are unreliable (Thea Rachinsky, pers. comm. March 2026)

esc <- readxl::read_excel(
  file.path("data", "Phillips", "fsar-sog-cn-phillips.xlsx"),
  sheet = "Data") %>%
  rename(year = "Analysis Year") %>%
  rename(esc.x = "Max Estimate") %>%
  mutate(nat_spawners =
           ifelse(is.na(`Natural Spawners Adult`), `Natural Spawners Total`, `Natural Spawners Adult`)) %>%
  summarise(
    escapement = sum(esc.x),
    nat_spawners = sum(nat_spawners),
    .by = c(year)
  ) %>%
  select (year, escapement, nat_spawners) %>%
  filter(year>=2012)   %>% #changed from 2009
  mutate(p_spawn = nat_spawners/escapement)

# Fill in NAs with the first value in the time-series
esc$p_spawn[is.na(esc$p_spawn)] <- na.omit(esc$p_spawn)[1]



# Phillips - CWT releases 1989-2020 (note, recoveries up to 2024). CWT smolt releases only up to 2019

rel <-  readxl::read_excel(
    file.path("data", "Phillips", "2026-03-10-PhillipsReleases_AllYears.xlsx"),
    sheet = "Actual Release"
  )

# Phillips - CWT releases of smolt0+ and seapen0+ 1989-2019
# Should I include 'Fed Fry' released in 2020? (now excluded)
# 2020 is the only year Fed Fry are released.No smolt0+/seapen0+ released in 2020
cwt_rel <- rel %>%
  filter(RELEASE_STAGE_NAME %in% c("Smolt 0+", "Seapen 0+", "Fed Fry")) %>%
  # Use only fed fry in Release Year 2020, as these are sized like smolts 0+
  filter(!(RELEASE_STAGE_NAME == "Fed Fry" & RELEASE_YEAR<=2019)) %>%
  summarise(n_CWT = sum(TaggedNum) - sum(ShedTagNum),
            .by = c(RELEASE_YEAR)) %>%
  arrange(RELEASE_YEAR)

lastCWTreleaseYr <- max(cwt_rel %>% filter(n_CWT>0) %>% pull(RELEASE_YEAR))


# Releases aligned by BY
g <- ggplot(cwt_rel, aes(RELEASE_YEAR - 1, n_CWT)) +
  geom_point() +
  geom_line() +
  labs(x = "Brood Year", y = "Phillips River CWT releases")

# Full matrix of ages (1-5) and years (2009 - 2023), which extends to the last
# CWT release year of smolt0+ in 2019 plus 4 years (for 5 year old return)
# Note, there was one 6-year old recovery from 2019 releases in 2024, but this
# is not included
full_matrix <- expand.grid(
  RELEASE_YEAR = min(esc$year):(lastCWTreleaseYr + 4), #max(rel$RELEASE_YEAR),
  Age = 1:5
) %>%
  as.data.frame()

full_year <- data.frame(RELEASE_YEAR = min(esc$year):(lastCWTreleaseYr + 4))#max(rel$RELEASE_YEAR))

# Get tag codes for CWT releases, including smolt0+, seapen0+, excl. Fed Fry

cwt_rel_tags <- rel %>%
  filter(RELEASE_STAGE_NAME %in% c("Smolt 0+", "Seapen 0+", "Fed Fry")) %>%
  # Use only fed fry in Release Year 2020, as these are sized like smolts 0+
  filter(!(RELEASE_STAGE_NAME == "Fed Fry" & RELEASE_YEAR<=2019)) %>%
  filter(RELEASE_YEAR >= min(full_year) & RELEASE_YEAR <= max(full_year)) %>%
  summarise(n_CWT = sum(TaggedNum) - sum(ShedTagNum), .by = c(RELEASE_YEAR, MRP_TAGCODE)) %>%
  filter(n_CWT > 0) %>%
  arrange(RELEASE_YEAR) %>%
  rename(tag_code = "MRP_TAGCODE") %>%
  select(tag_code, RELEASE_YEAR)

cwt_rel_tags$tag_code <- as.numeric(cwt_rel_tags$tag_code)

# CWT recoveries
cwt_dat <- readr::read_csv("data/Phillips/PHI_camp_recovery_wfisheries.csv") %>%
  rename(Age = "age")
#   mutate(RELEASE_YEAR = run_year - age + 1)

# problems(cwt_dat)
# cwt_dat[1220,13]

# Only recoveries from release type smolt 0+ and seapen 0+ included by matching
# tag code in the recoveries with tagcodes from releases for smolts0+ and
# seapen 0+ (~300 recoveries removed)
# In release year 2020, only fed fry are released and are excluded for now

# Used RELEASE_YEAR from release file (cwt_rel_tags) instead of
# RELEASE_YEAR = runyear - age + 1 in recoveries (cwt_dat), as there are some
# mistakes in aging (e.g., recovery_id 624505 which I think should be age 4, not 5)

cwt_dat_subset <- inner_join(cwt_rel_tags, cwt_dat, by=c("tag_code"))

cwt_esc <- cwt_dat_subset %>%
  filter(fishery_type == "escapement", Coarse_description %in% c("Escapement", "Subsistence")) %>%
  summarise(n = sum(adjusted_estimated_number), .by = c(RELEASE_YEAR, Age)) %>%
  right_join(full_matrix, by = c("RELEASE_YEAR", "Age")) %>%
  reshape2::acast(list("RELEASE_YEAR", "Age"), value.var = "n", fill = 0)

# Preterminal CWT
cwt_pt <- cwt_dat_subset %>%
  filter(fishery_type == "pre-terminal", Age < 7) %>%
  summarise(n = sum(adjusted_estimated_number), .by = c(RELEASE_YEAR, Age)) %>%
  right_join(full_matrix, by = c("RELEASE_YEAR", "Age")) %>%
  reshape2::acast(list("RELEASE_YEAR", "Age"), value.var = "n", fill = 0)

# cwt_t <- cwt_dat_subset %>%
#   filter(fishery_type == "terminal") %>%
#   summarise(n = sum(adjusted_estimated_number), .by = c(RELEASE_YEAR, Age)) %>%
#   right_join(full_matrix, by = c("RELEASE_YEAR", "Age")) %>%
#   reshape2::acast(list("RELEASE_YEAR", "Age"), value.var = "n", fill = 0)
#


#### Process data ----

esc <- esc %>%
  right_join(
    full_matrix %>% filter(Age == 1) %>% select(RELEASE_YEAR),
    by = c("year" = "RELEASE_YEAR")
  ) %>%
  arrange(year)

cwt_rel <- cwt_rel %>%
  right_join(
    full_matrix %>% filter(Age == 1) %>% select(RELEASE_YEAR),
    by = c("RELEASE_YEAR")
  ) %>%
  arrange(RELEASE_YEAR)

cwt_rel$n_CWT[which(is.na(cwt_rel$n_CWT))] <- 0

#All releases

# Total hatchery releases from Phillips: all release types
# File includes releases from:
# Big Bay, Fanny Bay/JNST, Phillips Lk, Phillips R, Quatam Est
rel_total <- rel %>%
  summarise(n_rel = sum(TotalRelease), .by = c(RELEASE_YEAR)) %>%
  arrange(RELEASE_YEAR)  %>%
  right_join(
    full_matrix %>% filter(Age == 1) %>% select(RELEASE_YEAR),
    by = c("RELEASE_YEAR")
  ) %>%
  arrange(RELEASE_YEAR)

rel_total$n_rel[which(is.na(rel_total$n_rel))] <- 0

# Plot CWT data
if (FALSE) {
  cwt_plot <- rbind(
    cwt_esc %>% reshape2::melt() %>% mutate(type = "Escapement"),
    cwt_pt %>% reshape2::melt() %>% mutate(type = "Preterminal catch")#,
    # cwt_t %>% reshape2::melt() %>% mutate(type = "Terminal catch")
  ) %>%
    rename(RelYear = Var1, Age = Var2, n = value) %>%
    mutate(p = n/sum(n, na.rm = TRUE), .by = c(RelYear, type))

  g <- cwt_plot %>%
    filter(!is.na(p)) %>%
    ggplot(aes(RelYear - 1, p, fill = factor(Age, levels = 5:2))) +
    facet_wrap(vars(type), ncol = 1) +
    geom_col(width = 1, colour = "grey40", linewidth = 0.1) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Brood Year", y = "Proportion", fill = "Age", title = "Phillips R CWT") +
    coord_cartesian(expand = FALSE)
  ggsave(paste0("figures/Phillips_CWT_proportion.png"), g, height = 6, width = 6)

  g2 <- g +
    coord_cartesian(expand = FALSE, xlim = c(2008.5, 2020.5))
  ggsave("figures/Phillips_CWT_prop2.png", g2, height = 4, width = 3)

  g <- cwt_plot %>%
    filter(!is.na(p)) %>%
    ggplot(aes(RelYear - 1, n, fill = factor(Age, levels = 5:2))) +
    facet_wrap(vars(type), ncol = 1) +
    geom_col(width = 1, colour = "grey40", linewidth = 0.1) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Brood Year", y = "Catch", fill = "Age", title = "Phillips R CWT") +
    coord_cartesian(expand = FALSE)
  ggsave(paste0("figures/Phillips_CWT_catch.png"), g, height = 6, width = 6)

  g2 <- g +
    coord_cartesian(expand = FALSE, xlim = c(2008.5, 2020.5))
  ggsave("figures/Phillips_CWT_esc_prop2.png", g2, height = 4, width = 6)
}

# Data object for model
Ldyr <- nrow(cwt_esc)
Nages <- 5

mat <- c(0, 0.01, 0.05, 0.2, 1) # Need to tune this vector for initial abundance #from WCVI = c(0, 0.1, 0.4, 0.95, 1)
vulPT <- c(0, 0.075, 0.9, 0.9, 1)
vulT <- vulPT

M_CTC <- -log(1 - c(0.9, 0.3, 0.2, 0.1, 0.1)) # CTC 23-06 p.9; CWT Exploitation Rate analyses
M_CTC[1] <- 4 # Need to tune this value for initial abundance

fec_Quinsam <- c(0, 0, 800, 2000, 2500) # Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.


# Model assumption of catch expansion factor
# Use alternative values to change data weighting of CWT (re-adjust numbers accordingly)
cwtExp <- 1

d <- list(
  Nages = Nages,
  Ldyr = Ldyr,
  lht = 1,
  n_r = 1,
  s_enroute = 1, # For Upper SoG 0.855 = 10% return migration M followed by 5% terminal ER
  cwtrelease = cwt_rel$n_CWT,
  cwtesc = array(round(cwt_esc/cwtExp), c(Ldyr, Nages, 1)), #Aligned by RY
  cwtcatPT = array(round(cwt_pt/cwtExp), c(Ldyr, Nages, 1)), #Aligned by RY
  cwtcatT = NULL, #array(round(cwt_t/cwtExp), c(Ldyr, Nages, 1)), #Aligned by RY
  bvulPT = vulPT,
  # bvulT = vulT,
  RelRegFPT = rep(1, Ldyr),
  RelRegFT = rep(1, Ldyr),
  bmatt = mat,
  mobase = M_CTC,
  hatchsurv = 0.8,#From M. Clarke life-cycle table. Walters and Korman (2024) used 0.5; 1 used for WCVI Chinook
  gamma = 0.8,
  ssum = 1,
  fec = fec_Quinsam*0.95,
  obsescape = esc$escapement,
  propwildspawn = round(esc$p_spawn, 2),
  hatchrelease = c(rel_total$n_rel, 0),
  s_enroute = 1,
  finitPT = 0.8,
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

start <- list(log_so = log(2 * max(d$obsescape)))

# Fit with sampling rate = 1
fit <- fit_CM(d, start = start, map = map, do_fit = TRUE)
samp <- sample_CM(fit, chains = 4, cores = 4, iter = 10000, thin = 5,
                  control=list(adapt_delta = 0.999, stepsize = 0.01,
                               max_treedepth = 20))
saveRDS(samp, file = paste0("CM/Phillips_CM_05.06.26.rds"))

samp <- readRDS(file = "CM/Phillips_CM_05.06.26.rds")
report <- salmonMSE:::get_report(samp)
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
#shinystan::launch_shinystan(samp)

rs_names <- c("Smolt 0+")
salmonMSE::report_CM(
  samp,
  rs_names = rs_names, name = "Phillips", year = unique(full_matrix$RELEASE_YEAR),
  dir = "CM", filename = "Phillips_05.06"
)


if (FALSE) { # Diagnostic figures do not run when sourcing file

  # Check fits quickly
  report <- salmonMSE::get_report(samp)
  d <- salmonMSE::get_CMdata(samp@.MISC$CMfit)
  CM_fit_esc(report, d)

  # Compare brood year pHOS when not fitted
  year <- unique(full_year$RelYear)
  salmonMSE:::.CM_ts(report, year1 = min(year), var = "pHOScensus_brood", ci = TRUE, ylab = "pHOScensus")
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

}
