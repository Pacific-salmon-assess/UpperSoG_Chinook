# Run projections for Quinsam/Campbell

#libraries
library(tidyverse)
library(readxl)
library(salmonMSE)

# First, define brood function
f_brood <- function(NO, HO, stray, m = 0, pmax_esc = 0.7) {

  NOB <- array(0, dim(NO))
  HOB_marked <- HOB_unmarked <- array(0, dim(HO))
  HOB_stray <- array(0, dim(stray))

  # This function will take the maximum amount of brood (pmax_esc = 0.7, 0.8, etc of escapement)
  # of age 4 and 5, not selective for brood origin
  # However, salmonMSE will return brood that exceeds release target
  # Assume there are no strays in the system
  ptake <- c(0, 0, 0, rep(pmax_esc, 2))

  NOB[] <- ptake * NO
  HOB_unmarked[] <- ptake * (1 - m) * HO
  HOB_marked[] <- ptake * m * HO

  list(NOB = NOB, HOB_marked = HOB_marked, HOB_unmarked = HOB_unmarked, HOB_stray = HOB_stray)
}

# Define rule for non-brood removals. Not used
premove_HOS <- function(NO, HO, m = 1, p.x=0.102) {
  # Assumes marking, such that only marked fish are removed for CWT sampling
  # p.x is the ppn of total escapement that are taken for CWT sampling, estimated from Quinsam time-series.
  # However, only marked fish are sampled, so ppn of HO escapement that is taken (p) is larger than p.x
  p <- 0
  if (sum(HO)) {
    pHOS <- sum(HO)/sum(NO, HO)
    p <- p.x / (pHOS)
  }
  # return(p)
  return(p)
}



# Get hatchery releases
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
           str_starts(RELEASE_SITE_NAME, "Orange") |
           str_starts(RELEASE_SITE_NAME, "Deepwater") |
           str_starts(RELEASE_SITE_NAME, "Drew") |
           str_starts(RELEASE_SITE_NAME, "Taku") ) %>%
  filter(RELEASE_STAGE_NAME %in% c("Smolt 0+", "Seapen 0+")) %>%
  summarise(n_rel = sum(TotalRelease), .by = c(RELEASE_YEAR)) %>%
  arrange(RELEASE_YEAR)

# Average hatchery releases from Q/C over last 5 years
hatch_rel <- rel_Quinsam %>%
  filter(RELEASE_YEAR %in% c(max(RELEASE_YEAR) - seq(4, 0))) %>%
  pull(n_rel) %>%
  mean()


# # Run salmonMSE - single instance

plot.base.case <- FALSE

if(plot.base.case == TRUE){
  SOM <- readRDS(file.path("SOM", "SOM_QC_low.rds"))

    # Add brood rule and premove_HOS rules to SOM
  SOM@Hatchery@f_brood <- f_brood
  # SOM@Hatchery@premove_HOS <- premove_HOS

  SMSE <- salmonMSE(SOM)
  report(SMSE , dir = "SMSE")
  saveRDS(SMSE , file = file.path("SMSE", paste0("QC_08.18.26low.rds")))

  png(here("figures", "QC_low_projecions_ts.png"),
      width = 6, height = 8, units = "in", res = 300)
  par(mfrow = c(3, 2))

  plot_statevar_ts(SMSE, var = "NOS", s = 1, figure = TRUE,
                   xlab = "Projection Year", quant = TRUE)
  mtext("(a)",  side = 3, adj = 0, line = 0.5, font = 1, cex = 1)
  plot_statevar_ts(SMSE, var = "HOS", s = 1, figure = TRUE,
                   xlab = "Projection Year", quant = TRUE)
  mtext("(b)",  side = 3, adj = 0,line = 0.5, font = 1, cex = 1)
  plot_statevar_ts(SMSE, var = "pHOS_census", s = 1, figure = TRUE,
                   xlab = "Projection Year", quant = TRUE, ylab="pHOS")
  mtext("(c)",  side = 3, adj = 0,line = 0.5, font = 1, cex = 1)
  plot_statevar_ts(SMSE, var = "NOB", s = 1, figure = TRUE,
                   xlab = "Projection Year", quant = TRUE, ylab="NOB")
  mtext("(d)",  side = 3, adj = 0,line = 0.5, font = 1, cex = 1)
  plot_statevar_ts(SMSE, var = "pNOB", s = 1, figure = TRUE,
                   xlab = "Projection Year", quant = TRUE, ylab="pNOB")
  mtext("(e)",  side = 3, adj = 0,line = 0.5, font = 1, cex = 1)
  plot_statevar_ts(SMSE, var = "PNI", s = 1, figure = TRUE,
                   xlab = "Projection Year", quant = TRUE, ylab="PNI")
  mtext("(f)",  side = 3, adj = 0,line = 0.5, font = 1, cex = 1)
  dev.off()
}


# SMSE <- readRDS(file.path("SMSE", paste0("QC_07.20.26.rds")))
# SMSE <- readRDS(file.path("SMSE", paste0("QC", 1, ".rds")))
# report(SMSE , dir = "SMSE")

#------------------------------------------------------------------------------
# Run over multiple SOMs

# Create grid of u_preterminal and n_yearling
# Default is u_preterminal is mean(UPT)=0.449 n_yearling in hatch_rel=3561051


  g <- expand.grid(
    u_preterminal = seq(0, 0.5, 0.02),# Add in FMSY
    n_yearling = hatch_rel * seq(0.5, 1.5, 0.05)
  )

  nOM <- nrow(g)

  folder_path <- "SMSE/QC"
  if (!dir.exists(here::here(folder_path))) {
    dir.create(here::here(folder_path), recursive = TRUE)
  }


  library(snowfall)
  library(parallel)
  sfInit(TRUE, nOM, cpus=(parallel::detectCores()-1))
  # Define brood function for parallel process
  f_brood <- function(NO, HO, stray, m = 0, pmax_esc = 0.7) {

    NOB <- array(0, dim(NO))
    HOB_marked <- HOB_unmarked <- array(0, dim(HO))
    HOB_stray <- array(0, dim(stray))

    # This function will take the maximum amount of brood (pmax_esc = 0.7, 0.8, etc of escapement)
    # of age 4 and 5, not selective for brood origin
    # However, salmonMSE will return brood that exceeds release target
    # Assume there are no strays in the system
    ptake <- c(0, 0, 0, rep(pmax_esc, 2))

    NOB[] <- ptake * NO
    HOB_unmarked[] <- ptake * (1 - m) * HO
    HOB_marked[] <- ptake * m * HO

    list(NOB = NOB, HOB_marked = HOB_marked, HOB_unmarked = HOB_unmarked, HOB_stray = HOB_stray)
  }

  sfExport(list = "f_brood")
  tictoc::tic()
  SMSE_list <- sfLapply(1:nrow(g), function(i, g) {
    require(salmonMSE)

    SOM <- readRDS(file.path("SOM", "SOM_QC_high.rds"))
    SOM@Hatchery@f_brood <- f_brood
    SOM@Hatchery@n_yearling <- g$n_yearling[i]
    SOM@Harvest@u_preterminal <- g$u_preterminal[i]
    SMSE <- salmonMSE(SOM)

    saveRDS(SMSE, file = file.path("SMSE", "QC", paste0("QC", i, "high.rds")))
    # saveRDS(SMSE, file = file.path("SMSE", "QC", paste0("QC_basecase.rds")))

    invisible()

  }, g = g)
  tictoc::toc()

#------------------------------------------------------------------------------
