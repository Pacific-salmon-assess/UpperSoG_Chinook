# Run projections for Adam River
# Libraries
library(tidyverse)
library(salmonMSE)
library(readxl)
library(here)

source("03-Adam-OM.R")

#------------------------------------------------------------------------------
# Run over multiple SOMs

# Create grid of u_preterminal and n_yearling
# Default is u_preterminal is mean(UPT)=0.449 n_yearling in hatch_rel=3561051

hatch_rel <- 0

folder_path <- "SMSE/Adam"
if (!dir.exists(here::here(folder_path))) {
  dir.create(here::here(folder_path), recursive = TRUE)
}

  g <- expand.grid(
    u_preterminal = seq(0, 0.5, 0.02),# Add in FMSY
    n_yearling = hatch_rel * seq(0.5, 1.5, 0.05)
  )

  nOM <- nrow(g)


  # # Run salmonMSE - single instance

  plot.base.case <- FALSE

  if(plot.base.case == TRUE){
    SOM <- readRDS(file.path("SOM", "SOM_Adam_low.rds"))
    # SOM <- readRDS(file.path("SOM", "SOM_Adam_base.rds"))

    # Add brood rule and premove_HOS rules to SOM
    # SOM@Hatchery@f_brood <- f_brood
    # SOM@Hatchery@premove_HOS <- premove_HOS

    SMSE <- salmonMSE(SOM)
    report(SMSE , dir = "SMSE")
    saveRDS(SMSE , file = file.path("SMSE", paste0("Adam_08.19.26low.rds")))

    png(here("figures", "Adam_low_projecions_ts.png"),
        width = 6, height = 8, units = "in", res = 300)
    par(mfrow = c(3, 2))

    plot_statevar_ts(SMSE, var = "NOS", s = 1, figure = TRUE,
                     xlab = "Projection Year", quant = TRUE)
    mtext("(a)",  side = 3, adj = 0, line = 0.5, font = 1, cex = 1)
    dev.off()
  }



  library(snowfall)
  library(parallel)
  sfInit(TRUE, nOM, cpus=(parallel::detectCores()-1))


  tictoc::tic()
  SMSE_list <- sfLapply(1:nrow(g), function(i, g) {
    require(salmonMSE)

    SOM <- readRDS(file.path("SOM", "SOM_Adam_low.rds"))
    # SOM@Hatchery@n_yearling <- g$n_yearling[i]
    SOM@Harvest@u_preterminal <- g$u_preterminal[i]
    SMSE <- salmonMSE(SOM)

    saveRDS(SMSE, file = file.path("SMSE", "Adam", paste0("Adam", i, "low.rds")))

    invisible()

  }, g = g)
  tictoc::toc()

#------------------------------------------------------------------------------
