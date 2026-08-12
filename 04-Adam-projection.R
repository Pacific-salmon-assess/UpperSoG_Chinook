# Run projections for Quinsam/Campbell


# # Run salmonMSE - single instance
#
# SMSE <- salmonMSE(SOM)
# report(SMSE , dir = "SMSE")
# saveRDS(SMSE , file = file.path("SMSE", paste0("QC_07.20.26.rds")))

# SMSE <- readRDS(file.path("SMSE", paste0("QC_07.20.26.rds")))
# SMSE <- readRDS(file.path("SMSE", paste0("QC", 1, ".rds")))
# report(SMSE , dir = "SMSE")

#------------------------------------------------------------------------------
# Run over multiple SOMs

# Create grid of u_preterminal and n_yearling
# Default is u_preterminal is mean(UPT)=0.449 n_yearling in hatch_rel=3561051

folder_path <- "SMSE/Adam"
if (!dir.exists(here::here(folder_path))) {
  dir.create(here::here(folder_path), recursive = TRUE)
}

  g <- expand.grid(
    u_preterminal = seq(0, 0.5, 0.02),# Add in FMSY
    n_yearling = hatch_rel * seq(0.5, 1.5, 0.05)
  )

  nOM <- nrow(g)



  library(snowfall)
  library(parallel)
  sfInit(TRUE, nOM, cpus=(parallel::detectCores()-1))


  tictoc::tic()
  SMSE_list <- sfLapply(1:nrow(g), function(i, g) {
    require(salmonMSE)

    SOM <- readRDS(file.path("SOM", "SOM_Adam.rds"))
    # SOM@Hatchery@n_yearling <- g$n_yearling[i]
    SOM@Harvest@u_preterminal <- g$u_preterminal[i]
    SMSE <- salmonMSE(SOM)

    saveRDS(SMSE, file = file.path("SMSE", "Adam", paste0("Adam", i, ".rds")))

    invisible()

  }, g = g)
  tictoc::toc()

#------------------------------------------------------------------------------
