# Run projections for Quinsam/Campbell


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

folder_path <- "SMSE/Nimpkish"
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
  sfExport("f_brood")


  tictoc::tic()
  SMSE_list <- sfLapply(1:nrow(g), function(i, g) {
    require(salmonMSE)

    SOM <- readRDS(file.path("SOM", "SOM_Nimpkish_base.rds"))
    SOM@Hatchery@f_brood <- f_brood
    SOM@Hatchery@n_yearling <- g$n_yearling[i]
    SOM@Harvest@u_preterminal <- g$u_preterminal[i]
    SMSE <- salmonMSE(SOM)

    saveRDS(SMSE, file = file.path("SMSE", "Nimpkish", paste0("Nimpkish", i, ".rds")))

    invisible()

  }, g = g)
  tictoc::toc()

#------------------------------------------------------------------------------
