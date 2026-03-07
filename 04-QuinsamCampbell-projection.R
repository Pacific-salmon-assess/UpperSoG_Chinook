

# First, define brood function
f_brood <- function(NO, HO, stray, m = 0, pmax_esc = 0.7) {

  NOB <- array(0, dim(NO))
  HOB_marked <- HOB_unmarked <- array(0, dim(HO))
  HOB_stray <- array(0, dim(stray))

  # This function will take the maximum amount of brood (70% of escapement) indiscriminately of brood origin
  # However, salmonMSE will return brood that exceeds release target
  # Brood - assume there are no strays in the system

  max_brood <- pmax_esc * sum(NO, HO)

  ptake <- max_brood/sum(NO, HO) # ptake <- pmax_esc

  NOB[] <- ptake * NO
  HOB_unmarked[] <- ptake * HO
  HOB_marked[] <- 0
  HOB_unmarked[] <- ptake * (1 - m) * HO
  HOB_marked[] <- ptake * m * HO

  list(NOB = NOB, HOB_marked = HOB_marked, HOB_unmarked = HOB_unmarked, HOB_stray = HOB_stray)
}

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


# Add brood rule and premove_HOS rules to SOM
SOM@Hatchery@f_brood <- f_brood
SOM@Hatchery@premove_HOS <- premove_HOS

# # Add in alternative (higher) fecundity
# fec_QC <- c(0, 0, 800, 2000, 2500) # Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.
#
# SOM@Bio@fec <- fec_QC
# SOM@Hatchery@fec_brood <- fec_QC


# Run salmonMSE

SMSE <- salmonMSE(SOM)
# SMSE_highfec <- salmonMSE(SOM)
report(SMSE , dir = "SMSE")
saveRDS(SMSE , file = file.path("SMSE", paste0("QC_07.03.26.rds")))


# #------------------------------------------------------------------------------
# # Run over multiple SOMs
#
# # Create grid of u_preterminal and n_yearling
# # Default is u_preterminal is mean(UPT)=0.449 n_yearling in hatch_rel=3561051
#
#
#   g <- expand.grid(
#     u_preterminal = c(0, 0.1, 0.2, 0.3, 0.4, mean(UPT), 0.5),# Add in FMSY
#     n_yearling = c(1E6, 2E6, 3E6, 4E6, 5E6)
#   )
#
#   nOM <- nrow(g)
#
#
#
#   library(snowfall)
#   sfInit(TRUE, nOM)
#
#
#   tictoc::tic()
#   SMSE_list <- sfLapply(1:nrow(g), function(i, g) {
#     require(salmonMSE)
#
#     SOM <- readRDS(file.path("SOM", "SOM_base.rds"))
#     SOM@Hatchery@n_yearling <- g$n_yearling[i]
#     SOM@Harvest@u_preterminal <- g$u_preterminal[i]
#     SMSE <- salmonMSE(SOM)
#
#     saveRDS(SMSE, file = file.path("SMSE", paste0("QC", i, ".rds")))
#
#     invisible()
#
#   }, g = g)
#   tictoc::toc()
#
# #------------------------------------------------------------------------------
