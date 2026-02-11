

# First, define brood function
# Code assumes zero strays and mark rate = 1!
# f_brood <- function(NO, HO, stray, m = 0, ptarget_NOB = 0, pmax_NOB = 0.5, pmax_esc = 0.7, min_esc = 0) {
#
#   if (!m %in% c(0, 1)) stop("Brood function assumes mark rate of either zero or one.")
#
#   NOB <- array(0, dim(NO))
#   HOB_marked <- HOB_unmarked <- array(0, dim(HO))
#   HOB_stray <- array(0, dim(stray))
#
#   # Rule 1: no brood if fewer than 600 returns
#   if (sum(NO, HO) > min_esc) {
#     # Rule 2: Brood <= 33% of in-river return
#     max_brood <- pmax_esc * sum(NO, HO)
#
#     if (m == 1) {
#       # # Rule 3: total brood cannot exceed twice the natural brood available.
#       # # This means when natural fish are scarce, total brood scales down.
#       # brood_total_cap <- min(max_brood, 2 * sum(NO))
#
#       # Rule 3 Revised: take the min of 50% (pmax_NOB) of sum(NO) and total broodcap
#       # This stops brood from taking all/most of the NO fish
#       # Take as many natural-origin fish as possible (up to the brood cap)
#       brood_total_cap <- max_brood
#       NOB_total <- min(pmax_NOB * sum(NO), brood_total_cap)
#
#       # Fill remaining brood with hatchery fish to reach ptarget_NOB
#       # ptarget = NOB/(NOB + HOB) --> HOB = (NOB - ptarget * NOB)/ptarget
#       HOB_total <- min((NOB_total - ptarget_NOB * NOB_total)/ptarget_NOB, brood_total_cap - NOB_total)
#
#       # Safety: ensure hatchery never exceeds natural due to numerical jitter
#       #pNOB <- NOB_total/(NOB_total + HOB_total)
#       #if (pNOB < ptarget_NOB) {
#       #  excess <- brood_hatchery - brood_natural
#       #  brood_hatchery <- brood_hatchery - excess
#       #}
#
#       if (sum(NO)) {
#         ptake_NOB <- NOB_total/sum(NO)
#         NOB[] <- ptake_NOB * NO
#       }
#
#       if (sum(HO)) {
#         ptake_HOB <- HOB_total/sum(HO)
#         HOB_marked[] <- ptake_HOB * HO
#       }
#     } else if (m == 0) {
#
#       # Rule 3 does not apply (use Rule 2)
#       brood_total_cap <- max_brood
#
#       pHOS <- sum(HO)/sum(NO, HO)
#       pNOB <- 1 - pHOS
#
#       NOB_total <- pNOB * brood_total_cap
#       HOB_total <- (1 - pNOB) * brood_total_cap
#
#       if (sum(NO)) {
#         ptake_NOB <- NOB_total/sum(NO)
#         NOB[] <- ptake_NOB * NO
#       }
#
#       if (sum(HO)) {
#         ptake_HOB <- HOB_total/sum(HO)
#         HOB_unmarked[] <- ptake_HOB * HO
#       }
#     } else {
#       stop("Brood rule only accommodates mark rate of zero or 1")
#     }
#   }
#
#   list(NOB = NOB, HOB_marked = HOB_marked, HOB_unmarked = HOB_unmarked, HOB_stray = HOB_stray)
# }

# Add brood rule to SOM
# SOM@Hatchery@f_brood <- f_brood
#
# # Add in alternative (higher) fecundity
# fec_QC <- c(0, 0, 800, 2000, 2500) # Walters and Korman (2024) removing age6=3000; Filipovic et al. (in revision) RPA.
#
# SOM@Bio@fec <- fec_QC
# SOM@Hatchery@fec_brood <- fec_QC


# Run salmonMSE

SMSE <- salmonMSE(SOM)
# SMSE_highfec <- salmonMSE(SOM)
report(SMSE , dir = "SMSE")
saveRDS(SMSE , file = file.path("SMSE", paste0("QC_11.02.26.rds")))


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
