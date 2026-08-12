

# First, create grid of pNOB target and ESSR exploitation rate target
# Default is Catch target of 1000 and pNOB = 0.5
g <- expand.grid(
  ESSR_ER = c(0.5, 0.75, 1),
  pNOB_target = c(0.5, 0.75, 1)
)

nOM <- nrow(g)

# Additional runs
gadd <- data.frame(
  name = c("Traditionals", "Fed Fry",
           "HighSurv", "HighSurvNoHarvestNoHatchery", "NoHarvestNoHatchery",
           "DensityDependence", "DeclineMaturity")
) %>%
  mutate(i = 1:nrow(.), OM = nOM + i)

library(snowfall)
sfInit(TRUE, nOM)
sfExport(f_brood)


tictoc::tic()
SMSE_list <- sfLapply(1:nrow(g), function(i, g) {
  require(salmonMSE)

  SOM <- readRDS(file.path("SOM", "SOM_base.rds"))
  SOM@Hatchery@ptarget_NOB <- g$pNOB_target[i]
  SOM@Hatchery@phatchery <- g$ESSR_ER[i]
  SMSE <- salmonMSE(SOM)

  saveRDS(SMSE, file = file.path("SMSE", paste0("Sarita", i, ".rds")))

  invisible()

}, g = g)
tictoc::toc()

tictoc::tic()
SMSE_list <- sfLapply(1:nrow(gadd), function(i, g) {
  require(salmonMSE)

  if (grepl("HighSurv", g$name[i])) {
    SOM <- readRDS(file.path("SOM", "SOM_highsurv.rds"))
  } else {
    SOM <- readRDS(file.path("SOM", "SOM_base.rds"))
  }

  if (g$name[i] == "Traditionals") {
    SOM@Hatchery@n_yearling <- c(0.1, 0.9) * 500000
  } else if (g$name[i] == "Fed Fry") {
    SOM@Hatchery@n_yearling <- c(0.9, 0.1) * 500000
  } else if (grepl("NoHarvestNoHatchery", g$name[i])) {
    SOM@Hatchery@n_yearling[] <- 0
    SOM@Hatchery@stray_external[] <- 0

    SOM@Harvest@u_preterminal <-
      SOM@Harvest@u_terminal <- 1e-8

  } else if (g$name[i] == "DensityDependence") {

    SOM@Habitat@prespawn_rel <- "HS"
    SOM@Habitat@prespawn_prod <- 1
    SOM@Habitat@prespawn_capacity <- 1300

  } else if (g$name[i] == "DeclineMaturity") {
    SOM <- readRDS(file.path("SOM", "SOM_declinemat.rds"))
  }

  SMSE <- salmonMSE(SOM)
  saveRDS(SMSE, file = file.path("SMSE", paste0("Sarita", g$OM[i], ".rds")))

  invisible()
}, g = gadd)


sfStop()
tictoc::toc()



#SOM <- readRDS("SOM/SOM_base_3sim.rds")
#SMSE <- salmonMSE(SOM)
#saveRDS(SMSE, file = "SMSE/Sarita_base2.rds")

#lhg_name <- c("Early Smalls", "Late Larges")
#rs_name <- c("Fed Fry", "Smolt 0+")
#report(SMSE, dir = "SMSE", filename = "Sarita_base2")
