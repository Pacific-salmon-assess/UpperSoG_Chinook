#-------------------------------------------------------------------------------
# Code to plot projection outputs: contour (ice-cream) plots and trade-off plots
# For Salmon River
#-------------------------------------------------------------------------------
# Libraries
library(salmonMSE)
library(tidyverse)
library(reshape2)
library(here)
library(patchwork)
library(ggplot2)

pop <- "Salmon"

# Get CM output for current PT and T ERs
samp <-  readRDS("CM/Salmon_10.22.26.base.rds")
# year1 <- 2002 # For Salmon
# scenario <- "base"
# if (scenario !="base") {file_suffix <- paste0("_", scenario)} else {file_suffix <- ""}

folder_path <- "figures/SMSE/Salmon/NoHatch"

if (!dir.exists(here::here(folder_path))) {
  dir.create(here::here(folder_path), recursive = TRUE)
}

# # Get UMSY from the Conditioning model = 28%
# report <- salmonMSE:::get_report(samp)
# d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
#
# UMSY_s <- salmonMSE:::.CM_MSY(report, d, simple = TRUE, type = "u")
# UMSY_avg <- salmonMSE:::.CM_MSY(report, d, simple = TRUE,
#                                 mean_bio = TRUE, type = "u")
# quantile(UMSY_avg, c(0.05, 0.25, 0.5, 0.75, 0.95))

# Now plot catches and spawner abundances across a range of pre-terminal ERs (scenarios), with no hatchery relesases

# Identify scenarios and management options ----
gr <- expand.grid(
  u_preterminal = seq(0, 0.5, 0.02),# Add in FMSY
  n_yearling = 0
)
gr$n <- 1:length(gr$u_preterminal)


# Grab all model output ----
# scenario_unique <- unique(gr$Option_name) # Represented by individual table

SMSE_list <- lapply(gr$n, function(i) {
  SMSE <- readRDS(file.path("SMSE",pop, "NoHatch", paste0(pop, i, ".rds")))
  return(SMSE)
})

# Performance metrics (all simulations at end of projection) ----
y <- 29



NS <- sapply(SMSE_list, function(x) {
  NS_a <- x@NOS[, 1, , y] + x@HOS[, 1, , y]
  rowSums(NS_a)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, `Natural Spawners` = value, n = Var2)


Aggcatch <- sapply(SMSE_list, function(x) {#Aggregate catch
  tc_HOS <- apply(x@KT_HOS, c(1, 2, 4), sum)
  tc_NOS <- apply(x@KT_NOS, c(1, 2, 4), sum)
  ptc_HOS <- apply(x@KPT_HOS, c(1, 2, 4), sum)
  ptc_NOS <- apply(x@KPT_NOS, c(1, 2, 4), sum)
  tc_NOS[, ,y]  + ptc_NOS[, ,y] + tc_HOS[, ,y]+ ptc_HOS[, ,y]
  }) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Aggcatch = value, n = Var2)

# State variables for all simulations in year y in one data frame
val_sim_all <- list(NS, Aggcatch) %>%
Reduce(left_join, .) %>%
  left_join(gr %>% select( n), by = "n") %>%
  reshape2::melt(id.vars = c("n", "Simulation"))


# Median and range for state variables across simulations
val_sim <- val_sim_all %>%
  summarise(m = mean(value),
            median = median(value, na.rm = TRUE),
            lwr = quantile(value, 0.25, na.rm = TRUE),
            upr = quantile(value, 0.75, na.rm = TRUE),
            .by = c( "n", "variable"))

#### Trade-off figures


# Spawners vs total catch
Sp_Catch <- val_sim %>%
  filter(variable %in% c("Natural Spawners", "Aggcatch")) %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  reshape2::dcast(u_preterminal + n_yearling ~ variable, value.var = "median")


# Get abundances and catch of Ricker UMSY
df <- val_sim %>%
  filter(variable %in% c("Natural Spawners", "Aggcatch")) %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  reshape2::dcast(u_preterminal + n_yearling ~ variable, value.var = "median")
current <- df %>% filter(u_preterminal == round(UPTmed,2) ) %>%
  rename(y = 'Natural Spawners', x = 'Aggcatch') %>% select(x, y)

log_cr <- rstan::extract(samp)[["log_cr"]]
umsy <- calc_Umsy_Ricker(log_cr) %>% median() - 0.0275
# Round umsy to an even number
umsy <- round(umsy * 100 /2) * 2 /100
rickerMSY <- df %>% filter(u_preterminal == umsy) %>%
  rename(y = 'Natural Spawners', x = 'Aggcatch') %>% select(x, y)

# Get abundances and catch of salmonMSE UMSY
# source('99-RefPts-UpperSoG.R')
df <- val_sim %>%
  filter(variable %in% c("Natural Spawners", "Aggcatch")) %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  reshape2::dcast(u_preterminal + n_yearling ~ variable, value.var = "median")
umsy_sMSE <- round(SalmonRP$salmonMSE[['UMSY (preterminal)']] * 100 /2) * 2 /100
sMSEMSY <- df %>% filter(u_preterminal == umsy_sMSE) %>%
  rename(y = 'Natural Spawners', x = 'Aggcatch') %>% select(x, y)

gto3 <- plot_tradeoff(
  pm1 = Sp_Catch$`Aggcatch`,
  pm2 = Sp_Catch$'Natural Spawners',
  x1 = factor(Sp_Catch$u_preterminal),
   x2 = Sp_Catch$n_yearling,
  xlab = "Aggregate catch (median)",
  ylab = "Natural spawners (median)",
  x1lab = "Pre-terminal\nexploitation\nrate",
  x2lab = "Relative\nhatchery\nproduction"
) +
  geom_point(size = 3) +
  scale_shape_manual(values = 16) +
  # annotate("point", x = current$x, y = current$y, shape = 18, size = 5)
  annotate("point", x = rickerMSY$x, y = rickerMSY$y, shape = 8, size = 3) +
  annotate("point", x = sMSEMSY$x, y = sMSEMSY$y, shape = 6, size = 3)

gto3

ggsave(file.path("figures", "SMSE", pop, "NoHatch",  paste0("gto3_NoHatch_", pop, file_suffix, ".png")), gto3,
       height = 7, width = 7)
#

