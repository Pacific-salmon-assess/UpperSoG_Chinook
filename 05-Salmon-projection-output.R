# Libraries
library(salmonMSE)
library(tidyverse)
library(reshape2)
library(here)
library(patchwork)
library(ggplot2)

pop <- "Salmon"
#---------------------------------------------------------------------------
# Note: Update Sgen and 85%SMSY values below if the CM changes
#---------------------------------------------------------------------------

# Functions for plotting
source("92-decision-table-plots.R") #for alternative version of decision tables

# Get CM output for current PT and T ERs
samp <-  readRDS("CM/Salmon_08.18.26.prior.rds")
year1 <- 2002 # For Salmon
scenario <- "base"
if (scenario !="base") {file_suffix <- paste0("_", scenario)} else {file_suffix <- ""}

folder_path <- "figures/SMSE/Salmon"

if (!dir.exists(here::here(folder_path))) {
  dir.create(here::here(folder_path), recursive = TRUE)
}


# Identify scenarios and management options ----
gr <- readr::read_csv(file.path("tables", "scenarios.csv"))
nOM <- nrow(gr)

# Grab all model output ----
# scenario_unique <- unique(gr$Option_name) # Represented by individual table

SMSE_list <- lapply(gr$n, function(i) {
  SMSE <- readRDS(file.path("SMSE",pop, paste0(pop, i, ".rds")))

  # Update PNI = 1 when there is no brood & pHOS = 0
  Brood <- SMSE@HOB[,,5,] + SMSE@NOB[,,5,]
  # Brood is age-selective. Take age-5s for purposes of setting PNI to 1
  NoBrood <- Brood < 0.001
  pHOS_zero <- SMSE@pHOS_effective < 0.001
  # if(any(NoBrood)) print(i)
  # if(any(pHOS_zero)) print(i)
  SMSE@PNI[pHOS_zero] <- 1
  SMSE@PNI[NoBrood] <- 1

  return(SMSE)
})

# Performance metrics (all simulations at end of projection) ----
y <- 29

PNI <- sapply(SMSE_list, function(x) x@PNI[, 1, y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, PNI = value, n = Var2)

NS <- sapply(SMSE_list, function(x) {
  NS_a <- x@NOS[, 1, , y] + x@HOS[, 1, , y]
  rowSums(NS_a)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, `Natural Spawners` = value, n = Var2)

Ret <- sapply(SMSE_list, function(x) {
  Rec_a <- x@Return_NOS[, 1, , y] + x@Return_HOS[, 1, , y]
  rowSums(Rec_a)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, `Returns` = value, n = Var2)


NOS <- sapply(SMSE_list, function(x) {
  NS_a <- x@NOS[, 1, , y]
  rowSums(NS_a)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, `Natural-Origin Spawners` = value, n = Var2)


# MA <- sapply(SMSE_list, function(x) {
#   NS_a <- x@NOS[, 1, , y] + x@HOS[, 1, , y]
#   MA <- apply(NS_a, 1, function(w) weighted.mean(x = 1:5, w = w))
#   return(MA)
# }) %>%
#   reshape2::melt() %>%
#   rename(Simulation = Var1, `Mean age` = value, n = Var2)

p_wild <- sapply(SMSE_list, function(x) x@p_wild[, , y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, pWILD = value, n = Var2)

pHOSeff <- sapply(SMSE_list, function(x) x@pHOS_effective[, , y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, pHOSeff = value, n = Var2)

pNOBeff <- sapply(SMSE_list, function(x) x@pNOB[, , y] + x@pNOB[, , y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, pNOBeff = value, n = Var2)

Brood <- sapply(SMSE_list, function(x) x@NOB[, , 4, y] + x@NOB[, , 5, y] +
                  x@HOB[, , 4, y] + x@HOB[, , 5, y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Brood = value, n = Var2)

# IRR <- sapply(SMSE_list, function(x) {
#   apply(x@Escapement_NOS[, , , y] + x@Escapement_HOS[, , , y], 1, sum)
# }) %>%
#   reshape2::melt() %>%
#   rename(Simulation = Var1, IR_Return = value, n = Var2)


Tcatch <- sapply(SMSE_list, function(x) { #Terminal catch
  tc_HOS <- apply(x@KT_HOS, c(1, 2, 4), sum) #Sum over ages
  tc_NOS <- apply(x@KT_NOS, c(1, 2, 4), sum)
  tc_NOS[, ,y] + tc_HOS[,, y]
  }) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Tcatch = value, n = Var2)

PTcatch <- sapply(SMSE_list, function(x) {#Pre-terminal catch
  ptc_HOS <- apply(x@KPT_HOS, c(1, 2, 4), sum)
  ptc_NOS <- apply(x@KPT_NOS, c(1, 2, 4), sum)
  ptc_NOS[, ,y] + ptc_HOS[, ,y]
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, PTcatch = value, n = Var2)

Aggcatch <- sapply(SMSE_list, function(x) {#Aggregate catch
  tc_HOS <- apply(x@KT_HOS, c(1, 2, 4), sum)
  tc_NOS <- apply(x@KT_NOS, c(1, 2, 4), sum)
  ptc_HOS <- apply(x@KPT_HOS, c(1, 2, 4), sum)
  ptc_NOS <- apply(x@KPT_NOS, c(1, 2, 4), sum)
  tc_NOS[, ,y]  + ptc_NOS[, ,y] + tc_HOS[, ,y]+ ptc_HOS[, ,y]
  }) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Aggcatch = value, n = Var2)

AggHOcatch <- sapply(SMSE_list, function(x) {#Aggregate hatchery-origin catch
  tc_HOS <- apply(x@KT_HOS, c(1, 2, 4), sum)
  ptc_HOS <- apply(x@KPT_HOS, c(1, 2, 4), sum)
  tc_HOS[, ,y] + ptc_HOS[, ,y]
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, AggHOcatch = value, n = Var2)

AggNOcatch <- sapply(SMSE_list, function(x) {#Aggregate natural-origin catch
  tc_NOS <- apply(x@KT_NOS, c(1, 2, 4), sum)
  ptc_NOS <- apply(x@KPT_NOS, c(1, 2, 4), sum)
  tc_NOS[, ,y] + ptc_NOS[, ,y]
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, AggNOcatch = value, n = Var2)


Egg <- sapply(SMSE_list, function(x) x@Egg_NOS[, 1, y] + x@Egg_HOS[, 1, y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Egg = value, n = Var2) %>%
  mutate(Egg = Egg/1e6)

Rel <- sapply(SMSE_list, function(x) x@Smolt_Rel[, 1, y]) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Releases = value, n = Var2) %>%
  mutate(Releases = Releases/1e5)

#Option_name <- data.frame(
#  Option = gr$Option[1:9]
#) %>%
#  mutate(scenario = paste0("(", 1:9, ") ", Option))

# State variables for all simulations in year y in one data frame
val_sim_all <- list(PNI, NS, NOS, pNOBeff, pHOSeff, p_wild, #MA,
                    Brood, PTcatch, Tcatch, Aggcatch, AggHOcatch, Ret,
                    AggNOcatch, Egg, Rel) %>%
Reduce(left_join, .) %>%
  left_join(gr %>% select( n), by = "n") %>%#Option_name,
  # rename(Option = Option_name) %>%
  reshape2::melt(id.vars = c("n", "Simulation"))#, "Option"
readr::write_csv(val_sim_all, file = paste0("tables/", pop, "_outcomes_sim", file_suffix, ".csv")) # Save for Slick object

# Median and range for state variables across simulations
val_sim <- val_sim_all %>%
  summarise(m = mean(value),
            median = median(value, na.rm = TRUE),
            lwr = quantile(value, 0.25, na.rm = TRUE),
            upr = quantile(value, 0.75, na.rm = TRUE),
            .by = c( "n", "variable")) #%>%  #"Option",
#left_join(Option_name, by = "Option") %>%
#mutate(scenario = factor(scenario, rev(Option_name$scenario)))

# Probability (end of projection) ----
  SMSY85 <- read.csv( here(
    "data",
    "BenchmarksSimple_NatPops.csv")) %>%
    filter(Benchmark == "85% SMSY") %>%
    filter(Pop == pop) %>% select(med) %>% pull()
  Sgen <- read.csv( here(
    "data",
    "BenchmarksSimple_NatPops.csv")) %>%
    filter(Benchmark == "Sgen") %>%
    filter(Pop == pop) %>% select(med) %>% pull()

#Sgen = 730
  P_Sgen_NOS <- sapply(SMSE_list, function(x, Sgen = 129) {
    val <- rowSums(x@NOS[, 1, , y]) >= 129
    mean(val)
  })
#SMSY85 =1744
  P_Smsy85_NOS <- sapply(SMSE_list, function(x, SMSY = 251) {
    val <- rowSums(x@NOS[, 1, , y]) >= 251
    mean(val)
  })

  P_Sgen <- sapply(SMSE_list, function(x, Sgen = 129) {
    val <- rowSums(x@NOS[, 1, , y] + x@HOS[, 1, , y]) >= 129
    mean(val)
  })

  P_Smsy85 <- sapply(SMSE_list, function(x, SMSY = 251) {
    val <- rowSums(x@NOS[, 1, , y] + x@HOS[, 1, , y]) >=  251
    mean(val)
  })


P_1500 <- sapply(SMSE_list, function(x) {
  val <- rowSums(x@NOS[, 1, , y] + x@HOS[, 1, , y]) >= 1500
  mean(val)
})

P_PNI50 <- sapply(SMSE_list, function(x) {
  PNI <- x@PNI[, 1, y]
  mean(PNI >= 0.5, na.rm = TRUE)
})


# Check proportion of simulations where PNI is not defined
# (No brood, but intermittently due to brood rule requiring > 600 spawners)
# This comes out to between 0-5%, not a big deal
if (FALSE) {
  PNI_NA <- sapply(SMSE_list, function(x) {
    PNI <- x@PNI[, 1, y]
    mean(is.na(PNI))
  })
  PNI_NA
}

  val_prob <- data.frame(n = 1:nrow(gr)) %>%
    mutate(
      P_PNI50 = P_PNI50,
      # P_Sgen_NOS = P_Sgen_NOS,
      # P_Smsy85_NOS = P_Smsy85_NOS,
      P_Sgen_NS = P_Sgen,
      P_Smsy85_NS = P_Smsy85,
      P_1500_NS = P_1500
    ) %>%
    reshape2::melt(id.vars = "n") %>%
    left_join(select(gr,  n)) #%>%  #Option_name,
    # rename(Option = Option_name)






readr::write_csv(val_prob, file = paste0("tables/", pop, "_outcomes_prob", file_suffix, ".csv") ) # Save for Slick object

### Big data frame of state variables for each simulation and year (for shinysalmon app)
### CSV file is 200 MB, but we'll have to convert the data frame to arrays and save as an R object (.rds) later to reduce disk space
# Not implemetned
if(FALSE){
  state_var2 <- c("Egg_NOS", "Egg_HOS", "Smolt_Rel", "Smolt_NOS", "Smolt_HOS", "KPT_NOS", "KPT_HOS", "Return_NOS", "Return_HOS",
                  "KT_NOS", "KT_HOS",
                  "Escapement_NOS", "Escapement_HOS", "NOB", "HOB",  "NOS", "HOS", "pNOB", "PNI", "pHOS_effective")

  df <- lapply(1:length(state_var2), function(j) {
    lapply(1:length(SMSE_list), function(i) {
      .ts_fn(SMSE_list[[i]], var = state_var2[j], all_sims = TRUE) %>%
        reshape2::melt() %>%
        rename(Simulation = Var1, Year = Var2) %>%
        mutate(variable = state_var2[j], n = gr$n[i])
    }) %>%
      bind_rows()
  }) %>%
    bind_rows() %>%
    left_join(select(gr,  n), by = "n")# %>% Option_name,
  # rename(Option = Option_name)
  readr::write_csv(df, file = paste0("tables/Quinsam_outcomes_sim_year_app", file_suffix, ".csv")) # Save for Slick object
  rm(df)
}

#### Run loop over each scenario to create performance metric tables ----

# Time series medians
#for (i in 1:length(scenario_unique)) {
#
#  ind <- gr$Scenario == scenario_unique[i]
#
#  dir <- file.path("figures", "SMSE", paste0("Set", i, "_"))
#
#  #### Time series barplots (annual medians for by scenario for each management option)
#  ind2 <- which(ind)
#
#  png(paste0(dir, "spawners_prop.png"), height = 6, width = 6, units = "in", res = 400)
#  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
#  for (j in 1:nrow(Option_name)) {
#    jj <- ind2[j]
#    plot_spawners(SMSE_list[[jj]])
#    box()
#    title(Option_name$scenario[j], cex.main = 0.75)
#  }
#  dev.off()
#
#  png(paste0(dir, "spawners.png"), height = 6, width = 6, units = "in", res = 400)
#  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
#  for (j in 1:nrow(Option_name)) {
#    jj <- ind2[j]
#    plot_spawners(SMSE_list[[jj]], prop = FALSE, ylim = c(0, 4000))
#    box()
#    title(Option_name$scenario[j], cex.main = 0.75)
#  }
#  dev.off()
#
#  png(paste0(dir, "p_brood.png"), height = 6, width = 6, units = "in", res = 400)
#  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
#  for (j in 1:nrow(Option_name)) {
#    jj <- ind2[j]
#    plot_escapement(SMSE_list[[jj]], ylim = c(0, 1))
#    box()
#    title(Option_name$scenario[j], cex.main = 0.75)
#  }
#  dev.off()
#
#  png(paste0(dir, "pHOS_fitness.png"), height = 6, width = 6, units = "in", res = 400)
#  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
#  for (j in 1:nrow(Option_name)) {
#    jj <- ind2[j]
#    plot_fitness(SMSE_list[[jj]], ylim = c(0, 1))
#    title(Option_name$scenario[j], cex.main = 0.75)
#  }
#  dev.off()
#
#  png(paste0(dir, "RS_HOS.png"), height = 6, width = 6, units = "in", res = 400)
#  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
#  for (j in 1:nrow(Option_name)) {
#    jj <- ind2[j]
#    plot_RS(SMSE_list[[jj]], var = "HOS", type = "abs",
#            name = c("Fed Fry", "Traditionals"), ylim = c(0, 4000))
#    title(Option_name$scenario[j], cex.main = 0.75)
#  }
#  dev.off()
#
#  png(paste0(dir, "RS_Rel.png"), height = 6, width = 6, units = "in", res = 400)
#  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
#  for (j in 1:nrow(Option_name)) {
#    jj <- ind2[j]
#    plot_RS(SMSE_list[[jj]], var = "Smolt", type = "abs", name = c("Fed Fry", "Traditionals"),
#            ylab = "Hatchery releases")
#    title(Option_name$scenario[j], cex.main = 0.75)
#  }
#  dev.off()
#
#  png(paste0(dir, "RS_Esc.png"), height = 6, width = 6, units = "in", res = 400)
#  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
#  for (j in 1:nrow(Option_name)) {
#    jj <- ind2[j]
#    plot_RS(SMSE_list[[jj]], var = "Esc", type = "abs", name = c("Fed Fry", "Traditionals"),
#            ylab = "HO Escapement",
#            ylim = c(0, 5000))
#    title(Option_name$scenario[j], cex.main = 0.75)
#  }
#  dev.off()
#
#  png(paste0(dir, "LHG_NOS.png"), height = 6, width = 6, units = "in", res = 400)
#  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
#  for (j in 1:nrow(Option_name)) {
#    jj <- ind2[j]
#    plot_LHG(SMSE_list[[jj]], var = "NOS", type = "abs", name = c("Early Smalls", "Late Larges"),
#             ylab = "NOS",
#             ylim = c(0, 700))
#    title(Option_name$scenario[j], cex.main = 0.75)
#  }
#  dev.off()
#
#  png(paste0(dir, "LHG_Esc.png"), height = 6, width = 6, units = "in", res = 400)
#  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
#  for (j in 1:nrow(Option_name)) {
#    jj <- ind2[j]
#    plot_LHG(SMSE_list[[jj]], var = "Esc", type = "abs", name = c("Early Smalls", "Late Larges"),
#             ylab = "NO Escapement",
#             ylim = c(0, 900))
#    title(Option_name$scenario[j], cex.main = 0.75)
#  }
#  dev.off()
#
#  png(paste0(dir, "LHG_Smolt.png"), height = 6, width = 6, units = "in", res = 400)
#  par(mfrow = c(3, 3), mar = c(5, 4, 1, 1))
#  for (j in 1:nrow(Option_name)) {
#    jj <- ind2[j]
#    plot_LHG(SMSE_list[[jj]], var = "Smolt", type = "abs", name = c("Early Smalls", "Late Larges"),
#             ylab = "NO outmigrating juveniles",
#             ylim = c(0, 4e5))
#    title(Option_name$scenario[j], cex.main = 0.75)
#  }
#  dev.off()
#
#}


### Tables of performance metrics by freshwater survival
# Not implemented
if(FALSE){
  pm_primary <- c("PNI", "Natural Spawners", "P_PNI50", "P_1500_NS")
  pm_ancillary <- c("IR_Return", "Brood", "Egg", "Releases",
                    "pNOBeff", "pHOSeff", "pWILD", "P_Sgen_NOS", "P_Smsy85_NOS", "P_Sgen_NS", "P_Smsy85_NS")

  for (i in LETTERS[1:3]) {

    # Make figure of performance metrics (all simulations at end of projection) ----
    #g <- val_sim %>%
    #  filter(Scenario == scenario_unique[i]) %>%
    #  left_join(gr) %>%
    #  #filter(variable %in% pm_primary) %>%
    #  plot_dotplot() +
    #  scale_shape_manual(values = c(1, 4, 16)) +
    #  #theme(strip.placement = "outside") +
    #  theme(legend.position = "bottom") +
    #  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
    #  labs(x = NULL, y = NULL, shape = "IRER 1300", colour = "pNOB target") +
    #  ggtitle(scenario_unique[i]) +
    #  scale_x_discrete(labels = bold_scenario) +
    #  geom_vline(xintercept = c(3, 6) + 0.5)
    #ggsave(paste0(dir, "performance_metrics.png"), g, width = 6, height = 7)


    # Full performance metrics (medians at end of the projection) ----
    d <- rbind(
      val_sim %>% select(Option, Scenario, variable, median),
      val_prob %>% select(Option, Scenario, variable, value) %>% rename(median = value)
    ) %>%
      filter(grepl(i, Scenario)) %>%
      filter(variable %in% c(pm_primary, pm_ancillary)) %>%
      mutate(variable = factor(variable, c(pm_primary, pm_ancillary)))

    g <- plot_table(d, ncol = 1) +
      geom_vline(xintercept = length(pm_primary) + 0.5, linewidth = 1, linetype = 2) +
      scale_y_discrete(labels = font_fn, limits = rev) +
      geom_hline(yintercept = 3.5, linewidth = 1)
    ggsave(file.path("figures", "SMSE", paste0("performance_table_", i, file_suffix, ".png")), g, width = 5.5, height = 7)
  }

  # Full performance table all scenarios
  glist <- lapply(LETTERS[1:3], function(i) {
    d <- rbind(
      val_sim %>% select(Option, Scenario, variable, median),
      val_prob %>% select(Option, Scenario, variable, value) %>% rename(median = value)
    ) %>%
      filter(grepl(i, Scenario)) %>%
      filter(variable %in% c(pm_primary, pm_ancillary)) %>%
      mutate(variable = factor(variable, c(pm_primary, pm_ancillary)))

    g <- plot_table(d, ncol = 1) +
      geom_vline(xintercept = length(pm_primary) + 0.5, linewidth = 1, linetype = 2) +
      geom_hline(yintercept = 3.5, linewidth = 1)

    if (i == "A") {
      g <- g + scale_y_discrete(labels = font_fn, limits = rev)
    } else {
      g <- g + scale_y_discrete(labels = NULL, limits = rev)
    }

    g
  })

  g <- ggpubr::ggarrange(plotlist = glist, ncol = 3, widths = c(3, 2, 2))
  ggsave(paste0("figures/SMSE/performance_table_full", file_suffix, ".png"), g, width = 17, height = 8)

  # Short performance metrics
  glist <- lapply(LETTERS[1:3], function(i) {
    d <- rbind(
      val_sim %>% select(Option, Scenario, variable, median),
      val_prob %>% select(Option, Scenario, variable, value) %>% rename(median = value)
    ) %>%
      filter(grepl(i, Scenario)) %>%
      filter(variable %in% pm_primary) %>%
      mutate(variable = factor(variable, pm_primary))

    g <- plot_table(d, ncol = 1) +
      geom_vline(xintercept = c(2.5, 4.5), linewidth = 1, linetype = 2) +
      geom_hline(yintercept = 3.5, linewidth = 1)

    if (i == "A") {
      g <- g + scale_y_discrete(labels = font_fn, limits = rev)
    } else {
      g <- g + scale_y_discrete(labels = NULL, limits = rev)
    }

    g
  })

  g <- ggpubr::ggarrange(plotlist = glist, ncol = 3, widths = c(3, 2, 2))
  ggsave(paste0("figures/SMSE/performance_table_short", file_suffix, ".png"), g, width = 8, height = 7)

}

#### Decision table for all scenarios ----

### PNI decision tables
# Option of icecream plot with legend:
Ccols <- c("0-0.10" = '#d01c8b', "0.11-0.33" = '#f1b6da',
           "0.34-0.65" = '#f7f7f7',
           "0.66-0.89" = '#b8e186', "0.90-1" = '#4dac26')

CcolSp <- c("0-200" = '#d01c8b', "201-400" = '#f1b6da',
            "401-600" = '#f7f7f7', "601-800" = '#b8e186',
            ">800" = '#4dac26')


CcolRel <- c("0-1" = '#d01c8b',
            "1.1-2" = '#f1b6da',
            "2.1-3" = '#f7f7f7', "3.1-4" = '#b8e186',
            ">4" = '#4dac26')

CcolPNI <- c("0-0.30" = '#d01c8b', "0.31-0.50" = '#f1b6da',
             "0.51-0.70" = '#f7f7f7', "0.71-0.79" = '#b8e186',
             "0.80-1.0"= '#4dac26')

CcolTc <- c("0-10" = '#d01c8b', "11-20" = '#f1b6da',
             "21-30" = '#f7f7f7', "31-40" = '#b8e186',
             ">40" = '#4dac26')


Ccolc <- c("0-100" = '#d01c8b', "101-200" = '#f1b6da',
             "201-300" = '#f7f7f7', "301-400" = '#b8e186',
             ">400" = '#4dac26')

Ccol400 <- c("0-100" = '#d01c8b', "101-200" = '#f1b6da',
           "201-300" = '#f7f7f7', "301-400" = '#b8e186',
           ">400" = '#4dac26')

Ccol100 <- c("0-25" = '#d01c8b', "26-50" = '#f1b6da',
             "51-75" = '#f7f7f7', "76-100" = '#b8e186',
             ">100" = '#4dac26')


Ccol40 <- c("0-10" = '#d01c8b', "11-20" = '#f1b6da',
            "21-30" = '#f7f7f7', "31-40" = '#b8e186',
            ">40" = '#4dac26')

Ccolhoc <- c("0-50" = '#d01c8b', "51-100" = '#f1b6da',
           "101-150" = '#f7f7f7', "151-200" = '#b8e186',
           ">200" = '#4dac26')


CcolTc <- c("0-5" = '#d01c8b', "6-10" = '#f1b6da',
           "11-15" = '#f7f7f7', "16-20" = '#b8e186',
           ">20" = '#4dac26')


CcolRet <- c("0-400" = '#d01c8b', "401-800" = '#f1b6da',
            "801-1200" = '#f7f7f7', "1201-1600" = '#b8e186',
            ">1600" = '#4dac26')


# Get UMSY and current U to label on plots
report <- salmonMSE:::get_report(samp)
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)

year <- year1 + seq(1, d$Ldyr) - 1
year_borrow <- seq(max(year) - 9, max(year) - 5)

UT <- salmonMSE:::.CM_ER(report, brood = FALSE, type = "T",  r = 1,
                         index_AEQ = match(year_borrow, year)) %>% t()
UPT <- salmonMSE:::.CM_ER(report, brood = FALSE, type = "PT",  r = 1,
                          index_AEQ = match(year_borrow, year)) %>% t()
UTrange <-  UT %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975), na.rm=T) %>%
  t()
UTmed <- UTrange[seq((d$Ldyr - 4),d$Ldyr), "50%"] %>%
  mean()
UPTrange <-  UPT %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975), na.rm=T) %>%
  t()
UPTmed <- UPTrange[seq((d$Ldyr - 4),d$Ldyr), "50%"] %>%
  mean()

Umsy <- read.csv(  here(
  "data",
  "Equilibrium trade-off analysis",
  "R-OUT_SMU_ref-pt_values_eq-trade-off.csv"
) ) %>%
  filter(variable == "Umsy") %>%
  pull(mid)

g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "PNI") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value <= 0.3, "0-0.30",
                        ifelse(value > 0.30 & value <= 0.5, "0.31-0.50",
                               ifelse(value > 0.5 & value <= 0.7, "0.51-0.70",
                                      ifelse(value > 0.7 & value < 0.80, "0.71-0.79",
                                             ifelse(value >= 0.8 & value <= 1, "0.80-1.0", NA))))))

gPNI <- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value, z = value)) +
  geom_raster() +
  scale_fill_manual(values = CcolPNI,
                    name = "PNI") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x = "Proportional change in\nhatchery releases", y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", pop, paste0("icecream_PNI_", pop, file_suffix, ".png")), gPNI, width = 5, height = 5)


# # Option of icecream plot without legend and numbers on the plot:
#
# g <- val_sim %>%
#   left_join(select(gr, u_preterminal, n_yearling, n)) %>%
#   filter(variable == "PNI") %>%
#   select(u_preterminal, n_yearling, median, n) %>%
#   rename(value = median) %>%
#   decision_table_grid(
#     ncol = 3,
#     title = "Median PNI",
#     fill_scheme =
#       scale_fill_gradientn(
#         values = c(0, 0.7, 1),
#         colours = c("deeppink", "lightgreen", "green4")
#       )
#   )
# g$facet$params$free$y <- TRUE
# g$facet$params$free$x <- TRUE
# # ggsave(file.path("figures", "SMSE", "decisiontable_PNI_QC.png"), g, width = 5, height = 5)
#
#
# Option of icecream plot with legend:

g <- val_prob %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "P_PNI50") %>%
  select(u_preterminal, n_yearling, value, n) %>%
  mutate(value = ifelse(value < 0.11, "0-0.10",
                          ifelse(value >= 0.11 & value <= 0.33, "0.11-0.33",
                                 ifelse(value > 0.33 & value < 0.66, "0.34-0.65",
                                        ifelse(value >= 0.66 & value < 0.9, "0.66-0.89",
                                               "0.90-1")))))

gPNIprob <- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value,
                             z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccols,
                    name = "Probability\nof PNI > 0.5") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x = "Proportional change in\nhatchery releases", y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")

ggsave(file.path("figures", "SMSE", pop, paste0("icecream_probPNI_", pop,
                                                file_suffix, ".png")),
                 gPNIprob, width = 5, height = 5)


# # Option of icecream plot without legend and numbers on the plot:
#
# val <- seq(0, 1, 0.01)
# cols <- ifelse(val >= 0.5, "lightgreen", "deeppink")
# g <- val_prob %>%
#   left_join(select(gr, u_preterminal, n_yearling, n)) %>%
#   filter(variable == "P_PNI50") %>%
#   select(u_preterminal, n_yearling, value, n) %>%
#   decision_table_grid(
#     ncol = 3,
#     title = "Probability PNI > 0.5",
#     fill_scheme =
#       scale_fill_gradientn(
#         values = val,
#         colours = cols
#       )
#   )
# g$facet$params$free$y <- TRUE
# g$facet$params$free$x <- TRUE
# # ggsave(file.path("figures", "SMSE", "decisiontable_PNI50_QC.png"), g, width = 5, height = 5)

### Natural spawners decision tables
# Option of icecream plot with legend:

g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Natural Spawners") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median)

if(scenario == "base" | scenario =="high"){
  g <- g %>%
    mutate(value = ifelse(value < 200, "0-200",
                          ifelse(value >= 200 & value <= 400, "201-400",
                                 ifelse(value > 400 & value < 600, "401-600",
                                        ifelse(value >= 600 & value < 800, "601-800",
                                               ">800")))))
  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-200", "201-400", "401-600",
                                        "601-800", ">800"))
  col <- CcolSp
}

if(scenario == "low"){
  g <- g %>%
    mutate(value = ifelse(value <= 100, "0-100",
                          ifelse(value > 100 & value <= 200, "101-200",
                                 ifelse(value > 200 & value <= 300, "201-300",
                                        ifelse(value > 300 & value <= 400, "301-400",
                                               ">400")))))
  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-100", "101-200", "201-300",
                                        "301-400", ">400"))
  col <- Ccol400
}

gSp<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = col,
                    name = "Natural spawners") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x = "Proportional change in\nhatchery releases", y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", pop, paste0("icecream_NatSp_", pop,
                                                file_suffix, ".png")),
       gSp, width = 5, height = 5)


# option for icream plot with numbers on plot (no legend)
g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Natural Spawners") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  mutate(value = round(median)) %>%
  decision_table_grid(ncol = 3, title = "Natural spawners",
                      fill_scheme = scale_fill_gradient2(low = "deeppink",
                                                         high = "green4",
                                                         mid = "white",
                                                         midpoint = 1500))
g$facet$params$free$y <- TRUE
g$facet$params$free$x <- TRUE
# ggsave(file.path("figures", "SMSE", "decisiontable_sp_QC.png"), g, width = 5, height = 5)

### NS <1500 prob decision tables
# Option of icecream plot without legend:
g <- val_prob %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "P_1500_NS") %>%
  select(u_preterminal, n_yearling, value, n) %>%
  decision_table_grid(
    ncol = 3,
    title = "Probabilty > 1500 natural spawners",
    fill_scheme =
      scale_fill_gradientn(
        values = c(0, 0.7, 1),
        colours = c("deeppink", "lightgreen", "green4")
      )
  )
g$facet$params$free$y <- TRUE
g$facet$params$free$x <- TRUE
# ggsave(file.path("figures", "SMSE", "decisiontable_P_1500_QC.png"), g, width = 5, height = 5)

# Option of icecream plot with legend:
g <- val_prob %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "P_1500_NS") %>%
  select(u_preterminal, n_yearling, value, n) %>%
  mutate(value = ifelse(value < 0.11, "0-0.10",
                        ifelse(value >= 0.11 & value <= 0.33, "0.11-0.33",
                               ifelse(value > 0.33 & value < 0.66, "0.34-0.65",
                                      ifelse(value >= 0.66 & value < 0.9, "0.66-0.89",
                                             "0.90-1")))))

gNSprob<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value,
                           z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccols,
                    name = "Probability natural\nspawners > 1500") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x ="Proportional change in\nhatchery releases" , y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", paste0("icecream_NOSProb1500_", pop,
                                           file_suffix, ".png")),
       gNSprob, width = 5, height = 5)

# NS relative to Sgen
g <- val_prob %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "P_Sgen_NS") %>%
  select(u_preterminal, n_yearling, value, n) %>%
  mutate(value = ifelse(value < 0.11, "0-0.10",
                        ifelse(value >= 0.11 & value <= 0.33, "0.11-0.33",
                               ifelse(value > 0.33 & value < 0.66, "0.34-0.65",
                                      ifelse(value >= 0.66 & value < 0.9, "0.66-0.89",
                                             "0.90-1")))))

gSgen<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value,
                         z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccols,
                    name = "Probability natural\nspawners > Sgen") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x ="Proportional change in\nhatchery releases" , y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", pop, paste0("icecream_Sgen_", pop,
                                                file_suffix, ".png")),
       gSgen, width = 5, height = 5)

# NS relative to 85%SMSY
g <- val_prob %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "P_Smsy85_NS") %>%
  select(u_preterminal, n_yearling, value, n) %>%
  mutate(value = ifelse(value < 0.11, "0-0.10",
                        ifelse(value >= 0.11 & value <= 0.33, "0.11-0.33",
                               ifelse(value > 0.33 & value < 0.66, "0.34-0.65",
                                      ifelse(value >= 0.66 & value < 0.9, "0.66-0.89",
                                             "0.90-1")))))

gSmsy85<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value,
                           z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccols,
                    name = "Probability natural\nspawners > 85%Smsy") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x ="Proportional change in\nhatchery releases" , y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", pop, paste0("icecream_Smsy_", pop,
                                                file_suffix, ".png")),
       gSmsy85, width = 5, height = 5)

### Returns decision tables
# Option of icecream plot with legend:

g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Returns") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median)

if(scenario == "base" | scenario =="high") {
  g <- g %>%
    mutate(value = ifelse(value < 200, "0-200",
                          ifelse(value >= 200 & value <= 400, "201-400",
                                 ifelse(value > 400 & value < 600, "401-600",
                                        ifelse(value >= 600 & value < 800, "601-800",
                                               ">800")))))
  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-200", "201-400", "401-600",
                                        "601-800", ">800"))
  col <- CcolSp
}
if(scenario == "low"){
  g <- g %>%
    mutate(value = ifelse(value <= 100, "0-100",
                          ifelse(value > 100 & value <= 200, "101-200",
                                 ifelse(value > 200 & value <= 300, "201-300",
                                        ifelse(value > 300 & value <= 400, "301-400",
                                               ">400")))))
  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-100", "101-200", "201-300",
                                        "301-400", ">400"))
  col <- Ccol400
}

gRet<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value,
                        z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = col,
                    name = "Returns") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x = "Proportional change in\nhatchery releases", y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", pop, paste0("icecream_Ret_", pop,
                                                file_suffix, ".png")),
       gRet, width = 5, height = 5)


### Hatchery releases
# Option of icecream plot without legend:

g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Releases") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  decision_table_grid(ncol = 3, "Hatchery releases\n(100,000s)")
g$facet$params$free$y <- TRUE
g$facet$params$free$x <- TRUE
#ggsave(file.path("figures", "SMSE", "decisiontable_rel.png"), g, width = 5, height = 5)


# Option of icecream plot with legend:
g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Releases") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value <= 1, "0-1",
                        ifelse(value > 1 & value <= 2, "1.1-2",
                               ifelse(value > 2 & value < 3, "2.1-3",
                                      ifelse(value > 2 & value < 3, "3.1-4",
                                             ">4")))) )
# Specify order of legend elements:
g$value <- factor(g$value, levels = c("0-1", "1.1-2", "2.1-3",
                                      "3.1-4", ">4"))

gRel<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value,
                        z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = CcolRel,
                    name = "Hatchery releases\n(100,000s)") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x = "Proportional change in\nhatchery releases", y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")

ggsave(file.path("figures", "SMSE", pop, paste0("icecream_Rel_", pop,
                                                file_suffix, ".png")),
       gRel, width = 5, height = 5)


### Total Catches (preterminal + terminal)
# Option of icecream plot with legend:

g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Aggcatch") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median)

if(scenario == "base" | scenario == "high"){
  g <- g %>%
    mutate(value = ifelse(value < 100, "0-100",
                          ifelse(value >= 100 & value <= 200, "101-200",
                                 ifelse(value > 200 & value < 300, "201-300",
                                        ifelse(value >= 300 & value < 400, "301-400",
                                               ">400")))))

  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-100", "101-200", "201-300",
                                        "301-400", ">400"))
  col <- Ccolc

}
if(scenario == "low"){
  g <- g %>%
    mutate(value = ifelse(value <= 25, "0-25",
                          ifelse(value > 25 & value <= 50, "26-50",
                                 ifelse(value > 50 & value <= 75, "51-75",
                                        ifelse(value > 75 & value <= 100, "76-100",
                                               ">100")))))

  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-25", "26-50", "51-75",
                                        "76-100", ">100"))
  col <- Ccol100

}

gc<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value,
                      z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = col,
                    name = "Total catch") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x = "Proportional change in\nhatchery releases", y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", pop, paste0("icecream_catch_", pop,
                                                file_suffix, ".png")),
       gc, width = 5, height = 5)


### Total Catches (preterminal + terminal), hatchery-orgin only
# Option of icecream plot with legend:

g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "AggHOcatch") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median)
if(scenario == "base" | scenario == "high"){
  g <- g %>%
    mutate(value = ifelse(value < 50, "0-50",
                          ifelse(value >= 50 & value <= 100, "51-100",
                                 ifelse(value > 100 & value < 150, "101-150",
                                        ifelse(value >= 150 & value < 200, "151-200",
                                               ">200")))))

  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-50", "51-100", "101-150",
                                        "151-200", ">200"))
  col <- Ccolhoc
}

if(scenario == "low"){
  g <- g %>%
    mutate(value = ifelse(value <= 10, "0-10",
                          ifelse(value > 10 & value <= 20, "11-20",
                                 ifelse(value > 20 & value <= 30, "21-30",
                                        ifelse(value > 30 & value <= 40, "31-40",
                                               ">40")))))

  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-10", "11-20", "21-30",
                                        "31-40", ">40"))
  col <- Ccol40
}

ghoc<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value,
                        z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = col,
                    name = "Total hatchery-\norigin catch") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x = "Proportional change in\nhatchery releases", y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", pop, paste0("icecream_HOcatch_", pop,
                                                file_suffix, ".png")),
       ghoc, width = 5, height = 5)


### Total Catches (preterminal + terminal), natural-origin only
# Option of icecream plot with legend:

g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "AggNOcatch") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median)
if(scenario == "base" | scenario =="high"){
  g <- g %>%
    mutate(value = ifelse(value < 50, "0-50",
                          ifelse(value >= 50 & value <= 100, "51-100",
                                 ifelse(value > 100 & value < 150, "101-150",
                                        ifelse(value >= 150 & value < 200, "151-200",
                                               ">200")))))

  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-50", "51-100", "101-150",
                                        "151-200", ">200"))
  col <- Ccolhoc

}
if(scenario == "low"){
  g <- g %>%
    mutate(value = ifelse(value <= 10, "0-10",
                          ifelse(value > 10 & value <= 20, "11-20",
                                 ifelse(value > 20 & value <= 30, "21-30",
                                        ifelse(value > 30 & value <= 40, "31-40",
                                               ">40")))))

  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-10", "11-20", "21-30",
                                        "31-40", ">40"))
  col <- Ccol40

}
gnoc<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value,
                        z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = col,
                    name = "Total natural-\norigin catch") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x = "Proportional change in\nhatchery releases", y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", pop, paste0("icecream_NOcatch_", pop,
                                                file_suffix, ".png")),
       gnoc, width = 5, height = 5)

### Pre-terminal Catches
# Option of icecream plot with legend:
g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "PTcatch") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median)

if(scenario == "base" | scenario == "high"){
  g <- g %>%
    mutate(value = ifelse(value < 100, "0-100",
                          ifelse(value >= 100 & value <= 200, "101-200",
                                 ifelse(value > 200 & value < 300, "201-300",
                                        ifelse(value >= 300 & value < 400, "301-400",
                                               ">400")))))

  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-100", "101-200", "201-300",
                                        "301-400", ">400"))
  col <- Ccolc
}

if(scenario == "low"){
  g <- g %>%
    mutate(value = ifelse(value <= 25, "0-25",
                          ifelse(value > 25 & value <= 50, "26-50",
                                 ifelse(value > 50 & value <= 75, "51-75",
                                        ifelse(value > 75 & value <= 100, "76-100",
                                               ">100")))))

  # Specify order of legend elements:
  g$value <- factor(g$value, levels = c("0-25", "26-50", "51-75",
                                        "76-100", ">100"))
  col <- Ccol100

}


gptc<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal,fill = value,
                        z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = col,
                    name = "Pre-terminal\ncatch") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x = "Proportional change in\nhatchery releases", y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", pop, paste0("icecream_ptcatch_", pop,
                                                file_suffix, ".png")),
       gptc, width = 5, height = 5)


### Terminal Catches
# Option of icecream plot with legend:

g <- val_sim %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Tcatch") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value < 5, "0-5",
                        ifelse(value >= 5 & value <= 10, "6-10",
                               ifelse(value > 10 & value < 15, "11-15",
                                      ifelse(value >= 15 & value < 20, "16-20",
                                             ">20")))))

# Specify order of legend elements:
g$value <- factor(g$value, levels = c("0-5", "6-10", "11-15",
                                      "16-20", ">20"))

gtc<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value,
                       z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = CcolTc,
                    name = "Terminal\ncatch") +
  theme(legend.position = "top", legend.text = element_text(size = 10),
        legend.title = element_text(size = 13)) +
  guides(fill = guide_legend(ncol = 2), title.position = "top") +
  labs(x = "Proportional change in\nhatchery releases", y = "Pre-terminal ER") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.key = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  geom_hline(yintercept = (Umsy - UTmed), col="grey40") +
  geom_hline(yintercept = (UPTmed), col="grey40", linetype = "dashed") +
  geom_vline(xintercept = 1, col="grey40", linetype="dashed") +
  # annotate(geom="text", x=(Umsy - UTmed), y=Inf, hjust= -0.1, vjust= - 1,
  #          label="Umsy (proxy)",
  #          colour="grey40", size=3) +
  coord_cartesian(clip = "off")
ggsave(file.path("figures", "SMSE", pop, paste0("icecream_tcatch_", pop,
                                                file_suffix, ".png")),
       gtc, width = 5, height = 5)



# Combine ice-cream plots

gic1 <- (gSp + gRet)/
  (gSgen + gSmsy85)

gic2 <- (gPNI + gc)/
  (ghoc + gnoc)

gic3 <- (gptc + gtc)/
  plot_spacer() +
  plot_layout(heights = c(1, 1))




ggsave(file.path("figures", "SMSE", pop, paste0("gic1_", pop, file_suffix,
                                                ".png")), gic1,
       height = 9, width = 7)
ggsave(file.path("figures", "SMSE", pop, paste0("gic2_", pop, file_suffix,
                                                ".png")), gic2,
       height = 9, width = 7)
ggsave(file.path("figures", "SMSE", pop, paste0("gic3_", pop, file_suffix,
                                                ".png")), gic3,
       height = 7, width = 7)




#### Trade-off figures
# PNI vs natural spawners
PNI_NS <- val_sim %>%
  filter(variable %in% c("PNI", "Natural Spawners")) %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  reshape2::dcast(u_preterminal + n_yearling ~ variable, value.var = "median") %>%
  mutate(
    bin = cut(
      n_yearling, breaks = seq(0.5, 1.5, length.out = 5),
      labels = c("50-75%", "75-100%", "100-125%", "125-150%"),
      include.lowest = TRUE, right = TRUE
    )
  )
gto <- plot_tradeoff(
  pm1 = PNI_NS$`Natural Spawners`,
  pm2 = PNI_NS$PNI,
  x1 = PNI_NS$u_preterminal,
  x2 = PNI_NS$bin,
  xlab = "Natural spawners (median)",
  ylab = "PNI (median)",
  x1lab = "Exploitation\nrate",
  x2lab = "Relative\nhatchery\nproduction"
)
ggsave(file.path("figures", "SMSE", pop, paste0("gto_", pop, file_suffix, ".png")), gto,
       height = 9, width = 7)



PNI_NS <- val_sim %>%
  filter(variable %in% c("PNI", "Natural Spawners")) %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  reshape2::dcast(u_preterminal + n_yearling ~ variable, value.var = "median") %>%
  filter(u_preterminal %in% c(0.00, 0.10, 0.20, 0.30, 0.40, 0.50)) %>%
  filter(n_yearling %in% c(0.50, 0.75, 1.00, 1.25, 1.50))
  # mutate(
  #   bin = cut(
  #     g$n_yearling, breaks = seq(0.5, 1.5, length.out = 5),
  #     labels = c("50-75%", "75-100%", "100-125%", "125-150%"),
  #     include.lowest = TRUE, right = TRUE
  #   )
  # )
gto <- plot_tradeoff(
  pm1 = PNI_NS$`Natural Spawners`,
  pm2 = PNI_NS$PNI,
  x1 = PNI_NS$u_preterminal,
  x2 = PNI_NS$n_yearling,
  xlab = "Natural spawners (median)",
  ylab = "PNI (median)",
  x1lab = "Exploitation\nrate",
  x2lab = "Relative\nhatchery\nproduction"
)
ggsave(file.path("figures", "SMSE", pop, paste0("gto_", pop, file_suffix, "v2.png")), gto,
       height = 9, width = 7)



# #### Plot Simulated CYER- Quang's archived plots
# CYER_PT <- calc_CYER(SMSE_list[[1]], PT = TRUE)
# CYER_T <- calc_CYER(SMSE_list[[1]], PT = FALSE)
#
# png("simulated_CYER.png", height = 6, width = 5, units = "in", res = 400)
# par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
# matplot(t(CYER_PT), type = "l", lty = 1, xlab = "Projection Year", ylab = "Preterminal CYER")
# matplot(t(CYER_T), type = "l", lty = 1, xlab = "Projection Year", ylab = "Terminal CYER")
# dev.off()


# #### Plot maturity and exploitation rate at age ----
# png("figures/SMSE/maturity_ER.png", height = 6, width = 3, units = "in", res = 400)
# par(mfrow = c(3, 1), mar = c(5, 4, 1, 1))
#
# SOM <- SMSE_list[[1]]@Misc$SOM
#
# salmonMSE:::plot_Mjuv_RS(SOM@Hatchery[[1]]@p_mature_HOS[, , 1, ],
#                          RS_names = c("Fed Fry", "Traditionals"), ylab = "Proportion mature")
#
# SOM@Harvest[[1]]@vulPT <- SMSE_list[[1]]@ExPT_NOS[, 1, , 30]
#
# salmonMSE:::plot_SOM(SOM@Harvest[[1]], "vulPT",
#                      type = "age", nsim = SOM@nsim, maxage = SOM@Bio[[1]]@maxage,
#                      proyears = SOM@proyears,
#                      ylab = "Preterminal exploitation rate")
#
# SOM@Harvest[[1]]@vulT <- SMSE_list[[1]]@ExT_NOS[, 1, , 30]
# salmonMSE:::plot_SOM(SOM@Harvest[[1]], "vulT",
#                      type = "age", nsim = SOM@nsim, maxage = SOM@Bio[[1]]@maxage,
#                      proyears = SOM@proyears,
#                      ylab = "Terminal exploitation rate")
#
# dev.off()
#
#
