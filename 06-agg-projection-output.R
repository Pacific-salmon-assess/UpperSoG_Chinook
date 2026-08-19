# Libraries
library(salmonMSE)
library(tidyverse)
library(reshape2)
library(here)
library(patchwork)
library(ggplot2)

# Functions for plotting
source("92-decision-table-plots.R") #for alternative version of decision tables

save.files <- TRUE

folder_path <- "figures/SMSE/Aggregate"

if (!dir.exists(here::here(folder_path))) {
  dir.create(here::here(folder_path), recursive = TRUE)
}


# Identify scenarios and management options ----
gr <- readr::read_csv(file.path("tables", "scenarios.csv"))
nOM <- nrow(gr)


# pop <- "QC"

for (pop in c("QC", "Nimpkish", "Salmon", "Adam")){

if(pop == "QC") {
  samp <- readRDS("CM/QuinsamCampbell_07.29.26.rds")
  year1 <- 1984 # For Quinsam/Campbell
}
if(pop == "Nimpkish") {
  samp <- readRDS("CM/Woss_07.22.26.prior.rds")
  year1 <- 2001
}
if(pop == "Salmon") {
  samp <- readRDS("CM/Salmon_08.18.26.prior.rds")
  year1 <- 2002
}
if(pop == "Adam") {
  samp <- readRDS("CM/Adam_08.06.26.prior.rds")
  year1 <- 2002
  }


# Grab all model output ----
# scenario_unique <- unique(gr$Option_name) # Represented by individual table

SMSE_list <- lapply(gr$n, function(i) {
  SMSE <- readRDS(file.path("SMSE", pop, paste0(pop, i, ".rds")))

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

NS <- sapply(SMSE_list, function(x) {
  NS_a <- x@NOS[, 1, , y] + x@HOS[, 1, , y]
  rowSums(NS_a)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, `Natural Spawners` = value, n = Var2)


assign(paste0("NS_", pop), NS)


Ret <- sapply(SMSE_list, function(x) {
  Rec_a <- x@Return_NOS[, 1, , y] + x@Return_HOS[, 1, , y]
  rowSums(Rec_a)
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, `Returns` = value, n = Var2)

assign(paste0("Ret_", pop), Ret)

Tcatch <- sapply(SMSE_list, function(x) { #Terminal catch
  tc_HOS <- apply(x@KT_HOS, c(1, 2, 4), sum) #Sum over ages
  tc_NOS <- apply(x@KT_NOS, c(1, 2, 4), sum)
  tc_HOS[,, y] + tc_NOS[, ,y]
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Tcatch = value, n = Var2)

assign(paste0("Tcatch_", pop), Tcatch)

PTcatch <- sapply(SMSE_list, function(x) {#Pre-terminal catch
  ptc_HOS <- apply(x@KPT_HOS, c(1, 2, 4), sum)
  ptc_NOS <- apply(x@KPT_NOS, c(1, 2, 4), sum)
  ptc_HOS[, ,y] + ptc_NOS[, ,y]
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, PTcatch = value, n = Var2)

assign(paste0("PTcatch_", pop), PTcatch)

Aggcatch <- sapply(SMSE_list, function(x) {#Aggregate catch
  tc_HOS <- apply(x@KT_HOS, c(1, 2, 4), sum)
  tc_NOS <- apply(x@KT_NOS, c(1, 2, 4), sum)
  ptc_HOS <- apply(x@KPT_HOS, c(1, 2, 4), sum)
  ptc_NOS <- apply(x@KPT_NOS, c(1, 2, 4), sum)
  tc_HOS[, ,y] + tc_NOS[, ,y] + ptc_HOS[, ,y] + ptc_NOS[, ,y]
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, Aggcatch = value, n = Var2)

assign(paste0("Aggcatch_", pop), Aggcatch)

AggHOcatch <- sapply(SMSE_list, function(x) {#Aggregate hatchery-origin catch
  tc_HOS <- apply(x@KT_HOS, c(1, 2, 4), sum)
  ptc_HOS <- apply(x@KPT_HOS, c(1, 2, 4), sum)
  tc_HOS[, ,y] + ptc_HOS[, ,y]
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, AggHOcatch = value, n = Var2)

assign(paste0("AggHOcatch_", pop), AggHOcatch)

AggNOcatch <- sapply(SMSE_list, function(x) {#Aggregate natural-origin catch
  tc_NOS <- apply(x@KT_NOS, c(1, 2, 4), sum)
  ptc_NOS <- apply(x@KPT_NOS, c(1, 2, 4), sum)
  tc_NOS[, ,y] + ptc_NOS[, ,y]
}) %>%
  reshape2::melt() %>%
  rename(Simulation = Var1, AggNOcatch = value, n = Var2)

assign(paste0("AggNOcatch_", pop), AggNOcatch)

# Get PT ER for icecream plots
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
assign(paste0("UPTmed_", pop), UPTmed)
assign(paste0("UTmed_", pop), UTmed)


}

# Sum population-specific metrics to natural systems and total systems
# Input aggregate Sgen and %SMSY
# Calc prob of NS > Sgen and 85%SMSY
# icecream plots: Spanwers, Ret, Catch, PT Catch, HO catch, NO Catch
# Add dashed lines for average PT ER, hatch rel =1. no UMSY? or actual total UMSY with * that not comparable to PT ER Ineed  to account for pop specific T er

# Sum variables (third column) among populations
NS_tot <- NS_QC
NS_tot[,3] <- NS_QC[,3] + NS_Nimpkish[,3] + NS_Salmon[,3] + NS_Adam[,3]
NS_nat <- NS_Nimpkish
NS_nat[,3] <- NS_Nimpkish[,3] + NS_Salmon[,3] + NS_Adam[,3]

Ret_tot <- Ret_QC
Ret_tot[,3] <- Ret_QC[,3] + Ret_Nimpkish[,3] + Ret_Salmon[,3] + Ret_Adam[,3]
Ret_nat <- Ret_Nimpkish
Ret_nat[,3] <- Ret_Nimpkish[,3] + Ret_Salmon[,3] + Ret_Adam[,3]

Tcatch_tot <- Tcatch_QC
Tcatch_tot[,3] <- Tcatch_QC[,3] + Tcatch_Nimpkish[,3] + Tcatch_Salmon[,3] +
  Tcatch_Adam[,3]
Tcatch_nat <- Tcatch_Nimpkish
Tcatch_nat[,3] <- Tcatch_Nimpkish[,3] + Tcatch_Salmon[,3] + Tcatch_Adam[,3]

PTcatch_tot <- PTcatch_QC
PTcatch_tot[,3] <- PTcatch_QC[,3] + PTcatch_Nimpkish[,3] + PTcatch_Salmon[,3] +
  PTcatch_Adam[,3]
PTcatch_nat <- PTcatch_Nimpkish
PTcatch_nat[,3] <- PTcatch_Nimpkish[,3] + PTcatch_Salmon[,3] + PTcatch_Adam[,3]

Aggcatch_tot <- Aggcatch_QC
Aggcatch_tot[,3] <- Aggcatch_QC[,3] + Aggcatch_Nimpkish[,3] +
  Aggcatch_Salmon[,3] + Aggcatch_Adam[,3]
Aggcatch_nat <- Aggcatch_Nimpkish
Aggcatch_nat[,3] <- Aggcatch_Nimpkish[,3] + Aggcatch_Salmon[,3] +
  Aggcatch_Adam[,3]

AggHOcatch_tot <- AggHOcatch_QC
AggHOcatch_tot[,3] <- AggHOcatch_QC[,3] + AggHOcatch_Nimpkish[,3] +
  AggHOcatch_Salmon[,3] + AggHOcatch_Adam[,3]
AggHOcatch_nat <- AggHOcatch_Nimpkish
AggHOcatch_nat[,3] <- AggHOcatch_Nimpkish[,3] + AggHOcatch_Salmon[,3] +
  AggHOcatch_Adam[,3]

AggNOcatch_tot <- AggNOcatch_QC
AggNOcatch_tot[,3] <- AggNOcatch_QC[,3] + AggNOcatch_Nimpkish[,3] +
  AggNOcatch_Salmon[,3] + AggNOcatch_Adam[,3]
AggNOcatch_nat <- AggNOcatch_Nimpkish
AggNOcatch_nat[,3] <- AggNOcatch_Nimpkish[,3] + AggNOcatch_Salmon[,3] +
  AggNOcatch_Adam[,3]

# State variables for all simulations in year y in one data frame
val_sim_all_tot <- list(NS_tot, NS_nat, Ret_tot, Ret_nat, Tcatch_tot, Tcatch_nat,
                    PTcatch_tot, PTcatch_nat, Aggcatch_tot, Aggcatch_nat,
                    AggHOcatch_tot, AggHOcatch_nat, AggNOcatch_tot,
                    AggNOcatch_nat) %>%
  Reduce(left_join, .) %>%
  left_join(gr %>% select( n), by = "n") %>%#Option_name,
  # rename(Option = Option_name) %>%
  reshape2::melt(id.vars = c("n", "Simulation"))#, "Option"
if(save.files){
  readr::write_csv(val_sim_all, file = "tables/Agg_outcomes_tot_sim.csv") # Save for Slick object
}

# Median and range for state variables across simulations
val_sim_tot <- val_sim_all_tot %>%
  summarise(m = mean(value),
            median = median(value, na.rm = TRUE),
            lwr = quantile(value, 0.25, na.rm = TRUE),
            upr = quantile(value, 0.75, na.rm = TRUE),
            .by = c( "n", "variable"))


# State variables for all simulations in year y in one data frame
val_sim_all_nat <- list(NS_nat, Ret_nat, Tcatch_nat,
                        PTcatch_nat, Aggcatch_nat,
                        AggHOcatch_nat, AggNOcatch_nat) %>%
  Reduce(left_join, .) %>%
  left_join(gr %>% select( n), by = "n") %>%#Option_name,
  # rename(Option = Option_name) %>%
  reshape2::melt(id.vars = c("n", "Simulation"))#, "Option"
if(save.files){
  readr::write_csv(val_sim_all, file = "tables/Agg_outcomes_nat_sim.csv") # Save for Slick object
}

# Median and range for state variables across simulations
val_sim_nat <- val_sim_all_nat %>%
  summarise(m = mean(value),
            median = median(value, na.rm = TRUE),
            lwr = quantile(value, 0.25, na.rm = TRUE),
            upr = quantile(value, 0.75, na.rm = TRUE),
            .by = c( "n", "variable"))


# MSY ref points from conditioning model

CM_85SMSY <- read.csv( here(
  "data",
  "BenchmarksSimple_NatPops.csv")) %>%
  filter(Benchmark == "85% SMSY") %>%
  summarize (RP = sum(med)) %>%
  pull()

CM_Sgen <- read.csv( here(
  "data",
  "BenchmarksSimple_NatPops.csv")) %>%
  filter(Benchmark == "Sgen") %>%
  summarize (RP = sum(med)) %>%
  pull()


# Calculate probability that natural spawners (from natural-domincated systems)
# are greater than benchmarks
val1 <- val2 <- NA
for(i in 1:dim(gr)[1]){
  NS.x <- NS_nat %>% filter(n == i)
  val1[i] <- mean(NS.x[,3] >= CM_Sgen)
  val2[i]<- mean(NS.x[,3] >= CM_85SMSY)
}
P_NS_Sgen <- data.frame(n = 1:length(val1), P_NS_Sgen = val1)
P_NS_85SMSY <- data.frame(n = 1:length(val2), P_NS_85SMSY = val2)

# val_prob <- data.frame(n = 1:nrow(gr)) %>%
#   mutate(
#     P_NS_Sgen = P_NS_Sgen,
#     P_NS_85SMSY = P_NS_85SMSY
#   ) %>%
#   reshape2::melt(id.vars = "n") %>%
#   left_join(select(gr,  n))
val_prob <- left_join(P_NS_Sgen, P_NS_85SMSY) %>%
  reshape2::melt(id.vars = "n")


if(save.files) {
  readr::write_csv(val_prob, file = "tables/Agg_outcomes_prob.csv") # Save for Slick object
}


# Get UMSY and current U to label on plots

# Average UPT among populations, excluding Adam as PT is not estimated
UPTmed <- (UPTmed_QC + UPTmed_Nimpkish + UPTmed_Salmon)/3

Umsy <- read.csv(  here(
  "data",
  "Equilibrium trade-off analysis",
  "R-OUT_SMU_ref-pt_values_eq-trade-off.csv"
) ) %>%
  filter(variable == "Umsy") %>%
  pull(mid)


CcolSp <- c("0-5,000" = '#d01c8b', "5,001-10,000" = '#f1b6da',
            "10,001-15,000" = '#f7f7f7',
            "15,001-20,000" = '#b8e186', ">20,000" = '#4dac26')

Ccols <- c("0-0.10" = '#d01c8b', "0.11-0.33" = '#f1b6da',
           "0.34-0.65" = '#f7f7f7',
           "0.66-0.89" = '#b8e186', "0.90-1" = '#4dac26')

Ccolc <- c("0-2,000" = '#d01c8b', "2,001-4,000" = '#f1b6da',
           "4,001-6,000" = '#f7f7f7', "6,001-8,000" = '#b8e186',
           ">8,000" = '#4dac26')


Ccolhoc <- c("0-1,000" = '#d01c8b', "1,001-2,000" = '#f1b6da',
           "2,001-3,000" = '#f7f7f7', "3,001-4,000" = '#b8e186',
           ">4,000" = '#4dac26')


Ccolnoc <- c("0-1,500" = '#d01c8b', "1,501-3,000" = '#f1b6da',
             "3,001-4,500" = '#f7f7f7', "4,501-6,000" = '#b8e186',
             ">6,200" = '#4dac26')

CcolTc <- c("0-500" = '#d01c8b', "501-1,000" = '#f1b6da',
            "1,001-1,500" = '#f7f7f7', "1,501-2,000" = '#b8e186',
            ">2,000" = '#4dac26')


g <- val_sim_tot %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Natural Spawners") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value < 5000, "0-5,000",
                        ifelse(value >= 5000 & value <= 10000, "5,001-10,000",
                               ifelse(value > 10000 & value < 15000, "10,001-15,000",
                                      ifelse(value >= 15000 & value < 20000, "15,001-20,000",
                                             ">20,000")))))
# Specify order of legend elements:
g$value <- factor(g$value, levels = c("0-5,000", "5,001-10,000", "10,001-15,000",
                                      "15,001-20,000", ">20,000"))

gSp_tot<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal,  fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = CcolSp,
                    name = "Natural spawners\n(all pops)") +
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

if(save.files){
  ggsave(file.path("figures", "SMSE", "Aggregate", "icecream_NatSp_Agg_tot.png"),
         gSp_tot, width = 7, height = 5)
}

g <- val_sim_nat %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Natural Spawners") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value < 2000, "0-2,000",
                        ifelse(value >= 2000 & value <= 4000, "2,001-4,000",
                               ifelse(value > 4000 & value < 6000, "4,001-6,000",
                                      ifelse(value >= 6000 & value < 8000, "6,001-8,000",
                                             ">8,000")))))
# Specify order of legend elements:
g$value <- factor(g$value, levels = c("0-2,000", "2,001-4,000", "4,001-6,000",
                                      "6,001-8,000", ">8,000"))

gSp_nat<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal,  fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccolc,
                    name = "Natural spawners\n(natural-dominated pops)") +
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

if(save.files){
  ggsave(file.path("figures", "SMSE", "Aggregate", "icecream_NatSp_Agg_nat.png"),
         gSp_nat, width = 7, height = 5)
}

# NS relative to Sgen
g <- val_prob %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "P_NS_Sgen") %>%
  select(u_preterminal, n_yearling, value, n) %>%
  mutate(value = ifelse(value < 0.11, "0-0.10",
                        ifelse(value >= 0.11 & value <= 0.33, "0.11-0.33",
                               ifelse(value > 0.33 & value < 0.66, "0.34-0.65",
                                      ifelse(value >= 0.66 & value < 0.9, "0.66-0.89",
                                             "0.90-1")))))

gSgen<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal,  fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccols,
                    name = "Probability natural\nspawners > Sgen\n(natural-dominated pops)") +
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

ggsave(file.path("figures", "SMSE", "Aggregate", "icecream_Sgen_Agg.png"),
       gSgen, width = 7, height = 5)

# NS relative to 85%SMSY
g <- val_prob %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "P_NS_85SMSY") %>%
  select(u_preterminal, n_yearling, value, n) %>%
  mutate(value = ifelse(value < 0.11, "0-0.10",
                        ifelse(value >= 0.11 & value <= 0.33, "0.11-0.33",
                               ifelse(value > 0.33 & value < 0.66, "0.34-0.65",
                                      ifelse(value >= 0.66 & value < 0.9, "0.66-0.89",
                                             "0.90-1")))))

gSmsy85<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal,  fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccols,
                    name = "Probability natural\nspawners > 85%Smsy\n(natural-dominated pops)") +
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

ggsave(file.path("figures", "SMSE", "Aggregate", "icecream_85Smsy_Agg.png"),
       gSmsy85, width = 7, height = 5)

### Returns decision tables
# Option of icecream plot with legend:

g <- val_sim_tot %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Returns") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value < 5000, "0-5,000",
                        ifelse(value >= 5000 & value <= 10000, "5,001-10,000",
                               ifelse(value > 10000 & value < 15000, "10,001-15,000",
                                      ifelse(value >= 15000 & value < 20000, "15,001-20,000",
                                             ">20,000")))))
# Specify order of legend elements:
g$value <- factor(g$value, levels = c("0-5,000", "5,001-10,000", "10,001-15,000",
                                      "15,001-20,000", ">20,000"))

gRet_tot<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal,  fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = CcolSp,
                    name = "Returns\n(all pops)") +
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

if(save.files){
  ggsave(file.path("figures", "SMSE", "Aggregate", "icecream_Ret_Agg_tot.png"),
         gRet_tot, width = 7, height = 5)
}

g <- val_sim_tot %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Aggcatch") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value < 2000, "0-2,000",
                        ifelse(value >= 2000 & value <= 4000, "2,001-4,000",
                               ifelse(value > 4000 & value < 6000, "4,001-6,000",
                                      ifelse(value >= 6000 & value < 8000, "6,001-8,000",
                                             ">8,000")))))

# Specify order of legend elements:
g$value <- factor(g$value, levels = c("0-2,000", "2,001-4,000", "4,001-6,000",
                                      "6,001-8,000", ">8,000"))

gc_tot<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal,  fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccolc,
                    name = "Total catch\n(all pops)") +
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

if(save.files){
  ggsave(file.path("figures", "SMSE", "Aggregate", "icecream_catch_Agg_tot.png"),
         gc_tot, width = 7, height = 5)
}


g <- val_sim_tot %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "AggHOcatch") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value < 1000, "0-1,000",
                        ifelse(value >= 1000 & value <= 2000, "1,001-2,000",
                               ifelse(value > 2000 & value < 3000, "2,001-3,000",
                                      ifelse(value >= 3000 & value < 4000, "3,001-4,000",
                                             ">4,000")))))

# Specify order of legend elements:
g$value <- factor(g$value, levels = c("0-1,000", "1,001-2,000", "2,001-3,000",
                                      "3,001-4,000", ">4,000"))

ghoc_tot<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal,  fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccolhoc,
                    name = "Total hatchery-\norigin catch\n(all pops)") +
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

if(save.files){
  ggsave(file.path("figures", "SMSE", "Aggregate", "icecream_HOcatch_Agg_tot.png"),
         ghoc_tot, width = 7, height = 5)
}

### Total Catches (preterminal + terminal), natural-origin only
# Option of icecream plot with legend:

g <- val_sim_tot %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "AggNOcatch") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value < 1500, "0-1,500",
                        ifelse(value >= 1500 & value <= 3000, "1,501-3,000",
                               ifelse(value > 3000 & value < 4500, "3,001-4,500",
                                      ifelse(value >= 4500 & value < 6000, "4,501-6,000",
                                             ">6,000")))))

# Specify order of legend elements:
g$value <- factor(g$value, levels = c("0-1,500", "1,501-3,000", "3,001-4,500",
                                      "4,501-6,000", ">6,000"))

gnoc_tot<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal,  fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccolnoc,
                    name = "Total natural-\norigin catch\n(all pops)") +
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

if(save.files){
  ggsave(file.path("figures", "SMSE", "Aggregate",  "icecream_NOcatch_Agg_tot.png"),
         gnoc_tot, width = 7, height = 5)
}

### Pre-terminal Catches
# Option of icecream plot with legend:
g <- val_sim_tot %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "PTcatch") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value < 2000, "0-2,000",
                        ifelse(value >= 2000 & value <= 4000, "2,001-4,000",
                               ifelse(value > 4000 & value < 6000, "4,001-6,000",
                                      ifelse(value >= 6000 & value < 8000, "6,001-8,000",
                                             ">8,000")))))
# Specify order of legend elements:
g$value <- factor(g$value, levels = c("0-2,000", "2,001-4,000", "4,001-6,000",
                                      "6,001-8,000", ">8,000"))

gptc_tot<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal, fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = Ccolc,
                    name = "Pre-terminal\ncatch\n(all pops)") +
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

if(save.files){
  ggsave(file.path("figures", "SMSE", "Aggregate",  "icecream_ptcatch_Agg_tot.png"),
         gptc_tot, width = 7, height = 5)
}

### Terminal Catches
# Option of icecream plot with legend:

g <- val_sim_tot %>%
  left_join(select(gr, u_preterminal, n_yearling, n)) %>%
  filter(variable == "Tcatch") %>%
  select(u_preterminal, n_yearling, median, n) %>%
  rename(value = median) %>%
  mutate(value = ifelse(value <= 500, "0-500",
                        ifelse(value > 500 & value <= 1000, "501-1,000",
                               ifelse(value > 1000 & value <= 1500, "1,001-1,500",
                                      ifelse(value > 1500 & value <= 2000, "1,501-2,000",
                                             ">2,000")))))
# Specify order of legend elements:
g$value <- factor(g$value, levels = c("0-500", "501-1,000", "1,001-1,500",
                                      "1,501-2,000", ">2,000"))

gtc_tot<- g %>% ggplot(aes(x = n_yearling, y = u_preterminal,  fill = value, z = value)) +
  geom_raster() +
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = CcolTc,
                    name = "Terminal\ncatch\n(all pops)") +
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

if(save.files){
  ggsave(file.path("figures", "SMSE", "Aggregate", "icecream_tcatch_Agg_tot.png"),
         gtc_tot, width = 7, height = 5)
}

# Combine ice-cream plots

gic1 <- (gSp_tot + gSp_nat)/
  (gSgen + gSmsy85)

gic2 <- (gRet_tot  + gc_tot) /
  (gptc_tot + gtc_tot)

gic3 <- (gnoc_tot + ghoc_tot)/
  plot_spacer() +
  plot_layout(heights = c(1, 1))


if(save.files){
  ggsave(file.path("figures", "SMSE", "Aggregate", "gic1_Agg.png"),
         gic1, height = 9, width = 7)
  ggsave(file.path("figures", "SMSE", "Aggregate", "gic2_Agg.png"),
         gic2, height = 9, width = 7)
  ggsave(file.path("figures", "SMSE", "Aggregate", "gic3_Agg.png"),
         gic3, height = 9, width = 7)

}
