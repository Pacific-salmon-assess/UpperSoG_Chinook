# Reference points

# Libraries
library(patchwork)
library(ggplot2)
library(tidyverse)
library(salmonMSE)
library(here)


# Equilibrium trade-off ref points
EqTradeOff_RP <- read.csv( here(
  "data",
  "Equilibrium trade-off analysis",
  "R-OUT_SMU_ref-pt_values_eq-trade-off.csv")) %>%
  filter(variable == "Smsy") %>%
  select(mid) %>%
  pull %>%
  round()

# Empirical alternative to USR
Emp_RP <- read.csv( here(
  "data",
  "EmpiricalUSR.csv"))

# MSY ref points from conditioning model

CM_RP <- read.csv( here(
  "data",
  "Benchmarks_NatPops.csv")) %>%
  filter(Benchmark == "85% SMSY") %>%
  summarize (RP = sum(med)) %>%
  pull()


# Habitat based ref points

Hab_RP <-  read.csv( here(
  "data",
  "UpperSoGChinook_out_posteriorpredictive.csv"))

Agg_40Smax <- Hab_RP %>%
  summarise(RP = 0.4*sum(SMAX_median)) %>%
  pull(RP) %>%
  round()
Nat_40Smax <- Hab_RP %>%
  filter(Stock %in% c("Nimpkish", "Adam/Eve", "Salmon")) %>%
  summarise(RP = 0.4*sum(SMAX_median)) %>%
  pull(RP) %>%
  round()
Agg_20Smax <- Hab_RP %>%
  summarise(RP = 0.2*sum(SMAX_median)) %>%
  pull(RP) %>%
  round()
Nat_20Smax <- Hab_RP %>%
  filter(Stock %in% c("Nimpkish", "Adam/Eve", "Salmon")) %>%
  summarise(RP = 0.2*sum(SMAX_median)) %>%
  pull(RP) %>%
  round()

df <- data.frame(Population = c(rep("Natural-dominated populations",3),
                 rep("All populations", 3)),
                 Category = c("Data-rich", "Data-rich",
                              "Data-limited", "Data-limited",
                              "Alternative", "Alternative"),
                 Type = c("Sum of 85% SMSY",
                          "Aggregate 85%SMSY (eq. trade-off)",
                          "Sum of 40% SMAX",
                          "Sum of 40% SMAX",
                          "50th percentile",
                          "Average of base period 2002-2011"),
                 Value = c(CM_RP,
                           EqTradeOff_RP,
                           Nat_40Smax,
                           Agg_40Smax,
                           round(Emp_RP$medspawners),
                           round(Emp_RP$avespawners_bp)
                           )
                 )

write.csv(df, here(
  "data",
  "USR.csv"),
row.names = FALSE
)

# UMSY

Umsy <- read.csv( here(
  "data",
  "Equilibrium trade-off analysis",
  "R-OUT_SMU_ref-pt_values_eq-trade-off.csv")) %>%
  filter(variable == "Umsy") %>%
  select(mid) %>%
  pull %>%
  round(2)
