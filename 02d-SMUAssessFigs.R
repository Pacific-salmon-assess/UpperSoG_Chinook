#-------------------------------------------------------------------------------
# Plot SMU level assessment figures for FSAR
#-------------------------------------------------------------------------------

# INPUTS:
# Conditioning model results

# Population productivity (logalpha): base, low or high

# Aggregated benchmarks for natural-dominated populations
# UMSY: "tables/R-OUT_SMU_ref-pt_values_eq-trade-off.csv"

# Aggregated lower benchmarks:   'medlBench_Nat'
# Aggregated upper benchmarks from Equilibrium Trade-off Analysis in file:
#   "tables/R-OUT_SMU_ref-pt_values_eq-trade-off.csv"

# OUTPUTS:
# Empirical threshold for aggregate abundances
# File names:
# "tables/EmpiricalUSR.low.csv"
# "tables/EmpiricalUSR.base.csv"
# "tables/EmpiricalUSR.high.csv"

# Four-panel figure for FSAR
# File names:
# "figures/Fourpanel.base.png",
# "figures/Fourpanel.low.png",
# "figures/Fourpanel.high.png"
# time-series of empirical values for projection model
# Files names: "Woss_timeseries.csv", "Salmon_timeseries.csv", "Adam_timeseries.csv"
  # in folders: data/CMoutput/base, data/CMoutput/low, or data/CMoutput/high

#-------------------------------------------------------------------------------

# Libraries
library(patchwork)
library(ggplot2)
library(tidyverse)
library(salmonMSE)
library(here)

#-------------------------------------------------------------------------------
# INPUTS:

prod <- "ensemble"#"base","low", "high", "ensemble"

# Get aggregated lower benchmark:   'medlBench_Nat'
# source(here("02-PopBenchFig.R"))

# Get aggregated upper benchmark:
meduBench_Nat <- read.csv(
  here(
  "tables",
  paste0("R-OUT_SMU_ref-pt_values_eq-trade-off.", prod, ".csv")
  )
) %>% filter(variable == "Smsy") %>% pull(mid)



# Input Conditioning Model results
if (prod!="ensemble"){
  ERM_QC <- readRDS(paste0("CM/QuinsamCampbell_09.22.26.", prod, ".rds"))
  report_QC <- salmonMSE:::get_report(ERM_QC)

  ERM_Adam <- readRDS(paste0("CM/Adam_09.26.26.", prod, ".rds"))
  report_Adam <- salmonMSE:::get_report(ERM_Adam)

  ERM_Salmon <- readRDS(paste0("CM/Salmon_09.26.26.", prod, ".rds"))
  report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

  ERM_Woss <- readRDS(paste0("CM/Woss_09.26.26.", prod, ".rds"))
  report_Woss <- salmonMSE:::get_report(ERM_Woss)
}

if(prod == "ensemble"){
  # Get number of MC trials (len)
  ERM_len <- readRDS(paste0("CM/QuinsamCampbell_09.22.26.base.rds"))
  report_len <- salmonMSE:::get_report(ERM_len)
  len <- length(report_len)

  # Pull reports for individual runs, for low, base, and high prod scenarios
  ERM_QC_low <- readRDS(paste0("CM/QuinsamCampbell_09.22.26.low.rds"))
  ERM_QC_base <- readRDS(paste0("CM/QuinsamCampbell_09.22.26.base.rds"))
  ERM_QC_high <- readRDS(paste0("CM/QuinsamCampbell_09.22.26.high.rds"))
  # Combine reports, 25% low, 50% base, and 25% high productivity
  report_QC <- c(sample( get_report( ERM_QC_low), len * 0.25),
                 sample( get_report( ERM_QC_base), len * 0.5),
                 sample( get_report( ERM_QC_high), len * 0.25)
  )
  ERM_QC <- ERM_QC_base

  # Pull reports for individual runs, for low, base, and high prod scenarios
  ERM_Woss_low <- readRDS(paste0("CM/Woss_09.26.26.low.rds"))
  ERM_Woss_base <- readRDS(paste0("CM/Woss_09.26.26.base.rds"))
  ERM_Woss_high <- readRDS(paste0("CM/Woss_09.26.26.high.rds"))
  # Combine reports, 25% low, 50% base, and 25% high productivity
  report_Woss <- c(sample( get_report( ERM_Woss_low), len * 0.25),
                   sample( get_report( ERM_Woss_base), len * 0.5),
                   sample( get_report( ERM_Woss_high), len * 0.25)
  )
  ERM_Woss <- ERM_Woss_base

  # Pull reports for individual runs, for low, base, and high prod scenarios
  ERM_Salmon_low <- readRDS(paste0("CM/Salmon_09.26.26.low.rds"))
  ERM_Salmon_base <- readRDS(paste0("CM/Salmon_09.26.26.base.rds"))
  ERM_Salmon_high <- readRDS(paste0("CM/Salmon_09.26.26.high.rds"))
  # Combine reports, 25% low, 50% base, and 25% high productivity
  report_Salmon <- c(sample( get_report( ERM_Salmon_low), len * 0.25),
                     sample( get_report( ERM_Salmon_base), len * 0.5),
                     sample( get_report( ERM_Salmon_high), len * 0.25)
  )
  ERM_Salmon <- ERM_Salmon_base

  # Pull reports for individual runs, for low, base, and high prod scenarios
  ERM_Adam_low <- readRDS(paste0("CM/Adam_09.26.26.low.rds"))
  ERM_Adam_base <- readRDS(paste0("CM/Adam_09.26.26.base.rds"))
  ERM_Adam_high <- readRDS(paste0("CM/Adam_09.26.26.high.rds"))
  # Combine reports, 25% low, 50% base, and 25% high productivity
  report_Adam <- c(sample( get_report( ERM_Adam_low), len * 0.25),
                   sample( get_report( ERM_Adam_base), len * 0.5),
                   sample( get_report( ERM_Adam_high), len * 0.25)
  )
  ERM_Adam <- ERM_Adam_base

}


folder_path <- paste0("data/CMoutput/", prod)

if (!dir.exists(here::here(folder_path))) {
  dir.create(here::here(folder_path), recursive = TRUE)
}

# Make folder for CM results, and put csv file there...


for (pop in c("QC", "Woss", "Adam", "Salmon")){

  samp <- get(paste0("ERM_", pop))
  d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
  # brood <- FALSE
  # type <- "T"# "PT" #
  if(pop == "QC") year1 <- 1984
  if(pop == "Woss") year1 <- 2001
  if(pop == "Salmon") year1 <- 2002
  if(pop == "Adam") year1 <- 2002

  report <- get(paste0("report_", pop))

  year <- year1 + seq(1, d$Ldyr) - 1
  year_borrow <- seq(max(year) - 9, max(year) - 5)


  catchPT <- sapply(report, getElement, "catchPT") %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    t()
  catchT <- sapply(report, getElement, "catchT") %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    t()


  #spawners
  esc <- sapply(report, getElement, "spawners") %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    t()


  ER <- salmonMSE:::.CM_ER(report, type = "all", r = 1, brood = FALSE,
                           index_AEQ = match(year_borrow, year)) %>%
    apply(1,quantile, probs = c(0.025, 0.5, 0.975), na.rm=T) %>%
    t()


  recr.out <- sapply(report, getElement, "recr")
  # Extracted as vector across years, ages and NO vs HO
  # First create array for every Posterior draw
  splitrow <- function(x){array(x, dim=c(length(year), d$Nages, 2))}

  recr <- seq_len(ncol(recr.out)) %>%
    map(~ splitrow(recr.out[,.x])) %>%
    map (~ .x[,,1] + .x[,,2]) %>% # Sum HO and NO
    map(~ .x[,1] + .x[,2] + .x[,3] + .x[,4] + .x[,5] ) %>% # Sum over ages
    simplify2array() %>%
    apply(FUN = quantile, MARGIN = 1, probs = c(0.025, 0.5, 0.975)) %>% # Extract quantiles
    t()

  df <- data.frame(year = rep(year, 3),
                   catchPT = c(catchPT[,'2.5%'], catchPT[,'50%'], catchPT[,'97.5%']),
                   catchT = c(catchT[,'2.5%'], catchT[,'50%'], catchT[,'97.5%']),
                   esc = c(esc[,'2.5%'], esc[,'50%'], esc[,'97.5%']),
                   ER = c(ER[,'2.5%'], ER[,'50%'], ER[,'97.5%']),
                   recr = c(recr[,'2.5%'], recr[,'50%'], recr[,'97.5%']),
                   label = c( rep("2.5%", length(year)),
                              rep("50%", length(year)),
                              rep("97.5%", length(year))
                   )
  )
  write.table(df, file =
                here::here(paste0("data/CMoutput/", prod, "/", pop, "_timeseries.csv")))

  df_withpop <- df %>% mutate(Pop = pop)
  if(pop == "Woss") df_withpop <- df %>% mutate (Pop = "Nimpkish")
  assign(paste0("Timeseries_", pop), df_withpop)

}

ts <- rbind(Timeseries_QC, Timeseries_Adam, Timeseries_Salmon, Timeseries_Woss)


# Dataframe for catches

# #To get population-specific plots:
# ts <- ts %>% filter(Pop == "Woss")

ts_catchPT <- ts %>%
  select(c(year, catchPT, label, Pop)) %>%
  pivot_wider(names_from=c(Pop,label), values_from=c(catchPT)) %>%
  filter(year >= 2002) %>%
  mutate(lwr =  rowSums(select(., contains("2.5%")), na.rm = TRUE)) %>%
  mutate(med =  rowSums(select(., contains("50%")), na.rm = TRUE)) %>%
  mutate(upr =  rowSums(select(., contains("97.5%")), na.rm = TRUE)) %>%
  select(c(year, lwr, med, upr)) %>%
  mutate(label="PreterminalCatch")

ts_catchT <- ts %>%
  select(c(year, catchT, label, Pop)) %>%
  pivot_wider(names_from=c(Pop,label), values_from=c(catchT)) %>%
  filter(year >= 2002) %>%
  mutate(lwr =  rowSums(select(., contains("2.5%")), na.rm = TRUE)) %>%
  mutate(med =  rowSums(select(., contains("50%")), na.rm = TRUE)) %>%
  mutate(upr =  rowSums(select(., contains("97.5%")), na.rm = TRUE)) %>%
  select(c(year, lwr, med, upr)) %>%
  mutate(label="TerminalCatch")

ts_catchAll <- ts %>% mutate(catchAll = catchPT + catchT) %>%
  select(c(year, catchAll, label, Pop)) %>%
  pivot_wider(names_from=c(Pop,label), values_from=c(catchAll)) %>%
  filter(year >= 2002) %>%
  mutate(lwr =  rowSums(select(., contains("2.5%")), na.rm = TRUE)) %>%
  mutate(med =  rowSums(select(., contains("50%")), na.rm = TRUE)) %>%
  mutate(upr =  rowSums(select(., contains("97.5%")), na.rm = TRUE)) %>%
  select(c(year, lwr, med, upr)) %>%
  mutate(label="TotalCatch")

catchPT_bp <- ts_catchPT %>%
  filter(year>=2002 & year <=2011) %>%
  pull(med) %>%
  mean()

catchT_bp <- ts_catchT %>%
  filter(year>=2002 & year <=2011) %>%
  pull(med) %>%
  mean()

catchAll_bp <- ts_catchAll %>%
  filter(year>=2002 & year <=2011) %>%
  pull(med) %>%
  mean()

catch <- rbind(ts_catchPT, ts_catchT, ts_catchAll) %>%
  mutate(med=med/1000, upr=upr/1000, lwr=lwr/1000)

# Dataframe for spawners

# For population specific spawner plots
# ts <- ts %>% filter(Pop == "Woss")
# ts_Nsp <- NULL

ts_Tsp <- ts %>% select (c(year, esc, label, Pop)) %>%
  pivot_wider(names_from=c(Pop,label), values_from=c(esc)) %>%
  filter(year >= 2002) %>%
  mutate(lwr =  rowSums(select(., contains("2.5%")), na.rm = TRUE)) %>%
  mutate(med =  rowSums(select(., contains("50%")), na.rm = TRUE)) %>%
  mutate(upr =  rowSums(select(., contains("97.5%")), na.rm = TRUE)) %>%
  select(c(year, lwr, med, upr)) %>%
  mutate(label="Total")

# Spawners from natural-dominated systems
ts_Nsp <- ts %>% select (c(year, esc, label, Pop)) %>%
  filter(Pop != "QC") %>%
  pivot_wider(names_from=c(Pop,label), values_from=c(esc)) %>%
  filter(year >= 2002) %>%
  mutate(lwr =  rowSums(select(., contains("2.5%")), na.rm = TRUE)) %>%
  mutate(med =  rowSums(select(., contains("50%")), na.rm = TRUE)) %>%
  mutate(upr =  rowSums(select(., contains("97.5%")), na.rm = TRUE)) %>%
  select(c(year, lwr, med, upr)) %>%
  mutate(label="Natural")


spawners <- rbind(ts_Tsp, ts_Nsp) %>%
  mutate(med=med/1000, upr=upr/1000, lwr=lwr/1000)

# Calculated empirical alternatives to USR

medspawners <- spawners %>%
  filter(label == "Total") %>%
  pull(med) %>%
  median()*1000

spawners_bp <- spawners %>%
  filter(label == "Total") %>%
  filter(year < 2012 & year > 2001) %>%
  select(med)
avespawners_bp <- mean(spawners_bp$med)*1000

empUSR <- data.frame(medspawners = medspawners, avespawners_bp = avespawners_bp)

write.csv(
  empUSR,
  here(
    "tables",
    paste0("EmpiricalUSR.", prod, ".csv")
  ),
  row.names = FALSE
)



# Dataframe for ER

ts_ER <- ts %>% select (c(year, ER, label, Pop)) %>%
  pivot_wider(names_from=c(label), values_from=c(ER))
ts_ER <- ts_ER %>%
  mutate(lwr =  ts_ER$'2.5%') %>%
  mutate(med =  ts_ER$'50%') %>%
  mutate(upr =  ts_ER$'97.5%') %>%
  select(c(year, Pop, lwr, med, upr)) %>%
  mutate(
    Pop = if_else(Pop == "QC", "Quinsam/Campbell", Pop)
  )

# Dataframe for recruitment
# For population specific recruitment plots
# ts <- ts %>% filter(Pop == "Woss")
# ts_Nrec <- NULL

ts_Trec <- ts %>% select (c(year, recr, label, Pop)) %>%
  pivot_wider(names_from=c(Pop,label), values_from=c(recr)) %>%
  filter(year >= 2002) %>%
  mutate(lwr =  rowSums(select(., contains("2.5%")), na.rm = TRUE)) %>%
  mutate(med =  rowSums(select(., contains("50%")), na.rm = TRUE)) %>%
  mutate(upr =  rowSums(select(., contains("97.5%")), na.rm = TRUE)) %>%
  select(c(year, lwr, med, upr)) %>%
  mutate(label="Total")

# Recruitment from natural-dominated systems
ts_Nrec <- ts %>% select (c(year, recr, label, Pop)) %>%
  filter(Pop != "QC") %>%
  pivot_wider(names_from=c(Pop,label), values_from=c(recr)) %>%
  filter(year >= 2002) %>%
  mutate(lwr =  rowSums(select(., contains("2.5%")), na.rm = TRUE)) %>%
  mutate(med =  rowSums(select(., contains("50%")), na.rm = TRUE)) %>%
  mutate(upr =  rowSums(select(., contains("97.5%")), na.rm = TRUE)) %>%
  select(c(year, lwr, med, upr)) %>%
  mutate(label="Natural")


recruits <- rbind(ts_Trec, ts_Nrec) %>%
  mutate(med=med/1000, upr=upr/1000, lwr=lwr/1000)

#------------------------------------------------------------------------------
# Plots

# Catch
gcatch <- ggplot(catch, aes(x = year, y = med, group = label)) +
         geom_line(aes(colour = label), linewidth = 0.8) +
         geom_ribbon(aes(ymin = lwr, ymax = upr, fill = label), alpha = 0.2, colour=NA) +
        # geom_ribbon(aes(ymin = lwr, ymax = upr, colour = label), alpha = 0.2, colour=NA) +
  ylab("Catch ('000s)") +
  xlab("Year") +
  xlab("Year") +
  annotate("text", x=-Inf, y=Inf , label = "a", hjust= 3, vjust= 1.2, size =5) +
  coord_cartesian(clip = "off") +
  # guides(linewidth = "none") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(
    axis.text = element_text(colour = "grey40"),
    axis.title = element_text(colour = "grey40"),
    axis.line = element_line(colour = "grey40"),
    axis.ticks = element_line(colour = "grey40"),
    panel.border = element_rect(colour = "grey40")
  ) +
  theme(panel.border = element_rect(
    colour = "grey40")) +
  theme(
    legend.position = c(0.05, 0.65),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  guides(colour = guide_legend(reverse = TRUE),
         fill = guide_legend(reverse = TRUE)) #+
  # geom_hline(yintercept= catchAll_bp/1000,
  #            colour= "#619CFF",
  #            linetype= "dashed",
  #            linewidth=0.3) +
  # geom_hline(yintercept= catchPT_bp/1000,
  #            colour= "#F8766D",
  #            linetype= "dashed",
  #            linewidth=0.3) +
  # geom_hline(yintercept= catchT_bp/1000,
  #            colour= "#00BA38",
  #            linetype= "dashed",
  #            linewidth=0.3)

gcatch

# Black and White is too hard to interpret with 3 lines and 3 bands
if (FALSE){
  ggplot(catch, aes(x = year, y = med, group = label)) +
    geom_line(aes(linetype = label, linewidth = label)) +
    scale_linetype_manual(
      values = c(
        TotalCatch = "solid",
        PreterminalCatch = "dashed",
        TerminalCatch = "dotted"
      )
    ) +
    scale_linewidth_manual(
      values = c(
        TotalCatch = 1,
        PreterminalCatch = 0.5,
        TerminalCatch = 0.5
      )
    ) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = "black") +
    # geom_ribbon(aes(ymin = lwr, ymax = upr, colour = label), alpha = 0.2, colour=NA) +
    ylab("Catch ('000s)") +
    xlab("Year") +
    theme_classic() +
    theme(legend.title = element_blank())
}

# Spawners
gspawners <- ggplot(spawners, aes(x = year, y = med, group = label)) +
  geom_line(aes(linewidth = label, linetype = label)) +
  scale_linetype_manual(
    values = c(
      Total = "solid",
      Natural = "dashed"
    )
  ) +
  scale_linewidth_manual(
    values = c(
      Total = 1,
      Natural = 0.75
    )
  ) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, #colour = "black",
              linetype = "solid", linewidth = 0.1) +
  ylab("Spawners ('000s)") +
  xlab("Year") +
   annotate("text", x=-Inf, y=Inf , label = "b", hjust= 3, vjust= 1.2, size =5) +
  coord_cartesian(clip = "off") +
  guides(linewidth = "none") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(
    axis.text = element_text(colour = "grey40"),
    axis.title = element_text(colour = "grey40"),
    axis.line = element_line(colour = "grey40"),
    axis.ticks = element_line(colour = "grey40"),
    panel.border = element_rect(colour = "grey40")
  ) +
  theme(panel.border = element_rect(
    colour = "grey40")) +
  theme(
    legend.position = c(0.05, 0.75),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  guides(
    linetype = guide_legend(
      override.aes = list(linewidth = 0.75),
      reverse = TRUE
      )
  ) +
  geom_hline(yintercept= medlBench_Nat/1000, colour= "darkred", linetype= "dotted", linewidth=0.7) +
  geom_hline(yintercept= meduBench_Nat/1000, colour= "darkgreen", linetype= "dashed", linewidth=0.7) +
  geom_hline(yintercept= medspawners/1000, colour= "navyblue", linetype= "solid", linewidth=0.5, alpha=0.2)



  # + theme(legend.position = "none")

gspawners

# Add 85%SMY and Sgen for natural dominated systems

# ER

UMSY_out <- read.csv( here(
  "tables",
  paste0("R-OUT_SMU_ref-pt_values_eq-trade-off.", prod, ".csv")))
UMSY_med <- UMSY_out %>% filter(variable == "Umsy") %>% pull(mid)
UMSY_lwr <- UMSY_out %>% filter(variable == "Umsy") %>% pull(lwr)
UMSY_upr <- UMSY_out %>% filter(variable == "Umsy") %>% pull(upr)

gER <- ggplot(ts_ER, aes(x = year, y = med, group = Pop, colour = Pop,
                            fill = Pop)) +
  geom_line(size = 1) +
  # geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, colour=NA) +
  # geom_ribbon(aes(ymin = lwr, ymax = upr, colour = label), alpha = 0.2, colour=NA) +
  ylab("Exploitation Rate") +
  xlab("Year") +
  annotate("text", x=-Inf, y=Inf , label = "c", hjust= 4, vjust= 1.2, size =5) +
  coord_cartesian(clip = "off", xlim = c(1987, 2025)) +
  guides(linewidth = "none") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(
    axis.text = element_text(colour = "grey40"),
    axis.title = element_text(colour = "grey40"),
    axis.line = element_line(colour = "grey40"),
    axis.ticks = element_line(colour = "grey40"),
    panel.border = element_rect(colour = "grey40")
  ) +
  theme(panel.border = element_rect(
    colour = "grey40")) +
  theme(
    legend.position = c(0.05, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  geom_hline(yintercept = UMSY_med, colour = "grey40", linetype = "dashed",
             linewidth = 0.6)+
  annotate(
    "rect",
    xmin = 1988, xmax = 2024,
    ymin = UMSY_lwr, ymax = UMSY_upr,
    fill = "darkgrey",
    color = NA,
    alpha = 0.2
  )


# + theme(legend.position = "none")
gER
 # Add UMSY for natural dominated systems

# Recruitment (Total, natural domianated systems)

# Spawners
grec <- ggplot(recruits, aes(x = year, y = med, group = label)) +
  geom_line(aes(linewidth = label, linetype = label)) +
  scale_linetype_manual(
    values = c(
      Total = "solid",
      Natural = "dashed"
    )
  ) +
  scale_linewidth_manual(
    values = c(
      Total = 1,
      Natural = 0.75
    )
  ) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, #colour = "black",
              linetype = "solid", linewidth = 0.1) +
  ylab("Recruits ('000s)") +
  xlab("Year") +
  annotate("text", x=-Inf, y=Inf , label = "d", hjust= 3, vjust= 1.2, size =5) +
  coord_cartesian(clip = "off") +
  guides(linewidth = "none") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(
    axis.text = element_text(colour = "grey40"),
    axis.title = element_text(colour = "grey40"),
    axis.line = element_line(colour = "grey40"),
    axis.ticks = element_line(colour = "grey40"),
    panel.border = element_rect(colour = "grey40")
  ) +
  theme(panel.border = element_rect(
    colour = "grey40")) +
  theme(
    legend.position = c(0.05, 0.75),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank()
  ) +
  guides(
    linetype = guide_legend(
      override.aes = list(linewidth = 0.75),
      reverse = TRUE
    )
  )

# + theme(legend.position = "none")

grec


g4panel <- (gcatch + gspawners)/
  (gER + grec)

ggsave(paste0("figures/Fourpanel.", prod, ".png"), g4panel, height = 6, width = 9)




