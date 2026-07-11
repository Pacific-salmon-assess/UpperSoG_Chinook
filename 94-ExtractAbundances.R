# Extract Abundances: Catches, Spawners, Recruitment

# Libraries
library(patchwork)
library(ggplot2)
library(tidyverse)
library(salmonMSE)
library(here)


# Pull results from conditioning model
ERM_QC <- readRDS("CM/QuinsamCampbell_06.15.26.rds")
report_QC <- salmonMSE:::get_report(ERM_QC)

ERM_Adam <- readRDS("CM/Adam_06.22.26.prior.rds")#readRDS("CM/Adam_06.01.26.JSpt.rds")
report_Adam <- salmonMSE:::get_report(ERM_Adam)

ERM_Salmon <- readRDS("CM/Salmon_06.22.26.prior.rds") #file:///C:/github/UpperSoG_Chinook/CM/Salmon_06.10.JSesc.NGTSt.html
report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

ERM_Woss <- readRDS("CM/Woss_06.22.26.prior.rds")
report_Woss <- salmonMSE:::get_report(ERM_Woss)



folder_path <- "data/CMoutput"

if (!dir.exists(here::here(folder_path))) {
  dir.create(here::here(folder_path), recursive = TRUE)
}

# Make folder for CM results, and put csv file there...

# pop <- "QC"
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
                here::here(paste0("data/CMoutput/", pop, "_timeseries.csv")))

  df_withpop <- df %>% mutate(Pop = pop)
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
medspawners <- median(spawners$med)*1000
spawners_bp <- spawners %>% filter(y<2012 & year>2001) %>% select(med)
avespawners_bp <- mean(spawners_bp$med)*1000

empUSR <- data.frame(medspawners = medspawners, avespawners_bp = avespawners_bp)

write.csv(
  empUSR,
  here(
    "data",
    "EmpiricalUSR.csv"
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
         fill = guide_legend(reverse = TRUE))

# + theme(legend.position = "none")
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
  geom_hline(yintercept= medlBench_Nat/1000, colour= "darkred", linetype= "dashed") +
  geom_hline(yintercept= meduBench_Nat/1000, colour= "darkgreen", linetype= "dashed")


  # + theme(legend.position = "none")

gspawners

# Add 85%SMY and Sgen for natural dominated systems

# ER

UMSY_out <- read.csv( here(
  "data",
  "Equilibrium trade-off analysis",
  "R-OUT_SMU_ref-pt_values_eq-trade-off.csv"))
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
  geom_hline(yintercept = UMSY_med, colour = "darkgrey", linetype = "dashed")+
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

ggsave(paste("figures/Fourpanel.png", sep=""), g4panel, height = 6, width = 9)
