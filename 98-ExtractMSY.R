# Extract MSY values from CM

# Libraries
library(salmonMSE)
library(tidyverse)
library(ggplot2)
library(here)
library(knitr)


# Input Conditioning Model results
ERM_QC <- readRDS("CM/QuinsamCampbell_06.15.26.rds")
report_QC <- salmonMSE:::get_report(ERM_QC)

ERM_Adam <- readRDS("CM/Adam_06.22.26.prior.rds")
report_Adam <- salmonMSE:::get_report(ERM_Adam)

ERM_Salmon <- readRDS("CM/Salmon_06.22.26.prior.rds")
report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

ERM_Woss <- readRDS("CM/Woss_06.22.26.prior.rds")
report_Woss <- salmonMSE:::get_report(ERM_Woss)

# Set up population
# pop <- "Adam"#"Woss"#"Salmon"#"Adam"#QC"

for (pop in c("Adam", "Salmon", "Woss")){
  samp <- get(paste0("ERM_", pop))
  d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
  if(pop == "QC") year1 <- 1984
  if(pop == "Nimpkish") year1 <- 2001
  if(pop == "Salmon") year1 <- 2002
  if(pop == "Adam") year1 <- 2002

  report <- get(paste0("report_", pop))

  year <- year1 + seq(1, d$Ldyr) - 1
  year_borrow <- seq(max(year) - 9, max(year) - 5)

  write_figures <- FALSE #Toggle to true to write figures to file

  theme_set(theme_bw())

  #------------------------------------------------------------------------------
  # MSY calculated with Ricker (lambert) equations

  alpha <- sapply(report, getElement, "alpha")# for egg-smolt rel
  beta <- sapply(report, getElement, "beta")# for egg-smolt rel
  # hist(alpha)

  # Time-varying reference points
  Sgen_s <- salmonMSE:::.CM_MSY(report, d, simple = TRUE, type = "Sgen")
  SMSY_s <- salmonMSE:::.CM_MSY(report, d, simple = TRUE, type = "spawner")
  UMSY_s <- salmonMSE:::.CM_MSY(report, d, simple = TRUE, type = "u")

  # MSY reference points where natural mortality and maturity are averaged across a range of years
  Sgen_avg <- salmonMSE:::.CM_MSY(report, d, simple = TRUE, index = match(year1:2024, year), mean_bio = TRUE, type = "Sgen")
  SMSY_avg <- salmonMSE:::.CM_MSY(report, d, simple = TRUE, index = match(year1:2024, year), mean_bio = TRUE, type = "spawner")



  # Get quantilies of time-varying SMSY
  SMSY_s_q <-  apply(SMSY_s, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = FALSE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1")) %>%
    rename(median='50%', lower='2.5%', upper='97.5%') %>%
    mutate(label="SMSY")

  # Plot 1: All SMSY values, including negative values
  g1 <- SMSY_s_q %>%
    ggplot(aes(Year, .data$median, colour= label, fill=label)) +
    geom_line() +
    scale_fill_manual(values = c("SMSY" = "chartreuse4")) +
    scale_colour_manual(values = c("SMSY" = "chartreuse4")) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color= NA) +
    labs(x = "Calendar Year", y = ylab)+
    theme(legend.title = element_blank()) +
    ylab("Spawners")

  # ggsave(paste("figures/", pop, "SMSY_calc_v1.png", sep=""), g1, height = 3.5, width = 6)

  # Get spawner time-series and bind to SMSY data frame to plot SMSY with spawners

  # Create df for spawners

  spawners <- sapply(report, getElement, "spawners") %>%
    apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
    t() %>%
    as.data.frame()

  Spawners_q <- data.frame(Year = year,
                           lower = spawners$'2.5%',#rep(NA,length(year)),
                           median = spawners$'50%', #d$obsescape,
                           upper = spawners$'97.5%',#rep(NA, length(year)),
                           label = "Spawners")
  SMSY_s_q <- rbind(SMSY_s_q, Spawners_q)

  # Plot 2: Add spawner times series to plot 1
  g2 <- SMSY_s_q %>%
    ggplot(aes(Year, .data$median, colour= label, fill=label)) +
    geom_line() +
    scale_fill_manual(values = c("SMSY" = "chartreuse4", "Spawners" = NA)) +
    scale_colour_manual(values = c("SMSY" = "chartreuse4", "Spawners" = "black")) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color= NA) +
    labs(x = "Calendar Year", y = ylab)+
    coord_cartesian(ylim = c(-30000,50000)) + #c(-3000,3000))+#c(-8000,14000)) + #c(-50000,50000)) + #c(-6000,4500)) + #c(-2500,3000))+#
    theme(legend.title = element_blank()) +
    ylab("Spawners")

  # ggsave(paste("figures/", pop, "SMSY_calc_v2.png", sep=""), g2, height = 3.5, width = 6)

  # Plot 3: Remove negative SMSY values
  na.rm <- TRUE
  if (na.rm) SMSY_s[SMSY_s <= 0] <- 0#NA_real_

  SMSY_s_q <- apply(SMSY_s, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = na.rm) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1")) %>%
    rename(median='50%', lower='2.5%', upper='97.5%') %>%
    mutate(label="SMSY")

  SMSY_s_q <- rbind(SMSY_s_q, Spawners_q)

  g3 <- SMSY_s_q %>%
    ggplot(aes(Year, .data$median, colour= label, fill=label)) +
    geom_line() +
    scale_fill_manual(values = c("SMSY" = "chartreuse4", "Spawners" = NA)) +
    scale_colour_manual(values = c("SMSY" = "chartreuse4", "Spawners" = "black")) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color= NA) +
    labs(x = "Calendar Year", y = ylab)+
    coord_cartesian(ylim = c(0, 15000)) + #c(0,3500)) + #c(0, 100000))  + #c(0,3500)) + #c(0, 50000)) +#c(0, 3500)) +# c(0,4500)) +
    theme(legend.title = element_blank()) +
    ylab("Spawners")

  # ggsave(paste("figures/", pop, "SMSY_calc_v3.png", sep=""), g3, height = 3.5, width = 6)

  # Plot 4: Add Sgen, and removing negative values from calculated Sgen values
  na.rm <- TRUE
  if (na.rm) Sgen_s[Sgen_s <= 0] <- 0#NA_real_

  Sgen_s_q <- apply(Sgen_s, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1")) %>%
    rename(median='50%', lower='2.5%', upper='97.5%')
  Sgen_s_q$label <- "Sgen"

  Sgen_s_q <- rbind(SMSY_s_q, Sgen_s_q, Spawners_q)

  g4 <- Sgen_s_q %>%
    ggplot(aes(Year, .data$median, colour= label, fill=label)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color= NA) +
    scale_fill_manual(values = c("Sgen" = "darkorange", "SMSY" = "chartreuse4", "Spawners" = NA)) +
    scale_colour_manual(values = c("Sgen" = "darkorange", "SMSY" = "chartreuse4", "Spawners" = "black")) +
    labs(x = "Calendar Year", y = ylab)+
    coord_cartesian(ylim = c(0,8000)) +#c(0,8000)) +#c(0,3000)) +# c(0,800)) + # c(0, 6000)) +#c(0,100000)) + #c(0,6000)) + #c(0, 50000)) + #c(0, 3500)) + #c(0,4500)) +
    ylab("Spawners") +
    theme(legend.title = element_blank())

  # ggsave(paste("figures/", pop, "SMSY_calc_v4.png", sep=""), g4, height = 3.5, width = 6)

  # Get average SMSY and Sgen values
  medianSMSY <- median(SMSY_avg)#SMSY_s_q %>% filter(label== "SMSY") %>% summarize(mean=mean(median)) %>% pull(mean)
  lwrSMSY <- quantile(SMSY_avg, prob=0.05)
  if(lwrSMSY < 0) lwrSMSY <- 0
  uprSMSY <- quantile(SMSY_avg, prob=0.95)
  medianSgen <- median(Sgen_avg)#Sgen_s_q %>% filter(label== "Sgen") %>% summarize(mean=mean(median)) %>% pull(mean)
  lwrSgen <- quantile(Sgen_avg, prob=0.05)
  if(lwrSgen < 0) lwrSgen <- 0
  uprSgen <- quantile(Sgen_avg, prob=0.95)
  cat("median SMSY", medianSMSY)
  cat("median Sgen", medianSgen)

  assign(paste0("medSMSY_", pop), medianSMSY)
  assign(paste0("medSgen_", pop), medianSgen)
  assign(paste0("lwrSMSY_", pop), lwrSMSY)
  assign(paste0("uprSMSY_", pop), uprSMSY)
  assign(paste0("lwrSgen_", pop), lwrSgen)
  assign(paste0("uprSgen_", pop), uprSgen)

  # meanSMSY <- SMSY_s_q %>% filter(label== "SMSY") %>% summarize(mean=mean(median)) %>% pull(mean)
  # meanSgen <- Sgen_s_q %>% filter(label== "Sgen") %>% summarize(mean=mean(median)) %>% pull(mean)

  # Plot 5: Add average SMSY and Sgen values
  g5 <- g4 +
    geom_hline(yintercept = medianSMSY, lty="dashed", colour = "chartreuse4") +
    geom_hline(yintercept = medianSgen, lty="dashed", colour = "darkorange")


  # Plot 6: Add uncertainty intervals in average SMSY and Sgen values, and dotted
  # lines for time-varying values. Include numerical averages values to the plot.
  # This plot is used for presentation

  na.rm <- TRUE
  if (na.rm) Sgen_s[Sgen_s <= 0] <- 0#NA_real_
  if (na.rm) SMSY_s[SMSY_s <= 0] <- 0#NA_real_

  SMSY_s_q <- apply(SMSY_s, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1")) %>%
    rename(median='50%', lower='2.5%', upper='97.5%') %>%
    mutate(lower = rep(NA, length(Year)), upper = rep(NA, length(Year))) %>%
    mutate(label="SMSY")#
  SMSYAve_s_q <- data.frame(Year = year, median = median(SMSY_avg),
                            lower = as.numeric(apply(SMSY_avg,1,
                                                     quantile,probs=c( 0.025),
                                                     na.rm=T)),
                            upper = as.numeric(apply(SMSY_avg,1,
                                                     quantile,probs=c( 0.975),
                                                     na.rm=T)),
                            label="SMSY-ave")


  SMSY_s_q <- rbind(SMSY_s_q, SMSYAve_s_q)

  Sgen_s_q <- apply(Sgen_s, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1")) %>%
    rename(median='50%', lower='2.5%', upper='97.5%') %>%
    mutate(lower = rep(NA, length(Year)), upper = rep(NA, length(Year))) %>%
    mutate(label="Sgen")#
  SgenAve_s_q <- data.frame(Year = year, median = median(Sgen_avg),
                            lower = apply(Sgen_avg,1,
                                          quantile,probs=c( 0.025),
                                          na.rm=T),
                            upper = apply(Sgen_avg,1,
                                          quantile,probs=c( 0.975),
                                          na.rm=T),
                            label = "Sgen-ave")
  Sgen_s_q <- rbind(Sgen_s_q, SgenAve_s_q)

  # Create df for spawners

  SMSY_s_q <- rbind(SMSY_s_q, Spawners_q)

  SMSY_s_q <- rbind(SMSY_s_q, Sgen_s_q, Spawners_q)
  SMSY_s_q <- SMSY_s_q %>%
    mutate(lower = ifelse(lower < 0,0,lower))



  g6 <- SMSY_s_q %>%
    ggplot(aes(Year, .data$median, colour= label, fill=label, group=label)) +
    geom_line(aes(linetype = label)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color= NA,
                show.legend = FALSE) +
    scale_fill_manual(values = c("Sgen-ave" = "darkorange",
                                 "SMSY-ave" = "chartreuse4", "Spawners" = "black")) +
    scale_colour_manual(values = c("Sgen" = "darkorange",
                                   "Sgen-ave" = "darkorange",
                                   "SMSY" = "chartreuse4",
                                   "SMSY-ave" = "chartreuse4",
                                   "Spawners" = "black")) +
    scale_linetype_manual(values = c("dashed", "solid", "dashed", "solid", "solid")) +
    labs(x = "Calendar Year", y = ylab) +
    coord_cartesian(ylim =
                      c(0,max(c(Spawners_q$median, SMSY_s_q$median, SMSY_s_q$upper),
                                na.rm=T))) +#c(0,8000)) +#c(0,3000)) +# c(0,800))
    ylab("Spawners") +
    theme(legend.title = element_blank())+
    annotate(geom="text", x=2023, y=median(SMSY_avg)+median(SMSY_avg)/3,
             label=c(round (median(SMSY_avg))), size=2, colour = "chartreuse4") +
    annotate(geom="text", x=2023, y=median(Sgen_avg)+median(Sgen_avg)/3,
             label=c(round (median(Sgen_avg))), size=2, colour = "darkorange")

  g6
  if(write_figures) ggsave(paste("figures/", pop, "QC_SMSY_calc_v6.png", sep=""), g6, height = 3.5, width = 6)

  #---------------------------------------------------------------------------
  # UMSY plots

  #Set negative UMSY to zero (where prod<1)
  UMSY_s[UMSY_s <= 0] <- 0#NA_real_

  UMSY_avg <- salmonMSE:::.CM_MSY(report, d, simple = TRUE,
                                  mean_bio = TRUE, type = "u")
  quantile(UMSY_avg, c(0.05, 0.25, 0.5, 0.75, 0.95))


  UMSY_s_q <-  apply(UMSY_s, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1")) %>%
    rename(median='50%', lower='2.5%', upper='97.5%') %>%
    mutate(lower = rep(NA, length(Year)), upper = rep(NA, length(Year))) %>%
    mutate(label="UMSY")

  UMSYAve_s_q <- data.frame(Year = year, median = median(UMSY_avg),
                            lower = apply(UMSY_avg,1,
                                          quantile,probs=c( 0.025),
                                          na.rm=T),
                            upper = apply(UMSY_avg,1,
                                          quantile,probs=c( 0.975),
                                          na.rm=T),
                            label = "UMSY-ave")
  UMSY_s_q <- rbind(UMSY_s_q, UMSYAve_s_q)
  UMSY_s_q <- UMSY_s_q %>%
    mutate(lower = ifelse(lower < 0,0,lower))


  # Get ERS (see also file 98-ExtractERs.R )
  df.all <- salmonMSE:::.CM_ER(report, type = "all", r =1 , brood = FALSE,
                               index_AEQ = match(year_borrow, year))
  ER.all <- apply(df.all, 1, function(row)
    quantile(row, probs = c(0.025, 0.5, 0.975), na.rm=T))
  ER <- as.data.frame(t(ER.all))
  ER <- ER %>% rename("lower"="2.5%", "median"="50%" ,"upper"="97.5%") %>%
    mutate (label="ER", Year = seq(from = year1, length.out = length(ER[,1])))
  ER$label <- "ER (AEQ)"

  gUMSY <- rbind(UMSY_s_q, ER) %>%
    ggplot(aes(Year, .data$median, colour= label, fill=label, group=label)) +
    geom_line(aes(linetype = label)) +
    scale_fill_manual(values = c("UMSY" = "black", "UMSY-ave" = "black",
                                 "ER (AEQ)" = "maroon")) +
    scale_colour_manual(values = c("UMSY" = "black", "UMSY-ave" = "black",
                                   "ER (AEQ)" = "maroon")) +
    scale_linetype_manual(values = c( "solid", "dashed","solid")) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color= NA,
                show.legend = FALSE) +
    labs(x = "Calendar Year", y = ylab)+
    theme(legend.title = element_blank()) +
    ylab("UMSY") +
    annotate(geom="text", x=2023, y=median(UMSY_avg)+0.02,
             label=c(round (median(UMSY_avg),2)), size=2) +
    ylim(0,1)

  gUMSY

  if(write_figures) ggsave(paste("figures/", pop, "QC_UMSY_calc.png", sep=""), gUMSY, height = 3.5, width = 6)


}


# -----------------------------------------------------------------------------
# Sum benchmarks
# Need to assign medSMSY and medSgen values above first

  meduBench_Nat <- 0.85* (medSMSY_Adam + medSMSY_Salmon + medSMSY_Woss)
  medlBench_Nat <- (medSgen_Adam + medSgen_Salmon + medSgen_Woss)

  df_bench_SR <- data.frame(Pop = c(rep("Adam",2),
                                 rep("Salmon", 2),
                                 rep("Nimpkish", 2)),
                         Benchmark = rep(c( "Sgen", "85% SMSY"), 3),
                         med = c(round(medSgen_Adam,0),
                                 round(0.85 * medSMSY_Adam),
                                 round(medSgen_Salmon),
                                 round(0.85 *medSMSY_Salmon),
                                 round(medSgen_Woss),
                                 round(0.85 * medSMSY_Woss)),
                         lwr = c(round(lwrSgen_Adam),
                                 round(0.85 * lwrSMSY_Adam),
                                 round(lwrSgen_Salmon),
                                 round(0.85 * lwrSMSY_Salmon),
                                 round(lwrSgen_Woss),
                                 round(0.85 * lwrSMSY_Woss)),
                         upr = c(round(uprSgen_Adam),
                                 round(0.85 * uprSMSY_Adam),
                                 round(uprSgen_Salmon),
                                 round(0.85 * uprSMSY_Salmon),
                                 round(uprSgen_Woss),
                                 round(0.85 * uprSMSY_Woss))
                         )



# Get Habitat based estimates of Smax
Smax_hab <- as.data.frame( read.csv(
  here("data", "UpperSoGChinook_out_posteriorpredictive.csv"))
)

# Smax_hab[ Smax_hab[,'Stock'] =="Adam/Eve",] $ SMAX_median * 0.2
# Smax_hab %>% filter(Stock== "Adam/Eve") %>% pull(SMAX_median) * 0.2

df_bench_hab <- data.frame(Pop = c(rep("Adam",2),
                               rep("Salmon", 2),
                               rep("Nimpkish", 2)),
                       Benchmark = rep(c("20% SMAX", "40% SMAX"), 3),
                       med = c(Smax_hab %>% filter(Stock== "Adam/Eve") %>% pull(SMAX_median) * 0.2,
                               Smax_hab %>% filter(Stock== "Adam/Eve") %>% pull(SMAX_median) * 0.4,
                               Smax_hab %>% filter(Stock== "Salmon") %>% pull(SMAX_median) * 0.2,
                               Smax_hab %>% filter(Stock== "Salmon") %>% pull(SMAX_median) * 0.4,
                               Smax_hab %>% filter(Stock== "Nimpkish") %>% pull(SMAX_median) * 0.2,
                               Smax_hab %>% filter(Stock== "Nimpkish") %>% pull(SMAX_median) * 0.4),
                       lwr = c(Smax_hab %>% filter(Stock== "Adam/Eve") %>% pull(SMAX_lwr5) * 0.2,
                               Smax_hab %>% filter(Stock== "Adam/Eve") %>% pull(SMAX_lwr5) * 0.4,
                               Smax_hab %>% filter(Stock== "Salmon") %>% pull(SMAX_lwr5) * 0.2,
                               Smax_hab %>% filter(Stock== "Salmon") %>% pull(SMAX_lwr5) * 0.4,
                               Smax_hab %>% filter(Stock== "Nimpkish") %>% pull(SMAX_lwr5) * 0.2,
                               Smax_hab %>% filter(Stock== "Nimpkish") %>% pull(SMAX_lwr5) * 0.4),
                       upr = c(Smax_hab %>% filter(Stock== "Adam/Eve") %>% pull(SMAX_upr95) * 0.2,
                               Smax_hab %>% filter(Stock== "Adam/Eve") %>% pull(SMAX_upr95) * 0.4,
                               Smax_hab %>% filter(Stock== "Salmon") %>% pull(SMAX_upr95) * 0.2,
                               Smax_hab %>% filter(Stock== "Salmon") %>% pull(SMAX_upr95) * 0.4,
                               Smax_hab %>% filter(Stock== "Nimpkish") %>% pull(SMAX_upr95) * 0.2,
                               Smax_hab %>% filter(Stock== "Nimpkish") %>% pull(SMAX_upr95) * 0.4))
df_bench_hab$med <- round(df_bench_hab$med)
df_bench_hab$lwr <- round(df_bench_hab$lwr)
df_bench_hab$upr <- round(df_bench_hab$upr)

order.pop <- c("Adam", "Salmon", "Nimpkish")
order.bench <- c("Sgen", "85% SMSY", "20% SMAX", "40%SMAX")

df_bench <- rbind(df_bench_SR, df_bench_hab) %>%
  arrange(match(Pop, order.pop), match(Benchmark, order.bench))
kable(df_bench)

if(write_figures) {
  write.csv(df_bench, here("data", "Benchmarks_NatPops.csv"),
            row.names = FALSE)
}# read.csv(file = here("data", "Benchmarks_NatPops.csv"))









#------------------------------------------------------------------------------
# Not implemented
# Get numerical estimates of MSY values

if (FALSE){
  SMSY <- salmonMSE:::.CM_MSY(report, d, type = "spawner", ncores = 10) # year x simulation
  #EMSY <- salmonMSE:::.CM_MSY(report, d, type = "egg", ncores = 10)
  UMSY <- salmonMSE:::.CM_MSY(report, d, type = "u", ncores = 10)
  Sgen <- salmonMSE:::.CM_Sgen(report, d, ncores = 10)

  # gSMSY <- salmonMSE:::CM_MSY(report, d, type = "spawner", ncores = 7, year1=2010)
  # gSgen <- salmonMSE:::CM_Sgen(report, d, ncores = 7 , year1=2010)
  gUMSY <- salmonMSE:::CM_MSY(report, d, type = "u", ncores = 7 , year1=2010)
  ggsave(paste("figures/", pop, "UMSY.png", sep=""), gUMSY, height = 3.5, width = 6)

  #------------------------------------------------------------------------------
  # SMSY estimated numerically: Plots
  MSY <- SMSY
  type <- "spawner"
  na.rm <- FALSE
  if (na.rm) MSY[MSY <= 0] <- NA_real_

  MSY_q <- apply(MSY, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = na.rm) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1"))


  MSY_q$label <- "SMSY"
  MSY_q <- MSY_q %>% rename(median='50%', lower='2.5%', upper='97.5%')

  g1 <- MSY_q %>%
    ggplot(aes(Year, .data$median, colour= label, fill=label)) +
    geom_line() +
    scale_fill_manual(values = c("SMSY" = "chartreuse4")) +
    scale_colour_manual(values = c("SMSY" = "chartreuse4")) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color= NA) +
    labs(x = "Calendar Year", y = ylab)+
    theme(legend.title = element_blank()) +
    ylab("Spawners")

  g1
  ggsave(paste("figures/", pop, "SMSY_v1.png", sep=""), g1, height = 3.5, width = 6)

  Spawners_q <- data.frame(Year=MSY_q$Year,
                           lower=rep(NA,length(MSY_q$Year)),
                           median=d$obsescape,
                           upper=rep(NA, length(MSY_q$Year)),
                           label="Spawners")
  MSY_q <- rbind(MSY_q, Spawners_q)

  g2 <- MSY_q %>%
    ggplot(aes(Year, .data$median, colour= label, fill=label)) +
    geom_line() +
    scale_fill_manual(values = c("SMSY" = "chartreuse4", "Spawners" = NA)) +
    scale_colour_manual(values = c("SMSY" = "chartreuse4", "Spawners" = "black")) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color= NA) +
    labs(x = "Calendar Year", y = ylab)+
    coord_cartesian(ylim = c(-20000, 20000)) + # c(-6000,4500)) +
    theme(legend.title = element_blank()) +
    ylab("Spawners")

  g2
  ggsave(paste("figures/", pop, "SMSY_v2.png", sep=""), g2, height = 3.5, width = 6)

  # SMSY: removing negative values
  na.rm <- TRUE
  if (na.rm) MSY[MSY <= 0] <- NA_real_

  MSY_q <- apply(MSY, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = na.rm) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1")) %>%
    rename(median='50%', lower='2.5%', upper='97.5%') %>%
    mutate(label="SMSY")


  MSY_q <- rbind(MSY_q, Spawners_q)

  g3 <- MSY_q %>%
    ggplot(aes(Year, .data$median, colour= label, fill=label)) +
    geom_line() +
    scale_fill_manual(values = c("SMSY" = "chartreuse4", "Spawners" = NA)) +
    scale_colour_manual(values = c("SMSY" = "chartreuse4", "Spawners" = "black")) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color= NA) +
    labs(x = "Calendar Year", y = ylab)+
    coord_cartesian(ylim = c(0,22000)) + #c(0,4500)) +
    theme(legend.title = element_blank()) +
    ylab("Spawners")

  g3
  ggsave(paste("figures/", pop, "SMSY_v3.png", sep=""), g3, height = 3.5, width = 6)


  # Adding Sgen
  na.rm <- TRUE
  Sgen_q <- apply(Sgen, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1")) %>%
    rename(median='50%', lower='2.5%', upper='97.5%')
  Sgen_q$label <- "Sgen"

  Sgen_q <- rbind(MSY_q, Sgen_q, Spawners_q)

  g4 <- Sgen_q %>%
    ggplot(aes(Year, .data$median, colour= label, fill=label)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color= NA) +
    scale_fill_manual(values = c("Sgen" = "darkorange", "SMSY" = "chartreuse4", "Spawners" = NA)) +
    scale_colour_manual(values = c("Sgen" = "darkorange", "SMSY" = "chartreuse4", "Spawners" = "black")) +
    labs(x = "Calendar Year", y = ylab)+
    coord_cartesian(ylim = c(0,6000)) + #c(0, 100000)) + #c(0, 22000)) + #c(0,4500)) +
    ylab("Spawners") +
    theme(legend.title = element_blank())

  g4
  ggsave(paste("figures/", pop, "SMSY_v4.png", sep=""), g4, height = 3.5, width = 6)


  meanSMSY <- MSY_q %>% filter(label== "SMSY") %>% summarize(mean=mean(median)) %>% pull(mean)
  meanSgen <- Sgen_q %>% filter(label== "Sgen") %>% summarize(mean=mean(median)) %>% pull(mean)


  g5 <- g4 +
    geom_hline(yintercept = meanSMSY, lty="dashed", colour = "chartreuse4") +
    geom_hline(yintercept = meanSgen, lty="dashed", colour = "darkorange")
  g5

  ggsave(paste("figures/", pop, "SMSY_v5.png", sep=""), g5, height = 3.5, width = 6)




  EMSY_q <- apply(EMSY, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = FALSE) %>%
    reshape2::melt() %>%
    mutate(Year = Var2 + year1 - 1) %>%
    reshape2::dcast(list("Year", "Var1")) %>%
    rename(median='50%', lower='2.5%', upper='97.5%') %>%
    mutate(label = "EMSY")
  mean(EMSY_q$median)



  # Replace negative SMSY values with NA
  SMSY_good <- SMSY
  SMSY_good[SMSY_good <= 0] <- NA_real_
  year <- 2010:2024
  median_posSMSY <- apply(SMSY_good, 1, median, na.rm = TRUE)
  plot(year, median_posSMSY, type = 'o')
  prop_bad <- apply(SMSY_good, 1, function(x) mean(is.na(x))) # Annual proportion of MCMC samples with SMSY < 0
  plot(year, prop_bad) # Better-behaved SMSY values once we exclude the most recent brood years

  mean(median_SMSY[1:11], na.rm = T)
  mean(median_posSMSY[1:11], na.rm = T)

  median_allSMSY <- apply(SMSY,1,median, na.rm = TRUE)
  plot(year, median_allSMSY, type='o')
  abline(h=0)

  median_UMSY <- apply(UMSY, 1, median, na.rm = TRUE)
  plot(year, median_UMSY, type='o')


}


# Additional Plots Not implemented
if (FALSE) {
  # Figure of SMSY
  # SMSY has some very skewed values for some stocks
  g <- CM_MSY(report, d, simple = TRUE, index = match(year1:2024, year), mean_bio = TRUE, type = "spawner")
  g
  g <- CM_MSY(report, d, simple = TRUE, index = match(year1:2024, year), mean_bio = TRUE, type = "u")
  g

  # Figure of time-varying productiity
  CM_prod(report, d, year1 = year1, index = match(year1:2024, year), mean_bio = TRUE)
  CM_prod(report, d, year1 = year1, index = match(year1:2024, year), mean_bio = FALSE)

  # Vulnerability and maturity in a single year, histograms of numerical vs. Lambert SMSY and UMSY
  par(mfrow = c(3, 2))

  vulPT <- sapply(report, getElement, "vulPT")
  matplot(vulPT, type = 'l', xlab = "Age", ylab = "Preterminal vulnerability")

  matt <- sapply(report, getElement, "matt", simplify = "array")
  y <- 10 # Select year
  matplot(matt[y, , , ], type = "l", xlab = "Age", ylab = "Proportion mature",
          main = paste("Year", year1 + y - 1))

  SMSY_y <- SMSY[y, ]
  UMSY_y <- UMSY[y, ]
  hist(SMSY_y[SMSY_y < 1e5], main = paste("Median =", median(SMSY_y) |> round()),
       xlab = paste("Numerical SMSY", year1 + y - 1))
  median(SMSY[y, ])
  hist(SMSY_s[SMSY_s < 1e5 & SMSY_s > 0],
       main = paste("Median =", median(SMSY_s[SMSY_s > 0], na.rm = TRUE) |> round()),
       xlab = "Lambert SMSY")

  hist(UMSY[y, ], main = paste("Median =", median(UMSY_y) |> round(2)),
       xlab = paste("Numerical UMSY", year1 + y - 1))
  median(UMSY[y, ])
  hist(UMSY_s[UMSY_s > 0],
       main = paste("Median =", median(UMSY_s[UMSY_s], na.rm = TRUE) |> round(2)),
       xlab = "Lambert UMSY")
}

#------------------------------------------------------------------------------
# Not implemented
# Compare reference points- numerical and estimated

if(FALSE){
  require(salmonMSE)
  compare_ref <- function(a = 3, # Units of recruits/spawner
                          b = 1/1000, # Units of reciprocal spawners
                          maxage = 5,
                          rel_F = c(0, 1), # e_1 and e_2
                          p_female = 1,
                          vulPT = rep(0, maxage),
                          vulT = c(0, 0, 0.2, 0.5, 1),
                          M = c(1, 0.8, 0.6, 0.4),
                          fec = c(0, 1500, 3000, 3200, 3500),
                          p_mature = c(0, 0.1, 0.2, 0.3, 1)) {


    # Calculate eggs per smolt
    phi <- salmonMSE:::calc_phi(
      Mjuv = M,
      p_mature = p_mature,
      p_female = p_female,
      fec = fec,
      s_enroute = 1
    )

    # To convert Smax from units of spawners to eggs:
    # 1. Calculate unfished spawners from Ricker SRR, set R = S --> S0 = log(a)/b
    # 2. Calculate spawners per juvenile (spro)
    # 3. Calculate smolt per egg (phi)
    # 4. unfished eggs (eo) = unfished spawners * juvenile per spawner * egg per smolt
    # 5. Emax = log(a) / unfished eggs
    spro <- salmonMSE:::calc_phi(
      Mjuv = M,
      p_mature = p_mature,
      p_female = p_female,
      fec = fec,
      s_enroute = 1,
      output = "spawner"
    )

    Smax <- 1/b
    so <- Smax * log(a)
    eo <- so / spro * phi
    beta_eggs <- log(a)/eo
    Emax_eggs <- 1/beta_eggs

    SRRpars <- data.frame(
      kappa = a,
      Smax = Smax,
      Emax = Emax_eggs,
      tau = spro,
      phi = phi,
      SRrel = "Ricker"
    )

    # salmonMSE ref
    ref <- calc_MSY(
      M, fec, p_female, rel_F, vulPT, vulT, p_mature,
      s_enroute = 1, SRRpars = SRRpars
    )

    ref["Sgen"] <- calc_Sgen(
      M, fec, p_female, rel_F, vulPT, vulT, p_mature,
      s_enroute = 1, SRRpars = SRRpars, SMSY = ref["Spawners_MSY"]
    )
    ref["Catch/Return"] <- ref["KT_MSY"]/ref["Return_MSY"]

    ref_salmonMSE <- structure(
      ref[c("UPT_MSY", "UT_MSY", "Catch/Return",  "Spawners_MSY", "Sgen")],
      names = c("UMSY (preterminal)", "UMSY (terminal)", "Terminal Catch/Return",  "SMSY", "Sgen")
    )

    # Ricker SRR calcs
    ref_Ricker <- local({
      umsy <- calc_Umsy_Ricker(log(a))
      smsy <- calc_Smsy_Ricker(log(a), b)
      sgen <- calc_Sgen_Ricker(log(a), b)
      structure(c(umsy, smsy, sgen), names = c("UMSY", "SMSY", "Sgen"))
    })

    output <- list(
      salmonMSE = ref_salmonMSE,
      Ricker = ref_Ricker
    )
    return(output)
  }

  compare_ref(
    a = 3,
    b = 1/1000,
    maxage = 5,
    rel_F = c(1, 0),
    p_female = 1,
    vulPT = c(0, 0.1, 0.2, 0.4, 1),
    vulT = c(0, 0, 0, 0, 0),
    M = c(1, 0.3, 0.2, 0.1),
    fec = rep(1, 5),
    # fec = c(0, 1000, 2000, 3000, 3500),
    p_mature = c(0, 0.1, 0.2, 0.3, 1)
  )

}

