# Extract MSY values from CM

# Libraries
library(salmonMSE)
library(tidyverse)
library(ggplot2)
library(here)
library(knitr)
library(patchwork)


# Input Conditioning Model results
ERM_QC <- readRDS("CM/QuinsamCampbell_07.29.26.rds")
report_QC <- salmonMSE:::get_report(ERM_QC)

ERM_Adam <- readRDS("CM/Adam_08.06.26.prior.rds")
report_Adam <- salmonMSE:::get_report(ERM_Adam)

ERM_Salmon <- readRDS("CM/Salmon_08.06.26.prior.rds")
report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

ERM_Woss <- readRDS("CM/Woss_07.22.26.prior.rds")
report_Woss <- salmonMSE:::get_report(ERM_Woss)

# Set up population
# pop <- "Adam"#"Woss"#"Salmon"#"Adam"#QC"

for (pop in c("Adam", "Salmon", "Woss")){
  samp <- get(paste0("ERM_", pop))
  d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
  if(pop == "QC") year1 <- 1984
  if(pop == "Woss") year1 <- 2001
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


  #  SMAX values where where natural mortality and maturity are averaged across a range of years
  get_beta_s <- function(report, samp){
    d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
    alpha <- sapply(report, getElement, "alpha")# for egg-smolt rel
    beta <- sapply(report, getElement, "beta")# for egg-smolt rel
    alpha_s <- salmonMSE:::.CM_prod(report, d, mean_bio = TRUE) # Ricker alpha, per spawner
    epro <- t(alpha_s)/alpha # s, y

    spro <- sapply(1:length(report), function(x) { # vector
      mo <- apply(report[[x]]$mo[, , drop = FALSE], 2, mean)
      matt <- apply(report[[x]]$matt[, , d$r_matt, drop = FALSE], 2, mean)

      lo <- salmonMSE:::calc_survival(mo, matt) # smolt survival at replacement
      spro <- sum(lo * d$ssum * matt)
      return(spro)
    })

    beta_s <- beta * epro/spro # Ricker beta, per spawner
    Srep <- log(t(alpha_s))/beta_s
    return(as.vector(beta_s))
  }
  beta_avg <- get_beta_s(report, samp)
  SMAX_avg <- 1/beta_avg



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

  # Get average SMSY, Sgen and Smax values
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
  medianSMAX <- median(SMAX_avg)
  lwrSMAX <- quantile(SMAX_avg, prob=0.05)
  if(lwrSMAX < 0) lwrSMAX <- 0
  uprSMAX <- quantile(SMAX_avg, prob=0.95)


  assign(paste0("medSMSY_", pop), medianSMSY)
  assign(paste0("medSgen_", pop), medianSgen)
  assign(paste0("lwrSMSY_", pop), lwrSMSY)
  assign(paste0("uprSMSY_", pop), uprSMSY)
  assign(paste0("lwrSgen_", pop), lwrSgen)
  assign(paste0("uprSgen_", pop), uprSgen)

  assign(paste0("medSMAX_", pop), medianSMAX)
  assign(paste0("lwrSMAX_", pop), lwrSMAX)
  assign(paste0("uprSMAX_", pop), uprSMAX)
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
             label=c(round (median(Sgen_avg))), size=2, colour = "darkorange") +
    annotate( "text",  x = -Inf, y = Inf,
              label = ifelse(pop =="Adam", "f", "i"),
              hjust = -1, vjust = 1.8, size = 6, fontface = "bold")


  g6
  if(write_figures) ggsave(paste("figures/", pop, "_SMSY_calc_v6.png", sep=""), g6, height = 3.5, width = 6)
  assign(paste0("gSp_", pop), g6)


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
    ylab("Exploitation rate (AEQ)") +
    annotate(geom="text", x=2023, y=median(UMSY_avg)+0.02,
             label=c(round (median(UMSY_avg),2)), size=2) +
    ylim(0,1) +
    annotate( "text",  x = -Inf, y = Inf,
              label = ifelse(pop =="Adam", "g", "j"),
              hjust = -1, vjust = 1.8, size = 6, fontface = "bold")


  gUMSY

  if(write_figures) ggsave(paste("figures/", pop, "_UMSY_calc.png", sep=""), gUMSY, height = 3.5, width = 6)
  assign(paste0("gU_", pop), gUMSY)

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


  df_Smaxbench_SR <- data.frame(Pop = c(rep("Adam", 2),
                                        rep("Salmon", 2),
                                        rep("Nimpkish", 2)),
                            Benchmark = rep(c("20% SMAX","40% SMAX"), 3),
                            med = c(round(medSMAX_Adam) * 0.2,
                                    round(medSMAX_Adam) * 0.4,
                                    round(medSMAX_Salmon) * 0.2,
                                    round(medSMAX_Salmon) * 0.4,
                                    round(medSMAX_Woss) * 0.2,
                                    round(medSMAX_Woss) * 0.4),
                            lwr = c(round(lwrSMAX_Adam) * 0.2,
                                    round(lwrSMAX_Adam) * 0.4,
                                    round(lwrSMAX_Salmon) * 0.2,
                                    round(lwrSMAX_Salmon) * 0.4,
                                    round(lwrSMAX_Woss) * 0.2,
                                    round(lwrSMAX_Woss) * 0.4),
                            upr = c(round(uprSMAX_Adam) * 0.2,
                                    round(uprSMAX_Adam) * 0.4,
                                    round(uprSMAX_Salmon) * 0.2,
                                    round(uprSMAX_Salmon) * 0.4,
                                    round(uprSMAX_Woss) * 0.2,
                                    round(uprSMAX_Woss) * 0.4)
                            )

  order.pop <- c("Adam", "Salmon", "Nimpkish")
  order.bench <- c("Sgen", "85% SMSY", "20% SMAX", "40% SMAX")

  df_bench_full <- rbind(df_bench_SR, df_Smaxbench_SR) %>%
    arrange(match(Pop, order.pop), match(Benchmark, order.bench))
  kable(df_bench_full)

  df_bench_simple <- df_bench_SR %>%
    arrange(match(Pop, order.pop), match(Benchmark, order.bench))
  kable(df_bench_simple)

  if(write_figures) {
    write.csv(df_bench_full, here("data", "BenchmarksFull_NatPops.csv"),
              row.names = FALSE)
    write.csv(df_bench_simple, here("data", "BenchmarksSimple_NatPops.csv"),
              row.names = FALSE)
  }# read.csv(file = here("data", "Benchmarks_NatPops.csv"))



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
  write.csv(df_bench, here("data", "BenchmarksHab_NatPops.csv"),
            row.names = FALSE)
}# read.csv(file = here("data", "Benchmarks_NatPops.csv"))




#------------------------------------------------------------------------------
# Get additional outputs from CM for population assessment
pop <- "Adam"# "Woss", "Salmon" "Adam" "QC"

samp <- get(paste0("ERM_", pop))
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
if(pop == "QC") year1 <- 1984
if(pop == "Woss") year1 <- 2001
if(pop == "Salmon") year1 <- 2002
if(pop == "Adam") year1 <- 2002
report <- get(paste0("report_", pop))
brood <- FALSE
ci <-  TRUE
r <-  1
year <- year1 + seq(1, d$Ldyr) - 1
year_borrow <- seq(max(year) - 9, max(year) - 5)
theme_set(theme_bw() +
            theme( panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
)

# SRR plot
# gSRR <- CM_SRR(report, year1)

egg <- sapply(report, getElement, "egg") %>% apply(1, median)
smolt <- sapply(report, function(x) x$N[-1, 1, 1]) %>% apply(1, median)
epred <- seq(0, 1.1 * max(egg), length.out = 50)
spred <- sapply(report, function(x) x$alpha * epred * exp(-x$beta * epred)) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

year <- year1 + seq(1, length(egg)) - 1

df <- data.frame(
  year = year,
  egg = egg,
  smolt = smolt
)

df_med <- data.frame(
  egg = epred,
  smolt = spred[2, ]
)

df_poly1 <- data.frame(
  egg = epred,
  smolt = spred[1, ]
)

df_poly2 <- data.frame(
  egg = rev(epred),
  smolt = rev(spred[3, ])
)

if(pop == "QC") {
  gSRR <- ggplot(df, aes(.data$egg, .data$smolt)) +
    geom_point(shape = 1) +
    # geom_line(data = df_med) +
    # geom_polygon(data = rbind(df_poly1, df_poly2), fill = "grey", alpha = 0.5) +
    labs(x = "Egg production", y = "Smolt production") +
    expand_limits(x = 0, y = 0) +
    geom_hline(yintercept = median(smolt), linetype = "dashed") +
    annotate( "text",  x = -Inf, y = Inf, label = "a",  hjust = -1,
              vjust = 1.8, size = 6, fontface = "bold")


}

if(pop != "QC"){
  gSRR <- ggplot(df, aes(.data$egg, .data$smolt)) +
    geom_point(shape = 1) +
    geom_line(data = df_med) +
    geom_polygon(data = rbind(df_poly1, df_poly2), fill = "grey", alpha = 0.5) +
    labs(x = "Egg production", y = "Smolt production") +
    expand_limits(x = 0, y = 0) +
    annotate( "text",  x = -Inf, y = Inf, label = "a",  hjust = -0.3,
              vjust = 1.3, size = 6, fontface = "bold")

}

if (requireNamespace("ggrepel", quietly = TRUE)) {
  gSRR <- gSRR + ggrepel::geom_text_repel(aes(label = .data$year))
}



# Marine survival plot to age 2
# gSurv <- CM_surv2(report, year1)

NO <- exp(-sapply(report, function(x) x$mo[, 1])) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(type = "Natural", Year = year1 + 1:n() - 1)

HO <- exp(-sapply(report, getElement, "moplot")) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(type = "Hatchery", Year = year1 + 1:n() - 1)

df <- rbind(NO, HO) %>%
  mutate(BroodYear = .data$Year - 1)

gSurv <- ggplot(df, aes(.data$BroodYear, .data$`50%`)) +
  geom_line(aes(colour = .data$type)) +
  labs(x = "Brood Year", y = "Marine survival to age 2", fill = "Origin", colour = "Origin") +
  expand_limits(y = 0) +
  theme(legend.position = "right") +
  annotate( "text",  x = -Inf, y = Inf, label = "b",  hjust = -1, vjust = 1.8,
            size = 6, fontface = "bold")

gSurv <- gSurv + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`, fill = .data$type), alpha = 0.25)


# Maturity plot
r <- d$r_matt
brood <- TRUE
rs_names <- "smolt0+"
n_r <- d$n_r
if (missing(rs_names)) rs_names <- seq(1, n_r)

# bmatt <- data.frame(Age = 1:d$Nages, value = d$bmatt) %>%
#   dplyr::filter(Age > 1)
matt <- sapply(report, function(i)
  salmonMSE:::CY2BY(i[["matt"]][, , r]), simplify = 'array')
matt_q <- apply(matt, 1:2, quantile,
                probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 + year1 - 1) %>%
  rename(Age = Var3) %>%
  reshape2::dcast(Age + Var2 + Year ~ Var1) %>%
  filter(!is.na(`50%`))

matt_q$Year <- matt_q$Year - 1 # Release year offset by 1

gMat <- matt_q %>%
  dplyr::filter(Age > 1) %>%
  ggplot(aes(Year, .data$`50%`, fill = factor(.data$Age),
             colour = factor(.data$Age))) +
  geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
  # geom_hline(data = bmatt, linetype = 2, aes(yintercept = .data$value,
  #                                            colour = factor(.data$Age))) +
  labs(x = ifelse(brood, "Brood Year", "Return Year"), y = "Proportion mature",
       colour = "Age", fill = "Age") +
  annotate( "text",  x = -Inf, y = Inf, label = "c",  hjust = -1,
            vjust = 1.8, size = 6, fontface = "bold")



# # For Q/C take >2010 values (approximately last 10 years)
# matt_q <- rename(matt_q, median='50%')#|> filter(Year>2010)
# out <- matt_q |> group_by(Age) |> summarise(mean_ppn = mean(median, na.rm=T))
# matt_q
# out
#
# age3 <- out$mean_ppn[3]-out$mean_ppn[2]
# age4 <- out$mean_ppn[4]-out$mean_ppn[3]
# age5 <- out$mean_ppn[5]-out$mean_ppn[4]
# age3
# age4
# age5

# Plot productivity with horizontal line at zero
prod <- salmonMSE:::.CM_prod(report, d)

prod_q <- apply(prod, 1, quantile, probs = c(0.025, 0.5, 0.975),
                na.rm = TRUE) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 + year1 - 1) %>%
  reshape2::dcast(list("Year", "Var1"))

# mean(prod_q[,3])
assign(paste0("prod_", pop), prod_q)

gProd <- prod_q %>%
  ggplot(aes(Year, .data$`50%`)) +
  geom_line() +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
  labs(x = "Calendar Year", y = "Productivity") +
  geom_hline(yintercept = 1) +
  annotate( "text",  x = -Inf, y = Inf, label = "d",  hjust = -1,
            vjust = 1.8, size = 6, fontface = "bold")


# Plot ERs
type <- "PT"
brood <- FALSE
ER_list <- salmonMSE:::.CM_ER(report, type = type,
                              index_AEQ = match(year_borrow, year),
                              brood = FALSE, r = 1, simplify = FALSE)

gERPT <- salmonMSE:::.CM_ts(
  ER_list,
  year1 = if (brood) year1 - 1 else year1,
  ci,
  var = "ER",
  xlab = ifelse(brood, "Brood Year", "Calendar Year"),
  ylab = switch(
    type,
    "PT" = "Preterminal exploitation rate (AEQ)",
    "T" = "Terminal exploitation rate (AEQ)",
    "all" = "Total exploitation rate (AEQ)")
  ) +
  annotate( "text",  x = -Inf, y = Inf, label = "e", hjust = -1,
            vjust = 1.8, size = 6, fontface = "bold")


type <- "T"
brood <- FALSE
ER_list <- salmonMSE:::.CM_ER(report, type = type,
                              index_AEQ = match(year_borrow, year),
                              brood = FALSE, r = 1, simplify = FALSE)

gERT <- salmonMSE:::.CM_ts(
  ER_list,
  year1 = if (brood) year1 - 1 else year1,
  ci,
  var = "ER",
  xlab = ifelse(brood, "Brood Year", "Calendar Year"),
  ylab = switch(
    type,
    "PT" = "Preterminal exploitation rate (AEQ)",
    "T" = "Terminal exploitation rate (AEQ)",
    "all" = "Total exploitation rate (AEQ)")) +
  annotate( "text",  x = -Inf, y = Inf, label = "f",  hjust = -1,
            vjust = 1.8, size = 6, fontface = "bold")



# gERPT <- CM_ER(report, brood = FALSE, type = "PT", year1, ci = TRUE, r = 1,
#                at_age = FALSE, index_AEQ = match(year_borrow, year))
#
# gERT <- CM_ER(report, brood = FALSE, type = "T", year1, ci = TRUE, r = 1,
#               at_age = FALSE, index_AEQ = match(year_borrow, year))


# Get pHOS
pHOSeff <- sapply(report, getElement, "pHOSeff") %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 + year1 - 1) %>%
  reshape2::dcast(list("Year", "Var1"), value.var = "value")

gHOS <- ggplot(pHOSeff, aes(.data$Year, .data$`50%`)) +
  geom_line() +
  labs(x = "Year", y = "pHOSeff") +
  expand_limits(y = 0) +
  ylim(0, 1) +
  annotate( "text",  x = -Inf, y = Inf, label = "g",
            hjust = -1, vjust = 1.8, size = 6, fontface = "bold")

gHOS <- gHOS + geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), fill = alpha("grey", 0.5))

# Get PNI
# brood <- sapply(report, getElement, "brood")
pNOB <- sapply(report, getElement, "pNOB")
pHOSeff <- sapply(report, getElement, "pHOSeff")
PNIx <-  pNOB/(pNOB + pHOSeff)

PNI <- PNIx %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 + year1 - 1) %>%
  reshape2::dcast(list("Year", "Var1"), value.var = "value")

gPNI <- ggplot(PNI, aes(.data$Year, .data$`50%`)) +
  geom_line() +
  labs(x = "Year", y = "PNI") +
  expand_limits(y = 0) +
  ylim(0, 1) +
  annotate( "text",  x = -Inf, y = Inf, label = "h",
            hjust = -1, vjust = 1.8, size = 6, fontface = "bold")


gPNI <- gPNI +
  geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`),
              fill = alpha("grey", 0.5))
if(pop!= "QC") {
  gPNI <- gPNI +
    geom_point() +
    geom_errorbar(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`), width = 0.1)
}

if(pop!="QC") {
  gSp <- get(paste0("gSp_", pop))
  gU <- get(paste0("gU_", pop))
 }

gCM1 <- (gSRR + gSurv)/
  (gMat + gProd)
if(pop != "Adam"){
  gCM2 <- (gERPT + gERT)/
  (gHOS + gPNI)
}

if(pop == "Adam"){
  gCM2 <- (gERPT + plot_spacer())/
    (gSp + gU) +
    plot_layout(heights = c(1,1))
}

if(pop != "QC"){
  gCM3 <- (gSp + gU)/
    plot_spacer() +
    plot_layout(heights = c(1, 1))
}


if(write_figures) ggsave(paste("figures/", pop, "_gCM1.png", sep=""), gCM1, height = 7, width = 10)

if(write_figures) ggsave(paste("figures/", pop, "_gCM2.png", sep=""), gCM2, height = 7, width = 10)

if(!pop %in% c("QC", "Adam") & write_figures == TRUE) ggsave(paste("figures/", pop, "_gCM3.png", sep=""), gCM3, height = 7, width = 10)





if (all( vapply(
  c("prod_QC", "prod_Woss", "prod_Salmon", "prod_Adam"), exists, logical(1)))
  ) {
  prod_QC$pop <- "Quinsam/Campbell"
  prod_Woss$pop <- "Nimpkish"
  prod_Adam$pop <- "Adam"
  prod_Salmon$pop <- "Salmon"
  df <- rbind(prod_QC, prod_Woss, prod_Adam, prod_Salmon)
  gProdAll <- ggplot(df, aes(.data$Year, .data$`50%`, group=pop, colour=pop)) +
    geom_line() +
    labs(x = "Year", y = "Productivity (recruits/spawner)") +
    expand_limits(y = 0) +
    labs(colour = NULL)

  if(write_figures) ggsave(here("figures/ProdAll.png"), gProdAll, height = 5, width = 8)
}

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

