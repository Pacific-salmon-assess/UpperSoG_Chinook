#-------------------------------------------------------------------------------
# Plot Priors and Posterior of spawner abundances at equilibrium
#-------------------------------------------------------------------------------

# INPUTS:

# Conditioning models outputs from population-specific "01" files.
# See "Input Conditioning Model results" section below.
#
# Habitat-based prior on unfished spawners at equilibrium
# data/UpperSoGChinook_posteriorpredictive_Sept25_iter1.csv

# Population productivity (logalpha): base, low, or high

# OUTPUTS:
# Figures of priors and posteriors, and of posterior estimate of productivity
# In folders: figures/CM (for base prod), figures/CM/low, or figures/CM/high for
  #low and high productivity
# File names:
# Woss_priorpost.png, Salmon_priorpost.png, Adam_priorpost.png
# Woss_prod.png, Salmon_prod.png, Adam_prod.png

#-------------------------------------------------------------------------------
# Libraries
library(scales)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Plot priors and posteriors from all prod scenarios together, otherwise
# plotted for a single population and single scenario
meta.plot <- FALSE

theme_set(
  theme_bw() +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
)




#-------------------------------------------------------------------------------
# Set up scenarios

for (prod in c("low", "base", "high")){
  # Input Conditioning Model results
  ERM_QC <- readRDS(paste0("CM/QuinsamCampbell_09.29.26.", prod, ".rds"))
  report_QC <- salmonMSE:::get_report(ERM_QC)

  ERM_Adam <- readRDS(paste0("CM/Adam_09.29.26.", prod, ".rds"))
  report_Adam <- salmonMSE:::get_report(ERM_Adam)

  ERM_Salmon <- readRDS(paste0("CM/Salmon_09.29.26.", prod, ".rds"))
  report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

  ERM_Woss <- readRDS(paste0("CM/Woss_09.29.26.", prod, ".rds"))
  report_Woss <- salmonMSE:::get_report(ERM_Woss)


  # Create folders for output figures
  if(prod == "base") folder_path <- "figures/CM/"
  if(prod == "low") folder_path <- "figures/CM/low/"
  if(prod == "high") folder_path <- "figures/CM/high/"

  if (!dir.exists(here::here(folder_path))) {
    dir.create(here::here(folder_path), recursive = TRUE)
  }

  for (pop in c("Adam", "Salmon", "Woss")){
    pop.prior <- pop
    if(pop=="Adam") pop.prior <- "Adam/Eve"
    if(pop=="Woss") pop.prior <- "Nimpkish"
    report <- get(paste0("report_", pop))
    samp <- get(paste0("ERM_", pop))
    d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
    if(pop == "QC") year1 <- 1984
    if(pop == "Woss") year1 <- 2001
    if(pop == "Salmon") year1 <- 2002
    if(pop == "Adam") year1 <- 2002

    if(pop == "Woss") {
      if (prod == "low") lab <- "a"
      if (prod == "base") lab <- "b"
      if (prod == "high") lab <- "c"
    }
    if(pop == "Salmon") {
      if (prod == "low") lab <- "d"
      if (prod == "base") lab <- "e"
      if (prod == "high") lab <- "f"
    }
    if(pop == "Adam") {
      if (prod == "low") lab <- "g"
      if (prod == "base") lab <- "h"
      if (prod == "high") lab <- "i"
    }

    # Smax prior plots
    data_Smax_prior <- as.data.frame( read.csv(
      ("data/UpperSoGChinook_posteriorpredictive_Sept25_iter1.csv"))
       # ("data/UpperSoGChinook_posteriorpredictive_Sept25_iter3.csv"))
    )
    mean_Smax_prior <- data_Smax_prior %>% filter(Stock == pop.prior) %>% pull(logSMAX_mean) %>% exp()
    logSmax_prior_sd <- data_Smax_prior %>% filter(Stock == pop.prior) %>%
      mutate(sigma=(SMAX_upr975-logSMAX_median)/2) %>%
      pull(sigma)
    if (pop == "QC") {
      mean_Smax_prior <- exp(d$smax_mu)
      logSmax_prior_sd <- d$smax_sd
    }

    post_log_smax <- log(sapply(report, getElement, "smax"))# for egg-smolt rel
    len <- length(post_log_smax)
    set.seed(234)
    prior_log_smax <- rnorm(len, mean(log(mean_Smax_prior)), logSmax_prior_sd)

    df <- data.frame(smax = c(prior_log_smax, post_log_smax),
                     label = c(rep("Prior", len), rep("Posterior", len)))

    if(pop=="Woss") xlimits <- c(log(2000), log(40000))
    if(pop=="Salmon") xlimits <- c(log(500), log(10000))
    if(pop=="Adam") xlimits <- c(log(500), log(10000))
    if(pop=="QC") xlimits <- c(log(500), log(100000))

    p <- ggplot(df, aes(x=smax, colour=label, fill = label)) +
      geom_density()
    dens_data <- ggplot_build(p)$data[[1]]
    ylab_coord <- max(dens_data$y)


    gSmaxPP <- ggplot(df, aes(x=smax, colour=label, fill=label)) +
      geom_density(alpha = 0.2) +
      theme(legend.position = "none") +
      # scale_x_continuous(trans = scales::exp_trans()) +
      # scale_x_continuous(labels = function(x) round(exp(x), 0)) +
      scale_x_continuous(
        breaks = log(c(200, 500, 1000, 2000, 5000,10000, 20000)),  # choose meaningful raw values
        labels = c("200", "500", "1000", "2000", "5000", "10,000", "20,000")
      ) +
      theme(legend.title = element_blank()) +
      coord_cartesian(xlim = xlimits ) +
      geom_vline(xintercept = median(prior_log_smax), colour = "#00BFC4", linetype="dashed") +
      geom_vline(xintercept = median(post_log_smax), colour = "#F8766D" , linetype="dashed") +
      labs( x= NULL, y = NULL)

    if(!meta.plot){
      gSmaxPP <- gSmaxPP +
        xlab("Smax") +
        ylab("Density") +
        annotate( "text",  x = xlimits[1], y = ylab_coord, label = "a",
                  hjust = 0, vjust = 0, size = 4)
    }
    if(meta.plot){
      gSmaxPP <- gSmaxPP +
        annotate( "text",  x = xlimits[1] , y = ylab_coord, label = lab,
                  hjust = 0, vjust = 0, size = 4)

      if(lab == "h"){
        gSmaxPP <- gSmaxPP +
          xlab("Smax")
      }
      if(lab == "d"){
        gSmaxPP <- gSmaxPP +
          ylab("Density")
      }
      if(lab == "a"){
        gSmaxPP <- gSmaxPP +
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = "low",
            vjust = - 0.4,
            hjust = + 6
          ) +
          coord_cartesian(clip = "off")
      }
      if(lab == "b"){
        gSmaxPP <- gSmaxPP +
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = "base",
            vjust = - 0.4,
            hjust = + 4
          ) +
          coord_cartesian(clip = "off")
      }
      if(lab == "c"){
        gSmaxPP <- gSmaxPP +
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = "high",
            vjust = - 0.4,
            hjust = + 4.5
          ) +
          coord_cartesian(clip = "off")
      }
    }

    gSmaxPP
    # gWoss <- g + theme(legend.position = "none") + xlab(NULL) + ylab(NULL)
    # gSalmon <- g + theme(legend.position = "none") + xlab(NULL) + ylab(NULL)
    # gAdam <- g + theme(legend.position = "none") + xlab("Smax") + ylab("Density")


    # ggsave(paste0(folder_path, pop, "_priorpost.png"), g, height = 3.5, width = 6)

    # Show productivity
    # Extract full life cycle alpha values from conditioning model output
    get_alpha_s <- function(report, samp){
      d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
      alpha_s <- salmonMSE:::.CM_prod(report, d, mean_bio = TRUE) # Ricker alpha, per spawner
      return(as.vector(alpha_s))
    }
    logalpha_post <- log(get_alpha_s(report, samp))
    logalpha_post <- logalpha_post[logalpha_post > 0] # remove negative values

    logalpha_prior <- rnorm(length(logalpha_post)*2, d$cr_mu, d$cr_sd)
    logalpha_prior <- logalpha_prior[logalpha_prior > 0] # remove negative values

    df <- data.frame(prod = c(logalpha_post, logalpha_prior),
                     label = c(rep("Posterior", length(logalpha_post)),
                               rep("Prior", length(logalpha_prior))
                     )
    )

    p <- ggplot(df, aes(x=prod, colour=label, fill = label)) +
      geom_density()
    dens_data <- ggplot_build(p)$data[[1]]
    ylab_coord <- max(dens_data$y)

    gProdPP <- ggplot(df, aes(x=prod, colour=label, fill = label)) +
      geom_density(alpha = 0.2) +
      theme(legend.position = "none") +
      # theme(legend.title = element_blank()) +
       geom_vline(xintercept = median(logalpha_prior), colour = "#00BFC4",
                 linetype="dashed" ) +
      geom_vline(xintercept = median(logalpha_post), colour = "#F8766D",
                 linetype="dashed" ) +
      coord_cartesian(xlim = c(0, 2.5)) +
      labs( x= NULL, y = NULL)



    if(!meta.plot){
      gProdPP <- gProdPP +
        xlab("log productivity (log(recruits/spawner))") +
        ylab("Density") +
        annotate( "text",  x = min(df$prod), y = ylab_coord, label = "b",
                  hjust =  0, vjust = 0, size = 4) +
        theme(legend.title = element_blank())

    }
    if(meta.plot){
      gProdPP <- gProdPP +
        annotate( "text",  x =  min(df$prod) , y = ylab_coord, label = lab,
                  hjust = 0, vjust = 0, size = 4)

      if(lab == "h"){
        gProdPP <- gProdPP +
          xlab("log productivity (log(recruits/spawner))")
      }
      if(lab == "d"){
        gProdPP <- gProdPP +
          ylab("Density")
      }
      if(lab == "a"){
        gProdPP <- gProdPP +
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = "low",
            vjust = - 0.4,
            hjust = + 6
          ) +
          coord_cartesian(clip = "off")
      }
      if(lab == "b"){
        gProdPP <- gProdPP +
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = "base",
            vjust = - 0.4,
            hjust = + 4
          ) +
          coord_cartesian(clip = "off")
      }
      if(lab == "c"){
        gProdPP <- gProdPP +
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = "high",
            vjust = - 0.4,
            hjust = + 4.5
          ) +
          coord_cartesian(clip = "off")
      }
    }



    gProdPP
    gPP <- (gSmaxPP + gProdPP)


    assign(paste0("gProdPP_", pop, prod), gProdPP)
    assign(paste0("gSmaxPP_", pop, prod), gSmaxPP)

    # gWoss_prod <- g +  xlab(NULL) + ylab(NULL)
    # gSalmon_prod <- g + xlab(NULL) + ylab(NULL)
    # gAdam_prod <- g +  xlab("Productivity (log(recruits/spawner))") + ylab(NULL)

    # library(patchwork)
    # gComb <- (gWoss + gWoss_prod) /
    #         (gSalmon + gSalmon_prod)/
    #         (gAdam + gAdam_prod)
    # ggsave(paste("figures/CM/postprior_combined.png", sep=""), gComb, height = 8, width = 8)


    if(!meta.plot) ggsave(paste0(folder_path, pop, "_prod.png"), gPP, height = 4, width = 10)

  }

}

if(!meta.plot){
  gSmaxPP_meta <-
    (gSmaxPP_Wosslow + gSmaxPP_Wossbase + gSmaxPP_Wosshigh)/
    (gSmaxPP_Salmonlow + gSmaxPP_Salmonbase + gSmaxPP_Salmonhigh)/
    (gSmaxPP_Adamlow + gSmaxPP_Adambase + gSmaxPP_Adamhigh)

  gSmaxPP_meta

  gProdPP_meta <-
    (gProdPP_Wosslow + gProdPP_Wossbase + gProdPP_Wosshigh)/
    (gProdPP_Salmonlow + gProdPP_Salmonbase + gProdPP_Salmonhigh)/
    (gProdPP_Adamlow + gProdPP_Adambase + gProdPP_Adamhigh)

  gProdPP_meta
  ggsave(here("figures", "Smax_PP_meta.png"), gSmaxPP_meta, height = 8, width = 8)
  ggsave(here("figures", "Prod_PP_meta.png"), gProdPP_meta, height = 8, width = 8)

}

#-------------------------------------------------------------------------------
# Plot assumptions about productivity
#-------------------------------------------------------------------------------

if(FALSE){
  # Get base case productivity from Salmon River
  ERM <- readRDS(paste0("CM/Salmon_09.29.26.base.rds"))
  d <- salmonMSE:::get_CMdata(ERM@.MISC$CMfit)


  logalpha_prior <- rnorm(50000, d$cr_mu, d$cr_sd)
  logalpha_prior <- logalpha_prior[logalpha_prior > 0] # remove negative values


  cr_mu_lo <- log(exp(d$cr_mu) - 0.5)
  cr_mu_hi <- log(exp(d$cr_mu) + 0.5)
  logalpha_prior_lo <- rnorm(50000, cr_mu_lo, d$cr_sd)
  logalpha_prior_hi <- rnorm(50000, cr_mu_hi, d$cr_sd)


  df <- data.frame(prod = c(logalpha_prior_lo, logalpha_prior, logalpha_prior_hi),
                   label = c(rep("low", length(logalpha_prior_lo)),
                             rep("base", length(logalpha_prior)),
                             rep("high", length(logalpha_prior_hi))
                   )
  )


  gLogProdSA <- ggplot(df, aes(x=(prod), colour=label, fill = label)) +
    geom_density(alpha = 0.2) +
    theme(legend.position = "none") +
     xlab("log productivity (log(recruits/spawner))") +
    ylab("Density") +
    theme(legend.title = element_blank()) +
    geom_vline(xintercept = mean((logalpha_prior)), colour = "#F8766D",
              linetype="dashed" ) +
    geom_vline(xintercept = mean((logalpha_prior_hi)), colour = "#00BA38",
             linetype="dashed" ) +
    geom_vline(xintercept = mean((logalpha_prior_lo)), colour = "#619CFF",
               linetype="dashed" ) +
    annotate( "text",  x = -Inf, y = Inf, label = "a",
              hjust = -1, vjust = 1.8, size = 3)

    gRawProdSA <- ggplot(df, aes(x=exp(prod), colour=label, fill = label)) +
      geom_density(alpha = 0.2) +
      xlab("productivity (recruits/spawner)") +
      ylab("Density") +
      theme(legend.title = element_blank()) +
      geom_vline(xintercept = mean(exp(logalpha_prior)), colour = "#F8766D",
                 linetype="dashed" ) +
      geom_vline(xintercept = mean(exp(logalpha_prior_hi)), colour = "#00BA38",
                 linetype="dashed" ) +
      geom_vline(xintercept = mean(exp(logalpha_prior_lo)), colour = "#619CFF",
                 linetype="dashed" ) +
     annotate( "text",  x = -Inf, y = Inf, label = "b",
              hjust = -1, vjust = 1.8, size = 3) +
      xlim(0,12)

    gProdSA <- gLogProdSA + gRawProdSA
    gProdSA

  ggsave(here("figures", "prod_SensAnal.png"), gProdSA, height = 3, width = 6)


}

#-------------------------------------------------------------------------------
# Plot assumption about SR curve to derive Srep-prior from Smax-prior and alpha
#-------------------------------------------------------------------------------
if(FALSE){

# Set up population
pop <- "Salmon" #"Adam" #"Salmon"#"Woss"
pop.prior <- pop
if(pop=="Adam") pop.prior <- "Adam/Eve"
if(pop=="Woss") pop.prior <- "Nimpkish"
report <- get(paste0("report_", pop))
samp <- get(paste0("ERM_", pop))
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
if(pop == "QC") year1 <- 1984
if(pop == "Woss") year1 <- 2001
if(pop == "Salmon") year1 <- 2002
if(pop == "Adam") year1 <- 2002


# # Srep prior
# data_Srep_prior <- as.data.frame( read.csv(
#   ("data/UpperSoGChinook_out_posteriorpredictive_NEWWArev.csv"))
#   # ("data/UpperSoGChinook_out_posteriorpredictive_UPDATEDAWA_Aug18.csv"))
# )
# Srep_prior <- data_Srep_prior %>% filter(Stock==pop.prior) %>% pull(SREP_median)
# logSrep_prior_sd <- data_Srep_prior %>% filter(Stock==pop.prior) %>%
#   mutate(sigma=(log(SREP_upr95)-log(SREP_median))/2) %>%
#   pull(sigma)


data_Smax_prior <- as.data.frame( read.csv(
  ("data/UpperSoGChinook_posteriorpredictive_Sept25_iter1.csv"))
  # ("data/UpperSoGChinook_out_posteriorpredictive_NEWWArev.csv"))
  # ("data/UpperSoGChinook_out_posteriorpredictive_UPDATEDAWA_Aug18.csv"))
)
med_Smax_prior <- data_Smax_prior %>% filter(Stock == pop.prior) %>% pull(SMAX_median)
mean_Smax_prior <- data_Smax_prior %>% filter(Stock == pop.prior) %>% pull(SMAX_mean)
logSmax_prior_sd <- data_Smax_prior %>% filter(Stock == pop.prior) %>%
  mutate(sigma=(log(SMAX_upr95)-log(SMAX_median))/2) %>%
  pull(sigma)

# Productivity for Cowichan Chinook (= mean prod for Fraser river, 3 stocks)
# (Greenberg et al. in prep, Table S4)
mean_logalpha <- 0.87
logalpha_sig <- 0.23
mean_logalpha_IWAM <- 1.82910



logalpha <- rnorm(2000, mean_logalpha, logalpha_sig)
logalpha_IWAM <- rnorm(2000, mean_logalpha_IWAM, logalpha_sig)
logSmax <- rnorm(2000, log(med_Smax_prior), logSmax_prior_sd)
# logSmax <- rnorm(2000, log(mean_Smax_prior) - logSmax_prior_sd^2/2, logSmax_prior_sd)

Srep_prior <-  logalpha * exp(logSmax) # Srep = log(alpha)/beta, beta = 1/Smax
Srep_prior_IWAM <-  logalpha_IWAM * exp(logSmax) # Srep = log(alpha)/beta, beta = 1/Smax
logSrep_prior_sd <- sd(log(Srep_prior))

print(median(Srep_prior))

ylim <- mean(Srep_prior_IWAM) * 1.1

S <- seq(0, ylim, 100)
# R <- S * exp(mean_logalpha - (1/mean_Smax_prior) * S)

Rsamp <- matrix(nrow= length(S), ncol = 2000)
Rsamp_IWAM <- matrix(nrow= length(S), ncol = 2000)

for (i in 1:2000){
  Rsamp[,i] <- S * exp(logalpha[i] - (1/exp(logSmax[i])) * S)
  Rsamp_IWAM[,i] <- S * exp(logalpha_IWAM[i] - (1/exp(logSmax[i])) * S)
}
R <- apply(Rsamp, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t() %>% as.data.frame() %>%
  rename(Recruits = '50%') %>%
  rename(lower = '2.5%') %>%
  rename(upper = '97.5%')
R_IWAM <- apply(Rsamp_IWAM, 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = T) %>%
  t() %>% as.data.frame() %>%
  rename(Recruits = '50%') %>%
  rename(lower = '2.5%') %>%
  rename(upper = '97.5%')

R$Spawners <- S
R$Productivity <- "Local-recent"
R_IWAM$Spawners <- S
R_IWAM$Productivity <- "Regional-historical"

R <- rbind(R, R_IWAM)

g <- ggplot(R, aes(x=Spawners, y= Recruits, group = Productivity,
              colour = Productivity, fill = Productivity)) +
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "solid",
              linewidth = 1.2,
              col= grey(0.6)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = median(Srep_prior),
             linewidth = 1.2,
             linetype = "dashed",
             col = "#F8766D") +
  geom_vline(xintercept = median(med_Smax_prior),
             linewidth = 1.2,
             linetype = "dotted",
             col = grey(0.6)) +
  geom_vline(xintercept = median(Srep_prior_IWAM),
             linewidth = 1.2,
             linetype = "dashed",
             col = "#00BFC4") +

  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.3
  )


ggsave(paste("figures/SR_priors_SalmonRiver.png", sep=""), g, height = 4, width = 6)

# Save to ResDoc repository
# ggsave("C:/github/UpperSoG_Chinook_ResDoc/figures/SR_priors_SalmonRiver.png", g, height = 4, width = 6)
}

