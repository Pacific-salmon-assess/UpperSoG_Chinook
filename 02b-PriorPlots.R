#-------------------------------------------------------------------------------
# Plot Priors and Posterior of spawner abundances at equilibrium
#-------------------------------------------------------------------------------

# INPUTS:

# Conditioning models outputs from population-specific "01" files.
# See "Input Conditioning Model results" section below.
#
# Habitat-based prior on unfished spawners at equilibrium
# data/UpperSoGChinook_posteriorpredictive_Sept25_iter1.csv

# Population productivity (logalpha): base, low or high

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

prod <- "high"#"base","low", "high #

theme_set(
  theme_bw() +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
)

# Create folders for output figures
if(prod == "base") folder_path <- "figures/CM/"
if(prod == "low") folder_path <- "figures/CM/low/"
if(prod == "high") folder_path <- "figures/CM/high/"

if (!dir.exists(here::here(folder_path))) {
  dir.create(here::here(folder_path), recursive = TRUE)
}

#-------------------------------------------------------------------------------
# Input Conditioning Model results
ERM_QC <- readRDS(paste0("CM/QuinsamCampbell_09.22.26.", prod, ".rds"))
report_QC <- salmonMSE:::get_report(ERM_QC)

ERM_Adam <- readRDS(paste0("CM/Adam_09.26.26.", prod, ".rds"))
report_Adam <- salmonMSE:::get_report(ERM_Adam)

ERM_Salmon <- readRDS(paste0("CM/Salmon_09.26.26.", prod, ".rds"))
report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

ERM_Woss <- readRDS(paste0("CM/Woss_09.26.26.", prod, ".rds"))
report_Woss <- salmonMSE:::get_report(ERM_Woss)


#-------------------------------------------------------------------------------
# Set up population
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


  # # Srep prior (use this for high productivity assumption from IWAM)
  # data_Srep_prior <- as.data.frame( read.csv(
  #   ("data/UpperSoGChinook_out_posteriorpredictive_NEWWArev.csv"))
  #   # ("data/UpperSoGChinook_out_posteriorpredictive_UPDATEDAWA_Aug18.csv"))
  # )
  # Srep_prior <- data_Srep_prior %>% filter(Stock==pop.prior) %>% pull(SREP_median)
  # logSrep_prior_sd <- data_Srep_prior %>% filter(Stock==pop.prior) %>%
  #   mutate(sigma=(log(SREP_upr95)-log(SREP_median))/2) %>%
  #   pull(sigma)


#   # # Srep prior (use this for low productivity assumption from Greendberg et al)
#   data_Smax_prior <- as.data.frame( read.csv(
#     ("data/UpperSoGChinook_out_posteriorpredictive_NEWWArev.csv"))
#     #("data/UpperSoGChinook_out_posteriorpredictive_UPDATEDAWA_Aug18.csv"))
#   )
#   mean_Smax_prior <- data_Smax_prior %>% filter(Stock == pop.prior) %>% pull(SMAX_median)
#   logSmax_prior_sd <- data_Smax_prior %>% filter(Stock == pop.prior) %>%
#     mutate(sigma=(log(SMAX_upr95)-log(SMAX_median))/2) %>%
#     pull(sigma)
#
#   # Productivity for Cowichan Chinook (= mean prod for Fraser river, 3 stocks)
#   # (Greenberg et al. in prep, Table S4)
#   mean_logalpha <- 0.87
#   logalpha_sig <- 0.23
#
#   logalpha <- rnorm(2000, mean_logalpha, logalpha_sig)
#   logSmax <- rnorm(2000, log(mean_Smax_prior), logSmax_prior_sd)
#
#   Srep_prior <-  logalpha * exp(logSmax) # Srep = log(alpha)/beta, beta = 1/Smax
#   logSrep_prior_sd <- sd(log(Srep_prior))
#
#
#
#   post_log_so <- log(sapply(report, getElement, "so"))# for egg-smolt rel
#   len <- length(post_log_so)
#   set.seed(234)
#   prior_log_so <- rnorm(len, mean(log(Srep_prior)), logSrep_prior_sd)
#
#   df <- data.frame(so = c(post_log_so, prior_log_so),
#                    label = c(rep("Posterior", len), rep("Habitat-based\n(as prior)", len)))
#
#   g <- ggplot(df, aes(x=so, colour=label, fill=label)) +
#     geom_density(alpha = 0.2) +
#     # scale_x_continuous(trans = scales::exp_trans()) +
#     # scale_x_continuous(labels = function(x) round(exp(x), 0)) +
#     scale_x_continuous(
#       breaks = log(c(1000, 2000, 5000,10000, 20000, 50000)),  # choose meaningful raw values
#       labels = c("1000", "2000", "5000", "10,000", "20,000", "50,000")
#     ) +
#     theme(legend.title = element_blank()) +
#     geom_vline(xintercept = mean(post_log_so), colour = "#00BFC4", linetype="dashed" ) +
#     geom_vline(xintercept = mean(log(Srep_prior)), colour = "#F8766D" , linetype="dashed") +
#     xlab("Unfished equilibrium spawner abundances") +
#     ylab("Density") #+
#     # annotate( "text",  x = -Inf, y = Inf, label = "d", # for combined plot only
#     #            hjust = -1, vjust = 1.8, size = 6)
#
#   g
# #
# # gWoss <- gWoss + theme(legend.position = "none") +    xlab("Srep")
# # gWoss_srep <- gWoss_srep+ theme(legend.position = "none")+    xlab("Srep")
# # gAdam <- gAdam + theme(legend.position = "none")+    xlab("Srep")
# # gAdam_srep <- gAdam_srep + theme(legend.position = "none")+    xlab("Srep")
# # gSalmon <- g + theme(legend.position = "none")+    xlab("Srep")
# # gSalmon_srep <- gSalmon_srep + theme(legend.position = "none")+    xlab("Srep")
# #
# # library(patchwork)
# # gComb <- (gWoss_srep + gWoss) /
# #         (gSalmon_srep + gSalmon)/
# #         (gAdam_srep + gAdam)
# # ggsave(paste("figures/CM/postprior_combined.png", sep=""), gComb, height = 8, width = 6)
# #
#
#   ggsave(paste("figures/CM/", pop, "_priorpost.png", sep=""), g, height = 3.5, width = 6)

  # # Plot Smax from IWAM vs posterior
  # # Srep prior
  # Smax_IWAM <- data_Srep_prior %>% filter(Stock==pop.prior) %>% pull(SMAX_median)
  # logSmax_IWAM_sd <- data_Srep_prior %>% filter(Stock==pop.prior) %>%
  #   mutate(sigma=(log(SMAX_upr95)-log(SMAX_median))/2) %>%
  #   pull(sigma)
  # set.seed(234)
  # logSmax_IWAM <- rnorm(len, log(Smax_IWAM), logSmax_IWAM_sd)
  #
  # # Posterior Smax
  #
  # get_beta_s <- function(report, samp){
  #   d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
  #   alpha <- sapply(report, getElement, "alpha")# for egg-smolt rel
  #   beta <- sapply(report, getElement, "beta")# for egg-smolt rel
  #   alpha_s <- salmonMSE:::.CM_prod(report, d, mean_bio = TRUE) # Ricker alpha, per spawner
  #   epro <- t(alpha_s)/alpha # s, y
  #
  #   spro <- sapply(1:length(report), function(x) { # vector
  #     mo <- apply(report[[x]]$mo[, , drop = FALSE], 2, mean)
  #     matt <- apply(report[[x]]$matt[, , d$r_matt, drop = FALSE], 2, mean)
  #
  #     lo <- salmonMSE:::calc_survival(mo, matt) # smolt survival at replacement
  #     spro <- sum(lo * d$ssum * matt)
  #     return(spro)
  #   })
  #
  #   beta_s <- beta * epro/spro # Ricker beta, per spawner
  #   Srep <- log(t(alpha_s))/beta_s
  #   return(as.vector(beta_s))
  # }
  # beta_post <- get_beta_s(report, samp)
  # Smax_post <- 1/beta_post
  #
  #
  #
  # df <- data.frame(so = c(log(Smax_post), logSmax_IWAM),
  #                  label = c(rep("Posterior", len), rep("Habitat-based", len)))
  #
  # if(pop=="Woss"){
  #   g <- ggplot(df, aes(x=so, colour=label, fill=label)) +
  #     geom_density(alpha = 0.2) +
  #     # scale_x_continuous(trans = scales::exp_trans()) +
  #     # scale_x_continuous(labels = function(x) round(exp(x), 0)) +
  #     scale_x_continuous(
  #       breaks = log(c(1000, 2000, 5000,10000, 20000, 50000)),  # choose meaningful raw values
  #       labels = c("1000", "2000", "5000", "10,000", "20,000", "50,000")
  #     ) +
  #     theme(legend.title = element_blank()) +
  #     geom_vline(xintercept = median(log(Smax_post)), colour = "#00BFC4" , linetype="dashed") +
  #     geom_vline(xintercept = median(logSmax_IWAM), colour = "#F8766D", linetype="dashed" ) +
  #     xlab("Smax") +
  #     ylab("Density")
  #   g
  #
  # }
  # if(pop=="Salmon"|pop=="Adam"){
  #   g <- ggplot(df, aes(x=so, colour=label, fill=label)) +
  #     geom_density(alpha = 0.2) +
  #     # scale_x_continuous(trans = scales::exp_trans()) +
  #     # scale_x_continuous(labels = function(x) round(exp(x), 0)) +
  #     scale_x_continuous(
  #       breaks = log(c(200, 500, 1000, 2000, 5000,10000, 20000)),  # choose meaningful raw values
  #       labels = c("200", "500", "1000", "2000", "5000", "10,000", "20,000")
  #     ) +
  #     theme(legend.title = element_blank()) +
  #     geom_vline(xintercept = median(log(Smax_post)), colour = "#00BFC4" , linetype="dashed") +
  #     geom_vline(xintercept = median(logSmax_IWAM), colour = "#F8766D", linetype="dashed" ) +
  #     xlab("Smax") +
  #     ylab("Density") +
  #     coord_cartesian(xlim = c(log(200), log(20000)))
  #   g
  #
  # }
  #
  # ggsave(paste("figures/CM/", pop, "_Smax_IWAM-post.png", sep=""), g, height = 3.5, width = 6)
  #
  # #


  # Smax prior plots
  data_Smax_prior <- as.data.frame( read.csv(
    ("data/UpperSoGChinook_posteriorpredictive_Sept25_iter1.csv"))
    # ("data/UpperSoGChinook_posteriorpredictive_Sept25_iter3.csv"))
  # ("data/UpperSoGChinook_out_posteriorpredictive_NEWWArev.csv"))
    # ("data/UpperSoGChinook_out_posteriorpredictive_UPDATEDAWA_Aug18.csv"))
  )
  mean_Smax_prior <- data_Smax_prior %>% filter(Stock == pop.prior) %>% pull(logSMAX_mean) %>% exp()
  logSmax_prior_sd <- data_Smax_prior %>% filter(Stock == pop.prior) %>%
    mutate(sigma=(SMAX_upr975-logSMAX_median)/2) %>%
    pull(sigma)


  post_log_smax <- log(sapply(report, getElement, "smax"))# for egg-smolt rel
  len <- length(post_log_smax)
  set.seed(234)
  prior_log_smax <- rnorm(len, mean(log(mean_Smax_prior)), logSmax_prior_sd)

  df <- data.frame(smax = c(post_log_smax, prior_log_smax),
                   label = c(rep("Posterior", len), rep("Habitat-based\n(as prior)", len)))

  if(pop=="Woss") xlimits <- c(log(2000), log(40000))
  if(pop=="Salmon") xlimits <- c(log(500), log(10000))
  if(pop=="Adam") xlimits <- c(log(500), log(10000))
  g <- ggplot(df, aes(x=smax, colour=label, fill=label)) +
    geom_density(alpha = 0.2) +
    # scale_x_continuous(trans = scales::exp_trans()) +
    # scale_x_continuous(labels = function(x) round(exp(x), 0)) +
    scale_x_continuous(
      breaks = log(c(200, 500, 1000, 2000, 5000,10000, 20000)),  # choose meaningful raw values
      labels = c("200", "500", "1000", "2000", "5000", "10,000", "20,000")
    ) +
    theme(legend.title = element_blank()) +
    xlab("Smax") +
    ylab("Density") +
    coord_cartesian(xlim = xlimits ) +
    geom_vline(xintercept = median(post_log_smax), colour = "#00BFC4" , linetype="dashed") +
    geom_vline(xintercept = median(prior_log_smax), colour = "#F8766D", linetype="dashed" ) #+
    # annotate( "text",  x = -Inf, y = Inf, label = "e", # for combined plot only
    #           hjust = -1, vjust = 1.8, size = 6)

  g
  # gWoss <- g + theme(legend.position = "none") + xlab(NULL) + ylab(NULL)
  # gSalmon <- g + theme(legend.position = "none") + xlab(NULL) + ylab(NULL)
  # gAdam <- g + theme(legend.position = "none") + xlab("Smax") + ylab("Density")


  ggsave(paste0(folder_path, pop, "_priorpost.png"), g, height = 3.5, width = 6)

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


  g <- ggplot(df, aes(x=prod, colour=label, fill = label)) +
    geom_density(alpha = 0.2) +
    # theme(legend.position = "none") +
    xlab("log productivity (log(recruits/spawner))") +
    ylab("Density") +
    geom_vline(xintercept = median(logalpha_post), colour = "#F8766D",
               linetype="dashed" ) +
    geom_vline(xintercept = median(logalpha_prior), colour = "#00BFC4",
               linetype="dashed" ) +
    coord_cartesian(xlim = c(0, 2.5)) #+
    # annotate( "text",  x = -Inf, y = Inf, label = "d", # for combined plot only
    #           hjust = -1, vjust = 1.8, size = 6)

  g
  # gWoss_prod <- g +  xlab(NULL) + ylab(NULL)
  # gSalmon_prod <- g + xlab(NULL) + ylab(NULL)
  # gAdam_prod <- g +  xlab("Productivity (log(recruits/spawner))") + ylab(NULL)

  # library(patchwork)
  # gComb <- (gWoss + gWoss_prod) /
  #         (gSalmon + gSalmon_prod)/
  #         (gAdam + gAdam_prod)
  # ggsave(paste("figures/CM/postprior_combined.png", sep=""), gComb, height = 8, width = 8)


  ggsave(paste0(folder_path, pop, "_prod.png"), g, height = 3.5, width = 6)

}


#-------------------------------------------------------------------------------
# Plot assumptions about productivity
#-------------------------------------------------------------------------------

if(FALSE){
  pop <- "Salmon"
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


  logalpha_prior <- rnorm(10000, d$cr_mu, d$cr_sd)
  logalpha_prior <- logalpha_prior[logalpha_prior > 0] # remove negative values

  sqrt(exp(d$cr_sd^2) - 1)

  cr_mu_lo <- log(exp(d$cr_mu) - 0.5)
  cr_mu_hi <- log(exp(d$cr_mu) + 0.5)
  logalpha_prior_lo <- rnorm(10000, cr_mu_lo, d$cr_sd)
  logalpha_prior_hi <- rnorm(10000, cr_mu_hi, d$cr_sd)


  # sig_raw <- sd(exp(logalpha_prior_hi))
  # sig_raw
  # sqrt(exp(sig_raw^2-1))


  df <- data.frame(prod = c(logalpha_prior_lo, logalpha_prior, logalpha_prior_hi),
                   label = c(rep("low", length(logalpha_prior_lo)),
                             rep("base", length(logalpha_prior)),
                             rep("high", length(logalpha_prior_hi))
                   )
  )


  g <- ggplot(df, aes(x=(prod), colour=label, fill = label)) + #exp(prod)
    geom_density(alpha = 0.2) +
    # theme(legend.position = "none") +
    # xlab("productivity (recruits/spawner)") +
    xlab("log productivity (log(recruits/spawner))") +
    ylab("Density") +
    theme(legend.title = element_blank()) +
    geom_vline(xintercept = median((logalpha_prior)), colour = "#F8766D",
              linetype="dashed" ) +
    geom_vline(xintercept = median((logalpha_prior_hi)), colour = "#00BA38",
             linetype="dashed" ) +
    geom_vline(xintercept = median((logalpha_prior_lo)), colour = "#619CFF",
               linetype="dashed" )
  # annotate( "text",  x = -Inf, y = Inf, label = "d", # for combined plot only
  #           hjust = -1, vjust = 1.8, size = 6)

  g

  ggsave(paste("figures/CM/logprod_sensanal.png", sep=""), g, height = 3.5, width = 6)


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

