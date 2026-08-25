# Prior and posterior plots

library(scales)
library(tidyverse)
library(ggplot2)


theme_set(
  theme_bw() +
    theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
)

# Input results
ERM_QuinsamCampbell <- readRDS("CM/QuinsamCampbell_07.29.26.rds")
report_QC <- salmonMSE:::get_report(ERM_QuinsamCampbell)

ERM_Adam <- readRDS("CM/Adam_08.06.26.prior.rds")#readRDS("CM/Adam_06.01.26.JSpt.rds")
report_Adam <- salmonMSE:::get_report(ERM_Adam)

ERM_Salmon <- readRDS("CM/Salmon_08.18.26.prior.rds") #file:///C:/github/UpperSoG_Chinook/CM/Salmon_06.10.JSesc.NGTSt.html
# ERM_Salmon <- readRDS("CM/Salmon_08.22.26.altprior.rds") #file:///C:/github/UpperSoG_Chinook/CM/Salmon_06.10.JSesc.NGTSt.html
report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

ERM_Woss <- readRDS("CM/Woss_07.22.26.prior.rds")
report_Woss <- salmonMSE:::get_report(ERM_Woss)

# Set up population
pop <- "Salmon" #"Adam" #"Salmon"#"Woss"
pop.prior <- pop
if(pop=="Adam") pop.prior <- "Adam/Eve"
if(pop=="Woss") pop.prior <- "Nimpkish"

report <- report_Salmon#report_QC#report_SalmonPhillips#report_AdamPhillips#report_WossPhillips
year1 <- 2002#2002#2001#1984#2010#1984
samp <- ERM_Salmon#readRDS(paste("CM/",pop,"_06.01.26.JSpt.rds", sep=""))#paste("ERM_", pop, "Phillips", sep="")
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)

# Srep prior
data_Srep_prior <- as.data.frame( read.csv(
  ("data/UpperSoGChinook_out_posteriorpredictive_NEWWArev.csv"))
  # ("data/UpperSoGChinook_out_posteriorpredictive_UPDATEDAWA_Aug18.csv"))
)
Srep_prior <- data_Srep_prior %>% filter(Stock==pop.prior) %>% pull(SREP_median)
logSrep_prior_sd <- data_Srep_prior %>% filter(Stock==pop.prior) %>%
  mutate(sigma=(log(SREP_upr95)-log(SREP_median))/2) %>%
  pull(sigma)


post_log_so <- log(sapply(report, getElement, "so"))# for egg-smolt rel
len <- length(post_log_so)
set.seed(234)
prior_log_so <- rnorm(len, log(Srep_prior), logSrep_prior_sd)

df <- data.frame(so = c(post_log_so, prior_log_so),
                 label = c(rep("Posterior", len), rep("Habitat-based\n(as prior)", len)))

g <- ggplot(df, aes(x=so, colour=label, fill=label)) +
  geom_density(alpha = 0.2) +
  # scale_x_continuous(trans = scales::exp_trans()) +
  # scale_x_continuous(labels = function(x) round(exp(x), 0)) +
  scale_x_continuous(
    breaks = log(c(1000, 2000, 5000,10000, 20000, 50000)),  # choose meaningful raw values
    labels = c("1000", "2000", "5000", "10,000", "20,000", "50,000")
  ) +
  theme(legend.title = element_blank()) +
  geom_vline(xintercept = mean(post_log_so), colour = "#00BFC4", linetype="dashed" ) +
  geom_vline(xintercept = log(Srep_prior), colour = "#F8766D" , linetype="dashed") +
  xlab("Unfished equilibrium spawner abundances") +
  ylab("Density")
g

ggsave(paste("figures/", pop, "_priorpost.png", sep=""), g, height = 3.5, width = 6)

# Plot Smax from IWAM vs posterior
# Srep prior
Smax_IWAM <- data_Srep_prior %>% filter(Stock==pop.prior) %>% pull(SMAX_median)
logSmax_IWAM_sd <- data_Srep_prior %>% filter(Stock==pop.prior) %>%
  mutate(sigma=(log(SMAX_upr95)-log(SMAX_median))/2) %>%
  pull(sigma)
set.seed(234)
logSmax_IWAM <- rnorm(len, log(Smax_IWAM), logSmax_IWAM_sd)

# Posterior Smax

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
beta_post <- get_beta_s(report, samp)
Smax_post <- 1/beta_post



df <- data.frame(so = c(log(Smax_post), logSmax_IWAM),
                 label = c(rep("Posterior", len), rep("Habitat-based", len)))

if(pop=="Woss"){
  g <- ggplot(df, aes(x=so, colour=label, fill=label)) +
    geom_density(alpha = 0.2) +
    # scale_x_continuous(trans = scales::exp_trans()) +
    # scale_x_continuous(labels = function(x) round(exp(x), 0)) +
    scale_x_continuous(
      breaks = log(c(1000, 2000, 5000,10000, 20000, 50000)),  # choose meaningful raw values
      labels = c("1000", "2000", "5000", "10,000", "20,000", "50,000")
    ) +
    theme(legend.title = element_blank()) +
    geom_vline(xintercept = median(log(Smax_post)), colour = "#00BFC4" , linetype="dashed") +
    geom_vline(xintercept = median(logSmax_IWAM), colour = "#F8766D", linetype="dashed" ) +
    xlab("Smax") +
    ylab("Density")
  g

}
if(pop=="Salmon"|pop=="Adam"){
  g <- ggplot(df, aes(x=so, colour=label, fill=label)) +
    geom_density(alpha = 0.2) +
    # scale_x_continuous(trans = scales::exp_trans()) +
    # scale_x_continuous(labels = function(x) round(exp(x), 0)) +
    scale_x_continuous(
      breaks = log(c(200, 500, 1000, 2000, 5000,10000, 20000)),  # choose meaningful raw values
      labels = c("200", "500", "1000", "2000", "5000", "10,000", "20,000")
    ) +
    theme(legend.title = element_blank()) +
    geom_vline(xintercept = median(log(Smax_post)), colour = "#00BFC4" , linetype="dashed") +
    geom_vline(xintercept = median(logSmax_IWAM), colour = "#F8766D", linetype="dashed" ) +
    xlab("Smax") +
    ylab("Density") +
    coord_cartesian(xlim = c(log(200), log(20000)))
  g

}

ggsave(paste("figures/", pop, "_Smax_IWAM-post.png", sep=""), g, height = 3.5, width = 6)


# Show productivity
# Extract full life cycle alpha values from conditioning model output
get_alpha_s <- function(report, samp){
  d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
  alpha_s <- salmonMSE:::.CM_prod(report, d, mean_bio = TRUE) # Ricker alpha, per spawner
  return(as.vector(alpha_s))
}
alpha_post <- get_alpha_s(report, samp)

df <- data.frame(prod = alpha_post,
                 label = "productivity")

g <- ggplot(df, aes(x=alpha_post, colour=label, fill = label)) +
  geom_density(alpha = 0.2) +
  theme(legend.position = "none") +
  xlab("Productivity (recruits/spawner)") +
  ylab("Density") +
  geom_vline(xintercept = median(alpha_post), colour = "#F8766D",
             linetype="dashed" ) +
  coord_cartesian(xlim = c(0, 15))

g

ggsave(paste("figures/", pop, "_prod.png", sep=""), g, height = 3.5, width = 6)

