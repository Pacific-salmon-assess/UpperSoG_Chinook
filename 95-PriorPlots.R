# Prior and posterior plots

library(scales)

# Input results
ERM_QuinsamCampbell <- readRDS("CM/QuinsamCampbell_06.15.26.rds")
report_QC <- salmonMSE:::get_report(ERM_QuinsamCampbell)

ERM_Adam <- readRDS("CM/Adam_06.22.26.prior.rds")#readRDS("CM/Adam_06.01.26.JSpt.rds")
report_Adam <- salmonMSE:::get_report(ERM_Adam)

ERM_Salmon <- readRDS("CM/Salmon_06.22.26.prior.rds") #file:///C:/github/UpperSoG_Chinook/CM/Salmon_06.10.JSesc.NGTSt.html
report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

ERM_Woss <- readRDS("CM/Woss_06.22.26.prior.rds")
report_Woss <- salmonMSE:::get_report(ERM_Woss)

# Set up population
pop <- "Woss" #"Adam" #"Salmon"#"Woss"
pop.prior <- pop
if(pop=="Adam") pop.prior <- "Adam/Eve"
if(pop=="Woss") pop.prior <- "Nimpkish"

report <- report_Woss#report_QC#report_SalmonPhillips#report_AdamPhillips#report_WossPhillips
year1 <- 2001#2002#2001#1984#2010#1984
samp <- ERM_Woss#readRDS(paste("CM/",pop,"_06.01.26.JSpt.rds", sep=""))#paste("ERM_", pop, "Phillips", sep="")
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)

# Smax prior
data_Smax_prior <- as.data.frame( read.csv(
  ("data/UpperSoGChinook_out_posteriorpredictive.csv"))
)
Smax_prior <- data_Smax_prior %>% filter(Stock==pop.prior) %>% pull(SREP_median)
logSmax_prior_sd <- data_Smax_prior %>% filter(Stock==pop.prior) %>%
  mutate(sigma=(log(SREP_upr95)-log(SREP_median))/2) %>%
  pull(sigma)


post_log_so <- log(sapply(report, getElement, "so"))# for egg-smolt rel
len <- length(post_log_so)
prior_log_so <- rnorm(len, log(Smax_prior), logSmax_prior_sd)

df <- data.frame(so = c(post_log_so, prior_log_so),
                 label = c(rep("Posterior", len), rep("Prior", len)))

g <- ggplot(df, aes(x=so, colour=label, fill=label)) +
  geom_density(alpha = 0.2) +
  # scale_x_continuous(trans = scales::exp_trans()) +
  # scale_x_continuous(labels = function(x) round(exp(x), 0)) +
  scale_x_continuous(
    breaks = log(c(1000, 2000, 5000,10000, 20000, 50000)),  # choose meaningful raw values
    labels = c("1000", "2000", "5000", "10,000", "20,000", "50,000")
  ) +
  theme(legend.title = element_blank()) +
  geom_vline(xintercept = log(Smax_prior), colour = "#00BFC4" , linetype="dashed") +
  geom_vline(xintercept = mean(post_log_so), colour = "#F8766D", linetype="dashed" ) +
  xlab("Unfished equilibrium spawner abundances") +
  ylab("Density")
g

ggsave(paste("figures/", pop, "_priorpost.png", sep=""), g, height = 3.5, width = 6)

