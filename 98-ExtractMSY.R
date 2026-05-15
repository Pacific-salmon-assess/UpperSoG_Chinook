# Extract MSY values from CM

# Input results
ERM_QuinsamCampbell <- readRDS("CM/QuinsamCampbell_05.01.26.rds")
report_QC <- salmonMSE:::get_report(ERM_QuinsamCampbell)

ERM_AdamPhillips <- readRDS("CM/AdamPhillips_CM_05.09.26.rds")
report_AdamPhillips <- salmonMSE:::get_report(ERM_AdamPhillips)

ERM_SalmonPhillips <- readRDS("CM/SalmonPhillips_CM_05.09.26.rds")
report_SalmonPhillips <- salmonMSE:::get_report(ERM_SalmonPhillips)

ERM_WossPhillips <- readRDS("CM/WossPhillips_CM_05.09.26.rds")
report_WossPhillips <- salmonMSE:::get_report(ERM_WossPhillips)

# Set up population
pop <- "Adam"#"Woss"#"Salmon"#"Adam"
report <- report_AdamPhillips#report_SalmonPhillips#report_AdamPhillips
year1 <- 2010#1981
samp <- readRDS(paste("CM/",pop,"Phillips_CM_05.09.26.rds", sep=""))#paste("ERM_", pop, "Phillips", sep="")
d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)


SMSY <- salmonMSE:::.CM_MSY(report, d, type = "spawner", ncores = 10) # year x simulation
#EMSY <- salmonMSE:::.CM_MSY(report, d, type = "egg", ncores = 10)
UMSY <- salmonMSE:::.CM_MSY(report, d, type = "u", ncores = 10)
Sgen <- salmonMSE:::.CM_Sgen(report, d, ncores = 10)

# gSMSY <- salmonMSE:::CM_MSY(report, d, type = "spawner", ncores = 7, year1=2010)
# gSgen <- salmonMSE:::CM_Sgen(report, d, ncores = 7 , year1=2010)
gUMSY <- salmonMSE:::CM_MSY(report, d, type = "u", ncores = 7 , year1=2010)
ggsave(paste("figures/", pop, "UMSY.png", sep=""), gUMSY, height = 3.5, width = 6)

# alpha <- sapply(report, getElement, "alpha")# for egg-smolt rel
# hist(alpha)
# mean(alpha)
# median(alpha)

#------------------------------------------------------------------------------
# MSY calculated with Ricker (lambert) equations

alpha_s <- salmonMSE:::.CM_prod(report, d) # Ricker alpha, per spawner
epro <- t(alpha_s)/alpha # egg per smolt: s, y
spro <- sapply(1:d$Ldyr, function(y) {
  sapply(1:length(report), function(x) {
    mo <- report[[x]]$mo[y, ]
    matt <- report[[x]]$matt[y, , d$r_matt]
    lo <- salmonMSE:::calc_survival(mo, matt) # smolt survival at replacement
    spro <- sum(lo * d$ssum * matt)
    return(spro)
  })
})
beta_s <- beta * epro/spro # Ricker beta, per spawner
SMSY_s <-  (1 - gsl:::lambert_W0(exp(1 - log(t(alpha_s)))))/beta_s
UMSY_s <-  (1 - gsl::lambert_W0(exp(1 - log(t(alpha_s)))))

calc_Sgen_Ricker <- function(loga, b){
  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
  a <- exp(loga)

  return(-1/b*gsl::lambert_W0(-b*sMSY/a))
}
Sgen_s <- calc_Sgen_Ricker(log(t(alpha_s)), beta_s)

SMSY_s_q <-  apply(t(SMSY_s), 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = FALSE) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 + year1 - 1) %>%
  reshape2::dcast(list("Year", "Var1")) %>%
  rename(median='50%', lower='2.5%', upper='97.5%') %>%
  mutate(label="SMSY")

# Plots
g1 <- SMSY_s_q %>%
  ggplot(aes(Year, .data$median, colour= label, fill=label)) +
  geom_line() +
  scale_fill_manual(values = c("SMSY" = "chartreuse4")) +
  scale_colour_manual(values = c("SMSY" = "chartreuse4")) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color= NA) +
  labs(x = "Calendar Year", y = ylab)+
  theme(legend.title = element_blank()) +
  ylab("Spawners")

g1
ggsave(paste("figures/", pop, "SMSY_calc_v1.png", sep=""), g1, height = 3.5, width = 6)

Spanwers_q <- data.frame(Year=MSY_q$Year,
                         lower=rep(NA,length(MSY_q$Year)),
                         median=d$obsescape,
                         upper=rep(NA, length(MSY_q$Year)),
                         label="Spawners")
SMSY_s_q <- rbind(SMSY_s_q, Spawners_q)

g2 <- SMSY_s_q %>%
  ggplot(aes(Year, .data$median, colour= label, fill=label)) +
  geom_line() +
  scale_fill_manual(values = c("SMSY" = "chartreuse4", "Spawners" = NA)) +
  scale_colour_manual(values = c("SMSY" = "chartreuse4", "Spawners" = "black")) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color= NA) +
  labs(x = "Calendar Year", y = ylab)+
  coord_cartesian(ylim = c(-2500,3000))+#c(-6000,4500)) +
  theme(legend.title = element_blank()) +
  ylab("Spawners")

g2
ggsave(paste("figures/", pop, "SMSY_calc_v2.png", sep=""), g2, height = 3.5, width = 6)

# SMSY: removing negative values
na.rm <- TRUE
if (na.rm) SMSY_s[SMSY_s <= 0] <- NA_real_

SMSY_s_q <- apply(t(SMSY_s), 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = na.rm) %>%
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
  coord_cartesian(ylim = c(0, 3500)) +# c(0,4500)) +
  theme(legend.title = element_blank()) +
  ylab("Spawners")

g3
ggsave(paste("figures/", pop, "SMSY_calc_v3.png", sep=""), g3, height = 3.5, width = 6)



# Adding Sgen, and removing negative values from calculated Sgen values
na.rm <- TRUE
if (na.rm) Sgen_s[Sgen_s <= 0] <- NA_real_

Sgen_s_q <- apply(t(Sgen_s), 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
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
  coord_cartesian(ylim = c(0, 3500)) + #c(0,4500)) +
  ylab("Spawners") +
  theme(legend.title = element_blank())

g4
ggsave(paste("figures/", pop, "SMSY_calc_v4.png", sep=""), g4, height = 3.5, width = 6)

meanSMSY <- SMSY_s_q %>% filter(label== "SMSY") %>% summarize(mean=mean(median)) %>% pull(mean)
meanSgen <- Sgen_s_q %>% filter(label== "Sgen") %>% summarize(mean=mean(median)) %>% pull(mean)

g5 <- g4 +
  geom_hline(yintercept = meanSMSY, lty="dashed", colour = "chartreuse4") +
  geom_hline(yintercept = meanSgen, lty="dashed", colour = "darkorange")
g5

ggsave(paste("figures/", pop, "SMSY_calc_v5.png", sep=""), g5, height = 3.5, width = 6)

#Remove negative UMSY from calcualted values (where prod<1)
UMSY_s[UMSY_s <= 0] <- NA_real_
UMSY_s_q <-  apply(t(UMSY_s), 1, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE) %>%
  reshape2::melt() %>%
  mutate(Year = Var2 + year1 - 1) %>%
  reshape2::dcast(list("Year", "Var1")) %>%
  rename(median='50%', lower='2.5%', upper='97.5%') %>%
  mutate(label="UMSY")

gUMSY <- UMSY_s_q %>%
  ggplot(aes(Year, .data$median, colour= label, fill=label)) +
  geom_line() +
  scale_fill_manual(values = c("UMSY" = "black")) +
  scale_colour_manual(values = c("UMSY" = "black")) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color= NA) +
  labs(x = "Calendar Year", y = ylab)+
  theme(legend.title = element_blank()) +
  ylab("UMSY")
gUMSY
ggsave(paste("figures/", pop, "UMSY_calc.png", sep=""), gUMSY, height = 3.5, width = 6)


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

Spanwers_q <- data.frame(Year=MSY_q$Year,
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
  coord_cartesian(ylim = c(0, 70000)) + #c(0, 22000)) + #c(0,4500)) +
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

#------------------------------------------------------------------------------
# Misc.
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





# Code comparing reference points
  require(salmonMSE)

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
