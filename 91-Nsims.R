# How many simulations needed to stabilize results?
library(ggplot2)
library(here)

# See SMSE output for base case of Q/C
SMSE <- readRDS(file.path("SMSE","QC", paste0("QC_basecase.rds")))



# look at time series of PNI as this is most variable,
dim(SMSE@PNI)


# First Update PNI = 1 when there is no brood & pHOS = 0
Brood <- SMSE@HOB[,,5,] + SMSE@NOB[,,5,]
# Brood is age-selective. Take age-5s for purposes of setting PNI to 1
NoBrood <- Brood < 0.001
pHOS_zero <- SMSE@pHOS_effective < 0.001
# if(any(NoBrood)) print(i)
# if(any(pHOS_zero)) print(i)
SMSE@PNI[pHOS_zero] <- 1
SMSE@PNI[NoBrood] <- 1


# Use the last year of the time-series as the performance metric

x <- SMSE@PNI[,1,30]

cum_median <- sapply(seq_along(x), function(i) {
  median(x[1:i])
})

df <- data.frame(PNI=cum_median, Nsims=1:1000)
g <- ggplot(df, aes(x=Nsims, y=PNI)) +
  geom_point() +
  theme_set(theme_bw()) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_vline(xintercept = 500, linetype="dashed")

ggsave(g, file = here::here("figures", "Nsims.png"),  width = 7, height = 5)
cum_median[500]

