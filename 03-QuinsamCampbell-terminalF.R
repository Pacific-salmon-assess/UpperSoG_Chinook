
# Get Fraser Run Size time-series, and detrend and standardize

data <- as.data.frame(read.csv("data/Fraser Sockeye Run Size_2026-02-26.csv"))
# Data extracted from here: https://psc1.shinyapps.io/PSC_Annual_Fraser/ 26 Feb 2026

dat <- data %>% filter(Year>1980&Year<2024) %>%
  group_by(Year) %>% summarize(Run=sum(Run.Size))

trend <-  predict(lm(Run ~ Year, data = dat))
dat$detrended <-  dat$Run - trend
plot(dat$Run~dat$Year, type="l")
plot(dat$detrended~dat$Year, type="l")


dat <- dat %>%
  mutate(mean.Run=mean(detrended)) %>%
  mutate(sd.Run=sd(detrended)) %>%
  mutate(std.Run=(detrended-mean.Run)/sd.Run)


plot(dat$std.Run~dat$Year, type="l", xlab="Year", ylab="Standardized, detrended Run Size Fraser Sockeye")

# Create random time-series with 4 year cycle that is similar to Fraser Run size

set.seed(123)
nsims <- 3
t <- 1:43                                # time index
sd.peak <- 0.5#0.2
mean.peak <- 0.2
sd.noise <- 0.5#0.5
peak <- cycle <- noise <- y <- matrix(NA, nrow=nsims, ncol=length(t))

for (n in 1:nsims){
  phi <- runif(1, 0, 2*pi) #initial point in period
  peak[n,] <- rlnorm(43, meanlog = mean.peak - 0.5 * sd.peak^2, sdlog = sd.peak)

  cycle[n,] <- peak[n,] * sin(2 * pi * t / 4 + phi) - 1        # 4‑period cycle

  noise[n,] <- rlnorm(43, meanlog= -0.5 * sd.noise^2, sdlog = sd.noise)
  y[n,] <- cycle[n,] + noise[n,]


}
# plot(t, y[1,], type = "l", col="green")
# lines(t, y[2,], type = "l", col="blue")
# lines(t, y[3,], type = "l", col="red")


plot(dat$std.Run~dat$Year, type="n", xlab="Year", ylab="Standardized, detrended Run Size Fraser Sockeye")
lines(dat$Year, y[1,], type = "l", col="grey")
lines(dat$Year, y[2,], type = "l", col="grey")
lines(dat$Year, y[3,], type = "l", col="grey")
lines(dat$Year, dat$std.Run, type="l", col="black")

# multiply standardized, detrended randome values (y) by a proposed SD for terminal ERs, added to mean of 15%

ER <- y * 0.02 + 0.15 # Mean of 0.15 with sd=0.01

plot(dat$Year, ER[3,], type="n", xlab="Year", ylab="Stochastic terminal ERs with 4-year cycle from Fraser run size")
lines(dat$Year, ER[1,], type = "l", col="grey")
lines(dat$Year, ER[2,], type = "l", col="grey")
lines(dat$Year, ER[3,], type = "l", col="black")

# Do these look plausible?

