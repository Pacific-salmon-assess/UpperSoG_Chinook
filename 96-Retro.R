# Run retrospective peels

#See "01-QuinsamCampbell-camp-CM.R"
# Changed final year in full_matrix and full_year to Ryears, instead of 2024
# Changed output file to include final year in the name

library(tidyverse)
library(readxl)
library(salmonMSE)

Ryears <- c(2019, 2020, 2021, 2022, 2023)
RunRetro_QC(Ryears = Ryears)
year1 <- 1984


samp1 <- readRDS(file = paste("CM/QuinsamCampbell_06.19.26.", Ryears[1], ".rds", sep=""))
samp2 <- readRDS(file = paste("CM/QuinsamCampbell_06.19.26.", Ryears[2], ".rds", sep=""))
samp3 <- readRDS(file = paste("CM/QuinsamCampbell_06.19.26.", Ryears[3], ".rds", sep=""))
samp4 <- readRDS(file = paste("CM/QuinsamCampbell_06.19.26.", Ryears[4], ".rds", sep=""))
samp5 <- readRDS(file = paste("CM/QuinsamCampbell_06.19.26.", Ryears[5], ".rds", sep=""))
samp <- readRDS(file = paste("CM/QuinsamCampbell_06.15.26.rds", sep=""))
#samp.x <- rep(0,length(Ryears))
for (i in (1:length(Ryears) ) ){
   # samp.x <- setNames(samp.x, paste0("samp", seq_along(i)))

  # setNames(values, paste0("element", seq_along(values)))

  samp.x <- get(paste0("samp", i))



  report <- salmonMSE:::get_report(samp.x)

  d <- salmonMSE:::get_CMdata(samp.x@.MISC$CMfit)
  # Repeat to fill in AEQs for incomplete brood years
  year <- year1 + seq(1, d$Ldyr) - 1
  year_borrow <- seq(max(year) - 9, max(year) - 5)


  # Get matrices of individual AEQ ER values by year x MCMC simulation.
  df.all <- salmonMSE:::.CM_ER(report, type = "all", r=1, brood=FALSE, index_AEQ = match(year_borrow, year))

  quantile(df.all, probs = c(0.025, 0.5, 0.975), na.rm=T)

  ER.all <- apply(df.all, 1, function(row) quantile(row, probs = c(0.025, 0.5, 0.975), na.rm=T))

  ER.x <- as.data.frame(t(ER.all))
  ER.x <- ER.x %>% rename("lower"="2.5%", "median"="50%" ,"upper"="97.5%") %>%
    mutate (label="ER", Year = seq(from = year1, length.out = length(ER.x[,1]))) %>%
    mutate(FinalYear = Ryears[i])


  assign(paste0("ER", Ryears[i]), ER.x)


}



# Get full data ER
report <- salmonMSE:::get_report(samp)
d <- salmonMSE:::get_CMdata(samp.x@.MISC$CMfit)
year <- year1 + seq(1, d$Ldyr) - 1
year_borrow <- seq(max(year) - 9, max(year) - 5)


df.all <- salmonMSE:::.CM_ER(report, type = "all", r=1, brood=FALSE, index_AEQ = match(year_borrow, year))
ER.all <- apply(df.all, 1, function(row) quantile(row, probs = c(0.025, 0.5, 0.975), na.rm=T))
ER2024 <- as.data.frame(t(ER.all))
ER2024 <- ER2024 %>% rename("lower"="2.5%", "median"="50%" ,"upper"="97.5%") %>%
  mutate (label="ER", Year = seq(from = year1, length.out = length(ER2024[,1]))) %>%
  mutate(FinalYear = 2024)


ER.retro <- rbind(ER2019, ER2020 , ER2021, ER2022, ER2023, ER2024)

# Look at individual years
ER.retro %>% filter(Year == "2020")
ER.retro %>% filter(Year == "2023")

#ggplot

g <- ER.retro %>% ggplot(aes(Year, median,  group=as.character(FinalYear))) +
  geom_line() +
  labs(x = "Calendar Year", y = "ER-total")+
  theme(legend.title = element_blank())
g

ggsave(paste("figures/QC_retro.jpg", sep=""), g, height = 3.5, width = 6)

