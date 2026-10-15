

library(tidyverse)
library(here)

# To update with Upper SoG escapement (Quinsam, Campbell, Adam, Salmon, Nimpkish) data once available


# Load the Excel file
dfa <- readxl::read_excel(
  file.path("data", "SOG_N_Escapement-Salmon_Adam_Nimpkish.xlsx"),
  sheet = "Data")
# Define the prefixes to match
prefixes <- c("Adam", "Nimpkish", "Salmon")

# Filter rows where Description starts with any of the prefixes (case-insensitive)
filtered_dfa <- dfa %>%
  filter(str_starts(Description, prefixes[1]) |
           str_starts(Description, prefixes[2]) |
           str_starts(Description, prefixes[3])) %>%
  mutate(Description = case_when( str_starts(Description, prefixes[1]) ~ prefixes[1],
    str_starts(Description, prefixes[2]) ~ prefixes[2],
    str_starts(Description, prefixes[3]) ~ prefixes[3])) %>%
  rename(year = "Analysis Year") %>%
  rename(escapement="Max Estimate") %>%
  rename(pop="Description") %>%
  mutate(IsUsed = ifelse(year < 2005, "no", "yes")) %>%
  select(year, escapement, pop, IsUsed)

# Add missing years to Nimpkish, where monitoring was not possible (2021,2022)

nimpkish_dfa <- filtered_dfa %>% filter(pop == "Nimpkish")
full_years_nimpkish <- data.frame(year = seq(min(nimpkish_dfa$year), max(nimpkish_dfa$year)))
df_complete_N <- left_join(full_years_nimpkish,nimpkish_dfa) %>%
  mutate(IsUsed = ifelse(year < 2005, "no", "yes")) %>%
  mutate(pop = "Nimpkish")



filtered_dfa_complete <- filtered_dfa %>%
  # Remove rows for Nimpkish
  filter(!pop %in% df_complete_N$pop) %>%
  # Add the updated rows
  bind_rows(df_complete_N)


# Load the Excel file
dfb <- readxl::read_excel(
  file.path("data", "Quinsam", "fsar-sog-cn-cq-nuseds.xlsx"),
  sheet = "Data")
# Define the prefixes to match
prefixes <- c("Quinsam", "Campbell")

# Filter rows where Description starts with any of the prefixes (case-insensitive)
filtered_dfb <- dfb %>%
  rename(WaterbodyName = "Waterbody Name") %>%
  filter(str_starts(str_to_lower(WaterbodyName), str_to_lower(prefixes[1])) |
           str_starts(str_to_lower(WaterbodyName), str_to_lower(prefixes[2]))) %>%
  mutate(WaterbodyName = case_when( str_starts(str_to_lower(WaterbodyName), str_to_lower(prefixes[1])) ~ prefixes[1],
                                  str_starts(str_to_lower(WaterbodyName), str_to_lower(prefixes[2])) ~ prefixes[2])) %>%
  rename(year = "Analysis Year") %>%
  rename(escapement="Max Estimate") %>%
  rename(pop="WaterbodyName") %>%
  mutate(IsUsed = ifelse(year < 2005, "no", "yes")) %>%
  select(year, escapement, pop, IsUsed)

filtered_df <- bind_rows(filtered_dfa_complete, filtered_dfb)
  # right_join(
  #   full_table %>% filter(Age == 1) %>% select(BROOD_YEAR),
  #   by = c("year" = "BROOD_YEAR")
  # )

g <- filtered_df %>%
  ggplot(aes(year, escapement)) +
  geom_line() +
  labs(x = "Year", y = "Total escapement") +
  #   scale_shape_manual(values = c(16, 1)) +
  facet_wrap(vars(pop), scales = "free_y", ncol = 3) +
  expand_limits(y = 0) +
  xlim(1985,2024) +
  geom_point(aes(col = IsUsed)) +
  theme(legend.position = 'bottom')

ggsave("figures/total_escapement.png", g, height = 6, width = 7)

# esc <- readr::read_csv("data/R-OUT_infilled_indicators_escapement_timeseries.csv")
#
# g <- esc %>%
#   mutate(river = gsub("_", " ", river) %>% stringr::str_to_title()) %>%
#   mutate(`Infill?` = ifelse(is.na(infill), "No", "Yes")) %>%
#   ggplot(aes(year, escapement)) +
#   geom_line() +
#   geom_point(aes(shape = `Infill?`)) +
#   labs(x = "Year", y = "Total escapement") +
#   scale_shape_manual(values = c(16, 1)) +
#   facet_wrap(vars(river), scales = "free_y", ncol = 3) +
#   expand_limits(y = 0) +
#   theme(legend.position = 'bottom')
# ggsave("figures/total_escapement.png", g, height = 9, width = 7)


#-------------------------------------------------------------------------------
# Plots of observed escapement and posterior estimates
# Quinsam Campbell
# Load conditioning model results - Quinsam/Campbell


png(here("figures", "escapements.png"))
par(mfrow = c(2, 2))

ERM <- readRDS("CM/QuinsamCampbell_09.29.26.base.rds")
report <- salmonMSE:::get_report(ERM)
d <- salmonMSE:::get_CMdata(ERM@.MISC$CMfit)
year1 <- 1984 # from conditioning model
year <- year1 + seq(1, d$Ldyr) - 1
obs <- d$obsescape
esc <- exp(sapply(report, getElement, "logpredesc")) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975))
plot(NULL, NULL, type = "n", pch = 16, xlim = range(year), ylim = c(0, 1.1) * range(esc, obs, na.rm = TRUE), xlab = "", ylab = "")
polygon(c(rev(year), year), c(rev(esc[1, ]), esc[3, ]), col = alpha("grey", 0.5), border = NA)
lines(year, esc[2, ], lwd = 1.5)
points(year, obs, pch = 16)
lines(year, obs, lty = 3)
mtext("(a) Quinsam/Cambpell Rivers", side = 3, adj = 0, line = 0.5)

ERM <- readRDS("CM/Woss_09.29.26.base.rds")
report <- salmonMSE:::get_report(ERM)
d <- salmonMSE:::get_CMdata(ERM@.MISC$CMfit)
year1 <- 2001 # from conditioning model
year <- year1 + seq(1, d$Ldyr) - 1
obs <- d$obsescape
esc <- exp(sapply(report, getElement, "logpredesc")) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975))
plot(NULL, NULL, type = "n", pch = 16, xlim = range(year), ylim = c(0, 1.1) * range(esc, obs, na.rm = TRUE), xlab = "", ylab = "")
polygon(c(rev(year), year), c(rev(esc[1, ]), esc[3, ]), col = alpha("grey", 0.5), border = NA)
lines(year, esc[2, ], lwd = 1.5)
points(year, obs, pch = 16)
lines(year, obs, lty = 3)
mtext("(b) Nimpkish River", side = 3, adj = 0, line = 0.5)

ERM <- readRDS("CM/Salmon_09.29.26.base.rds")
report <- salmonMSE:::get_report(ERM)
d <- salmonMSE:::get_CMdata(ERM@.MISC$CMfit)
year1 <- 2002 # from conditioning model
year <- year1 + seq(1, d$Ldyr) - 1
obs <- d$obsescape
esc <- exp(sapply(report, getElement, "logpredesc")) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975))
plot(NULL, NULL, type = "n", pch = 16, xlim = range(year), ylim = c(0, 1.1) * range(esc, obs, na.rm = TRUE), xlab = "", ylab = "")
polygon(c(rev(year), year), c(rev(esc[1, ]), esc[3, ]), col = alpha("grey", 0.5), border = NA)
lines(year, esc[2, ], lwd = 1.5)
points(year, obs, pch = 16)
lines(year, obs, lty = 3)
mtext("(c) Salmon River", side = 3, adj = 0, line = 0.5)

ERM <- readRDS("CM/Adam_09.29.26.base.rds")
report <- salmonMSE:::get_report(ERM)
d <- salmonMSE:::get_CMdata(ERM@.MISC$CMfit)
year1 <- 2002 # from conditioning model
year <- year1 + seq(1, d$Ldyr) - 1
obs <- d$obsescape
esc <- exp(sapply(report, getElement, "logpredesc")) %>%
  apply(1, quantile, probs = c(0.025, 0.5, 0.975))
plot(NULL, NULL, type = "n", pch = 16, xlim = range(year), ylim = c(0, 1.1) * range(esc, obs, na.rm = TRUE), xlab = "", ylab = "")
polygon(c(rev(year), year), c(rev(esc[1, ]), esc[3, ]), col = alpha("grey", 0.5), border = NA)
lines(year, esc[2, ], lwd = 1.5)
points(year, obs, pch = 16)
lines(year, obs, lty = 3)
mtext("(d) Adam/Eve Rivers", side = 3, adj = 0, line = 0.5)

mtext("Return Year", side = 1, outer = TRUE, line = -1)
mtext("Escapement", side = 2, outer = TRUE, line = -1)



dev.off()
