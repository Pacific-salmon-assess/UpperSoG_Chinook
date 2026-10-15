library(tidyverse)
library(readxl)

# Quinsam CWT recovery
rec <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-02-17-QuinsamChinook_Analyses_2005-2024.xlsx"),
  sheet = "Expanded"
) %>%
  mutate(is_catch = TotCatch > 0, is_esc = Escape > 0)

# CWT by release strategy
cwt_rs <- rec %>%
  mutate(RS = ifelse(RELEASE_STAGE_NAME == "Fed Fry", "Fed Fry", "Seapen/Smolt 0+")) %>%
  summarise(
    n_catch = sum(TotCatch),
    n_esc = sum(Escape),
    .by = c(Age, BROOD_YEAR, RS)
  )

g <- cwt_rs %>%
  filter(!is.na(Age), Age %in% seq(2, 5)) %>%
  mutate(Age = paste("Age", Age)) %>%
  ggplot(aes(BROOD_YEAR, n_catch, colour = RS)) +
  facet_grid(vars(RS), vars(Age), scales = "free_y") +
  geom_line() +
  geom_point() +
  labs(x = "Brood Year", y = "CWT catch") +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")
ggsave("figures/Quinsam_CWT_catch_rs.png", g, height = 3.5, width = 6)

g <- cwt_rs %>%
  filter(RS == "Seapen/Smolt 0+") %>%
  filter(!is.na(Age), Age %in% seq(2, 5)) %>%
  mutate(Age = paste("Age", Age)) %>%
  ggplot(aes(BROOD_YEAR, n_catch)) +
  facet_grid(vars(Age), scales = "free_y", rows=vars( )) +
  geom_line() +
  geom_point() +
  labs(x = "Brood Year", y = "CWT catch (trad)") +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")
ggsave("figures/Quinsam_CWT_catch.png", g, height = 3.5, width = 6)

g <- cwt_rs %>%
  filter(!is.na(Age), Age %in% seq(2, 5)) %>%
  mutate(Age = paste("Age", Age)) %>%
  ggplot(aes(BROOD_YEAR, n_esc)) +
  facet_grid(vars(RS), vars(Age), scales = "free_y") +
  geom_line() +
  geom_point() +
  labs(x = "Brood Year", y = "CWT escapement", colour = "Release Strategy") +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")
ggsave("figures/Quinsam_CWT_esc_rs.png", g, height = 3.5, width = 6)

g <- cwt_rs %>%
  filter(!is.na(Age), Age %in% seq(2, 5)) %>%
  mutate(Age = paste("Age", Age)) %>%
  filter(RS == "Seapen/Smolt 0+") %>%
  ggplot(aes(BROOD_YEAR, n_esc)) +
  facet_grid(vars(Age), scales = "free_y", rows=vars( )) +
  geom_line() +
  geom_point() +
  labs(x = "Brood Year", y = "CWT escapement (trad)") +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")
ggsave("figures/Quinsam_CWT_esc.png", g, height = 3.5, width = 6)

cwt_esc_annual <- cwt_rs %>%
  filter(RS == "Seapen/Smolt 0+") %>%
  summarise(n = sum(n_esc, na.rm=T), .by = BROOD_YEAR)
g <- cwt_esc_annual %>%
  ggplot(aes(BROOD_YEAR, n)) +
  geom_line() +
  geom_point() +
  labs(y = "CWT escapement (trad)")
ggsave("figures/Quinsam_CWT_esc_total.png", g, width = 4, height = 3)


# Full matrix of ages (1-5) and years (1979 - 2023)
full_matrix <- expand.grid(
  BroodYear = 2005:2021,
  Age = 1:5
) %>%
  as.data.frame()

g <- cwt_rs %>%
  filter(RS == "Seapen/Smolt 0+") %>%
  summarise(n = sum(n_esc), .by = c(BROOD_YEAR, Age)) %>%
  rename(BroodYear = BROOD_YEAR) %>%
  right_join(full_matrix, by = c("BroodYear", "Age")) %>%
  arrange(Age, BroodYear) %>%
  mutate(p = n/sum(n, na.rm = TRUE), .by = BroodYear) %>%
  filter(!is.na(p)) %>%
  ggplot(aes(BroodYear, p, fill = factor(Age, levels = 5:2))) +
  geom_col(width = 1, colour = "grey40") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Brood Year", y = "Proportion", fill = "Age", title = "CWT escapement (trad)") +
  coord_cartesian(expand = FALSE)
ggsave("figures/Quinsam_CWT_esc_prop.png", g, height = 4, width = 6)


# Quinsam - CWT releases
rel <- readxl::read_excel(
  file.path("data", "Quinsam", "2025-02-17-QuinsamChinook_Analyses_2005-2024.xlsx"),
  sheet = "Releases"
)

rel_rs <- rel %>%
  filter(RELEASE_STAGE_NAME %in% c("Fed Fry", "Smolt 0+", "Seapen 0+")) %>%
  mutate(RS = ifelse(RELEASE_STAGE_NAME == "Fed Fry", "Fed Fry", "Seapen/Smolt 0+")) %>%
  summarise(n_CWT = sum(TaggedNum) - sum(ShedTagNum), .by = c(BROOD_YEAR, RS))

g <- ggplot(rel_rs, aes(BROOD_YEAR, n_CWT)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(RS), ncol = 1, scales = "free_y") +
  labs(x = "Brood year", y = "CWT releases") +
  theme(legend.position = "bottom") +
  expand_limits(y = 0)
ggsave("figures/Quinsam_CWT_rel_bystrategy.png", g, height = 5, width = 4)

rel_trad <- rel %>%
  filter(RELEASE_STAGE_NAME %in% c("Smolt 0+", "Seapen 0+")) %>%
  mutate(RS = "Seapen/Smolt 0+") %>%
  summarise(n_CWT = sum(TaggedNum) - sum(ShedTagNum), .by = c(BROOD_YEAR, RS))

g <- ggplot(rel_trad, aes(BROOD_YEAR, n_CWT)) +
  geom_point() +
  geom_line() +
  labs(x = "Brood year", y = "CWT releases (trad)") +
  theme(legend.position = "bottom") +
  expand_limits(y = 0)
ggsave("figures/Quinsam_CWT_rel.png", g, height = 5, width = 4)



