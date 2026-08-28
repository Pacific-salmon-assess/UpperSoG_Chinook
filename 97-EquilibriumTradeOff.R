# Code to produce equilibrium trade-off plots

# Code adapted from https://github.com/SCA-stock-assess/WCVI-CN-ResDoc/blob/main/2.%20R%20code/Equilibrium%20trade-off%20analysis/equilibrium_trade-off.R
# Brown et al. (2026)
# https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-5/Fs70-5-2026-020-eng.pdf

# Packages ----------------------------------------------------------------

pkgs <- c("tidyverse", "here", "readxl", "janitor", "gsl")
#install.packages(pkgs)

library(here)
library(tidyverse)
library(readxl)
library(janitor)
library(gsl) # For Lambert's W function used to calculate Umsy
library(salmonMSE)

# Disable scientific notation in outputs
options(scipen = 999)


# Get data

ERM_Woss <- readRDS("CM/Woss_07.22.26.prior.rds")
report_Woss <- salmonMSE:::get_report(ERM_Woss)
ERM_Adam <- readRDS("CM/Adam_08.06.26.prior.rds")
report_Adam <- salmonMSE:::get_report(ERM_Adam)
ERM_Salmon <- readRDS("CM/Salmon_08.18.26.prior.rds")
report_Salmon <- salmonMSE:::get_report(ERM_Salmon)

# Define functions for calculating Umsy, Srep, and Heq, and extracting alpha
# and beta values from conditioning model output ----------------------------


# Umsy using Scheuerell (2016) explicit solution, where alpha=loga
Umsy <- function(alpha) {

  Umsy = 1 - gsl::lambert_W0(exp(1-alpha))
  return(Umsy)

}


# Seq - equilibrium population escapement at fixed harvest rate
Seq <- function(alpha, beta, U) {

  Seq = (alpha - (-log(1-U)))/beta

  return(Seq)
}


# Heq - equilibrium harvest
Heq <- function(alpha, beta, U) {

  Seq = (alpha - (-log(1-U)))/beta

  Heq = (Seq * exp(alpha-(beta*Seq))) - Seq

  return(Heq)
}

# Sgen using Scheuerell (2016) explicit solution for SMSY, where alpha = loga
# Taken from salmonMSE R package
Sgen <- function(alpha, beta){

  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - alpha) ) ) / beta

  a <- exp(alpha)

  return( -1 / beta * gsl::lambert_W0( -beta * sMSY / a ) )
}

# Ugen,  where alpha = loga
Ugen <- function(alpha, beta){
  sgen <- Sgen(alpha, beta)
  fgen <-  alpha - beta * sgen # slope of diagnoal line intersecting curve at sgen
  ugen <- 1 - exp(-fgen)
}


# Extract full life cycle alpha and beta values from conditioning model output
get_alpha_s <- function(report, samp){
  d <- salmonMSE:::get_CMdata(samp@.MISC$CMfit)
  alpha_s <- salmonMSE:::.CM_prod(report, d, mean_bio = TRUE) # Ricker alpha, per spawner
  return(as.vector(alpha_s))
}

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

# report <- report_Woss
# samp <- ERM_Woss
# salmonMSE:::get_CMdata(samp@.MISC$CMfit)
#
# summary(as.vector(salmonMSE:::.CM_SMSY(report,d, mean_bio=TRUE, type="spawner")))
# summary(as.vector(calc_Smsy_Ricker(log(salmonMSE:::.CM_prod(report, d, mean_bio = TRUE)), get_beta_s(report,samp))))

# Calculate equilibrium harvest curves -------------------------------------


# Define range of harvest rates
U <- tibble(U = seq(0, 0.99, by = 0.01))

# Parameter distributions for each population ----------------------------------

pop_params <- tribble(
  ~pop, ~a_dat, ~beta_dat,
  "Adam", get_alpha_s(report_Adam, ERM_Adam),  get_beta_s(report_Adam, ERM_Adam),
  "Salmon", get_alpha_s(report_Salmon, ERM_Salmon),  get_beta_s(report_Salmon, ERM_Salmon),
  "Nimpkish", get_alpha_s(report_Woss, ERM_Woss),  get_beta_s(report_Woss, ERM_Woss))



pop_params <-
  pop_params |>
  rowwise() |>
  mutate(
    data = list(
      tibble(
        log_a =log(a_dat),
        beta = beta_dat
      ) |>
        crossing(U) |>
        mutate(
          Seq = Seq(alpha = log_a, beta = beta, U = U),
          Heq = Heq(alpha = log_a, beta = beta, U = U),
          Sgen = Sgen(alpha = log_a, beta =beta)
        )
    )
  )


# Calculate CU-level Umsy
umsy_vals <- pop_params |>
  unnest(data) |>
  distinct(pop, log_a, beta) |>
  mutate(umsy = Umsy(log_a)) |>
  mutate(ugen = Ugen(log_a, beta))|>
  summarise(
    .by = pop,
    umsy_mid = median(umsy),
    umsy_lwr = quantile(umsy, 0.25),
    umsy_upr = quantile(umsy, 0.75),
    ugen_mid = median(ugen),
    ugen_lwr = quantile(ugen, 0.25),
    ugen_upr = quantile(ugen, 0.75)
  )


# Unpack and summarize data
eq_sum_data <- pop_params |>
  unnest(data) |>
  filter(!if_any(Seq:Heq, ~(is.na(.x) | is.infinite(.x)))) |>
  mutate(
    .by = pop,
    id = row_number(),
    #across(Seq:Heq, ~if_else(.x < 0, 0, .x))
  ) |>
  summarize(
    .by = c(id, U),
    across(Seq:Heq, sum)
  ) |>
  summarize(
    .by = U,
    Seq = median(Seq),
    Heq_mid = median(Heq),
    Heq_lwr = quantile(Heq, 0.25),
    Heq_upr = quantile(Heq, 0.75)
  ) |>
  mutate(
    umsy_mid = list(umsy_vals$umsy_mid),
    umsy_lwr = list(umsy_vals$umsy_lwr),
    umsy_upr = list(umsy_vals$umsy_upr),
    across(
      contains("umsy"),
      ~map2(
        .x = U,
        .y = .,
        ~ifelse(.x > .y, 1, 0)
      ) |>
        map(sum) |>
        unlist(),
      .names = "num_pop_exceed_{.col}"
    )
  ) |>
  mutate(
    ugen_mid = list(umsy_vals$ugen_mid),
    ugen_lwr = list(umsy_vals$ugen_lwr),
    ugen_upr = list(umsy_vals$ugen_upr),
    across(
      contains("ugen"),
      ~map2(
        .x = U,
        .y = .,
        ~ifelse(.x > .y, 1, 0)
      ) |>
        map(sum) |>
        unlist(),
      .names = "num_pop_exceed_{.col}"
    )
  )




# Calculate secondary y-axis transformation
ratio <- max(eq_sum_data$Heq_upr)*1.05/max(eq_sum_data$num_pop_exceed_umsy_upr)
ratio2 <- max(eq_sum_data$Heq_upr)*1.05/max(eq_sum_data$num_pop_exceed_ugen_upr)


# Calculate Smsy median, upper, and lower values
mid_Heq_Smsy = max(eq_sum_data$Heq_mid)
lwr_Heq_Smsy = max(eq_sum_data$Heq_lwr)
upr_Heq_Smsy = max(eq_sum_data$Heq_upr)


# label Smsy
smsy_lab <- tibble(
  Smsy_mid = eq_sum_data$Seq[eq_sum_data$Heq_mid == mid_Heq_Smsy],
  # Heq and Smsy lower and upper bounds are reversed
  Smsy_upr = eq_sum_data$Seq[eq_sum_data$Heq_lwr == lwr_Heq_Smsy],
  Smsy_lwr = eq_sum_data$Seq[eq_sum_data$Heq_upr == upr_Heq_Smsy],
  Heq_mid = mid_Heq_Smsy
) |>
  mutate(label = Smsy_mid)


# U labels at 5% increments
U_lab <- eq_sum_data |>
  mutate(U= as.character(U)) |>
  select(U, Seq) |>
  right_join(tibble(U = as.character(seq(0, 0.95, by = 0.05)))) |>
  distinct(Seq, .keep_all = TRUE) |>
  mutate(U = as.numeric(U))


# Plot summarized data
(eq_plot <- ggplot(eq_sum_data, aes(x = Seq, y = Heq_mid)) +
    # Add vertical lines showing harvest rate steps
    geom_vline(
      xintercept = U_lab$Seq,
      colour = "grey90"
    ) +
    # Label harvest rate steps
    annotate(
      "text",
      x = U_lab$Seq,
      y = upr_Heq_Smsy,
      label = scales::percent(U_lab$U),
      # Make text vertical and offset to the right of the lines
      angle = 270,
      vjust = -0.5,
      hjust = 0,
      colour = "grey75"
    ) +
    # Add stepped line showing # CUs where agg ER exceeds Umsy
    geom_step(
      aes(y = num_pop_exceed_umsy_mid*ratio),
      colour = "red",
      linewidth = 1.25
    ) +
    # # Add stepped line showing # CUs where agg ER exceeds Ugen
    # geom_step(
    #   aes(y = num_pop_exceed_ugen_mid*ratio),
    #   colour = "darkgreen",
    #   linewidth = 1.25
    # )+
    # Add stepped CI corresponding to the stepped line for UMSY
    geom_rect(
      aes(
        xmin = Seq,
        xmax = lead(Seq),
        ymin = num_pop_exceed_umsy_lwr*ratio,
        ymax = num_pop_exceed_umsy_upr*ratio
      ),
      fill = "red",
      alpha = 0.2
    ) +
    # # Add stepped CI corresponding to the stepped line, for Ugen
    # geom_rect(
    #   aes(
    #     xmin = Seq,
    #     xmax = lead(Seq),
    #     ymin = num_pop_exceed_ugen_lwr*ratio,
    #     ymax = num_pop_exceed_ugen_upr*ratio
    #   ),
    #   fill = "darkgreen",
    #   alpha = 0.2
    # ) +
    geom_line(
      colour = "blue",
      linewidth = 1
    ) +
    geom_ribbon(
      aes(ymin = Heq_lwr, ymax = Heq_upr),
      fill = "blue",
      alpha = 0.2
    ) +
    # Pointrange showing the estimated Smsy
    geom_pointrange(
      data = smsy_lab,
      aes(x = Smsy_mid, xmin = Smsy_lwr, xmax = Smsy_upr),
      size = 1,
      linewidth = 1
    ) +
    # Add text annotation that states midpoint and IQR of Smsy
    annotate(
      "text",
      x = smsy_lab$label,
      y = smsy_lab$Heq_mid,
      label = paste0(
        "S[MSY] ==~",
        round(smsy_lab$label, 0),
        "~(IQR:~",
        round(smsy_lab$Smsy_lwr),
        "-",
        round(smsy_lab$Smsy_upr),
        ")"
      ),
      hjust = 0.1,
      vjust = -0.4,#-1
      parse = TRUE
    ) +
    scale_x_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      labels = scales::comma,
      expand = expansion(mult = c(0, 0)),
      sec.axis = sec_axis(
        transform = ~./ratio,
        name = "Number of pops where agg. ER > Umsy"
      )
    ) +
    coord_cartesian(
      xlim = c(0, max(eq_sum_data$Seq)),
      ylim = c(0, max(eq_sum_data$Heq_upr*1.05))
    ) +
    labs(
      y = "Aggregate equilibrium harvest",
      x = "Aggregate equilibrium spawners"
    ) +
    theme_classic() +
    # Match axis title colour to corresponding line colour
    theme(
      axis.title.y.left = element_text(colour = "blue"),
      axis.title.y.right = element_text(colour = "red"),
    ) +
    #Hard coded order of populations
    annotate(geom="text", x=5000, y=0.92*ratio,
             label="Salmon",
             colour="red", size=3) +
    annotate(geom="text", x=5000, y=1.8*ratio,
             label="Adam", colour="red",
             size=3) +
    annotate(geom="text", x=3000, y=2.92*ratio,
             label="Nimpkish", colour="red",
             size=3)

)


# Save plot
ggsave(
  eq_plot,
  filename = here(
    "figures",
    "R-PLOT_equilibrium_harvest_curves.png"
    # "R-PLOT_equilibrium_harvest_curves_Ugen.png"
  ),
  height = 4.5,
  width = 7, units = "in"
)




# Compile and save table of key output values -----------------------------


# Smsy values are already saved in a dataframe for the trade-off plot
eq_outputs <- smsy_lab |>
  select(-label) |>
  pivot_longer(everything()) |>
  # Manually add Umsy values and Heq lwr/upr
  add_row(
    name = c(
      "Heq_lwr",
      "Heq_upr",
      "Umsy_mid",
      "Umsy_lwr",
      "Umsy_upr"
    ),
    value = c(
      lwr_Heq_Smsy,
      upr_Heq_Smsy,
      eq_sum_data$U[eq_sum_data$Heq_mid == mid_Heq_Smsy],
      eq_sum_data$U[eq_sum_data$Heq_lwr == lwr_Heq_Smsy],
      eq_sum_data$U[eq_sum_data$Heq_upr == upr_Heq_Smsy]
    )
  ) |>
  separate(
    name,
    c("variable", "increment"),
    sep = "_"
  ) |>
  pivot_wider(
    names_from = increment,
    values_from = value
  )


# Save the output table of key values
write.csv(
  eq_outputs,
  here(
    "data",
    "Equilibrium trade-off analysis",
    "R-OUT_SMU_ref-pt_values_eq-trade-off.csv"
  ),
  row.names = FALSE
)
