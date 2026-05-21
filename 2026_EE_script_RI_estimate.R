############################################################
## REPRODUCTIVE ISOLATION ESTIMATES 
############################################################

############################################################
## 0. PACKAGES AND GLOBAL SETTINGS
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(grid)
library(ggpubr)
library(patchwork)

set.seed(1830)
nboot_main <- 10000

# File paths
path_ri_model      <- "data/07a_2026_EE_vis_ri_calculation1.csv"
path_ri_strengths  <- "data/07b_2026_EE_vis_ri_calculation2.csv"
path_egg           <- "data/02_2026_EE_vis_egg.csv"
path_larval        <- "data/03_2026_EE_vis_larval_performance.csv"
path_fidelity      <- "data/04_2026_EE_data_fidelity.csv"
path_choice_mating <- "data/06_2026_EE_data_choice_mating.csv"

# Shared factor levels and colours
barrier_levels <- c(
  "seasonal", "habitat", "sexual", "prehatching", "csp",
  "hybrid inviability", "total"
)

species_levels <- c(
  "Hdk (M x L)", "Hni x Hya", "Hdk (M x D)",
  "Hvm x Hpt", "Hse x Hpl", "Hvp x Hen"
)

barrier_cols <- c(
  "seasonal" = "#b3de69",
  "habitat" = "#8dd3c7",
  "sexual" = "#fccde5",
  "prehatching" = "#fdb462",
  "csp" = "#fb8072",
  "hybrid inviability" = "#80b1d3",
  "total" = "black"
)

############################################################
## 1. HELPER FUNCTIONS
############################################################

# General RI formula
ri_from_C_H <- function(C, H) {
  1 - 2 * (H / (H + C))
}

# RI from binary data: 1 = heterospecific, 0 = conspecific
calc_RI_binary <- function(x) {
  H <- sum(x == 1)
  C <- sum(x == 0)
  1 - 2 * (H / (H + C))
}

bootstrap_RI_binary <- function(x, nboot = 10000, seed = 1830) {
  set.seed(seed)
  boot_vals <- replicate(nboot, {
    samp <- sample(x, size = length(x), replace = TRUE)
    calc_RI_binary(samp)
  })
  list(
    point = calc_RI_binary(x),
    mean = mean(boot_vals),
    se = sd(boot_vals),
    ci = quantile(boot_vals, c(0.025, 0.975)),
    boot = boot_vals
  )
}

# Exact binomial-based RI, used for habitat isolation
ri_from_binom <- function(H, N) {
  bt <- binom.test(H, N)
  p_hat <- H / N
  RI_hat <- 1 - 2 * p_hat

  p_low <- bt$conf.int[1]
  p_high <- bt$conf.int[2]

  RI_low <- 1 - 2 * p_high
  RI_high <- 1 - 2 * p_low

  list(
    p_hat = p_hat,
    p_CI = c(p_low, p_high),
    RI_hat = RI_hat,
    RI_CI = c(RI_low, RI_high)
  )
}

# Approximate bootstrap for RI from published probabilities
bootstrap_RI_from_probs <- function(pC, pH, nC, nH, nboot = 10000, seed = 123) {
  set.seed(seed)

  boot_vals <- replicate(nboot, {
    C_boot <- rbinom(1, size = nC, prob = pC) / nC
    H_boot <- rbinom(1, size = nH, prob = pH) / nH
    ri_from_C_H(C_boot, H_boot)
  })

  list(
    point = ri_from_C_H(pC, pH),
    mean = mean(boot_vals),
    se = sd(boot_vals),
    ci = quantile(boot_vals, c(0.025, 0.975)),
    boot = boot_vals
  )
}

# Approximate RI bootstrap from published summary data
# Use when only proportions and sample sizes are available.
# This treats successes as binomial outcomes and propagates uncertainty
# from the published sample size and rate.
bootstrap_RI_from_published_binomial <- function(pC, pH, nC, nH,
                                                 nboot = 10000,
                                                 seed = 1830) {
  bootstrap_RI_from_probs(
    pC = pC,
    pH = pH,
    nC = nC,
    nH = nH,
    nboot = nboot,
    seed = seed
  )
}

# Sequential absolute RI contribution model
compute_ri <- function(df) {
  df <- df %>% arrange(order)

  e_values <- matrix(0, nrow = nrow(df), ncol = 1001)
  colnames(e_values) <- sprintf("e%.3f", seq(0, 1, by = 0.001))

  for (i in seq_along(colnames(e_values))) {
    e_col <- colnames(e_values)[i]
    e_index <- as.numeric(gsub("e", "", e_col))

    for (j in 1:nrow(df)) {
      if (df$barrier[j] == "habitat") {
        e_values[j, i] <- df$strength[j] * e_index
      } else {
        e_values[j, i] <- df$strength[j] * (1 - sum(e_values[1:j, i]))
      }
    }
  }

  df <- cbind(df, e_values)
  return(df)
}

# Weighted and unweighted egg hatchability helpers
egg_weighted_hatch <- function(df) {
  sum(df$htc) / sum(df$prod)
}

calc_RI_egg_weighted <- function(data, conspecific_pair, hetero_pair) {
  C_df <- data %>% filter(pair == conspecific_pair)
  H_df <- data %>% filter(pair == hetero_pair)

  if (nrow(C_df) == 0 || nrow(H_df) == 0) {
    stop("Missing conspecific or heterospecific rows for egg data.")
  }

  C <- egg_weighted_hatch(C_df)
  H <- egg_weighted_hatch(H_df)
  RI <- ri_from_C_H(C, H)

  list(
    C = C,
    H = H,
    RI = RI,
    n_conspecific = nrow(C_df),
    n_heterospecific = nrow(H_df)
  )
}

bootstrap_RI_egg_weighted <- function(data, conspecific_pair, hetero_pair,
                                      nboot = 10000, seed = 1830) {
  set.seed(seed)

  C_df <- data %>% filter(pair == conspecific_pair)
  H_df <- data %>% filter(pair == hetero_pair)

  if (nrow(C_df) == 0 || nrow(H_df) == 0) {
    stop("Missing conspecific or heterospecific rows for egg data.")
  }

  boot_vals <- replicate(nboot, {
    C_boot <- C_df %>% slice_sample(n = nrow(C_df), replace = TRUE)
    H_boot <- H_df %>% slice_sample(n = nrow(H_df), replace = TRUE)

    C <- egg_weighted_hatch(C_boot)
    H <- egg_weighted_hatch(H_boot)

    ri_from_C_H(C, H)
  })

  list(
    point = calc_RI_egg_weighted(data, conspecific_pair, hetero_pair)$RI,
    mean = mean(boot_vals),
    se = sd(boot_vals),
    ci = quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE),
    boot = boot_vals
  )
}

calc_RI_egg_unweighted <- function(data, conspecific_pair, hetero_pair) {
  C_df <- data %>% filter(pair == conspecific_pair)
  H_df <- data %>% filter(pair == hetero_pair)

  if (nrow(C_df) == 0 || nrow(H_df) == 0) {
    stop("Missing conspecific or heterospecific rows for egg data.")
  }

  C <- mean(C_df$hatch_rate, na.rm = TRUE)
  H <- mean(H_df$hatch_rate, na.rm = TRUE)
  RI <- ri_from_C_H(C, H)

  list(
    C = C,
    H = H,
    RI = RI,
    n_conspecific = nrow(C_df),
    n_heterospecific = nrow(H_df)
  )
}

bootstrap_RI_egg_unweighted <- function(data, conspecific_pair, hetero_pair,
                                        nboot = 10000, seed = 1830) {
  set.seed(seed)

  C_df <- data %>% filter(pair == conspecific_pair)
  H_df <- data %>% filter(pair == hetero_pair)

  if (nrow(C_df) == 0 || nrow(H_df) == 0) {
    stop("Missing conspecific or heterospecific rows for egg data.")
  }

  boot_vals <- replicate(nboot, {
    C_boot <- C_df %>% slice_sample(n = nrow(C_df), replace = TRUE)
    H_boot <- H_df %>% slice_sample(n = nrow(H_df), replace = TRUE)

    C <- mean(C_boot$hatch_rate, na.rm = TRUE)
    H <- mean(H_boot$hatch_rate, na.rm = TRUE)

    ri_from_C_H(C, H)
  })

  list(
    point = calc_RI_egg_unweighted(data, conspecific_pair, hetero_pair)$RI,
    mean = mean(boot_vals),
    se = sd(boot_vals),
    ci = quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE),
    boot = boot_vals
  )
}

# Larval performance helpers
larv_weighted_survival <- function(df, succ_col, tot_col) {
  sum(df[[succ_col]], na.rm = TRUE) / sum(df[[tot_col]], na.rm = TRUE)
}

calc_RI_larval_weighted <- function(data, conspecific_pair, hetero_pair, succ_col, tot_col) {
  C_df <- data %>% filter(pair == conspecific_pair)
  H_df <- data %>% filter(pair == hetero_pair)

  if (nrow(C_df) == 0 || nrow(H_df) == 0) {
    stop("Missing conspecific or heterospecific rows for larval data.")
  }

  C <- larv_weighted_survival(C_df, succ_col, tot_col)
  H <- larv_weighted_survival(H_df, succ_col, tot_col)
  RI <- ri_from_C_H(C, H)

  list(
    C = C,
    H = H,
    RI = RI,
    n_conspecific = nrow(C_df),
    n_heterospecific = nrow(H_df)
  )
}

bootstrap_RI_larval_weighted <- function(data, conspecific_pair, hetero_pair,
                                         succ_col, tot_col,
                                         nboot = 10000, seed = 1830) {
  set.seed(seed)

  C_df <- data %>% filter(pair == conspecific_pair)
  H_df <- data %>% filter(pair == hetero_pair)

  if (nrow(C_df) == 0 || nrow(H_df) == 0) {
    stop("Missing conspecific or heterospecific rows for larval data.")
  }

  boot_vals <- replicate(nboot, {
    C_boot <- C_df %>% slice_sample(n = nrow(C_df), replace = TRUE)
    H_boot <- H_df %>% slice_sample(n = nrow(H_df), replace = TRUE)

    C <- larv_weighted_survival(C_boot, succ_col, tot_col)
    H <- larv_weighted_survival(H_boot, succ_col, tot_col)

    ri_from_C_H(C, H)
  })

  list(
    point = calc_RI_larval_weighted(data, conspecific_pair, hetero_pair, succ_col, tot_col)$RI,
    mean = mean(boot_vals, na.rm = TRUE),
    se = sd(boot_vals, na.rm = TRUE),
    ci = quantile(boot_vals, c(0.025, 0.975), na.rm = TRUE),
    boot = boot_vals
  )
}

# Replace RI point and CI values for one barrier in a plotting dataframe
update_barrier <- function(df, barrier_name, x, y, xlow, xhigh, ylow, yhigh) {
  df %>%
    mutate(
      strength.1 = ifelse(barrier == barrier_name, x, strength.1),
      strength.2 = ifelse(barrier == barrier_name, y, strength.2),
      x_low      = ifelse(barrier == barrier_name, xlow, x_low),
      x_high     = ifelse(barrier == barrier_name, xhigh, x_high),
      y_low      = ifelse(barrier == barrier_name, ylow, y_low),
      y_high     = ifelse(barrier == barrier_name, yhigh, y_high)
    )
}

# Print a compact RI summary line
print_ri_line <- function(label, point, low = NA, high = NA, digits = 6) {
  if (is.na(low) || is.na(high)) {
    cat(label, " point:", round(point, digits), "\n")
  } else {
    cat(label, " point:", round(point, digits), " CI:", round(low, digits), round(high, digits), "\n")
  }
}

# Convert one comparison into a long-format RI table

make_ri_estimate_table <- function(lineage, axis_1, axis_2,
                                   barriers,
                                   strength_1, strength_2,
                                   low_1, high_1,
                                   low_2, high_2,
                                   method, notes = NA_character_) {
  tibble(
    lineage = lineage,
    barrier = barriers,
    axis_1 = axis_1,
    RI_axis_1 = strength_1,
    CI_low_axis_1 = low_1,
    CI_high_axis_1 = high_1,
    axis_2 = axis_2,
    RI_axis_2 = strength_2,
    CI_low_axis_2 = low_2,
    CI_high_axis_2 = high_2,
    method = method,
    notes = notes
  )
}

# Bootstrap RI from aggregated conspecific and heterospecific counts.
bootstrap_RI_from_counts <- function(C, H, nboot = 10000, seed = 1830) {
  set.seed(seed)

  x <- c(rep(0, C), rep(1, H))

  if (length(x) == 0) {
    return(list(
      point = NA_real_,
      mean = NA_real_,
      se = NA_real_,
      ci = c(`2.5%` = NA_real_, `97.5%` = NA_real_),
      boot = rep(NA_real_, nboot)
    ))
  }

  bootstrap_RI_binary(x, nboot = nboot, seed = seed)
}

# Bootstrap realised habitat RI from host-sighting counts.
# Net/off-host sightings are excluded from RI but retained in the table.
bootstrap_realised_habitat_counts <- function(C, H, nboot = 10000, seed = 1830) {
  bootstrap_RI_from_counts(C = C, H = H, nboot = nboot, seed = seed)
}

# Bootstrap realised sexual RI from mating-pair counts.
# Null/no-mating observations are excluded from RI but retained for mating rate.
bootstrap_realised_sexual_counts <- function(par, hyb, nboot = 10000, seed = 1830) {
  bootstrap_RI_from_counts(C = par, H = hyb, nboot = nboot, seed = seed)
}

############################################################
## 2. PLOT HELPER FUNCTIONS
############################################################

make_ri_model_plot <- function(data, species_name, linetype_values, linetype_labels) {
  ggplot(data[data$species == species_name, ], aes(x = e_value, y = RI, colour = barrier)) +
    geom_line(size = 0.5, alpha = 0.7, aes(linetype = pair)) +
    scale_colour_manual(values = barrier_cols) +
    scale_linetype_manual(values = linetype_values, labels = linetype_labels) +
    scale_y_continuous(
      breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1),
      labels = c("", "0", "0.25", "0.5", "0.75", "1")
    ) +
    scale_x_continuous(
      expand = c(0.01, 0.01),
      trans = "reverse",
      breaks = c(1, 0.75, 0.5, 0.25, 0),
      labels = c("1", "0.75", "0.5", "0.25", "0")
    ) +
    xlab("\nEnvironmental stability") +
    ylab(expression(atop("\n", italic(RI) ~ " absolute contribution"))) +
    theme_bw() +
    theme(
      legend.position = c(0.70, 0.95),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(face = "italic", size = 6),
      legend.key.size = unit(0.01, "cm"),
      legend.box.spacing = unit(0.5, "cm"),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      panel.spacing = unit(1, "lines"),
      panel.grid = element_blank()
    ) +
    guides(
      colour = "none",
      fill = "none",
      linetype = guide_legend(keywidth = unit(0.6, "cm"))
    )
}

make_ri_strength_plot <- function(plot_df, x_lab, y_lab, title_expr, show_legend = TRUE) {
  legend_theme <- if (show_legend) {
    theme(
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.01, "cm"),
      legend.box.spacing = unit(0.5, "cm")
    )
  } else {
    theme(legend.position = "none")
  }

  ggplot(plot_df, aes(x = strength.1, y = strength.2, colour = barrier)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_segment(
      aes(x = strength.1, xend = strength.1, y = y_low, yend = y_high, color = barrier),
      linewidth = 1
    ) +
    geom_segment(
      aes(y = strength.2, yend = strength.2, x = x_low, xend = x_high, color = barrier),
      linewidth = 1
    ) +
    geom_point(aes(fill = barrier), size = 4, shape = 21, stroke = 0.7, colour = "black") +
    scale_colour_manual(values = barrier_cols, drop = FALSE) +
    scale_fill_manual(values = barrier_cols, drop = FALSE) +
    coord_cartesian(xlim = c(-0.1, 1), ylim = c(-0.1, 1)) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    xlab(x_lab) +
    ylab(y_lab) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 16),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    legend_theme +
    guides(
      colour = guide_legend(override.aes = list(size = 2)),
      fill = guide_legend(override.aes = list(size = 2))
    ) +
    ggtitle(title_expr)
}

############################################################
## 3. ABSOLUTE RI CONTRIBUTION MODEL
############################################################

ri_data <- read.csv(path_ri_model)
ri_data$barrier <- factor(ri_data$barrier, levels = barrier_levels[barrier_levels != "total"])

ri_mod1 <- ri_data %>%
  filter(!is.na(order)) %>%
  group_by(species, pair) %>%
  do(compute_ri(.))

ri_mod2 <- ri_mod1 %>%
  filter(!is.na(order)) %>%
  group_by(species, pair) %>%
  summarise(across(starts_with("e"), sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(barrier = "total")

ri_mod <- bind_rows(ri_mod1, ri_mod2) %>%
  filter(species %in% c("Hdk (M x D)", "Hdk (M x L)", "Hvm x Hpt")) %>%
  pivot_longer(cols = starts_with("e"), names_to = "e_value", values_to = "RI") %>%
  mutate(e_value = as.numeric(sub("e", "", e_value)))

ri_mod$barrier <- factor(ri_mod$barrier, levels = barrier_levels)
ri_mod$species <- factor(ri_mod$species, levels = species_levels)
levels(ri_mod$species) <- gsub(" x ", " × ", levels(ri_mod$species))

plot.ri.model3 <- make_ri_model_plot(
  data = ri_mod,
  species_name = "Hdk (M × D)",
  linetype_values = c(
    "H. diekei D-race" = "solid",
    "H. diekei M-race" = "dashed"
  ),
  linetype_labels = c(
    expression(italic("H. diekei") ~ italic(D) ~ "-race"),
    expression(italic("H. diekei") ~ italic(M) ~ "-race")
  )
)

plot.ri.model4 <- make_ri_model_plot(
  data = ri_mod,
  species_name = "Hdk (M × L)",
  linetype_values = c(
    "H. diekei L-race" = "solid",
    "H. diekei M-race" = "dashed"
  ),
  linetype_labels = c(
    expression(italic("H. diekei") ~ italic(L) ~ "-race"),
    expression(italic("H. diekei") ~ italic(M) ~ "-race")
  )
)

############################################################
## 4. H. diekei M x D: RI ESTIMATES
############################################################

# Sexual isolation
sex_data_M <- c(rep(1, 4), rep(0, 18))
sex_data_D <- c(rep(1, 16), rep(0, 24))

sex_M_res <- bootstrap_RI_binary(sex_data_M, nboot = nboot_main, seed = 1830)
sex_D_res <- bootstrap_RI_binary(sex_data_D, nboot = nboot_main, seed = 1830)

# Habitat isolation
hab_M_res <- ri_from_binom(H = 0, N = 37)
hab_D_res <- ri_from_binom(H = 0, N = 92)

# Egg hatchability
egg2 <- read.csv(path_egg, na.strings = "na")

egg_clean <- egg2 %>%
  mutate(
    htc = ifelse(htc %in% c("na", "NA", ""), NA, htc),
    htc = as.numeric(htc),
    prod = as.numeric(prod),
    value = as.numeric(value),
    fem = as.character(fem),
    mle = as.character(mle),
    pair = as.character(pair)
  ) %>%
  filter(!is.na(htc), !is.na(prod), !is.na(value), value > 0, prod > 0)

egg_expanded <- egg_clean %>%
  uncount(weights = value) %>%
  mutate(hatch_rate = htc / prod)

egg_data_D <- egg_expanded %>% filter(fem == "drace")
egg_data_M <- egg_expanded %>% filter(fem == "mrace")

egg_D_point_weighted <- calc_RI_egg_weighted(egg_data_D, "dd", "dm")
egg_M_point_weighted <- calc_RI_egg_weighted(egg_data_M, "mm", "md")

egg_D_boot_weighted <- bootstrap_RI_egg_weighted(egg_data_D, "dd", "dm", nboot = nboot_main, seed = 1830)
egg_M_boot_weighted <- bootstrap_RI_egg_weighted(egg_data_M, "mm", "md", nboot = nboot_main, seed = 1830)

egg_D_point_unweighted <- calc_RI_egg_unweighted(egg_data_D, "dd", "dm")
egg_M_point_unweighted <- calc_RI_egg_unweighted(egg_data_M, "mm", "md")

egg_D_boot_unweighted <- bootstrap_RI_egg_unweighted(egg_data_D, "dd", "dm", nboot = nboot_main, seed = 1830)
egg_M_boot_unweighted <- bootstrap_RI_egg_unweighted(egg_data_M, "mm", "md", nboot = nboot_main, seed = 1830)

# Larval performance
lp2 <- read.csv(path_larval)

lp2$pair <- factor(lp2$pair, levels = c("mm", "md", "dd", "dm"))
lp2$stg <- factor(lp2$stg, levels = c("Larval acceptance", "Survival to the second instar", "Reaching adulthood"))

larv_clean <- lp2 %>%
  rename(
    succ_mik  = succ.mik,
    fail_mik  = fail.mik,
    tot_mik   = tot.mik,
    rsucc_mik = rsucc.mik,
    rfail_mik = rfail.mik,
    succ_dic  = succ.dic,
    fail_dic  = fail.dic,
    tot_dic   = tot.dic,
    rsucc_dic = rsucc.dic,
    rfail_dic = rfail.dic,
    stg       = stg
  ) %>%
  mutate(
    pair = as.character(pair),
    fem  = as.character(fem),
    mle  = as.character(mle),
    fam  = as.character(fam),
    stg  = as.character(stg)
  )

larv_adult <- larv_clean %>%
  filter(stg == "Reaching adulthood")

larv_data_M <- larv_adult %>%
  filter(fem == "mrace", pair %in% c("mm", "md"))

larv_data_D <- larv_adult %>%
  filter(fem == "drace", pair %in% c("dd", "dm"))

larv_M_point_weighted <- calc_RI_larval_weighted(
  data = larv_data_M,
  conspecific_pair = "mm",
  hetero_pair = "md",
  succ_col = "succ_mik",
  tot_col = "tot_mik"
)

larv_D_point_weighted <- calc_RI_larval_weighted(
  data = larv_data_D,
  conspecific_pair = "dd",
  hetero_pair = "dm",
  succ_col = "succ_dic",
  tot_col = "tot_dic"
)

larv_M_boot_weighted <- bootstrap_RI_larval_weighted(
  data = larv_data_M,
  conspecific_pair = "mm",
  hetero_pair = "md",
  succ_col = "succ_mik",
  tot_col = "tot_mik",
  nboot = nboot_main,
  seed = 1830
)

larv_D_boot_weighted <- bootstrap_RI_larval_weighted(
  data = larv_data_D,
  conspecific_pair = "dd",
  hetero_pair = "dm",
  succ_col = "succ_dic",
  tot_col = "tot_dic",
  nboot = nboot_main,
  seed = 1830
)

############################################################
## 5. OPTIONAL PRINT CHECKS: H. diekei M x D
############################################################

cat("\n===== SEXUAL =====\n")
print_ri_line("M", sex_M_res$point, sex_M_res$ci[1], sex_M_res$ci[2])
print_ri_line("D", sex_D_res$point, sex_D_res$ci[1], sex_D_res$ci[2])

cat("\n===== HABITAT =====\n")
print_ri_line("M", hab_M_res$RI_hat, hab_M_res$RI_CI[1], hab_M_res$RI_CI[2])
print_ri_line("D", hab_D_res$RI_hat, hab_D_res$RI_CI[1], hab_D_res$RI_CI[2])

cat("\n===== EGG (WEIGHTED) =====\n")
print_ri_line("M", egg_M_point_weighted$RI, egg_M_boot_weighted$ci[1], egg_M_boot_weighted$ci[2])
print_ri_line("D", egg_D_point_weighted$RI, egg_D_boot_weighted$ci[1], egg_D_boot_weighted$ci[2])

cat("\n===== EGG (UNWEIGHTED) =====\n")
print_ri_line("M", egg_M_point_unweighted$RI, egg_M_boot_unweighted$ci[1], egg_M_boot_unweighted$ci[2])
print_ri_line("D", egg_D_point_unweighted$RI, egg_D_boot_unweighted$ci[1], egg_D_boot_unweighted$ci[2])

cat("\n===== LARVAL PERFORMANCE =====\n")
print_ri_line("M", larv_M_point_weighted$RI, larv_M_boot_weighted$ci[1], larv_M_boot_weighted$ci[2])
print_ri_line("D", larv_D_point_weighted$RI, larv_D_boot_weighted$ci[1], larv_D_boot_weighted$ci[2])

############################################################
## 6. H. diekei M x D: BUILD PLOTTING DATA
############################################################

ri_data2 <- read.csv(path_ri_strengths)
ri_data2$barrier <- factor(ri_data2$barrier, levels = barrier_levels)
ri_data2$species <- factor(ri_data2$species)
levels(ri_data2$species) <- gsub(" x ", " × ", levels(ri_data2$species))

plot_df <- ri_data2 %>%
  filter(species == "Hdk (M × D)") %>%
  mutate(
    x_low  = strength.1,
    x_high = strength.1,
    y_low  = strength.2,
    y_high = strength.2
  ) %>%
  update_barrier(
    barrier_name = "habitat",
    x = hab_M_res$RI_hat,
    y = hab_D_res$RI_hat,
    xlow = hab_M_res$RI_CI[1],
    xhigh = hab_M_res$RI_CI[2],
    ylow = hab_D_res$RI_CI[1],
    yhigh = hab_D_res$RI_CI[2]
  ) %>%
  update_barrier(
    barrier_name = "sexual",
    x = sex_M_res$point,
    y = sex_D_res$point,
    xlow = sex_M_res$ci[1],
    xhigh = sex_M_res$ci[2],
    ylow = sex_D_res$ci[1],
    yhigh = sex_D_res$ci[2]
  ) %>%
  update_barrier(
    barrier_name = "prehatching",
    x = egg_M_point_weighted$RI,
    y = egg_D_point_weighted$RI,
    xlow = egg_M_boot_weighted$ci[1],
    xhigh = egg_M_boot_weighted$ci[2],
    ylow = egg_D_boot_weighted$ci[1],
    yhigh = egg_D_boot_weighted$ci[2]
  ) %>%
  update_barrier(
    barrier_name = "hybrid inviability",
    x = larv_M_point_weighted$RI,
    y = larv_D_point_weighted$RI,
    xlow = larv_M_boot_weighted$ci[1],
    xhigh = larv_M_boot_weighted$ci[2],
    ylow = larv_D_boot_weighted$ci[1],
    yhigh = larv_D_boot_weighted$ci[2]
  )

plot_ri2 <- make_ri_strength_plot(
  plot_df = plot_df,
  x_lab = expression(atop("\n", italic(RI) ~ "strength" ~ italic(M) ~ "-race")),
  y_lab = expression(atop("\n", italic(RI) ~ "strength" ~ italic(D) ~ "-race")),
  title_expr = expression(italic(F[ST]) == 1.0000 ~ ";" ~ italic(D[XY]) == 0.0079),
  show_legend = TRUE
)

############################################################
## 7. H. diekei M x L: PUBLISHED RI ESTIMATES
############################################################

# Habitat isolation
hab_ML_M_hab <- ri_from_binom(H = 0, N = 86)
hab_ML_L_hab <- ri_from_binom(H = 0, N = 68)

# Sexual isolation: approximate bootstrap from published proportions
sex_ML_M_sex <- bootstrap_RI_from_probs(
  pC = 0.63, pH = 0.74, nC = 59, nH = 64,
  nboot = nboot_main, seed = 1830
)

sex_ML_L_sex <- bootstrap_RI_from_probs(
  pC = 0.79, pH = 0.76, nC = 64, nH = 59,
  nboot = nboot_main, seed = 1830
)

# Egg hatchability: approximate bootstrap from published Table 5
# The SD in Table 5 refers to number of eggs per batch, not hatchability.
# Therefore, uncertainty is approximated using total number of eggs and
# published hatching rate as a binomial outcome.
egg_ML_published <- tibble(
  female_race = c("M-race", "M-race", "L-race", "L-race"),
  male_race   = c("M-race", "L-race", "L-race", "M-race"),
  pair        = c("MM", "ML", "LL", "LM"),
  n_females   = c(11, 10, 13, 10),
  n_batches   = c(66, 60, 73, 56),
  n_eggs      = c(1302, 1130, 1075, 978),
  eggs_per_batch_mean = c(20.2, 17.6, 19.0, 17.9),
  eggs_per_batch_sd   = c(2.7, 3.9, 2.0, 3.7),
  hatch_rate  = c(0.723, 0.702, 0.643, 0.678)
) %>%
  mutate(
    n_hatched = round(n_eggs * hatch_rate),
    n_unhatched = n_eggs - n_hatched
  )

# M female: conspecific = MM, heterospecific = ML
egg_ML_M_boot <- bootstrap_RI_from_published_binomial(
  pC = egg_ML_published$hatch_rate[egg_ML_published$pair == "MM"],
  pH = egg_ML_published$hatch_rate[egg_ML_published$pair == "ML"],
  nC = egg_ML_published$n_eggs[egg_ML_published$pair == "MM"],
  nH = egg_ML_published$n_eggs[egg_ML_published$pair == "ML"],
  nboot = nboot_main,
  seed = 1830
)

# L female: conspecific = LL, heterospecific = LM
egg_ML_L_boot <- bootstrap_RI_from_published_binomial(
  pC = egg_ML_published$hatch_rate[egg_ML_published$pair == "LL"],
  pH = egg_ML_published$hatch_rate[egg_ML_published$pair == "LM"],
  nC = egg_ML_published$n_eggs[egg_ML_published$pair == "LL"],
  nH = egg_ML_published$n_eggs[egg_ML_published$pair == "LM"],
  nboot = nboot_main,
  seed = 1830
)

egg_ML_M_egg_RI <- egg_ML_M_boot$point
egg_ML_L_egg_RI <- egg_ML_L_boot$point

# Hybrid inviability / larval performance: approximate bootstrap from published Table 6
# Here we use survivorship to adulthood on the maternal host plant,
# matching the original RI calculation:
# M-race axis: C = Mikania race on Mikania, H = F1ML on Mikania
# L-race axis: C = Leucas race on Leucas, H = F1LM on Leucas
larv_ML_published <- tibble(
  comparison = c("M maternal host", "M maternal host", "L maternal host", "L maternal host"),
  larvae_type = c("Mikania race", "F1ML", "Leucas race", "F1LM"),
  role = c("C", "H", "C", "H"),
  host = c("Mikania", "Mikania", "Leucas", "Leucas"),
  n_families = c(10, 6, 10, 7),
  n_larvae = c(100, 60, 100, 70),
  acceptance_first_instar = c(1.00, 0.93, 0.96, 0.97),
  survival_second_instar = c(0.92, 0.77, 0.90, 0.90),
  adult_survival = c(0.64, 0.65, 0.81, 0.79)
) %>%
  mutate(
    n_adult = round(n_larvae * adult_survival),
    n_not_adult = n_larvae - n_adult
  )

# M maternal host: C = Mikania race, H = F1ML
larv_ML_M_boot <- bootstrap_RI_from_published_binomial(
  pC = larv_ML_published$adult_survival[
    larv_ML_published$comparison == "M maternal host" &
      larv_ML_published$role == "C"
  ],
  pH = larv_ML_published$adult_survival[
    larv_ML_published$comparison == "M maternal host" &
      larv_ML_published$role == "H"
  ],
  nC = larv_ML_published$n_larvae[
    larv_ML_published$comparison == "M maternal host" &
      larv_ML_published$role == "C"
  ],
  nH = larv_ML_published$n_larvae[
    larv_ML_published$comparison == "M maternal host" &
      larv_ML_published$role == "H"
  ],
  nboot = nboot_main,
  seed = 1830
)

# L maternal host: C = Leucas race, H = F1LM
larv_ML_L_boot <- bootstrap_RI_from_published_binomial(
  pC = larv_ML_published$adult_survival[
    larv_ML_published$comparison == "L maternal host" &
      larv_ML_published$role == "C"
  ],
  pH = larv_ML_published$adult_survival[
    larv_ML_published$comparison == "L maternal host" &
      larv_ML_published$role == "H"
  ],
  nC = larv_ML_published$n_larvae[
    larv_ML_published$comparison == "L maternal host" &
      larv_ML_published$role == "C"
  ],
  nH = larv_ML_published$n_larvae[
    larv_ML_published$comparison == "L maternal host" &
      larv_ML_published$role == "H"
  ],
  nboot = nboot_main,
  seed = 1830
)

larv_ML_M_larv_RI <- larv_ML_M_boot$point
larv_ML_L_larv_RI <- larv_ML_L_boot$point

cat("\n===== H. diekei M x L: HABITAT ISOLATION =====\n")
print_ri_line("Mikania race", hab_ML_M_hab$RI_hat, hab_ML_M_hab$RI_CI[1], hab_ML_M_hab$RI_CI[2])
print_ri_line("Leucas race", hab_ML_L_hab$RI_hat, hab_ML_L_hab$RI_CI[1], hab_ML_L_hab$RI_CI[2])

cat("\n===== H. diekei M x L: SEXUAL ISOLATION =====\n")
print_ri_line("Mikania race", sex_ML_M_sex$point, sex_ML_M_sex$ci[1], sex_ML_M_sex$ci[2])
print_ri_line("Leucas race", sex_ML_L_sex$point, sex_ML_L_sex$ci[1], sex_ML_L_sex$ci[2])

cat("\n===== H. diekei M x L: EGG HATCHABILITY =====\n")
print_ri_line("Mikania race", egg_ML_M_boot$point, egg_ML_M_boot$ci[1], egg_ML_M_boot$ci[2])
print_ri_line("Leucas race", egg_ML_L_boot$point, egg_ML_L_boot$ci[1], egg_ML_L_boot$ci[2])

cat("\n===== H. diekei M x L: HYBRID INVIABILITY / LARVAL PERFORMANCE =====\n")
print_ri_line("Mikania race", larv_ML_M_boot$point, larv_ML_M_boot$ci[1], larv_ML_M_boot$ci[2])
print_ri_line("Leucas race", larv_ML_L_boot$point, larv_ML_L_boot$ci[1], larv_ML_L_boot$ci[2])

ri_ML_results <- tibble(
  lineage = "H. diekei (Mikania x Leucas)",
  barrier = c("habitat", "sexual", "prehatching", "hybrid inviability"),
  strength_M = c(
    hab_ML_M_hab$RI_hat,
    sex_ML_M_sex$point,
    egg_ML_M_boot$point,
    larv_ML_M_boot$point
  ),
  strength_L = c(
    hab_ML_L_hab$RI_hat,
    sex_ML_L_sex$point,
    egg_ML_L_boot$point,
    larv_ML_L_boot$point
  ),
  x_low = c(
    hab_ML_M_hab$RI_CI[1],
    sex_ML_M_sex$ci[1],
    egg_ML_M_boot$ci[1],
    larv_ML_M_boot$ci[1]
  ),
  x_high = c(
    hab_ML_M_hab$RI_CI[2],
    sex_ML_M_sex$ci[2],
    egg_ML_M_boot$ci[2],
    larv_ML_M_boot$ci[2]
  ),
  y_low = c(
    hab_ML_L_hab$RI_CI[1],
    sex_ML_L_sex$ci[1],
    egg_ML_L_boot$ci[1],
    larv_ML_L_boot$ci[1]
  ),
  y_high = c(
    hab_ML_L_hab$RI_CI[2],
    sex_ML_L_sex$ci[2],
    egg_ML_L_boot$ci[2],
    larv_ML_L_boot$ci[2]
  ),
  notes = c(
    "female feeding-choice data; exact binomial CI",
    "approximate bootstrap from published mating success and total attempts",
    "approximate binomial bootstrap from published hatching rate and total number of eggs",
    "approximate binomial bootstrap from published adulthood survivorship and number of larvae"
  )
)

print(ri_ML_results)

plot_ML_df <- ri_ML_results %>%
  mutate(
    barrier = factor(
      barrier,
      levels = c("habitat", "sexual", "prehatching", "hybrid inviability", "total")
    ),
    strength.1 = strength_M,
    strength.2 = strength_L
  )

plot_ri_ML <- make_ri_strength_plot(
  plot_df = plot_ML_df,
  x_lab = expression(atop("\n", italic(RI) ~ "strength" ~ italic(M) ~ "-race")),
  y_lab = expression(atop("\n", italic(RI) ~ "strength" ~ italic(L) ~ "-race")),
  title_expr = expression(italic(F[ST]) == 0.367 ~ ";" ~ italic(D[XY]) == 0.0049),
  show_legend = FALSE
)

############################################################
## 8. REALISED FIELD-CAGE RI ESTIMATES
############################################################

############################################################
## 8A. H. diekei M x D: realised habitat isolation
## From field-cage host sightings.
## M-race: C = Mikania sightings, H = Dicliptera sightings.
## D-race: C = Dicliptera sightings, H = Mikania sightings.
## Net sightings are excluded from RI but retained as off-host use.
############################################################

fidelity_MD <- read.csv(path_fidelity)

realised_habitat_MD_point <- fidelity_MD %>%
  filter(arr %in% c("sep", "check")) %>%
  group_by(arr, race) %>%
  summarise(
    mik = sum(mik, na.rm = TRUE),
    dic = sum(dic, na.rm = TRUE),
    net = sum(net, na.rm = TRUE),
    total_sightings = mik + dic + net,
    .groups = "drop"
  ) %>%
  mutate(
    axis = case_when(
      race == "mrace" ~ "M-race",
      race == "drace" ~ "D-race",
      TRUE ~ as.character(race)
    ),
    C = case_when(
      race == "mrace" ~ mik,
      race == "drace" ~ dic,
      TRUE ~ NA_real_
    ),
    H = case_when(
      race == "mrace" ~ dic,
      race == "drace" ~ mik,
      TRUE ~ NA_real_
    ),
    net_prop = net / total_sightings,
    RI = ri_from_C_H(C, H)
  )

realised_habitat_MD_boot <- realised_habitat_MD_point %>%
  rowwise() %>%
  mutate(
    boot = list(bootstrap_realised_habitat_counts(C, H, nboot = nboot_main, seed = 1830)),
    RI_boot_mean = boot$mean,
    RI_boot_se = boot$se,
    RI_low = boot$ci[1],
    RI_high = boot$ci[2]
  ) %>%
  ungroup() %>%
  dplyr::select(-boot)

realised_habitat_MD_table <- realised_habitat_MD_boot %>%
  transmute(
    lineage = "H. diekei (Mikania x Dicliptera)",
    arrangement = arr,
    estimate = "realised habitat isolation",
    axis = axis,
    C,
    H,
    off_host_net = net,
    total_observed = total_sightings,
    off_host_net_prop = net_prop,
    RI,
    RI_boot_mean,
    RI_boot_se,
    RI_low,
    RI_high,
    method = "bootstrap from field-cage host-sighting counts",
    notes = "C and H are conspecific-host and heterospecific-host sightings; net sightings excluded from RI."
  )

############################################################
## 8B. H. diekei M x D: realised sexual isolation
## From field-cage mating observations.
## par = conspecific mating pair; hyb = heterospecific mating pair;
## null = no mating. Keep sep and check only.
############################################################

choice_mating_MD <- read.csv(path_choice_mating)

realised_sexual_MD_point <- choice_mating_MD %>%
  filter(arr %in% c("sep", "check")) %>%
  group_by(arr, male) %>%
  summarise(
    par = sum(par, na.rm = TRUE),
    hyb = sum(hyb, na.rm = TRUE),
    null = sum(null, na.rm = TRUE),
    mating_observed = par + hyb,
    total_observed = par + hyb + null,
    mating_rate = mating_observed / total_observed,
    .groups = "drop"
  ) %>%
  mutate(
    axis = case_when(
      male == "mrace" ~ "M-race male",
      male == "drace" ~ "D-race male",
      TRUE ~ as.character(male)
    ),
    RI = ri_from_C_H(par, hyb)
  )

realised_sexual_MD_boot <- realised_sexual_MD_point %>%
  rowwise() %>%
  mutate(
    boot = list(bootstrap_realised_sexual_counts(par, hyb, nboot = nboot_main, seed = 1830)),
    RI_boot_mean = boot$mean,
    RI_boot_se = boot$se,
    RI_low = boot$ci[1],
    RI_high = boot$ci[2]
  ) %>%
  ungroup() %>%
  dplyr::select(-boot)

realised_sexual_MD_table <- realised_sexual_MD_boot %>%
  transmute(
    lineage = "H. diekei (Mikania x Dicliptera)",
    arrangement = arr,
    estimate = "realised sexual isolation",
    axis = axis,
    C = par,
    H = hyb,
    null,
    mating_observed,
    total_observed,
    mating_rate,
    RI,
    RI_boot_mean,
    RI_boot_se,
    RI_low,
    RI_high,
    method = "bootstrap from field-cage mating-pair counts",
    notes = "C = par/conspecific mating; H = hyb/heterospecific mating; null excluded from RI but retained for mating rate."
  )

############################################################
## 8C. H. diekei M x L: realised sexual isolation
## From published field-cage Table 2.
## S-cage = separate arrangement; C-cage = checkerboard arrangement.
############################################################

realised_sexual_ML_counts <- tibble(
  arrangement = c("sep", "check"),
  cage = c("S-cage", "C-cage"),
  design = c("Separate", "Checkerboard"),
  LL = c(16, 15),
  MM = c(12, 16),
  LM = c(0, 2),  # L female x M male
  ML = c(0, 1)   # M female x L male
) %>%
  mutate(
    total_matings = LL + MM + LM + ML,
    intra_racial_matings = LL + MM,
    inter_racial_matings = LM + ML
  )

# Female-specific axes:
# L-race female: C = LL, H = LM
# M-race female: C = MM, H = ML
realised_sexual_ML_long <- bind_rows(
  realised_sexual_ML_counts %>%
    transmute(
      lineage = "H. diekei (Mikania x Leucas)",
      arrangement,
      estimate = "realised sexual isolation",
      axis = "L-race female",
      C = LL,
      H = LM,
      total_matings
    ),
  realised_sexual_ML_counts %>%
    transmute(
      lineage = "H. diekei (Mikania x Leucas)",
      arrangement,
      estimate = "realised sexual isolation",
      axis = "M-race female",
      C = MM,
      H = ML,
      total_matings
    )
)

realised_sexual_ML_table <- realised_sexual_ML_long %>%
  rowwise() %>%
  mutate(
    boot = list(bootstrap_realised_sexual_counts(C, H, nboot = nboot_main, seed = 1830)),
    RI = boot$point,
    RI_boot_mean = boot$mean,
    RI_boot_se = boot$se,
    RI_low = boot$ci[1],
    RI_high = boot$ci[2]
  ) %>%
  ungroup() %>%
  dplyr::select(-boot) %>%
  mutate(
    method = "bootstrap from published field-cage mating-pair counts",
    notes = "C = intra-racial mating for that female race; H = inter-racial mating for that female race."
  )

############################################################
## 8D. H. diekei M x L: realised habitat isolation
## From published field-cage Table 3.
## Uses original-host vs alternative-host sightings.
## Net sightings are excluded from RI but retained as off-host use.
##
## For realised RI comparable to the two-way field cages, the main values
## are S-cage and C-cage. L-cage and M-cage are retained in the output
## as single-host cage controls.
############################################################

realised_habitat_ML_counts <- tibble(
  cage = c(
    "L-cage", "L-cage", "L-cage",
    "M-cage", "M-cage", "M-cage",
    "S-cage", "S-cage", "S-cage", "S-cage",
    "C-cage", "C-cage", "C-cage", "C-cage"
  ),
  race = c(
    "Leucas", "Mikania", "Leucas",
    "Leucas", "Mikania", "Mikania",
    "Leucas", "Leucas", "Mikania", "Mikania",
    "Leucas", "Leucas", "Mikania", "Mikania"
  ),
  sex = c(
    "Female", "Female", "Male",
    "Female", "Female", "Male",
    "Female", "Male", "Female", "Male",
    "Female", "Male", "Female", "Male"
  ),
  n_collected_adults = c(13, 11, 13, 13, 13, 14, 11, 14, 8, 11, 10, 11, 15, 13),
  total_sightings = c(183, 132, 181, 200, 157, 155, 158, 167, 175, 196, 148, 157, 208, 195),
  original_host = c(173, 132, 173, 196, 154, 153, 156, 166, 167, 192, 148, 157, 206, 195),
  alternative_host = c(3, 0, 2, 0, 0, 2, 2, 1, 8, 5, 0, 0, 2, 0),
  net = c(7, 0, 6, 4, 3, 0, 1, 4, 1, 7, 3, 2, 0, 0),
  dispersal_average = c(1.12, 0.32, 1.37, 0.55, 0.91, 0.93, 0.83, 0.78, 0.63, 0.41, 1.27, 0.98, 0.77, 0.15),
  dispersal_se = c(0.11, 0.11, 0.23, 0.17, 0.20, 0.20, 0.18, 0.13, 0.16, 0.12, 0.20, 0.17, 0.17, 0.05),
  intra_host_migration = c(0.36, 0.12, 0.43, 0.20, 0.33, 0.39, 0.43, 0.43, 0.32, 0.24, 0.56, 0.47, 0.35, 0.30),
  inter_host_migration = c(0.03, 0, 0.03, 0, 0, 0.03, 0.03, 0.01, 0.06, 0.04, 0.00, 0.00, 0.01, 0.01)
) %>%
  mutate(
    arrangement = case_when(
      cage == "S-cage" ~ "sep",
      cage == "C-cage" ~ "check",
      cage == "L-cage" ~ "L-cage control",
      cage == "M-cage" ~ "M-cage control"
    ),
    axis = paste0(race, " ", sex),
    C = original_host,
    H = alternative_host,
    RI = ri_from_C_H(C, H),
    net_prop = net / total_sightings
  )

realised_habitat_ML_table <- realised_habitat_ML_counts %>%
  rowwise() %>%
  mutate(
    boot = list(bootstrap_realised_habitat_counts(C, H, nboot = nboot_main, seed = 1830)),
    RI_boot_mean = boot$mean,
    RI_boot_se = boot$se,
    RI_low = boot$ci[1],
    RI_high = boot$ci[2]
  ) %>%
  ungroup() %>%
  dplyr::select(-boot) %>%
  transmute(
    lineage = "H. diekei (Mikania x Leucas)",
    arrangement,
    cage,
    estimate = "realised habitat isolation",
    axis,
    C,
    H,
    off_host_net = net,
    total_observed = total_sightings,
    off_host_net_prop = net_prop,
    RI,
    RI_boot_mean,
    RI_boot_se,
    RI_low,
    RI_high,
    method = "bootstrap from published field-cage host-sighting counts",
    notes = "C = original-host sightings; H = alternative-host sightings; net excluded from RI. S-cage and C-cage are the two-way field-cage treatments."
  )

############################################################
## 8E. Save realised field-cage RI outputs
############################################################

realised_field_cage_RI_summary <- bind_rows(
  realised_habitat_MD_table,
  realised_sexual_MD_table,
  realised_habitat_ML_table,
  realised_sexual_ML_table
)

realised_field_cage_RI_summary_rounded <- realised_field_cage_RI_summary %>%
  mutate(
    across(
      any_of(c(
        "RI", "RI_boot_mean", "RI_boot_se", "RI_low", "RI_high",
        "off_host_net_prop", "mating_rate"
      )),
      ~ round(.x, 4)
    )
  )

print(realised_field_cage_RI_summary)

write.csv(
  realised_field_cage_RI_summary,
  file = "data/07c_2026_EE_output_ri_estimates_realised_summary.csv",
  row.names = FALSE
)


############################################################
## 9. EXPORT RI ESTIMATE TABLES
############################################################

## H. diekei M x D

ri_MD_results <- make_ri_estimate_table(
  lineage = "H. diekei (Mikania x Dicliptera)",
  axis_1 = "M-race",
  axis_2 = "D-race",
  barriers = c("habitat", "sexual", "prehatching_weighted", "prehatching_unweighted", "hybrid inviability"),
  strength_1 = c(
    hab_M_res$RI_hat,
    sex_M_res$point,
    egg_M_point_weighted$RI,
    egg_M_point_unweighted$RI,
    larv_M_point_weighted$RI
  ),
  strength_2 = c(
    hab_D_res$RI_hat,
    sex_D_res$point,
    egg_D_point_weighted$RI,
    egg_D_point_unweighted$RI,
    larv_D_point_weighted$RI
  ),
  low_1 = c(
    hab_M_res$RI_CI[1],
    sex_M_res$ci[1],
    egg_M_boot_weighted$ci[1],
    egg_M_boot_unweighted$ci[1],
    larv_M_boot_weighted$ci[1]
  ),
  high_1 = c(
    hab_M_res$RI_CI[2],
    sex_M_res$ci[2],
    egg_M_boot_weighted$ci[2],
    egg_M_boot_unweighted$ci[2],
    larv_M_boot_weighted$ci[2]
  ),
  low_2 = c(
    hab_D_res$RI_CI[1],
    sex_D_res$ci[1],
    egg_D_boot_weighted$ci[1],
    egg_D_boot_unweighted$ci[1],
    larv_D_boot_weighted$ci[1]
  ),
  high_2 = c(
    hab_D_res$RI_CI[2],
    sex_D_res$ci[2],
    egg_D_boot_weighted$ci[2],
    egg_D_boot_unweighted$ci[2],
    larv_D_boot_weighted$ci[2]
  ),
  method = c(
    "exact binomial CI",
    "bootstrap from binary mating data",
    "bootstrap, weighted hatchability",
    "bootstrap, unweighted hatchability",
    "bootstrap, weighted larval survival"
  ),
  notes = c(
    "Habitat isolation",
    "Sexual isolation",
    "Postmating prehatching isolation; weighted estimate used in plot",
    "Postmating prehatching isolation; sensitivity estimate",
    "Hybrid inviability / larval performance"
  )
)

## H. diekei M x L

ri_ML_results_export <- make_ri_estimate_table(
  lineage = "H. diekei (Mikania x Leucas)",
  axis_1 = "M-race",
  axis_2 = "L-race",
  barriers = ri_ML_results$barrier,
  strength_1 = ri_ML_results$strength_M,
  strength_2 = ri_ML_results$strength_L,
  low_1 = ri_ML_results$x_low,
  high_1 = ri_ML_results$x_high,
  low_2 = ri_ML_results$y_low,
  high_2 = ri_ML_results$y_high,
  method = c(
    "exact binomial CI",
    "approximate bootstrap from published mating success",
    "approximate binomial bootstrap from published hatching rate and total number of eggs",
    "approximate binomial bootstrap from published adulthood survivorship and number of larvae"
  ),
  notes = ri_ML_results$notes
)


## Combined table and output

ri_estimates_table <- bind_rows(
  ri_MD_results,
  ri_ML_results_export
)

print(ri_estimates_table)

write.csv(
  ri_estimates_table,
  file = "data/07c_2026_EE_output_ri_estimates_summary2.csv",
  row.names = FALSE
)


############################################################
## 10. COMBINE AND SAVE PLOTS
############################################################

plot_ri_ML2 <- plot_ri_ML +
  theme(axis.title.y = element_text(size = 12))

plot_ri22 <- plot_ri2 +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

plot_model42 <- plot.ri.model4 +
  theme(axis.title.y = element_text(size = 12))

plot_model32 <- plot.ri.model3 +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

plot.ri.comb.4 <- ggarrange(
  plot_ri2,
  plot.ri.model3,
  labels = c("a", "b"),
  heights = c(3, 3),
  widths = c(3, 3),
  ncol = 1,
  nrow = 2,
  align = "v"
)

final_plot <- (plot_ri_ML2 | plot_ri22) /
  (plot_model42 | plot_model32) +
  plot_annotation(tag_levels = list(c("a", "", "b", "")))

final_plot

path_plot_base     <- "figures/diekei_dic_ri"

ggsave(filename = paste0(path_plot_base, ".png"), width = 6.5, height = 6, device = "png", dpi = 1200)
ggsave(filename = paste0(path_plot_base, ".pdf"), width = 6.5, height = 6, device = "pdf")


############################################################
## 11. PLOT REALISED RI + LAB EXPERIMENT RI
############################################################

############################################################
## Helper: plot HI and SI only
############################################################

barrier_cols_HI_SI <- c(
  "habitat" = "#8dd3c7",
  "sexual"  = "#fccde5"
)

context_shapes <- c(
  "without host plant" = 21,
  "with host plant parapatry" = 22,
  "with host plant sympatry" = 24
)

make_HI_SI_context_plot <- function(plot_df, x_lab, y_lab, title_expr,
                                    show_legend = TRUE) {
  
  legend_theme <- if (show_legend) {
    theme(
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.25, "cm"),
      legend.box.spacing = unit(0.5, "cm")
    )
  } else {
    theme(legend.position = "none")
  }
  
  ggplot(plot_df, aes(x = strength.1, y = strength.2)) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "black", linewidth = 0.5
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed", color = "black", linewidth = 0.5
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed", color = "black", linewidth = 0.5
    ) +
    
    geom_segment(
      aes(
        x = strength.1, xend = strength.1,
        y = y_low, yend = y_high,
        colour = barrier
      ),
      linewidth = 1
    ) +
    
    geom_segment(
      aes(
        y = strength.2, yend = strength.2,
        x = x_low, xend = x_high,
        colour = barrier
      ),
      linewidth = 1
    ) +
    
    geom_point(
      aes(fill = barrier, shape = context),
      size = 4,
      stroke = 0.7,
      colour = "black"
    ) +
    
    scale_colour_manual(
      values = barrier_cols_HI_SI,
      limits = c("habitat", "sexual"),
      labels = c("habitat isolation", "sexual isolation"),
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = barrier_cols_HI_SI,
      limits = c("habitat", "sexual"),
      labels = c("habitat isolation", "sexual isolation"),
      drop = FALSE
    ) +
    
    scale_shape_manual(
      values = context_shapes,
      limits = c(
        "without host plant",
        "with host plant parapatry",
        "with host plant sympatry"
      ),
      drop = FALSE
    ) +
    
    coord_cartesian(xlim = c(-0.1, 1), ylim = c(-0.1, 1)) +
    
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +
    
    xlab(x_lab) +
    ylab(y_lab) +
    
    theme_bw() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 16),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    legend_theme +
    
    # Important:
    # - hide colour legend so CI lines do not create a second black-looking legend
    # - use fill legend for coloured barrier points
    guides(
      colour = "none",
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          size = 4,
          colour = "black"
        ),
        order = 1
      ),
      shape = guide_legend(
        override.aes = list(
          fill = "grey80",
          colour = "black",
          size = 4
        ),
        order = 2
      )
    ) +
    
    ggtitle(title_expr)
}

############################################################
## H. diekei M x D: build HI + SI plot data
############################################################

# Lab experiment / no host plant
lab_MD_plot_df <- tibble(
  context = "without host plant",
  barrier = c("habitat", "sexual"),
  strength.1 = c(
    hab_M_res$RI_hat,
    sex_M_res$point
  ),
  strength.2 = c(
    hab_D_res$RI_hat,
    sex_D_res$point
  ),
  x_low = c(
    hab_M_res$RI_CI[1],
    sex_M_res$ci[1]
  ),
  x_high = c(
    hab_M_res$RI_CI[2],
    sex_M_res$ci[2]
  ),
  y_low = c(
    hab_D_res$RI_CI[1],
    sex_D_res$ci[1]
  ),
  y_high = c(
    hab_D_res$RI_CI[2],
    sex_D_res$ci[2]
  )
)

# Field cage / with host plants: habitat
realised_habitat_MD_plot <- realised_habitat_MD_table %>%
  filter(arrangement %in% c("sep", "check")) %>%
  mutate(
    context = case_when(
      arrangement == "sep" ~ "with host plant parapatry",
      arrangement == "check" ~ "with host plant sympatry"
    ),
    barrier = "habitat",
    axis_plot = case_when(
      axis == "M-race" ~ "M",
      axis == "D-race" ~ "D",
      TRUE ~ axis
    )
  ) %>%
  dplyr::select(context, barrier, axis_plot, RI, RI_low, RI_high) %>%
  pivot_wider(
    names_from = axis_plot,
    values_from = c(RI, RI_low, RI_high)
  ) %>%
  transmute(
    context,
    barrier,
    strength.1 = RI_M,
    strength.2 = RI_D,
    x_low = RI_low_M,
    x_high = RI_high_M,
    y_low = RI_low_D,
    y_high = RI_high_D
  )

# Field cage / with host plants: sexual
realised_sexual_MD_plot <- realised_sexual_MD_table %>%
  filter(arrangement %in% c("sep", "check")) %>%
  mutate(
    context = case_when(
      arrangement == "sep" ~ "with host plant parapatry",
      arrangement == "check" ~ "with host plant sympatry"
    ),
    barrier = "sexual",
    axis_plot = case_when(
      axis == "M-race male" ~ "M",
      axis == "D-race male" ~ "D",
      TRUE ~ axis
    )
  ) %>%
  dplyr::select(context, barrier, axis_plot, RI, RI_low, RI_high) %>%
  pivot_wider(
    names_from = axis_plot,
    values_from = c(RI, RI_low, RI_high)
  ) %>%
  transmute(
    context,
    barrier,
    strength.1 = RI_M,
    strength.2 = RI_D,
    x_low = RI_low_M,
    x_high = RI_high_M,
    y_low = RI_low_D,
    y_high = RI_high_D
  )

HI_SI_MD_plot_df <- bind_rows(
  lab_MD_plot_df,
  realised_habitat_MD_plot,
  realised_sexual_MD_plot
) %>%
  mutate(
    context = factor(
      context,
      levels = c(
        "without host plant",
        "with host plant parapatry",
        "with host plant sympatry"
      )
    ),
    barrier = factor(barrier, levels = c("habitat", "sexual"))
  )

plot_HI_SI_MD <- make_HI_SI_context_plot(
  plot_df = HI_SI_MD_plot_df,
  x_lab = expression(atop("\n", italic(RI) ~ "strength" ~ italic(M) ~ "-race")),
  y_lab = expression(atop("\n", italic(RI) ~ "strength" ~ italic(D) ~ "-race")),
  title_expr = expression(italic(F[ST]) == 1.0000 ~ ";" ~ italic(D[XY]) == 0.0079),
  show_legend = TRUE
)

############################################################
## H. diekei M x L: build HI + SI plot data
############################################################

# Lab experiment / no host plant
lab_ML_plot_df <- tibble(
  context = "without host plant",
  barrier = c("habitat", "sexual"),
  strength.1 = c(
    hab_ML_M_hab$RI_hat,
    sex_ML_M_sex$point
  ),
  strength.2 = c(
    hab_ML_L_hab$RI_hat,
    sex_ML_L_sex$point
  ),
  x_low = c(
    hab_ML_M_hab$RI_CI[1],
    sex_ML_M_sex$ci[1]
  ),
  x_high = c(
    hab_ML_M_hab$RI_CI[2],
    sex_ML_M_sex$ci[2]
  ),
  y_low = c(
    hab_ML_L_hab$RI_CI[1],
    sex_ML_L_sex$ci[1]
  ),
  y_high = c(
    hab_ML_L_hab$RI_CI[2],
    sex_ML_L_sex$ci[2]
  )
)

# Field cage / with host plants: habitat
# Pool female and male sightings within race and arrangement.
realised_habitat_ML_pooled <- realised_habitat_ML_counts %>%
  filter(arrangement %in% c("sep", "check")) %>%
  group_by(arrangement, race) %>%
  summarise(
    C = sum(original_host, na.rm = TRUE),
    H = sum(alternative_host, na.rm = TRUE),
    off_host_net = sum(net, na.rm = TRUE),
    total_observed = sum(total_sightings, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    boot = list(bootstrap_realised_habitat_counts(C, H, nboot = nboot_main, seed = 1830)),
    RI = boot$point,
    RI_low = boot$ci[1],
    RI_high = boot$ci[2],
    axis_plot = case_when(
      race == "Mikania" ~ "M",
      race == "Leucas" ~ "L",
      TRUE ~ race
    )
  ) %>%
  ungroup() %>%
  dplyr::select(-boot)

realised_habitat_ML_plot <- realised_habitat_ML_pooled %>%
  mutate(
    context = case_when(
      arrangement == "sep" ~ "with host plant parapatry",
      arrangement == "check" ~ "with host plant sympatry"
    ),
    barrier = "habitat"
  ) %>%
  dplyr::select(context, barrier, axis_plot, RI, RI_low, RI_high) %>%
  pivot_wider(
    names_from = axis_plot,
    values_from = c(RI, RI_low, RI_high)
  ) %>%
  transmute(
    context,
    barrier,
    strength.1 = RI_M,
    strength.2 = RI_L,
    x_low = RI_low_M,
    x_high = RI_high_M,
    y_low = RI_low_L,
    y_high = RI_high_L
  )

# Field cage / with host plants: sexual
realised_sexual_ML_plot <- realised_sexual_ML_table %>%
  filter(arrangement %in% c("sep", "check")) %>%
  mutate(
    context = case_when(
      arrangement == "sep" ~ "with host plant parapatry",
      arrangement == "check" ~ "with host plant sympatry"
    ),
    barrier = "sexual",
    axis_plot = case_when(
      axis == "M-race female" ~ "M",
      axis == "L-race female" ~ "L",
      TRUE ~ axis
    )
  ) %>%
  dplyr::select(context, barrier, axis_plot, RI, RI_low, RI_high) %>%
  pivot_wider(
    names_from = axis_plot,
    values_from = c(RI, RI_low, RI_high)
  ) %>%
  transmute(
    context,
    barrier,
    strength.1 = RI_M,
    strength.2 = RI_L,
    x_low = RI_low_M,
    x_high = RI_high_M,
    y_low = RI_low_L,
    y_high = RI_high_L
  )

HI_SI_ML_plot_df <- bind_rows(
  lab_ML_plot_df,
  realised_habitat_ML_plot,
  realised_sexual_ML_plot
) %>%
  mutate(
    context = factor(
      context,
      levels = c(
        "without host plant",
        "with host plant parapatry",
        "with host plant sympatry"
      )
    ),
    barrier = factor(barrier, levels = c("habitat", "sexual"))
  )

plot_HI_SI_ML <- make_HI_SI_context_plot(
  plot_df = HI_SI_ML_plot_df,
  x_lab = expression(atop("\n", italic(RI) ~ "strength" ~ italic(M) ~ "-race")),
  y_lab = expression(atop("\n", italic(RI) ~ "strength" ~ italic(L) ~ "-race")),
  title_expr = expression(italic(F[ST]) == 0.367 ~ ";" ~ italic(D[XY]) == 0.0049),
  show_legend = FALSE
)

############################################################
## Combined HI + SI figure
############################################################

final_HI_SI_context_plot <- plot_HI_SI_ML | plot_HI_SI_MD

final_HI_SI_context_plot

ggsave(
  filename = "figures/diekei_dic_ri_context.png",
  plot = final_HI_SI_context_plot,
  width = 7,
  height = 3,
  device = "png",
  dpi = 1200
)

ggsave(
  filename = "figures/diekei_dic_ri_context.pdf",
  plot = final_HI_SI_context_plot,
  width = 7,
  height = 3,
  device = "pdf"
)


############################################################
## Save plot data
############################################################

write.csv(
  HI_SI_ML_plot_df,
  file = "data/07c_2026_EE_output_ri_estimates_realised_summary_ml.csv",
  row.names = FALSE
)

write.csv(
  HI_SI_MD_plot_df,
  file = "data/07c_2026_EE_output_ri_estimates_realised_summary_md.csv",
  row.names = FALSE
)
