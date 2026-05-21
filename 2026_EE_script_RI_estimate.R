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
## 4. H. diekei M x D: UPDATED RI ESTIMATES
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

# Egg hatchability: point estimate only
egg_ML_M_egg_RI <- ri_from_C_H(C = 0.723, H = 0.702)
egg_ML_L_egg_RI <- ri_from_C_H(C = 0.643, H = 0.678)

# Hybrid inviability / larval performance: point estimate only
larv_ML_M_larv_RI <- ri_from_C_H(C = 0.64, H = 0.65)
larv_ML_L_larv_RI <- ri_from_C_H(C = 0.81, H = 0.79)

cat("\n===== H. diekei M x L: HABITAT ISOLATION =====\n")
print_ri_line("Mikania race", hab_ML_M_hab$RI_hat, hab_ML_M_hab$RI_CI[1], hab_ML_M_hab$RI_CI[2])
print_ri_line("Leucas race", hab_ML_L_hab$RI_hat, hab_ML_L_hab$RI_CI[1], hab_ML_L_hab$RI_CI[2])

cat("\n===== H. diekei M x L: SEXUAL ISOLATION =====\n")
print_ri_line("Mikania race", sex_ML_M_sex$point, sex_ML_M_sex$ci[1], sex_ML_M_sex$ci[2])
print_ri_line("Leucas race", sex_ML_L_sex$point, sex_ML_L_sex$ci[1], sex_ML_L_sex$ci[2])

cat("\n===== H. diekei M x L: EGG HATCHABILITY =====\n")
print_ri_line("Mikania race", egg_ML_M_egg_RI)
print_ri_line("Leucas race", egg_ML_L_egg_RI)

cat("\n===== H. diekei M x L: HYBRID INVIABILITY / LARVAL PERFORMANCE =====\n")
print_ri_line("Mikania race", larv_ML_M_larv_RI)
print_ri_line("Leucas race", larv_ML_L_larv_RI)

ri_ML_results <- tibble(
  lineage = "H. diekei (Mikania x Leucas)",
  barrier = c("habitat", "sexual", "prehatching", "hybrid inviability"),
  strength_M = c(
    hab_ML_M_hab$RI_hat,
    sex_ML_M_sex$point,
    egg_ML_M_egg_RI,
    larv_ML_M_larv_RI
  ),
  strength_L = c(
    hab_ML_L_hab$RI_hat,
    sex_ML_L_sex$point,
    egg_ML_L_egg_RI,
    larv_ML_L_larv_RI
  ),
  x_low = c(
    hab_ML_M_hab$RI_CI[1],
    sex_ML_M_sex$ci[1],
    egg_ML_M_egg_RI,
    larv_ML_M_larv_RI
  ),
  x_high = c(
    hab_ML_M_hab$RI_CI[2],
    sex_ML_M_sex$ci[2],
    egg_ML_M_egg_RI,
    larv_ML_M_larv_RI
  ),
  y_low = c(
    hab_ML_L_hab$RI_CI[1],
    sex_ML_L_sex$ci[1],
    egg_ML_L_egg_RI,
    larv_ML_L_larv_RI
  ),
  y_high = c(
    hab_ML_L_hab$RI_CI[2],
    sex_ML_L_sex$ci[2],
    egg_ML_L_egg_RI,
    larv_ML_L_larv_RI
  ),
  notes = c(
    "female feeding-choice data; exact binomial CI",
    "approximate bootstrap from published mating success and total attempts",
    "point estimate only from published hatchability",
    "point estimate only from published adulthood survivorship on maternal host"
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
## 8. COMBINE AND SAVE PLOTS
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
