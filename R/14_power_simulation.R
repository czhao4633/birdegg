# 14_power_simulation.R
# How many independent terrestrial flight losses would be needed to test the
# HWI-sized reversion?
#
# Design (Montgomerie dataset, MCC tree, model M1):
#   * Residuals are simulated from the fitted Pagel's-lambda model of each
#     trait: e ~ N(0, sigma2 * [lambda * C + (1 - lambda) * diag(C)]), with
#     sigma2, lambda and the mass / egg-size slopes taken from the real fit.
#   * The 5 real terrestrial origins (12 species) are kept; k - 5 extra
#     origins are created by relabelling randomly chosen volant species as
#     terrestrial. Extra origins are single species, like most extant
#     flightless birds without egg data (mostly island rails).
#   * Two scenarios share the same residuals and the same extra origins:
#       H0-like : true terrestrial effect = 0
#       H1-like : true effect = -benchmark (HWI-predicted reversion, from
#                 output/tables/equivalence_tests.csv)
#   * Each simulated data set is refitted with M1. Recorded:
#       power   = P(one-sided p < 0.05 for a negative effect | reversion)
#       exclude = P(delta_min < benchmark | no effect), i.e. probability that
#                 a one-sided equivalence test rules out the benchmark.
#
# Caveats: MCC tree only; extra origins placed at random volant tips (real
# candidates are concentrated in rails and ducks); real flightless species
# would also bring new covariate values. Treat results as order-of-magnitude.
#
# Each replicate has its own seed (SEED + offset), so results do not depend
# on the number of cores or the scheduling of parallel workers.
# Outputs: output/tables/power_simulation.csv,
#          data/processed/power_simulation_reps.rds (all replicates),
#          output/figures/figS1_power_simulation.{pdf,png}

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(parallel)
})

N_REP   <- 1000L
K_GRID  <- c(5L, 8L, 10L, 15L, 20L, 30L, 45L, 60L)  # terrestrial origins
TRAITS  <- c("elongation", "pointedness", "polar_asym")
N_CORES <- max(1L, min(8L, parallel::detectCores() - 1L))

mcc <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
d <- as.data.frame(readRDS(file.path(PATHS$processed, "analysis_montgomerie.rds")))
rownames(d) <- d$birdtree_name
tr <- keep.tip(mcc, d$birdtree_name)
d <- d[tr$tip.label, ]
stopifnot(is.ultrametric(tr))
n <- nrow(d)
volant_tips <- which(d$flight_status == "volant")
C <- vcv(tr)

bench <- read.csv(file.path(PATHS$tables, "equivalence_tests.csv")) |>
  filter(model == "M1_status") |>
  select(trait, benchmark_sd)

# ---- Fitted generating model per trait ----------------------------------------
gen <- lapply(setNames(TRAITS, TRAITS), function(trait) {
  f <- phylolm(as.formula(paste(trait, "~ flight_status + log_mass + log_egg_size")),
               data = d, phy = tr, model = "lambda")
  lam <- f$optpar
  V <- f$sigma2 * (lam * C + (1 - lam) * diag(diag(C)))
  list(b0 = coef(f)[["(Intercept)"]], b_mass = coef(f)[["log_mass"]],
       b_egg = coef(f)[["log_egg_size"]], b_div = coef(f)[["flight_statusdiving"]],
       L = t(chol(V)),
       sd_vol = sd(d[[trait]][d$flight_status == "volant"]),
       bench = bench$benchmark_sd[bench$trait == trait])
})

one_rep <- function(ti, ki, r) {
  trait <- TRAITS[ti]; k <- K_GRID[ki]; g <- gen[[trait]]
  set.seed(SEED + r + 10000L * ki + 1000000L * ti)
  status <- as.character(d$flight_status)
  if (k > 5) status[sample(volant_tips, k - 5)] <- "terrestrial"
  status <- factor(status, levels = c("volant", "terrestrial", "diving"))
  e <- as.vector(g$L %*% rnorm(n))
  mu <- g$b0 + g$b_mass * d$log_mass + g$b_egg * d$log_egg_size +
    g$b_div * (status == "diving")
  dd <- data.frame(status = status, log_mass = d$log_mass,
                   log_egg_size = d$log_egg_size, row.names = rownames(d))
  out <- list()
  for (sc in c("no_effect", "reversion")) {
    eff <- if (sc == "reversion") -g$bench * g$sd_vol else 0
    dd$y <- mu + eff * (status == "terrestrial") + e
    s <- summary(phylolm(y ~ status + log_mass + log_egg_size, data = dd, phy = tr,
                         model = "lambda"))$coefficients
    est <- s["statusterrestrial", "Estimate"]; se <- s["statusterrestrial", "StdErr"]
    dfree <- n - nrow(s)
    out[[sc]] <- data.frame(
      trait = trait, n_origins = k, rep = r, scenario = sc,
      std_est = est / g$sd_vol,
      p_one_sided_neg = pt(est / se, dfree),
      delta_min = (qt(0.95, dfree) * se - est) / g$sd_vol)
  }
  do.call(rbind, out)
}

grid <- expand.grid(r = seq_len(N_REP), ki = seq_along(K_GRID), ti = seq_along(TRAITS))
# The ~10-min simulation is cached; the cache is reused only if it was made
# with the same settings (seed, replicates, grid, traits, benchmarks).
cache_file <- file.path(PATHS$processed, "power_simulation_reps.rds")
settings <- list(SEED = SEED, N_REP = N_REP, K_GRID = K_GRID, TRAITS = TRAITS,
                 bench = bench$benchmark_sd)
sims <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
if (is.null(sims) || !identical(attr(sims, "settings"), settings)) {
  log_msg("Simulating ", nrow(grid), " replicates x 2 scenarios on ", N_CORES, " cores ...")
  t0 <- Sys.time()
  sims <- bind_rows(mclapply(seq_len(nrow(grid)), function(j)
    one_rep(grid$ti[j], grid$ki[j], grid$r[j]), mc.cores = N_CORES))
  log_msg("Done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  attr(sims, "settings") <- settings
  saveRDS(sims, cache_file)
} else {
  log_msg("Using cached simulation: ", cache_file)
}
stopifnot(nrow(sims) == 2 * nrow(grid))

summ <- sims |>
  left_join(bench, by = "trait") |>
  group_by(trait, n_origins) |>
  summarise(power = mean(p_one_sided_neg[scenario == "reversion"] < 0.05),
            exclude = mean(delta_min[scenario == "no_effect"] < first(benchmark_sd)),
            median_delta_min = median(delta_min[scenario == "no_effect"]),
            mean_est_reversion = mean(std_est[scenario == "reversion"]),
            benchmark_sd = first(benchmark_sd),
            .groups = "drop") |>
  relocate(benchmark_sd, .after = n_origins) |>
  mutate(mc_se_power = sqrt(power * (1 - power) / N_REP),
         mc_se_exclude = sqrt(exclude * (1 - exclude) / N_REP),
         n_rep = N_REP) |>
  mutate(across(where(is.double), ~ signif(.x, 3)))
write_table(summ, "power_simulation")
print(summ, n = Inf, width = Inf)

# ---- Figure S1 -----------------------------------------------------------------
TRAIT_LAB <- c(elongation = "Elongation", pointedness = "Pointedness",
               polar_asym = "Polar asymmetry")
TRAIT_COL <- c(elongation = "#E69F00", pointedness = "#009E73", polar_asym = "#CC79A7")
theme_pub <- function(base = 8) {
  theme_classic(base_size = base) +
    theme(axis.line = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
          axis.text = element_text(colour = "grey20"),
          legend.key.size = unit(3, "mm"),
          plot.tag = element_text(face = "bold"))
}
# Ceiling from R/13: current origins + extant terrestrial origins without egg
# data (certain), and + those whose locomotion still needs review (maximum)
inv <- read.csv(file.path(PATHS$tables, "origin_inventory.csv"))
k_now <- sum(inv$has_egg_data & inv$locomotion == "terrestrial")
k_sure <- k_now + sum(!inv$has_egg_data & inv$locomotion == "terrestrial")
k_max <- k_sure + sum(!inv$has_egg_data & inv$locomotion == "review")
log_msg("Terrestrial origins: now ", k_now, ", extant certain ", k_sure, ", extant max ", k_max)
marks <- data.frame(x = c(k_now, k_sure, k_max),
                    lab = c(paste0("now (", k_now, ")"),
                            paste0("all extant\nterrestrial (", k_sure, ")"),
                            paste0("+ unclassified\n(", k_max, ")")))
panel <- function(col, ylab) {
  ggplot(summ, aes(n_origins, .data[[col]], colour = trait)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    geom_vline(data = marks, aes(xintercept = x), linewidth = 0.3, colour = "grey75") +
    geom_text(data = marks, aes(x = x + 0.8, y = 0.02, label = lab), inherit.aes = FALSE,
              hjust = 0, vjust = 0, size = 2.1, colour = "grey40", lineheight = 0.9) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.4) +
    scale_colour_manual(values = TRAIT_COL, labels = TRAIT_LAB, name = NULL) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(breaks = K_GRID) +
    labs(x = "Independent terrestrial flight losses", y = ylab) +
    theme_pub()
}
figS1 <- panel("power", "Power to detect HWI-sized reversion") +
  panel("exclude", "Probability of excluding it if absent") +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
for (ext in c("pdf", "png")) {
  ggsave(file.path(PATHS$figures, paste0("figS1_power_simulation.", ext)), figS1,
         width = 180, height = 80, units = "mm", dpi = 300,
         device = if (ext == "pdf") cairo_pdf else "png", bg = "white")
}
log_msg("Wrote figS1_power_simulation (.pdf, .png)")
