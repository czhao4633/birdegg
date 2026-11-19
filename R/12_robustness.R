# 12_robustness.R
# Three robustness checks:
#
#   A. Covariate dependence. M1 includes log egg size, which may lie on the
#      causal path from flight loss to egg shape (if flight loss changes egg
#      size relative to the female, adjusting for it removes part of the
#      effect). Refit the flight-status model with
#        M0_status_only : trait ~ flight_status
#        M6_status_mass : trait ~ flight_status + log_mass
#      in both datasets, on the MCC tree and the 100 posterior trees.
#   B. Equivalence (one-sided, directional). Reversion predicts a negative
#      terrestrial effect. For each model we report the smallest reversion
#      (in volant-SD units) that can be rejected at alpha = 0.05,
#        delta_min = (t_0.95,df * SE - estimate) / SD_volant,
#      and test against an a-priori benchmark derived from Stoddard et al.'s
#      own mechanism: the HWI slope (fitted on volant species only) times the
#      observed drop in HWI from each terrestrial species' volant sister group
#      to the species itself, averaged over terrestrial species.
#   C. Extrapolation diagnostics for the Montgomerie dataset (the same checks
#      that R/11 runs for Stoddard): ratites removed, mass ranges, and raw
#      ratite vs tinamou means.
#
# Outputs (output/tables/):
#   robustness_covariates.csv, equivalence_tests.csv,
#   montgomerie_diagnostics.csv
# No random numbers are drawn; results are deterministic.

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  library(dplyr)
  library(purrr)
  library(parallel)
})

mcc   <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
post  <- readRDS(file.path(PATHS$processed, "trees_posterior_100.rds"))
pairs <- readRDS(file.path(PATHS$processed, "sister_pairs.rds"))

DATA <- list(
  montgomerie = list(d = readRDS(file.path(PATHS$processed, "analysis_montgomerie.rds")),
                     traits = c("elongation", "pointedness", "polar_asym")),
  stoddard    = list(d = readRDS(file.path(PATHS$processed, "analysis_stoddard.rds")),
                     traits = c("ellipticity", "asymmetry"))
)

FORMULAS <- list(
  M1_status      = "~ flight_status + log_mass + log_egg_size",
  M0_status_only = "~ flight_status",
  M6_status_mass = "~ flight_status + log_mass"
)
RATITE_RE <- "^(Struthio|Rhea|Dromaius|Casuarius|Apteryx)_"
TERMS <- c("flight_statusterrestrial", "flight_statusdiving")
N_CORES <- max(1L, min(8L, parallel::detectCores() - 1L))

#' Fit one PGLS (Pagel's lambda); return the flight-status (or hwi) rows with
#' SE, residual df and the trait SD among volant species.
fit_terms <- function(d, tree, trait, rhs, terms = TERMS) {
  d <- as.data.frame(d)
  rownames(d) <- d$birdtree_name
  vars <- all.vars(as.formula(paste(trait, rhs)))
  d <- d[complete.cases(d[, vars]), , drop = FALSE]
  tr <- keep.tip(tree, d$birdtree_name)
  d <- d[tr$tip.label, , drop = FALSE]
  fit <- phylolm(as.formula(paste(trait, rhs)), data = d, phy = tr, model = "lambda")
  s <- summary(fit)$coefficients
  terms <- intersect(terms, rownames(s))
  sd_vol <- sd(d[[trait]][d$flight_status == "volant"])
  tibble(term = terms, estimate = s[terms, "Estimate"], se = s[terms, "StdErr"],
         p = s[terms, "p.value"], df = nrow(d) - nrow(s), n = nrow(d),
         lambda = fit$optpar, trait_sd_volant = sd_vol,
         std_estimate = s[terms, "Estimate"] / sd_vol)
}

#' Smallest reversion (SD units, predicted negative direction) rejected at
#' one-sided alpha = 0.05, and the one-sided p-value for H0: beta <= -delta.
delta_min <- function(est, se, df, sdv) (qt(0.95, df) * se - est) / sdv
p_noninf  <- function(est, se, df, sdv, delta) pt((est + delta * sdv) / se, df, lower.tail = FALSE)

# ---- A. Covariate dependence: MCC tree ----------------------------------------
mcc_fits <- list()
for (dn in names(DATA)) for (tr in DATA[[dn]]$traits) for (m in names(FORMULAS)) {
  mcc_fits[[length(mcc_fits) + 1]] <- fit_terms(DATA[[dn]]$d, mcc, tr, FORMULAS[[m]]) |>
    mutate(dataset = dn, trait = tr, model = m)
}
# Montgomerie without ratites (C; terrestrial = kagu, kakapo, nativehen, rail)
for (tr in DATA$montgomerie$traits) for (m in c("M1_status", "M0_status_only")) {
  d_nr <- DATA$montgomerie$d |> filter(!grepl(RATITE_RE, birdtree_name))
  mcc_fits[[length(mcc_fits) + 1]] <- fit_terms(d_nr, mcc, tr, FORMULAS[[m]]) |>
    mutate(dataset = "montgomerie", trait = tr, model = paste0(m, "_no_ratites"))
}
mcc_fits <- bind_rows(mcc_fits)

# ---- A. Covariate dependence: posterior trees ------------------------------------
# M1 posterior fits already exist (R/07); refit only the two new models.
log_msg("Refitting M0/M6 on ", length(post), " posterior trees (", N_CORES, " cores) ...")
jobs <- expand.grid(dataset = names(DATA), model = c("M0_status_only", "M6_status_mass"),
                    stringsAsFactors = FALSE)
post_new <- bind_rows(lapply(seq_len(nrow(jobs)), function(j) {
  dn <- jobs$dataset[j]; m <- jobs$model[j]
  bind_rows(mclapply(seq_along(post), function(i) {
    bind_rows(lapply(DATA[[dn]]$traits, function(tr)
      fit_terms(DATA[[dn]]$d, post[[i]], tr, FORMULAS[[m]]) |>
        mutate(trait = tr, tree = i)))
  }, mc.cores = N_CORES)) |> mutate(dataset = dn, model = m)
}))
post_m1 <- readRDS(file.path(PATHS$processed, "pgls_posterior.rds")) |>
  filter(model == "M1_status", term %in% TERMS) |>
  mutate(df = n - 5L)  # intercept + 2 status + log_mass + log_egg_size
post_all <- bind_rows(post_m1 |> select(any_of(names(post_new))), post_new)

post_sum <- post_all |>
  group_by(dataset, trait, model, term) |>
  summarise(post_std_median = median(std_estimate),
            post_prop_p_lt_005 = mean(p < 0.05),
            post_delta_min_median = median(delta_min(estimate, se, df, trait_sd_volant)),
            post_delta_min_max = max(delta_min(estimate, se, df, trait_sd_volant)),
            .groups = "drop")

robust <- mcc_fits |>
  left_join(post_sum, by = c("dataset", "trait", "model", "term")) |>
  mutate(delta_min_sd = delta_min(estimate, se, df, trait_sd_volant)) |>
  relocate(dataset, trait, model, term) |>
  mutate(across(where(is.double), ~ signif(.x, 4)))
write_table(robust, "robustness_covariates")
print(robust |> filter(term == "flight_statusterrestrial") |>
        select(dataset, trait, model, std_estimate, p, delta_min_sd,
               post_prop_p_lt_005, post_delta_min_max), n = Inf, width = Inf)

# ---- B. Equivalence against the HWI benchmark (Montgomerie) ---------------------
d_m <- DATA$montgomerie$d
# HWI slope among volant species only (avoids using flightless species to
# predict themselves); same covariates as Stoddard et al.'s model (M2).
volant <- d_m |> filter(flight_status == "volant")
hwi_slope <- map_dfr(DATA$montgomerie$traits, function(tr) {
  f <- fit_terms(volant, mcc, tr, "~ hwi + log_mass + log_egg_size", terms = "hwi")
  tibble(trait = tr, hwi_std_slope = f$estimate / f$trait_sd_volant)
})
# Observed HWI drop: sister-group mean HWI minus each terrestrial species' HWI
terr_pairs <- pairs |> filter(dataset == "montgomerie", locomotion == "terrestrial")
hwi_of <- setNames(d_m$hwi, d_m$birdtree_name)
hwi_drop <- map_dfr(seq_len(nrow(terr_pairs)), function(k) {
  fl  <- strsplit(terr_pairs$flightless_taxa[k], ", ")[[1]]
  sis <- strsplit(terr_pairs$sister_taxa[k], ", ")[[1]]
  tibble(origin = terr_pairs$origin_id[k], species = fl,
         drop = mean(hwi_of[sis], na.rm = TRUE) - hwi_of[fl])
})
mean_drop <- mean(hwi_drop$drop)
log_msg("Mean HWI drop over ", nrow(hwi_drop), " terrestrial species: ", round(mean_drop, 2))

equiv <- robust |>
  filter(dataset == "montgomerie", term == "flight_statusterrestrial",
         model %in% names(FORMULAS)) |>
  left_join(hwi_slope, by = "trait") |>
  mutate(mean_hwi_drop = mean_drop,
         benchmark_sd = hwi_std_slope * mean_hwi_drop,
         p_vs_benchmark = p_noninf(estimate, se, df, trait_sd_volant, benchmark_sd),
         p_vs_0.5sd = p_noninf(estimate, se, df, trait_sd_volant, 0.5),
         benchmark_excluded = delta_min_sd < benchmark_sd) |>
  select(trait, model, std_estimate, se, df, delta_min_sd, post_delta_min_median,
         post_delta_min_max, hwi_std_slope, mean_hwi_drop, benchmark_sd,
         p_vs_benchmark, p_vs_0.5sd, benchmark_excluded) |>
  mutate(across(where(is.double), ~ signif(.x, 4)))
write_table(equiv, "equivalence_tests")
print(equiv, width = Inf)

# ---- C. Montgomerie extrapolation diagnostics -----------------------------------
is_rat <- grepl(RATITE_RE, d_m$birdtree_name)
max_vol <- max(d_m$mass_g[d_m$flight_status == "volant"])
rat_pair <- terr_pairs |> filter(origin_id == "Ratites")
diag <- tibble(
  check = c("Heaviest volant species, mass (g)",
            "Lightest non-kiwi ratite, mass (g)",
            "Ratites heavier than heaviest volant species",
            "Ratites in dataset",
            paste("Mean", DATA$montgomerie$traits, "ratites"),
            paste("Mean", DATA$montgomerie$traits, "sister group (tinamous)")),
  value = c(max_vol,
            min(d_m$mass_g[is_rat & !grepl("^Apteryx", d_m$birdtree_name)]),
            sum(d_m$mass_g[is_rat] > max_vol),
            sum(is_rat),
            unlist(rat_pair[paste0("fl_", DATA$montgomerie$traits)]),
            unlist(rat_pair[paste0("sis_", DATA$montgomerie$traits)]))
) |> mutate(value = signif(value, 4))
write_table(diag, "montgomerie_diagnostics")
print(diag, width = Inf)
