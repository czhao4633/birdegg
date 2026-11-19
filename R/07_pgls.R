# 07_pgls.R
# Global phylogenetic regressions (phylolm, Pagel's lambda).
#
# Models, fitted for every shape trait in each dataset:
#   M1 status   : trait ~ flight_status + log_mass + log_egg_size
#   M2 hwi      : trait ~ hwi + log_mass + log_egg_size   (as in Stoddard 2017)
#   M3 confound : M1 + log_clutch + inc_site              (Montgomerie only)
#   M4 volume   : M1 with log(egg volume) as egg size     (Montgomerie only)
#   M5 no kiwi  : M1 without Apteryx (very large eggs relative to body mass)
# flight_status reference level = volant; the terrestrial-vs-diving contrast
# (H2) is computed from the M1 coefficients and their covariance.
#
# Robustness: M1 and M2 refitted on the 100 posterior trees.
# Outputs: output/tables/pgls_main.csv, pgls_posterior.csv

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  library(dplyr)
  library(purrr)
})

mcc  <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
post <- readRDS(file.path(PATHS$processed, "trees_posterior_100.rds"))

DATA <- list(
  montgomerie = list(d = readRDS(file.path(PATHS$processed, "analysis_montgomerie.rds")),
                     traits = c("elongation", "pointedness", "polar_asym")),
  stoddard    = list(d = readRDS(file.path(PATHS$processed, "analysis_stoddard.rds")),
                     traits = c("ellipticity", "asymmetry"))
)

FORMULAS <- list(
  M1_status   = "~ flight_status + log_mass + log_egg_size",
  M2_hwi      = "~ hwi + log_mass + log_egg_size",
  M3_confound = "~ flight_status + log_mass + log_egg_size + log_clutch + inc_site",
  M4_volume   = "~ flight_status + log_mass + log_egg_volume"
)

#' Fit one PGLS and return a tidy table (one row per coefficient), plus the
#' derived terrestrial - diving contrast when both levels are present.
fit_pgls <- function(d, tree, trait, rhs) {
  d <- as.data.frame(d)
  rownames(d) <- d$birdtree_name
  vars <- all.vars(as.formula(paste(trait, rhs)))
  d <- d[complete.cases(d[, vars]), , drop = FALSE]
  tr <- keep.tip(tree, d$birdtree_name)
  d <- d[tr$tip.label, , drop = FALSE]
  fit <- phylolm(as.formula(paste(trait, rhs)), data = d, phy = tr, model = "lambda")
  s <- summary(fit)$coefficients
  out <- tibble(term = rownames(s), estimate = s[, "Estimate"], se = s[, "StdErr"],
                t = s[, "t.value"], p = s[, "p.value"])
  b <- coef(fit); V <- vcov(fit)
  if (all(c("flight_statusterrestrial", "flight_statusdiving") %in% names(b))) {
    k <- c("flight_statusterrestrial", "flight_statusdiving")
    est <- b[k[1]] - b[k[2]]
    se <- sqrt(V[k[1], k[1]] + V[k[2], k[2]] - 2 * V[k[1], k[2]])
    dfree <- nrow(d) - length(b)
    out <- bind_rows(out, tibble(term = "terrestrial_minus_diving", estimate = est,
                                 se = se, t = est / se, p = 2 * pt(-abs(est / se), dfree)))
  }
  out |> mutate(lambda = fit$optpar, n = nrow(d),
                # effect in SD units of the trait among volant species
                trait_sd_volant = sd(d[[trait]][d$flight_status == "volant"]),
                std_estimate = estimate / trait_sd_volant)
}

# ---- Main fits on the MCC tree ------------------------------------------------
main <- list()
for (dn in names(DATA)) {
  for (tr in DATA[[dn]]$traits) {
    for (m in names(FORMULAS)) {
      if (dn == "stoddard" && m %in% c("M3_confound", "M4_volume")) next
      main[[length(main) + 1]] <- fit_pgls(DATA[[dn]]$d, mcc, tr, FORMULAS[[m]]) |>
        mutate(dataset = dn, trait = tr, model = m)
    }
    # M5: drop kiwi
    no_kiwi <- DATA[[dn]]$d |> filter(!grepl("^Apteryx_", birdtree_name))
    main[[length(main) + 1]] <- fit_pgls(no_kiwi, mcc, tr, FORMULAS$M1_status) |>
      mutate(dataset = dn, trait = tr, model = "M5_no_kiwi")
  }
}
main <- bind_rows(main) |>
  relocate(dataset, trait, model) |>
  mutate(across(c(estimate, se, std_estimate), ~ signif(.x, 4)),
         across(c(t, lambda), ~ round(.x, 3)), p = signif(p, 3))
write_table(main, "pgls_main")

focus <- c("flight_statusterrestrial", "flight_statusdiving",
           "terrestrial_minus_diving", "hwi")
print(main |> filter(term %in% focus) |>
        select(dataset, trait, model, term, estimate, std_estimate, se, p, lambda, n),
      n = Inf, width = Inf)

# ---- Posterior robustness -------------------------------------------------------
log_msg("Refitting M1 and M2 on ", length(post), " posterior trees ...")
post_res <- list()
for (dn in names(DATA)) {
  for (tr in DATA[[dn]]$traits) {
    for (m in c("M1_status", "M2_hwi")) {
      res <- map_dfr(seq_along(post), function(i) {
        fit_pgls(DATA[[dn]]$d, post[[i]], tr, FORMULAS[[m]]) |>
          filter(term %in% focus) |>
          mutate(tree = i)
      })
      post_res[[length(post_res) + 1]] <- res |> mutate(dataset = dn, trait = tr, model = m)
    }
    log_msg("  done: ", dn, " / ", tr)
  }
}
post_res <- bind_rows(post_res) |> relocate(dataset, trait, model, tree)
saveRDS(post_res, file.path(PATHS$processed, "pgls_posterior.rds"))

post_summary <- post_res |>
  group_by(dataset, trait, model, term) |>
  summarise(est_median = median(estimate),
            est_lo95 = quantile(estimate, 0.025), est_hi95 = quantile(estimate, 0.975),
            std_median = median(std_estimate),
            prop_p_lt_005 = mean(p < 0.05),
            prop_same_sign_as_median = mean(sign(estimate) == sign(median(estimate))),
            lambda_median = median(lambda),
            .groups = "drop") |>
  mutate(across(where(is.double), ~ signif(.x, 4)))
write_table(post_summary, "pgls_posterior_summary")
print(post_summary, n = Inf, width = Inf)
