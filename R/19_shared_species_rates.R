# 19_shared_species_rates.R
# Shared-species check of the rate results: in R/09 the flightless rate
# ratio for Montgomerie elongation (and pointedness) is > 1, but for the
# closely related Stoddard ellipticity it is < 1. Neither survived the
# simulation null of R/15, but the opposite directions need explaining.
#
# Three candidate causes:
#   1. species composition (the two datasets contain different flightless
#      lineages);
#   2. the metrics themselves (elongation and ellipticity are correlated,
#      r = 0.87, but not identical);
#   3. noise (ratios estimated from very few, short regimes).
#
# Approach, all on the MCC tree with the regime painting of R/09:
#   A. refit R2 / R3 on the species shared by both datasets, so that the
#      tree, species and regimes are identical and only the metric differs;
#   B. leave one flightless origin out at a time (shared species, R2);
#   C. lambda null (as in R/15) on the shared species, to see whether the
#      shared-species ratios are distinguishable from a single rate.
# Outputs: output/tables/shared_rates.csv, shared_rates_loo.csv,
#          shared_rates_null.csv, shared_traits_cor.csv;
#          data/processed/shared_rates_null_cache/*.rds

# Single-threaded BLAS in parallel workers (same reason as in R/09).
if (Sys.getenv("OPENBLAS_NUM_THREADS") != "1") {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  status <- system2(file.path(R.home("bin"), "Rscript"), script,
                    env = c("OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1"))
  quit(save = "no", status = status)
}

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(phylolm)
  library(dplyr)
  library(purrr)
  library(parallel)
})

N_REP <- as.integer(Sys.getenv("RATE_NULL_NREP", "100"))  # env var only for smoke tests
N_CORES <- max(1L, parallel::detectCores() - 1L)

mcc  <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
mont <- readRDS(file.path(PATHS$processed, "analysis_montgomerie.rds"))
stod <- readRDS(file.path(PATHS$processed, "analysis_stoddard.rds"))

# ---- Shared species ------------------------------------------------------------
shared <- inner_join(
  mont |> select(birdtree_name, flight_status_m = flight_status, flightless_m = flightless,
                 elongation, pointedness, polar_asym),
  stod |> select(birdtree_name, flight_status_s = flight_status, flightless_s = flightless,
                 ellipticity, asymmetry),
  by = "birdtree_name")
# Both datasets take flight status from the same source, so it must agree.
stopifnot(identical(as.character(shared$flight_status_m), as.character(shared$flight_status_s)))
shared <- shared |>
  rename(flight_status = flight_status_m, flightless = flightless_m) |>
  select(-flight_status_s, -flightless_s)
log_msg("Shared species: ", nrow(shared), " (", sum(shared$flightless), " flightless)")

TRAITS <- tribble(
  ~dataset,      ~trait,
  "montgomerie", "elongation",
  "montgomerie", "pointedness",
  "stoddard",    "ellipticity"
)

# ---- Helpers (same regime painting and fit as R/09) --------------------------
paint_regimes <- function(tree, d, three = FALSE) {
  tree <- keep.tip(tree, d$birdtree_name)
  tree$node.label <- NULL
  st <- paintSubTree(tree, node = Ntip(tree) + 1, state = "volant")
  fl_tips <- d$birdtree_name[d$flightless]
  for (cl in flightless_clades(tree, fl_tips)) {
    node <- if (length(cl) == 1) which(tree$tip.label == cl) else getMRCA(tree, cl)
    state <- if (three) as.character(d$flight_status[match(cl[1], d$birdtree_name)]) else "flightless"
    st <- paintSubTree(st, node = node, state = state, stem = 0.5)
  }
  st
}

fit_one <- function(d, trait, three) {
  st <- paint_regimes(mcc, d, three)
  x <- setNames(as.numeric(scale(d[[trait]])), d$birdtree_name)[st$tip.label]
  f <- brownie.lite(st, x)
  r <- f$sig2.multiple
  tibble(regimes = if (three) "R3" else "R2",
         ratio_flightless  = if (!three) r["flightless"] / r["volant"] else NA_real_,
         ratio_terrestrial = if (three) r["terrestrial"] / r["volant"] else NA_real_,
         ratio_diving      = if (three) r["diving"] / r["volant"] else NA_real_,
         lrt_p = f$P.chisq,
         converged = f$convergence == "Optimization has converged.")
}

# ---- A. Full data vs shared species ---------------------------------------------
full_d <- list(montgomerie = mont, stoddard = stod)
fits <- pmap_dfr(TRAITS, function(dataset, trait) {
  bind_rows(
    bind_rows(fit_one(full_d[[dataset]], trait, FALSE), fit_one(full_d[[dataset]], trait, TRUE)) |>
      mutate(species_set = "full", n_species = nrow(full_d[[dataset]])),
    bind_rows(fit_one(shared, trait, FALSE), fit_one(shared, trait, TRUE)) |>
      mutate(species_set = "shared", n_species = nrow(shared))
  ) |> mutate(dataset = dataset, trait = trait, .before = 1)
})
write_table(fits |> mutate(across(where(is.double), ~ signif(.x, 3))), "shared_rates")
print(fits, width = Inf)

# Correlation of the metrics on the shared species
cor_tbl <- bind_rows(
  tibble(subset = "all shared", n = nrow(shared),
         r_elong_ellip = cor(shared$elongation, shared$ellipticity),
         r_point_ellip = cor(shared$pointedness, shared$ellipticity)),
  tibble(subset = "flightless shared", n = sum(shared$flightless),
         r_elong_ellip = cor(shared$elongation[shared$flightless], shared$ellipticity[shared$flightless]),
         r_point_ellip = cor(shared$pointedness[shared$flightless], shared$ellipticity[shared$flightless]))
) |> mutate(across(where(is.double), ~ signif(.x, 3)))
write_table(cor_tbl, "shared_traits_cor")
print(cor_tbl)

# ---- B. Leave one flightless origin out (shared species, R2) ---------------------
tree_sh <- keep.tip(mcc, shared$birdtree_name)
origins <- flightless_clades(tree_sh, shared$birdtree_name[shared$flightless])
origin_label <- function(cl) {
  if (length(cl) == 1) return(cl)
  paste0(cl[1], " + ", length(cl) - 1, " more")
}
loo <- map_dfr(seq_along(origins), function(k) {
  d_k <- shared |> filter(!birdtree_name %in% origins[[k]])
  pmap_dfr(TRAITS, function(dataset, trait) {
    fit_one(d_k, trait, FALSE) |>
      mutate(dataset = dataset, trait = trait, .before = 1)
  }) |>
    mutate(dropped_origin = origin_label(origins[[k]]),
           dropped_n = length(origins[[k]]),
           dropped_status = as.character(shared$flight_status[match(origins[[k]][1], shared$birdtree_name)]),
           .before = 1)
})
write_table(loo |> select(-ratio_terrestrial, -ratio_diving) |>
              mutate(across(where(is.double), ~ signif(.x, 3))), "shared_rates_loo")
print(loo |> select(dropped_origin, dropped_status, trait, ratio_flightless, converged), n = Inf)

# ---- C. Lambda null on the shared species (R2) -------------------------------------
lambda_tree <- function(tree, lambda) {   # as in R/15
  H <- max(node.depth.edgelength(tree))
  is_tip <- tree$edge[, 2] <= Ntip(tree)
  tree$edge.length[!is_tip] <- lambda * tree$edge.length[!is_tip]
  tree$edge.length[is_tip] <- lambda * tree$edge.length[is_tip] + (1 - lambda) * H
  tree
}

cache_dir <- file.path(PATHS$processed, "shared_rates_null_cache")
dir.create(cache_dir, showWarnings = FALSE)
st2 <- paint_regimes(mcc, shared, three = FALSE)
tree_sh$node.label <- NULL

null_res <- map_dfr(seq_len(nrow(TRAITS)), function(j) {
  tr <- TRAITS$trait[j]
  f <- file.path(cache_dir, sprintf("%s_lambda_n%d.rds", tr, N_REP))
  if (!file.exists(f)) {
    x <- setNames(as.numeric(scale(shared[[tr]])), shared$birdtree_name)[tree_sh$tip.label]
    lam <- phylolm(y ~ 1, data = data.frame(y = x, row.names = names(x)),
                   phy = tree_sh, model = "lambda")
    gtree <- lambda_tree(tree_sh, lam$optpar)
    log_msg("Shared-species lambda null: ", tr, " (lambda = ", signif(lam$optpar, 3),
            "), ", N_REP, " replicates")
    seed_off <- 1900000L + 100000L * j
    reps <- mclapply(seq_len(N_REP), function(r) {
      set.seed(SEED + seed_off + r)
      y <- fastBM(gtree, sig2 = lam$sigma2)[st2$tip.label]
      f2 <- brownie.lite(st2, y)
      data.frame(rep = r,
                 ratio_flightless = f2$sig2.multiple["flightless"] / f2$sig2.multiple["volant"],
                 lrt_p = f2$P.chisq,
                 converged = f2$convergence == "Optimization has converged.",
                 lambda = lam$optpar)
    }, mc.cores = N_CORES)
    saveRDS(bind_rows(reps) |> mutate(trait = tr), f)
  }
  readRDS(f)
})

obs_sh <- fits |> filter(species_set == "shared", regimes == "R2") |>
  select(trait, obs_ratio = ratio_flightless)
null_summ <- null_res |>
  filter(converged) |>
  left_join(obs_sh, by = "trait") |>
  group_by(trait) |>
  summarise(n_converged = n(),
            lambda = first(lambda),
            obs_ratio = first(obs_ratio),
            null_median = median(ratio_flightless),
            null_lo95 = quantile(ratio_flightless, 0.025),
            null_hi95 = quantile(ratio_flightless, 0.975),
            p_le_obs = mean(ratio_flightless <= first(obs_ratio)),
            p_ge_obs = mean(ratio_flightless >= first(obs_ratio)),
            null_lrt_false_pos = mean(lrt_p < 0.05),
            .groups = "drop") |>
  mutate(across(where(is.double), ~ signif(.x, 3)))
write_table(null_summ, "shared_rates_null")
print(null_summ, width = Inf)
