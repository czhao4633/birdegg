# 09_rates.R
# H3: do egg-shape traits evolve faster in flightless lineages?
# Multi-rate Brownian motion with phytools::brownie.lite (O'Meara et al. 2006).
#
# Regime mapping: each flightless origin (maximal flightless clade) is painted
# from the midpoint of its stem branch (stem = 0.5), because flight was lost
# at an unknown point along the stem; everything else is "volant".
# Models compared with the single-rate model by likelihood-ratio test:
#   R2 : volant vs flightless
#   R3 : volant vs terrestrial vs diving
# Fitted on the MCC tree and on the 100 posterior trees, for both datasets.
# Traits are standardised (mean 0, SD 1) before fitting; rates are therefore
# in SD^2 per Myr, and only rate *ratios* are interpreted.
#
# Caveat: the flightless regime is dominated by the
# two species-rich, old clades (penguins, ratites); single-species origins
# contribute half a stem branch plus a tip branch each.

# OpenBLAS reads its thread count only at start-up. Parallel workers with
# multithreaded BLAS oversubscribe the CPU, and
# RhpcBLASctl cannot control Homebrew's OpenBLAS on macOS. So, unless the
# limit is already set, rerun this script in a child R with it set.
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

#' Paint flightless regimes on a tree pruned to the dataset's species.
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

fit_rates <- function(tree, d, trait, three) {
  st <- paint_regimes(tree, d, three)
  # Standardise the trait: rates are tiny on the raw scale and the optimiser
  # then hits a singular covariance matrix. Rate ratios are scale-invariant.
  x <- setNames(as.numeric(scale(d[[trait]])), d$birdtree_name)[st$tip.label]
  f <- brownie.lite(st, x)
  rates <- f$sig2.multiple
  tibble(
    regimes = if (three) "R3" else "R2",
    sig2_single = f$sig2.single,
    sig2_volant = rates["volant"],
    sig2_flightless = if (!three) rates["flightless"] else NA_real_,
    sig2_terrestrial = if (three) rates["terrestrial"] else NA_real_,
    sig2_diving = if (three) rates["diving"] else NA_real_,
    logL_single = f$logL1, logL_multi = f$logL.multiple,
    lrt_p = f$P.chisq, converged = f$convergence == "Optimization has converged."
  )
}

run_all <- function(tree) {
  map_dfr(names(DATA), function(dn) {
    map_dfr(DATA[[dn]]$traits, function(tr) {
      bind_rows(fit_rates(tree, DATA[[dn]]$d, tr, three = FALSE),
                fit_rates(tree, DATA[[dn]]$d, tr, three = TRUE)) |>
        mutate(dataset = dn, trait = tr, .before = 1)
    })
  }) |>
    mutate(ratio_flightless = sig2_flightless / sig2_volant,
           ratio_terrestrial = sig2_terrestrial / sig2_volant,
           ratio_diving = sig2_diving / sig2_volant)
}

main <- run_all(mcc)
write_table(main |> mutate(across(where(is.double), ~ signif(.x, 4))), "rates_main")
print(main |> select(dataset, trait, regimes, ratio_flightless, ratio_terrestrial,
                     ratio_diving, lrt_p, converged), width = Inf)

log_msg("Refitting on ", length(post), " posterior trees ...")
# Parallel over trees (fork-based; deterministic: no random numbers are drawn).
# BLAS is single-threaded here (see top of script). Results are cached per
# tree so an interrupted run resumes.
cache_dir <- file.path(PATHS$processed, "rates_posterior_cache")
dir.create(cache_dir, showWarnings = FALSE)
progress_file <- file.path(cache_dir, "progress.log")
n_cores <- max(1L, parallel::detectCores() - 1L)
invisible(parallel::mclapply(seq_along(post), function(i) {
  f <- file.path(cache_dir, sprintf("tree_%03d.rds", i))
  if (file.exists(f)) return(NULL)
  saveRDS(run_all(post[[i]]) |> mutate(tree = i), f)
  cat(format(Sys.time(), "%H:%M:%S"), "tree", i, "done\n",
      file = progress_file, append = TRUE)
  NULL
}, mc.cores = n_cores, mc.preschedule = FALSE))
post_res <- bind_rows(lapply(
  file.path(cache_dir, sprintf("tree_%03d.rds", seq_along(post))), readRDS))
stopifnot(n_distinct(post_res$tree) == length(post))
saveRDS(post_res, file.path(PATHS$processed, "rates_posterior.rds"))

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
# Only fits whose optimiser converged are summarised (2-8 % per model failed).
post_summary <- post_res |>
  group_by(dataset, trait, regimes) |>
  mutate(n_converged = sum(converged)) |>
  filter(converged) |>
  summarise(
    n_converged = first(n_converged),
    across(c(ratio_flightless, ratio_terrestrial, ratio_diving),
           list(median = ~ median(.x, na.rm = TRUE), lo95 = ~ q(.x, 0.025),
                hi95 = ~ q(.x, 0.975))),
    prop_lrt_p_lt_005 = mean(lrt_p < 0.05),
    .groups = "drop"
  ) |>
  mutate(across(where(is.double), ~ signif(.x, 3)))
write_table(post_summary, "rates_posterior_summary")
print(post_summary, width = Inf)
