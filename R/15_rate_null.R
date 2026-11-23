# 15_rate_null.R
# Is "asymmetry evolves more slowly in flightless lineages" (H3, R/09) an
# artefact of estimating a rate from very few lineages?
#
# For each trait, data are simulated on the MCC tree with ONE rate for all
# lineages (no true difference), then refitted exactly as in R/09
# (phytools::brownie.lite, flightless regimes painted from stem midpoints,
# R2 = volant vs flightless, R3 = volant / terrestrial / diving).
# Two null generators:
#   BM     : single-rate Brownian motion with the fitted single-rate sigma2
#   lambda : Pagel's-lambda model fitted to the trait (intercept only), i.e.
#            BM plus extra tip variance, which brownie.lite does not model;
#            checks whether that misspecification alone produces low ratios.
# Reported per trait: where the observed MCC ratio falls in the null
# distribution (one-sided p = share of null ratios <= observed), the null
# median ratio, and the null false-positive rate of the LRT.
#
# Seeds are per replicate (SEED + offset). Results are cached per
# trait x generator so an interrupted run resumes.
# Outputs: output/tables/rate_null.csv,
#          data/processed/rate_null_cache/*.rds

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
  library(parallel)
})

N_REP <- as.integer(Sys.getenv("RATE_NULL_NREP", "100"))  # env var only for smoke tests
N_CORES <- max(1L, parallel::detectCores() - 1L)

mcc <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
DATA <- list(
  montgomerie = list(d = readRDS(file.path(PATHS$processed, "analysis_montgomerie.rds")),
                     traits = c("elongation", "pointedness", "polar_asym")),
  stoddard    = list(d = readRDS(file.path(PATHS$processed, "analysis_stoddard.rds")),
                     traits = c("ellipticity", "asymmetry"))
)
observed <- read.csv(file.path(PATHS$tables, "rates_main.csv"))

# Same regime painting as R/09
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

#' Pagel's lambda transform of an ultrametric tree: internal branches x lambda,
#' tip branches extended so every tip keeps its original depth.
lambda_tree <- function(tree, lambda) {
  H <- max(node.depth.edgelength(tree))
  is_tip <- tree$edge[, 2] <= Ntip(tree)
  tree$edge.length[!is_tip] <- lambda * tree$edge.length[!is_tip]
  tree$edge.length[is_tip] <- lambda * tree$edge.length[is_tip] + (1 - lambda) * H
  tree
}

cache_dir <- file.path(PATHS$processed, "rate_null_cache")
dir.create(cache_dir, showWarnings = FALSE)

jobs <- do.call(rbind, lapply(names(DATA), function(dn)
  expand.grid(dataset = dn, trait = DATA[[dn]]$traits, generator = c("BM", "lambda"),
              stringsAsFactors = FALSE)))

res <- list()
for (j in seq_len(nrow(jobs))) {
  dn <- jobs$dataset[j]; tr <- jobs$trait[j]; gen <- jobs$generator[j]
  f <- file.path(cache_dir, sprintf("%s_%s_%s_n%d.rds", dn, tr, gen, N_REP))
  if (!file.exists(f)) {
    d <- DATA[[dn]]$d
    st2 <- paint_regimes(mcc, d, three = FALSE)
    st3 <- paint_regimes(mcc, d, three = TRUE)
    tree <- keep.tip(mcc, d$birdtree_name); tree$node.label <- NULL
    x <- setNames(as.numeric(scale(d[[tr]])), d$birdtree_name)[tree$tip.label]
    # Generating model fitted to the observed (standardised) trait
    sig2 <- brownie.lite(st2, x[st2$tip.label])$sig2.single
    gtree <- tree
    if (gen == "lambda") {
      lam <- phylolm(y ~ 1, data = data.frame(y = x, row.names = names(x)),
                     phy = tree, model = "lambda")
      gtree <- lambda_tree(tree, lam$optpar)
      sig2 <- lam$sigma2
    }
    log_msg("Null ", dn, " / ", tr, " / ", gen, ": ", N_REP, " replicates")
    seed_off <- 100000L * j
    reps <- mclapply(seq_len(N_REP), function(r) {
      set.seed(SEED + seed_off + r)
      y <- fastBM(gtree, sig2 = sig2)[st2$tip.label]
      f2 <- brownie.lite(st2, y); f3 <- brownie.lite(st3, y[st3$tip.label])
      data.frame(rep = r,
                 ratio_flightless = f2$sig2.multiple["flightless"] / f2$sig2.multiple["volant"],
                 lrt_p_R2 = f2$P.chisq, conv_R2 = f2$convergence == "Optimization has converged.",
                 ratio_terrestrial = f3$sig2.multiple["terrestrial"] / f3$sig2.multiple["volant"],
                 ratio_diving = f3$sig2.multiple["diving"] / f3$sig2.multiple["volant"],
                 conv_R3 = f3$convergence == "Optimization has converged.")
    }, mc.cores = N_CORES)
    saveRDS(bind_rows(reps) |> mutate(dataset = dn, trait = tr, generator = gen), f)
  }
  res[[j]] <- readRDS(f)
}
sims <- bind_rows(res)

obs <- observed |>
  select(dataset, trait, regimes, ratio_flightless, ratio_terrestrial, ratio_diving)
obs2 <- obs |> filter(regimes == "R2") |> select(dataset, trait, obs_flightless = ratio_flightless)
obs3 <- obs |> filter(regimes == "R3") |>
  select(dataset, trait, obs_terrestrial = ratio_terrestrial, obs_diving = ratio_diving)

summ <- sims |>
  left_join(obs2, by = c("dataset", "trait")) |>
  left_join(obs3, by = c("dataset", "trait")) |>
  group_by(dataset, trait, generator) |>
  summarise(
    n_rep = n(),
    obs_flightless = first(obs_flightless),
    null_flightless_median = median(ratio_flightless[conv_R2]),
    null_flightless_lo95 = quantile(ratio_flightless[conv_R2], 0.025),
    p_flightless_le_obs = mean(ratio_flightless[conv_R2] <= first(obs_flightless)),
    null_lrt_false_pos = mean(lrt_p_R2[conv_R2] < 0.05),
    obs_terrestrial = first(obs_terrestrial),
    null_terrestrial_median = median(ratio_terrestrial[conv_R3]),
    p_terrestrial_le_obs = mean(ratio_terrestrial[conv_R3] <= first(obs_terrestrial)),
    n_converged_R2 = sum(conv_R2), n_converged_R3 = sum(conv_R3),
    .groups = "drop") |>
  mutate(across(where(is.double), ~ signif(.x, 3)))
write_table(summ, "rate_null")
print(summ, width = Inf, n = Inf)
