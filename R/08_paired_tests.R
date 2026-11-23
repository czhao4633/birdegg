# 08_paired_tests.R
# Sister-pair tests on the Montgomerie pairs (H1, H2, H4).
#
# Difference per pair = mean(flightless clade) - mean(volant sister group),
# also expressed in SD units of the trait among volant species.
# H1 predicts terrestrial differences < 0 (rounder = lower elongation, more
# symmetric = lower pointedness / polar asymmetry). H2 predicts ~0 for diving.
#
# Tests (terrestrial pairs, n = 4; diving pairs, n = 2, are descriptive only):
#   - one-sample t-test and exact Wilcoxon signed-rank test (one-sided, < 0)
#   - exact sign-flip permutation test of the mean (2^n assignments)
#     NOTE: with n = 4 the smallest attainable one-sided p is 1/16 = 0.0625.
#   - age-matched phylogenetic randomisation: each observed difference is
#     compared with differences between all-volant sister clades on the MCC
#     tree whose age lies within a factor of 2 of the pair's divergence time;
#     the null distribution of the mean standardised difference is built by
#     drawing one matched volant contrast per pair (10,000 draws).
# H4: standardised difference ~ divergence time (exploratory; n = 4 and 6).

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
})
set.seed(SEED)

pairs <- readRDS(file.path(PATHS$processed, "sister_pairs.rds")) |>
  filter(include_paired_test)
d <- readRDS(file.path(PATHS$processed, "analysis_montgomerie.rds"))
traits <- c("elongation", "pointedness", "polar_asym")
sd_volant <- sapply(traits, function(t) sd(d[[t]][d$flight_status == "volant"]))

# ---- Pair differences ------------------------------------------------------------
diffs <- pairs |>
  select(origin_id, locomotion, divergence_myr, divergence_source,
         n_flightless, n_sister, starts_with("fl_"), starts_with("sis_")) |>
  pivot_longer(c(starts_with("fl_"), starts_with("sis_")),
               names_to = c("side", "trait"), names_pattern = "^(fl|sis)_(.*)$") |>
  pivot_wider(names_from = side, values_from = value) |>
  filter(trait %in% traits) |>   # drop Stoddard columns (NA for these rows)
  mutate(diff = fl - sis, diff_std = diff / sd_volant[trait])
stopifnot(!anyNA(diffs$diff))
write_table(diffs |> mutate(across(c(fl, sis, diff, diff_std), ~ signif(.x, 4))),
            "paired_differences")

# ---- Classical and exact tests -------------------------------------------------
sign_flip_p <- function(x) {
  # exact one-sided p for mean(x) < 0 over all 2^n sign assignments
  n <- length(x)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  null <- as.vector(signs %*% abs(x)) / n
  mean(null <= mean(x) + 1e-12)
}

tests <- diffs |>
  group_by(locomotion, trait) |>
  summarise(
    n_pairs = n(),
    mean_diff = mean(diff), mean_diff_std = mean(diff_std),
    n_negative = sum(diff < 0),
    t_p_one_sided = if (n() >= 3) t.test(diff, alternative = "less")$p.value else NA,
    wilcox_p_one_sided = if (n() >= 3)
      wilcox.test(diff, alternative = "less", exact = TRUE)$p.value else NA,
    signflip_p_one_sided = if (n() >= 3) sign_flip_p(diff) else NA,
    min_attainable_p = 1 / 2^n(),
    .groups = "drop"
  )

# ---- Age-matched phylogenetic randomisation -----------------------------------
mcc <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
tree <- keep.tip(mcc, d$birdtree_name)
vol <- d$birdtree_name[d$flight_status == "volant"]
depth <- node.depth.edgelength(tree)
height <- max(depth[seq_len(Ntip(tree))])
desc_tips <- function(node) {
  if (node <= Ntip(tree)) tree$tip.label[node] else extract.clade(tree, node)$tip.label
}
# Candidate contrasts: internal nodes whose two child clades are all volant
cand <- list()
for (node in (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode)) {
  ch <- tree$edge[tree$edge[, 1] == node, 2]
  if (length(ch) != 2) next
  a <- desc_tips(ch[1]); b <- desc_tips(ch[2])
  if (!all(c(a, b) %in% vol)) next
  cand[[length(cand) + 1]] <- list(age = height - depth[node], a = a, b = b)
}
cand_age <- sapply(cand, `[[`, "age")
log_msg("Volant sister-clade contrasts available: ", length(cand))

clade_mean <- function(tips, trait) mean(d[[trait]][match(tips, d$birdtree_name)])
null_contrast <- function(k, trait) {
  # random orientation makes the null symmetric around 0
  s <- sample(c(-1, 1), 1)
  s * (clade_mean(cand[[k]]$a, trait) - clade_mean(cand[[k]]$b, trait)) / sd_volant[trait]
}

N_RAND <- 10000L
phylo_rand <- list()
for (tr in traits) {
  for (loc in c("terrestrial", "diving", "all")) {
    sub <- diffs |> filter(trait == tr, loc == "all" | locomotion == loc)
    pools <- lapply(sub$divergence_myr, function(a) which(cand_age >= a / 2 & cand_age <= a * 2))
    # precompute every contrast in each pool once, then resample
    pool_vals <- lapply(pools, function(p) sapply(p, function(k) null_contrast(k, tr)))
    null_means <- replicate(N_RAND, mean(sapply(pool_vals, function(v) {
      sample(c(-1, 1), 1) * abs(v[sample.int(length(v), 1)])
    })))
    obs <- mean(sub$diff_std)
    phylo_rand[[length(phylo_rand) + 1]] <- tibble(
      trait = tr, locomotion = loc, n_pairs = nrow(sub),
      obs_mean_diff_std = obs,
      pool_sizes = paste(sapply(pools, length), collapse = "/"),
      p_one_sided_less = mean(null_means <= obs),
      p_two_sided = mean(abs(null_means) >= abs(obs))
    )
  }
}
phylo_rand <- bind_rows(phylo_rand)

# Per-pair percentile within its own age-matched null
pair_pct <- diffs |>
  rowwise() |>
  mutate(null_percentile = {
    p <- which(cand_age >= divergence_myr / 2 & cand_age <= divergence_myr * 2)
    v <- abs(sapply(p, function(k) null_contrast(k, trait)))
    # proportion of matched volant contrasts with |diff| >= |observed|
    mean(v >= abs(diff_std))
  },
  # 95th percentile of |difference| among matched volant contrasts: the null
  # is symmetric, so +/- this value is its central 90 % range (used in Fig. 3)
  null_abs_q90 = {
    p <- which(cand_age >= divergence_myr / 2 & cand_age <= divergence_myr * 2)
    unname(quantile(abs(sapply(p, function(k) null_contrast(k, trait))), 0.90))
  },
  n_null = length(which(cand_age >= divergence_myr / 2 & cand_age <= divergence_myr * 2))) |>
  ungroup() |>
  select(origin_id, locomotion, trait, divergence_myr, diff_std, n_null,
         null_percentile, null_abs_q90)

# ---- H4: difference vs divergence time --------------------------------------------
h4 <- diffs |>
  group_by(trait) |>
  group_modify(function(x, ...) {
    bind_rows(lapply(c("terrestrial", "all"), function(loc) {
      y <- if (loc == "all") x else filter(x, locomotion == loc)
      f <- lm(diff_std ~ divergence_myr, data = y)
      tibble(subset = loc, n = nrow(y),
             slope_per_10myr = 10 * coef(f)[2],
             slope_p = summary(f)$coefficients[2, 4],
             spearman_rho = cor(y$diff_std, y$divergence_myr, method = "spearman"))
    }))
  }) |>
  ungroup()

r4 <- function(df) mutate(df, across(where(is.double), ~ signif(.x, 3)))
write_table(r4(tests), "paired_tests")
write_table(r4(phylo_rand), "paired_phylo_randomisation")
write_table(r4(pair_pct), "paired_null_percentiles")
write_table(r4(h4), "h4_divergence_time")
print(r4(diffs |> select(origin_id, locomotion, trait, diff_std, divergence_myr)), n = Inf)
print(r4(tests), width = Inf); print(r4(phylo_rand), width = Inf)
print(r4(pair_pct), n = Inf); print(r4(h4), width = Inf)
