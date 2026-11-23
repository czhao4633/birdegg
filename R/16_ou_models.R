# 16_ou_models.R
# Does flight loss shift the OPTIMUM of egg shape?
# Reversion, stated as an Ornstein-Uhlenbeck model: flightless lineages are
# pulled towards a different (rounder) optimum theta than volant lineages.
#
# Models (mvMORPH::mvOU / mvBM, MCC tree, regimes painted as in R/09):
#   BM1  : one-rate Brownian motion
#   OU1  : one optimum for all lineages
#   OUM2 : separate optima for volant and flightless lineages
#   OUM3 : separate optima for volant, terrestrial and diving lineages
# Traits are standardised (mean 0, SD 1), so differences between optima are
# in SD units of the whole dataset. Reversion predicts
#   theta_terrestrial - theta_volant < 0.
# Also reported: alpha as a phylogenetic half-life (ln 2 / alpha, Myr) and
# AICc weights. Deterministic (no random numbers).
# Output: output/tables/ou_models.csv

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
  library(mvMORPH)
  library(dplyr)
  library(parallel)
})

mcc <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
DATA <- list(
  montgomerie = list(d = readRDS(file.path(PATHS$processed, "analysis_montgomerie.rds")),
                     traits = c("elongation", "pointedness", "polar_asym")),
  stoddard    = list(d = readRDS(file.path(PATHS$processed, "analysis_stoddard.rds")),
                     traits = c("ellipticity", "asymmetry"))
)

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

fit_one <- function(dn, tr) {
  d <- DATA[[dn]]$d
  st2 <- paint_regimes(mcc, d, FALSE)
  st3 <- paint_regimes(mcc, d, TRUE)
  x <- setNames(as.numeric(scale(d[[tr]])), d$birdtree_name)[st2$tip.label]
  q <- list(echo = FALSE, diagnostic = FALSE)
  fits <- list(
    BM1  = do.call(mvBM,  c(list(tree = st2, data = x, model = "BM1"), q)),
    OU1  = do.call(mvOU,  c(list(tree = st2, data = x, model = "OU1"), q)),
    OUM2 = do.call(mvOU,  c(list(tree = st2, data = x, model = "OUM",
                                 param = list(root = FALSE)), q)),
    OUM3 = do.call(mvOU,  c(list(tree = st3, data = x[st3$tip.label], model = "OUM",
                                 param = list(root = FALSE)), q))
  )
  aicc <- sapply(fits, function(f) f$AICc)
  w <- exp(-0.5 * (aicc - min(aicc))); w <- w / sum(w)
  th2 <- setNames(as.vector(fits$OUM2$theta), rownames(fits$OUM2$theta))
  th3 <- setNames(as.vector(fits$OUM3$theta), rownames(fits$OUM3$theta))
  a2 <- as.numeric(fits$OUM2$alpha); a3 <- as.numeric(fits$OUM3$alpha)
  tibble(
    dataset = dn, trait = tr,
    aicc_BM1 = aicc["BM1"], aicc_OU1 = aicc["OU1"],
    aicc_OUM2 = aicc["OUM2"], aicc_OUM3 = aicc["OUM3"],
    w_BM1 = w["BM1"], w_OU1 = w["OU1"], w_OUM2 = w["OUM2"], w_OUM3 = w["OUM3"],
    lrt_p_OUM2_vs_OU1 = pchisq(2 * (fits$OUM2$LogLik - fits$OU1$LogLik), 1, lower.tail = FALSE),
    lrt_p_OUM3_vs_OU1 = pchisq(2 * (fits$OUM3$LogLik - fits$OU1$LogLik), 2, lower.tail = FALSE),
    theta_flightless_minus_volant = th2["flightless"] - th2["volant"],
    theta_terrestrial_minus_volant = th3["terrestrial"] - th3["volant"],
    theta_diving_minus_volant = th3["diving"] - th3["volant"],
    half_life_myr_OUM2 = log(2) / a2, half_life_myr_OUM3 = log(2) / a3,
    tree_height_myr = max(node.depth.edgelength(st2))
  )
}

jobs <- do.call(rbind, lapply(names(DATA), function(dn)
  data.frame(dataset = dn, trait = DATA[[dn]]$traits, stringsAsFactors = FALSE)))
res <- bind_rows(mclapply(seq_len(nrow(jobs)), function(j)
  fit_one(jobs$dataset[j], jobs$trait[j]),
  mc.cores = min(nrow(jobs), max(1L, parallel::detectCores() - 1L))))
res <- res |> mutate(across(where(is.double), ~ signif(.x, 4)))
write_table(res, "ou_models")
print(as.data.frame(res))
