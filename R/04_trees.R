# 04_trees.R
# Prepare phylogenies from the BirdTree Hackett Stage 2 posterior (trees 1-1000):
#   1. prune every tree to the species present in either egg dataset
#   2. draw N_POSTERIOR_TREES trees at random (fixed seed) for robustness runs
#   3. compute a maximum clade credibility (MCC) tree from all 1,000 pruned
#      trees for the main analyses (BirdTree distributes no consensus tree)
# Outputs (data/processed/): trees_pruned_1000.rds, trees_posterior_100.rds,
# tree_mcc.rds

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
})

egg  <- readRDS(file.path(PATHS$processed, "egg_shape_birdtree.rds"))
mont <- readRDS(file.path(PATHS$processed, "egg_montgomerie_birdtree.rds"))
keep <- union(egg$birdtree_name, mont$birdtree_name)
log_msg("Species kept (union of egg datasets): ", length(keep))

tre_file <- file.path(PATHS$raw, "birdtree", "HackettStage2_0001_1000", "AllBirdsHackett1.tre")

# ---- 1. Prune (stream the 460 MB file one tree per line) --------------------
out_1000 <- file.path(PATHS$processed, "trees_pruned_1000.rds")
if (!file.exists(out_1000)) {
  con <- file(tre_file, "r")
  trees <- vector("list", 1000L)
  i <- 0L
  while (length(line <- readLines(con, n = 1L)) > 0) {
    if (!nzchar(line)) next
    i <- i + 1L
    trees[[i]] <- keep.tip(read.tree(text = line), keep)
    if (i %% 100 == 0) log_msg("Pruned ", i, " trees")
  }
  close(con)
  trees <- trees[seq_len(i)]
  class(trees) <- "multiPhylo"
  stopifnot(length(trees) == 1000L)
  saveRDS(trees, out_1000)
} else {
  trees <- readRDS(out_1000)
}

# ---- 2. Random posterior sample ---------------------------------------------
set.seed(SEED)
idx <- sort(sample(length(trees), N_POSTERIOR_TREES))
post <- trees[idx]
attr(post, "source_index") <- idx
saveRDS(post, file.path(PATHS$processed, "trees_posterior_100.rds"))
write.csv(data.frame(sampled_tree_index = idx),
          file.path(PATHS$tables, "posterior_tree_indices.csv"), row.names = FALSE)

# ---- 3. MCC tree ---------------------------------------------------------------
# phangorn::maxCladeCred returns the sampled tree with the highest product of
# clade frequencies; its branch lengths are those of that sampled tree, so the
# MCC tree stays ultrametric.
mcc <- maxCladeCred(trees)
cred <- maxCladeCred(trees, tree = FALSE)
attr(mcc, "mcc_index") <- which.max(cred)
log_msg("MCC tree = posterior tree #", which.max(cred),
        "; ultrametric: ", is.ultrametric(mcc), "; tips: ", Ntip(mcc))
saveRDS(mcc, file.path(PATHS$processed, "tree_mcc.rds"))
