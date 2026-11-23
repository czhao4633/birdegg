# 05_pairs.R
# For every independent loss of flight, find its sister group of
# volant species with egg data, and estimate the divergence time between them
# (an upper bound on how long the lineage has been flightless).
#
# Run separately for each egg dataset, because each prunes the tree
# differently and so can yield a different sister group.
#
# Definitions
#   origin        : maximal all-flightless clade on the MCC tree (00_setup.R)
#   sister group  : all other descendants of the origin's parent node, i.e.
#                   the closest relatives *among species with egg data*
#   stem age      : age of that parent node (Myr) = divergence time
#   crown age     : age of the MRCA of the flightless clade (if > 1 species)
# Robustness across the 100 posterior trees: monophyly of each flightless
# clade, identity of its sister group, distribution of stem ages.

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(readr)
})

mcc  <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
post <- readRDS(file.path(PATHS$processed, "trees_posterior_100.rds"))
fl   <- readRDS(file.path(PATHS$processed, "flightless.rds")) |>
  filter(locomotion %in% c("terrestrial", "diving"))

datasets <- list(
  montgomerie = list(
    data = readRDS(file.path(PATHS$processed, "egg_montgomerie_birdtree.rds")),
    traits = c("elongation", "pointedness", "polar_asym")
  ),
  stoddard = list(
    data = readRDS(file.path(PATHS$processed, "egg_shape_birdtree.rds")),
    traits = c("ellipticity", "asymmetry")
  )
)

node_age <- function(tree, node) {
  # Age of a node in an ultrametric tree (Myr before present). Computed from
  # root-to-node depths rather than branching.times(), whose names switch to
  # node labels when a tree carries them (the MCC tree does).
  depth <- node.depth.edgelength(tree)
  max(depth[seq_len(Ntip(tree))]) - depth[node]
}

#' Origin + sister group for one flightless clade on one tree.
pair_on_tree <- function(tree, clade_tips, flightless_tips) {
  node <- if (length(clade_tips) == 1) which(tree$tip.label == clade_tips) else
    getMRCA(tree, clade_tips)
  parent <- tree$edge[tree$edge[, 2] == node, 1]
  sib_tips <- setdiff(extract.clade(tree, parent)$tip.label, clade_tips)
  list(
    monophyletic = length(clade_tips) == 1 || is.monophyletic(tree, clade_tips),
    sister_tips = setdiff(sib_tips, flightless_tips),
    sister_has_flightless = any(sib_tips %in% flightless_tips),
    stem_age = node_age(tree, parent),
    crown_age = if (length(clade_tips) > 1) node_age(tree, node) else 0
  )
}

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# Species with genetic data = tips of the BirdTree Stage 1 trees. Stage 2
# species absent from Stage 1 were placed by taxonomic constraints only, so
# their sister relationships and divergence times are not data-informed.
stage1_file <- file.path(PATHS$raw, "birdtree", "HackettStage1_0001_1000",
                         "HackettStage1Full_1.tre")
genetic_tips <- read.tree(text = readLines(stage1_file, n = 1))$tip.label
n_imputed <- function(tips) sum(!tips %in% genetic_tips)

all_pairs <- list()
for (dname in names(datasets)) {
  d <- datasets[[dname]]
  spp <- d$data$birdtree_name
  tree <- keep.tip(mcc, spp)
  fl_d <- fl |> filter(birdtree_name %in% spp)
  clades <- flightless_clades(tree, fl_d$birdtree_name)
  log_msg(dname, ": ", length(clades), " origins on the MCC tree")

  rows <- lapply(seq_along(clades), function(k) {
    ct <- clades[[k]]
    loco <- unique(fl_d$locomotion[match(ct, fl_d$birdtree_name)])
    stopifnot(length(loco) == 1)  # an origin must not mix locomotion types
    p <- pair_on_tree(tree, ct, fl_d$birdtree_name)

    # Posterior robustness
    pp <- lapply(post, function(tr) {
      tr <- keep.tip(tr, spp)
      if (!(length(ct) == 1 || is.monophyletic(tr, ct))) {
        return(c(mono = 0, same_sister = NA, jac = NA, stem = NA))
      }
      q <- pair_on_tree(tr, ct, fl_d$birdtree_name)
      c(mono = 1, same_sister = setequal(q$sister_tips, p$sister_tips),
        jac = jaccard(q$sister_tips, p$sister_tips), stem = q$stem_age)
    })
    pp <- do.call(rbind, pp)

    trait_means <- function(tips) {
      x <- d$data[match(tips, d$data$birdtree_name), d$traits, drop = FALSE]
      colMeans(x, na.rm = TRUE)
    }
    tibble(
      dataset = dname,
      origin_id = NA_character_,
      locomotion = loco,
      n_flightless = length(ct),
      flightless_taxa = paste(sort(ct), collapse = ", "),
      n_sister = length(p$sister_tips),
      sister_taxa = paste(sort(p$sister_tips), collapse = ", "),
      sister_has_flightless = p$sister_has_flightless,
      n_flightless_imputed = n_imputed(ct),
      n_sister_imputed = n_imputed(p$sister_tips),
      stem_age_mcc = round(p$stem_age, 2),
      crown_age_mcc = round(p$crown_age, 2),
      stem_age_median = round(median(pp[, "stem"], na.rm = TRUE), 2),
      stem_age_lo95 = round(quantile(pp[, "stem"], 0.025, na.rm = TRUE), 2),
      stem_age_hi95 = round(quantile(pp[, "stem"], 0.975, na.rm = TRUE), 2),
      prop_monophyletic = mean(pp[, "mono"]),
      prop_same_sister = mean(pp[, "same_sister"], na.rm = TRUE),
      mean_sister_jaccard = round(mean(pp[, "jac"], na.rm = TRUE), 2),
      !!!setNames(as.list(trait_means(ct)), paste0("fl_", d$traits)),
      !!!setNames(as.list(trait_means(p$sister_tips)), paste0("sis_", d$traits))
    )
  })
  all_pairs[[dname]] <- bind_rows(rows)
}

# ---- Human-readable origin labels --------------------------------------------
label_origin <- function(taxa) {
  case_when(
    grepl("Struthio|Rhea|Dromaius|Casuarius|Apteryx", taxa) ~ "Ratites",
    grepl("Spheniscus|Pygoscelis|Eudyptes|Aptenodytes", taxa) ~ "Penguins",
    TRUE ~ sub(",.*$", "", taxa)
  )
}
pairs <- bind_rows(all_pairs) |>
  mutate(origin_id = label_origin(flightless_taxa))

# ---- Comparison with Saitta et al. pairs -------------------------------------
sp <- read.csv(file.path(PATHS$raw, "saitta2025", "extracted", "sisterpairs.csv"),
               fileEncoding = "UTF-8-BOM")
fixes <- read_csv(file.path(PATHS$curated, "name_fixes.csv"),
                  show_col_types = FALSE, comment = "#") |>
  filter(dataset == "saitta2025")
sayol <- readRDS(file.path(PATHS$processed, "sayol_birdtree.rds"))
sp <- sp |>
  mutate(
    binomial = to_birdtree_name(sub("^(\\S+)\\s+(\\S+).*$", "\\1 \\2",
                                    gsub("\\(=[^)]*\\)\\s*", "", Taxon))),
    binomial = coalesce(to_birdtree_name(fixes$birdtree_name)[
      match(binomial, to_birdtree_name(fixes$original_name))], binomial),
    birdtree_name = coalesce(sayol$birdtree_name[
      match(binomial, to_birdtree_name(sayol$species))], binomial),
    max_age_myr = max_age / 1e6
  )
saitta_fl <- sp |> filter(flighted == 0)
saitta_partner <- sp |> filter(flighted == 1, !is.na(pair)) |>
  select(pair, saitta_partner = birdtree_name)

pairs$saitta_flightless <- NA_character_
pairs$saitta_partner <- NA_character_
pairs$saitta_max_age_myr <- NA_real_
for (i in seq_len(nrow(pairs))) {
  hit <- saitta_fl |>
    filter(birdtree_name %in% strsplit(pairs$flightless_taxa[i], ", ")[[1]])
  if (nrow(hit) == 0) next
  pairs$saitta_flightless[i] <- paste(hit$birdtree_name, collapse = ", ")
  partners <- saitta_partner$saitta_partner[match(hit$pair, saitta_partner$pair)]
  pairs$saitta_partner[i] <- paste(unique(na.omit(partners)), collapse = ", ")
  pairs$saitta_max_age_myr[i] <- max(hit$max_age_myr, na.rm = TRUE)
}
pairs$saitta_max_age_myr[is.infinite(pairs$saitta_max_age_myr)] <- NA

# ---- Literature overrides for taxonomy-imputed species -----------------------
pairs$pair_basis <- "BirdTree MCC"
ov <- read_csv(file.path(PATHS$curated, "pair_overrides.csv"),
               show_col_types = FALSE, comment = "#")
for (k in seq_len(nrow(ov))) {
  i <- which(pairs$dataset == ov$dataset[k] & pairs$flightless_taxa == ov$flightless[k])
  stopifnot(length(i) == 1)
  d <- datasets[[ov$dataset[k]]]
  stopifnot(ov$partner[k] %in% d$data$birdtree_name)
  x <- d$data[d$data$birdtree_name == ov$partner[k], d$traits]
  pairs[i, paste0("sis_", d$traits)] <- as.list(x)
  pairs$sister_taxa[i] <- ov$partner[k]
  pairs$n_sister[i] <- 1L
  pairs$n_sister_imputed[i] <- n_imputed(ov$partner[k])
  pairs$sister_has_flightless[i] <- FALSE
  # BirdTree ages are uninformative for imputed taxa: use the literature value
  pairs$stem_age_mcc[i] <- ov$divergence_myr[k]
  pairs[i, c("stem_age_median", "stem_age_lo95", "stem_age_hi95",
             "prop_same_sister", "mean_sister_jaccard")] <- NA
  pairs$pair_basis[i] <- paste0("Literature: ", ov$source[k])
}

# Divergence time used for H4: posterior median for BirdTree pairs (the MCC
# tree is a single sample, e.g. ratites 65.9 Myr vs median 80.7 Myr), the
# literature value for overridden pairs.
pairs <- pairs |>
  mutate(divergence_myr = if_else(pair_basis == "BirdTree MCC",
                                  stem_age_median, stem_age_mcc),
         divergence_source = if_else(pair_basis == "BirdTree MCC",
                                     "BirdTree posterior median (100 trees)",
                                     "Saitta et al. 2025"))

# ---- Which pairs enter the paired test ---------------------------------------
# Excluded: imputed flightless species with no literature partner in the data
# (Atlantisia), meaningless sister groups (Stoddard kagu: 42 Strisores spp.),
# and every Stoddard pair, because 2 valid pairs are too few for a test.
pairs <- pairs |>
  mutate(
    include_paired_test = case_when(
      dataset == "stoddard" ~ FALSE,
      n_flightless_imputed > 0 & pair_basis == "BirdTree MCC" ~ FALSE,
      TRUE ~ TRUE
    ),
    exclusion_reason = case_when(
      dataset == "stoddard" & origin_id == "Rhynochetos_jubatus" ~
        "Sister group is 42 Strisores spp. (no sunbittern in data)",
      dataset == "stoddard" ~ "Stoddard: only 2 valid pairs, no paired test",
      !include_paired_test ~ "Imputed in BirdTree; no literature partner with egg data",
      TRUE ~ ""
    )
  )

saveRDS(pairs, file.path(PATHS$processed, "sister_pairs.rds"))
write_table(pairs, "stage2_sister_pairs")
print(pairs |> select(dataset, origin_id, locomotion, n_flightless, n_sister,
                      stem_age_mcc, stem_age_median, stem_age_lo95, stem_age_hi95,
                      prop_monophyletic, prop_same_sister, saitta_max_age_myr),
      n = Inf, width = Inf)
