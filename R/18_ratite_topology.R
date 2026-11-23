# 18_ratite_topology.R
# Sensitivity to the ratite topology.
# BirdTree's Hackett backbone has ratites monophyletic, sister to tinamous, so
# ratites count as ONE flight loss. Nuclear-genome studies instead place
# tinamous inside the ratites, with the ostrich sister to all other
# palaeognaths (Harshman et al. 2008, PNAS 105:13462, doi:10.1073/pnas.0803242105;
# Mitchell et al. 2014). Relationships among rheas, tinamous and the
# Australasian ratites (kiwi, emu, cassowaries) are still unresolved.
#
# Alternative tree (MCC topology otherwise unchanged):
#   (Struthio, (Rheidae, Tinamidae, Australasian ratites))
# with the three-way node left unresolved (a polytomy). Ages: the palaeognath
# root keeps its MCC age (65.8 Myr); the new polytomy is placed halfway
# between the root and the tinamou crown (both MCC ages), because it must be
# older than the tinamou crown and no single published age fits this
# unresolved node. Crown ages of all subclades are unchanged.
# Kiwi + emu + cassowaries stay one clade, as on BirdTree (both published
# alternatives keep them together or unresolved), so ratites become THREE
# independent losses: ostrich, rheas, Australasian ratites.
#
# Re-run on this tree: PGLS M1 + one-sided equivalence bound (R/12); paired
# sister contrasts with an exact sign-flip test; multi-rate BM (R2, R/09).
# Note: all three ratite origins have the tinamous as their volant sister
# group, so their paired contrasts are not independent of each other.
# Output: output/tables/ratite_topology_pgls.csv, ratite_topology_pairs.csv,
#         ratite_topology_rates.csv. Deterministic.

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  library(phytools)
  library(dplyr)
})

mcc <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
DATA <- list(
  montgomerie = list(d = readRDS(file.path(PATHS$processed, "analysis_montgomerie.rds")),
                     traits = c("elongation", "pointedness", "polar_asym")),
  stoddard    = list(d = readRDS(file.path(PATHS$processed, "analysis_stoddard.rds")),
                     traits = c("ellipticity", "asymmetry"))
)
RATITE_RE <- "^(Struthio|Rhea|Dromaius|Casuarius|Apteryx)_"
TINAMOU_FAMILY <- "Tinamidae"

node_age <- function(tree, node) max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[node]
clade_newick <- function(tree, tips, stem) {
  # Newick of the clade of `tips` (or a single tip) with a stem of length `stem`
  if (length(tips) == 1) return(paste0(tips, ":", stem))
  sub <- extract.clade(tree, getMRCA(tree, tips))
  paste0(sub(";$", "", write.tree(sub)), ":", stem)
}

#' Rebuild the palaeognath clade of `tree` as (Struthio,(Rheas,Tinamous,Australasians)).
alt_ratite_tree <- function(tree, d) {
  tree$node.label <- NULL
  tips <- tree$tip.label
  fam <- setNames(d$family, d$birdtree_name)
  tin <- tips[fam[tips] == TINAMOU_FAMILY]
  ost <- grep("^Struthio_", tips, value = TRUE)
  rhea <- grep("^Rhea_", tips, value = TRUE)
  aus <- grep("^(Dromaius|Casuarius|Apteryx)_", tips, value = TRUE)
  stopifnot(length(ost) == 1, length(tin) > 0, length(rhea) > 0, length(aus) > 0)
  pal <- c(ost, rhea, aus, tin)
  stopifnot(is.monophyletic(tree, pal))
  root_age <- node_age(tree, getMRCA(tree, pal))
  crown <- function(x) if (length(x) == 1) 0 else node_age(tree, getMRCA(tree, x))
  t_poly <- (root_age + crown(tin)) / 2
  stopifnot(t_poly > max(crown(tin), crown(rhea), crown(aus)))
  new_pal <- read.tree(text = paste0(
    "(", ost, ":", root_age, ",(",
    clade_newick(tree, rhea, t_poly - crown(rhea)), ",",
    clade_newick(tree, tin, t_poly - crown(tin)), ",",
    clade_newick(tree, aus, t_poly - crown(aus)), "):", root_age - t_poly, ");"))
  # Replace the palaeognath clade: keep one placeholder tip, cut its branch at
  # the old crown node, and graft the new clade there.
  x <- drop.tip(tree, setdiff(pal, ost))
  e <- which(x$edge[, 2] == which(x$tip.label == ost))
  x$edge.length[e] <- x$edge.length[e] - root_age
  x$tip.label[x$tip.label == ost] <- "PLACEHOLDER"
  out <- bind.tree(x, new_pal, where = which(x$tip.label == "PLACEHOLDER"))
  stopifnot(setequal(out$tip.label, tips), is.ultrametric(out, tol = 1e-6))
  out
}

fit_m1 <- function(d, tree, trait) {
  d <- as.data.frame(d); rownames(d) <- d$birdtree_name
  d <- d[tree$tip.label, ]
  fit <- phylolm(as.formula(paste(trait, "~ flight_status + log_mass + log_egg_size")),
                 data = d, phy = tree, model = "lambda")
  s <- summary(fit)$coefficients
  sd_vol <- sd(d[[trait]][d$flight_status == "volant"])
  tm <- "flight_statusterrestrial"
  dfree <- nrow(d) - nrow(s)
  tibble(std_estimate = s[tm, "Estimate"] / sd_vol, p = s[tm, "p.value"],
         delta_min_sd = (qt(0.95, dfree) * s[tm, "StdErr"] - s[tm, "Estimate"]) / sd_vol)
}

signflip_p <- function(x) {
  # exact one-sided sign-flip test of mean(x) < 0
  s <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(x))))
  mean(as.vector(s %*% abs(x)) / length(x) <= mean(x) + 1e-12)
}

pgls_rows <- list(); pair_rows <- list(); rate_rows <- list()
for (dn in names(DATA)) {
  d <- DATA[[dn]]$d
  base <- keep.tip(mcc, d$birdtree_name); base$node.label <- NULL
  alt <- alt_ratite_tree(base, d)
  fl_tips <- d$birdtree_name[d$flightless]
  origins <- flightless_clades(alt, fl_tips)
  log_msg(dn, ": ", length(flightless_clades(base, fl_tips)), " origins on MCC, ",
          length(origins), " on the alternative ratite tree")
  for (tr in DATA[[dn]]$traits) {
    pgls_rows[[length(pgls_rows) + 1]] <- bind_rows(
      fit_m1(d, base, tr) |> mutate(tree = "MCC (ratites one origin)"),
      fit_m1(d, alt, tr) |> mutate(tree = "alternative (ratites three origins)")) |>
      mutate(dataset = dn, trait = tr, .before = 1)
  }
  # Paired contrasts: each terrestrial origin vs volant descendants of its parent node
  sd_vol <- sapply(DATA[[dn]]$traits, function(tr) sd(d[[tr]][d$flight_status == "volant"]))
  for (cl in origins) {
    loco <- unique(as.character(d$flight_status[match(cl, d$birdtree_name)]))
    node <- if (length(cl) == 1) which(alt$tip.label == cl) else getMRCA(alt, cl)
    parent <- alt$edge[alt$edge[, 2] == node, 1]
    sis <- setdiff(extract.clade(alt, parent)$tip.label, c(cl, fl_tips))
    if (length(sis) == 0) next
    for (tr in DATA[[dn]]$traits) {
      fl_mean <- mean(d[[tr]][match(cl, d$birdtree_name)])
      sis_mean <- mean(d[[tr]][match(sis, d$birdtree_name)])
      pair_rows[[length(pair_rows) + 1]] <- tibble(
        dataset = dn, trait = tr, locomotion = loco,
        origin = if (length(cl) > 3) paste0(cl[1], " + ", length(cl) - 1, " spp.") else paste(cl, collapse = ", "),
        n_sister = length(sis), sister_is_tinamous = all(grepl("Tinam|Crypturellus|Eudromia|Nothura|Rhynchotus|Nothocercus|Nothoprocta|Tinamotis", sis)),
        diff_std = (fl_mean - sis_mean) / sd_vol[[tr]])
    }
  }
  # Multi-rate BM, R2 (volant vs flightless), as in R/09
  st_alt <- paintSubTree(alt, node = Ntip(alt) + 1, state = "volant")
  for (cl in origins) {
    node <- if (length(cl) == 1) which(alt$tip.label == cl) else getMRCA(alt, cl)
    st_alt <- paintSubTree(st_alt, node = node, state = "flightless", stem = 0.5)
  }
  for (tr in DATA[[dn]]$traits) {
    x <- setNames(as.numeric(scale(d[[tr]])), d$birdtree_name)[st_alt$tip.label]
    f <- brownie.lite(st_alt, x)
    rate_rows[[length(rate_rows) + 1]] <- tibble(
      dataset = dn, trait = tr,
      ratio_flightless = f$sig2.multiple["flightless"] / f$sig2.multiple["volant"],
      lrt_p = f$P.chisq, converged = f$convergence == "Optimization has converged.")
  }
}

pgls <- bind_rows(pgls_rows) |> mutate(across(where(is.double), ~ signif(.x, 4)))
pairs <- bind_rows(pair_rows)
# Test set (Montgomerie, terrestrial): the three ratite contrasts from the
# alternative tree, plus the non-ratite terrestrial pairs exactly as used in
# R/08 (same partners, including the literature partner of the Tasmanian
# nativehen; the Inaccessible Island rail stays excluded, D2.x).
orig <- read.csv(file.path(PATHS$tables, "paired_differences.csv")) |>
  filter(locomotion == "terrestrial", origin_id != "Ratites") |>
  transmute(dataset = "montgomerie", trait, locomotion, origin = origin_id,
            n_sister, sister_is_tinamous = FALSE, diff_std, source = "R/08 pair")
alt_rat <- pairs |>
  filter(dataset == "montgomerie", locomotion == "terrestrial",
         grepl(RATITE_RE, sub(" .*$", "", origin))) |>
  mutate(source = "alternative tree")
test_set <- bind_rows(alt_rat, orig)
tests <- test_set |>
  group_by(dataset, trait) |>
  summarise(n_pairs = n(), mean_diff_std = mean(diff_std),
            n_negative = sum(diff_std < 0),
            signflip_p_one_sided = signflip_p(diff_std), .groups = "drop")
pairs <- test_set
rates <- bind_rows(rate_rows)
write_table(pgls, "ratite_topology_pgls")
write_table(bind_rows(pairs |> mutate(across(where(is.double), ~ signif(.x, 4))),
                      tests |> mutate(origin = "TEST (terrestrial pairs)",
                                      across(where(is.double), ~ signif(.x, 4)))),
            "ratite_topology_pairs")
write_table(rates |> mutate(across(where(is.double), ~ signif(.x, 4))), "ratite_topology_rates")
print(as.data.frame(pgls)); print(as.data.frame(pairs)); print(as.data.frame(tests)); print(as.data.frame(rates))
