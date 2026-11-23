# 13_origin_inventory.R
# How many independent losses of flight exist among extant birds, and how many
# lack egg-shape data? Sets the ceiling for the power analysis in
# R/14_power_simulation.R.
#
# Uses the first tree of the full (9,993-species) BirdTree Hackett Stage 2
# set, i.e. without pruning to species with egg data, and the 57 extant
# flightless species of Sayol et al. (2020) as matched in R/03. Origins are
# maximal all-flightless clades (same definition as R/05). Many of the extra
# rails are placed in BirdTree by taxonomy alone, so the count is approximate.
# Output: output/tables/origin_inventory.csv

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(ape))

tre_file <- file.path(PATHS$raw, "birdtree", "HackettStage2_0001_1000", "AllBirdsHackett1.tre")
tree <- read.tree(text = readLines(tre_file, n = 1))  # first of the 1,000 trees
fl <- readRDS(file.path(PATHS$processed, "flightless.rds"))
stopifnot(all(fl$birdtree_name %in% tree$tip.label))

clades <- flightless_clades(tree, fl$birdtree_name)
loco <- setNames(fl$locomotion, fl$birdtree_name)
has_egg <- setNames(fl$in_egg_data, fl$birdtree_name)
inv <- do.call(rbind, lapply(clades, function(cl) data.frame(
  n_species = length(cl),
  locomotion = paste(sort(unique(loco[cl])), collapse = "/"),
  has_egg_data = any(has_egg[cl]),
  species = paste(sort(cl), collapse = ", ")
)))
inv <- inv[order(!inv$has_egg_data, inv$locomotion, inv$species), ]
write_table(inv, "origin_inventory")

log_msg(nrow(inv), " origins on the full tree; ", sum(inv$has_egg_data), " with egg data")
print(table(locomotion = inv$locomotion, egg_data = inv$has_egg_data))
