# 00_setup.R
# Shared paths, constants and helpers. Sourced at the top of every numbered
# script so that each script can also be run on its own.

suppressPackageStartupMessages({
  library(here)
})

# ---- Paths -----------------------------------------------------------------
PATHS <- list(
  raw       = here::here("data", "raw"),
  curated   = here::here("data", "curated"),   # hand-curated inputs (tracked in git)
  processed = here::here("data", "processed"),
  figures   = here::here("output", "figures"),
  tables    = here::here("output", "tables"),
  report    = here::here("report")
)
invisible(lapply(PATHS, dir.create, showWarnings = FALSE, recursive = TRUE))

# ---- Reproducibility --------------------------------------------------------
# One global seed; scripts that draw random numbers call set.seed(SEED) again
# locally so that their output does not depend on what ran before them.
SEED <- 1213L
set.seed(SEED)

# Number of posterior trees sampled from the BirdTree 1,000-tree set
N_POSTERIOR_TREES <- 100L

# ---- Helpers ----------------------------------------------------------------

#' Normalise a binomial to BirdTree style ("Genus_species").
to_birdtree_name <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  # Capitalise genus, lower-case epithet
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}

#' Write a CSV table to output/tables with a consistent format.
write_table <- function(df, name) {
  path <- file.path(PATHS$tables, paste0(name, ".csv"))
  readr::write_csv(df, path, na = "")
  message("Wrote ", path)
  invisible(path)
}

#' Split flightless tips into independent origins.
#' An origin is a maximal clade whose tips are all flightless (i.e. the tips
#' that a single loss of flight on the stem branch would explain).
#' Returns a list of character vectors of tip labels, one per origin.
flightless_clades <- function(tree, flightless_tips) {
  flightless_tips <- intersect(flightless_tips, tree$tip.label)
  clades <- list()
  remaining <- flightless_tips
  while (length(remaining) > 0) {
    node <- which(tree$tip.label == remaining[1])
    repeat {
      parent <- tree$edge[tree$edge[, 2] == node, 1]
      if (length(parent) == 0) break
      if (all(ape::extract.clade(tree, parent)$tip.label %in% flightless_tips)) {
        node <- parent
      } else break
    }
    clade <- if (node <= ape::Ntip(tree)) tree$tip.label[node] else
      ape::extract.clade(tree, node)$tip.label
    clades[[length(clades) + 1]] <- clade
    remaining <- setdiff(remaining, clade)
  }
  clades
}

#' Log a line to the console with a timestamp.
log_msg <- function(...) {
  message(format(Sys.time(), "[%H:%M:%S] "), ...)
}
