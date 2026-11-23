# 03_flightless.R
# Build the flightless-species list, classify locomotion (terrestrial vs
# diving), intersect with the egg-shape data and count independent origins.
#
# Inputs
#   Flight ability : Sayol et al. 2020 Data S2 ("Volancy"), via 02_clean.R
#   Locomotion     : Saitta et al. 2025 Dataset S1, swimming-mode columns
#                    (each taxon cites its own natural-history source)
#   Overrides      : data/curated/flightless_overrides.csv (manual classifications
#                    of flagged species; may be empty)
#
# Locomotion rule
#   diving      = Saitta mode 3 (underwater hind-feet rowing), 4 or 5
#                 (underwater wing-propelled flight) recorded
#   terrestrial = no aquatic mode recorded (modes 1-5 all empty)
#   review      = surface swimming only (modes 1-2), or species not scored by
#                 Saitta and its family has no unanimous scored members.
#   Unscored species inherit the family label only if every Saitta-scored
#   flightless member of that family has the same label.

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ape)
})

sayol  <- readRDS(file.path(PATHS$processed, "sayol_birdtree.rds"))
egg    <- readRDS(file.path(PATHS$processed, "egg_shape_birdtree.rds"))
mont   <- readRDS(file.path(PATHS$processed, "egg_montgomerie_birdtree.rds"))

# ---- Saitta swimming modes ---------------------------------------------------
saitta <- read.csv(file.path(PATHS$raw, "saitta2025", "extracted", "Dataset.csv"),
                   check.names = FALSE, fileEncoding = "UTF-8-BOM")
mode_cols <- names(saitta)[7:11]          # modes 1-5, see column headers
has_mode <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))
saitta_loco <- saitta |>
  transmute(
    saitta_taxon = Taxon,
    # binomial = first two words after removing "(=...)" synonyms
    binomial = to_birdtree_name(sub("^(\\S+)\\s+(\\S+).*$", "\\1 \\2",
                                     gsub("\\(=[^)]*\\)\\s*", "", Taxon))),
    surface    = has_mode(.data[[mode_cols[1]]]) | has_mode(.data[[mode_cols[2]]]),
    underwater = has_mode(.data[[mode_cols[3]]]) | has_mode(.data[[mode_cols[4]]]) |
                 has_mode(.data[[mode_cols[5]]]),
    saitta_source = .data[[names(saitta)[13]]],
    saitta_volant = trimws(.data[[names(saitta)[2]]])  # Y / N / N* / P / ? ...
  ) |>
  distinct(binomial, .keep_all = TRUE) |>
  mutate(saitta_class = case_when(
    underwater ~ "diving",
    surface    ~ "review",
    TRUE       ~ "terrestrial"
  ))

# Saitta uses some older names; map them to BirdTree via the Sayol table
# (Sayol original name -> BirdTree) where needed.
name_map <- sayol |> transmute(q = to_birdtree_name(species), birdtree_name)
saitta_fixes <- read_csv(file.path(PATHS$curated, "name_fixes.csv"),
                         show_col_types = FALSE, comment = "#") |>
  filter(dataset == "saitta2025") |>
  transmute(q = to_birdtree_name(original_name), to = to_birdtree_name(birdtree_name))
saitta_loco <- saitta_loco |>
  mutate(binomial = coalesce(saitta_fixes$to[match(binomial, saitta_fixes$q)], binomial),
         birdtree_name = coalesce(name_map$birdtree_name[match(binomial, name_map$q)],
                                  binomial))

# ---- Flightless list ---------------------------------------------------------
fl <- sayol |>
  filter(Volancy != "Volant") |>
  group_by(birdtree_name) |>
  summarise(sayol_names = paste(species, collapse = "; "),
            sayol_volancy = paste(unique(Volancy), collapse = "; "),
            order = first(Order), family = first(Family), .groups = "drop") |>
  left_join(saitta_loco |> select(birdtree_name, saitta_taxon, saitta_class,
                                  saitta_source, saitta_volant),
            by = "birdtree_name")

family_rule <- fl |>
  filter(!is.na(saitta_class)) |>
  group_by(family) |>
  summarise(family_class = if (n_distinct(saitta_class) == 1) first(saitta_class) else NA_character_,
            .groups = "drop")

fl <- fl |>
  left_join(family_rule, by = "family") |>
  mutate(
    locomotion = coalesce(saitta_class, family_class, "review"),
    locomotion_basis = case_when(
      !is.na(saitta_class) ~ paste0("Saitta et al. 2025 swim modes (", saitta_taxon, ")"),
      !is.na(family_class) ~ "Inherited: all Saitta-scored flightless members of family agree",
      TRUE ~ "Not scored by Saitta; family not unanimous"
    ),
    # Flight status itself is uncertain if Sayol is not "Flightless" or
    # Saitta records the species as "?" (unknown) or "P" (poor flier)
    volancy_disputed = sayol_volancy != "Flightless" |
      (!is.na(saitta_volant) & saitta_volant %in% c("?", "P")),
    locomotion = if_else(volancy_disputed, "review", locomotion),
    review_reason = case_when(
      sayol_volancy != "Flightless" ~ paste0("Sayol volancy = '", sayol_volancy, "'"),
      volancy_disputed ~ paste0("Sayol = Flightless but Saitta volant = '", saitta_volant, "'"),
      locomotion == "review" & !is.na(saitta_class) ~ "Surface swimmer only (no underwater mode)",
      locomotion == "review" ~ "No locomotion data",
      TRUE ~ ""
    )
  )

# Apply manual classifications of flagged species, if any
ov_path <- file.path(PATHS$curated, "flightless_overrides.csv")
if (file.exists(ov_path)) {
  ov <- read_csv(ov_path, show_col_types = FALSE, comment = "#")
  idx <- match(ov$birdtree_name, fl$birdtree_name)
  fl$locomotion[idx] <- ov$locomotion
  fl$locomotion_basis[idx] <- paste0("Manual: ", ov$reason)
  fl$review_reason[idx] <- ""
}

fl$in_stoddard    <- fl$birdtree_name %in% egg$birdtree_name
fl$in_montgomerie <- fl$birdtree_name %in% mont$birdtree_name
fl$in_egg_data    <- fl$in_stoddard | fl$in_montgomerie

write_csv(fl |> select(birdtree_name, order, family, sayol_names, sayol_volancy,
                       locomotion, locomotion_basis, review_reason, in_stoddard,
                       in_montgomerie,
                       saitta_source),
          file.path(PATHS$curated, "flightless_species.csv"), na = "")
saveRDS(fl, file.path(PATHS$processed, "flightless.rds"))

# ---- Independent origins -----------------------------------------------------
# Counted on BirdTree tree 1 (palaeognath topology is identical across the
# posterior: ratites monophyletic in 100/100 sampled trees).
# An origin = a maximal clade of flightless species among species with egg
# data (Fitch parsimony on the pruned tree, flightless = derived state).
tre_file <- file.path(PATHS$raw, "birdtree", "HackettStage2_0001_1000", "AllBirdsHackett1.tre")
tree <- read.tree(text = readLines(tre_file, n = 1))

# flightless_clades() is defined in 00_setup.R
count_origins <- function(tree, tips) length(flightless_clades(tree, tips))

# Count separately for each egg dataset (each on its own pruned tree) and for
# their union, by locomotion class.
datasets <- list(
  Stoddard    = egg$birdtree_name,
  Montgomerie = mont$birdtree_name,
  Union       = union(egg$birdtree_name, mont$birdtree_name)
)
summary_tbl <- bind_rows(lapply(names(datasets), function(d) {
  spp <- datasets[[d]]
  tr <- keep.tip(tree, spp)
  fl_d <- fl |> filter(birdtree_name %in% spp)
  bind_rows(lapply(c("terrestrial", "diving", "review"), function(cl) {
    tips <- fl_d$birdtree_name[fl_d$locomotion == cl]
    tibble(dataset = d, n_egg_species = length(spp), locomotion = cl,
           n_species = length(tips), origins_birdtree = count_origins(tr, tips),
           species = paste(sort(tips), collapse = ", "))
  }))
}))
print(summary_tbl |> select(-species), n = Inf)
write_table(summary_tbl, "stage1_flightless_in_egg_data")
write_table(fl |> filter(locomotion == "review") |>
              select(birdtree_name, family, sayol_volancy, review_reason,
                     in_stoddard, in_montgomerie),
            "flightless_needs_review")
