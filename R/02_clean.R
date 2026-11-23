# 02_clean.R
# Read the egg-shape (Stoddard et al. 2017), flight-ability (Sayol et al. 2020)
# and AVONET data, harmonise all species names to BirdTree (Jetz et al. 2012)
# taxonomy, and report match rates.
#
# Name-matching cascade (first hit wins):
#   1. exact match to a BirdTree tip label
#   2. AVONET BirdLife -> BirdTree crosswalk (only one-to-one targets)
#   3. AVONET eBird -> BirdLife -> BirdTree crosswalk (only one-to-one targets)
#   4. hand-curated fixes in data/curated/name_fixes.csv (each with a reason)
# Anything still unmatched is written to output/tables/unmatched_*.csv.

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(readr)
})

avonet_path <- file.path(PATHS$raw, "avonet", "AVONET_Supplementary_dataset_1.xlsx")
tax <- read_csv(file.path(PATHS$raw, "birdtree", "BLIOCPhyloMasterTax.csv"),
                show_col_types = FALSE)
bt_tips <- tax$TipLabel

# ---- Crosswalk tables -------------------------------------------------------
bl_bt <- read_excel(avonet_path, sheet = "BirdLife–BirdTree crosswalk", na = c("", "NA")) |>
  filter(!is.na(Species3)) |>
  mutate(from = to_birdtree_name(Species1), to = to_birdtree_name(Species3)) |>
  distinct(from, to)
# keep only names that map to a single BirdTree species
bl_bt <- bl_bt |> group_by(from) |> filter(n() == 1) |> ungroup()

eb_bl <- read_excel(avonet_path, sheet = "BirdLife-eBird crosswalk", na = c("", "NA")) |>
  filter(!is.na(Species2), !is.na(Species1)) |>
  mutate(from = to_birdtree_name(Species2), via = to_birdtree_name(Species1)) |>
  distinct(from, via) |>
  inner_join(bl_bt, by = c("via" = "from")) |>
  distinct(from, to) |>
  group_by(from) |> filter(n() == 1) |> ungroup()

fixes_path <- file.path(PATHS$curated, "name_fixes.csv")
fixes <- if (file.exists(fixes_path)) {
  read_csv(fixes_path, show_col_types = FALSE, comment = "#") |>
    mutate(from = to_birdtree_name(original_name), to = to_birdtree_name(birdtree_name))
} else {
  tibble(dataset = character(), from = character(), to = character(), reason = character())
}

#' Map a vector of names to BirdTree; returns a tibble with the method used.
match_birdtree <- function(names, dataset) {
  nm <- to_birdtree_name(names)
  fx <- fixes |> filter(.data$dataset == !!dataset)
  out <- tibble(original_name = names, query = nm) |>
    mutate(
      birdtree_name = case_when(
        query %in% bt_tips ~ query,
        query %in% bl_bt$from ~ bl_bt$to[match(query, bl_bt$from)],
        query %in% eb_bl$from ~ eb_bl$to[match(query, eb_bl$from)],
        query %in% fx$from ~ fx$to[match(query, fx$from)],
        TRUE ~ NA_character_
      ),
      match_method = case_when(
        query %in% bt_tips ~ "exact",
        query %in% bl_bt$from ~ "BirdLife crosswalk",
        query %in% eb_bl$from ~ "eBird crosswalk",
        query %in% fx$from ~ "manual fix",
        TRUE ~ "unmatched"
      )
    )
  # Guard: every mapped name must be a real tip
  stopifnot(all(is.na(out$birdtree_name) | out$birdtree_name %in% bt_tips))
  out
}

# ---- Egg shape (Stoddard et al. 2017) --------------------------------------
egg_raw <- read_excel(
  file.path(PATHS$raw, "stoddard2017", "aaj1945_datas1_egg_shape_by_species_v2.xlsx"),
  sheet = 1
) |>
  rename(asymmetry = Asymmetry, ellipticity = Ellipticity,
         egg_length_cm = `AvgLength (cm)`, n_images = `Number of images`,
         n_eggs = `Number of eggs`)

log_msg("Stoddard rows: ", nrow(egg_raw),
        "; rows with missing shape: ", sum(is.na(egg_raw$asymmetry)))

egg_raw <- egg_raw |> filter(!is.na(asymmetry), !is.na(ellipticity))
egg_match <- match_birdtree(egg_raw$Species, "stoddard2017")
egg_raw <- bind_cols(egg_raw, egg_match |> select(birdtree_name, match_method))

# Several source names can collapse onto one BirdTree species (duplicates in
# the source, or lumping). Combine them with an egg-count-weighted mean.
egg <- egg_raw |>
  filter(!is.na(birdtree_name)) |>
  group_by(birdtree_name) |>
  summarise(
    n_source_rows = n(),
    source_names = paste(unique(Species), collapse = "; "),
    stoddard_order = first(Order),
    asymmetry = weighted.mean(asymmetry, n_eggs),
    ellipticity = weighted.mean(ellipticity, n_eggs),
    egg_length_cm = weighted.mean(egg_length_cm, n_eggs),
    n_eggs = sum(n_eggs),
    .groups = "drop"
  )

# ---- Egg shape (Montgomerie et al.; Biggins et al. 2018 metrics) -----------
# The authors supply their own BirdTree name (scinameBT); we still run it
# through the cascade so that the same rules apply to every dataset.
mont_raw <- read_csv(file.path(PATHS$raw, "montgomerie2022", "DATAeggshape.csv"),
                     show_col_types = FALSE)
mont_match <- match_birdtree(mont_raw$scinameBT, "montgomerie2022")
mont_raw <- bind_cols(mont_raw, mont_match |> select(birdtree_name, match_method))

# Species-level file has no per-species egg counts, so rows collapsing onto one
# BirdTree species are combined with an unweighted mean.
mont <- mont_raw |>
  filter(!is.na(birdtree_name)) |>
  group_by(birdtree_name) |>
  summarise(
    n_source_rows = n(),
    source_names = paste(unique(sciname), collapse = "; "),
    elongation = mean(avEL), pointedness = mean(avPT), polar_asym = mean(avPA),
    egg_length_mm = mean(avLENGTH), egg_width_mm = mean(avWIDTH),
    egg_volume_ml = mean(avVOL),
    female_mass_g = mean(ADmass), clutch_size = mean(CSegg),
    inc_site = first(IncSite), loco_montgomerie = first(LOCO),
    devel = first(devel), hwi_montgomerie = mean(HWI),
    .groups = "drop"
  )

# ---- Flight ability (Sayol et al. 2020) ------------------------------------
sayol_raw <- read_excel(file.path(PATHS$raw, "sayol2020", "abb6095_data_file_s2.xlsx"))
sayol_match <- match_birdtree(sayol_raw$species, "sayol2020")
sayol <- bind_cols(sayol_raw, sayol_match |> select(birdtree_name, match_method))

# ---- AVONET (already in BirdTree taxonomy) ---------------------------------
avonet <- read_excel(avonet_path, sheet = "AVONET3_BirdTree") |>
  transmute(
    birdtree_name = to_birdtree_name(Species3),
    family = Family3, order = Order3,
    mass_g = Mass, hwi = `Hand-Wing.Index`,
    habitat = Habitat, migration = Migration,
    trophic_niche = Trophic.Niche, primary_lifestyle = Primary.Lifestyle
  )
stopifnot(all(avonet$birdtree_name %in% bt_tips))

# ---- Reports ----------------------------------------------------------------
match_summary <- bind_rows(
  egg_match |> mutate(dataset = "Stoddard 2017 egg shape"),
  mont_match |> mutate(dataset = "Montgomerie egg shape"),
  sayol_match |> mutate(dataset = "Sayol 2020 flight ability")
) |>
  count(dataset, match_method) |>
  group_by(dataset) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  ungroup()
print(match_summary, n = Inf)
write_table(match_summary, "stage1_name_match_summary")

write_table(
  egg_raw |> filter(is.na(birdtree_name)) |>
    select(Order, Family, Species, n_eggs),
  "unmatched_stoddard2017"
)
write_table(
  sayol |> filter(is.na(birdtree_name)) |>
    select(Order, Family, species, Volancy),
  "unmatched_sayol2020"
)
write_table(
  egg |> filter(n_source_rows > 1) |> select(birdtree_name, source_names, n_source_rows),
  "stoddard_merged_rows"
)
write_table(
  mont_raw |> filter(is.na(birdtree_name)) |>
    select(Order, Family, sciname, scinameBT, common, LOCO),
  "unmatched_montgomerie2022"
)
write_table(
  mont |> filter(n_source_rows > 1) |> select(birdtree_name, source_names, n_source_rows),
  "montgomerie_merged_rows"
)

# ---- Save -------------------------------------------------------------------
saveRDS(egg,    file.path(PATHS$processed, "egg_shape_birdtree.rds"))
saveRDS(mont,   file.path(PATHS$processed, "egg_montgomerie_birdtree.rds"))
saveRDS(sayol,  file.path(PATHS$processed, "sayol_birdtree.rds"))
saveRDS(avonet, file.path(PATHS$processed, "avonet_birdtree.rds"))
log_msg("Egg-shape species after harmonisation: ", nrow(egg))
