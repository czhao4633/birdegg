# 06_analysis_data.R
# Assemble one analysis table per egg dataset: egg shape + flight status +
# AVONET covariates (+ clutch size and incubation site for Montgomerie).
#
# Flight status (factor, reference = volant):
#   volant, terrestrial (flightless), diving (flightless)
# Covariates:
#   log_mass     = log(AVONET body mass, g)            -- both datasets
#   log_egg_size = log(egg length, mm)                 -- both datasets
#                  (Stoddard gives no egg width, so length is the only size
#                   measure available in both; Montgomerie volume is used in
#                   a sensitivity model)
#   hwi          = AVONET hand-wing index

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(dplyr))

avonet <- readRDS(file.path(PATHS$processed, "avonet_birdtree.rds"))
fl     <- readRDS(file.path(PATHS$processed, "flightless.rds")) |>
  filter(locomotion %in% c("terrestrial", "diving"))

add_common <- function(d) {
  d |>
    left_join(avonet, by = "birdtree_name") |>
    mutate(
      flight_status = factor(
        coalesce(fl$locomotion[match(birdtree_name, fl$birdtree_name)], "volant"),
        levels = c("volant", "terrestrial", "diving")),
      flightless = flight_status != "volant",
      log_mass = log(mass_g)
    )
}

mont <- readRDS(file.path(PATHS$processed, "egg_montgomerie_birdtree.rds")) |>
  add_common() |>
  mutate(log_egg_size = log(egg_length_mm),
         log_egg_volume = log(egg_volume_ml),
         log_clutch = log(clutch_size))

stod <- readRDS(file.path(PATHS$processed, "egg_shape_birdtree.rds")) |>
  add_common() |>
  mutate(log_egg_size = log(egg_length_cm * 10))

report <- function(d, name, traits) {
  vars <- c(traits, "log_mass", "log_egg_size", "hwi")
  log_msg(name, ": ", nrow(d), " spp.; complete cases for core models: ",
          sum(complete.cases(d[, vars])))
  print(table(d$flight_status))
  miss <- colSums(is.na(d[, vars]))
  if (any(miss > 0)) print(miss[miss > 0])
}
report(mont, "Montgomerie", c("elongation", "pointedness", "polar_asym"))
report(stod, "Stoddard", c("ellipticity", "asymmetry"))
log_msg("Montgomerie missing clutch size: ", sum(is.na(mont$clutch_size)),
        "; missing incubation site: ", sum(is.na(mont$inc_site)))

saveRDS(mont, file.path(PATHS$processed, "analysis_montgomerie.rds"))
saveRDS(stod, file.path(PATHS$processed, "analysis_stoddard.rds"))
