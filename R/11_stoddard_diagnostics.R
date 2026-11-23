# 11_stoddard_diagnostics.R
# Reproduces the diagnosis of the one "significant" flight-status effect
# (Stoddard ellipticity, terrestrial):
#   1. M1 without body-mass / egg-size covariates
#   2. M1 with ratites removed (terrestrial group = kagu only)
#   3. Body-mass ranges of ratites vs the heaviest volant species
#   4. Raw mean ellipticity of ratites vs their sister group (tinamous)
# Output: output/tables/stoddard_diagnostics.csv

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  library(dplyr)
})

mcc  <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
d    <- readRDS(file.path(PATHS$processed, "analysis_stoddard.rds"))
pairs <- readRDS(file.path(PATHS$processed, "sister_pairs.rds"))
is_ratite <- grepl("^(Struthio|Rhea|Dromaius|Casuarius|Apteryx)_", d$birdtree_name)

fit_term <- function(dd, rhs, term = "flight_statusterrestrial") {
  dd <- as.data.frame(dd)
  rownames(dd) <- dd$birdtree_name
  tr <- keep.tip(mcc, dd$birdtree_name)
  dd <- dd[tr$tip.label, ]
  s <- summary(phylolm(as.formula(paste("ellipticity", rhs)), data = dd, phy = tr,
                       model = "lambda"))$coefficients
  sd_vol <- sd(dd$ellipticity[dd$flight_status == "volant"])
  tibble(estimate = s[term, "Estimate"], std_estimate = s[term, "Estimate"] / sd_vol,
         p = s[term, "p.value"], n = nrow(dd),
         n_terrestrial = sum(dd$flight_status == "terrestrial"))
}

rat_pair <- pairs |> filter(dataset == "stoddard", origin_id == "Ratites")
sister <- strsplit(rat_pair$sister_taxa, ", ")[[1]]

out <- bind_rows(
  fit_term(d, "~ flight_status + log_mass + log_egg_size") |>
    mutate(check = "M1 as fitted (MCC tree)"),
  fit_term(d, "~ flight_status") |>
    mutate(check = "M1 without mass and egg-size covariates"),
  fit_term(d[!is_ratite, ], "~ flight_status + log_mass + log_egg_size") |>
    mutate(check = "M1 without ratites (terrestrial = kagu only)")
) |> relocate(check)

mass <- tibble(
  check = c("Heaviest volant species, mass (g)", "Lightest ratite, mass (g)",
            "Ratites heavier than heaviest volant species",
            "Mean ellipticity, ratites", "Mean ellipticity, sister group (tinamous)"),
  value = c(max(d$mass_g[d$flight_status == "volant"]),
            min(d$mass_g[is_ratite & !grepl("^Apteryx", d$birdtree_name)]),
            sum(d$mass_g[is_ratite] > max(d$mass_g[d$flight_status == "volant"])),
            mean(d$ellipticity[is_ratite]),
            mean(d$ellipticity[match(sister, d$birdtree_name)]))
)
print(out, width = Inf); print(mass, width = Inf)
log_msg("Ratites in Stoddard: ", sum(is_ratite), " (",
        paste(sort(d$birdtree_name[is_ratite]), collapse = ", "), ")")
log_msg("Ratite sister group: ", paste(sister, collapse = ", "))
write_table(bind_rows(out |> mutate(across(where(is.double), ~ signif(.x, 4))),
                      mass |> mutate(value = signif(value, 4))),
            "stoddard_diagnostics")
