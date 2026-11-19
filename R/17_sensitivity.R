# 17_sensitivity.R
# Does the PGLS result (no terrestrial flight-status effect, M1) depend on a
# single origin, on two judgement calls in the flightless classification, or
# on how body mass is modelled?
#
#   A. Leave-one-origin-out: drop each flightless origin in turn.
#   B. Classification: kagu treated as volant (it is listed as a weak flyer);
#      Tasmanian nativehen removed (terrestrial vs surface-swimmer call).
#   C. Body mass: (i) natural-spline terms (3 df) for log mass and log egg
#      size instead of straight lines; (ii) only flightless species inside
#      the volant mass range (removes extrapolation entirely).
# MCC tree, Pagel's lambda, both datasets. For each fit: terrestrial and
# diving effects in volant-SD units, p, and the one-sided equivalence bound
# delta_min of R/12. Deterministic (no random numbers).
# Output: output/tables/sensitivity_pgls.csv

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  library(dplyr)
  library(splines)
})

mcc <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))
pairs <- readRDS(file.path(PATHS$processed, "sister_pairs.rds"))
DATA <- list(
  montgomerie = list(d = readRDS(file.path(PATHS$processed, "analysis_montgomerie.rds")),
                     traits = c("elongation", "pointedness", "polar_asym")),
  stoddard    = list(d = readRDS(file.path(PATHS$processed, "analysis_stoddard.rds")),
                     traits = c("ellipticity", "asymmetry"))
)
M1 <- "~ flight_status + log_mass + log_egg_size"
SPLINE <- "~ flight_status + ns(log_mass, df = 3) + ns(log_egg_size, df = 3)"
TERMS <- c("flight_statusterrestrial", "flight_statusdiving")

fit_terms <- function(d, trait, rhs) {
  d <- as.data.frame(d)
  d$flight_status <- droplevels(factor(d$flight_status, levels = c("volant", "terrestrial", "diving")))
  rownames(d) <- d$birdtree_name
  tr <- keep.tip(mcc, d$birdtree_name)
  d <- d[tr$tip.label, , drop = FALSE]
  fit <- phylolm(as.formula(paste(trait, rhs)), data = d, phy = tr, model = "lambda")
  s <- summary(fit)$coefficients
  sd_vol <- sd(d[[trait]][d$flight_status == "volant"])
  dfree <- nrow(d) - nrow(s)
  out <- lapply(TERMS, function(tm) {
    if (!tm %in% rownames(s)) return(NULL)
    est <- s[tm, "Estimate"]; se <- s[tm, "StdErr"]
    tibble(term = tm, std_estimate = est / sd_vol, p = s[tm, "p.value"],
           delta_min_sd = (qt(0.95, dfree) * se - est) / sd_vol)
  })
  bind_rows(out) |> mutate(n = nrow(d),
                           n_terrestrial = sum(d$flight_status == "terrestrial"),
                           n_diving = sum(d$flight_status == "diving"))
}

res <- list()
add <- function(dn, tr, analysis, variant, d, rhs = M1) {
  res[[length(res) + 1]] <<- fit_terms(d, tr, rhs) |>
    mutate(dataset = dn, trait = tr, analysis = analysis, variant = variant, .before = 1)
}

for (dn in names(DATA)) {
  d0 <- DATA[[dn]]$d
  origins <- pairs |> filter(dataset == dn)
  vol_max <- max(d0$mass_g[d0$flight_status == "volant"])
  for (tr in DATA[[dn]]$traits) {
    add(dn, tr, "baseline", "M1, all species", d0)
    # A. leave one origin out
    for (k in seq_len(nrow(origins))) {
      drop <- strsplit(origins$flightless_taxa[k], ", ")[[1]]
      add(dn, tr, "leave_one_origin_out", paste("without", origins$origin_id[k]),
          d0 |> filter(!birdtree_name %in% drop))
    }
    # B. classification
    if ("Rhynochetos_jubatus" %in% d0$birdtree_name) {
      add(dn, tr, "classification", "kagu as volant",
          d0 |> mutate(flight_status = replace(flight_status,
                                               birdtree_name == "Rhynochetos_jubatus", "volant")))
    }
    if ("Gallinula_mortierii" %in% d0$birdtree_name) {
      add(dn, tr, "classification", "without Tasmanian nativehen",
          d0 |> filter(birdtree_name != "Gallinula_mortierii"))
    }
    # C. body mass
    add(dn, tr, "body_mass", "spline mass and egg size (3 df)", d0, SPLINE)
    add(dn, tr, "body_mass", "flightless within volant mass range",
        d0 |> filter(flight_status == "volant" | mass_g <= vol_max))
  }
}
res <- bind_rows(res) |> mutate(across(where(is.double), ~ signif(.x, 4)))
write_table(res, "sensitivity_pgls")
print(as.data.frame(res |> filter(term == "flight_statusterrestrial") |>
                      select(dataset, trait, variant, std_estimate, p, delta_min_sd, n_terrestrial)))
