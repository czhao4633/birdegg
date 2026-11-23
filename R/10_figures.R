# 10_figures.R
# Publication figures (PDF + 300-dpi PNG in output/figures).
#   Fig. 1  egg-shape space, coloured by flight status (both datasets)
#   Fig. 2  MCC tree (Montgomerie species) + trait heatmap, flight losses marked
#   Fig. 3  sister-pair comparisons: (a) one line per pair; (b) differences vs
#           age-matched volant sister-clade contrasts
#   Fig. 4  robustness across 100 posterior trees: PGLS effects and rate ratios
#
# Colours: Okabe & Ito (2008) colour-blind-safe set for flight status
# (volant grey, terrestrial vermillion, diving blue); heatmap uses the
# ColorBrewer PRGn diverging scale (colour-blind safe) with grey midpoint.
# Deterministic: no random numbers except the fixed-seed jitter in Fig. 4.

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtree)
  library(ggrepel)
  library(patchwork)
})
set.seed(SEED)

STATUS_COL <- c(volant = "#9A9A9A", terrestrial = "#D55E00", diving = "#0072B2")
TRAIT_LAB <- c(
  elongation  = "Elongation",
  pointedness = "Pointedness",
  polar_asym  = "Polar asymmetry",
  ellipticity = "Ellipticity",
  asymmetry   = "Asymmetry"
)
# Short English labels for each flightless origin and its sister group
ORIGIN_LAB <- c(
  Ratites               = "Ratites vs tinamous",
  Penguins              = "Penguins vs tubenoses",
  Strigops_habroptila   = "Kakapo vs other parrots",
  Rhynochetos_jubatus   = "Kagu vs sunbittern",
  Gallinula_mortierii   = "Tasmanian vs black-tailed nativehen",
  Phalacrocorax_harrisi = "Flightless vs double-crested cormorant"
)

theme_pub <- function(base = 8) {
  theme_classic(base_size = base) +
    theme(axis.line = element_line(linewidth = 0.3, colour = "grey30"),
          axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
          axis.text = element_text(colour = "grey20"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0),
          legend.key.size = unit(3, "mm"),
          plot.tag = element_text(face = "bold"))
}

save_fig <- function(p, name, w_mm, h_mm) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(PATHS$figures, paste0(name, ".", ext)), p,
           width = w_mm, height = h_mm, units = "mm", dpi = 300,
           device = if (ext == "pdf") cairo_pdf else "png", bg = "white")
  }
  log_msg("Wrote ", name, " (.pdf, .png)")
}

mont <- readRDS(file.path(PATHS$processed, "analysis_montgomerie.rds"))
stod <- readRDS(file.path(PATHS$processed, "analysis_stoddard.rds"))
mcc  <- readRDS(file.path(PATHS$processed, "tree_mcc.rds"))

# Common names for the flightless origins (single-species origins by species)
COMMON <- c(Ratites = "Ratites", Penguins = "Penguins",
            Strigops_habroptila = "Kakapo", Rhynochetos_jubatus = "Kagu",
            Gallinula_mortierii = "Tasmanian nativehen",
            Phalacrocorax_harrisi = "Gal\u00e1pagos cormorant",
            Atlantisia_rogersi = "Inaccessible Island rail")

# Label for each flightless species = the origin it belongs to
origin_labels <- function(d) {
  tree <- keep.tip(mcc, d$birdtree_name)
  cl <- flightless_clades(tree, d$birdtree_name[d$flightless])
  lab <- sapply(cl, function(x) {
    if (any(grepl("^(Struthio|Rhea|Dromaius|Casuarius|Apteryx)_", x))) "Ratites"
    else if (any(grepl("^Spheniscus|^Eudyptes|^Aptenodytes|^Pygoscelis", x))) "Penguins"
    else x[1]
  })
  lab <- unname(COMMON[lab])
  stopifnot(!anyNA(lab))
  tibble(birdtree_name = unlist(cl), origin = rep(lab, lengths(cl)))
}

# ---- Fig. 1: egg-shape space ------------------------------------------------------
shape_panel <- function(d, x, y, title, label = TRUE) {
  d <- d |> left_join(origin_labels(d), by = "birdtree_name") |>
    arrange(flight_status != "volant")          # flightless drawn on top
  cent <- d |> filter(flightless) |>
    group_by(origin, flight_status) |>
    summarise(across(all_of(c(x, y)), mean), n = n(), .groups = "drop") |>
    mutate(lab = if_else(n > 1, paste0(origin, " (", n, ")"), origin))
  ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_point(data = filter(d, !flightless), aes(fill = flight_status),
               shape = 21, colour = STATUS_COL[["volant"]], stroke = 0.1,
               size = 0.8, alpha = 0.45) +
    geom_point(data = filter(d, flightless), aes(fill = flight_status),
               shape = 21, size = 1.9, colour = "white", stroke = 0.3) +
    { if (label) geom_text_repel(data = cent, aes(label = lab), size = 2.2,
                    colour = "grey15", min.segment.length = 0, segment.size = 0.25,
                    segment.colour = "grey40", box.padding = 0.45, force = 3,
                    seed = SEED, max.overlaps = Inf) } +
    scale_fill_manual(values = STATUS_COL, drop = FALSE, name = NULL,
                      labels = c(volant = "Volant", terrestrial = "Flightless, terrestrial",
                                 diving = "Flightless, diving")) +
    labs(x = TRAIT_LAB[x], y = TRAIT_LAB[y], subtitle = title) +
    theme_pub()
}
legend_fix <- guides(fill = guide_legend(override.aes = list(
  size = 2.2, alpha = 1, colour = "white")))
f1a <- shape_panel(mont, "elongation", "pointedness",
                   paste0("Montgomerie et al. (", nrow(mont), " spp.)"))
f1b <- shape_panel(mont, "elongation", "polar_asym",
                   paste0("Montgomerie et al. (", nrow(mont), " spp.)"), label = FALSE)
f1c <- shape_panel(stod, "asymmetry", "ellipticity",
                   paste0("Stoddard et al. 2017 (", nrow(stod), " spp.)"))
fig1 <- (f1a + f1b + f1c) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") & legend_fix &
  theme(legend.position = "bottom")
save_fig(fig1, "fig1_egg_shape_space", 180, 75)

# ---- Fig. 2: tree + heatmap ---------------------------------------------------------
tree <- keep.tip(mcc, mont$birdtree_name)
tree$node.label <- NULL
fl_tips <- mont$birdtree_name[mont$flightless]
clades <- flightless_clades(tree, fl_tips)
origin_node <- sapply(clades, function(cl)
  if (length(cl) == 1) which(tree$tip.label == cl) else getMRCA(tree, cl))
origin_status <- as.character(mont$flight_status[match(sapply(clades, `[`, 1),
                                                       mont$birdtree_name)])
# Branch status: every edge inside (and the stem of) a flightless clade
edge_status <- rep("volant", Ntip(tree) + tree$Nnode)
for (k in seq_along(clades)) {
  nodes <- c(origin_node[k], if (origin_node[k] > Ntip(tree))
    phangorn::Descendants(tree, origin_node[k], "all"))
  edge_status[nodes] <- origin_status[k]
}
node_df <- tibble(node = seq_along(edge_status),
                  status = factor(edge_status, levels = names(STATUS_COL)))
origin_df <- tibble(node = origin_node, status = factor(origin_status, names(STATUS_COL)),
                    lab = origin_labels(mont)$origin[match(sapply(clades, `[`, 1),
                                                           origin_labels(mont)$birdtree_name)])

heat <- mont |>
  transmute(birdtree_name,
            across(c(elongation, pointedness, polar_asym),
                   ~ (.x - mean(.x[!flightless])) / sd(.x[!flightless]))) |>
  as.data.frame()
rownames(heat) <- heat$birdtree_name
heat <- heat[, -1]
colnames(heat) <- c("EL", "PT", "PA")
heat[] <- lapply(heat, function(v) pmax(pmin(v, 3), -3))   # clamp to +/- 3 SD

p2 <- ggtree(tree, layout = "circular", linewidth = 0.15) %<+% node_df +
  aes(colour = status) +
  scale_colour_manual(values = STATUS_COL, guide = "none")
# Origin markers: a point at the start of each flightless stem branch
p2 <- p2 + geom_point2(data = function(x) {
  x |> inner_join(origin_df |> select(node), by = "node") |>
    mutate(x = x - branch.length / 2)
}, aes(x = x, y = y), colour = "black", fill = "white", shape = 21, size = 1.6,
stroke = 0.4, inherit.aes = FALSE)
p2 <- gheatmap(p2, heat, offset = 2, width = 0.22, colnames_angle = 0,
               colnames_offset_y = 0, font.size = 2, hjust = 0.5, color = NA) +
  scale_fill_distiller(palette = "PRGn", direction = 1, limits = c(-3, 3),
                       name = "SD from\nvolant mean", na.value = "grey90")
# Clade labels for each origin, outside the heatmap
# The two rail origins (Atlantisia, Gallinula mortierii) sit next to each other
# on the circle; they share one label to avoid overprinting.
rails <- origin_df$lab %in% c("Inaccessible Island rail", "Tasmanian nativehen")
origin_df$lab[rails] <- c("Two rails", "")[seq_len(sum(rails))]
for (k in seq_along(origin_node)) {
  p2 <- p2 + geom_cladelab(node = origin_node[k], label = origin_df$lab[k],
                           offset = 36, barsize = 1.2, fontsize = 2.4,
                           barcolour = STATUS_COL[origin_status[k]],
                           textcolour = "grey15", angle = "auto", hjust = 0,
                           offset.text = 2, extend = 2)
}
# extra radial room so the outer clade labels are not clipped
p2 <- p2 + xlim(NA, max(p2$data$x, na.rm = TRUE) * 1.5) +
  theme(legend.position = c(0.06, 0.08),
                 legend.title = element_text(size = 6),
                 legend.text = element_text(size = 6),
                 legend.key.size = unit(3, "mm"))
save_fig(p2, "fig2_tree_heatmap", 180, 180)

# ---- Fig. 3: sister pairs -------------------------------------------------------------
pd <- read.csv(file.path(PATHS$tables, "paired_differences.csv"))
nulls <- read.csv(file.path(PATHS$tables, "paired_null_percentiles.csv"))
sd_vol <- sapply(c("elongation", "pointedness", "polar_asym"),
                 function(t) sd(mont[[t]][!mont$flightless]))
mean_vol <- sapply(names(sd_vol), function(t) mean(mont[[t]][!mont$flightless]))
pair_order <- pd |> distinct(origin_id, divergence_myr) |>
  arrange(divergence_myr) |> pull(origin_id)
pd <- pd |>
  mutate(fl_z = (fl - mean_vol[trait]) / sd_vol[trait],
         sis_z = (sis - mean_vol[trait]) / sd_vol[trait],
         pair = factor(ORIGIN_LAB[origin_id], ORIGIN_LAB[pair_order]),
         # very short labels: the slope panels are narrow
         short = recode(origin_id, Gallinula_mortierii = "Nativehen",
                        Phalacrocorax_harrisi = "Cormorant",
                        Strigops_habroptila = "Kakapo", Rhynochetos_jubatus = "Kagu"),
         pair_age = factor(sprintf("%s (%.0f Myr)", ORIGIN_LAB[origin_id], divergence_myr),
                           unique(sprintf("%s (%.0f Myr)", ORIGIN_LAB[pair_order],
                             divergence_myr[match(pair_order, origin_id)]))),
         trait_lab = factor(TRAIT_LAB[trait], TRAIT_LAB[1:3]))

slope <- pd |>
  select(pair, short, locomotion, trait_lab, `Volant sister` = sis_z, Flightless = fl_z) |>
  pivot_longer(c(`Volant sister`, Flightless), names_to = "side", values_to = "z") |>
  mutate(side = factor(side, c("Volant sister", "Flightless")))
f3a <- ggplot(slope, aes(side, z, group = pair, colour = locomotion)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.6) +
  geom_text_repel(data = filter(slope, side == "Flightless"),
                  aes(label = short), size = 2.2, direction = "y", hjust = 0,
                  nudge_x = 0.25, segment.size = 0.2, segment.colour = "grey60",
                  colour = "grey20", seed = SEED, xlim = c(2.1, NA)) +
  facet_wrap(~ trait_lab, nrow = 1) +
  scale_colour_manual(values = STATUS_COL[c("terrestrial", "diving")],
                      breaks = c("terrestrial", "diving"),
                      labels = c(terrestrial = "Terrestrial", diving = "Diving"),
                      name = "Flightless type") +
  scale_x_discrete(expand = expansion(add = c(0.25, 1.1))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Trait value (SD units, relative to volant mean)") +
  theme_pub()

forest <- pd |>
  left_join(nulls |> select(origin_id, trait, null_abs_q90, n_null),
            by = c("origin_id", "trait"))
f3b <- ggplot(forest, aes(diff_std, pair_age)) +
  geom_linerange(aes(xmin = -null_abs_q90, xmax = null_abs_q90),
                 linewidth = 2.2, colour = "grey88") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_point(aes(colour = locomotion), size = 2) +
  facet_wrap(~ trait_lab, nrow = 1) +
  scale_colour_manual(values = STATUS_COL[c("terrestrial", "diving")], guide = "none") +
  labs(x = paste0("Flightless \u2212 sister difference (SD units); H1 predicts < 0 for terrestrial pairs\n",
                  "Grey bars: central 90% of age-matched contrasts between volant sister clades"),
       y = NULL) +
  theme_pub()
fig3 <- f3a / f3b + plot_layout(heights = c(1.15, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
save_fig(fig3, "fig3_sister_pairs", 180, 150)

# ---- Fig. 4: posterior robustness ------------------------------------------------------
pg <- readRDS(file.path(PATHS$processed, "pgls_posterior.rds")) |>
  filter(model == "M1_status",
         term %in% c("flight_statusterrestrial", "flight_statusdiving")) |>
  mutate(status = factor(if_else(term == "flight_statusterrestrial", "terrestrial", "diving"),
                         c("terrestrial", "diving")),
         trait_lab = factor(paste0(TRAIT_LAB[trait], "\n(",
                                   if_else(dataset == "montgomerie", "M", "S"), ")")))
lvl <- levels(pg$trait_lab)
lvl <- c(grep("\\(M\\)", lvl, value = TRUE), grep("\\(S\\)", lvl, value = TRUE))
pg$trait_lab <- factor(pg$trait_lab, lvl)
pgm <- read.csv(file.path(PATHS$tables, "pgls_main.csv")) |>
  filter(model == "M1_status",
         term %in% c("flight_statusterrestrial", "flight_statusdiving")) |>
  mutate(status = factor(if_else(term == "flight_statusterrestrial", "terrestrial", "diving"),
                         c("terrestrial", "diving")),
         trait_lab = factor(paste0(TRAIT_LAB[trait], "\n(",
                                   if_else(dataset == "montgomerie", "M", "S"), ")"), lvl))
# share of posterior trees with p < 0.05, printed above each group
sig <- pg |> group_by(trait_lab, status) |>
  summarise(lab = sprintf("%d%%", round(100 * mean(p < 0.05))), .groups = "drop")
dodge <- position_dodge(width = 0.6)
f4a <- ggplot(pg, aes(trait_lab, as.numeric(std_estimate), colour = status)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_vline(xintercept = 3.5, linewidth = 0.25, linetype = "dashed", colour = "grey60") +
  geom_text(data = sig, aes(y = 0.75, label = lab), position = dodge, size = 2.1,
            show.legend = FALSE) +
  annotate("text", x = 0.45, y = 0.92, hjust = 0, size = 2.1, colour = "grey30",
           label = "Share of trees with p < 0.05:") +
  annotate("text", x = 5.45, y = -2.05, hjust = 1, size = 2.1, colour = "grey30",
           label = "S terrestrial effects: covariate-dependent, few species (see text)") +
  geom_point(position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.6,
                                             seed = SEED),
             size = 0.6, alpha = 0.5, stroke = 0) +
  geom_point(data = pgm, aes(y = std_estimate), position = dodge, shape = 23,
             size = 2, fill = "white", stroke = 0.6) +
  scale_colour_manual(values = STATUS_COL[c("terrestrial", "diving")],
                      labels = c(terrestrial = "Terrestrial", diving = "Diving"),
                      name = "Flightless type") +
  labs(x = NULL, y = "PGLS effect vs volant (SD units)",
       subtitle = "Model M1, 100 posterior trees (dots) and MCC tree (diamonds)") +
  theme_pub()

rt <- readRDS(file.path(PATHS$processed, "rates_posterior.rds")) |>
  filter(converged) |>
  transmute(dataset, trait, regimes,
            `Flightless (R2)` = unname(ratio_flightless),
            `Terrestrial (R3)` = unname(ratio_terrestrial),
            `Diving (R3)` = unname(ratio_diving)) |>
  pivot_longer(-c(dataset, trait, regimes), names_to = "regime", values_to = "ratio") |>
  filter(!is.na(ratio)) |>
  mutate(regime = factor(regime, c("Flightless (R2)", "Terrestrial (R3)", "Diving (R3)")),
         trait_lab = factor(paste0(TRAIT_LAB[trait], "\n(",
                                   if_else(dataset == "montgomerie", "M", "S"), ")"), lvl))
REG_COL <- c(`Flightless (R2)` = "grey25", `Terrestrial (R3)` = STATUS_COL[["terrestrial"]],
             `Diving (R3)` = STATUS_COL[["diving"]])
f4b <- ggplot(rt, aes(trait_lab, ratio, colour = regime)) +
  geom_hline(yintercept = 1, linewidth = 0.25, colour = "grey50") +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.shape = NA,
               linewidth = 0.3, fill = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75,
                                             seed = SEED),
             size = 0.45, alpha = 0.4, stroke = 0) +
  scale_y_log10(breaks = c(0.03, 0.1, 0.3, 1, 3, 10, 30), labels = c(0.03, 0.1, 0.3, 1, 3, 10, 30)) +
  scale_colour_manual(values = REG_COL, name = "Regime") +
  labs(x = NULL, y = "Rate ratio vs volant (log scale)",
       subtitle = "Multi-rate Brownian motion, converged fits on 100 posterior trees") +
  theme_pub()
fig4 <- f4a / f4b + plot_annotation(tag_levels = "a",
  caption = "M = Montgomerie et al. dataset (Biggins metrics); S = Stoddard et al. 2017") &
  theme(legend.position = "right", plot.caption = element_text(size = 6, colour = "grey30"))
save_fig(fig4, "fig4_posterior_robustness", 180, 140)
