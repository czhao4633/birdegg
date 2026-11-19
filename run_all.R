# run_all.R
# Reproduce the full analysis: Rscript run_all.R
# Each numbered script in R/ is self-contained and is run in a fresh R
# process so that no state leaks between steps.

scripts <- sort(list.files("R", pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
scripts <- setdiff(scripts, file.path("R", "00_setup.R"))

for (s in scripts) {
  message("\n==== Running ", s, " ====")
  status <- system2("Rscript", s)
  if (status != 0) stop("Script failed: ", s)
}

# Render the report if Quarto is available
if (nzchar(Sys.which("quarto"))) {
  system2("quarto", c("render", file.path("report", "report.qmd")))
} else {
  message("Quarto not found; skipping report rendering.")
}
