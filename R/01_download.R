# 01_download.R
# Download the source datasets into data/raw/ and verify their checksums.
#
# Sources (URLs verified 2026-11-17):
#   AVONET (Tobias et al. 2022, Ecol. Lett. 25:581-597, doi:10.1111/ele.13898)
#     figshare article 16586228, file "AVONET Supplementary dataset 1.xlsx"
#   BirdTree (Jetz et al. 2012, Nature 491:444-448, doi:10.1038/nature11631)
#     data.vertlife.org GCS bucket: Hackett Stage 2 and Stage 1, trees 1-1000,
#     plus the master taxonomy BLIOCPhyloMasterTax.csv
#   Egg shape (Montgomerie et al.), Zenodo record 5826243 (mirror of Dryad)
#   Egg shape (Stoddard et al. 2017, Science 356:1249-1254,
#     doi:10.1126/science.aaj1945), Data S1 v2.
#     science.org refuses scripted downloads, so this file must be placed in
#     data/raw/stoddard2017/ by hand (see README). This script only checks it.

source(here::here("R", "00_setup.R"))

options(timeout = 3600)  # tree archives are ~200 MB

GCS <- "https://storage.googleapis.com/data.vertlife.org/birdtree/"
# Montgomerie, Birkhead, Thompson & Hemmings, Zenodo/Dryad
# doi:10.5061/dryad.sj3tx9648 (CC0): egg shape (Biggins et al. 2018 metrics)
# for 955 species, with clutch size and incubation-site structure.
ZEN <- "https://zenodo.org/api/records/5826243/files/"

# expected_md5: hex MD5 published by the host (figshare: hex; GCS: base64,
# converted to hex below). NA = no published checksum.
b64_to_hex <- function(b64) {
  paste(as.character(jsonlite::base64_dec(b64)), collapse = "")
}

files <- data.frame(
  dest = c(
    "avonet/AVONET_Supplementary_dataset_1.xlsx",
    "birdtree/BLIOCPhyloMasterTax.csv",
    "birdtree/HackettStage2_0001_1000.zip",
    "birdtree/HackettStage1_0001_1000.zip",
    "montgomerie2022/README.txt",
    "montgomerie2022/DATAeggshape.csv",
    "montgomerie2022/DATAraweggs.csv"
  ),
  url = c(
    "https://ndownloader.figshare.com/files/34480856",
    paste0(GCS, "BLIOCPhyloMasterTax.csv"),
    paste0(GCS, "Stage2/HackettStage2_0001_1000.zip"),
    paste0(GCS, "Stage1/HackettStage1_0001_1000.zip"),
    paste0(ZEN, "README.txt/content"),
    paste0(ZEN, "DATAeggshape.csv/content"),
    paste0(ZEN, "DATAraweggs.csv/content")
  ),
  expected_md5 = c(
    "1445afdcfb6df784010c2ca034544bc8",
    b64_to_hex("EERTMO6UNBSibb4IB+9zLQ=="),
    b64_to_hex("JRS61ik+7z9UHQ2cvYRVhQ=="),
    b64_to_hex("e8Urmzc46tM4ROck0dIV0A=="),
    "c110b8a55328f0bcb6639cb24e8726ae",
    "c6169eb05d7567aef84215de0cc4e6d8",
    "2139b639fda518afa5029723c05ee9b0"
  ),
  stringsAsFactors = FALSE
)

download_one <- function(dest, url, expected_md5) {
  path <- file.path(PATHS$raw, dest)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(path) && unname(tools::md5sum(path)) == expected_md5) {
    log_msg("Already present and verified: ", dest)
  } else {
    log_msg("Downloading ", dest)
    utils::download.file(url, path, mode = "wb", quiet = TRUE)
  }
  md5 <- unname(tools::md5sum(path))
  if (!is.na(expected_md5) && md5 != expected_md5) {
    stop("MD5 mismatch for ", dest, ": got ", md5, ", expected ", expected_md5)
  }
  data.frame(file = dest, url = url, md5 = md5,
             bytes = file.size(path), retrieved = as.character(Sys.Date()))
}

manifest <- do.call(rbind, Map(download_one, files$dest, files$url, files$expected_md5))

# Unzip the tree archives (each contains one NEXUS file with 1,000 trees)
for (z in grep("\\.zip$", files$dest, value = TRUE)) {
  zp <- file.path(PATHS$raw, z)
  exdir <- sub("\\.zip$", "", zp)
  # junkpaths: the archives store files under a long absolute server path
  if (!dir.exists(exdir)) utils::unzip(zp, exdir = exdir, junkpaths = TRUE)
}

# ---- Manually supplied file -------------------------------------------------
stoddard_dir <- file.path(PATHS$raw, "stoddard2017")
# Create every manual-download folder so users know where files go
for (d in c("stoddard2017", "sayol2020", "saitta2025")) {
  dir.create(file.path(PATHS$raw, d), showWarnings = FALSE)
}
stoddard_file <- file.path(stoddard_dir, "aaj1945_datas1_egg_shape_by_species_v2.xlsx")
if (file.exists(stoddard_file)) {
  manifest <- rbind(manifest, data.frame(
    file = "stoddard2017/aaj1945_datas1_egg_shape_by_species_v2.xlsx",
    url = "https://www.science.org/doi/10.1126/science.aaj1945 (Supplementary Materials, Data S1 v2; manual download)",
    md5 = unname(tools::md5sum(stoddard_file)),
    bytes = file.size(stoddard_file), retrieved = as.character(as.Date(file.mtime(stoddard_file)))
  ))
} else {
  warning(
    "Stoddard et al. 2017 Data S1 not found. Download ",
    "'aaj1945_datas1_egg_shape_by_species_v2.xlsx' from the Supplementary ",
    "Materials of doi:10.1126/science.aaj1945 and save it to ", stoddard_file
  )
}

# Other manual downloads (publisher / Dryad block scripted access)
manual <- list(
  list(file = "sayol2020/abb6095_data_file_s2.xlsx",
       url = "https://doi.org/10.1126/sciadv.abb6095 (Supplementary Materials, Data file S2; manual download)"),
  list(file = "saitta2025/Supplementary_Data_S1-S5.zip",
       url = "https://doi.org/10.5061/dryad.qrfj6q5s7 (manual download)")
)
for (m in manual) {
  f <- file.path(PATHS$raw, m$file)
  if (!file.exists(f)) {
    warning("Manual download missing: ", f, " (source: ", m$url, "; see README)")
    next
  }
  manifest <- rbind(manifest, data.frame(
    file = m$file, url = m$url, md5 = unname(tools::md5sum(f)),
    bytes = file.size(f), retrieved = as.character(as.Date(file.mtime(f)))))
}

# Extract the two Saitta et al. files used downstream (03, 05). The archive
# also contains macOS "__MACOSX/._*" resource forks, which are skipped.
saitta_zip <- file.path(PATHS$raw, "saitta2025", "Supplementary_Data_S1-S5.zip")
if (file.exists(saitta_zip)) {
  members <- utils::unzip(saitta_zip, list = TRUE)$Name
  want <- members[!grepl("__MACOSX", members) &
                    grepl("(S1_Dataset/Dataset\\.csv|data/sisterpairs\\.csv)$", members)]
  stopifnot(length(want) == 2)
  utils::unzip(saitta_zip, files = want, junkpaths = TRUE,
               exdir = file.path(PATHS$raw, "saitta2025", "extracted"))
}

readr::write_csv(manifest, file.path(PATHS$raw, "MANIFEST.csv"))
print(manifest[, c("file", "bytes", "md5")])
