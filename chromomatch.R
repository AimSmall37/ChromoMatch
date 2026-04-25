#!/usr/bin/env Rscript
###############################################################################
# ChromoMatch
#
# Compares two autosomal DNA samples and finds shared IBD
# (identical-by-descent) segments reported in centimorgans (cM).
#
# Input sources:
#   1. Raw DNA files from 23andMe, AncestryDNA, FTDNA
#   2. PLINK binary datasets (.bed/.bim/.fam) — no external tools required
#
# Genetic maps: PLINK (.map), HapMap recombination maps, or auto-download
# Reference build: GRCh37 / hg19 (default)
#
# Usage (command line — raw files):
#   Rscript chromomatch.R <file1> <file2> [options]
#
# Usage (command line — PLINK dataset):
#   Rscript chromomatch.R --bfile <prefix> --sample1 <ID> --sample2 <ID>
#   Rscript chromomatch.R --bfile <prefix>          # interactive sample picker
#
# Usage (RStudio):
#   Uncomment and set variables in the CONFIGURATION section below,
#   then select all and click Run (or source the file).
#
# Options:
#   --genetic-map <dir>    Genetic map directory (default: auto-download)
#   --min-cm <num>         Minimum segment size in cM (default: 7)
#   --min-snps <num>       Minimum SNPs per segment (default: 500)
#   --mismatch-bunch <num> Sliding window size in SNPs (default: 250)
#   --seed-threshold <num> Max mismatch rate to seed IBD (default: 0.02)
#   --output <path>        Output CSV path (default: comparison_results.csv)
#   --plot / --no-plot     Generate chromosome plot (default: yes)
#   --plot-file <path>     Plot file path (default: comparison_plot.png)
#   --bfile <prefix>       PLINK binary prefix (.bed/.bim/.fam)
#   --sample1 <ID>         First sample IID (from .fam file)
#   --sample2 <ID>         Second sample IID (from .fam file)
#
# License: GPL-3.0
###############################################################################

# =============================================================================
# DEPENDENCIES
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table", repos = "https://cloud.r-project.org")
  }
  library(data.table)
})

# =============================================================================
# CONFIGURATION (RStudio Interactive Mode)
# =============================================================================
# Uncomment and edit the file paths below to run interactively in RStudio.
# Use EITHER file1/file2 (raw DNA) OR bfile (PLINK). Not both.
#
# --- Raw DNA files ---
# file1           <- "C:/path/to/person1_raw_dna.txt"
# file2           <- "C:/path/to/person2_raw_dna.txt"
#
# --- PLINK dataset ---
# bfile           <- "C:/path/to/plink_dataset"   # prefix for .bed/.bim/.fam
# sample1         <- "SampleA"                     # IID from .fam (or leave
# sample2         <- "SampleB"                     #   unset for interactive)
#
# --- Shared options ---
# genetic_map_dir <- "./genetic_map_grch37"
# min_cm          <- 7
# min_snps        <- 500
# mismatch_bunch  <- 250
# seed_threshold  <- 0.02
# output_file     <- "comparison_results.csv"
# plot_output     <- TRUE
# plot_file       <- "comparison_plot.png"

# =============================================================================
# DEFAULT PARAMETERS
# =============================================================================

defaults <- list(
  min_cm          = 7,       # Minimum segment cM to report
  min_snps        = 500,     # Minimum SNPs per segment
  mismatch_bunch  = 250,     # Sliding window size (SNPs)
  seed_threshold  = 0.02,    # Max mismatch rate to seed an IBD region
  genetic_map_dir = NULL,    # NULL = attempt auto-download
  output_file     = "comparison_results.csv",
  plot            = TRUE,
  plot_file       = "comparison_plot.png",
)

# =============================================================================
# ARGUMENT HANDLING
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
input_mode <- "none"   # will be set to "raw", "plink", or trigger usage

if (length(args) >= 1) {
  # ---- Command-line mode: scan for --bfile first ----
  bfile_idx <- match("--bfile", args)

  if (!is.na(bfile_idx)) {
    # PLINK mode
    input_mode <- "plink"
    plink_prefix <- args[bfile_idx + 1]
    plink_sample1 <- NULL
    plink_sample2 <- NULL

    # Parse remaining flags (skip the positional --bfile value)
    i <- 1
    while (i <= length(args)) {
      if (i == bfile_idx) { i <- i + 2; next }
      switch(args[i],
        "--sample1"        = { plink_sample1 <- args[i+1]; i <- i+2 },
        "--sample2"        = { plink_sample2 <- args[i+1]; i <- i+2 },
        "--genetic-map"    = { defaults$genetic_map_dir <- args[i+1]; i <- i+2 },
        "--min-cm"         = { defaults$min_cm <- as.numeric(args[i+1]); i <- i+2 },
        "--min-snps"       = { defaults$min_snps <- as.integer(args[i+1]); i <- i+2 },
        "--mismatch-bunch" = { defaults$mismatch_bunch <- as.integer(args[i+1]); i <- i+2 },
        "--seed-threshold" = { defaults$seed_threshold <- as.numeric(args[i+1]); i <- i+2 },
        "--output"         = { defaults$output_file <- args[i+1]; i <- i+2 },
        "--plot"           = { defaults$plot <- TRUE; i <- i+1 },
        "--no-plot"        = { defaults$plot <- FALSE; i <- i+1 },
        "--plot-file"      = { defaults$plot_file <- args[i+1]; i <- i+2 },
        { cat("Unknown option:", args[i], "\n"); i <- i+1 }
      )
    }

  } else if (length(args) >= 2 && !startsWith(args[1], "--")) {
    # Raw file mode (two positional arguments)
    input_mode <- "raw"
    file1 <- args[1]
    file2 <- args[2]
    i <- 3
    while (i <= length(args)) {
      switch(args[i],
        "--genetic-map"    = { defaults$genetic_map_dir <- args[i+1]; i <- i+2 },
        "--min-cm"         = { defaults$min_cm <- as.numeric(args[i+1]); i <- i+2 },
        "--min-snps"       = { defaults$min_snps <- as.integer(args[i+1]); i <- i+2 },
        "--mismatch-bunch" = { defaults$mismatch_bunch <- as.integer(args[i+1]); i <- i+2 },
        "--seed-threshold" = { defaults$seed_threshold <- as.numeric(args[i+1]); i <- i+2 },
        "--output"         = { defaults$output_file <- args[i+1]; i <- i+2 },
        "--plot"           = { defaults$plot <- TRUE; i <- i+1 },
        "--no-plot"        = { defaults$plot <- FALSE; i <- i+1 },
        "--plot-file"      = { defaults$plot_file <- args[i+1]; i <- i+2 },
        { cat("Unknown option:", args[i], "\n"); i <- i+1 }
      )
    }
  }
}

# ---- RStudio interactive mode fallback ----
if (input_mode == "none") {
  if (exists("bfile") && !is.null(get("bfile", envir = .GlobalEnv))) {
    input_mode <- "plink"
    plink_prefix <- get("bfile", envir = .GlobalEnv)
    plink_sample1 <- if (exists("sample1")) get("sample1", envir = .GlobalEnv) else NULL
    plink_sample2 <- if (exists("sample2")) get("sample2", envir = .GlobalEnv) else NULL
  } else if (exists("file1") && exists("file2")) {
    input_mode <- "raw"
    file1 <- get("file1", envir = .GlobalEnv)
    file2 <- get("file2", envir = .GlobalEnv)
  }

  # Pick up shared options from global env
  if (input_mode != "none") {
    if (exists("genetic_map_dir")) defaults$genetic_map_dir <- genetic_map_dir
    if (exists("min_cm"))          defaults$min_cm          <- min_cm
    if (exists("min_snps"))        defaults$min_snps        <- min_snps
    if (exists("mismatch_bunch"))  defaults$mismatch_bunch  <- mismatch_bunch
    if (exists("seed_threshold"))  defaults$seed_threshold  <- seed_threshold
    if (exists("output_file"))     defaults$output_file     <- output_file
    if (exists("plot_output"))     defaults$plot            <- plot_output
    if (exists("plot_file"))       defaults$plot_file       <- plot_file
  }
}

if (input_mode == "none") {
  message("
============================================================
ChromoMatch
============================================================

RStudio: Uncomment and set file1/file2 (or bfile) near the
         top of the script, then select all and click Run.

Command line (raw DNA files):
  Rscript chromomatch.R <file1> <file2> [options]

Command line (PLINK dataset):
  Rscript chromomatch.R --bfile <prefix> --sample1 ID1 --sample2 ID2
  Rscript chromomatch.R --bfile <prefix>   (interactive sample picker)

Options:
  --genetic-map DIR       Genetic map directory
  --min-cm NUM            Min segment cM (default: 7)
  --min-snps NUM          Min SNPs per segment (default: 500)
  --mismatch-bunch NUM    Sliding window size (default: 250)
  --seed-threshold NUM    IBD seed mismatch rate (default: 0.02)
  --output FILE           Output CSV (default: comparison_results.csv)
  --plot / --no-plot      Chromosome plot (default: yes)
  --plot-file FILE        Plot path (default: comparison_plot.png)
  --bfile PREFIX          PLINK binary dataset (.bed/.bim/.fam)
  --sample1 ID            First sample IID
  --sample2 ID            Second sample IID

Supported raw formats: 23andMe, AncestryDNA, FTDNA
============================================================")
  stop("Please provide input files before running.", call. = FALSE)
}

# =============================================================================
# RAW DNA FILE FORMAT DETECTION AND PARSING
# =============================================================================

detect_and_read_raw <- function(filepath) {
  con <- file(filepath, "r")
  header_lines <- c()
  first_data <- NULL

  while (TRUE) {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#")) {
      header_lines <- c(header_lines, line)
    } else {
      first_data <- line
      break
    }
  }
  close(con)

  fmt <- "unknown"
  build <- "37"
  header_text <- paste(header_lines, collapse = "\n")

  if (grepl("23andMe", header_text, ignore.case = TRUE)) {
    fmt <- "23andme"
    if (grepl("build 38|GRCh38", header_text)) build <- "38"
  } else if (grepl("AncestryDNA", header_text, ignore.case = TRUE)) {
    fmt <- "ancestry"
    if (grepl("build 38|GRCh38", header_text)) build <- "38"
  } else if (grepl('^"?RSID"?,', first_data, ignore.case = TRUE) ||
             grepl('^"rs', first_data)) {
    fmt <- "ftdna"
  } else if (grepl("^rsid\t", first_data, ignore.case = TRUE)) {
    fields <- strsplit(first_data, "\t")[[1]]
    if (length(fields) == 5) fmt <- "ancestry"
    else if (length(fields) == 4) fmt <- "23andme"
  }

  cat("  Format:", fmt, "| Build:", build, "\n")

  if (fmt == "23andme") {
    dt <- fread(filepath, sep = "\t", header = FALSE,
                skip = length(header_lines),
                colClasses = rep("character", 4), showProgress = FALSE)
    setnames(dt, c("rsid", "chromosome", "position", "genotype"))
    if (dt[1, rsid] == "rsid") dt <- dt[-1]
    dt[, position := as.integer(position)]

  } else if (fmt == "ancestry") {
    dt <- fread(filepath, sep = "\t", header = FALSE,
                skip = length(header_lines),
                colClasses = rep("character", 5), showProgress = FALSE)
    setnames(dt, c("rsid", "chromosome", "position", "allele1", "allele2"))
    if (dt[1, rsid] == "rsid") dt <- dt[-1]
    dt[, position := as.integer(position)]
    dt[, genotype := paste0(allele1, allele2)]
    dt[, c("allele1", "allele2") := NULL]

  } else if (fmt == "ftdna") {
    dt <- fread(filepath, sep = ",", header = FALSE,
                skip = max(1, length(header_lines)),
                colClasses = rep("character", 4), showProgress = FALSE)
    for (col in names(dt)) dt[, (col) := gsub('"', '', get(col))]
    if (toupper(dt[[1]][1]) == "RSID") dt <- dt[-1]
    setnames(dt, c("rsid", "chromosome", "position", "genotype"))
    dt[, position := as.integer(position)]

  } else {
    dt <- fread(filepath, sep = "auto", header = TRUE, showProgress = FALSE)
    if (ncol(dt) == 5) {
      setnames(dt, c("rsid", "chromosome", "position", "allele1", "allele2"))
      dt[, genotype := paste0(allele1, allele2)]
      dt[, c("allele1", "allele2") := NULL]
    } else if (ncol(dt) >= 4) {
      setnames(dt, 1:4, c("rsid", "chromosome", "position", "genotype"))
    }
  }

  # Standardize
  dt[, chromosome := toupper(gsub("chr", "", chromosome, ignore.case = TRUE))]
  dt[, position := as.integer(position)]
  dt <- dt[chromosome %in% as.character(1:22)]
  dt[, chromosome := as.integer(chromosome)]
  dt <- dt[genotype != "--" & genotype != "00" & genotype != "  " &
           !is.na(position) & nchar(genotype) == 2]
  setkey(dt, chromosome, position)

  cat("  Loaded", format(nrow(dt), big.mark = ","), "autosomal SNPs\n")

  fmt_display <- switch(fmt,
    "23andme"  = "23andMe",
    "ancestry" = "AncestryDNA",
    "ftdna"    = "FTDNA",
    "unknown"
  )

  return(list(data = dt, format = fmt, format_display = fmt_display,
              build = build, snp_count = nrow(dt)))
}

# =============================================================================
# PLINK BINARY (.bed/.bim/.fam) READER — pure R, no external tools
# =============================================================================

#' Read the .fam file and return a data.table of sample information
read_plink_fam <- function(fam_path) {
  fam <- fread(fam_path, header = FALSE, colClasses = rep("character", 6),
               showProgress = FALSE)
  setnames(fam, c("FID", "IID", "father", "mother", "sex", "pheno"))
  return(fam)
}

#' Read the .bim file and return a data.table of variant information
read_plink_bim <- function(bim_path) {
  bim <- fread(bim_path, header = FALSE, showProgress = FALSE)
  setnames(bim, c("chromosome", "rsid", "cM", "position", "allele1", "allele2"))
  bim[, chromosome := as.character(chromosome)]
  return(bim)
}

#' Resolve a sample identifier string against a .fam data.table.
#' Accepts "FID IID" (compound), bare "IID" (if unique), or row number.
#' Returns the 1-based row index into the fam table.
resolve_sample_id <- function(id_string, fam, label = "sample") {
  id_string <- trimws(id_string)

  # 1. Try numeric row index
  num <- suppressWarnings(as.integer(id_string))
  if (!is.na(num) && num >= 1 && num <= nrow(fam)) return(num)

  # 2. Try "FID IID" compound match (space-separated)
  parts <- strsplit(id_string, "\\s+")[[1]]
  if (length(parts) >= 2) {
    fid <- parts[1]
    iid <- paste(parts[-1], collapse = " ")  # IID may itself contain spaces
    idx <- which(fam$FID == fid & fam$IID == iid)
    if (length(idx) == 1) return(idx)
    if (length(idx) > 1) {
      cat("  WARNING: Multiple rows match FID='", fid, "' IID='", iid,
          "' — using first\n", sep = "")
      return(idx[1])
    }
    # If two-word input didn't match as FID+IID, fall through and try as bare IID
  }

  # 3. Try bare IID match
  idx <- which(fam$IID == id_string)
  if (length(idx) == 1) return(idx)
  if (length(idx) > 1) {
    cat("  ERROR: IID '", id_string, "' matches multiple samples:\n", sep = "")
    for (j in idx) cat("    ", fam$FID[j], " ", fam$IID[j], "\n")
    stop("Ambiguous ", label, ". Use 'FID IID' format, e.g. --", label,
         " \"", fam$FID[idx[1]], " ", id_string, "\"", call. = FALSE)
  }

  # 4. Try case-insensitive partial match on "FID IID" display string
  display <- paste(fam$FID, fam$IID)
  idx <- grep(id_string, display, ignore.case = TRUE)
  if (length(idx) == 1) return(idx)
  if (length(idx) > 1) {
    cat("  Multiple partial matches for '", id_string, "':\n", sep = "")
    for (j in idx) cat("    ", fam$FID[j], " ", fam$IID[j], "\n")
    stop("Ambiguous ", label, ". Please be more specific.", call. = FALSE)
  }

  stop(label, " not found in .fam: ", id_string, call. = FALSE)
}

#' Extract genotypes for exactly two samples from a PLINK .bed file.
#' Returns a list with two data.tables in the same format as detect_and_read_raw().
#'
#' sample1_id and sample2_id can be "FID IID", bare "IID" (if unique), or
#' row numbers. They are resolved via resolve_sample_id().
#'
#' BED format (variant-major, PLINK 1.9+):
#'   - 3-byte magic header: 0x6C 0x1B 0x01
#'   - Then for each variant: ceil(n_samples/4) bytes
#'   - Each byte encodes 4 genotypes in 2-bit pairs (low bits first):
#'       00 = homozygous allele1 (A1/A1)
#'       01 = missing
#'       10 = heterozygous (A1/A2)
#'       11 = homozygous allele2 (A2/A2)
read_plink_bed <- function(bfile_prefix, sample1_id, sample2_id) {
  bed_path <- paste0(bfile_prefix, ".bed")
  bim_path <- paste0(bfile_prefix, ".bim")
  fam_path <- paste0(bfile_prefix, ".fam")

  for (f in c(bed_path, bim_path, fam_path)) {
    if (!file.exists(f)) stop("File not found: ", f, call. = FALSE)
  }

  # Read metadata
  fam <- read_plink_fam(fam_path)
  bim <- read_plink_bim(bim_path)

  n_samples <- nrow(fam)
  n_variants <- nrow(bim)
  cat("  PLINK dataset:", format(n_variants, big.mark = ","), "variants,",
      n_samples, "samples\n")

  # Resolve sample identifiers to row indices
  idx1 <- resolve_sample_id(sample1_id, fam, "sample1")
  idx2 <- resolve_sample_id(sample2_id, fam, "sample2")
  sample1_label <- paste(fam$FID[idx1], fam$IID[idx1])
  sample2_label <- paste(fam$FID[idx2], fam$IID[idx2])
  cat("  Sample 1: [", idx1, "] ", sample1_label, "\n", sep = "")
  cat("  Sample 2: [", idx2, "] ", sample2_label, "\n", sep = "")

  # Compute byte positions for these two samples within each variant block
  bytes_per_variant <- as.integer(ceiling(n_samples / 4))

  # For sample at 1-based index idx:
  #   byte offset within variant block = floor((idx-1) / 4)
  #   bit shift within that byte       = ((idx-1) %% 4) * 2
  byte_offset1 <- as.integer(floor((idx1 - 1) / 4))
  bit_shift1   <- as.integer(((idx1 - 1) %% 4) * 2)
  byte_offset2 <- as.integer(floor((idx2 - 1) / 4))
  bit_shift2   <- as.integer(((idx2 - 1) %% 4) * 2)

  # Read entire BED file as raw bytes
  cat("  Reading BED file...\n")
  raw_bytes <- readBin(bed_path, what = "raw", n = file.info(bed_path)$size)

  # Validate magic number
  if (length(raw_bytes) < 3 ||
      raw_bytes[1] != as.raw(0x6C) ||
      raw_bytes[2] != as.raw(0x1B)) {
    stop("Not a valid PLINK BED file (bad magic number)", call. = FALSE)
  }
  if (raw_bytes[3] != as.raw(0x01)) {
    stop("Sample-major BED format not supported. ",
         "Convert with: plink --bfile ... --make-bed --out ...", call. = FALSE)
  }

  # Extract genotype codes for both samples across all variants
  # We work in integer arithmetic on the raw bytes
  header_size <- 3L
  geno1 <- integer(n_variants)
  geno2 <- integer(n_variants)

  # Vectorized extraction: compute the byte index for each variant
  variant_offsets <- seq(0L, by = bytes_per_variant, length.out = n_variants)

  # Byte indices (1-based) into raw_bytes for each sample at each variant
  byte_idx1 <- header_size + variant_offsets + byte_offset1 + 1L
  byte_idx2 <- header_size + variant_offsets + byte_offset2 + 1L

  # Extract and decode the 2-bit genotype codes
  geno1 <- bitwAnd(bitwShiftR(as.integer(raw_bytes[byte_idx1]), bit_shift1), 3L)
  geno2 <- bitwAnd(bitwShiftR(as.integer(raw_bytes[byte_idx2]), bit_shift2), 3L)

  # Decode: 0 = A1/A1, 1 = missing, 2 = het (A1/A2), 3 = A2/A2
  # Build genotype strings from allele1/allele2 in bim
  a1 <- bim$allele1  # "clear bits" allele
  a2 <- bim$allele2  # "set bits" allele

  genotype1 <- ifelse(geno1 == 0L, paste0(a1, a1),
               ifelse(geno1 == 2L, paste0(a1, a2),
               ifelse(geno1 == 3L, paste0(a2, a2), NA_character_)))
  genotype2 <- ifelse(geno2 == 0L, paste0(a1, a1),
               ifelse(geno2 == 2L, paste0(a1, a2),
               ifelse(geno2 == 3L, paste0(a2, a2), NA_character_)))

  # Build data.tables in the same format as raw DNA files
  dt1 <- data.table(
    rsid       = bim$rsid,
    chromosome = bim$chromosome,
    position   = bim$position,
    genotype   = genotype1
  )
  dt2 <- data.table(
    rsid       = bim$rsid,
    chromosome = bim$chromosome,
    position   = bim$position,
    genotype   = genotype2
  )

  # Standardize (same filtering as raw DNA path)
  standardize_plink_dt <- function(dt, label) {
    dt[, chromosome := toupper(gsub("chr", "", chromosome, ignore.case = TRUE))]
    dt[, position := as.integer(position)]
    dt <- dt[chromosome %in% as.character(1:22)]
    dt[, chromosome := as.integer(chromosome)]
    dt <- dt[!is.na(genotype) & genotype != "00" & nchar(genotype) == 2 &
             !is.na(position)]
    setkey(dt, chromosome, position)
    cat("  ", label, ":", format(nrow(dt), big.mark = ","), "autosomal SNPs\n")
    return(dt)
  }

  dt1 <- standardize_plink_dt(dt1, sample1_label)
  dt2 <- standardize_plink_dt(dt2, sample2_label)

  return(list(
    dt1 = list(data = dt1, format = "plink", format_display = "PLINK",
               build = "37", snp_count = nrow(dt1)),
    dt2 = list(data = dt2, format = "plink", format_display = "PLINK",
               build = "37", snp_count = nrow(dt2)),
    sample1_label = sample1_label,
    sample2_label = sample2_label,
    fam = fam
  ))
}

#' Interactive sample picker — presents a searchable/numbered list
#' Shows "FID IID" for each sample. User can select by number, exact
#' "FID IID", exact "IID" (if unique), or partial text match.
interactive_sample_selection <- function(fam) {
  display_ids <- paste(fam$FID, fam$IID)
  n <- length(display_ids)

  cat("\n============================================================\n")
  cat("PLINK Dataset — Sample Selection\n")
  cat("============================================================\n")
  cat("  ", n, "samples available:\n\n")

  # Display as a numbered list (FID + IID)
  col_width <- max(nchar(display_ids)) + 7  # room for index + padding
  cols <- max(1, floor(72 / col_width))
  for (i in seq_along(display_ids)) {
    entry <- sprintf("  %3d) %-*s", i, col_width - 7, display_ids[i])
    cat(entry)
    if (i %% cols == 0 || i == n) cat("\n")
  }

  cat("\n")

  pick_sample <- function(prompt) {
    repeat {
      cat(prompt)
      input <- trimws(readline())

      # Try numeric index
      num <- suppressWarnings(as.integer(input))
      if (!is.na(num) && num >= 1 && num <= n)
        return(display_ids[num])

      # Try exact "FID IID" match
      match_idx <- which(display_ids == input)
      if (length(match_idx) == 1) return(display_ids[match_idx])

      # Try bare IID match (if unique)
      iid_idx <- which(fam$IID == input)
      if (length(iid_idx) == 1) return(display_ids[iid_idx])

      # Try case-insensitive partial match
      match_idx <- grep(input, display_ids, ignore.case = TRUE)
      if (length(match_idx) == 1) return(display_ids[match_idx])
      if (length(match_idx) > 1) {
        cat("  Multiple matches:\n")
        for (j in match_idx)
          cat("    ", j, ") ", display_ids[j], "\n", sep = "")
        cat("  Please be more specific or enter a number.\n")
        next
      }

      cat("  Not found. Enter a number (1-", n,
          "), 'FID IID', or a search term.\n", sep = "")
    }
  }

  s1 <- pick_sample("  Select sample 1 (number, name, or search): ")
  cat("  -> ", s1, "\n")
  s2 <- pick_sample("  Select sample 2 (number, name, or search): ")
  cat("  -> ", s2, "\n")

  if (s1 == s2) stop("Sample 1 and Sample 2 must be different.", call. = FALSE)
  cat("\n")
  return(list(sample1 = s1, sample2 = s2))
}

# =============================================================================
# GENETIC MAP LOADING
# =============================================================================

load_genetic_map <- function(map_dir = NULL) {

  if (!is.null(map_dir) && dir.exists(map_dir)) {
    cat("  Source:", map_dir, "\n")
    map_list <- list()

    for (chr in 1:22) {
      patterns <- c(
        paste0("plink.chr", chr, ".GRCh37.map"),
        paste0("plink.chr", chr, ".GRCh38.map"),
        paste0("genetic_map_chr", chr, "_b37.txt"),
        paste0("genetic_map_chr", chr, "_b36.txt"),
        paste0("genetic_map_GRCh37_chr", chr, ".txt"),
        paste0("genetic_map_chr", chr, ".txt"),
        paste0("chr", chr, ".txt"),
        paste0("chr", chr, ".map")
      )

      map_file <- NULL
      for (p in patterns) {
        fp <- file.path(map_dir, p)
        if (file.exists(fp)) { map_file <- fp; break }
      }
      if (is.null(map_file)) {
        all_files <- list.files(map_dir, full.names = TRUE)
        chr_pat <- paste0("chr", chr, "([^0-9]|$)")
        candidates <- all_files[grepl(chr_pat, basename(all_files))]
        if (length(candidates) > 0) map_file <- candidates[1]
      }
      if (is.null(map_file)) next

      m <- tryCatch(fread(map_file, showProgress = FALSE, header = FALSE),
                    error = function(e) NULL)
      if (is.null(m) || nrow(m) == 0 || ncol(m) < 3) next

      position_vals <- NULL; cm_vals <- NULL
      nc <- ncol(m)

      if (nc == 4) {
        has_header <- any(grepl("pos|map|rate|cm|chrom", names(m), ignore.case = TRUE))
        if (has_header) {
          pos_col <- grep("pos|bp", names(m), ignore.case = TRUE, value = TRUE)
          cm_col  <- grep("^map|^cm", names(m), ignore.case = TRUE, value = TRUE)
          if (length(pos_col) > 0 && length(cm_col) > 0) {
            position_vals <- as.numeric(m[[pos_col[1]]])
            cm_vals       <- as.numeric(m[[cm_col[1]]])
          }
        }
        if (is.null(position_vals)) {
          col4 <- suppressWarnings(as.numeric(head(m[[4]], 20)))
          col4 <- col4[!is.na(col4)]
          col2 <- suppressWarnings(as.numeric(head(m[[2]], 20)))
          col2 <- col2[!is.na(col2)]
          col4_bp <- length(col4) > 0 && max(col4) > 10000
          col2_bp <- length(col2) > 0 && max(col2) > 10000
          if (col4_bp && !col2_bp) {
            position_vals <- as.numeric(m[[4]]); cm_vals <- as.numeric(m[[3]])
          } else if (col2_bp && !col4_bp) {
            position_vals <- as.numeric(m[[2]]); cm_vals <- as.numeric(m[[4]])
          } else if (col4_bp) {
            position_vals <- as.numeric(m[[4]]); cm_vals <- as.numeric(m[[3]])
          }
        }
      } else if (nc == 3) {
        position_vals <- as.numeric(m[[1]]); cm_vals <- as.numeric(m[[3]])
      } else if (nc >= 5) {
        pos_col <- grep("pos|bp", names(m), ignore.case = TRUE, value = TRUE)
        cm_col  <- grep("^map|^cm", names(m), ignore.case = TRUE, value = TRUE)
        if (length(pos_col) > 0 && length(cm_col) > 0) {
          position_vals <- as.numeric(m[[pos_col[1]]])
          cm_vals       <- as.numeric(m[[cm_col[1]]])
        }
      }

      if (!is.null(position_vals) && !is.null(cm_vals)) {
        map_list[[chr]] <- data.table(
          chromosome = chr, position = as.integer(position_vals), cM = cm_vals
        )
      }
    }

    if (length(map_list) > 0) {
      gmap <- rbindlist(map_list)
      gmap <- gmap[is.finite(position) & is.finite(cM) & position > 0]
      gmap <- unique(gmap, by = c("chromosome", "position"))
      setkey(gmap, chromosome, position)
      cat("  Loaded", format(nrow(gmap), big.mark = ","),
          "map points across", length(unique(gmap$chromosome)), "chromosomes\n")
      return(gmap)
    }
  }

  # Auto-download fallback
  cat("  Downloading Beagle PLINK GRCh37 genetic map...\n")
  map_url <- "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip"
  map_cache <- file.path(tempdir(), "genetic_map_grch37")
  if (!dir.exists(map_cache)) {
    dir.create(map_cache, recursive = TRUE)
    zip_file <- file.path(tempdir(), "genetic_map.zip")
    tryCatch({
      download.file(map_url, zip_file, quiet = TRUE, mode = "wb")
      unzip(zip_file, exdir = map_cache)
    }, error = function(e) {
      cat("  WARNING: Download failed. Using 1 Mb ~ 1 cM approximation.\n")
      return(NULL)
    })
  }
  if (dir.exists(map_cache) && length(list.files(map_cache)) > 0)
    return(load_genetic_map(map_cache))
  return(NULL)
}

# =============================================================================
# cM INTERPOLATION
# =============================================================================

interpolate_cM <- function(positions, chr, genetic_map) {
  if (is.null(genetic_map)) return(positions / 1e6)
  map_chr <- genetic_map[chromosome == chr]
  map_chr <- map_chr[is.finite(position) & is.finite(cM)]
  map_chr <- unique(map_chr, by = "position")
  setorder(map_chr, position)
  if (nrow(map_chr) < 2) return(positions / 1e6)
  result <- tryCatch(
    approx(x = map_chr$position, y = map_chr$cM, xout = positions, rule = 2)$y,
    error = function(e) positions / 1e6
  )
  na_idx <- is.na(result)
  if (any(na_idx)) result[na_idx] <- positions[na_idx] / 1e6
  return(result)
}

# =============================================================================
# HALF-IDENTICAL MATCH TEST
# =============================================================================

is_half_identical <- function(geno1, geno2) {
  a1 <- substr(geno1, 1, 1); b1 <- substr(geno1, 2, 2)
  a2 <- substr(geno2, 1, 1); b2 <- substr(geno2, 2, 2)
  return(a1 == a2 | a1 == b2 | b1 == a2 | b1 == b2)
}

# =============================================================================
# IBD SEGMENT DETECTION
# =============================================================================
# Seed-and-extend algorithm with sliding windows:
#   1. Slide a window of `mismatch_bunch` SNPs across the chromosome
#   2. Seed IBD where window mismatch rate <= seed_threshold (default 2%)
#   3. Extend seeds into adjacent windows below an adaptive threshold
#      (half the observed background mismatch rate, capped at 4%)
#   4. Trim segment edges to seed-quality window boundaries
#   5. Validate: overall mismatch rate, seed density, size thresholds
#   6. Stricter filtering for small segments (< 15 cM)

find_shared_segments <- function(dt1, dt2, genetic_map,
                                  min_cm, min_snps, mismatch_bunch,
                                  seed_threshold) {
  cat("\nComparing genotypes...\n")
  merged <- merge(dt1, dt2, by = c("rsid", "chromosome", "position"),
                  suffixes = c("_1", "_2"))
  setkey(merged, chromosome, position)
  total_overlapping <- nrow(merged)
  cat("  Overlapping SNPs:", format(total_overlapping, big.mark = ","), "\n")

  merged[, match := is_half_identical(genotype_1, genotype_2)]
  all_segments <- list()

  for (chr in sort(unique(merged$chromosome))) {
    chr_data <- merged[chromosome == chr]
    n_snps <- nrow(chr_data)
    if (n_snps < min_snps) next

    n_matches <- sum(chr_data$match)
    cat("  Chr", sprintf("%2d", chr), ":", format(n_snps, big.mark = ","),
        "SNPs |", format(n_matches, big.mark = ","), "matches",
        paste0("(", round(100 * n_matches / n_snps, 1), "%)"), "\n")

    chr_data[, cM := interpolate_cM(position, chr, genetic_map)]

    segs <- detect_segments_on_chr(
      chr_data$match, chr_data$position, chr_data$cM,
      n_snps, min_cm, min_snps, mismatch_bunch, seed_threshold, chr
    )
    if (length(segs) > 0) all_segments <- c(all_segments, segs)
  }

  if (length(all_segments) == 0) {
    cat("\nNo shared segments found above threshold.\n")
    return(list(
      segments = data.table(Chr = integer(), Start_Position = integer(),
                            End_Position = integer(), cM = numeric(),
                            SNPs = integer()),
      total_overlapping = total_overlapping
    ))
  }

  result <- rbindlist(all_segments)
  setorder(result, Chr, Start_Position)
  return(list(segments = result, total_overlapping = total_overlapping))
}

detect_segments_on_chr <- function(matches, positions, cm_positions, n_snps,
                                    min_cm, min_snps, window_size,
                                    seed_threshold, chr) {
  segments <- list()
  if (n_snps < window_size) return(segments)

  mismatch_vec <- as.integer(!matches)
  total_mismatch_rate <- sum(mismatch_vec) / n_snps

  # Sliding window mismatch rates via cumulative sum
  cum <- c(0L, cumsum(mismatch_vec))
  n_windows <- n_snps - window_size + 1
  window_mismatches <- cum[(window_size + 1):(n_snps + 1)] - cum[1:n_windows]
  window_rates <- window_mismatches / window_size

  # Thresholds
  extend_threshold <- min(total_mismatch_rate * 0.50, 0.04)
  is_seed <- window_rates <= seed_threshold
  is_extendable <- window_rates <= extend_threshold

  # Seed and extend (bidirectional)
  ibd_window <- logical(n_windows)
  ibd_window[is_seed] <- TRUE
  for (w in 2:n_windows)
    if (!ibd_window[w] && ibd_window[w-1] && is_extendable[w]) ibd_window[w] <- TRUE
  for (w in (n_windows - 1):1)
    if (!ibd_window[w] && ibd_window[w+1] && is_extendable[w]) ibd_window[w] <- TRUE

  if (!any(ibd_window)) return(segments)

  # Window runs -> SNP-level segments, trimmed to seed boundaries
  w_runs <- rle(ibd_window)
  w_ends <- cumsum(w_runs$lengths)
  w_starts <- c(1L, w_ends[-length(w_ends)] + 1L)

  ibd_snp <- logical(n_snps)
  for (r in seq_along(w_runs$lengths)) {
    if (!w_runs$values[r]) next
    seed_in_run <- which(is_seed[w_starts[r]:w_ends[r]]) + w_starts[r] - 1L
    if (length(seed_in_run) == 0) next
    s <- seed_in_run[1]
    e <- min(seed_in_run[length(seed_in_run)] + window_size - 1, n_snps)
    ibd_snp[s:e] <- TRUE
  }

  # Extract and validate segments
  runs <- rle(ibd_snp)
  if (!any(runs$values)) return(segments)
  run_ends <- cumsum(runs$lengths)
  run_starts <- c(1L, run_ends[-length(run_ends)] + 1L)

  for (r in seq_along(runs$lengths)) {
    if (!runs$values[r]) next
    s <- run_starts[r]; e <- run_ends[r]
    n <- e - s + 1
    cm <- cm_positions[e] - cm_positions[s]
    mm_rate <- sum(mismatch_vec[s:e]) / n

    if (mm_rate > 0.025) next                    # Overall mismatch gate
    if (cm < 15 && mm_rate > 0.015) next         # Stricter for small segments

    # Seed density check
    ws <- s; we <- min(e - window_size + 1, n_windows)
    if (we >= ws && sum(is_seed[ws:we]) / (we - ws + 1) < 0.30) next

    if (cm >= min_cm && n >= min_snps) {
      segments[[length(segments) + 1]] <- data.table(
        Chr = chr, Start_Position = positions[s],
        End_Position = positions[e], cM = round(cm, 2), SNPs = n
      )
    }
  }
  return(segments)
}

# =============================================================================
# CHROMOSOME VISUALIZATION — Modern Dark Theme
# =============================================================================

plot_chromosome_comparison <- function(result, plot_file, file_info = NULL) {
  segs <- result$segments
  if (nrow(segs) == 0) { cat("  No segments to plot.\n"); return(invisible(NULL)) }

  chr_lengths <- c(
    249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
    159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
    133275309, 107349540, 102531392, 90354753, 81195210, 78077248,
    59128983, 63025520, 48129895, 51304566
  )
  names(chr_lengths) <- 1:22
  max_len <- max(chr_lengths); n_chr <- 22

  # --- Color palette ---
  bg_color      <- "#1a1a2e"    # deep navy background
  panel_color   <- "#16213e"    # slightly lighter panel
  chr_fill      <- "#2a2a4a"    # chromosome bar fill
  chr_border    <- "#3a3a5c"    # chromosome bar border
  seg_fill      <- "#00d4aa"    # teal/mint for segments
  seg_border    <- "#00b894"    # slightly darker teal border
  seg_glow      <- "#00d4aa40"  # translucent glow behind segments
  text_primary  <- "#e8e8f0"    # bright text
  text_secondary <- "#8888aa"   # dimmer text
  accent        <- "#6c5ce7"    # purple accent
  grid_color    <- "#2a2a4a"    # subtle grid lines

  png(plot_file, width = 1600, height = 900, res = 110, bg = bg_color)

  par(mar = c(4.5, 4.5, 4, 19), xpd = FALSE, bg = bg_color,
      fg = text_primary, col.axis = text_secondary, col.lab = text_secondary,
      col.main = text_primary, family = "sans")

  plot(NULL, xlim = c(0, max_len), ylim = c(0.5, n_chr + 0.5),
       xlab = "", ylab = "", yaxt = "n", xaxt = "n",
       main = "", bty = "n")

  # Title with accent line
  mtext("ChromoMatch", side = 3, line = 2.2, cex = 1.6, font = 2,
        col = text_primary, adj = 0)
  mtext("Shared IBD Segments", side = 3, line = 0.8, cex = 1.0,
        col = text_secondary, adj = 0)

  # Subtle horizontal grid lines
  for (i in 1:n_chr) {
    y <- n_chr - i + 1
    abline(h = y, col = grid_color, lwd = 0.5)
  }

  # X-axis
  tick_at <- axTicks(1)
  axis(1, at = tick_at,
       labels = paste0(round(tick_at / 1e6), "M"),
       col = grid_color, col.ticks = text_secondary,
       cex.axis = 0.8, lwd = 0.5)
  mtext("Position (bp)", side = 1, line = 2.8, cex = 0.9, col = text_secondary)

  # Y-axis
  axis(2, at = 1:n_chr, labels = rev(1:n_chr), las = 2,
       col = grid_color, col.ticks = text_secondary,
       cex.axis = 0.85, lwd = 0.5)
  mtext("Chromosome", side = 2, line = 2.8, cex = 0.9, col = text_secondary)

  # Draw chromosome bars with rounded appearance
  for (i in 1:n_chr) {
    y <- n_chr - i + 1
    rect(0, y - 0.28, chr_lengths[i], y + 0.28,
         col = chr_fill, border = chr_border, lwd = 0.8)
  }

  # Draw segments: glow layer first, then solid
  for (i in 1:nrow(segs)) {
    y <- n_chr - segs$Chr[i] + 1
    # Glow
    rect(segs$Start_Position[i] - max_len * 0.002, y - 0.38,
         segs$End_Position[i] + max_len * 0.002, y + 0.38,
         col = seg_glow, border = NA)
    # Solid segment
    rect(segs$Start_Position[i], y - 0.28,
         segs$End_Position[i], y + 0.28,
         col = seg_fill, border = seg_border, lwd = 0.8)
  }

  # ---- Right-side summary panel ----
  par(xpd = TRUE)
  x <- max_len * 1.05
  yt <- n_chr

  # Summary stats
  text(x, yt, paste("Total cM:", round(sum(segs$cM), 2)),
       adj = 0, cex = 1.1, font = 2, col = seg_fill)
  text(x, yt - 1, paste("Overlapping SNPs:",
       format(result$total_overlapping, big.mark = ",")),
       adj = 0, cex = 0.9, col = text_primary)
  text(x, yt - 2, paste("Segments:", nrow(segs)),
       adj = 0, cex = 0.9, col = text_primary)

  # Estimated relationship
  y_after_stats <- yt - 2
  if (!is.null(file_info) && !is.null(file_info$relationship)) {
    y_after_stats <- yt - 3.2
    text(x, yt - 3.2, file_info$relationship,
         adj = 0, cex = 0.8, font = 3, col = "#f8b500")
  }

  # Accent divider
  segments(x, y_after_stats - 0.7, x + max_len * 0.18, y_after_stats - 0.7,
           col = accent, lwd = 2)

  # Dataset name (for PLINK mode)
  y_cursor <- y_after_stats - 1.0
  if (!is.null(file_info) && !is.null(file_info$dataset_name)) {
    text(x, y_cursor, paste("Dataset:", file_info$dataset_name),
         adj = 0, cex = 0.7, family = "mono", col = text_secondary)
    y_cursor <- y_cursor - 0.7
  }

  # File/sample information
  if (!is.null(file_info)) {
    y_cursor <- y_cursor - 0.1
    text(x, y_cursor, file_info$name1,
         adj = 0, cex = 0.65, family = "mono", col = seg_fill)
    text(x, y_cursor - 0.6, paste0("  ", format(file_info$snps1, big.mark = ","),
                              " SNPs (", file_info$fmt1, ")"),
         adj = 0, cex = 0.6, family = "mono", col = text_secondary)
    y_cursor <- y_cursor - 1.4
    text(x, y_cursor, file_info$name2,
         adj = 0, cex = 0.65, family = "mono", col = seg_fill)
    text(x, y_cursor - 0.6, paste0("  ", format(file_info$snps2, big.mark = ","),
                              " SNPs (", file_info$fmt2, ")"),
         adj = 0, cex = 0.6, family = "mono", col = text_secondary)
    y_cursor <- y_cursor - 1.4
  }

  # Segment table
  y_tbl <- y_cursor - 0.4
  segments(x, y_tbl + 0.3, x + max_len * 0.18, y_tbl + 0.3,
           col = accent, lwd = 1)
  text(x, y_tbl,
       sprintf("%-4s %12s %12s %6s %6s", "Chr", "Start", "End", "cM", "SNPs"),
       adj = 0, cex = 0.65, family = "mono", col = text_secondary)

  max_rows <- min(nrow(segs), 14)
  for (i in 1:max_rows) {
    row_col <- if (i %% 2 == 0) text_secondary else text_primary
    text(x, y_tbl - i * 0.68,
         sprintf("%-4d %12s %12s %6.2f %6s", segs$Chr[i],
                 format(segs$Start_Position[i], big.mark = ","),
                 format(segs$End_Position[i], big.mark = ","),
                 segs$cM[i], format(segs$SNPs[i], big.mark = ",")),
         adj = 0, cex = 0.6, family = "mono", col = row_col)
  }
  if (nrow(segs) > max_rows)
    text(x, y_tbl - (max_rows + 1) * 0.68,
         paste("... +", nrow(segs) - max_rows, "more"),
         adj = 0, cex = 0.65, font = 3, col = text_secondary)

  dev.off()
  cat("  Plot saved to:", plot_file, "\n")
}

# =============================================================================
# RELATIONSHIP ESTIMATION (Blaine Bettinger's Shared cM Project ranges)
# =============================================================================

estimate_relationship <- function(total_cm, silent = FALSE) {
  rel <- if (total_cm > 3400) "Parent/Child or Identical Twin (~3400 cM)"
    else if (total_cm > 2300) "Full Sibling (~2550 cM)"
    else if (total_cm > 1700) "Grandparent, Aunt/Uncle, Half-Sibling (~1750 cM)"
    else if (total_cm > 1100) "1st Cousin, Great-Grandparent (~1200 cM)"
    else if (total_cm > 500)  "1st Cousin 1x removed, Half-1C (~600-900 cM)"
    else if (total_cm > 200)  "2nd Cousin or 2C1R (~250 cM)"
    else if (total_cm > 90)   "3rd Cousin (~100 cM)"
    else if (total_cm > 20)   "4th-6th Cousin"
    else                      "Very distant or no detectable relationship"
  if (!silent) cat("  Estimated relationship:", rel, "\n")
  return(invisible(rel))
}

# =============================================================================
# RESULTS DISPLAY
# =============================================================================

print_results <- function(result) {
  segs <- result$segments
  cat("\n============================================================\n")
  cat("ChromoMatch RESULTS\n")
  cat("============================================================\n")
  cat("  Total shared cM: ", round(sum(segs$cM), 2), "\n")
  cat("  SNPs overlapping:", format(result$total_overlapping, big.mark = ","), "\n")
  cat("  Segments:        ", nrow(segs), "\n")
  cat("------------------------------------------------------------\n")
  cat(sprintf("  %-4s %14s %14s %8s %7s\n", "Chr", "Start", "End", "cM", "SNPs"))
  cat("  ", paste(rep("-", 51), collapse = ""), "\n")
  for (i in 1:nrow(segs)) {
    cat(sprintf("  %-4d %14s %14s %8.2f %7s\n",
                segs$Chr[i],
                format(segs$Start_Position[i], big.mark = ","),
                format(segs$End_Position[i], big.mark = ","),
                segs$cM[i],
                format(segs$SNPs[i], big.mark = ",")))
  }
  cat("============================================================\n")
  estimate_relationship(sum(segs$cM))
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  cat("============================================================\n")
  cat("ChromoMatch\n")
  cat("============================================================\n")

  if (input_mode == "plink") {
    # ------ PLINK binary input ------
    cat("  Mode:  PLINK binary dataset\n")
    cat("  Bfile:", plink_prefix, "\n")

    # Resolve sample IDs (interactive if not provided)
    s1 <- plink_sample1
    s2 <- plink_sample2
    if (is.null(s1) || is.null(s2)) {
      fam <- read_plink_fam(paste0(plink_prefix, ".fam"))
      picks <- interactive_sample_selection(fam)
      s1 <- picks$sample1
      s2 <- picks$sample2
    }

    cat("  Sample 1:", s1, "\n")
    cat("  Sample 2:", s2, "\n")
    cat("  Min cM:", defaults$min_cm, "| Min SNPs:", defaults$min_snps, "\n")
    cat("  Window:", defaults$mismatch_bunch,
        "| Seed threshold:", defaults$seed_threshold, "\n")
    cat("------------------------------------------------------------\n")

    cat("\nReading PLINK dataset...\n")
    plink_data <- read_plink_bed(plink_prefix, s1, s2)
    raw1 <- plink_data$dt1
    raw2 <- plink_data$dt2

    # Labels for plot (resolved FID+IID from the fam file)
    label1 <- plink_data$sample1_label
    label2 <- plink_data$sample2_label

  } else {
    # ------ Raw DNA file input ------
    cat("  Mode:  Raw DNA files\n")
    cat("  File 1:", file1, "\n")
    cat("  File 2:", file2, "\n")
    cat("  Min cM:", defaults$min_cm, "| Min SNPs:", defaults$min_snps, "\n")
    cat("  Window:", defaults$mismatch_bunch,
        "| Seed threshold:", defaults$seed_threshold, "\n")
    cat("------------------------------------------------------------\n")

    cat("\nReading files...\n")
    raw1 <- detect_and_read_raw(file1)
    raw2 <- detect_and_read_raw(file2)

    if (raw1$build != raw2$build)
      cat("\n  WARNING: Different builds (", raw1$build, "vs", raw2$build, ")\n")

    label1 <- basename(file1)
    label2 <- basename(file2)
  }

  cat("\nLoading genetic map...\n")
  genetic_map <- load_genetic_map(defaults$genetic_map_dir)

  result <- find_shared_segments(
    raw1$data, raw2$data, genetic_map,
    min_cm = defaults$min_cm, min_snps = defaults$min_snps,
    mismatch_bunch = defaults$mismatch_bunch,
    seed_threshold = defaults$seed_threshold
  )

  if (nrow(result$segments) > 0) {
    print_results(result)
    fwrite(result$segments, defaults$output_file)
    cat("\n  Results saved to:", defaults$output_file, "\n")
    if (defaults$plot) {
      cat("  Generating plot...\n")
      rel_str <- estimate_relationship(sum(result$segments$cM), silent = TRUE)
      file_info <- list(
        name1 = label1,
        snps1 = raw1$snp_count,
        fmt1  = raw1$format_display,
        name2 = label2,
        snps2 = raw2$snp_count,
        fmt2  = raw2$format_display,
        dataset_name = if (input_mode == "plink") basename(plink_prefix) else NULL,
        relationship = rel_str
      )
      plot_chromosome_comparison(result, defaults$plot_file, file_info)
    }
  }
  cat("\nDone.\n")
}

main()
