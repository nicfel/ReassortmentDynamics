#!/usr/bin/env Rscript
# =============================================================================
# buildCoalReXml.R
#
# Flexible builder for coalRe (BEAST 2) reassortment-network XMLs.
#
# Point it at a directory of per-segment FASTA files for a segmented virus and
# it will generate a ready-to-run analysis directory containing the XML(s) and
# LSF submission script(s), laid out like the example `dependent_Ne/` and
# `independent_Ne/` folders (one sub-folder per replicate, each holding an
# `.xml` and a `sub.sh`).
#
# It adapts automatically to the number of segments found in the FASTA
# directory (2, 3, 8, 10, ... whatever is there) by expanding / collapsing the
# per-segment blocks of the template.
#
# REQUIREMENTS on the input FASTA files (one per segment):
#   * All segments must share the SAME taxon headers (the intersection of taxa
#     is used; a warning is printed if they differ).
#   * Each file name must end in `_<LABEL>.fasta`; the LABEL is the string after
#     the LAST underscore (e.g. `HPAI_all_North_America_HA.fasta` -> "HA").
#   * Headers are `name|accession|DATE|...` style, pipe-separated, with the
#     sampling date in the 3rd field (configurable below).
#
# Only base R is used (no external packages required).
# =============================================================================

rm(list = ls())

# ============================== CONFIGURATION ================================
# --- Reassortment model: "dependent", "independent", or "both" ---------------
#   "dependent"   -> reassortment rate is a function of Ne   (isNeSkyline block)
#   "independent" -> reassortment rate is independent of Ne  (isISkyline block)
MODE <- "both"

# --- Input / output ----------------------------------------------------------
FASTA_DIR   <- "example_fasta"                 # dir with per-segment FASTA files
TEMPLATE    <- "inference_template_wgs_cr.xml"  # coalRe template (8-segment)
OUTPUT_ROOT <- "coalre_output"                  # dir that will hold <mode>_Ne/
NAME        <- "coalre"                         # analysis name (used in filenames)
N_REPS      <- 3                                 # number of replicate folders

# --- Optional explicit segment ordering (seg1, seg2, ...) --------------------
# Give the segment LABELS in the desired order, e.g. c("HA","NA","MP").
# Leave as NULL to order the segments alphabetically by label.
SEGMENT_ORDER <- NULL

# --- Molecular clock ---------------------------------------------------------
CLOCK_RATE <- 0.0035                             # substitutions/site/year

# --- Date parsing ------------------------------------------------------------
DATE_SEP    <- "\\|"    # regex separator used to split the taxon header
DATE_FIELD  <- 3        # which field (after splitting) holds the date
DATE_FORMAT <- "%Y-%m-%d"

# --- Rate-shift grid (Ne / reassortment-rate change points, in years) --------
# A fine grid over the sampled period plus a few coarse deep-time anchors.
N_FINE   <- 15      # # of evenly spaced shifts over the sampling window
N_COARSE <- 4       # # of coarse shifts beyond the sampling window
DEEP_MAX <- 30      # oldest coarse anchor (years before most recent sample)
# reassortment/Ne coupling is estimated independently only within the fine grid
INDEPENDENT_AFTER <- NULL   # NULL -> use N_FINE; or set an explicit integer

# --- Submission script -------------------------------------------------------
SCHEDULER    <- "LSF"                    # "LSF", "SLURM", or "PBS"
SUB_WALLTIME <- "720:00"                 # walltime as HH:MM (":00" appended for SLURM/PBS)
SUB_MEM      <- "10240"                  # memory (MB)
SUB_EMAIL    <- "user@example.com"       # notification email
SUB_MODULES  <- c("beast2", "beagle/4.0.0", "gcc/8.5.0")
SUB_BEAST    <- "beast -beagle -beagle_SSE"

# ============================================================================
#                          (implementation below)
# ============================================================================

# ---- tidy relative paths so the script can be run from anywhere -------------
FASTA_DIR   <- normalizePath(FASTA_DIR,  mustWork = TRUE)
TEMPLATE    <- normalizePath(TEMPLATE,   mustWork = TRUE)
dir.create(OUTPUT_ROOT, showWarnings = FALSE, recursive = TRUE)
OUTPUT_ROOT <- normalizePath(OUTPUT_ROOT, mustWork = TRUE)

# ---- small helpers ----------------------------------------------------------

# Remove quote-like characters and parentheses from a taxon header.
quote_chars <- "['‘’`´ʼʹʽˈˊ′‵‶‷]"
clean_name <- function(x) {
  x <- gsub(quote_chars, "", x, perl = TRUE)
  x <- gsub("[()]", "", x)
  trimws(x)
}

# Minimal FASTA reader -> named character vector (name = cleaned header).
read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines))]
  hdr   <- grepl("^>", lines)
  ids   <- clean_name(sub("^>", "", lines[hdr]))
  grp   <- cumsum(hdr)
  seqs  <- tapply(lines[!hdr], grp[!hdr], function(s) paste0(s, collapse = ""))
  seqs  <- seqs[order(as.integer(names(seqs)))]
  stats::setNames(as.character(seqs), ids)
}

write_fasta <- function(seqs, path) {
  con <- file(path, "w")
  on.exit(close(con))
  for (nm in names(seqs)) {
    writeLines(paste0(">", nm), con)
    writeLines(seqs[[nm]], con)
  }
}

# The segment label = string after the last underscore, before .fasta
label_from_file <- function(f) {
  base <- sub("\\.fasta$", "", basename(f), ignore.case = TRUE)
  sub("^.*_", "", base)
}

# -----------------------------------------------------------------------------
# Generic per-segment block expander.
#
# The template is written for `template_seg` segments (auto-detected, e.g. 8).
# Every per-segment line carries the token `seg<k>`, and structurally identical
# sub-blocks appear once per segment, in order (seg1, seg2, ..., seg8). This
# walks the template and, at every seg1 sub-block that begins a full
# seg1..seg<template_seg> run, re-emits the run for `n_seg` segments.
#
# It is tolerant of two real-world irregularities in the template:
#   * cosmetic blank lines that differ between sub-blocks (ignored for matching,
#     one blank re-inserted between emitted blocks);
#   * the "define once" idiom where the seg1 block DEFINES a shared object
#     (e.g. StrictClock.c) that later blocks only reference. Segment 1 is kept
#     verbatim; segment 2's block is used as the stencil for segments >= 2.
#
# Lone seg1 anchors that are not followed by a matching seg2 run, and anything
# inside an XML comment, are copied through unchanged.
# -----------------------------------------------------------------------------

seg_indices <- function(line) {
  toks <- regmatches(line, gregexpr("seg([0-9]+)", line))[[1]]
  if (!length(toks)) return(integer(0))
  as.integer(sub("seg", "", toks))
}
is_blank <- function(line) !nzchar(trimws(line))

expand_segments <- function(lines, n_seg) {
  n <- length(lines)
  template_seg <- max(c(0L, unlist(lapply(lines, seg_indices))))
  if (template_seg < 2L) return(lines)   # nothing segmented to expand

  out <- character(0)
  i <- 1L
  in_comment <- FALSE

  while (i <= n) {
    line <- lines[i]
    seg_i <- seg_indices(line)

    # A section can only start on a live (non-comment) line whose only segment
    # reference is seg1 and which does not itself open a comment.
    can_start <- !in_comment && !grepl("<!--", line, fixed = TRUE) &&
                 length(seg_i) >= 1L && all(seg_i == 1L)

    if (can_start) {
      # first following line that references seg2 -> start of block 2
      seg2_first <- NA_integer_
      j <- i + 1L
      while (j <= n) {
        sj <- seg_indices(lines[j])
        if (length(sj) && any(sj == 2L)) { seg2_first <- j; break }
        j <- j + 1L
      }

      section_ok <- FALSE
      if (!is.na(seg2_first) && seg2_first > i) {
        block1 <- lines[i:(seg2_first - 1L)]
        b1segs <- unlist(lapply(block1, seg_indices))
        # block 1 must reference only seg1
        if (length(b1segs) == 0L || all(b1segs == 1L)) {
          nb1 <- sum(!vapply(block1, is_blank, logical(1)))
          if (nb1 >= 1L) {
            # collect template_seg * nb1 non-blank lines starting at i,
            # remembering the source line index of each
            need <- template_seg * nb1
            idx <- integer(0)
            k <- i
            while (k <= n && length(idx) < need) {
              if (!is_blank(lines[k])) idx <- c(idx, k)
              k <- k + 1L
            }
            if (length(idx) == need) {
              groups <- lapply(1:template_seg, function(m)
                lines[idx[((m - 1L) * nb1 + 1L):(m * nb1)]])
              # each group m must reference only segment m
              section_ok <- all(vapply(1:template_seg, function(m) {
                gs <- unlist(lapply(groups[[m]], seg_indices))
                length(gs) == 0L || all(gs == m)
              }, logical(1)))
              last_consumed <- idx[need]
            }
          }
        }
      }

      if (section_ok) {
        for (kk in seq_len(n_seg)) {
          src   <- if (kk == 1L) 1L else 2L
          block <- gsub(paste0("seg", src), paste0("seg", kk), groups[[src]])
          out <- c(out, block)
          # readability: blank line only between multi-line blocks, not list items
          if (kk < n_seg && nb1 > 1L) out <- c(out, "")
        }
        i <- last_consumed + 1L
        next
      }
    }

    # verbatim copy + multi-line comment tracking
    out <- c(out, line)
    opens  <- length(regmatches(line, gregexpr("<!--", line, fixed = TRUE))[[1]])
    closes <- length(regmatches(line, gregexpr("-->",  line, fixed = TRUE))[[1]])
    if (opens > closes) in_comment <- TRUE
    else if (closes > 0) in_comment <- FALSE
    i <- i + 1L
  }
  out
}

# -----------------------------------------------------------------------------
# Fill scalar placeholders + apply the dependent/independent toggle.
# `mode` is "dependent" or "independent".
# -----------------------------------------------------------------------------
fill_template <- function(lines, mode, taxa, dates_chr, rateshifts,
                          indep_after, seq_lengths, labels, clock_rate,
                          n_seg) {

  # markers to comment-out (disable). Exactly one skyline block stays live.
  disable <- c("isconstant")
  if (mode == "dependent")   disable <- c(disable, "isISkyline")
  if (mode == "independent") disable <- c(disable, "isNeSkyline")

  tips_lines <- sprintf('\t\t<taxon spec="Taxon" id="%s"/>', taxa)
  heights_val <- paste(paste0(taxa, "=", dates_chr), collapse = ",")
  times_val   <- paste(rateshifts, collapse = " ")
  weights_val <- paste(seq_lengths, collapse = " ")

  out <- character(0)
  i <- 1L
  n <- length(lines)
  while (i <= n) {
    line <- lines[i]

    # --- disable non-selected model blocks (comment them out) ---------------
    for (mk in disable) {
      if (grepl(paste0(mk, "-->"), line, fixed = TRUE))
        line <- gsub(paste0(mk, "-->"), "", line, fixed = TRUE)
      if (grepl(paste0("<!--", mk), line, fixed = TRUE))
        line <- gsub(paste0("<!--", mk), "", line, fixed = TRUE)
    }

    # --- independent-of-Ne reassortment-rate prior swap ---------------------
    if (mode == "independent" &&
        grepl('spec="Prior" x="@reassortmentRate"', line, fixed = TRUE)) {
      out <- c(out,
        '                    <distribution spec="Prior">',
        '                        <x spec="coalre.dynamics.LogDifference" arg="@reassortmentRate"/>',
        '                        <distr spec="beast.base.inference.distribution.Normal" mean="0" sigma="1.0"/>')
      i <- i + 2L   # skip original opening line + the LogNormal <distr> line
      next
    }

    # --- multi-line list placeholders ---------------------------------------
    if (grepl("insert_tips", line)) {
      out <- c(out, tips_lines); i <- i + 1L; next
    }
    if (grepl("insert_heights", line)) {
      out <- c(out, gsub("insert_heights", heights_val, line)); i <- i + 1L; next
    }

    # --- scalar substitutions -----------------------------------------------
    line <- gsub("insert_times", times_val, line)
    line <- gsub("insert_clock_rate", format(clock_rate, scientific = FALSE), line)
    line <- gsub("insert_independent_after", indep_after, line)
    line <- gsub("insert_from_to_prior", "", line)

    # number of segments in the coalescent network
    line <- gsub('nSegments="[0-9]+"', sprintf('nSegments="%d"', n_seg), line)

    # per-segment fasta filenames -> ../<label>.fasta (relative to repX/)
    if (grepl("insert_seg[0-9]+\\.fasta", line)) {
      for (k in seq_len(n_seg)) {
        line <- gsub(sprintf("insert_seg%d\\.fasta", k),
                     sprintf("../%s.fasta", labels[k]), line)
      }
    }

    # mutation-rate delta-exchange weights = per-segment alignment lengths
    if (grepl('id="weightparameter"', line)) {
      line <- gsub('dimension="[0-9]+"', sprintf('dimension="%d"', n_seg), line)
      line <- sub(">[0-9 ]+</weightvector>",
                  sprintf(">%s</weightvector>", weights_val), line)
    }

    out <- c(out, line)
    i <- i + 1L
  }
  out
}

# -----------------------------------------------------------------------------
# Submission script (LSF / SLURM / PBS)
# -----------------------------------------------------------------------------

# SLURM/PBS want HH:MM:SS; if the user gave HH:MM (LSF style) append ":00".
walltime_hms <- function(w) if (grepl("^[0-9]+:[0-9]+$", w)) paste0(w, ":00") else w

make_sub_sh <- function(xml_file, job_name) {
  sched     <- toupper(SCHEDULER)
  modules   <- paste0("module load ", SUB_MODULES)
  run_line  <- sprintf("%s %s", SUB_BEAST, xml_file)

  if (sched == "LSF") {
    header <- c(
      "#!/bin/bash",
      sprintf("#BSUB -W %s # walltime", SUB_WALLTIME),
      sprintf("#BSUB -J %s # job name", job_name),
      sprintf("#BSUB -o %s_run.out", job_name),
      sprintf("#BSUB -e %s_run.err", job_name),
      sprintf("#BSUB -M %s", SUB_MEM),
      "#BSUB -N",
      sprintf("#BSUB -u %s    # sends email upon job completion", SUB_EMAIL))
    body <- c("", "", modules, "", run_line)

  } else if (sched == "SLURM") {
    header <- c(
      "#!/bin/bash",
      sprintf("#SBATCH --job-name=%s", job_name),
      sprintf("#SBATCH --output=%s_run.out", job_name),
      sprintf("#SBATCH --error=%s_run.err", job_name),
      sprintf("#SBATCH --time=%s", walltime_hms(SUB_WALLTIME)),
      sprintf("#SBATCH --mem=%s", SUB_MEM),
      "#SBATCH --nodes=1",
      "#SBATCH --ntasks=1",
      "#SBATCH --cpus-per-task=1",
      "#SBATCH --mail-type=END,FAIL",
      sprintf("#SBATCH --mail-user=%s", SUB_EMAIL))
    body <- c("", "", modules, "", run_line)

  } else if (sched == "PBS") {
    header <- c(
      "#!/bin/bash",
      sprintf("#PBS -N %s", job_name),
      sprintf("#PBS -o %s_run.out", job_name),
      sprintf("#PBS -e %s_run.err", job_name),
      sprintf("#PBS -l walltime=%s", walltime_hms(SUB_WALLTIME)),
      "#PBS -l nodes=1:ppn=1",
      sprintf("#PBS -l mem=%smb", SUB_MEM),
      "#PBS -m ae",
      sprintf("#PBS -M %s", SUB_EMAIL))
    # PBS starts the job in $HOME; move to the submission directory first.
    body <- c("", "cd $PBS_O_WORKDIR", "", modules, "", run_line)

  } else {
    stop("Unknown SCHEDULER '", SCHEDULER, "'. Use \"LSF\", \"SLURM\", or \"PBS\".")
  }

  c(header, body)
}

# =============================== DRIVER ======================================

build_one_mode <- function(mode, template_lines, taxa, dates_chr, rateshifts,
                           indep_after, seq_lengths, labels, n_seg) {

  message(sprintf("  building '%s' XML ...", mode))
  expanded <- expand_segments(template_lines, n_seg)
  filled   <- fill_template(expanded, mode, taxa, dates_chr, rateshifts,
                            indep_after, seq_lengths, labels, CLOCK_RATE, n_seg)

  # sanity checks
  leftover <- grep("insert_", filled, value = TRUE)
  if (length(leftover))
    warning(sprintf("[%s] unfilled placeholders remain:\n%s",
                    mode, paste(leftover, collapse = "\n")))

  out_dir <- file.path(OUTPUT_ROOT, paste0(mode, "_Ne"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # copy cleaned segment fastas once into the analysis dir
  # (done by caller via attribute) -- handled outside

  for (rep in seq_len(N_REPS)) {
    rep_dir <- file.path(out_dir, paste0("rep", rep))
    dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)
    xml_name <- sprintf("%s.%s.rep%d.xml", NAME, mode, rep)
    writeLines(filled, file.path(rep_dir, xml_name))
    job <- sprintf("%s_%s_rep%d", NAME, mode, rep)
    writeLines(make_sub_sh(xml_name, job), file.path(rep_dir, "sub.sh"))
  }
  message(sprintf("  -> %s (%d replicates)", out_dir, N_REPS))
  out_dir
}

main <- function() {
  # ---- discover segments ----------------------------------------------------
  files <- list.files(FASTA_DIR, pattern = "\\.fasta$", full.names = TRUE,
                       ignore.case = TRUE)
  if (length(files) < 1) stop("No .fasta files found in ", FASTA_DIR)
  labels <- vapply(files, label_from_file, character(1))
  names(files) <- labels

  # ordering
  if (!is.null(SEGMENT_ORDER)) {
    if (!setequal(SEGMENT_ORDER, labels))
      stop("SEGMENT_ORDER labels do not match FASTA labels: ",
           paste(setdiff(SEGMENT_ORDER, labels), collapse = ", "))
    labels <- SEGMENT_ORDER
  } else {
    labels <- sort(labels)
  }
  files <- files[labels]
  n_seg <- length(labels)
  message(sprintf("Found %d segment(s): %s", n_seg, paste(labels, collapse = ", ")))

  # ---- read + reconcile taxa ------------------------------------------------
  fastas <- lapply(files, read_fasta)
  taxa_sets <- lapply(fastas, names)
  common <- Reduce(intersect, taxa_sets)
  if (length(common) < 1) stop("No taxa shared across all segments.")
  for (k in seq_along(fastas)) {
    if (!setequal(taxa_sets[[k]], common))
      warning(sprintf("Segment '%s' taxa differ from the common set; using intersection (%d taxa).",
                      labels[k], length(common)))
  }
  # keep a stable taxon order (as in the first segment)
  taxa <- taxa_sets[[1]][taxa_sets[[1]] %in% common]

  # ---- dates ----------------------------------------------------------------
  dates_chr <- vapply(strsplit(taxa, DATE_SEP), `[`, character(1), DATE_FIELD)
  dates <- as.Date(dates_chr, format = DATE_FORMAT)
  if (any(is.na(dates)))
    stop("Could not parse dates for: ",
         paste(head(taxa[is.na(dates)]), collapse = ", "),
         "\nCheck DATE_SEP / DATE_FIELD / DATE_FORMAT.")

  # ---- rate-shift grid ------------------------------------------------------
  span <- as.numeric(max(dates) - min(dates)) / 365
  fine   <- seq(0, span, length.out = N_FINE)
  coarse <- seq(span * 1.5, DEEP_MAX, length.out = N_COARSE)
  rateshifts <- unique(round(c(fine, coarse), 6))
  indep_after <- if (is.null(INDEPENDENT_AFTER)) N_FINE else INDEPENDENT_AFTER
  message(sprintf("Sampling span: %.2f years; %d rate-shift points; independentAfter=%d",
                  span, length(rateshifts), indep_after))

  # ---- per-segment alignment lengths (delta-exchange weights) ---------------
  seq_lengths <- vapply(fastas, function(fa) nchar(fa[[which(names(fa) == taxa[1])]]),
                        integer(1))

  # ---- write cleaned/reconciled fastas into each analysis dir ---------------
  template_lines <- readLines(TEMPLATE, warn = FALSE)

  modes <- if (MODE == "both") c("dependent", "independent") else MODE
  for (mode in modes) {
    out_dir <- build_one_mode(mode, template_lines, taxa, dates_chr, rateshifts,
                              indep_after, seq_lengths, labels, n_seg)
    for (k in seq_len(n_seg)) {
      fa <- fastas[[k]][taxa]          # subset + reorder to common taxa
      write_fasta(fa, file.path(out_dir, paste0(labels[k], ".fasta")))
    }
  }
  message("Done.")
}

main()
