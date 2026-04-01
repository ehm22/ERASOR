library(DT)
library(shinythemes)
library(shiny)
library(shinydashboard)
library(GenomicFeatures)
library(GenomeInfoDb)   
library(AnnotationDbi)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(biomaRt)
library(Biostrings)
library(tidyverse)
library(cluster)
library(rlang)
library(dplyr)
library(purrr)
library(shinyBS)
library(openxlsx)
library(future)
library(future.apply)

####$$$####
library(stringr)
####$$$####

source("../tools/GGGenome_functions.R")
source("../tools/RNaseH_script.R")
source("../tools/Off_target_tissue.R")
source("../tools/Off_target_OMIM_Api.R")
source("../tools/Off_target_accessibility.R")

if (.Platform$OS.type == "windows") {
  plan(multisession, workers = 8)
} else {
  plan(multicore, workers = 8)
}
options(future.globals.maxSize = 6 * 1024^3, shiny.maxRequestSize = 8 * 1024^3)

# Helper function to track time of the run---------------------------------------
format_elapsed <- function(start_time, end_time = Sys.time()) {
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  hours   <- floor(elapsed / 3600)
  minutes <- floor((elapsed %% 3600) / 60)
  seconds <- floor(elapsed %% 60)
  
  parts <- c()
  if (hours > 0)   parts <- c(parts, sprintf("%dh", hours))
  if (minutes > 0) parts <- c(parts, sprintf("%dm", minutes))
  parts <- c(parts, sprintf("%ds", seconds))  
  
  paste(parts, collapse = " ")
}

# Helper function for renaming coumsn in RnaseH table

rename_rnaseh_cols <- function(df) {
  display_map <- c(
    position = "Position (nt)",
    average    = "RNase H score",
    window = "Nucleotide window"
  )
  
  nm <- names(df)
  nm2 <- ifelse(nm %in% names(display_map), display_map[nm], nm)
  names(df) <- make.unique(nm2, sep = " (dup) ")
  df
}

##############
# helper function to get output seqeunces from ggenome
render_subject_alignment_html <- function(subject_seq, match_string) {
  s_raw <- gsub("-", "", as.character(subject_seq))
  s <- strsplit(s_raw, "")[[1]]
  
  # Swap I and D for display so they match your intended biology
  ops_raw <- strsplit(as.character(match_string), "")[[1]]
  ops <- chartr("ID", "DI", ops_raw)
  
  out <- character(0)
  j <- 1L
  
  for (op in ops) {
    if (op == "=") {
      if (j <= length(s)) {
        out <- c(out, s[j])
        j <- j + 1L
      }
      
    } else if (op == "X") {
      if (j <= length(s)) {
        out <- c(out, paste0(
          "<span style='color:red;font-weight:bold;text-decoration:underline;'>",
          s[j],
          "</span>"
        ))
        j <- j + 1L
      }
      
    } else if (op == "I") {
      # query has a nucleotide that is missing in the subject
      out <- c(out,
               "<span style='color:black;font-weight:bold;'>-</span>"
      )
      
    } else if (op == "D") {
      # subject has an extra nucleotide compared to the query
      if (j <= length(s)) {
        out <- c(out, paste0(
          "<span style='color:#00cc00;font-weight:bold;text-decoration:underline;'>",
          s[j],
          "</span>"
        ))
        j <- j + 1L
      }
      
    } else {
      if (j <= length(s)) {
        out <- c(out, s[j])
        j <- j + 1L
      }
    }
  }
  
  paste0(out, collapse = "")
}
###################

# helper function for vcf and bcf file input for patient design
####$$$####

normalize_chr_style <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  ifelse(x %in% c("M", "MT"), "chrM", paste0("chr", x))
}

is_valid_region_string <- function(region) {
  grepl("^[^:]+:\\d+-\\d+$", trimws(region))
}

parse_region_string <- function(region) {
  region <- trimws(region)
  m <- regexec("^([^:]+):(\\d+)-(\\d+)$", region)
  hit <- regmatches(region, m)[[1]]
  if (length(hit) != 4) {
    stop("Region must have format chr:start-end")
  }
  list(
    chr = hit[2],
    start = as.integer(hit[3]),
    end = as.integer(hit[4])
  )
}

stage_patient_variant_files <- function(variant_input, index_input = NULL) {
  req(variant_input)
  
  workdir <- tempfile("patient_variant_job_")
  dir.create(workdir, recursive = TRUE)
  
  infile_name <- variant_input$name
  infile_path <- variant_input$datapath
  
  if (grepl("\\.vcf\\.gz$", infile_name, ignore.case = TRUE)) {
    local_variant <- file.path(workdir, "input.vcf.gz")
    variant_type <- "vcf.gz"
  } else if (grepl("\\.bcf$", infile_name, ignore.case = TRUE)) {
    local_variant <- file.path(workdir, "input.bcf")
    variant_type <- "bcf"
  } else if (grepl("\\.vcf$", infile_name, ignore.case = TRUE)) {
    local_variant <- file.path(workdir, "input.vcf")
    variant_type <- "vcf"
  } else {
    stop("Unsupported variant file type. Please upload .bcf, .vcf, or .vcf.gz")
  }
  
  ok <- file.copy(infile_path, local_variant, overwrite = TRUE)
  if (!ok) stop("Failed to stage uploaded variant file.")
  
  local_index <- NULL
  
  if (!is.null(index_input) && !is.null(index_input$datapath) && nzchar(index_input$datapath)) {
    idx_name <- index_input$name
    
    if (grepl("\\.tbi$", idx_name, ignore.case = TRUE)) {
      local_index <- paste0(local_variant, ".tbi")
    } else if (grepl("\\.csi$", idx_name, ignore.case = TRUE)) {
      local_index <- paste0(local_variant, ".csi")
    } else {
      stop("Index file must be .tbi or .csi")
    }
    
    ok_idx <- file.copy(index_input$datapath, local_index, overwrite = TRUE)
    if (!ok_idx) stop("Failed to stage uploaded index file.")
  }
  
  list(
    workdir = workdir,
    variant_path = local_variant,
    index_path = local_index,
    variant_type = variant_type,
    is_indexed = !is.null(local_index)
  )
}

get_vcf_samples <- function(variant_path) {
  bcftools_path <- get_bcftools_path()
  
  out <- system2(
    command = bcftools_path,
    args = c("query", "-l", variant_path),
    stdout = TRUE,
    stderr = TRUE
  )
  
  out <- trimws(out)
  out[nzchar(out)]
}

classify_variant_type <- function(ref, alt) {
  if (grepl(",", alt, fixed = TRUE)) {
    return("multiallelic")
  }
  
  if (grepl("^<", alt)) {
    return("symbolic")
  }
  
  if (nchar(ref) == 1 && nchar(alt) == 1) {
    return("SNV")
  }
  
  if (nchar(ref) < nchar(alt)) {
    return("insertion")
  }
  
  if (nchar(ref) > nchar(alt)) {
    return("deletion")
  }
  
  if (nchar(ref) == nchar(alt) && nchar(ref) > 1) {
    return("MNV_or_complex_substitution")
  }
  
  "other"
}

parse_gt_to_haplotypes <- function(ref, alt, gt) {
  gt <- trimws(as.character(gt))
  
  if (is.na(gt) || gt == "" || gt == "." || gt == "./." || gt == ".|.") {
    return(list(
      phased = FALSE,
      allele1_code = NA_character_,
      allele2_code = NA_character_,
      hap1_allele = NA_character_,
      hap2_allele = NA_character_,
      variant_on = "unknown",
      allele_note = "missing genotype"
    ))
  }
  
  phased <- grepl("\\|", gt)
  sep <- if (phased) "\\|" else "/"
  
  gt_parts <- strsplit(gt, sep)[[1]]
  if (length(gt_parts) == 1) {
    gt_parts <- c(gt_parts, NA_character_)
  }
  
  allele1_code <- gt_parts[1]
  allele2_code <- gt_parts[2]
  
  alt_parts <- strsplit(alt, ",", fixed = TRUE)[[1]]
  
  decode_allele <- function(code) {
    if (is.na(code) || code == ".") return(NA_character_)
    if (code == "0") return(ref)
    
    idx <- suppressWarnings(as.integer(code))
    if (is.na(idx) || idx < 1 || idx > length(alt_parts)) {
      return(NA_character_)
    }
    
    alt_parts[idx]
  }
  
  hap1_allele <- decode_allele(allele1_code)
  hap2_allele <- decode_allele(allele2_code)
  
  variant_on <- if (!phased) {
    "unphased"
  } else if (!is.na(allele1_code) && allele1_code != "0" &&
             !is.na(allele2_code) && allele2_code != "0") {
    "both_haplotypes"
  } else if (!is.na(allele1_code) && allele1_code != "0") {
    "haplotype_1"
  } else if (!is.na(allele2_code) && allele2_code != "0") {
    "haplotype_2"
  } else {
    "reference_only"
  }
  
  allele_note <- if (phased) {
    "phased genotype: haplotype assignment available"
  } else {
    "unphased genotype: haplotype assignment not known"
  }
  
  list(
    phased = phased,
    allele1_code = allele1_code,
    allele2_code = allele2_code,
    hap1_allele = hap1_allele,
    hap2_allele = hap2_allele,
    variant_on = variant_on,
    allele_note = allele_note
  )
}

read_patient_variants_bcftools <- function(variant_path, region_string, is_indexed = FALSE) {
  bcftools_path <- get_bcftools_path()
  
  region_flag <- if (isTRUE(is_indexed)) "-r" else "-t"
  
  out_file <- tempfile(fileext = ".vcf")
  err_file <- tempfile(fileext = ".log")
  
  args <- c(
    "view",
    region_flag, region_string,
    "-H",
    variant_path
  )
  
  status <- system2(
    command = bcftools_path,
    args = args,
    stdout = out_file,
    stderr = err_file
  )
  
  if (!identical(status, 0L)) {
    err_msg <- paste(readLines(err_file, warn = FALSE), collapse = "\n")
    stop(paste("bcftools failed:", err_msg))
  }
  
  if (!file.exists(out_file) || file.info(out_file)$size == 0) {
    return(tibble(
      chr = character(),
      pos = integer(),
      ref = character(),
      alt = character(),
      gt = character(),
      sample = character()
    ))
  }
  
  # Read raw VCF body lines
  lines <- readLines(out_file, warn = FALSE)
  if (length(lines) == 0) {
    return(tibble(
      chr = character(),
      pos = integer(),
      ref = character(),
      alt = character(),
      gt = character(),
      sample = character()
    ))
  }
  
  fields <- strsplit(lines, "\t", fixed = TRUE)
  
  df <- tibble(
    chr    = vapply(fields, `[`, character(1), 1),
    pos    = as.integer(vapply(fields, `[`, character(1), 2)),
    ref    = vapply(fields, `[`, character(1), 4),
    alt    = vapply(fields, `[`, character(1), 5),
    format = vapply(fields, function(x) if (length(x) >= 9) x[9] else NA_character_, character(1)),
    sample_field = vapply(fields, function(x) if (length(x) >= 10) x[10] else NA_character_, character(1))
  )
  
  extract_gt <- function(format_string, sample_string) {
    if (is.na(format_string) || is.na(sample_string)) return(NA_character_)
    
    fmt_parts <- strsplit(format_string, ":", fixed = TRUE)[[1]]
    samp_parts <- strsplit(sample_string, ":", fixed = TRUE)[[1]]
    
    gt_idx <- match("GT", fmt_parts)
    if (is.na(gt_idx) || gt_idx > length(samp_parts)) return(NA_character_)
    
    samp_parts[gt_idx]
  }
  
  df <- df %>%
    mutate(
      gt = purrr::map2_chr(format, sample_field, extract_gt),
      sample = NA_character_
    ) %>%
    select(chr, pos, ref, alt, gt, sample)
  
  gt_info <- purrr::pmap(
    list(df$ref, df$alt, df$gt),
    function(ref, alt, gt) parse_gt_to_haplotypes(ref, alt, gt)
  )
  
  gt_tbl <- tibble(
    phased = vapply(gt_info, `[[`, logical(1), "phased"),
    allele1_code = vapply(gt_info, `[[`, character(1), "allele1_code"),
    allele2_code = vapply(gt_info, `[[`, character(1), "allele2_code"),
    hap1_allele = vapply(gt_info, `[[`, character(1), "hap1_allele"),
    hap2_allele = vapply(gt_info, `[[`, character(1), "hap2_allele"),
    variant_on = vapply(gt_info, `[[`, character(1), "variant_on"),
    allele_note = vapply(gt_info, `[[`, character(1), "allele_note")
  )
  
  bind_cols(df, gt_tbl) %>%
    mutate(
      chr = normalize_chr_style(chr),
      variant_type = purrr::map2_chr(ref, alt, classify_variant_type)
    )
}

filter_snvs_to_gene <- function(variants_raw, gene_chr, gene_start, gene_end) {
  variants_raw %>%
    filter(
      chr == normalize_chr_style(gene_chr),
      pos >= gene_start,
      pos <= gene_end
    ) %>%
    distinct(chr, pos, ref, alt, .keep_all = TRUE)
}

make_patient_window <- function(gene_chr, gene_start, gene_end, snvs, flank = 0L, strand = "+") {
  flank <- as.integer(flank)
  if (nrow(snvs) == 0) {
    start_pos <- max(1L, gene_start - flank)
    end_pos   <- gene_end + flank
  } else {
    start_pos <- max(1L, min(snvs$pos) - flank)
    end_pos   <- max(snvs$pos) + flank
  }
  
  GenomicRanges::GRanges(
    seqnames = normalize_chr_style(gene_chr),
    ranges = IRanges::IRanges(start = start_pos, end = end_pos),
    strand = strand
  )
}

apply_snvs_to_window <- function(ref_seq, window_start, snvs) {
  seq_chars <- strsplit(as.character(ref_seq), "")[[1]]
  
  if (nrow(snvs) > 0) {
    for (i in seq_len(nrow(snvs))) {
      rel_pos <- snvs$pos[i] - window_start + 1L
      if (rel_pos >= 1L && rel_pos <= length(seq_chars)) {
        seq_chars[rel_pos] <- snvs$alt[i]
      }
    }
  }
  
  DNAString(paste(seq_chars, collapse = ""))
}

annotate_patient_snv_overlap <- function(target_annotation_patient, snvs_gene) {
  if (nrow(snvs_gene) == 0) {
    target_annotation_patient$patient_snv_count <- 0L
    target_annotation_patient$overlaps_patient_snv <- FALSE
    target_annotation_patient$distance_to_nearest_snv <- NA_integer_
    return(target_annotation_patient)
  }
  
  snv_pos <- snvs_gene$pos
  
  target_annotation_patient$patient_snv_count <- vapply(
    seq_len(nrow(target_annotation_patient)),
    function(i) {
      sum(
        snv_pos >= target_annotation_patient$chr_start[i] &
          snv_pos <= target_annotation_patient$chr_end[i]
      )
    },
    integer(1)
  )
  
  target_annotation_patient$overlaps_patient_snv <-
    target_annotation_patient$patient_snv_count > 0
  
  target_annotation_patient$distance_to_nearest_snv <- vapply(
    seq_len(nrow(target_annotation_patient)),
    function(i) {
      min(abs(snv_pos - target_annotation_patient$chr_start[i]),
          abs(snv_pos - target_annotation_patient$chr_end[i]))
    },
    integer(1)
  )
  
  target_annotation_patient
}
####$$$####

###$$$###
# helper with bcftools import
get_bcftools_path <- function() {
  bcftools_path <- Sys.which("bcftools")
  
  if (bcftools_path == "") {
    stop("bcftools is not installed in this Docker image or not available on PATH.")
  }
  
  bcftools_path
}
###$$$###


###$$$#####
# function to map the snvs to the reference sequence
apply_snvs_to_reference <- function(ref_seq, snv_df, region_start) {
  # Convert the DNAString reference sequence to a character vector of single bases
  seq_chars <- strsplit(as.character(ref_seq), "")[[1]]
  
  # If there are no SNVs, return the original sequence unchanged
  if (nrow(snv_df) == 0) {
    return(DNAString(paste(seq_chars, collapse = "")))
  }
  
  # Loop over each SNV and replace the reference base at that genomic position
  for (i in seq_len(nrow(snv_df))) {
    genomic_pos <- snv_df$pos[i]
    ref_base_vcf <- toupper(snv_df$ref[i])
    alt_base_vcf <- toupper(snv_df$alt[i])
    
    # Convert genomic position to position inside the selected sequence
    rel_pos <- genomic_pos - region_start + 1L
    
    # Skip anything outside the selected region
    if (rel_pos < 1L || rel_pos > length(seq_chars)) {
      next
    }
    
    # Optional safety check:
    # confirm the reference base in the genome matches the VCF REF allele
    ref_base_genome <- toupper(seq_chars[rel_pos])
    
    if (ref_base_genome != ref_base_vcf) {
      warning(
        paste0(
          "Reference mismatch at genomic position ", genomic_pos,
          ": genome has ", ref_base_genome,
          ", VCF REF is ", ref_base_vcf
        )
      )
      next
    }
    
    # Replace reference base with patient ALT base
    seq_chars[rel_pos] <- alt_base_vcf
  }
  
  DNAString(paste(seq_chars, collapse = ""))
}

#####$$$#####

# compute neurotox here for both tabs
# Make the tox score function

calculate_acute_neurotox <- function(xx) {
  # make sure input is in character format
  xx <- as.character(xx)
  
  # count the number of each nucleotide
  lf <- function(x) {
    x <- tolower(x)
    x <- strsplit(x, "")[[1]]
    x <- table(factor(x, levels = c("a", "c", "t", "g")))
    return(x)
  }
  cnt_nt <- as.data.frame(t(sapply(xx, lf)))
  # count number of nucleotides from the 3'-end untill it finds the first g
  gfree3 <- function(x) {
    x <- tolower(x)
    x <- strsplit(x, "")[[1]]
    tfg <- x == "g"
    if (sum(tfg) == 0) {
      l3 <- NA
    } else {
      posg <- c(1:length(x))[tfg]
      l3 <- length(x) - max(posg)
    }
    return(l3)
  }
  cnt_gfree3 <- sapply(xx, gfree3)
  cnt_gfree3[cnt_gfree3 > 20] <- 20      #Set max to 20
  cnt_gfree3[is.na(cnt_gfree3)] <- 20  #Set no g in ASO to 20
  ## Calculate final score based on trained parameters and return result
  calc_out <- round(
    136.0430 - 3.1263 * cnt_nt$a - 5.1100 * cnt_nt$c -
      4.7217 * cnt_nt$t - 10.1264 * cnt_nt$g + 1.3577 *
      cnt_gfree3,
    1
  )
  
  return(as.numeric(calc_out))
}

#-------------------------------------------------------------------------------

#################### test for sliding filter with numeric boxes ################
# --- Reusable range filter module (slider + numeric boxes) ----
rangeFilterServer <- function(id, min_allowed, max_allowed,
                              fixed = c("none", "left", "right")) {
  fixed <- match.arg(fixed)
  
  moduleServer(id, function(input, output, session) {
    
    if (fixed == "none") {
      
      # numeric -> slider
      observeEvent(c(input$min, input$max), ignoreInit = TRUE, {
        if (is.null(input$min) || is.null(input$max) || anyNA(c(input$min, input$max))) return()
        
        min_val <- max(min_allowed, min(input$min, input$max))
        max_val <- min(max_allowed, max(input$min, input$max))
        
        if (!identical(input$slider, c(min_val, max_val))) {
          updateSliderInput(session, "slider", value = c(min_val, max_val))
        }
      })
      
      # slider -> numeric
      observeEvent(input$slider, ignoreInit = TRUE, {
        rng <- input$slider
        if (is.null(rng) || length(rng) != 2 || anyNA(rng)) return()
        
        updateNumericInput(session, "min", value = rng[1])
        updateNumericInput(session, "max", value = rng[2])
      })
      
      reactive(input$slider)
      
    } else if (fixed == "left") {
      
      get_max <- function(x) {
        if (is.null(x) || anyNA(x)) return(NULL)
        if (length(x) == 2) x[2] else x[1]
      }
      
      observeEvent(input$slider, ignoreInit = TRUE, {
        v0 <- get_max(input$slider); if (is.null(v0)) return()
        v  <- min(max_allowed, max(min_allowed, v0))
        
        if (!identical(input$max, v))    updateNumericInput(session, "max", value = v)
        if (!identical(get_max(input$slider), v)) updateSliderInput(session, "slider", value = v)
      })
      
      observeEvent(input$max, ignoreInit = TRUE, {
        if (is.null(input$max) || anyNA(input$max)) return()
        v <- min(max_allowed, max(min_allowed, input$max))
        
        if (!identical(get_max(input$slider), v)) updateSliderInput(session, "slider", value = v)
        if (!identical(input$max, v))             updateNumericInput(session, "max", value = v)
      })
      
      reactive(c(min_allowed, get_max(input$slider)))
      
    } else { # fixed == "right"
      
      get_min <- function(x) {
        if (is.null(x) || anyNA(x)) return(NULL)
        x[1]
      }
      
      observeEvent(input$slider, ignoreInit = TRUE, {
        v0 <- get_min(input$slider); if (is.null(v0)) return()
        v  <- min(max_allowed, max(min_allowed, v0))
        
        if (!identical(input$min, v))    updateNumericInput(session, "min", value = v)
        if (!identical(get_min(input$slider), v)) updateSliderInput(session, "slider", value = v)
      })
      
      observeEvent(input$min, ignoreInit = TRUE, {
        if (is.null(input$min) || anyNA(input$min)) return()
        v <- min(max_allowed, max(min_allowed, input$min))
        
        if (!identical(get_min(input$slider), v)) updateSliderInput(session, "slider", value = v)
        if (!identical(input$min, v))             updateNumericInput(session, "min", value = v)
      })
      
      reactive(c(get_min(input$slider), max_allowed))
    }
  })
}
################################################################################

function(input, output, session) {
  
  # ----------------------------------- Notifications UI -----------------------
  # ----------------------------- GGGenome startup check -----------------------
  check_gggenome <- function() {
    tryCatch({
      all_offt("ACTGACTGACTGACTGACTG", mismatches_allowed = 0)
      TRUE
    }, error = function(e) {
      FALSE
    })
  }
  
  gggenome_available <- reactiveVal(NULL)
  
  observeEvent(TRUE, {
    available <- check_gggenome()
    gggenome_available(available)
    
    if (isTRUE(available)) {
      showNotification(
        "GGGenome is available.",
        type = "message",
        duration = 20
      )
    } else {
      showNotification(
        "GGGenome is currently unavailable. Running the script might lead to a crash/disconnect.",
        type = "error",
        duration = NULL
      )
    }
    
  }, once = TRUE)
  
  
  # ----------------------------- BioMart startup check -----------------------
  
  check_biomart <- function() {
    tryCatch({
      martHS <- biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl",
        host    = "https://www.ensembl.org"
      )
      
      martMM <- biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "mmusculus_gene_ensembl",
        host    = "https://www.ensembl.org"
      )
      
      # 1) human -> mouse ortholog query
      ortho <- biomaRt::getBM(
        attributes = "mmusculus_homolog_ensembl_gene",
        filters    = "ensembl_gene_id",
        values     = "ENSG00000139618",
        mart       = martHS,
        bmHeader   = FALSE
      )
      
      mm_id <- ortho$mmusculus_homolog_ensembl_gene
      mm_id <- mm_id[!is.na(mm_id) & nzchar(mm_id)]
      if (length(mm_id) == 0) return(FALSE)
      
      # 2) mouse sequence query
      mm_seq <- biomaRt::getBM(
        attributes = c("gene_exon_intron", "ensembl_gene_id"),
        filters    = "ensembl_gene_id",
        values     = mm_id[1],
        mart       = martMM
      )
      
      if (!is.data.frame(mm_seq) ||
          nrow(mm_seq) == 0 ||
          !"gene_exon_intron" %in% names(mm_seq)) {
        return(FALSE)
      }
      
      # 3) polymorphism query
      pm <- biomaRt::getBM(
        attributes = c("minor_allele_freq", "chromosome_start"),
        filters    = "ensembl_gene_id",
        values     = "ENSG00000139618",
        mart       = martHS
      )
      
      if (!is.data.frame(pm) ||
          !all(c("minor_allele_freq", "chromosome_start") %in% names(pm))) {
        return(FALSE)
      }
      
      TRUE
    }, error = function(e) {
      message("BioMart startup check failed: ", e$message)
      FALSE
    })
  }
  
  biomart_available <- reactiveVal(NULL)
  
  observeEvent(TRUE, {
    available <- check_biomart()
    biomart_available(available)
    
    if (isTRUE(available)) {
      showNotification(
        "BioMart is available.",
        type = "message",
        duration = 8
      )
    } else {
      showNotification(
        "connection to BioMart is unstable. Running the script might lead to a crash/disconnect.",
        type = "warning",
        duration = NULL
      )
    }
  }, once = TRUE)
  
  ############### ensembl autocomplete #############
  # ---- Ensembl autocomplete: build dropdown choices once at startup ----
  # Gene symbols are loaded from a precomputed local .rds file instead of
  # querying Ensembl live at startup. This avoids missing symbols in the
  # containerized user-facing Podman deployment.
  
  txdb_for_choices <- tryCatch({
    loadDb("/opt/ERASOR/txdb_hsa_biomart.db")
  }, error = function(e) {
    message("Primary txdb path not found for choices, trying local path...")
    loadDb("../txdb_hsa_biomart.db")
  })
  
  gdb_for_choices <- genes(txdb_for_choices)
  
  chr_to_keep <- c(as.character(1:22), "X", "Y", "MT")
  gdb_for_choices <- gdb_for_choices[seqnames(gdb_for_choices) %in% chr_to_keep]
  
  ens_ids <- names(gdb_for_choices)
  
  gene_map_file <- if (file.exists("/srv/shiny-server/ERASOR/gene_symbol_map.rds")) {
    "/srv/shiny-server/ERASOR/gene_symbol_map.rds"
  } else if (file.exists("/opt/ERASOR/gene_symbol_map.rds")) {
    "/opt/ERASOR/gene_symbol_map.rds"
  } else if (file.exists("../gene_symbol_map.rds")) {
    "../gene_symbol_map.rds"
  } else if (file.exists("gene_symbol_map.rds")) {
    "gene_symbol_map.rds"
  } else {
    NULL
  }
  
  gene_symbols <- rep(NA_character_, length(ens_ids))
  
  if (!is.null(gene_map_file)) {
    bm_map <- readRDS(gene_map_file)
    bm_map <- as.data.frame(bm_map, stringsAsFactors = FALSE)
    bm_map <- bm_map[!duplicated(bm_map$ensembl_gene_id), , drop = FALSE]
    
    idx <- match(ens_ids, bm_map$ensembl_gene_id)
    gene_symbols <- bm_map$hgnc_symbol[idx]
  } else {
    message("WARNING: gene_symbol_map.rds not found; falling back to Ensembl IDs only.")
  }
  
  gene_map <- data.frame(
    ensembl = ens_ids,
    symbol  = gene_symbols,
    stringsAsFactors = FALSE
  )
  
  gene_map$label <- ifelse(
    is.na(gene_map$symbol) | gene_map$symbol == "",
    gene_map$ensembl,
    paste0(gene_map$symbol, " (", gene_map$ensembl, ")")
  )
  
  gene_choices_df <- data.frame(
    label   = gene_map$label,
    value   = gene_map$ensembl,
    symbol  = ifelse(is.na(gene_map$symbol), "", gene_map$symbol),
    ensembl = gene_map$ensembl,
    stringsAsFactors = FALSE
  )
  
  observe({
    updateSelectizeInput(
      session,
      inputId  = "ensemble_id_input",
      choices  = gene_choices_df,
      selected = "ENSG00000174775",
      server   = TRUE,
      options  = list(
        valueField  = "value",
        labelField  = "label",
        searchField = c("label", "symbol", "ensembl"),
        placeholder = "Type gene symbol or Ensembl ID"
      )
    )
  })
  
  ######$$$########
  observe({
    updateSelectizeInput(
      session,
      inputId  = "ensemble_id_input_patient",
      choices  = gene_choices_df,
      selected = "ENSG00000174775",
      server   = TRUE,
      options  = list(
        valueField  = "value",
        labelField  = "label",
        searchField = c("label", "symbol", "ensembl"),
        placeholder = "Type gene symbol or Ensembl ID"
      )
    )
  })
  ######$$$########
  
  session$onFlushed(function() {
    shinyjs::hide("app_startup_loading")
    
    showNotification(
      "Finished loading",
      type = "default",
      duration = 10,
      closeButton = TRUE
    )
  }, once = TRUE)
  ######################################################################
  
  # ----------------------------- ASO parser preview -----------------------------
  parsed_asos <- reactive({
    req(input$single_aso_input)
    
    raw_text <- toupper(trimws(if (is.null(input$aso_seq_input)) "" else input$aso_seq_input))
    
    if (raw_text == "") {
      return(character(0))
    }
    
    asos <- regmatches(raw_text, gregexpr("[ACGT]+", raw_text))[[1]]
    asos <- unique(asos[nchar(asos) > 0])
    asos
  })
  
  output$parsed_asos <- renderText({
    req(input$single_aso_input)
    
    asos <- parsed_asos()
    
    if (length(asos) == 0) {
      return("No valid ASO sequences detected.")
    }
    
    paste0(seq_along(asos), ". ", asos, collapse = "\n")
  })
  # -------------------------------------------------------------------
  
  # Disable filters when selecting single ASO input
  
  observeEvent(input$single_aso_input, ignoreInit = TRUE, {
    updateCheckboxInput(session, "perfect_input",  value = FALSE)
    updateCheckboxInput(session, "mismatch_input", value = FALSE)
    updateCheckboxInput(session, "Poly_input", value = FALSE)
    updateCheckboxInput(session, "tox_input", value = FALSE)
    updateCheckboxInput(session, "Accessibility_input", value = FALSE)
    updateCheckboxInput(session, "gc_input", value = FALSE)
  })
  
  # Enable filters when disabling single ASO input
  
  observeEvent(input$single_aso_input, ignoreInit = TRUE, {
    if (!isTRUE(input$single_aso_input)) {
      updateCheckboxInput(session, "perfect_input",  value = TRUE)
      updateCheckboxInput(session, "mismatch_input", value = TRUE)
      updateCheckboxInput(session, "Poly_input", value = TRUE)
      updateCheckboxInput(session, "tox_input", value = TRUE)
      updateCheckboxInput(session, "Accessibility_input", value = TRUE)
      updateCheckboxInput(session, "gc_input", value = TRUE)
    }
  })
  
  # make output table collapsable
  show_all_cols <- reactiveVal(FALSE)
  
  observeEvent(input$toggle_cols, {
    current <- show_all_cols()
    show_all_cols(!current)
    
    # update button label
    if (current) {
      updateActionButton(inputId = "toggle_cols", label = "Show extended data")
    } else {
      updateActionButton(inputId = "toggle_cols", label = "Hide extended data ")
    }
  })
  
  
  # reset to standard values 
  
  observeEvent(input$reset_defaults, {
    shinyjs::reset("filters")
  })
  
  # sliding filter with numeric boxes
  tox_range <- rangeFilterServer("tox_score",    0, 136, fixed = "right")
  gc_range  <- rangeFilterServer("gc_content",   0, 100, fixed = "none")
  pm_range  <- rangeFilterServer("pm_freq",      0, 1,   fixed = "left")
  acc_range <- rangeFilterServer("accessibility",0, 1,   fixed = "none")
  perf_range <- rangeFilterServer("perfect_hits", 0, 200, fixed = "left")
  mm_range   <- rangeFilterServer("mismatch_hits",0, 500, fixed = "left")
  
  # starts timer, shows that the scirpt started and is loadin gthe progress bar
  
  observeEvent(input$run_button, {
    if (isTRUE(input$single_aso_input) && length(parsed_asos()) == 0) {
      showNotification(
        "No valid ASO sequences detected. Please enter sequences containing only A, C, G, and T.",
        type = "error",
        duration = 8
      )
      return(NULL)
    }
    start_time <- Sys.time()
    
    showNotification(
      "Script started",
      type = "default",
      duration = 10,
      closeButton = TRUE
    )
    
    
    withProgress(message = "Running analysis...", value = 0, {
      run_step <- function(label, amount, expr) {
        incProgress(0, detail = paste0(label, "..."))
        on.exit(
          incProgress(amount, detail = paste0(label)),
          add = TRUE
        )
        force(expr)
      }
      
      # ----------------------------------- Functions ------------------------------
      
      if (input$linux_input == TRUE) {
        # 2.5.1 Predict accessibility for an RNA target
        RNAplfold_R = function(seq.char,
                               L.in = 40,
                               W.in = 80,
                               u.in = 16) {
          cmmnd2 = paste("RNAplfold -L", L.in, "-W", W.in, "-u", u.in)
          seq.char = as.character(seq.char)
          cat(seq.char, file = paste("|", cmmnd2, sep = ""))
          acc.tx = read.delim(
            "plfold_lunp",
            as.is     = TRUE,
            skip      = 2,
            header    = FALSE,
            row.names = 1
          )
          acc.tx = acc.tx[, colSums(is.na(acc.tx)) != nrow(acc.tx)]
          colnames(acc.tx) = 1:ncol(acc.tx)
          file.remove("plfold_lunp")
          file.remove("plfold_dp.ps")
          return(acc.tx)
        }
        
        # 2.5.3 Predict duplex formation and self folding of oligonucleotides
        RNAduplex_R = function(seqs) {
          sys_cmd = system('RNAduplex',
                           input = c(seqs, seqs),
                           intern = TRUE)
          as.numeric(regmatches(sys_cmd, regexpr("-?\\d+\\.\\d+", sys_cmd)))
        }
        
        RNAselffold_R = function (seqs) {
          output = system('RNAfold --noPS',
                          input = c(seqs),
                          intern = T)
          output = unlist(strsplit(output[grepl('[0-9]', output)], '[(]'))
          as.double(gsub(' |[)]', '', output[grepl('[0-9]', output)]))
        }
      }
      
      # ----------------------------------- Data setup -----------------------------
      # Store all human pre-mRNA sequences
      
      txdb_hsa <- run_step("Retrieving sequences from Ensembl", 0.05, {
        tryCatch({
          loadDb("/opt/ERASOR/txdb_hsa_biomart.db")
        }, error = function(e) {
          message("Primary txdb path not found, trying local path...")
          loadDb("../txdb_hsa_biomart.db")
        })
      })
      
      # ----------------------------------- milestone 1 ----------------------------
      print("milestone1: loaded human database object")
      
      # Extract the genes
      gdb_hsa <- genes(txdb_hsa)
      
      # Harmonize seqname style between TxDb and BSgenome
      # (this fixes the 'invalid sequence name: 1' error)
      seqlevelsStyle(gdb_hsa) <- seqlevelsStyle(Hsapiens)
      
      # Keep only chromosomes present in BOTH objects
      common_chrs <- intersect(seqlevels(gdb_hsa), seqlevels(Hsapiens))
      gdb_hsa <- keepSeqlevels(gdb_hsa, common_chrs, pruning.mode = "coarse")
      
      # ----------------------------------- milestone 2 ----------------------------
      print("milestone2: Subsetted genes from specified chromosomes")
      
      # Get the sequences *
      HS <- getSeq(Hsapiens, gdb_hsa)
      
      # ----------------------------------- milestone 3 ----------------------------
      print("milestone3: Saved human gene sequences")
      # Target collect the pre-mRNA sequence
      
      # Define wanted Ensembl ID
      ensembl_ID = input$ensemble_id_input
      
      # Retrieve a specific RNA target using the Ensembl ID
      RNA_target = HS[names(HS) == ensembl_ID]
      
      
      # ----------------------------------- milestone 4 ----------------------------
      print("milestone4: Saved input RNA target (ENSEMBL) information: \n")
      print(RNA_target)
      #filters on ensembl ID
      target_ranges = gdb_hsa[names(gdb_hsa) == ensembl_ID]
      
      
      # Extracts the chromosome name for target range,extracts the start and end
      # Position of the genomic region and note from which strand it is.
      # Positive strand 1 ('+'), negative strand -1 ('-'), or unspecified 0 ('*').
      chr_coord = c(
        chr = as.numeric(as.character(seqnames(target_ranges))),
        start = start(target_ranges),
        end = end(target_ranges),
        strand = ifelse(strand(target_ranges) == "+", 1, -1)
      )
      
      # ----------------------------------- milestone 5 ----------------------------
      print("milestone5: Extracted chromosome coordinates of RNA target")
      
      # Set the width of rna_target as value L
      l = width(RNA_target)
      
      # Note length of subsequence
      oligo_lengths = input$oligo_length_range[1]:input$oligo_length_range[2]
      
      # Construct the dataframe
      target_annotation = lapply(oligo_lengths, function(i) {
        tibble(start = 1:(l - i + 1), length = i)
      }) %>%
        bind_rows() %>%
        mutate(end = start + length - 1)
      
      # Adds "name" sequence to the dataframe
      target_regions = DNAStringSet(RNA_target[[1]],
                                    start = target_annotation$start,
                                    width = target_annotation$length)
      names(target_regions) = as.character(target_regions)
      
      target_annotation$name = names(target_regions)
      
      # initliaze pm columns to avoi dcrash when pm calcualtation is disabled 
      target_annotation$PM_tot_freq <- NA_real_
      target_annotation$PM_max_freq <- NA_real_
      target_annotation$PM_count    <- NA_real_
      
      if (isTRUE(input$single_aso_input)) {
        aso_matches <- parsed_asos()
        
        if (length(aso_matches) > 0) {
          rev_comps <- as.character(
            reverseComplement(DNAStringSet(aso_matches))
          )
          
          target_annotation <- target_annotation[target_annotation$name %in% rev_comps, ]
          target_annotation$input_order <- match(target_annotation$name, rev_comps)
          target_annotation <- target_annotation[order(target_annotation$input_order), ]
        } else {
          target_annotation <- target_annotation[0, ]
        }
      }
      # ----------------------------------- milestone 6 ----------------------------
      print("milestone6: Enumerated all possible ASO target sequences")
      
      
      # ----------------------------------- milestone 7 ----------------------------
      print("milestone 7: Prefiltered Oligo sequences ending with G")
      
      if (input$linux_input == TRUE) {
        # 3.4 Estimate Transcript Accessibility for the RNA Target at Single-Nucleotide Resolution
        accessibility = RNAplfold_R(RNA_target, u.in = max(oligo_lengths)) %>%
          as_tibble() %>%
          mutate(end = 1:l) %>%
          gather(length, accessibility, -end) %>%
          mutate(length = as.double(length))
        
        target_annotation = left_join(target_annotation, accessibility, by =
                                        c('length', 'end'))
        
      }
      # ----------------------------------- milestone 8 ----------------------------
      print("milestone 8: Calculated ViennaRNA accessibility score and filtering")
      
      nucleobase_seq = reverseComplement(target_regions)
      
      # Voeg ze toe aan uni_tar gekoppeld aan name
      target_annotation$oligo_seq = as.character(nucleobase_seq[target_annotation$name])
      
      # Toxicity score acute neurotox score (desired >70)
      target_annotation$tox_score = calculate_acute_neurotox(target_annotation$oligo_seq)
      
      # ----------------------------------- milestone 8.5 --------------------------
      oligo_dna <- DNAStringSet(target_annotation$oligo_seq)
      gc_counts <- letterFrequency(oligo_dna, c("G", "C"))
      oligo_len <- width(oligo_dna)
      target_annotation$gc_content <- rowSums(gc_counts) / oligo_len * 100
      
      # ----------------------------------- milestone 9 ----------------------------
      print("milestone 9: Calculated toxicity score and filtering")
      
      # Define motif weights
      motif_weights <- c(
        CCAC = +0.3,
        TCCC = +0.3,
        ACTC = +0.2,
        GCCA = +0.2,
        CTCT = +0.1,
        GGGG = -0.2,
        ACTG = -0.2,
        TAA  = -0.2,
        CCGG = -0.1,
        AAA  = -0.1
      )
      
      seqs <- DNAStringSet(target_annotation$oligo_seq)
      
      motif_counts_mat <- sapply(
        names(motif_weights),
        function(m) vcountPattern(m, seqs),
        simplify = "matrix"
      )
      
      if (nrow(target_annotation) == 1) {
        motif_counts_mat <- t(motif_counts_mat)
      }
      
      motif_counts_df <- as.data.frame(motif_counts_mat, stringsAsFactors = FALSE)
      
      motif_scores <- as.matrix(motif_counts_mat) %*% as.numeric(motif_weights)
      motif_scores <- round(motif_scores[, 1], 6)
      
      target_annotation <- bind_cols(
        target_annotation,
        motif_counts_df,
        motif_cor_score = motif_scores
      )
      
      
      # ----------------------------------- milestone 9.5 --------------------------
      print("milestone 9.5: Calculated motif correlation score")
      
      # Bereken aantal "cg"
      target_annotation$CGs = run_step("Getting mouse ortholog data", 0.10, {(target_annotation$length -
                                                                                nchar(gsub('CG', '', target_annotation$name))) / 2
      })
      
      # ----------------------------------- milestone 10 ---------------------------
      print("milestone 10: Calculated CG motifs")
      
      options(timeout = 60)
      
      # Define the marts for mmusculus and hsapiens
      martHS <- useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl",
        host    = "https://www.ensembl.org"
      )
      
      martMM <- useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "mmusculus_gene_ensembl",
        host    = "https://www.ensembl.org"
      )
      
      # Get the orthologous Ensembl gene for the provided human Ensembl ID
      ortho_ENS = getBM(attributes = "mmusculus_homolog_ensembl_gene",
                        filters = "ensembl_gene_id",
                        values = ensembl_ID, mart = martHS,
                        bmHeader = FALSE)
      
      # ----------------------------------- milestone 11 ---------------------------
      print("milestone 11: Get mouse ortholog data for genes")
      
      RNA_target_mouse = DNAStringSet(
        getBM(attributes = c("gene_exon_intron","ensembl_gene_id"),
              filters = "ensembl_gene_id",
              values =
                ortho_ENS$mmusculus_homolog_ensembl_gene,
              mart = martMM)$gene_exon_intron)
      
      # ----------------------------------- milestone 12 ---------------------------
      print("milestone 12: ")
      
      # Obtain all human polymorphisms for the RNA target
      if (input$polymorphism_input == TRUE) {
        PMs = getBM(
          attributes = c("minor_allele_freq", "chromosome_start"),
          filters = "ensembl_gene_id",
          values = ensembl_ID,
          mart = martHS
        ) %>%
          as_tibble() %>%
          arrange(chromosome_start, desc(minor_allele_freq)) %>%
          filter(!is.na(minor_allele_freq),
                 !duplicated(chromosome_start)) %>%
          rename(chr_start = chromosome_start, PM_freq = minor_allele_freq)
      }
      
      ##If Ensembl is offline and still want to test -> load in manual test data.
      #PMs <- read.csv("~/PMs.csv")
      
      
      # ----------------------------------- milestone 13 ---------------------------
      print("milestone 13: Get human polymorfisms for RNA target")
      
      # Count Nucleobase sequence occurrences
      
      # Get the sequences
      tr = target_annotation$name
      
      # Count the sequences by making it a table
      replica = table(tr)
      
      # Save it as NoRepeats
      target_annotation$NoRepeats = as.vector(replica[tr])
      
      # ----------------------------------- milestone 14 ---------------------------
      print("milestone 14: Count amount of times ASO sequence is repeated in target gene")
      
      # High-Frequency Polymorphisms analysis
      
      # Correcting end and start cord based on direction
      if (chr_coord['strand'] == 1) {
        target_annotation$chr_start = chr_coord['start'] + target_annotation$start - 1
        target_annotation$chr_end = chr_coord['start'] + target_annotation$end - 1
      } else {
        target_annotation$chr_start = chr_coord['end'] - target_annotation$end + 1
        target_annotation$chr_end = chr_coord['end'] - target_annotation$start + 1
      }
      
      # ----------------------------------- milestone 15 ---------------------------
      print("milestone 15: Corrected start en end cord based on direction")
      
      # Keep unique names only and extract
      # Information base on chr_start from target.
      if (input$polymorphism_input == TRUE) {
        PM_freq = PMs %>%
          mutate(name = map(chr_start, function(X) {
            filter(target_annotation, chr_start <= X, chr_end >= X) %>%
              select(name, chr_start_anno = chr_start)
          })) %>%
          unnest_legacy() %>%
          rename(chr_start_PM = chr_start) %>%
          group_by(name, chr_start_anno) %>%
          summarise(
            PM_tot_freq = 1 - prod(1 - PM_freq),
            PM_max_freq = max(PM_freq),
            PM_count = n()
          )
        
        # ----------------------------------- milestone 16 -------------------------
        print("milestone 16: Retrieved point mutation frequency count and chance")
        
        # Combine data frames
        target_annotation = left_join(
          target_annotation,
          PM_freq,
          by = c("name" = "name", "chr_start" = "chr_start_anno")
        ) %>%
          mutate(
            PM_tot_freq = coalesce(PM_tot_freq.y, PM_tot_freq.x),
            PM_max_freq = coalesce(PM_max_freq.y, PM_max_freq.x),
            PM_count    = coalesce(PM_count.y, PM_count.x)
          ) %>%
          select(-PM_tot_freq.x, -PM_tot_freq.y,
                 -PM_max_freq.x, -PM_max_freq.y,
                 -PM_count.x, -PM_count.y)
      }
      
      
      # ----------------------------------- milestone 17 ---------------------------
      print("milestone 17: Joined PM freq table with target annotation")
      
      # Match RNA Target Regions to the Mouse Ortholog
      
      # Get length
      lm = width(RNA_target_mouse)
      
      # Make table of mouse information
      MM_tab = lapply(oligo_lengths, function(i) {
        tibble(st = 1:(lm - i + 1), w = i)
      }) %>%
        bind_rows()
      
      # ----------------------------------- milestone 18.1 -------------------------
      print("milestone 18.1: Match RNA target to mouse ortholog")
      
      # Makes DNAStringSet object with mouse info.
      RNAsitesMM = DNAStringSet(RNA_target_mouse[[1]],
                                start = MM_tab$st,
                                width = MM_tab$w)
      
      # Adds if conserved in mouse.
      
      target_annotation$conserved_in_mmusculus = run_step("Calculating secondary structure characteristics", 0.25, {target_annotation$name %in% RNAsitesMM
      })
      
      # ----------------------------------- milestone 18.2 -------------------------
      print("milestone 18.2: If selected: Matched mouse ortholog to target gene")
      
      if (input$linux_input == TRUE) {
        # Deze nog toevoegen wanneer viennaRNA werkt
        target_annotation$sec_energy = RNAselffold_R(target_annotation$oligo_seq)
        target_annotation$duplex_energy = RNAduplex_R(target_annotation$oligo_seq)
      }
      # ----------------------------------- milestone 19 ---------------------------
      print("milestone 19: Calculated secondary and duplex energy of ASO seq")
      
      unfiltered_total_data <- target_annotation
      
      nseq_prefilter <- nrow(target_annotation)
      
      nseq_ending_G <- nseq_prefilter - nrow(target_annotation %>%
                                               filter(!grepl("^C", name))
      )
      
      # nseq_toxscore <- nseq_prefilter - nrow(filter_function(target_annotation, input$numeric_input_e, "tox_score", input$dropdown_input_e))
      ## new with slidng
      tox_rng <- tox_range()  
      nseq_toxscore <- nseq_prefilter - nrow(
        target_annotation %>%
          dplyr::filter(
            !is.na(tox_score),
            tox_score >= tox_rng[1],
            tox_score <= tox_rng[2]
          )
      )
      
      gc_rng <- gc_range()
      nseq_gc <- nseq_prefilter - nrow(
        target_annotation %>%
          dplyr::filter(
            !is.na(gc_content),
            gc_content >= gc_rng[1],
            gc_content <= gc_rng[2]
          )
      )
      ##
      
      nseq_pmfreq <- if (isTRUE(input$polymorphism_input)) {
        pm_rng <- pm_range()
        
        nseq_prefilter - nrow(
          target_annotation %>%
            dplyr::mutate(PM_tot_freq = tidyr::replace_na(PM_tot_freq, 0)) %>%
            dplyr::filter(PM_tot_freq >= pm_rng[1], PM_tot_freq <= pm_rng[2])
        )
      } else {
        NA_integer_
      }
      
      
      nseq_accessible <- if (isTRUE(input$linux_input)) {
        acc_rng <- acc_range()
        nseq_prefilter - nrow(
          target_annotation %>%
            dplyr::filter(
              !is.na(accessibility),
              accessibility >= acc_rng[1],
              accessibility <= acc_rng[2]
            )
        )
      } else {
        NA_integer_
      }
      
      filter_numbers <- tibble(
        metric = c(
          "Prefiltered",
          "ASO ending with G",
          "Tox score",
          "GC content",
          "PM frequency",
          "Accessibility"
        ),
        count = c(
          nseq_prefilter,
          nseq_ending_G,
          nseq_toxscore,
          nseq_gc,
          nseq_pmfreq,
          nseq_accessible
        )
      )
      
      output$unfiltered_results_table <- renderDT ({
        datatable(
          filter_numbers,
          rownames = FALSE,
          options = list(dom = "t", ordering = FALSE),
          class = "compact stripe"
        )
      })
      
      # ----------------------------------- filtering ------------------------------
      
      apply_range <- function(df, col, rng) {
        df %>%
          dplyr::filter(
            !is.na(.data[[col]]),
            .data[[col]] >= rng[1],
            .data[[col]] <= rng[2]
          )
      }
      
      ta <- target_annotation  # startdataset
      print(nrow(target_annotation))
      # 1) ASO ending with G filter
      if (isTRUE(input$ASO_ending_G)) {
        ta_prev <- ta
        ta <- ta %>% filter(!grepl("^C", name))
        
        if (nrow(ta) == 0) {
          ta <- ta_prev
          message("Filter 'ASO ending with G' removed all rows; reverting to previous dataset.")
          showNotification("Filter 'ASO ending with G' removed all rows; reverting to previous dataset.", type = "warning")
        }
      }
      
      # 2) tox_score filter
      if (isTRUE(input$tox_input)) {
        ta_prev <- ta
        rng <- tox_range()  # c(min, max)
        
        ta <- ta %>%
          dplyr::filter(
            !is.na(tox_score),
            tox_score >= rng[1],
            tox_score <= rng[2]
          )
        
        if (nrow(ta) == 0) {
          ta <- ta_prev
          message("Filter 'tox_score' removed all rows; reverting to previous dataset.")
          showNotification("Filter 'tox_score' removed all rows; reverting to previous dataset.", type = "warning")
        }
      }
      # 2b) GC content filter 
      if (isTRUE(input$gc_input)) {   
        ta_prev <- ta
        rng <- gc_range()  # c(min, max)
        
        ta <- ta %>%
          dplyr::filter(
            !is.na(gc_content),
            gc_content >= rng[1],
            gc_content <= rng[2]
          )
        
        if (nrow(ta) == 0) {
          ta <- ta_prev
          message("Filter 'gc_content' removed all rows; reverting to previous dataset.")
          showNotification("Filter 'GC content' removed all rows; reverting to previous dataset.", type = "warning")
        }
      }
      ###
      
      # 3) accessibility filter (alleen als linux + checkbox aanstaan)
      if (isTRUE(input$linux_input) && isTRUE(input$Accessibility_input)) {
        ta_prev <- ta
        if (isTRUE(input$linux_input) && isTRUE(input$Accessibility_input)) {
          ta_prev <- ta
          rng <- acc_range()
          ta <- apply_range(ta, "accessibility", rng)
          
          if (nrow(ta) == 0) {
            ta <- ta_prev
            showNotification("Filter 'accessibility' removed all rows; reverting.", type = "warning")
          }
        }
        
        if (nrow(ta) == 0) {
          ta <- ta_prev
          message("Filter 'accessibility' removed all rows; reverting to previous dataset.")
          showNotification("Filter 'accessibility' removed all rows; reverting to previous dataset.", type = "warning")
        }
      }
      
      # 4) polymorphism NA->0 (geen filter; alleen transform)
      if (isTRUE(input$polymorphism_input)) {
        ta <- ta %>%
          mutate(across(c(PM_max_freq, PM_tot_freq, PM_count), ~ tidyr::replace_na(., 0.0)))
      }
      
      # 5) polymorphism frequency filter
      if (isTRUE(input$polymorphism_input) && isTRUE(input$Poly_input)) {
        ta_prev <- ta
        if (isTRUE(input$polymorphism_input)) {
          ta <- ta %>%
            mutate(across(c(PM_max_freq, PM_tot_freq, PM_count), ~ tidyr::replace_na(., 0.0)))
        }
        
        if (isTRUE(input$polymorphism_input) && isTRUE(input$Poly_input)) {
          ta_prev <- ta
          rng <- pm_range()  # c(min,max), min fixed at 0
          ta <- apply_range(ta, "PM_tot_freq", rng)
          
          if (nrow(ta) == 0) {
            ta <- ta_prev
            showNotification("Filter 'PM_tot_freq' removed all rows; reverting.", type = "warning")
          }
        }
        
        if (nrow(ta) == 0) {
          ta <- ta_prev
          message("Filter 'PM_tot_freq' removed all rows; reverting to previous dataset.")
          showNotification("Filter 'PM_tot_freq' removed all rows; reverting to previous dataset.", type = "warning")
        }
      }
      
      # 6) conserved filter
      if (isTRUE(input$Conserved_input)) {
        ta_prev <- ta
        ta <- ta %>% filter(conserved_in_mmusculus == TRUE)
        
        if (nrow(ta) == 0) {
          ta <- ta_prev
          message("Filter 'conserved_in_mmusculus' removed all rows; reverting to previous dataset.")
          showNotification("Filter 'conserved_in_mmusculus' removed all rows; reverting to previous dataset.", type = "warning")
        }
      }
      
      target_annotation <- ta
      
      
      # ----------------------------------- milestone 20 ---------------------------
      print("milestone 20: Calculated secondary and duplex energy of ASO seq")
      
      
      # count the number of pre-mRNA transcripts with a perfect match 
      
      # this part will take some time to run...
      uni_tar = dplyr::select(target_annotation, name, length)%>%
        unique() %>%
        split(.,.$length)
      
      uni_tar_list <- future_lapply(
        X = uni_tar,
        FUN = function(X, HS) {
          dict0 <- PDict(X$name, max.mismatch = 0)
          dict1 <- PDict(X$name, max.mismatch = 1)
          
          # perfect match count
          pm <- vwhichPDict(
            pdict = dict0, subject = HS,
            max.mismatch = 0, min.mismatch = 0
          )
          X$gene_hits_pm <- tabulate(unlist(pm), nbins = nrow(X))
          
          # single mismatch count (no indels)
          mm1 <- vwhichPDict(
            pdict = dict1, subject = HS,
            max.mismatch = 1, min.mismatch = 1
          )
          X$gene_hits_1mm <- tabulate(unlist(mm1), nbins = nrow(X))
          
          X
        },
        HS = HS,
        future.seed = TRUE
      )
      
      uni_tar <-run_step("Searching for off-targets.", 0.40, {bind_rows(uni_tar_list)
      })
      # ----------------------------------- milestone 21 ---------------------------
      print("milestone 21: Matched ASO sequences to potential off-targets (perfect match, one mismatch")
      
      target_annotation <- left_join(target_annotation, uni_tar, by = c("name", "length"))
      
      perform_offt <- TRUE
      
      if (isTRUE(input$perfect_input) && isTRUE(input$mismatch_input)) {
        
        prefilter <- nrow(target_annotation)
        
        pm_max <- perf_range()[2]
        mm_max <- mm_range()[2]
        
        target_annotation_filtered <- target_annotation %>%
          dplyr::filter(
            gene_hits_pm <= pm_max,
            gene_hits_1mm <= mm_max
          )
        
        postfilter <- nrow(target_annotation_filtered)
        removed <- prefilter - postfilter
        
        print(paste0(
          "Filtering Oligo sequences with perfect matches <= ", pm_max,
          " and 1 mismatch <= ", mm_max
        ))
        print(paste0("Rows before filtering: ", prefilter))
        print(paste0("Rows after filtering: ", postfilter))
        print(paste0("Filtering removed ", removed, " possible ASOs."))
        
        if (postfilter > 0) {
          target_annotation <- target_annotation_filtered
        } else {
          perform_offt <- FALSE
          showNotification("Filter 'off-targets' removed all rows; reverting.", type = "warning")
        }
      }
      # ----------------------------------- milestone 22 ---------------------------
      print("milestone 22: Filtered ASOs with too many off targets")
      
      if (isTRUE(input$single_aso_input) && "input_order" %in% names(target_annotation)) {
        target_annotation <- target_annotation[order(target_annotation$input_order), ]
      } else {
        target_annotation <- target_annotation[order(target_annotation$gene_hits_pm, target_annotation$gene_hits_1mm), ]
      }
      
      target_annotation <- head(target_annotation, 1000)
      
      if (isTRUE(perform_offt)) {
        
        ta <- target_annotation %>% 
          dplyr::select(name, length) %>% 
          dplyr::distinct()
        print(nrow(ta))
        
        # here it prints the number of ASOs for which it will look up the off targets
        # i could inlcude a break here, whcih shows the number of ASO candidates adn then the estimated 
        # time the off target search will take with the speed of 1 gene per 3.6 seconds or 1000 genes per hour
        # do i want there to be a notification or a break in code waiting for user input 
        # if user input then they will have to stay at the PC until the notification pops up, which is also not good
        # could also make an option at start which asks the max run time or max numebr of ASOs
        # quesion is if they put fewer ASOs, which ASOs from the filtering should be used for the off target search
        # cant just take the first 100 becuase i think they are either not sorted or they are and then it will jus tbe ASOs targeting th ebegnning of the gene
        
        summary_server <- tryCatch(
          {
            res_list <- future_lapply(
              X = seq_len(nrow(ta)),
              FUN = function(i) {
                
                cat("Worker started i=", i, "\n", file = "future_log.txt", append = TRUE)
                
                seq_i <- ta$name[[i]]
                len_i <- ta$length[[i]]
                
                # Retrieve the GGGenome information with the default setting of a maximum of 2 mismatches.
                df <- all_offt(seq_i, mismatches_allowed = 2)
                
                cat(
                  "Worker finished i=", i, 
                  " nrow=", if (is.data.frame(df)) nrow(df) else NA, "\n",
                  file = "future_log.txt", append = TRUE
                )
                
                if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
                
                df$name   <- seq_i
                df$length <- len_i
                df
              },
              future.seed = TRUE
            )
            
            out <- dplyr::bind_rows(Filter(Negate(is.null), res_list)) %>%
              dplyr::mutate(distance = mismatches + deletions + insertions) %>%
              dplyr::distinct(gene_name, match_string, query_seq, .keep_all = TRUE)
            
            out
          },
          error = function(e) {
            message("GGGenome is currently unavailable. Off-target features are disabled. ", e$message)
            NULL
          }
        )
      }
      # ----------------------------------- milestone 23 ----------------------------
      print("milestone 23: GGGenome searched for all ASO off-targets")
      if (!is.null(summary_server)) {
        tmp <- tempfile(fileext = ".bgz")
        
        # Download GnomAD lof metrics by gene
        download.file(
          url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", 
          destfile = tmp,
          mode = "wb"
        )
        gnomad_df <- read_tsv(tmp, show_col_types = FALSE,
                              col_select = c(gene, transcript, oe_lof))
        unlink(tmp)
        
        # ----------------------------------- milestone 24 ---------------------------
        print("milestone 24: Downloaded GnomAD lof metrics by gene to dataframe")
        
        off_targets_total = left_join(summary_server, gnomad_df, by = c('gene_name' = 'gene'))
        
        # ----------------------------------- milestone 25 ---------------------------
        print("milestone 25:")
        
        dist_counts <- off_targets_total %>%
          group_by(name, distance) %>%
          summarise(n_distance = n(), .groups = "drop") %>%
          complete(name, distance = 0:2, fill = list(n_distance = 0)) %>%
          mutate(distance = paste0("n_distance_", distance)) %>%  # bijv. n_distance_0, n_distance_1
          pivot_wider(
            names_from  = distance,
            values_from = n_distance,
            values_fill = 0
          )
        
        # 2. Per name de laagste oe_lof score bepalen
        oe_lof_min <- off_targets_total %>%
          filter(distance <= 1) %>%
          group_by(name) %>%
          summarise(min_oe_lof = min(oe_lof, na.rm = TRUE), .groups = "drop")
        
        # 3. Deze twee samenvattingen samenvoegen tot één tabel per name
        off_summary <- dist_counts %>%
          left_join(oe_lof_min, by = "name")
        
        # 4. Met left_join koppelen aan target_annotation op oligo_seq = name
        target_annotation <- target_annotation %>%
          left_join(off_summary, by = c("name" = "name"))
        
        # ----------------------------------- milestone 26 ---------------------------
        print("milestone 26:")
        
        target_annotation <- target_annotation %>%
          mutate(
            off_target_score =
              1    * n_distance_0 +
              0.3  * n_distance_1 +
              0.02 * n_distance_2
          )
      }
      
      # ----------------------------------- milestone 27 ---------------------------
      print("milestone 27")
      print(motif_weights)
      
      ## Split correlation motifs from dataframe
      motif_results_cols <- names(motif_scores)
      motif_results_cols <- append(motif_results_cols, c('name', 'oligo_seq'))
      motif_cols <- names(motif_weights)
      
      
      # Motif_results including motif count after filtering
      motif_results_filtered <- target_annotation %>%
        select(all_of(motif_results_cols))
      
      # Target annotation with just the correlation score
      target_annotation <- target_annotation %>%
        select(-all_of(motif_cols))
      
      # ----------------------------------- milestone 28 ---------------------------
      print("milestone28")
      
      # Check if the filtered result is empty
      if (nrow(target_annotation) == 0) {
        # If empty, return the original target_regions
        target_annotation <- unfiltered_total_data
        print("No results after filtering, exporting unfiltered list")
        showNotification(
          "No results after filtering, exporting unfiltered list",
          type = "default",
          duration = NULL,
          # Notification stays until clicked away
          closeButton = TRUE
        ) # Include a close button)
      }
      
      
      # Render the tables.
      ################ make download table identical to viewed table
      results1_lookup <- reactive({
        req(target_annotation)
        
        df <- target_annotation
        
        if (isTRUE(input$single_aso_input) && "input_order" %in% names(df)) {
          df <- df %>% dplyr::arrange(input_order)
        } else if ("off_target_score" %in% names(df)) {
          df <- df %>% dplyr::arrange(dplyr::coalesce(off_target_score, Inf))
        }
        
        df
      })
      
      results1_data <- reactive({
        df <- results1_lookup()
        
        column_order <- c(
          if (isTRUE(input$single_aso_input)) "input_order",
          "oligo_seq",
          "name",
          "length",
          "start",
          "end",
          "gc_content",
          "tox_score",
          "off_target_score",
          "n_distance_0",
          "n_distance_1",
          "n_distance_2",
          # "pli_score",
          "sec_energy",
          "duplex_energy",
          "motif_cor_score",
          "min_oe_lof",
          "conserved_in_mmusculus",
          "CGs",
          "chr_start",
          "chr_end",
          "PM_tot_freq",
          "PM_max_freq",
          "PM_count",
          "gene_hits_pm",
          "gene_hits_1mm",
          "NoRepeats",
          "accessibility"
        )
        
        df <- df %>%
          dplyr::select(dplyr::any_of(column_order), dplyr::everything())
        
        column_names <- c(
          input_order            = "Input order",
          oligo_seq              = "ASO sequence",
          name                   = "Target (DNA)",
          length                 = "ASO length (nt)",
          start                  = "Start position in gene",
          end                    = "End position in gene",
          gc_content             = "GC content (%)",
          tox_score              = "Acute neurotox score",
          off_target_score       = "Off-target score",
          n_distance_1           = "Off-targets 1 mismatch (GGGenome)",
          n_distance_2           = "Off-targets 2 mismatches (GGGenome)",
          n_distance_0           = "Perfect matches (GGGenome)",
          sec_energy             = "ASO self-folding energy",
          duplex_energy          = "ASO duplex energy",
          min_oe_lof             = "Min. off-target LoF",
          motif_cor_score        = "Motif correlation score",
          conserved_in_mmusculus = "Conserved in mouse",
          CGs                    = "Number of CpGs",
          chr_start              = "Chromosome start pos.",
          chr_end                = "Chromosome end pos.",
          PM_tot_freq            = "PM total freq.",
          PM_max_freq            = "PM max freq.",
          PM_count               = "PM count",
          gene_hits_pm           = "Perfect matches (Pedersen)",
          gene_hits_1mm          = "Off-targets 1 mismatch (Pedersen)",
          NoRepeats              = "Number of repeats",
          accessibility          = "Accessibility"
        )
        
        to_rename <- intersect(names(df), names(column_names))
        names(df)[match(to_rename, names(df))] <- column_names[to_rename]
        
        if (isTRUE(input$single_aso_input) && "Input order" %in% names(df)) {
          df <- df %>% dplyr::arrange(`Input order`)
        } else if ("Off-target score" %in% names(df)) {
          df <- df %>% dplyr::arrange(dplyr::coalesce(`Off-target score`, Inf))
        }
        
        df
      })
      #############
      
      output$results1 <- DT::renderDT({
        df <- results1_data()
        
        main_cols <- c(
          if (isTRUE(input$single_aso_input)) "Input order",
          "ASO sequence",
          "Target (DNA)",
          "ASO length (nt)",
          "Start position in gene",
          "End position in gene",
          "GC content (%)",
          "Acute neurotox score",
          "Off-target score",
          "Perfect matches (GGGenome)",
          "Off-targets 1 mismatch (GGGenome)",
          "Off-targets 2 mismatches (GGGenome)"
        )
        
        main_cols <- intersect(main_cols, names(df))
        main_idx  <- match(main_cols, names(df))
        extra_idx <- setdiff(seq_along(df), main_idx)
        
        column_defs <- list()
        if (!show_all_cols() && length(extra_idx) > 0) {
          column_defs[[1]] <- list(
            visible = FALSE,
            targets = extra_idx - 1L
          )
        }
        
        dt <- DT::datatable(
          df,
          rownames  = FALSE,
          selection = "single",
          options   = list(
            dom        = "tip",
            columnDefs = column_defs
          )
        )
        
        DT::formatRound(
          dt,
          columns = intersect(c("GC content (%)", "PM total freq.", "PM max freq."), names(df)),
          digits  = 0
        )
      })
      # ----------------------------- Offtarget analysis ---------------------------
      
      # Off-target information is stored in reactiveVals, so dependencies can be easily updated.
      current_seq <- reactiveVal(NULL)
      current_mismatch <- reactiveVal(2)
      current_offtargets <- reactiveVal(NULL)
      cached_results <- reactiveVal(list())
      
      # This function only runs when the application is running on Linux. 
      # It calculates the accessibility of off-targets with a distance of less than 2.
      # The 80-nt snippet is used in the RNAplfold_R function. 
      # The matrix from ViennaRNA is converted to long format to obtain the corresponding accessibility score for the off-target positions in the snippet.
      compute_offtarget_accessibility <- function(df) {
        if (input$linux_input != TRUE) return(df)
        
        for (i in seq_len(nrow(df))) {
          if (df$distance[i] >= 2){
            df$offtarget_accessibility[i] <- NA_real_
            next
          }
          
          offtarget_seq <- gsub("-", "", df$subject_seq[i])
          l_ot <- nchar(offtarget_seq)
          
          snip_result <- acc_snippet(
            begin_target  = df$start_target[i],
            end_target    = df$end_target[i],
            begin_snippet = df$snippet_start[i],
            end_snippet   = df$snippet_end[i],
            full_snippet  = df$snippet[i]
          )
          
          l_snipseq <- nchar(snip_result$snippet_seq)
          
          df$offtarget_accessibility[i] <-
            RNAplfold_R(
              snip_result$snippet_seq,
              u.in = l_ot
            ) %>%
            as_tibble() %>%
            mutate(end = 1:l_snipseq) %>%
            gather(length, accessibility, -end) %>%
            mutate(length = as.integer(length)) %>%
            filter(
              length == l_ot,
              end == snip_result$target_end_internal
            ) %>%
            pull(accessibility)
        }
        
        return(df)
      }
      
      observeEvent(input$apply_mismatch, {
        req(current_seq())
        
        mm  <- as.numeric(input$user_mismatch)
        seq <- toupper(current_seq())
        key <- paste0("mm", mm)
        
        # The selected subset of off-targets is checked to see if it's already stored in R. 
        # If not, the subset is generated by filtering the total off-target data frame. 
        # If the user selects a mismatch of 3, GGGenome is called again, but with a mismatch count of 3. 
        # The subsets are stored in R if they don't already exist, making them easier to retrieve. 
        # If the user is running Linux, the off-target accessibility information is also added to the subset of data frames.
        current_mismatch(mm)
        cache <- cached_results()
        
        if (!is.null(cache[[seq]]) && !is.null(cache[[seq]][[key]])) {
          subset_df <- cache[[seq]][[key]]
          
        } else {
          if (mm %in% c(0, 1, 2)) {
            subset_df <- off_targets_total %>%
              filter(toupper(name) == seq, `distance` <= mm)
            
          } else if (mm == 3) {
            showNotification(
              "Results loading",
              type = "default",
              duration = NULL,
              # Notification stays until clicked away
              closeButton = TRUE
            )
            
            # Retrieve the GGGenome information with the default setting of up to 3 mismatches if the user chooses this option.        new_res <- all_offt(seq, 3)
            new_res$name   <- seq
            new_res$length <- nchar(seq)
            new_res <- new_res %>%
              mutate(distance = mismatches + deletions + insertions,
                     gene_name = str_extract(line, "(?<=\\|)[^;]+")
              ) %>%
              distinct(gene_name, match_string, query_seq, .keep_all = TRUE)
            
            tmp <- tempfile(fileext = ".bgz")
            
            download.file(
              url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", 
              destfile = tmp,
              mode = "wb"
            )
            gnomad_df <- read_tsv(tmp, show_col_types = FALSE,
                                  col_select = c(gene, transcript, oe_lof))
            unlink(tmp)
            
            subset_df = left_join(new_res, gnomad_df, by = c('gene_name' = 'gene'))
          }
          
          subset_df <- compute_offtarget_accessibility(subset_df)
          
          if (is.null(cache[[seq]]))
            cache[[seq]] <- list()
          cache[[seq]][[key]] <- subset_df
          cached_results(cache)
        }
        
        if (input$linux_input == TRUE) {
          subset_df <- subset_df %>% 
            relocate(offtarget_accessibility, .after = insertions)
        }
        current_offtargets(subset_df)
        
        # The number of off-targets, off-target information in data frame form is displayed as output and the output can also be downloaded.
        output$numb_offtargets <- renderText({
          if (is.null(summary_server)) {
            return(NULL)
          }
          paste0("# off targets: ", nrow(subset_df))})
        
      })

      output$offtarget_results <- DT::renderDT({
        req(current_offtargets())
        df <- current_offtargets() %>%
          mutate(
            match_string_display = chartr("ID", "DI", match_string),
            subject_seq_display = purrr::map2_chr(
              subject_seq, match_string_display,
              render_subject_alignment_html
            ),
            off_target_nt_length = nchar(gsub("-", "", subject_seq))
          )
  
        
        df_view <- df %>%
          dplyr::select(-dplyr::any_of(c(
            "line",
            "subject_seq",
            "query_seq",
            "start_target",
            "end_target",
            "snippet",
            "snippet_start",
            "snippet_end",
            "name",
            "matches",
            "length",
            "match_string"
          ))) %>%
          dplyr::select(
            dplyr::any_of(c(
              "gene_name",
              "transcript",
              "match_string_display",
              "subject_seq_display",
              "off_target_nt_length",
              "distance",
              "mismatches",
              "deletions",
              "insertions",
              "offtarget_accessibility",
              "oe_lof"
            )),
            dplyr::everything()
          )
        
        display_map <- c(
          gene_name               = "Gene symbol",
          transcript              = "Transcript (Ensembl)",
          match_string_display    = "Match (alignment)",
          subject_seq_display     = "Off-target sequence",
          off_target_nt_length    = "Length (nt)",
          distance                = "Number of Mismatches/Indels",
          mismatches              = "Mismatches",
          deletions               = "Deletions",
          insertions              = "Insertions",
          offtarget_accessibility = "Accessibility (off-target)",
          oe_lof                  = "GnomAD oe_lof"
        )
        
        nm <- names(df_view)
        nm2 <- ifelse(nm %in% names(display_map), display_map[nm], nm)
        names(df_view) <- make.unique(nm2, sep = " (dup) ")
        
        dist_idx <- which(names(df_view) == "Number of Mismatches/Indels")
        if (length(dist_idx) == 0) {
          dist_idx <- which(grepl("^Number of Mismatches/Indels", names(df_view)))[1]
        }
        distance_col <- dist_idx - 1
        ####&*####
      #   DT::datatable(
      #     df_view,
      #     rownames = FALSE,
      #     options = list(order = list(list(distance_col, "asc")))
      #   )
      # })
        ####&*#### 
        
        dt <- DT::datatable(
          df_view,
          rownames = FALSE,
          escape = FALSE,
          options = list(order = list(list(distance_col, "asc")))
        )
      
        if ("Match (alignment)" %in% names(df_view)) {
          dt <- DT::formatStyle(
            dt,
            columns = "Match (alignment)",
            `font-family` = "monospace",
            `white-space` = "pre",
            `letter-spacing` = "1px"
          )
        }
        
        if ("Off-target sequence" %in% names(df_view)) {
          dt <- DT::formatStyle(
            dt,
            columns = "Off-target sequence",
            `font-family` = "monospace",
            `white-space` = "pre",
            `letter-spacing` = "1px"
          )
        }
        
        dt
      })
      
      output$download_offtarget <- downloadHandler(
        filename = function() {
          paste(
            'offtargets_',
            current_seq(),
            "_mismatches_",
            current_mismatch(),
            "_",
            Sys.Date(),
            '.csv'
          )
        },
        content = function(con) {
          if (is.null(summary_server)) {
            return(NULL)
          }
          data_offtarget <- current_offtargets()
          req(data_offtarget)
          write.csv2(data_offtarget, con, row.names = FALSE)
        }
      )
      
      # ---------------- Proteinatlas/OMIM off-target search -----------------------
      # Function for parsing OMIM and Protein atlast data
      
      PAtlas_cache <- reactiveVal(tibble())
      
      OMIM_cache <- reactiveVal(tibble())
      
      
      observeEvent(input$PAtlas_OMIM_search, {
        offtargets <- current_offtargets()
        req(offtargets)
        
        offtargets_d2 <- offtargets %>%
          filter(distance < 3)
        
        ens_ID <- unique(offtargets_d2$transcript)
        ens_ID <- ens_ID[!is.na(ens_ID) & ens_ID != ""]
        PAtlas_results <- tibble()
        pa_cached <- PAtlas_cache()
        cached_tx <- if (nrow(pa_cached) == 0 || !"transcript" %in% names(pa_cached)) character(0) else unique(pa_cached$transcript)
        new_tx <- setdiff(ens_ID, cached_tx)
        if (length(new_tx) > 0) {
          PAtlas_result_list <- lapply(new_tx, function(tx) {
            PA_parsed <- parse_PAtlas(tx, input$target_tissue)
            if (is.null(PA_parsed) || nrow(PA_parsed) == 0) return(NULL)
            mutate(PA_parsed, transcript = tx)
          })
          PAtlas_results <- bind_rows(Filter(Negate(is.null), PAtlas_result_list))
        }
        
        PAtlas_updated <- bind_rows(pa_cached, PAtlas_results) %>%
          distinct(transcript, .keep_all = TRUE)
        PAtlas_cache(PAtlas_updated)
        
        cols_to_add <- setdiff(names(PAtlas_updated), names(offtargets))
        cols_to_add <- union("transcript", cols_to_add)
        
        offtargets <- offtargets %>%
          left_join(PAtlas_updated %>%
                      select(all_of(cols_to_add)), 
                    by = c("transcript" = "transcript"))
        
        # OMIM, try/catch for when no OMIM key is in file
        tryCatch({
          gene_id <- unique(offtargets_d2$gene_name)
          gene_id <- gene_id[!is.na(gene_id) & gene_id != ""]
          
          omim_cached <- OMIM_cache()
          cached_gene <- if (nrow(omim_cached) == 0 || !"gene" %in% names(omim_cached)) character(0) else unique(omim_cached$gene)
          
          new_genes <- setdiff(gene_id, cached_gene)
          OMIM_results <- tibble()
          if (length(new_genes) > 0) {
            OMIM_results <- lapply(new_genes, OMIM_search) %>%
              bind_rows()}
          
          OMIM_updated <- bind_rows(omim_cached, OMIM_results) %>%
            distinct(gene, .keep_all = TRUE)
          
          OMIM_cache(OMIM_updated)
          
          cols_to_add <- setdiff(names(OMIM_updated), names(offtargets))
          cols_to_add <- union("gene", cols_to_add)
          
          offtargets <- offtargets %>%
            left_join(OMIM_updated %>%
                        select(all_of(cols_to_add)),
                      by = c("gene_name" = "gene"))
        }, error = function(e) {
          message("OMIM search skipped: No OMIM key found in file or OMIM API unavailable")
          offtargets
        })
        
        current_offtargets(offtargets)
        
      })
      
      
      # -------------------------------- RNase H analysis --------------------------
      
      # String reverse function
      reverse_string <- function(x) {
        paste0(rev(strsplit(x, "")[[1]]), collapse = "")
      }
      
      # A stored value for use in the second table and download.
      selected_target <- reactiveVal(NULL)
      oligo_sequence <- reactiveVal(NULL)
      rnaseh_stored <- reactiveVal(NULL)
      
      # Apply end modifications.
      observeEvent(input$add_mods, {
        row_data <- selected_target()
        if (is.null(row_data))
          return()
        
        # The used data collected from the selected row.
        rnaseh_data <- rnaseh_results(
          selected_row_name = row_data$name,
          oligo_seq = row_data$oligo_seq,
          mod_5prime = input$mod_5prime,
          mod_3prime = input$mod_3prime
        )
        
        # Stores the data is usable variable. 
        rnaseh_stored(rnaseh_data)
        
        # Renders table output for rnaseh_results. 
        output$rnaseh_results <- DT::renderDT({
          rnaseh_view <- rename_rnaseh_cols(rnaseh_data)
          DT::datatable(rnaseh_view, selection = list(mode = "single", selected = 1))
        })
        
        # Renders cleavage_visual div on rnaseh page.
        output$cleavage_visual <- renderUI(div())
      })
      
      # The second observer object gives a visual of the cleavage site on the target sequence.
      observeEvent(input$rnaseh_results_rows_selected, {
        row_number <- input$rnaseh_results_rows_selected
        if (length(row_number) == 0)
          return()
        
        row_data <- selected_target()
        if (is.null(row_data))
          return()
        
        rnaseh_data <- rnaseh_stored()
        if (is.null(rnaseh_data))
          return()
        
        # Selected row.
        selected_row <- rnaseh_data[row_number, ]
        
        # Oligo sequence.
        oligo_seq <- row_data$oligo_seq
        
        # Modification inputs to oligo sequence.
        mod5 <- input$mod_5prime
        mod3 <- input$mod_3prime
        
        oligo_len <- nchar(oligo_seq)
        
        # Modified regions of oligo sequence. 
        mod5_region <- substr(oligo_seq, 1, mod5)
        mod3_region <- substr(oligo_seq, oligo_len - mod3 + 1, oligo_len)
        
        # Unmodified mid-section of oligo sequence. 
        mid_region <- substr(oligo_seq, mod5 + 1, oligo_len - mod3)
        
        # Target mRNA sequence.
        position_string <- str_split_1(selected_row$position, " - ")
        start_pos <- as.numeric(position_string[1])
        
        # Getting both orientations of the target mRNA sequence. 
        target_seq_fw <- row_data$name
        target_seq_rv <- reverse_string(target_seq_fw)
        
        rna_len <- nchar(target_seq_fw)
        
        # Cleavage position/coordinates for the visualization. 
        cleavage_start_fw <- start_pos
        cleavage_pos_fw <- cleavage_start_fw + 6
        cleavage_pos_rv <- rna_len - cleavage_pos_fw
        cleavage_site_up <- cleavage_pos_rv - 1
        cleavage_site_down <- cleavage_pos_rv + 7
        
        # Creating all the substrings of the cleavage site for the visualization. 
        upstream <- substr(target_seq_rv, 1, cleavage_site_up - 1)
        site_start <- substr(target_seq_rv, cleavage_site_up, cleavage_pos_rv - 1)
        cut_site <- substr(target_seq_rv, cleavage_pos_rv, cleavage_pos_rv)
        site_down <- substr(target_seq_rv, cleavage_pos_rv + 1, cleavage_site_down)
        downstream <- substr(target_seq_rv, cleavage_site_down + 1, rna_len)
        
        # The oligo visualization with possible modifications in HTML-format that will be used on the RNase H page.
        oligo_visual_fw <- paste0(
          "<b style='color:darkorange;'>5'</b> ",
          "<span style='font-weight:bold; color:#90D5FF;'>",
          mod5_region,
          "</span>",
          mid_region,
          "<span style='font-weight:bold; color:#90D5FF;'>",
          mod3_region,
          "</span>",
          " <b style='color:darkorange;'>3'</b>"
        )
        
        # The target mRNA visualization with cleavage site in HTML-format that will be used on the RNase H page.
        rna_visual_rv <- paste0(
          "<b style='color:darkorange;'>3'</b> ",
          upstream,
          "<span style='background-color: lightblue; color: red; font-weight: bold;'>",
          site_start,
          cut_site,
          "|",
          site_down,
          "</span>",
          downstream,
          " <b style='color:darkorange;'>5'</b>"
        )
        
        # The structured div element with both oligo and target mRNA sequences in HTML-format being rendered. 
        output$cleavage_visual <- renderUI({
          HTML(
            paste0(
              "<h5>Oligo sequence (ASO): </h5>",
              "<div style='font-family: monospace; white-space: pre; font-size: 18px;'>",
              oligo_visual_fw,
              "</div>",
              "<h5>RNA sequence: </h5>",
              "<div style='font-family: monospace; white-space: pre; font-size: 18px;'>",
              rna_visual_rv,
              "</div>"
            )
          )
        })
      })
      
      # Download handler for RNase H results. 
      output$download_rnaseh <- downloadHandler(
        filename = function() {
          row_data <- selected_target()
          if (is.null(row_data)) {
            "RNaseH_results.xlsx"
          } else {
            paste0("RNaseH_results_", row_data$name, ".xlsx")
          }
        },
        content = function(file) {
          data <- rnaseh_stored()
          
          if (is.null(data))
            data <- data.frame(Message = "No avalable data")
          
          data[] <- lapply(data, as.character)
          
          wb <- createWorkbook()
          addWorksheet(wb, "RNaseH_results")
          writeData(wb, "RNaseH_results", data)
          saveWorkbook(wb, file, overwrite = TRUE)
        }
      )
      
      # ----------------------------- Table handlers -------------------------------
      
      # Function for table handler call
      handle_table_events <- function(input,
                                      output,
                                      session,
                                      table_id,
                                      table_data_reactive,
                                      off_targets_total,
                                      selected_target,
                                      oligo_sequence,
                                      rnaseh_stored,
                                      current_seq,
                                      current_mismatch,
                                      current_offtargets,
                                      cached_results) {
        observeEvent({
          list(input[[paste0(table_id, "_cell_clicked")]], input[[paste0(table_id, "_rows_selected")]])
        }, {
          cell <- input[[paste0(table_id, "_cell_clicked")]]
          row  <- input[[paste0(table_id, "_rows_selected")]]
          
          if (!is.null(cell) && !is.null(cell$row)) {
            session$sendCustomMessage("selectRow", list(table = table_id, row   = cell$row))
          }
          
          if (is.null(row) || length(row) == 0) return()
          other_table <- if (table_id == "results1")
            "results2"
          else
            "results1"
          proxy_other <- dataTableProxy(other_table)
          selectRows(proxy_other, NULL)
          
          df_current <- table_data_reactive()
          row_data <- df_current[row, , drop = FALSE]
          
          if (is.null(row_data$name) || nrow(row_data) == 0) return()
          
          seq <- toupper(row_data$name[[1]])
          if (!grepl("^[ACGT]+$", seq)) return()
          
          selected_target(row_data)
          oligo_sequence(row_data$oligo_seq[[1]])
          
          # Off-target functionality
          
          # It is checked whether the clicked line contains the correct information to continue.
          if (is.null(row_data$name)) return()
          seq <- toupper(row_data$name)
          if (!grepl("^[ACGT]+$", seq)) return()
          
          # RNaseH functionality
          selected_target(row_data)
          oligo_sequence(row_data$oligo_seq)
          
          rnaseh_data <- rnaseh_results(
            selected_row_name = row_data$name[[1]],
            oligo_seq = row_data$oligo_seq[[1]],
            mod_5prime = 5,
            mod_3prime = 5
          )
          
          output$rnaseh_title <- renderText(
            paste0("RNase H results for: ", row_data$name[[1]])
          )
          
          output$rnaseh_info <- renderUI({
            HTML(
              paste0(
                "length of sequence: ",
                row_data$length[[1]],
                "<br>",
                "Oligo sequence (ASO): ",
                oligo_sequence()
              )
            )
          })
          
          output$rnaseh_results <- DT::renderDT({
            rnaseh_view <- rename_rnaseh_cols(rnaseh_data)
            DT::datatable(rnaseh_view, selection = list(mode = "single", selected = 1))
          })
          
          updateTabsetPanel(session, "tabs_main_general", selected = "Off target results")
          
          # Off-target functionality
          if (!is.null(summary_server)) {
            
            # Once GGGenome has returned output, the reactiveVals are updated with the correct information for the selected row. 
            # First, R checks whether the correct data frame has already been saved; otherwise, it is created by filtering the off_targets_total data frame and saving it to easily access the correct data frame within the same session. 
            # If Linux is running, the off-target accessibility is also calculated and added to the data frame.
            
            current_seq(seq)
            current_mismatch(2)
            updateSelectInput(session, "user_mismatch", selected = 2)
            
            mm <- current_mismatch()
            
            key <- paste0("mm", mm)
            cache <- cached_results()
            
            if (!is.null(cache[[seq]]) && !is.null(cache[[seq]][[key]])) {
              default_subset <- cache[[seq]][[key]]
              
            } else {
              default_subset <- off_targets_total %>%
                filter(toupper(name) == seq, `distance` <= 2)
            }
            
            if (input$linux_input == TRUE) {
              default_subset <- compute_offtarget_accessibility(default_subset)
              default_subset <- default_subset %>% 
                relocate(offtarget_accessibility, .after = insertions)
            }
            
            # Cache
            cache[[seq]][[key]] <- default_subset
            cached_results(cache)
            
            # results
            default_subset[order(default_subset$distance), ]
            current_offtargets(default_subset)
            
            # The results are displayed as output.
            
            output$offtarget_title <- renderText({
              if (is.null(summary_server)) {
                return(NULL)
              }
              paste0("Off target results for: ", seq)
            })
            output$aso_seq <- renderText({
              if (is.null(summary_server)) {
                return(NULL)
              }
              paste0("ASO sequence: ", as.character(reverseComplement(DNAString(seq))))
            })
            output$numb_offtargets <- renderText({
              if (is.null(summary_server)) {
                return(NULL)
              }
              paste0("# off targets: ", nrow(default_subset))
            })
          }
        })
      }
      
      # Table one handler call
      handle_table_events(
        input = input,
        output = output,
        session = session,
        table_id = "results1",
        table_data = results1_lookup,
        off_targets_total = off_targets_total,
        selected_target = selected_target,
        oligo_sequence = oligo_sequence,
        rnaseh_stored = rnaseh_stored,
        current_seq = current_seq,
        current_mismatch = current_mismatch,
        current_offtargets = current_offtargets,
        cached_results = cached_results
      )
      
      # Table two handler call
      handle_table_events(
        input = input,
        output = output,
        session = session,
        table_id = "results2",
        table_data = nucleobase_select,
        off_targets_total = off_targets_total,
        selected_target = selected_target,
        oligo_sequence = oligo_sequence,
        rnaseh_stored = rnaseh_stored,
        current_seq = current_seq,
        current_mismatch = current_mismatch,
        current_offtargets = current_offtargets,
        cached_results = cached_results
      )
      
      # ----------------------------- End of script --------------------------------
      # unfiltered data downloader
      output$Download_unfiltered <- downloadHandler(
        
        filename = function() {
          paste0(
            "unfiltered_data_results_",
            format(Sys.time(), "%Y%m%d_%H%M%S"),
            ".csv"
          )
        },
        
        content = function(file) {
          write.csv(
            unfiltered_total_data,
            file,
            row.names = FALSE
          )
        }
      )
      
      target_annotation[order(target_annotation$off_target_score), ]
      # Filtered data downloader
      output$Download_filtered <- downloadHandler(
        
        filename = function() {
          paste0(
            "Filtered_data_results_",
            format(Sys.time(), "%Y%m%d_%H%M%S"),
            ".csv"
          )
        },
        
        content = function(file) {
          write.csv(
            results1_data(),
            file,
            row.names = FALSE
          )
        }
      )
    })
    
    end_time <- Sys.time()
    elapsed_str <- format_elapsed(start_time, end_time)
    
    showNotification(
      paste("Run completed in", elapsed_str),
      type = "default",
      duration = NULL,
      closeButton = TRUE
    )
  })
  
  ###############################################################
  ###############################################################
  ###########                Patient tab          ###############
  ###############################################################
  ###############################################################
  
  ######$$$#######
  observeEvent(input$ensemble_id_input_patient, {
    req(input$ensemble_id_input_patient)
    
    ensembl_ID <- input$ensemble_id_input_patient
    
    # use gene symbol annotation toupdate chromsome region
    target_ranges <- gdb_for_choices[names(gdb_for_choices) == ensembl_ID]
    
    if (length(target_ranges) != 1) {
      return()
    }
    
    gene_chr   <- normalize_chr_style(as.character(seqnames(target_ranges)))
    gene_start <- start(target_ranges)
    gene_end   <- end(target_ranges)
    
    region_string <- paste0(gene_chr, ":", gene_start, "-", gene_end)
    
    updateTextInput(
      session,
      inputId = "patient_region_input",
      value = region_string
    )
  }, ignoreInit = TRUE)
  
  observeEvent(input$run_button_patient, {
    
    req(input$ensemble_id_input_patient)
    req(input$patient_variant_file)
    
    withProgress(message = "Preparing patient-specific input...", value = 0, {
      
      incProgress(0.15, detail = "Loading gene annotation")
  
      # load gene annotation database
      txdb_hsa <- tryCatch({
        loadDb("/opt/ERASOR/txdb_hsa_biomart.db")
      }, error = function(e) {
        loadDb("../txdb_hsa_biomart.db")
      })
      
      gdb_hsa <- genes(txdb_hsa)
      seqlevelsStyle(gdb_hsa) <- seqlevelsStyle(Hsapiens)
      
      common_chrs <- intersect(seqlevels(gdb_hsa), seqlevels(Hsapiens))
      gdb_hsa <- keepSeqlevels(gdb_hsa, common_chrs, pruning.mode = "coarse")
      
      # extract chromosome region from gene
      ensembl_ID <- input$ensemble_id_input_patient
      target_ranges <- gdb_hsa[names(gdb_hsa) == ensembl_ID]
      
      if (length(target_ranges) != 1) {
        showNotification(
          "Could not uniquely identify selected gene.",
          type = "error",
          duration = NULL
        )
        return()
      }
      
      gene_chr    <- normalize_chr_style(as.character(seqnames(target_ranges)))
      gene_start  <- start(target_ranges)
      gene_end    <- end(target_ranges)
      gene_strand <- as.character(strand(target_ranges))
      
      incProgress(0.45, detail = "Checking uploaded files")
      
      # variant file and index file are copied into a temp dir
      
      staged <- tryCatch(
        stage_patient_variant_files(
          variant_input = input$patient_variant_file,
          index_input   = input$patient_variant_index_file
        ),
        error = function(e) {
          showNotification(
            paste("File preparation failed:", conditionMessage(e)),
            type = "error",
            duration = NULL
          )
          return(NULL)
        }
      )
      
      if (is.null(staged)) {
        return()
      }
      
      flank <- as.integer(input$patient_flank)
      
      # detect chromsome region
      region_string <- trimws(if (is.null(input$patient_region_input)) "" else input$patient_region_input)
      
      # validate chromosme region format
      if (nzchar(region_string)) {
        if (!is_valid_region_string(region_string)) {
          showNotification(
            "Region must have format chr:start-end",
            type = "error",
            duration = NULL
          )
          return()
        }
      } else {
        # no input means the script will take the gene chromsome region by default
        region_string <- paste0(
          gene_chr, ":",
          max(1L, gene_start - flank), "-",
          gene_end + flank
        )
      }
      # parse the chromosome region string into chr/start/end
      region_info <- parse_region_string(region_string)
      
      # extract patient variants
      incProgress(0.60, detail = "Reading patient variants")
      
      variants_raw <- tryCatch(
        read_patient_variants_bcftools(
          variant_path  = staged$variant_path,
          region_string = region_string,
          is_indexed    = staged$is_indexed
        ),
        error = function(e) {
          showNotification(
            paste("bcftools failed:", conditionMessage(e)),
            type = "error",
            duration = NULL
          )
          return(NULL)
        }
      )
      
      if (is.null(variants_raw)) {
        return()
      }
      
      target_chr <- normalize_chr_style(region_info$chr)
      
      variants_region <- variants_raw %>%
        dplyr::mutate(chr = normalize_chr_style(chr)) %>%
        dplyr::filter(
          chr == target_chr,
          pos >= region_info$start,
          pos <= region_info$end
        ) %>%
        dplyr::arrange(pos)
      
      incProgress(0.70, detail = "Getting reference sequence")
      
      region_gr <- GenomicRanges::GRanges(
        seqnames = region_info$chr,
        ranges   = IRanges::IRanges(
          start = region_info$start,
          end   = region_info$end
        ),
        strand = "+"
      )
      
      # Make chromosome naming match the BSgenome object
      GenomeInfoDb::seqlevelsStyle(region_gr) <- GenomeInfoDb::seqlevelsStyle(Hsapiens)
      
      ref_seq <- getSeq(Hsapiens, region_gr)[[1]]
      
      incProgress(0.80, detail = "Applying Variants to reference")
      
      patient_seq <- apply_snvs_to_reference(
        ref_seq = ref_seq,
        snv_df = variants_region %>%
          dplyr::filter(
            variant_type == "SNV",
            !grepl(",", alt)
          ),
        region_start = region_info$start
      )
      
      incProgress(0.75, detail = "Rendering upload summary")
      
      # generate a summary table for now, this can be replaced once this entire input process works
      patient_summary_df <- tibble(
        Metric = c(
          "Selected gene",
          "Gene chromosome",
          "Gene start",
          "Gene end",
          "Gene strand",
          "Variant file name",
          "Variant file type",
          "Variant file staged path",
          "Index file detected",
          "Index staged path",
          "Indexed input",
          "Region used",
          "Flank",
          "Reference sequence length",
          "Variants in selected region"
        ),
        Value = c(
          ensembl_ID,
          gene_chr,
          gene_start,
          gene_end,
          gene_strand,
          input$patient_variant_file$name,
          staged$variant_type,
          staged$variant_path,
          ifelse(is.null(staged$index_path), "no", "yes"),
          ifelse(is.null(staged$index_path), "", staged$index_path),
          ifelse(staged$is_indexed, "yes", "no"),
          region_string,
          flank,
          nchar(as.character(ref_seq)),
          nrow(variants_region)
        )
      )
      
      
      # render the data tables
      output$patient_summary <- renderTable({
        patient_summary_df
      }, striped = TRUE, bordered = TRUE, spacing = "s")
      
      
      # this output table is jsut to make sure things worked, can also be removed later on
      output$unfiltered_results_table_patient <- renderDT({
        datatable(
          tibble(
            step = c(
              "Gene selected",
              "Variant file uploaded",
              "Variant file staged",
              "Region resolved",
              "Reference sequence retrieved",
              "Patient variants read"
            ),
            status = c(
              "yes",
              "yes",
              "yes",
              "yes",
              "yes",
              ifelse(nrow(variants_region) > 0, "yes", "no variants found")
            )
          ),
          rownames = FALSE,
          options = list(dom = "t", ordering = FALSE),
          class = "compact stripe"
        )
      })
      
      
      # alos just a test table
      output$results1_patient <- renderDT({
        datatable(
          tibble(
            check = c(
              "Gene",
              "Variant file",
              "Index",
              "Region"
            ),
            value = c(
              ensembl_ID,
              input$patient_variant_file$name,
              ifelse(is.null(staged$index_path), "none", basename(staged$index_path)),
              region_string
            )
          ),
          rownames = FALSE,
          options = list(dom = "t", pageLength = 10)
        )
      })
      
      output$patient_snvs_table <- renderDT({
        datatable(
          variants_region %>%
            dplyr::select(
              chr,
              pos,
              ref,
              alt,
              variant_type,
              gt,
              phased,
              allele1_code,
              allele2_code,
              hap1_allele,
              hap2_allele,
              variant_on,
              allele_note,
              sample
            ),
          rownames = FALSE,
          options = list(pageLength = 10, scrollX = TRUE)
        )
      })
      
      output$patient_reference_sequence <- renderText({
        as.character(ref_seq)
      })
      
      # output$patient_modified_sequence <- renderText({
      #   as.character(patient_seq)
      # })
      
      # download buttons for output tables
      output$Download_unfiltered_patient <- downloadHandler(
        filename = function() {
          paste0("patient_input_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
        },
        content = function(file) {
          write.csv(patient_summary_df, file, row.names = FALSE)
        }
      )
      
      output$Download_filtered_patient <- downloadHandler(
        filename = function() {
          paste0("patient_input_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
        },
        content = function(file) {
          write.csv(patient_summary_df, file, row.names = FALSE)
        }
      )
      
      incProgress(1, detail = "Done")
      
      showNotification(
        "Patient input files were uploaded and staged successfully.",
        type = "message",
        duration = 8
      )
    })
  })
#######$$$############
  
  
}

# this code is run a single time, it does not need to be repeated
# it loads the biomart db for gene symbols once, it is then also imported with the docker
# podman struggles in the user version to connect to ensembl for the gene symbols so i rather have it access to it immediately
##################################################################
# library(GenomicFeatures)
# library(GenomeInfoDb)
# library(biomaRt)
# 
# txdb <- loadDb("/opt/ERASOR/txdb_hsa_biomart.db")
# 
# gdb <- genes(txdb)
# 
# chr_to_keep <- c(as.character(1:22), "X", "Y", "MT")
# gdb <- gdb[seqnames(gdb) %in% chr_to_keep]
# 
# ens_ids <- names(gdb)
# 
# mart <- biomaRt::useEnsembl(
#   biomart = "ENSEMBL_MART_ENSEMBL",
#   dataset = "hsapiens_gene_ensembl",
#   host    = "https://www.ensembl.org"
# )
# 
# bm_map <- biomaRt::getBM(
#   attributes = c("ensembl_gene_id", "hgnc_symbol"),
#   filters    = "ensembl_gene_id",
#   values     = ens_ids,
#   mart       = mart
# )
# 
# bm_map <- bm_map[!duplicated(bm_map$ensembl_gene_id), , drop = FALSE]
# 
# saveRDS(bm_map, "gene_symbol_map.rds")
##################################################################

# 1000 genes per hour
# 16 genes per min
# 0.277777 genese per second
# 1 gene per 3.6 seconds



