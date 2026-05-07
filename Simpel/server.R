dev_mode <- TRUE

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

# add docker?
library(R.utils)

#$$$#
library(stringr)
#$$$#

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

# Functions not included in computation ----------------------------------------
## Track time elapsed during the run -------------------------------------------
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

## Renaming RnaseH table columns -----------------------------------------------

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

## calc average and max rnaseh score for main tables
add_rnaseh_summary_columns <- function(df, mod_5prime = 0, mod_3prime = 0) {
  if (is.null(df) || nrow(df) == 0) {
    df$rnaseh_max_score <- numeric()
    df$rnaseh_mean_score <- numeric()
    return(df)
  }
  
  if (!all(c("name", "oligo_seq") %in% names(df))) {
    df$rnaseh_max_score <- NA_real_
    df$rnaseh_mean_score <- NA_real_
    return(df)
  }
  
  rnaseh_key_df <- df %>%
    dplyr::mutate(
      rnaseh_key = paste0(name, "__", oligo_seq)
    ) %>%
    dplyr::distinct(rnaseh_key, name, oligo_seq)
  
  rnaseh_summary <- purrr::pmap_dfr(
    list(
      rnaseh_key_df$rnaseh_key,
      rnaseh_key_df$name,
      rnaseh_key_df$oligo_seq
    ),
    function(rnaseh_key, target_seq, oligo_seq) {
      res <- tryCatch(
        rnaseh_results(
          selected_row_name = target_seq,
          oligo_seq = oligo_seq,
          mod_5prime = mod_5prime,
          mod_3prime = mod_3prime
        ),
        error = function(e) NULL
      )
      
      if (is.null(res) || nrow(res) == 0) {
        return(tibble::tibble(
          rnaseh_key = rnaseh_key,
          rnaseh_max_score = NA_real_,
          rnaseh_mean_score = NA_real_
        ))
      }
      
      score_col <- dplyr::case_when(
        "average" %in% names(res) ~ "average",
        "RNase H score" %in% names(res) ~ "RNase H score",
        TRUE ~ NA_character_
      )
      
      if (is.na(score_col)) {
        return(tibble::tibble(
          rnaseh_key = rnaseh_key,
          rnaseh_max_score = NA_real_,
          rnaseh_mean_score = NA_real_
        ))
      }
      
      scores <- suppressWarnings(as.numeric(res[[score_col]]))
      scores <- scores[is.finite(scores)]
      
      if (length(scores) == 0) {
        return(tibble::tibble(
          rnaseh_key = rnaseh_key,
          rnaseh_max_score = NA_real_,
          rnaseh_mean_score = NA_real_
        ))
      }
      
      tibble::tibble(
        rnaseh_key = rnaseh_key,
        rnaseh_max_score = max(scores, na.rm = TRUE),
        rnaseh_mean_score = mean(scores, na.rm = TRUE)
      )
    }
  )
  
  df %>%
    dplyr::mutate(
      rnaseh_key = paste0(name, "__", oligo_seq)
    ) %>%
    dplyr::left_join(rnaseh_summary, by = "rnaseh_key") %>%
    dplyr::select(-rnaseh_key)
}

## Hint boxes for columns in main table-----------------------------------------
col_with_tooltip <- function(label, tooltip) {
  paste0(
    "<span class='col-header-with-tooltip'>",
    "<span>", label, "</span>",
    "<span class='tooltip-icon'>",
    "<img src='questionmark.png' height='14px'>",
    "<span class='custom-tooltip'>",
    tooltip,
    "</span>",
    "</span>",
    "</span>"
  )
}

strip_header_html <- function(x) {
  x <- gsub("<span class='custom-tooltip'>.*?</span>", "", x, perl = TRUE)
  x <- gsub("<[^>]+>", "", x, perl = TRUE)
  trimws(x)
}

## Analysis mode helpers -------------------------------------------------------

has_reference_gene <- function(x) {
  !is.null(x) && nzchar(trimws(x))
}

build_target_annotation_from_asos <- function(asos) {
  asos <- unique(toupper(trimws(asos)))
  asos <- asos[grepl("^[ACGT]+$", asos)]
  
  if (length(asos) == 0) {
    return(tibble(
      oligo_seq = character(),
      name = character(),
      length = integer(),
      start = integer(),
      end = integer(),
      input_order = integer()
    ))
  }
  
  oligo_dna <- DNAStringSet(asos)
  target_dna <- reverseComplement(oligo_dna)
  
  tibble(
    oligo_seq = as.character(oligo_dna),
    name = as.character(target_dna),   # keep this so the rest of your code still works
    length = nchar(asos),
    start = NA_integer_,
    end = NA_integer_,
    input_order = seq_along(asos)
  )
}

add_empty_reference_columns <- function(df) {
  if (!"chr_start" %in% names(df)) df$chr_start <- NA_integer_
  if (!"chr_end" %in% names(df)) df$chr_end <- NA_integer_
  if (!"region_class" %in% names(df)) df$region_class <- NA_character_
  if (!"target_transcript" %in% names(df)) df$target_transcript <- NA_character_
  if (!"PM_tot_freq" %in% names(df)) df$PM_tot_freq <- NA_real_
  if (!"PM_max_freq" %in% names(df)) df$PM_max_freq <- NA_real_
  if (!"PM_count" %in% names(df)) df$PM_count <- NA_real_
  if (!"NoRepeats" %in% names(df)) df$NoRepeats <- NA_real_
  if (!"conserved_in_mmusculus" %in% names(df)) df$conserved_in_mmusculus <- NA
  df
}

#@@@#
# -------------------------------------------------------------------------
# Specific ASO input mode helpers
# -------------------------------------------------------------------------

single_aso_uses_reference <- function(mode) {
  mode %in% c("with_reference_full", "with_reference_consensus_cf")
}

single_aso_full_reference_mode <- function(mode) {
  identical(mode, "with_reference_full")
}

single_aso_consensus_cf_mode <- function(mode) {
  identical(mode, "with_reference_consensus_cf")
}

locate_specific_asos_in_reference_gene <- function(
    target_annotation,
    RNA_target,
    chr_coord
) {
  if (is.null(target_annotation) || nrow(target_annotation) == 0) {
    return(target_annotation)
  }
  
  if (is.null(RNA_target) || length(RNA_target) != 1 || is.null(chr_coord)) {
    return(add_empty_reference_columns(target_annotation))
  }
  
  out <- vector("list", 0)
  
  for (i in seq_len(nrow(target_annotation))) {
    row_i <- target_annotation[i, , drop = FALSE]
    target_seq <- as.character(row_i$name[[1]])
    
    hits <- Biostrings::matchPattern(
      pattern = Biostrings::DNAString(target_seq),
      subject = RNA_target[[1]],
      fixed = TRUE
    )
    
    if (length(hits) == 0) {
      row_i$start <- NA_integer_
      row_i$end <- NA_integer_
      row_i$chr_start <- NA_integer_
      row_i$chr_end <- NA_integer_
      row_i$reference_gene_match_count <- 0L
      row_i$reference_gene_match_status <- "No perfect target match in selected reference gene"
      
      out[[length(out) + 1L]] <- row_i
      next
    }
    
    for (h in seq_along(hits)) {
      row_h <- row_i
      
      gene_start_local <- BiocGenerics::start(hits)[h]
      gene_end_local   <- BiocGenerics::end(hits)[h]
      
      row_h$start <- as.integer(gene_start_local)
      row_h$end   <- as.integer(gene_end_local)
      
      if (chr_coord$strand == 1L) {
        row_h$chr_start <- chr_coord$start + row_h$start - 1L
        row_h$chr_end   <- chr_coord$start + row_h$end - 1L
      } else {
        row_h$chr_start <- chr_coord$end - row_h$end + 1L
        row_h$chr_end   <- chr_coord$end - row_h$start + 1L
      }
      
      row_h$reference_gene_match_count <- length(hits)
      row_h$reference_gene_match_status <- "Perfect target match in selected reference gene"
      
      out[[length(out) + 1L]] <- row_h
    }
  }
  
  dplyr::bind_rows(out) %>%
    add_empty_reference_columns()
}


build_consensus_with_posmap <- function(
    ref_seq,
    region_chr,
    region_start,
    variants_df,
    consensus_label = "reference"
) {
  ref_seq <- toupper(as.character(ref_seq))
  region_chr <- normalize_chr_style(region_chr)
  
  seq_chars <- strsplit(ref_seq, "", fixed = TRUE)[[1]]
  
  # Each consensus base maps back to a reference genomic coordinate.
  # Inserted bases get NA because they do not exist in the reference genome.
  pos_map <- seq.int(region_start, region_start + length(seq_chars) - 1L)
  
  if (is.null(variants_df) || nrow(variants_df) == 0) {
    return(list(
      label = consensus_label,
      seq = paste0(seq_chars, collapse = ""),
      pos_map = pos_map,
      variants_used = tibble()
    ))
  }
  
  variants_df <- variants_df %>%
    dplyr::mutate(
      chr = normalize_chr_style(chr),
      pos = as.integer(pos),
      ref = toupper(as.character(ref)),
      alt = toupper(as.character(alt))
    ) %>%
    dplyr::filter(
      chr == region_chr,
      pos >= region_start,
      pos <= region_start + length(seq_chars) - 1L,
      grepl("^[ACGT]+$", ref),
      grepl("^[ACGT]+$", alt)
    ) %>%
    dplyr::arrange(pos)
  
  if (nrow(variants_df) == 0) {
    return(list(
      label = consensus_label,
      seq = paste0(seq_chars, collapse = ""),
      pos_map = pos_map,
      variants_used = tibble()
    ))
  }
  
  offset <- 0L
  used_rows <- list()
  
  for (i in seq_len(nrow(variants_df))) {
    v <- variants_df[i, , drop = FALSE]
    
    ref <- v$ref[[1]]
    alt <- v$alt[[1]]
    ref_len <- nchar(ref)
    alt_len <- nchar(alt)
    
    rel_pos <- v$pos[[1]] - region_start + 1L + offset
    
    if (rel_pos < 1L || rel_pos > length(seq_chars)) {
      next
    }
    
    observed_ref <- paste0(
      seq_chars[rel_pos:(rel_pos + ref_len - 1L)],
      collapse = ""
    )
    
    if (!identical(observed_ref, ref)) {
      message(
        "Consensus REF mismatch: ",
        "chr=", v$chr[[1]],
        " pos=", v$pos[[1]],
        " expected=", ref,
        " observed=", observed_ref
      )
      next
    }
    
    left_idx <- if (rel_pos > 1L) seq_len(rel_pos - 1L) else integer()
    
    right_start <- rel_pos + ref_len
    
    right_idx <- if (right_start <= length(seq_chars)) {
      seq.int(right_start, length(seq_chars))
    } else {
      integer()
    }
    
    alt_chars <- strsplit(alt, "", fixed = TRUE)[[1]]
    
    if (alt_len == ref_len) {
      alt_pos_map <- seq.int(v$pos[[1]], v$pos[[1]] + alt_len - 1L)
      
    } else if (alt_len > ref_len) {
      # VCF insertion usually REF=A ALT=ATG.
      # Anchor/ref base maps to genomic coordinate; inserted bases get NA.
      alt_pos_map <- c(
        seq.int(v$pos[[1]], v$pos[[1]] + ref_len - 1L),
        rep(NA_integer_, alt_len - ref_len)
      )[seq_len(alt_len)]
      
    } else {
      # Deletion: only retained ALT bases have reference coordinates.
      alt_pos_map <- seq.int(v$pos[[1]], v$pos[[1]] + alt_len - 1L)
    }
    
    seq_chars <- c(
      seq_chars[left_idx],
      alt_chars,
      seq_chars[right_idx]
    )
    
    pos_map <- c(
      pos_map[left_idx],
      alt_pos_map,
      pos_map[right_idx]
    )
    
    offset <- offset + alt_len - ref_len
    
    used_rows[[length(used_rows) + 1L]] <- v
  }
  
  list(
    label = consensus_label,
    seq = paste0(seq_chars, collapse = ""),
    pos_map = pos_map,
    variants_used = dplyr::bind_rows(used_rows)
  )
}

orient_consensus_to_gene_strand <- function(cons, gene_strand = "+") {
  gene_strand <- as.character(gene_strand)[1]
  
  if (!identical(gene_strand, "-")) {
    return(cons)
  }
  
  cons$seq <- as.character(
    Biostrings::reverseComplement(
      Biostrings::DNAString(cons$seq)
    )
  )
  
  cons$pos_map <- rev(cons$pos_map)
  
  cons
}

build_patient_consensus_set <- function(
    ref_seq,
    region_chr,
    region_start,
    variants_df,
    gene_strand = "+",
    max_ambiguous_combinations = 32
) {
  consensus_list <- list()
  
  make_one_consensus <- function(vars, label) {
    cons <- build_consensus_with_posmap(
      ref_seq = ref_seq,
      region_chr = region_chr,
      region_start = region_start,
      variants_df = vars,
      consensus_label = label
    )
    
    orient_consensus_to_gene_strand(
      cons = cons,
      gene_strand = gene_strand
    )
  }
  
  # Always include reference.
  consensus_list[["reference"]] <- make_one_consensus(
    vars = tibble(),
    label = "reference"
  )
  
  if (is.null(variants_df) || nrow(variants_df) == 0) {
    return(consensus_list)
  }
  
  variants_df <- variants_df %>%
    dplyr::mutate(
      chr = normalize_chr_style(chr),
      pos = as.integer(pos),
      ref = toupper(as.character(ref)),
      alt = toupper(as.character(alt)),
      gt = as.character(gt),
      phased = as.logical(phased),
      variant_on = as.character(variant_on)
    ) %>%
    dplyr::arrange(pos)
  
  homozygous_alt <- variants_df %>%
    dplyr::filter(
      grepl("^([1-9][0-9]*)(/|\\|)\\1$", gt)
    ) %>%
    dplyr::mutate(
      consensus_variant_class = "homozygous_alt"
    )
  
  phased_vars <- variants_df %>%
    dplyr::filter(
      phased == TRUE,
      !grepl("^([1-9][0-9]*)\\|\\1$", gt)
    ) %>%
    dplyr::mutate(
      consensus_variant_class = "phased"
    )
  
  ambiguous_unphased_vars <- variants_df %>%
    dplyr::filter(
      phased == FALSE,
      !grepl("^([1-9][0-9]*)/\\1$", gt)
    ) %>%
    dplyr::mutate(
      consensus_variant_class = "ambiguous_unphased"
    )
  
  # -------------------------------------------------------------
  # Add ONE homozygous ALT consensus if applicable.
  # -------------------------------------------------------------
  if (nrow(homozygous_alt) > 0) {
    consensus_list[["homozygous_alt"]] <- make_one_consensus(
      vars = homozygous_alt,
      label = "homozygous_alt"
    )
  }
  
  # -------------------------------------------------------------
  # Add phased haplotype-specific consensuses only if phased variants exist.
  # Do NOT add homozygous_alt here, because it would duplicate rows.
  # -------------------------------------------------------------
  hap1_vars <- phased_vars %>%
    dplyr::filter(variant_on %in% c("haplotype_1", "both_haplotypes")) %>%
    dplyr::arrange(pos)
  
  hap2_vars <- phased_vars %>%
    dplyr::filter(variant_on %in% c("haplotype_2", "both_haplotypes")) %>%
    dplyr::arrange(pos)
  
  if (nrow(hap1_vars) > 0) {
    consensus_list[["haplotype_1"]] <- make_one_consensus(
      vars = hap1_vars,
      label = "haplotype_1"
    )
  }
  
  if (nrow(hap2_vars) > 0) {
    consensus_list[["haplotype_2"]] <- make_one_consensus(
      vars = hap2_vars,
      label = "haplotype_2"
    )
  }
  
  # -------------------------------------------------------------
  # Ambiguous unphased heterozygous variants.
  #
  # Important:
  # If homozygous ALT variants exist, include them in every ambiguous
  # possibility because they are part of both alleles.
  # This gives:
  #   reference
  #   homozygous_alt
  #   ambiguous_possibility_1
  #   ...
  # but avoids duplicate haplotype rows for 1/1.
  # -------------------------------------------------------------
  if (nrow(ambiguous_unphased_vars) > 0) {
    message(
      "Building capped ambiguous consensus possibilities: ",
      nrow(ambiguous_unphased_vars),
      " ambiguous unphased variants, max ",
      max_ambiguous_combinations,
      " consensus possibilities."
    )
    
    ambiguous_sets <- enumerate_unphased_variant_subsets(
      n = nrow(ambiguous_unphased_vars),
      max_combinations = max_ambiguous_combinations
    )
    
    for (j in seq_along(ambiguous_sets)) {
      idx <- ambiguous_sets[[j]]
      
      vars_j <- ambiguous_unphased_vars[idx, , drop = FALSE]
      
      vars_j <- dplyr::bind_rows(
        homozygous_alt,
        vars_j
      ) %>%
        dplyr::arrange(pos)
      
      consensus_list[[paste0("ambiguous_possibility_", j)]] <-
        make_one_consensus(
          vars = vars_j,
          label = paste0("ambiguous_possibility_", j)
        )
    }
  }
  
  attr(consensus_list, "homozygous_alt_variants") <- homozygous_alt
  attr(consensus_list, "ambiguous_unphased_variants") <- ambiguous_unphased_vars
  
  consensus_list
}

get_reference_seq_at_projected_hit <- function(
    chr,
    chr_start,
    chr_end,
    gene_strand = "+"
) {
  if (is.na(chr_start) || is.na(chr_end)) {
    return(NA_character_)
  }
  
  ref_seq <- tryCatch(
    get_seq_one_based(
      chr = chr,
      start_pos = chr_start,
      end_pos = chr_end
    ),
    error = function(e) NA_character_
  )
  
  if (is.na(ref_seq) || !nzchar(ref_seq)) {
    return(NA_character_)
  }
  
  if (identical(as.character(gene_strand), "-")) {
    ref_seq <- as.character(
      Biostrings::reverseComplement(
        Biostrings::DNAString(ref_seq)
      )
    )
  }
  
  toupper(ref_seq)
}


map_specific_asos_to_consensus_set <- function(
    target_annotation,
    consensus_list,
    gene_start,
    gene_end,
    gene_strand,
    gene_chr = NULL
) {
  if (is.null(target_annotation) || nrow(target_annotation) == 0) {
    return(target_annotation)
  }
  
  out <- vector("list", 0)
  
  for (i in seq_len(nrow(target_annotation))) {
    row_i <- target_annotation[i, , drop = FALSE]
    target_seq <- toupper(as.character(row_i$name[[1]]))
    
    row_hits <- vector("list", 0)
    
    for (consensus_name in names(consensus_list)) {
      cons <- consensus_list[[consensus_name]]
      
      hits <- Biostrings::matchPattern(
        pattern = Biostrings::DNAString(target_seq),
        subject = Biostrings::DNAString(cons$seq),
        fixed = TRUE
      )
      
      if (length(hits) == 0) {
        next
      }
      
      for (h in seq_along(hits)) {
        local_start <- BiocGenerics::start(hits)[h]
        local_end   <- BiocGenerics::end(hits)[h]
        
        hit_pos_map <- cons$pos_map[local_start:local_end]
        ref_positions <- hit_pos_map[!is.na(hit_pos_map)]
        
        row_h <- row_i
        
        row_h$consensus_source <- cons$label
        row_h$consensus_match_start <- local_start
        row_h$consensus_match_end <- local_end
        row_h$consensus_match_status <- "Matched patient/reference consensus"
        
        if (length(ref_positions) > 0) {
          row_h$chr_start <- min(ref_positions)
          row_h$chr_end <- max(ref_positions)
        } else {
          row_h$chr_start <- NA_integer_
          row_h$chr_end <- NA_integer_
        }
        
        if (!is.na(row_h$chr_start) && !is.na(row_h$chr_end)) {
          if (identical(as.character(gene_strand), "+")) {
            row_h$start <- row_h$chr_start - gene_start + 1L
            row_h$end   <- row_h$chr_end - gene_start + 1L
          } else {
            row_h$start <- gene_end - row_h$chr_end + 1L
            row_h$end   <- gene_end - row_h$chr_start + 1L
          }
        } else {
          row_h$start <- NA_integer_
          row_h$end <- NA_integer_
        }
        
        used_vars <- cons$variants_used
        
        if (!is.null(used_vars) &&
            nrow(used_vars) > 0 &&
            !is.na(row_h$chr_start) &&
            !is.na(row_h$chr_end)) {
          
          overlapping_vars <- used_vars %>%
            dplyr::filter(
              pos <= row_h$chr_end,
              pos + nchar(ref) - 1L >= row_h$chr_start
            )
          
          row_h$consensus_variant_overlap_count <- nrow(overlapping_vars)
          row_h$consensus_variant_specific_match <- nrow(overlapping_vars) > 0
          
          row_h$consensus_overlaps_ambiguous_variant <- if (
            nrow(overlapping_vars) > 0 &&
            "consensus_variant_class" %in% names(overlapping_vars)
          ) {
            any(overlapping_vars$consensus_variant_class == "ambiguous_unphased")
          } else {
            FALSE
          }
          
          row_h$consensus_overlaps_homozygous_alt_variant <- if (
            nrow(overlapping_vars) > 0 &&
            "consensus_variant_class" %in% names(overlapping_vars)
          ) {
            any(overlapping_vars$consensus_variant_class == "homozygous_alt")
          } else {
            FALSE
          }
          
          row_h$consensus_variant_notes <- if (nrow(overlapping_vars) > 0) {
            paste(
              purrr::pmap_chr(
                list(
                  overlapping_vars$chr,
                  overlapping_vars$pos,
                  overlapping_vars$ref,
                  overlapping_vars$alt
                ),
                make_hgvs_g
              ),
              collapse = "; "
            )
          } else {
            NA_character_
          }
          
        } else {
          row_h$consensus_variant_overlap_count <- 0L
          row_h$consensus_variant_specific_match <- FALSE
          row_h$consensus_variant_notes <- NA_character_
          row_h$consensus_overlaps_ambiguous_variant <- FALSE
          row_h$consensus_overlaps_homozygous_alt_variant <- FALSE
        }
        
        # -------------------------------------------------------------
        # Important correction:
        # A match to an altered consensus is NOT automatically a match
        # to the unchanged reference sequence at the projected coordinates.
        # For indels, chr_start/chr_end are only projected coordinates.
        # This checks the actual reference sequence at that locus.
        # -------------------------------------------------------------
        if (!is.null(gene_chr) &&
            !is.na(row_h$chr_start) &&
            !is.na(row_h$chr_end)) {
          
          row_h$reference_projected_target <- get_reference_seq_at_projected_hit(
            chr = gene_chr,
            chr_start = row_h$chr_start,
            chr_end = row_h$chr_end,
            gene_strand = gene_strand
          )
          
          row_h$reference_projected_match <- identical(
            toupper(as.character(row_h$reference_projected_target)),
            target_seq
          )
          
        } else {
          row_h$reference_projected_target <- NA_character_
          row_h$reference_projected_match <- NA
        }
        
        row_hits[[length(row_hits) + 1L]] <- row_h
      }
    }
    
    if (length(row_hits) == 0) {
      row_i$consensus_source <- NA_character_
      row_i$consensus_match_start <- NA_integer_
      row_i$consensus_match_end <- NA_integer_
      row_i$consensus_match_status <- "No match in reference/patient consensus"
      row_i$consensus_variant_overlap_count <- 0L
      row_i$consensus_variant_specific_match <- FALSE
      row_i$consensus_variant_notes <- NA_character_
      row_i$chr_start <- NA_integer_
      row_i$chr_end <- NA_integer_
      row_i$start <- NA_integer_
      row_i$end <- NA_integer_
      row_i$reference_projected_target <- NA_character_
      row_i$reference_projected_match <- NA
      row_i$consensus_overlaps_ambiguous_variant <- FALSE
      row_i$consensus_overlaps_homozygous_alt_variant <- FALSE
      
      out[[length(out) + 1L]] <- row_i
      
    } else {
      out <- c(out, row_hits)
    }
  }
  
  res <- dplyr::bind_rows(out)
  
  if (nrow(res) == 0) {
    return(res)
  }
  
  homozygous_alt_variants <- attr(consensus_list, "homozygous_alt_variants")
  
  if (is.null(homozygous_alt_variants)) {
    homozygous_alt_variants <- tibble()
  }
  
  if (nrow(homozygous_alt_variants) > 0) {
    homozygous_alt_variants <- homozygous_alt_variants %>%
      dplyr::mutate(
        chr = normalize_chr_style(chr),
        pos = as.integer(pos),
        ref = toupper(as.character(ref)),
        alt = toupper(as.character(alt)),
        ref_start = pos,
        ref_end = pos + nchar(ref) - 1L
      )
    
    res$overlaps_homozygous_alt_locus <- vapply(
      seq_len(nrow(res)),
      function(i) {
        if (is.na(res$chr_start[i]) || is.na(res$chr_end[i])) {
          return(FALSE)
        }
        
        any(
          homozygous_alt_variants$ref_start <= res$chr_end[i] &
            homozygous_alt_variants$ref_end >= res$chr_start[i]
        )
      },
      logical(1)
    )
  } else {
    res$overlaps_homozygous_alt_locus <- FALSE
  }
  
  # -------------------------------------------------------------
  # Remove misleading duplicate reference rows.
  #
  # If the same input ASO maps to the same projected locus in both:
  #   reference
  #   homozygous_alt / ambiguous / haplotype consensus
  #
  # and the reference sequence at that projected locus does NOT actually
  # equal the ASO target, then the reference row is not a real reference
  # target match and should be removed.
  # -------------------------------------------------------------
  res <- res %>%
    dplyr::group_by(input_order) %>%
    dplyr::mutate(
      input_has_homozygous_alt_match = any(
        consensus_source == "homozygous_alt" &
          consensus_variant_specific_match == TRUE,
        na.rm = TRUE
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      # Remove ambiguous consensus rows unless the ASO actually overlaps
      # an ambiguous/unphased variant.
      !(
        grepl("^ambiguous_possibility_", consensus_source) &
          consensus_overlaps_ambiguous_variant == FALSE
      ),
      
      # Remove reference rows at homozygous ALT loci.
      # In the patient, the reference allele is not present there.
      !(
        consensus_source == "reference" &
          overlaps_homozygous_alt_locus == TRUE &
          input_has_homozygous_alt_match == TRUE
      )
    ) %>%
    dplyr::select(-input_has_homozygous_alt_match)
  
  missing_input_orders <- setdiff(
    unique(target_annotation$input_order),
    unique(res$input_order)
  )
  
  if (length(missing_input_orders) > 0) {
    fallback_rows <- target_annotation %>%
      dplyr::filter(input_order %in% missing_input_orders) %>%
      dplyr::mutate(
        consensus_source = NA_character_,
        consensus_match_start = NA_integer_,
        consensus_match_end = NA_integer_,
        consensus_match_status = "No retained patient/reference consensus match after filtering",
        consensus_variant_overlap_count = 0L,
        consensus_variant_specific_match = FALSE,
        consensus_variant_notes = NA_character_,
        consensus_overlaps_ambiguous_variant = FALSE,
        consensus_overlaps_homozygous_alt_variant = FALSE,
        overlaps_homozygous_alt_locus = FALSE,
        reference_projected_target = NA_character_,
        reference_projected_match = NA
      )
    
    res <- dplyr::bind_rows(res, fallback_rows)
  }
  
  res <- res %>%
    dplyr::mutate(
      consensus_priority = dplyr::case_when(
        consensus_variant_specific_match == TRUE ~ 1L,
        consensus_source == "reference" & reference_projected_match == TRUE ~ 2L,
        consensus_source == "reference" ~ 3L,
        TRUE ~ 4L
      )
    ) %>%
    dplyr::arrange(
      input_order,
      start,
      end,
      consensus_priority
    ) %>%
    dplyr::group_by(
      input_order,
      name,
      start,
      end,
      chr_start,
      chr_end,
      consensus_source
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-consensus_priority)
  
  res %>%
    dplyr::group_by(input_order) %>%
    dplyr::mutate(
      consensus_match_count = sum(
        consensus_match_status == "Matched patient/reference consensus",
        na.rm = TRUE
      )
    ) %>%
    dplyr::ungroup()
}


add_specific_aso_reference_cf_annotations <- function(
    target_annotation,
    txdb_hsa,
    ensembl_ID,
    martHS,
    martMM,
    input_polymorphism = TRUE
) {
  target_annotation <- add_empty_reference_columns(target_annotation)
  
  target_annotation <- annotate_general_knockdown_region(
    target_annotation = target_annotation,
    txdb = txdb_hsa,
    ensembl_ID = ensembl_ID
  )
  
  # Polymorphism frequency
  PMs <- tibble(chr_start = integer(), PM_freq = numeric())
  
  if (isTRUE(input_polymorphism) && !is.null(martHS)) {
    PMs <- tryCatch(
      getBM(
        attributes = c("minor_allele_freq", "chromosome_start"),
        filters = "ensembl_gene_id",
        values = ensembl_ID,
        mart = martHS
      ) %>%
        as_tibble() %>%
        arrange(chromosome_start, desc(minor_allele_freq)) %>%
        filter(
          !is.na(minor_allele_freq),
          !duplicated(chromosome_start)
        ) %>%
        rename(
          chr_start = chromosome_start,
          PM_freq = minor_allele_freq
        ),
      error = function(e) tibble(
        chr_start = integer(),
        PM_freq = numeric()
      )
    )
  }
  
  if (isTRUE(input_polymorphism) && nrow(PMs) > 0) {
    PM_freq <- PMs %>%
      mutate(name = map(chr_start, function(X) {
        filter(
          target_annotation,
          !is.na(chr_start),
          !is.na(chr_end),
          chr_start <= X,
          chr_end >= X
        ) %>%
          select(name, chr_start_anno = chr_start)
      })) %>%
      unnest_legacy() %>%
      rename(chr_start_PM = chr_start) %>%
      group_by(name, chr_start_anno) %>%
      summarise(
        PM_tot_freq = 1 - prod(1 - PM_freq),
        PM_max_freq = max(PM_freq),
        PM_count = n(),
        .groups = "drop"
      )
    
    target_annotation <- left_join(
      target_annotation,
      PM_freq,
      by = c("name" = "name", "chr_start" = "chr_start_anno")
    ) %>%
      mutate(
        PM_tot_freq = coalesce(PM_tot_freq.y, PM_tot_freq.x),
        PM_max_freq = coalesce(PM_max_freq.y, PM_max_freq.x),
        PM_count    = coalesce(PM_count.y, PM_count.x)
      ) %>%
      select(
        -PM_tot_freq.x, -PM_tot_freq.y,
        -PM_max_freq.x, -PM_max_freq.y,
        -PM_count.x, -PM_count.y
      )
  }
  
  # Mouse conservation
  # Mouse conservation ------------------------------------------------------
  
  target_annotation$conserved_in_mmusculus <- NA
  target_annotation$mouse_conservation_status <- NA_character_
  target_annotation$mouse_ortholog_ids <- NA_character_
  
  if (is.null(martHS) || is.null(martMM)) {
    
    target_annotation$mouse_conservation_status <-
      "BioMart unavailable; mouse conservation not calculated"
    
  } else {
    
    ortho_ENS <- tryCatch(
      {
        biomaRt::getBM(
          attributes = "mmusculus_homolog_ensembl_gene",
          filters = "ensembl_gene_id",
          values = ensembl_ID,
          mart = martHS,
          bmHeader = FALSE
        )
      },
      error = function(e) {
        message("Mouse ortholog BioMart query failed: ", conditionMessage(e))
        data.frame(mmusculus_homolog_ensembl_gene = character())
      }
    )
    
    mouse_ids <- unique(ortho_ENS$mmusculus_homolog_ensembl_gene)
    mouse_ids <- mouse_ids[!is.na(mouse_ids) & nzchar(mouse_ids)]
    
    target_annotation$mouse_ortholog_ids <- if (length(mouse_ids) > 0) {
      paste(mouse_ids, collapse = "; ")
    } else {
      NA_character_
    }
    
    if (length(mouse_ids) == 0) {
      
      target_annotation$mouse_conservation_status <-
        "No mouse ortholog Ensembl ID returned by BioMart"
      
    } else {
      
      mouse_seq_df <- tryCatch(
        {
          biomaRt::getBM(
            attributes = c("ensembl_gene_id", "gene_exon_intron"),
            filters = "ensembl_gene_id",
            values = mouse_ids,
            mart = martMM
          )
        },
        error = function(e) {
          message("Mouse gene sequence BioMart query failed: ", conditionMessage(e))
          data.frame(
            ensembl_gene_id = character(),
            gene_exon_intron = character()
          )
        }
      )
      
      mouse_seq_df <- mouse_seq_df %>%
        dplyr::filter(
          !is.na(gene_exon_intron),
          nzchar(gene_exon_intron)
        )
      
      if (nrow(mouse_seq_df) == 0) {
        
        target_annotation$mouse_conservation_status <-
          "Mouse ortholog ID found, but no mouse gene_exon_intron sequence returned"
        
      } else {
        
        RNA_target_mouse <- Biostrings::DNAStringSet(
          mouse_seq_df$gene_exon_intron
        )
        
        oligo_lengths_current <- sort(unique(target_annotation$length))
        
        mouse_windows <- lapply(seq_along(RNA_target_mouse), function(j) {
          mouse_seq_j <- RNA_target_mouse[[j]]
          lm_j <- Biostrings::width(mouse_seq_j)
          
          one_mouse <- lapply(oligo_lengths_current, function(i) {
            i <- as.integer(i)
            
            if (lm_j < i) {
              return(Biostrings::DNAStringSet())
            }
            
            starts_j <- seq_len(lm_j - i + 1L)
            
            Biostrings::DNAStringSet(
              mouse_seq_j,
              start = starts_j,
              width = i
            )
          })
          
          do.call(c, one_mouse)
        })
        
        RNAsitesMM <- do.call(c, mouse_windows)
        
        if (length(RNAsitesMM) == 0) {
          
          target_annotation$conserved_in_mmusculus <- NA
          target_annotation$mouse_conservation_status <-
            "Mouse ortholog sequence found, but no mouse windows generated"
          
        } else {
          
          target_annotation$conserved_in_mmusculus <-
            target_annotation$name %in% as.character(RNAsitesMM)
          
          target_annotation$mouse_conservation_status <-
            paste0(
              "Mouse conservation calculated using ",
              length(mouse_ids),
              " mouse ortholog ID(s)"
            )
        }
      }
    }
  }
  
  target_annotation
}
#@@@#





reference_required_display_cols <- function() {
  c(
    "Start position in gene",
    "End position in gene",
    "Target region",
    "Target transcript(s)",
    "Chromosome start pos.",
    "Chromosome end pos.",
    "PM total freq.",
    "PM max freq.",
    "PM count",
    "Number of repeats",
    "Conserved in mouse",
    "Accessibility",
    "Perfect matches (Pedersen)",
    "Off-targets 1 mismatch (Pedersen)"
  )
}

## Creating off-target sequences with annotated indels/mismatches --------------
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

# target region and tracnript
annotate_general_knockdown_region <- function(target_annotation, txdb, ensembl_ID) {
  if (nrow(target_annotation) == 0) {
    return(target_annotation)
  }
  
  # ---------------- selected gene transcripts ----------------
  tx_df <- as.data.frame(
    transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
  )
  
  tx_df$gene_id_chr <- as.character(tx_df$gene_id)
  tx_df <- tx_df[tx_df$gene_id_chr == ensembl_ID, , drop = FALSE]
  
  if (nrow(tx_df) == 0) {
    target_annotation$region_class <- NA_character_
    target_annotation$target_transcript <- NA_character_
    return(target_annotation)
  }
  
  tx_df$tx_id_chr <- as.character(tx_df$tx_id)
  tx_df$tx_label <- ifelse(
    is.na(tx_df$tx_name) | tx_df$tx_name == "",
    tx_df$tx_id_chr,
    as.character(tx_df$tx_name)
  )
  
  tx_ids <- unique(tx_df$tx_id_chr)
  gene_chr <- as.character(tx_df$seqnames[1])
  gene_strand <- as.character(tx_df$strand[1])
  
  # ---------------- exons from selected gene transcripts ----------------
  ex_df <- as.data.frame(exons(txdb, columns = c("tx_id", "exon_rank")))
  
  if (nrow(ex_df) > 0) {
    ex_df$tx_id_chr <- as.character(ex_df$tx_id)
    ex_df <- ex_df[ex_df$tx_id_chr %in% tx_ids, , drop = FALSE]
    ex_df$tx_label <- tx_df$tx_label[match(ex_df$tx_id_chr, tx_df$tx_id_chr)]
  }
  
  # ---------------- CDS from selected gene transcripts ----------------
  cds_df <- as.data.frame(cds(txdb, columns = c("tx_id")))
  
  if (nrow(cds_df) > 0) {
    cds_df$tx_id_chr <- as.character(cds_df$tx_id)
    cds_df <- cds_df[cds_df$tx_id_chr %in% tx_ids, , drop = FALSE]
    cds_df$tx_label <- tx_df$tx_label[match(cds_df$tx_id_chr, tx_df$tx_id_chr)]
  }
  
  # ---------------- derive UTR intervals from exon minus CDS ----------------
  make_utr_df <- function(ex_df, cds_df, tx_df, strand_value) {
    if (nrow(ex_df) == 0 || nrow(cds_df) == 0) {
      return(data.frame())
    }
    
    out <- vector("list", 0)
    
    for (txid in unique(ex_df$tx_id_chr)) {
      ex_tx <- ex_df[ex_df$tx_id_chr == txid, , drop = FALSE]
      cds_tx <- cds_df[cds_df$tx_id_chr == txid, , drop = FALSE]
      tx_row <- tx_df[tx_df$tx_id_chr == txid, , drop = FALSE]
      
      if (nrow(ex_tx) == 0 || nrow(cds_tx) == 0 || nrow(tx_row) == 0) next
      
      tx_label <- tx_row$tx_label[1]
      cds_min <- min(cds_tx$start)
      cds_max <- max(cds_tx$end)
      
      for (k in seq_len(nrow(ex_tx))) {
        ex_start <- ex_tx$start[k]
        ex_end <- ex_tx$end[k]
        
        # left genomic side of CDS
        if (ex_start < cds_min) {
          s <- ex_start
          e <- min(ex_end, cds_min - 1)
          if (s <= e) {
            part <- if (strand_value == "+") "fiveUTR" else "threeUTR"
            out[[length(out) + 1]] <- data.frame(
              seqnames = as.character(ex_tx$seqnames[k]),
              start = s,
              end = e,
              strand = as.character(ex_tx$strand[k]),
              tx_id_chr = txid,
              tx_label = tx_label,
              utr_part = part,
              stringsAsFactors = FALSE
            )
          }
        }
        
        # right genomic side of CDS
        if (ex_end > cds_max) {
          s <- max(ex_start, cds_max + 1)
          e <- ex_end
          if (s <= e) {
            part <- if (strand_value == "+") "threeUTR" else "fiveUTR"
            out[[length(out) + 1]] <- data.frame(
              seqnames = as.character(ex_tx$seqnames[k]),
              start = s,
              end = e,
              strand = as.character(ex_tx$strand[k]),
              tx_id_chr = txid,
              tx_label = tx_label,
              utr_part = part,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    if (length(out) == 0) {
      return(data.frame())
    }
    
    dplyr::bind_rows(out)
  }
  
  utr_df <- make_utr_df(ex_df, cds_df, tx_df, gene_strand)
  
  # ---------------- annotate each ASO target site ----------------
  target_annotation$region_class <- NA_character_
  target_annotation$target_transcript <- NA_character_
  
  for (i in seq_len(nrow(target_annotation))) {
    q_start <- target_annotation$chr_start[i]
    q_end   <- target_annotation$chr_end[i]
    
    # 1) classify site in the selected gene:
    # full containment in UTR > full containment in exon > otherwise intron
    
    in_utr <- FALSE
    if (nrow(utr_df) > 0) {
      utr_hit <- utr_df[
        utr_df$seqnames == gene_chr &
          utr_df$start <= q_start &
          utr_df$end >= q_end,
        ,
        drop = FALSE
      ]
      in_utr <- nrow(utr_hit) > 0
    }
    
    in_exon <- FALSE
    if (nrow(ex_df) > 0) {
      ex_hit <- ex_df[
        as.character(ex_df$seqnames) == gene_chr &
          ex_df$start <= q_start &
          ex_df$end >= q_end,
        ,
        drop = FALSE
      ]
      in_exon <- nrow(ex_hit) > 0
    }
    
    if (in_utr) {
      target_annotation$region_class[i] <- "UTR"
    } else if (in_exon) {
      target_annotation$region_class[i] <- "Exonic"
    } else {
      target_annotation$region_class[i] <- "Intronic"
    }
    
    # 2) transcript presence:
    # only transcripts of THIS gene where the full target site is inside one exon
    tx_present <- character(0)
    
    if (nrow(ex_df) > 0) {
      ex_hit_tx <- ex_df[
        as.character(ex_df$seqnames) == gene_chr &
          ex_df$start <= q_start &
          ex_df$end >= q_end,
        ,
        drop = FALSE
      ]
      
      if (nrow(ex_hit_tx) > 0) {
        tx_present <- unique(ex_hit_tx$tx_label)
      }
    }
    
    target_annotation$target_transcript[i] <- if (length(tx_present) > 0) {
      paste(tx_present, collapse = "; ")
    } else {
      NA_character_
    }
  }
  
  target_annotation
}

# helper function for vcf and bcf file input for patient design
####$$$####

normalize_chr_style <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  ifelse(x %in% c("M", "MT"), "chrM", paste0("chr", x))
}

normalize_chr_style_nochr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  ifelse(x %in% c("M", "MT"), "MT", x)
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

flip_chr_style <- function(chr) {
  chr <- as.character(chr)
  if (grepl("^chr", chr, ignore.case = TRUE)) {
    sub("^chr", "", chr, ignore.case = TRUE)
  } else {
    paste0("chr", chr)
  }
}

flip_region_chr_style <- function(region) {
  x <- parse_region_string(region)
  paste0(flip_chr_style(x$chr), ":", x$start, "-", x$end)
}

empty_patient_variants <- function() {
  tibble(
    chr = character(),
    pos = integer(),
    ref = character(),
    alt = character(),
    gt = character(),
    sample = character(),
    phased = logical(),
    allele1_code = character(),
    allele2_code = character(),
    hap1_allele = character(),
    hap2_allele = character(),
    variant_on = character(),
    allele_note = character(),
    variant_type = character()
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
    return("Multiallelic")
  }
  
  if (grepl("^<", alt)) {
    return("symbolic")
  }
  
  if (nchar(ref) == 1 && nchar(alt) == 1) {
    return("SNV")
  }
  
  if (nchar(ref) < nchar(alt)) {
    return("Insertion")
  }
  
  if (nchar(ref) > nchar(alt)) {
    return("Deletion")
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

read_patient_variants_bcftools <- function(
    variant_path,
    region_string = NULL,
    is_indexed = FALSE,
    mode = c("auto", "region", "all")
) {
  mode <- match.arg(mode)
  bcftools_path <- get_bcftools_path()
  
  out_file <- tempfile(fileext = ".vcf")
  err_file <- tempfile(fileext = ".log")
  
  args <- c("view", "-H")
  
  if (mode == "region" && !is.null(region_string) && nzchar(region_string)) {
    region_flag <- if (isTRUE(is_indexed)) "-r" else "-t"
    args <- c("view", region_flag, region_string, "-H", variant_path)
  } else {
    args <- c("view", "-H", variant_path)
  }
  
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
    return(empty_patient_variants())
  }
  
  lines <- readLines(out_file, warn = FALSE)
  if (length(lines) == 0) {
    return(empty_patient_variants())
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
      chr = as.character(chr),
      variant_type = purrr::map2_chr(ref, alt, classify_variant_type)
    )
}

read_patient_variants_resilient <- function(
    variant_path,
    region_string,
    is_indexed = FALSE
) {
  stopifnot(is_valid_region_string(region_string))
  
  region_info <- parse_region_string(region_string)
  alt_region_string <- flip_region_chr_style(region_string)
  
  # 1. try requested region as-is
  v1 <- tryCatch(
    read_patient_variants_bcftools(
      variant_path  = variant_path,
      region_string = region_string,
      is_indexed    = is_indexed,
      mode          = "region"
    ),
    error = function(e) structure(empty_patient_variants(), error_msg = conditionMessage(e))
  )
  
  if (nrow(v1) > 0) {
    attr(v1, "query_region_used") <- region_string
    attr(v1, "query_mode_used") <- "region_exact"
    return(v1)
  }
  
  # 2. try flipped chr style
  v2 <- tryCatch(
    read_patient_variants_bcftools(
      variant_path  = variant_path,
      region_string = alt_region_string,
      is_indexed    = is_indexed,
      mode          = "region"
    ),
    error = function(e) structure(empty_patient_variants(), error_msg = conditionMessage(e))
  )
  
  if (nrow(v2) > 0) {
    attr(v2, "query_region_used") <- alt_region_string
    attr(v2, "query_mode_used") <- "region_flipped_chr_style"
    return(v2)
  }
  
  # 3. final fallback: read all rows, then filter in R
  vall <- tryCatch(
    read_patient_variants_bcftools(
      variant_path  = variant_path,
      region_string = NULL,
      is_indexed    = is_indexed,
      mode          = "all"
    ),
    error = function(e) {
      stop(paste(
        "No rows found by region query, and full-file fallback also failed:",
        conditionMessage(e)
      ))
    }
  )
  
  target_chr_chr <- normalize_chr_style(region_info$chr)
  target_chr_nochr <- normalize_chr_style_nochr(region_info$chr)
  
  vall_filt <- vall %>%
    mutate(chr_norm_chr = normalize_chr_style(chr),
           chr_norm_nochr = normalize_chr_style_nochr(chr)) %>%
    filter(
      (chr_norm_chr == target_chr_chr | chr_norm_nochr == target_chr_nochr),
      pos >= region_info$start,
      pos <= region_info$end
    ) %>%
    select(-chr_norm_chr, -chr_norm_nochr)
  
  attr(vall_filt, "query_region_used") <- region_string
  attr(vall_filt, "query_mode_used") <- "full_file_fallback_filtered_in_R"
  vall_filt
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


# to get 2 columsn with refernce genome seq and minus strand real geen seq for snv aso
revcomp_seq <- function(x) {
  x <- toupper(as.character(x))
  as.character(reverseComplement(DNAString(x)))
}

make_variant_label <- function(ref, alt) {
  if (nchar(ref) == 1 && nchar(alt) == 1) {
    return("SNV")
  }
  if (nchar(ref) < nchar(alt)) {
    return("Insertion")
  }
  if (nchar(ref) > nchar(alt)) {
    return("Deletion")
  }
  "Complex"
}

render_variant_on_genomic_seq_html <- function(
    target_seq,
    chr_start,
    chr_end,
    variant_pos,
    variant_ref,
    variant_alt
) {
  if (is.na(target_seq) || is.na(chr_start) || is.na(chr_end) ||
      is.na(variant_pos) || is.na(variant_ref) || is.na(variant_alt)) {
    return(as.character(target_seq))
  }
  
  chars <- strsplit(as.character(target_seq), "")[[1]]
  ref_disp <- toupper(as.character(variant_ref))
  alt_disp <- toupper(as.character(variant_alt))
  rel_pos <- variant_pos - chr_start + 1L
  
  if (rel_pos < 1L || rel_pos > length(chars)) {
    return(as.character(target_seq))
  }
  
  if (nchar(ref_disp) == 1 && nchar(alt_disp) == 1) {
    if (!identical(toupper(chars[rel_pos]), ref_disp)) {
      return(as.character(target_seq))
    }
    
    chars[rel_pos] <- paste0(
      "<span style='color:red;font-weight:bold;text-decoration:underline;'>",
      alt_disp,
      "</span>"
    )
    return(paste(chars, collapse = ""))
  }
  
  if (nchar(ref_disp) < nchar(alt_disp)) {
    if (substr(target_seq, rel_pos, rel_pos + nchar(ref_disp) - 1) != ref_disp) {
      return(as.character(target_seq))
    }
    
    inserted_part <- substr(alt_disp, nchar(ref_disp) + 1, nchar(alt_disp))
    green_insert <- paste0(
      "<span style='color:#00cc00;font-weight:bold;text-decoration:underline;'>",
      strsplit(inserted_part, "")[[1]],
      "</span>",
      collapse = ""
    )
    
    chars[rel_pos] <- paste0(chars[rel_pos], green_insert)
    return(paste(chars, collapse = ""))
  }
  
  if (nchar(ref_disp) > nchar(alt_disp)) {
    if (substr(target_seq, rel_pos, rel_pos + nchar(ref_disp) - 1) != ref_disp) {
      return(as.character(target_seq))
    }
    
    deleted_len <- nchar(ref_disp) - nchar(alt_disp)
    minus_html <- paste0(
      rep("<span style='color:black;font-weight:bold;'>-</span>", deleted_len),
      collapse = ""
    )
    
    chars[rel_pos] <- paste0(chars[rel_pos], minus_html)
    return(paste(chars, collapse = ""))
  }
  
  as.character(target_seq)
}

render_variant_on_gene_seq_html <- function(
    target_seq,
    chr_start,
    chr_end,
    variant_pos,
    variant_ref,
    variant_alt,
    gene_strand = "+"
) {
  if (is.na(target_seq) || is.na(chr_start) || is.na(chr_end) ||
      is.na(variant_pos) || is.na(variant_ref) || is.na(variant_alt)) {
    return(as.character(target_seq))
  }
  
  revcomp_seq <- function(x) {
    as.character(reverseComplement(DNAString(toupper(as.character(x)))))
  }
  
  chars <- strsplit(as.character(target_seq), "")[[1]]
  
  if (identical(gene_strand, "+")) {
    rel_pos <- variant_pos - chr_start + 1L
    ref_disp <- toupper(variant_ref)
    alt_disp <- toupper(variant_alt)
  } else {
    rel_pos <- chr_end - variant_pos + 1L
    ref_disp <- revcomp_seq(variant_ref)
    alt_disp <- revcomp_seq(variant_alt)
  }
  
  if (rel_pos < 1L || rel_pos > length(chars)) {
    return(as.character(target_seq))
  }
  
  if (nchar(ref_disp) == 1 && nchar(alt_disp) == 1) {
    if (!identical(toupper(chars[rel_pos]), ref_disp)) {
      return(as.character(target_seq))
    }
    
    chars[rel_pos] <- paste0(
      "<span style='color:red;font-weight:bold;text-decoration:underline;'>",
      alt_disp,
      "</span>"
    )
    return(paste(chars, collapse = ""))
  }
  
  if (nchar(ref_disp) < nchar(alt_disp)) {
    inserted_part <- substr(alt_disp, nchar(ref_disp) + 1, nchar(alt_disp))
    green_insert <- paste0(
      "<span style='color:#00cc00;font-weight:bold;text-decoration:underline;'>",
      strsplit(inserted_part, "")[[1]],
      "</span>",
      collapse = ""
    )
    
    chars[rel_pos] <- paste0(chars[rel_pos], green_insert)
    return(paste(chars, collapse = ""))
  }
  
  if (nchar(ref_disp) > nchar(alt_disp)) {
    deleted_len <- nchar(ref_disp) - nchar(alt_disp)
    minus_html <- paste0(
      rep("<span style='color:black;font-weight:bold;'>-</span>", deleted_len),
      collapse = ""
    )
    
    chars[rel_pos] <- paste0(chars[rel_pos], minus_html)
    return(paste(chars, collapse = ""))
  }
  
  as.character(target_seq)
}


html_base_span <- function(base, color, underline = TRUE) {
  deco <- if (isTRUE(underline)) "text-decoration:underline;" else ""
  paste0(
    "<span style='color:", color, ";font-weight:bold;", deco, "'>",
    base,
    "</span>"
  )
}

color_patient_alt_sequence_html <- function(
    seq,
    source_vec,
    deletion_marker_after = NA_integer_
) {
  bases <- split_dna(seq)
  
  if (length(bases) == 0 || length(source_vec) != length(bases)) {
    return(as.character(seq))
  }
  
  out <- character(0)
  
  for (i in seq_along(bases)) {
    src <- source_vec[[i]]
    
    rendered_base <- dplyr::case_when(
      src == "inserted" ~ html_base_span(bases[[i]], "#00cc00"),
      src %in% c("SNV", "Complex") ~ html_base_span(bases[[i]], "red"),
      TRUE ~ bases[[i]]
    )
    
    out <- c(out, rendered_base)
    
    if (!is.na(deletion_marker_after) && i == deletion_marker_after) {
      out <- c(out, html_base_span("-", "black"))
    }
  }
  
  paste0(out, collapse = "")
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


#$$$# 
# ths function could also jus tbe inthe patient tba code but if i want funcitons.R ill put it there
# creates the addiitonal columns for patient ASOs 

is_supported_patient_variant <- function(ref, alt) {
  !is.na(ref) &&
    !is.na(alt) &&
    grepl("^[ACGT]+$", toupper(ref)) &&
    grepl("^[ACGT]+$", toupper(alt))
}

label_patient_phasing_status <- function(gt, phased, variant_on) {
  gt <- as.character(gt)
  variant_on <- as.character(variant_on)
  
  if (isTRUE(phased)) {
    if (variant_on == "haplotype_1") return("Haplotype 1")
    if (variant_on == "haplotype_2") return("Haplotype 2")
    if (variant_on == "both_haplotypes") return("Both haplotypes")
    if (variant_on == "reference_only") return("Reference")
    return("Phased")
  }
  
  if (!is.na(gt) && grepl("^([1-9][0-9]*)/\\1$", gt)) return("Unphased homozygous alt")
  if (!is.na(gt) && grepl("^[0-9]+/[0-9]+$", gt)) return("Unphased")
  
  "Unphased"
}

build_patient_variant_aso_rows <- function(
    target_annotation_ref,
    variants_region,
    gene_strand,
    gene_chr,
    gene_start,
    gene_end,
    aso_lengths
) {
  resolved_rows <- list()
  ambiguous_rows <- list()
  
  if (is.null(variants_region) || nrow(variants_region) == 0) {
    return(list(
      resolved = tibble(),
      ambiguous = tibble()
    ))
  }
  
  vars <- variants_region %>%
    dplyr::filter(purrr::map2_lgl(ref, alt, is_supported_patient_variant)) %>%
    dplyr::mutate(
      pos = as.integer(pos),
      ref = toupper(ref),
      alt = toupper(alt),
      chr = normalize_chr_style(chr)
    ) %>%
    dplyr::arrange(pos)
  
  if (nrow(vars) == 0) {
    return(list(
      resolved = tibble(),
      ambiguous = tibble()
    ))
  }
  
  max_L <- max(as.integer(aso_lengths))
  
  # Conservative first implementation:
  # Generate ASOs from each single ALT allele independently.
  # Do not combine nearby unphased variants into one resolved ASO.
  for (i in seq_len(nrow(vars))) {
    v <- vars[i, , drop = FALSE]
    
    local_tbl <- build_single_variant_local_alt_table(
      chr = v$chr[[1]],
      variant_pos = v$pos[[1]],
      ref = v$ref[[1]],
      alt = v$alt[[1]],
      flank_left = max_L - 1L,
      flank_right = max_L - 1L
    )
    
    rows_i <- enumerate_single_variant_asos(
      local_tbl = local_tbl,
      aso_lengths = aso_lengths,
      gene_strand = gene_strand,
      variant_row = v,
      gene_start = gene_start,
      gene_end = gene_end,
      min_variant_bases = 1L
    )
    
    if (!is.null(rows_i) && nrow(rows_i) > 0) {
      resolved_rows[[length(resolved_rows) + 1L]] <- rows_i
    }
  }
  
  resolved <- dplyr::bind_rows(resolved_rows)
  
  if (nrow(resolved) > 0) {
    resolved <- resolved %>%
      dplyr::distinct(
        oligo_seq,
        variant_note,
        local_start,
        local_end,
        .keep_all = TRUE
      )
  }
  
  # Optional warning/ambiguous table:
  # Mark windows where multiple variants are close enough that a single ASO
  # could theoretically cover more than one variant, but phasing may be unknown.
  if (nrow(vars) > 1) {
    for (i in seq_len(nrow(vars))) {
      nearby <- vars %>%
        dplyr::filter(
          pos != vars$pos[[i]],
          abs(pos - vars$pos[[i]]) <= max_L
        )
      
      if (nrow(nearby) > 0) {
        this_group <- dplyr::bind_rows(vars[i, , drop = FALSE], nearby) %>%
          dplyr::arrange(pos)
        
        all_phased <- all(this_group$phased)
        
        if (!all_phased) {
          ambiguous_rows[[length(ambiguous_rows) + 1L]] <- tibble(
            allele_phasing_status = "Ambiguous",
            phase_status = ifelse(all(this_group$phased), "Phased", "Unphased"),
            genotype = paste(unique(this_group$gt), collapse = "; "),
            variant_overlap_count = nrow(this_group),
            variant_note = paste(
              vapply(
                seq_len(nrow(this_group)),
                function(j) make_hgvs_g(
                  this_group$chr[[j]],
                  this_group$pos[[j]],
                  this_group$ref[[j]],
                  this_group$alt[[j]]
                ),
                character(1)
              ),
              collapse = "; "
            ),
            ambiguity_reason = "Nearby unphased variants could occur in the same ASO window; not combined into one resolved ASO.",
            length = NA_integer_,
            start = NA_integer_,
            end = NA_integer_,
            region_class = NA_character_,
            target_transcript = NA_character_,
            chr_start = min(this_group$pos),
            chr_end = max(this_group$pos)
          )
        }
      }
    }
  }
  
  ambiguous <- dplyr::bind_rows(ambiguous_rows)
  
  if (nrow(ambiguous) > 0) {
    ambiguous <- ambiguous %>%
      dplyr::distinct(
        variant_note,
        ambiguity_reason,
        .keep_all = TRUE
      )
  }
  
  list(
    resolved = resolved,
    ambiguous = ambiguous
  )
}


# this one is fror making the dna target output row show the mismatch
render_variant_target_html <- function(
    target_seq,
    chr_start,
    chr_end,
    variant_pos,
    variant_ref,
    variant_alt,
    gene_strand = "+"
) {
  if (is.na(target_seq) || is.na(chr_start) || is.na(chr_end) ||
      is.na(variant_pos) || is.na(variant_ref) || is.na(variant_alt)) {
    return(as.character(target_seq))
  }
  
  revcomp_seq <- function(x) {
    x <- toupper(as.character(x))
    as.character(reverseComplement(DNAString(x)))
  }
  
  chars <- strsplit(as.character(target_seq), "")[[1]]
  strand_val <- as.character(gene_strand)[1]
  
  if (strand_val == "+") {
    rel_pos <- variant_pos - chr_start + 1L
    ref_disp <- toupper(variant_ref)
    alt_disp <- toupper(variant_alt)
  } else {
    rel_pos <- chr_end - variant_pos + 1L
    ref_disp <- revcomp_seq(variant_ref)
    alt_disp <- revcomp_seq(variant_alt)
  }
  
  if (rel_pos < 1L || rel_pos > length(chars)) {
    return(as.character(target_seq))
  }
  
  # ---------------- SNV ----------------
  if (nchar(ref_disp) == 1 && nchar(alt_disp) == 1) {
    if (!identical(toupper(chars[rel_pos]), alt_disp)) {
      return(as.character(target_seq))
    }
    
    chars[rel_pos] <- paste0(
      "<span style='color:red;font-weight:bold;text-decoration:underline;'>",
      chars[rel_pos],
      "</span>"
    )
    
    return(paste(chars, collapse = ""))
  }
  
  # ---------------- insertion ----------------
  # VCF style: ref=A alt=ATG  => inserted part = TG
  if (nchar(ref_disp) < nchar(alt_disp)) {
    inserted_part <- substr(alt_disp, nchar(ref_disp) + 1, nchar(alt_disp))
    
    green_insert <- paste0(
      "<span style='color:#00cc00;font-weight:bold;text-decoration:underline;'>",
      strsplit(inserted_part, "")[[1]],
      "</span>",
      collapse = ""
    )
    
    # place insertion after anchor base
    chars[rel_pos] <- paste0(chars[rel_pos], green_insert)
    return(paste(chars, collapse = ""))
  }
  
  # ---------------- deletion ----------------
  # VCF style: ref=ATG alt=A => deleted part = TG
  if (nchar(ref_disp) > nchar(alt_disp)) {
    deleted_len <- nchar(ref_disp) - nchar(alt_disp)
    
    minus_html <- paste0(
      rep("<span style='color:black;font-weight:bold;'>-</span>", deleted_len),
      collapse = ""
    )
    
    # place minus signs after anchor base
    chars[rel_pos] <- paste0(chars[rel_pos], minus_html)
    return(paste(chars, collapse = ""))
  }
  
  # fallback for complex substitutions
  as.character(target_seq)
}


# this is a but confusing becuase of variants on the minus strand, those are ereverse compleemtns so that vairant not actuallyneeeds to be changed

revcomp_chr <- function(x) {
  as.character(reverseComplement(DNAStringSet(as.character(x))))
}

get_genomic_window_seq <- function(chr, start_pos, end_pos) {
  gr <- GenomicRanges::GRanges(
    seqnames = normalize_chr_style(chr),
    ranges = IRanges::IRanges(start = start_pos, end = end_pos),
    strand = "+"
  )
  GenomeInfoDb::seqlevelsStyle(gr) <- GenomeInfoDb::seqlevelsStyle(Hsapiens)
  as.character(getSeq(Hsapiens, gr))
}



# this is for indels 
normalize_chr_for_bsgenome <- function(chr) {
  chr <- as.character(chr)
  hs_levels <- seqlevels(Hsapiens)
  
  if (any(grepl("^chr", hs_levels))) {
    normalize_chr_style(chr)
  } else {
    normalize_chr_style_nochr(chr)
  }
}

get_seq_one_based <- function(chr, start_pos, end_pos) {
  if (start_pos > end_pos) return("")
  
  gr <- GenomicRanges::GRanges(
    seqnames = normalize_chr_for_bsgenome(chr),
    ranges = IRanges::IRanges(start = start_pos, end = end_pos),
    strand = "+"
  )
  
  as.character(getSeq(Hsapiens, gr))
}


### this is all for indels only
split_dna <- function(x) {
  x <- toupper(as.character(x))
  
  if (is.na(x) || x == "" || nchar(x) == 0) {
    return(character(0))
  }
  
  strsplit(x, "", fixed = TRUE)[[1]]
}


safe_get_seq_one_based <- function(chr, start_pos, end_pos) {
  start_pos <- as.integer(start_pos)
  end_pos <- as.integer(end_pos)
  
  if (is.na(start_pos) || is.na(end_pos) || start_pos > end_pos) {
    return("")
  }
  
  get_seq_one_based(chr, start_pos, end_pos)
}


build_single_variant_local_alt_table <- function(
    chr,
    variant_pos,
    ref,
    alt,
    flank_left,
    flank_right
) {
  chr <- normalize_chr_style(chr)
  variant_pos <- as.integer(variant_pos)
  ref <- toupper(as.character(ref))
  alt <- toupper(as.character(alt))
  
  ref_len <- nchar(ref)
  alt_len <- nchar(alt)
  
  left_start <- max(1L, variant_pos - flank_left)
  left_end   <- variant_pos - 1L
  
  right_start <- variant_pos + ref_len
  right_end   <- variant_pos + ref_len + flank_right - 1L
  
  left_seq <- safe_get_seq_one_based(chr, left_start, left_end)
  right_seq <- safe_get_seq_one_based(chr, right_start, right_end)
  
  observed_ref <- safe_get_seq_one_based(
    chr,
    variant_pos,
    variant_pos + ref_len - 1L
  )
  
  if (!identical(toupper(observed_ref), ref)) {
    message(
      "REF mismatch in build_single_variant_local_alt_table: ",
      "chr=", chr,
      " pos=", variant_pos,
      " expected_ref=", ref,
      " observed_ref=", observed_ref
    )
    return(tibble())
  }
  
  variant_kind <- dplyr::case_when(
    ref_len == 1L && alt_len == 1L ~ "SNV",
    alt_len > ref_len ~ "Insertion",
    alt_len < ref_len ~ "Deletion",
    TRUE ~ "Complex"
  )
  
  left_bases <- split_dna(left_seq)
  
  left_tbl <- tibble(
    base = left_bases,
    ref_pos = if (length(left_bases) > 0) seq.int(left_start, left_end) else integer(),
    source = "left_flank",
    is_variant_base = FALSE
  )
  
  alt_bases <- split_dna(alt)
  
  if (variant_kind == "Insertion" &&
      substr(ref, 1, 1) == substr(alt, 1, 1)) {
    
    # VCF insertion: REF=A, ALT=A + inserted sequence.
    # The first ALT base is the anchor/reference base.
    # The remaining ALT bases are the inserted patient-specific bases.
    alt_ref_pos <- c(variant_pos, rep(NA_integer_, alt_len - 1L))
    alt_source <- c("anchor", rep("inserted", alt_len - 1L))
    alt_is_variant <- c(FALSE, rep(TRUE, alt_len - 1L))
    
  } else if (variant_kind == "Deletion" &&
             substr(ref, 1, 1) == substr(alt, 1, 1)) {
    
    # VCF deletion: REF=AGGT, ALT=A.
    # ALT contains only the retained anchor base.
    # The anchor base is real sequence and should NOT become "-".
    # We mark it as normal anchor, and later add one visual deletion marker
    # between the anchor and the right flank.
    alt_ref_pos <- seq.int(variant_pos, variant_pos + alt_len - 1L)
    alt_source <- rep("anchor", alt_len)
    alt_is_variant <- rep(TRUE, alt_len)
    
  } else if (ref_len == alt_len) {
    
    alt_ref_pos <- seq.int(variant_pos, variant_pos + alt_len - 1L)
    alt_source <- rep(variant_kind, alt_len)
    alt_is_variant <- alt_bases != split_dna(ref)
    
  } else {
    
    # Complex fallback.
    alt_ref_pos <- c(
      seq.int(variant_pos, variant_pos + min(ref_len, alt_len) - 1L),
      rep(NA_integer_, max(0L, alt_len - ref_len))
    )[seq_len(alt_len)]
    
    alt_source <- rep(variant_kind, alt_len)
    alt_is_variant <- rep(TRUE, alt_len)
  }
  
  alt_tbl <- tibble(
    base = alt_bases,
    ref_pos = alt_ref_pos,
    source = alt_source,
    is_variant_base = alt_is_variant
  )
  
  right_bases <- split_dna(right_seq)
  
  right_tbl <- tibble(
    base = right_bases,
    ref_pos = if (length(right_bases) > 0) seq.int(right_start, right_end) else integer(),
    source = "right_flank",
    is_variant_base = FALSE
  )
  
  dplyr::bind_rows(left_tbl, alt_tbl, right_tbl) %>%
    dplyr::mutate(local_pos = dplyr::row_number())
}


enumerate_single_variant_asos <- function(
    local_tbl,
    aso_lengths,
    gene_strand,
    variant_row,
    gene_start,
    gene_end,
    min_variant_bases = 1L
) {
  if (is.null(local_tbl) || nrow(local_tbl) == 0) {
    return(tibble())
  }
  
  local_seq <- paste0(local_tbl$base, collapse = "")
  
  if (!grepl("^[ACGT]+$", local_seq)) {
    return(tibble())
  }
  
  rows <- list()
  
  for (L in aso_lengths) {
    L <- as.integer(L)
    
    if (nchar(local_seq) < L) next
    
    for (local_start in seq_len(nchar(local_seq) - L + 1L)) {
      local_end <- local_start + L - 1L
      win <- local_tbl[local_start:local_end, , drop = FALSE]
      
      variant_base_count <- sum(win$is_variant_base, na.rm = TRUE)
      
      if (variant_base_count < min_variant_bases) next
      
      target_genomic_seq <- paste0(win$base, collapse = "")
      
      is_deletion <- make_variant_label(
        variant_row$ref[[1]],
        variant_row$alt[[1]]
      ) == "Deletion"
      
      deletion_marker_after_genomic <- NA_integer_
      
      if (is_deletion) {
        # For simple VCF deletions REF=AGGT, ALT=A,
        # put one visual marker after the retained anchor base.
        deletion_anchor_ref_pos <- variant_row$pos[[1]] + nchar(variant_row$alt[[1]]) - 1L
        
        marker_idx <- which(
          !is.na(win$ref_pos) &
            win$ref_pos == deletion_anchor_ref_pos
        )
        
        if (length(marker_idx) > 0) {
          deletion_marker_after_genomic <- marker_idx[[1]]
        }
      }
      
      target_gene_seq <- if (identical(as.character(gene_strand), "+")) {
        target_genomic_seq
      } else {
        as.character(reverseComplement(DNAString(target_genomic_seq)))
      }
      
      deletion_marker_after_gene <- deletion_marker_after_genomic
      
      if (is_deletion &&
          !is.na(deletion_marker_after_genomic) &&
          identical(as.character(gene_strand), "-")) {
        deletion_marker_after_gene <- nchar(target_gene_seq) - deletion_marker_after_genomic
      }
      
      oligo_seq <- as.character(reverseComplement(DNAString(target_gene_seq)))
      
      ref_positions <- win$ref_pos[!is.na(win$ref_pos)]
      
      chr_start <- if (length(ref_positions) > 0) {
        min(ref_positions)
      } else {
        variant_row$pos[[1]]
      }
      
      chr_end <- if (length(ref_positions) > 0) {
        max(ref_positions)
      } else {
        variant_row$pos[[1]]
      }
      
      gene_window_start <- if (identical(as.character(gene_strand), "+")) {
        chr_start - gene_start + 1L
      } else {
        gene_end - chr_end + 1L
      }
      
      gene_window_end <- if (identical(as.character(gene_strand), "+")) {
        chr_end - gene_start + 1L
      } else {
        gene_end - chr_start + 1L
      }
      
      rows[[length(rows) + 1L]] <- tibble(
        patient_specific_aso = "Yes",
        allele_phasing_status = label_patient_phasing_status(
          gt = variant_row$gt[[1]],
          phased = variant_row$phased[[1]],
          variant_on = variant_row$variant_on[[1]]
        ),
        phase_status = ifelse(isTRUE(variant_row$phased[[1]]), "Phased", "Unphased"),
        genotype = variant_row$gt[[1]],
        variant_overlap_count = 1L,
        variant_type_patient = make_variant_label(
          variant_row$ref[[1]],
          variant_row$alt[[1]]
        ),
        variant_note = make_hgvs_g(
          variant_row$chr[[1]],
          variant_row$pos[[1]],
          variant_row$ref[[1]],
          variant_row$alt[[1]]
        ),
        oligo_seq = oligo_seq,
        name = target_gene_seq,
        target_genomic_seq = target_genomic_seq,
        target_gene_seq = target_gene_seq,
        target_genomic_source = paste(win$source, collapse = ";"),
        target_gene_source = if (identical(as.character(gene_strand), "+")) {
          paste(win$source, collapse = ";")
        } else {
          paste(rev(win$source), collapse = ";")
        },
        deletion_marker_after_genomic = deletion_marker_after_genomic,
        deletion_marker_after_gene = deletion_marker_after_gene,
        length = L,
        local_start = local_start,
        local_end = local_end,
        start = gene_window_start,
        end = gene_window_end,
        chr = normalize_chr_style(variant_row$chr[[1]]),
        chr_start = chr_start,
        chr_end = chr_end,
        variant_pos = variant_row$pos[[1]],
        variant_ref = variant_row$ref[[1]],
        variant_alt = variant_row$alt[[1]],
        variant_pos_all = as.character(variant_row$pos[[1]]),
        variant_ref_all = variant_row$ref[[1]],
        variant_alt_all = variant_row$alt[[1]],
        variant_base_count = variant_base_count
      )
    }
  }
  
  dplyr::bind_rows(rows) %>%
    dplyr::distinct(
      oligo_seq,
      variant_note,
      local_start,
      local_end,
      .keep_all = TRUE
    )
}
###


### debug for overdropping
validate_variants_against_reference <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  obs_ref <- vapply(seq_len(nrow(df)), function(i) {
    get_seq_one_based(
      chr = df$chr[i],
      start_pos = df$pos[i],
      end_pos = df$pos[i] + nchar(df$ref[i]) - 1L
    )
  }, character(1))
  
  df %>%
    dplyr::mutate(
      ref_observed_genome = toupper(obs_ref),
      ref_matches_genome = toupper(ref) == toupper(ref_observed_genome)
    )
}


build_variant_haplotype_window <- function(chr, chr_start, chr_end, variant_pos, ref, alt) {
  ref <- toupper(ref)
  alt <- toupper(alt)
  
  L <- chr_end - chr_start + 1L
  ref_len <- nchar(ref)
  alt_len <- nchar(alt)
  
  # reference window on genomic forward strand
  ref_window <- get_seq_one_based(chr, chr_start, chr_end)
  if (nchar(ref_window) != L) return(NA_character_)
  
  rel_pos <- variant_pos - chr_start + 1L
  if (rel_pos < 1L || rel_pos > L) return(NA_character_)
  
  # VCF anchor check
  ref_from_window <- substr(ref_window, rel_pos, rel_pos + ref_len - 1L)
  if (!identical(ref_from_window, ref)) return(NA_character_)
  
  left_part  <- if (rel_pos > 1L) substr(ref_window, 1L, rel_pos - 1L) else ""
  right_start <- rel_pos + ref_len
  right_part <- if (right_start <= L) substr(ref_window, right_start, L) else ""
  
  hap <- paste0(left_part, alt, right_part)
  
  hap_len <- nchar(hap)
  
  if (hap_len == L) {
    return(hap)
  }
  
  if (hap_len > L) {
    # insertion: trim extra nt from the right side
    return(substr(hap, 1L, L))
  }
  
  # deletion: extend on the right side to restore length L
  need <- L - hap_len
  extra_right <- get_seq_one_based(chr, chr_end + 1L, chr_end + need)
  hap2 <- paste0(hap, extra_right)
  
  if (nchar(hap2) != L) return(NA_character_)
  hap2
}



# for multivaraint windows
build_multi_variant_haplotype_window <- function(chr, chr_start, chr_end, variants_df) {
  ref_window <- get_seq_one_based(chr, chr_start, chr_end)
  L <- nchar(ref_window)
  
  if (L == 0) {
    message("Empty ref window: chr=", chr, " start=", chr_start, " end=", chr_end)
    return(NA_character_)
  }
  
  variants_df <- variants_df %>% dplyr::arrange(pos)
  
  current_seq <- ref_window
  offset <- 0L
  
  for (i in seq_len(nrow(variants_df))) {
    vpos <- variants_df$pos[i]
    ref <- toupper(variants_df$ref[i])
    alt <- toupper(variants_df$alt[i])
    
    rel_pos <- (vpos - chr_start + 1L) + offset
    
    seq_ref <- substr(current_seq, rel_pos, rel_pos + nchar(ref) - 1L)
    
    if (!identical(seq_ref, ref)) {
      message(
        "Ref mismatch in build_multi_variant_haplotype_window: ",
        "chr=", chr,
        " window=", chr_start, "-", chr_end,
        " pos=", vpos,
        " rel_pos=", rel_pos,
        " expected_ref=", ref,
        " observed_ref=", seq_ref,
        " current_seq=", current_seq
      )
      return(NA_character_)
    }
    
    left  <- if (rel_pos > 1L) substr(current_seq, 1L, rel_pos - 1L) else ""
    right_start <- rel_pos + nchar(ref)
    right <- if (right_start <= nchar(current_seq)) substr(current_seq, right_start, nchar(current_seq)) else ""
    
    current_seq <- paste0(left, alt, right)
    offset <- offset + (nchar(alt) - nchar(ref))
  }
  
  if (nchar(current_seq) > L) {
    current_seq <- substr(current_seq, 1L, L)
  } else if (nchar(current_seq) < L) {
    need <- L - nchar(current_seq)
    extra_right <- get_seq_one_based(chr, chr_end + 1L, chr_end + need)
    current_seq <- paste0(current_seq, extra_right)
  }
  
  if (nchar(current_seq) != L) return(NA_character_)
  current_seq
}



# function to generate the maibguousASO table
enumerate_unphased_variant_subsets <- function(n, max_combinations = 32) {
  n <- as.integer(n)
  max_combinations <- as.integer(max_combinations)
  
  if (is.na(n) || n <= 0L || is.na(max_combinations) || max_combinations <= 0L) {
    return(list())
  }
  
  idx <- seq_len(n)
  out <- vector("list", 0)
  
  # First add single-variant possibilities.
  for (i in idx) {
    out[[length(out) + 1L]] <- i
    
    if (length(out) >= max_combinations) {
      return(out)
    }
  }
  
  # Then add nearby pair possibilities, capped.
  # This avoids trying to build all 2^n possible haplotypes.
  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      for (j in seq.int(i + 1L, n)) {
        out[[length(out) + 1L]] <- c(i, j)
        
        if (length(out) >= max_combinations) {
          return(out)
        }
      }
    }
  }
  
  out
}

#$$$# 

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



# function to help take multiallelic variants and pslit them into two entries
expand_multiallelic_patient_variants <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  out <- vector("list", 0)
  
  for (i in seq_len(nrow(df))) {
    row <- df[i, , drop = FALSE]
    
    alt_parts <- strsplit(as.character(row$alt), ",", fixed = TRUE)[[1]]
    gt <- as.character(row$gt)
    
    # split genotype into allele codes
    phased <- grepl("\\|", gt)
    sep <- if (phased) "\\|" else "/"
    gt_parts <- strsplit(gt, sep)[[1]]
    
    # keep only ALT allele indices actually present in genotype
    gt_nums <- suppressWarnings(as.integer(gt_parts))
    alt_idxs_present <- sort(unique(gt_nums[!is.na(gt_nums) & gt_nums > 0]))
    
    if (length(alt_idxs_present) == 0) next
    
    for (alt_idx in alt_idxs_present) {
      if (alt_idx > length(alt_parts)) next
      
      new_row <- row
      new_row$alt <- alt_parts[alt_idx]
      
      # allele-specific annotation
      new_row$multiallelic_original_alt <- row$alt
      new_row$alt_index <- alt_idx
      
      allele_on <- character(0)
      if (length(gt_parts) >= 1 && gt_parts[1] == as.character(alt_idx)) allele_on <- c(allele_on, "haplotype_1")
      if (length(gt_parts) >= 2 && gt_parts[2] == as.character(alt_idx)) allele_on <- c(allele_on, "haplotype_2")
      
      new_row$variant_on <- dplyr::case_when(
        length(allele_on) == 2 ~ "both_haplotypes",
        length(allele_on) == 1 ~ allele_on[1],
        TRUE ~ "unphased"
      )
      
      new_row$allele1_code <- if (length(gt_parts) >= 1) gt_parts[1] else NA_character_
      new_row$allele2_code <- if (length(gt_parts) >= 2) gt_parts[2] else NA_character_
      
      out[[length(out) + 1]] <- new_row
    }
  }
  
  dplyr::bind_rows(out)
}



# function for HGVS annotation
make_hgvs_g <- function(chr, pos, ref, alt) {
  chr <- as.character(chr)
  pos <- as.integer(pos)
  ref <- toupper(as.character(ref))
  alt <- toupper(as.character(alt))
  
  prefix <- paste0(chr, ":g.")
  
  # SNV
  if (nchar(ref) == 1 && nchar(alt) == 1) {
    return(paste0(prefix, pos, ref, ">", alt))
  }
  
  # deletion
  if (nchar(ref) > nchar(alt) && substr(ref, 1, 1) == substr(alt, 1, 1)) {
    deleted <- substr(ref, nchar(alt) + 1, nchar(ref))
    del_start <- pos + nchar(alt)
    del_end <- pos + nchar(ref) - 1
    
    if (nchar(deleted) == 1) {
      return(paste0(prefix, del_start, "del"))
    } else {
      return(paste0(prefix, del_start, "_", del_end, "del"))
    }
  }
  
  # insertion
  if (nchar(alt) > nchar(ref) && substr(ref, 1, 1) == substr(alt, 1, 1)) {
    inserted <- substr(alt, nchar(ref) + 1, nchar(alt))
    left_base <- pos
    right_base <- pos + 1
    return(paste0(prefix, left_base, "_", right_base, "ins", inserted))
  }
  
  # delins / complex substitution
  if (nchar(ref) > 1) {
    end_pos <- pos + nchar(ref) - 1
    return(paste0(prefix, pos, "_", end_pos, "delins", alt))
  }
  
  paste0(prefix, pos, "delins", alt)
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




#!!!!!!!!!!!!! dev mode
make_mock_aso_results <- function(n = 25) {
  random_dna <- function(len) {
    paste0(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = "")
  }
  
  lengths <- sample(18:20, n, replace = TRUE)
  targets <- vapply(lengths, random_dna, character(1))
  oligos <- vapply(targets, function(x) {
    as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
  }, character(1))
  
  tibble::tibble(
    input_order = seq_len(n),
    oligo_seq = oligos,
    name = targets,
    length = lengths,
    start = sample(1:2500, n),
    end = start + length - 1,
    gc_content = round(runif(n, 35, 70), 1),
    tox_score = round(runif(n, 45, 115), 1),
    off_target_score = round(runif(n, 0, 15), 2),
    n_distance_0 = sample(0:3, n, replace = TRUE),
    n_distance_1 = sample(0:20, n, replace = TRUE),
    n_distance_2 = sample(10:150, n, replace = TRUE),
    region_class = sample(c("Exonic", "Intronic", "UTR"), n, replace = TRUE),
    target_transcript = sample(
      c("ENST00000335137", "ENST00000413465", "ENST00000544455", NA),
      n,
      replace = TRUE
    ),
    motif_cor_score = round(runif(n, -0.6, 1.2), 2),
    max_pLI = round(runif(n, 0, 1), 3),
    min_LOEUF = round(runif(n, 0.1, 1.8), 3),
    conserved_in_mmusculus = sample(c(TRUE, FALSE, NA), n, replace = TRUE),
    CGs = sample(0:5, n, replace = TRUE),
    chr_start = sample(1000000:2000000, n),
    chr_end = chr_start + length - 1,
    PM_tot_freq = round(runif(n, 0, 0.08), 4),
    PM_max_freq = round(runif(n, 0, 0.05), 4),
    PM_count = sample(0:4, n, replace = TRUE),
    gene_hits_pm = sample(1:8, n, replace = TRUE),
    gene_hits_1mm = sample(5:80, n, replace = TRUE),
    NoRepeats = sample(1:3, n, replace = TRUE),
    sec_energy = round(runif(n, -8, 0), 2),
    duplex_energy = round(runif(n, -18, -4), 2),
    accessibility = signif(runif(n, 0, 0.001), 3)
  )
}

make_mock_offtargets <- function(n = 20) {
  random_dna <- function(len) {
    paste0(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = "")
  }
  
  tibble::tibble(
    gene_name = sample(c("TP53", "BRCA1", "DMD", "CFTR", "MYC", "APOE"), n, replace = TRUE),
    transcript = paste0("ENST", sample(10000000000:99999999999, n)),
    match_string_display = sample(c("==================", "====X=============", "=====I====D======="), n, replace = TRUE),
    subject_seq_display = vapply(sample(18:20, n, replace = TRUE), random_dna, character(1)),
    off_target_nt_length = sample(18:20, n, replace = TRUE),
    distance = sample(0:2, n, replace = TRUE),
    mismatches = sample(0:2, n, replace = TRUE),
    deletions = sample(0:1, n, replace = TRUE),
    insertions = sample(0:1, n, replace = TRUE),
    pLI = round(runif(n, 0, 1), 3),
    oe_lof_upper = round(runif(n, 0.1, 2), 3)
  )
}

make_mock_rnaseh <- function() {
  tibble::tibble(
    position = c("1 - 7", "2 - 8", "3 - 9", "4 - 10", "5 - 11"),
    average = round(runif(5, 0.2, 1.0), 3),
    window = c("ACTGCCA", "CTGCCAA", "TGCCAAT", "GCCAATG", "CCAATGG")
  )
}
#!!!!!!!!!!!!!

function(input, output, session) {
  
  
#!!!!!!!!!!!!!
  if (isFALSE(dev_mode)) {
  
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
  
  observeEvent(TRUE, {available <- check_gggenome()
    
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
      options(timeout = 60)
      
      martHS <- biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl"
      )
      
      test <- biomaRt::getBM(
        attributes = "ensembl_gene_id",
        filters    = "ensembl_gene_id",
        values     = "ENSG00000139618",
        mart       = martHS,
        bmHeader   = FALSE
      )
      
      is.data.frame(test) &&
        nrow(test) > 0 &&
        "ensembl_gene_id" %in% names(test)
      
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
  
  
    }
  #!!!!!!!!!!!! the checks in in dev mode rtight now
  
  use_reference_gene <- reactive({
    has_reference_gene(input$ensemble_id_input)
  })
  
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
      selected = last_reference_gene(),
      server   = TRUE,
      options  = list(
        valueField  = "value",
        labelField  = "label",
        searchField = c("label", "symbol", "ensembl"),
        placeholder = "Type gene symbol or Ensembl ID"
      )
    )
  })
  
  ## store sleected gene and then enable/disable it when selctin single aso input
  last_reference_gene <- reactiveVal("ENSG00000174775")
  observeEvent(input$ensemble_id_input, {
    if (!is.null(input$ensemble_id_input) &&
        nzchar(input$ensemble_id_input)) {
      last_reference_gene(input$ensemble_id_input)
    }
  }, ignoreInit = TRUE)
  
  ## lock single aso input after running since it messes up a lot of ui afterwards if zou change after the run
  
  ######$$$########
  observe({
    updateSelectizeInput(
      session,
      inputId  = "ensemble_id_input_patient",
      choices  = gene_choices_df,
      selected = last_reference_gene(),
      server   = TRUE,
      options  = list(
        valueField  = "value",
        labelField  = "label",
        searchField = c("label", "symbol", "ensembl"),
        placeholder = "Type gene symbol or Ensembl ID"
      )
    )
  })
  
  observe({
    updateSelectizeInput(
      session,
      inputId  = "ensemble_id_input_single_aso",
      choices  = gene_choices_df,
      selected = last_reference_gene(),
      server   = TRUE,
      options  = list(
        valueField  = "value",
        labelField  = "label",
        searchField = c("label", "symbol", "ensembl"),
        placeholder = "Type gene symbol or Ensembl ID"
      )
    )
  })
  
  observeEvent(input$ensemble_id_input_single_aso, {
    if (!is.null(input$ensemble_id_input_single_aso) &&
        nzchar(input$ensemble_id_input_single_aso)) {
      last_reference_gene(input$ensemble_id_input_single_aso)
    }
  }, ignoreInit = TRUE)
  
  # disable ASO generation when selction input summary onky
  observe({
    summary_only <- identical(input$patient_analysis_mode, "summary_only")
    
    if (summary_only) {
      shinyjs::hide("patient_optional_analysis_controls")
      shinyjs::hide("patient_filter_options_container")
    } else {
      shinyjs::show("patient_optional_analysis_controls")
      shinyjs::show("patient_filter_options_container")
    }
  })
  
  
  # helps with over dropping rows
  drop_overlapping_variants <- function(df) {
    if (nrow(df) <= 1) return(df)
    
    df <- df %>%
      dplyr::arrange(pos, dplyr::desc(nchar(ref)), dplyr::desc(nchar(alt)))
    
    keep <- rep(TRUE, nrow(df))
    occupied <- list()
    
    for (i in seq_len(nrow(df))) {
      this_start <- df$pos[i]
      this_end <- df$pos[i] + nchar(df$ref[i]) - 1L
      
      overlaps_existing <- FALSE
      
      for (j in seq_along(occupied)) {
        occ <- occupied[[j]]
        if (!(this_end < occ[1] || this_start > occ[2])) {
          overlaps_existing <- TRUE
          break
        }
      }
      
      if (overlaps_existing) {
        keep[i] <- FALSE
      } else {
        occupied[[length(occupied) + 1]] <- c(this_start, this_end)
      }
    }
    
    df[keep, , drop = FALSE]
  }
  
  
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
    raw_text <- toupper(trimws(
      if (is.null(input$aso_seq_input)) "" else input$aso_seq_input
    ))
    
    if (raw_text == "") {
      return(character(0))
    }
    
    asos <- regmatches(raw_text, gregexpr("[ACGT]+", raw_text))[[1]]
    asos <- unique(asos[nchar(asos) > 0])
    asos
  })
  
  output$parsed_asos <- renderText({
    asos <- parsed_asos()
    
    if (length(asos) == 0) {
      return("No valid ASO sequences detected.")
    }
    
    paste0(seq_along(asos), ". ", asos, collapse = "\n")
  })
  

  # -------------------------------------------------------------------
  
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
  
  
  #$$$#
  show_all_cols_patient <- reactiveVal(FALSE)
  show_all_cols_patient_snv <- reactiveVal(FALSE)
  
  observeEvent(input$toggle_cols_patient, {
    current <- show_all_cols_patient()
    show_all_cols_patient(!current)
    
    if (current) {
      updateActionButton(inputId = "toggle_cols_patient", label = "Show extended data")
    } else {
      updateActionButton(inputId = "toggle_cols_patient", label = "Hide extended data")
    }
  })
  
  observeEvent(input$toggle_cols_patient_snv, {
    current <- show_all_cols_patient_snv()
    show_all_cols_patient_snv(!current)
    
    if (current) {
      updateActionButton(
        inputId = "toggle_cols_patient_snv",
        label = "Show extended data"
      )
    } else {
      updateActionButton(
        inputId = "toggle_cols_patient_snv",
        label = "Hide extended data"
      )
    }
  })
  
  #$$$#
  
  # reset to standard values 
  
  observeEvent(input$reset_defaults, {
    shinyjs::reset("filters")
    updateCheckboxInput(session, "patient_filter_snv_asos", value = FALSE)
  })
  
  observeEvent(input$reset_defaults_patient, {
    updateCheckboxInput(session, "patient_ASO_ending_G", value = TRUE)
    updateCheckboxInput(session, "patient_Conserved_input", value = FALSE)
    updateCheckboxInput(session, "patient_polymorphism_input", value = TRUE)
    
    updateSliderInput(session, "patient_oligo_length_range", value = c(18, 20))
    
    updateCheckboxInput(session, "patient_gc_input", value = TRUE)
    updateNumericInput(session, "patient_gc_content-min", value = 40)
    updateNumericInput(session, "patient_gc_content-max", value = 60)
    updateSliderInput(session, "patient_gc_content-slider", value = c(40, 60))
    
    updateCheckboxInput(session, "patient_tox_input", value = TRUE)
    updateNumericInput(session, "patient_tox_score-min", value = 60)
    updateSliderInput(session, "patient_tox_score-slider", value = 60)
    
    updateCheckboxInput(session, "patient_perfect_input", value = TRUE)
    updateNumericInput(session, "patient_perfect_hits-max", value = 1)
    updateSliderInput(session, "patient_perfect_hits-slider", value = 1)
    
    updateCheckboxInput(session, "patient_mismatch_input", value = TRUE)
    updateNumericInput(session, "patient_mismatch_hits-max", value = 5)
    updateSliderInput(session, "patient_mismatch_hits-slider", value = 5)
    
    updateCheckboxInput(session, "patient_Poly_input", value = TRUE)
    updateNumericInput(session, "patient_pm_freq-min", value = 0)
    updateNumericInput(session, "patient_pm_freq-max", value = 0.05)
    updateSliderInput(session, "patient_pm_freq-slider", value = c(0, 0.05))
    
    updateCheckboxInput(session, "patient_filter_snv_asos", value = TRUE)
  })
  
  # Disable extend data button until table appears
  
  observe({
    shinyjs::disable("toggle_cols")
  })
  
  # sliding filter with numeric boxes
  tox_range <- rangeFilterServer("tox_score",    0, 136, fixed = "right")
  gc_range  <- rangeFilterServer("gc_content",   0, 100, fixed = "none")
  pm_range  <- rangeFilterServer("pm_freq",      0, 1,   fixed = "none")
  # acc_range <- rangeFilterServer("accessibility",0, 1,   fixed = "none")
  perf_range <- rangeFilterServer("perfect_hits", 0, 200, fixed = "left")
  mm_range   <- rangeFilterServer("mismatch_hits",0, 500, fixed = "left")
  
  
  #$$$#
  patient_tox_range  <- rangeFilterServer("patient_tox_score",    0, 136, fixed = "right")
  patient_gc_range   <- rangeFilterServer("patient_gc_content",   0, 100, fixed = "none")
  patient_pm_range   <- rangeFilterServer("patient_pm_freq",      0, 1,   fixed = "none")
  patient_perf_range <- rangeFilterServer("patient_perfect_hits", 0, 200, fixed = "left")
  patient_mm_range   <- rangeFilterServer("patient_mismatch_hits",0, 500, fixed = "left")
  #$$$#
  
  # starts timer, shows that the scirpt started and is loadin gthe progress bar
  
  observeEvent(input$run_button, {
    start_time <- Sys.time()
    
    showNotification(
      "Script started",
      type = "default",
      duration = 10,
      closeButton = TRUE
    )
    
    withProgress(message = "Running analysis...", value = 0, {

      
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
      
      txdb_hsa <- NULL
      gdb_hsa <- NULL
      HS <- NULL
      RNA_target <- NULL
      target_ranges <- NULL
      chr_coord <- NULL
      ensembl_ID <- NULL
      target_regions <- NULL
      
      need_txdb_hs <- TRUE
      
      
      
      if (need_txdb_hs) {
        txdb_hsa <- tryCatch({
            loadDb("/opt/ERASOR/txdb_hsa_biomart.db")
          }, error = function(e) {
            message("Primary txdb path not found, trying local path...")
            loadDb("../txdb_hsa_biomart.db")
          })
        
      incProgress(0.05, detail = "Loading sequences")
      
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
      }
      
      # ----------------------------------- milestone 2 ----------------------------
      print("milestone2: Subsetted genes from specified chromosomes")
      
      # Get the sequences *
      if (isTRUE(input$perfect_input) || isTRUE(input$mismatch_input) || isTRUE(use_reference_gene())) {
        HS <- getSeq(Hsapiens, gdb_hsa)}
        
      # ----------------------------------- milestone 3 ----------------------------
      print("milestone3: Saved human gene sequences")
      
      # Target collect the pre-mRNA sequence
      if (isTRUE(use_reference_gene())) {
      # Define wanted Ensembl ID
      ensembl_ID = input$ensemble_id_input
      
      # Retrieve a specific RNA target using the Ensembl ID
      RNA_target = HS[names(HS) == ensembl_ID]
      
      
      # ----------------------------------- milestone 4 ----------------------------
      print("milestone4: Saved input RNA target (ENSEMBL) information: \n")
      
      #filters on ensembl ID
      target_ranges = gdb_hsa[names(gdb_hsa) == ensembl_ID]
      
      if (length(RNA_target) != 1) {
        showNotification("Selected reference gene was not found.", type = "error", duration = 8)
        return(NULL)
      }
      
      if (length(target_ranges) != 1) {
        showNotification("Selected reference gene coordinates were not found.", type = "error", duration = 8)
        return(NULL)
      }
      
      # Extracts the chromosome name for target range,extracts the start and end
      # Position of the genomic region and note from which strand it is.
      # Positive strand 1 ('+'), negative strand -1 ('-'), or unspecified 0 ('*').
      chr_coord <- list(
        chr = as.character(seqnames(target_ranges)),
        start = as.integer(start(target_ranges)),
        end = as.integer(end(target_ranges)),
        strand = ifelse(as.character(strand(target_ranges)) == "+", 1L, -1L)
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
      
      # ----------------------------------- milestone 6 ----------------------------
      print("milestone6: Enumerated all possible ASO target sequences")
      
      } else {
        showNotification(
          "Please select a reference gene for General knockdown mode.",
          type = "error",
          duration = 8
        )
        return(NULL)
      }
      
      target_annotation <- add_empty_reference_columns(target_annotation)
      
      # ----------------------------------- milestone 7 ----------------------------
      print("milestone 7: Prefiltered Oligo sequences ending with G")
      
      # 3.4 Estimate Transcript Accessibility for the RNA Target at Single-Nucleotide Resolution
      if (isTRUE(use_reference_gene()) && isTRUE(input$linux_input)) {
        oligo_lengths <- sort(unique(target_annotation$length))
        
        accessibility <- RNAplfold_R(RNA_target, u.in = max(oligo_lengths)) %>%
          as_tibble() %>%
          mutate(end = 1:width(RNA_target)) %>%
          gather(length, accessibility, -end) %>%
          mutate(length = as.double(length))
        
        target_annotation <- left_join(
          target_annotation,
          accessibility,
          by = c("length", "end")
        )
      } else {
        target_annotation$accessibility <- NA_real_
      }
      # ----------------------------------- milestone 8 ----------------------------
      print("milestone 8: Calculated ViennaRNA accessibility score and filtering")
      
      if (isTRUE(use_reference_gene())) {
        nucleobase_seq <- reverseComplement(target_regions)
        target_annotation$oligo_seq <- as.character(nucleobase_seq[target_annotation$name])
      } else {
        target_annotation$oligo_seq <- toupper(parsed_asos())[target_annotation$input_order]
      }
      
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
      target_annotation$CGs = (target_annotation$length - nchar(gsub('CG', '', target_annotation$name))) / 2
      
      incProgress(0.1, detail = "Getting mouse ortholog data")
      
      # ----------------------------------- milestone 10 ---------------------------
      print("milestone 10: Calculated CG motifs")
      
      options(timeout = 60)
      
      if (isTRUE(use_reference_gene()) && isTRUE(biomart_available())) {
        
        martHS <- tryCatch(
          useEnsembl(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "hsapiens_gene_ensembl"
          ),
          error = function(e) useEnsembl(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "hsapiens_gene_ensembl",
            mirror  = "www"
          )
        )
        
        martMM <- tryCatch(
          useEnsembl(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "mmusculus_gene_ensembl"
          ),
          error = function(e) useEnsembl(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "mmusculus_gene_ensembl",
            mirror  = "www"
          )
        )
        
        ortho_ENS <- tryCatch(
          getBM(
            attributes = "mmusculus_homolog_ensembl_gene",
            filters = "ensembl_gene_id",
            values = ensembl_ID,
            mart = martHS,
            bmHeader = FALSE
          ),
          error = function(e) data.frame(
            mmusculus_homolog_ensembl_gene = character()
          )
        )
      
      # ----------------------------------- milestone 11 ---------------------------
      print("milestone 11: Get mouse ortholog data for genes")
        
        RNA_target_mouse <- tryCatch(
          DNAStringSet(
            getBM(
              attributes = c("gene_exon_intron", "ensembl_gene_id"),
              filters = "ensembl_gene_id",
              values = ortho_ENS$mmusculus_homolog_ensembl_gene,
              mart = martMM
            )$gene_exon_intron
          ),
          error = function(e) DNAStringSet()
        )
        
      } else {
        martHS <- NULL
        martMM <- NULL
        ortho_ENS <- NULL
        RNA_target_mouse <- DNAStringSet()
        target_annotation$conserved_in_mmusculus <- NA
      }
      
      # ----------------------------------- milestone 12 ---------------------------
      print("milestone 12: ")
      
      # Obtain all human polymorphisms for the RNA target
      if (isTRUE(use_reference_gene()) && 
          isTRUE(input$polymorphism_input) && 
          !is.null(martHS)) {
        
        PMs <- tryCatch(
          getBM(
            attributes = c("minor_allele_freq", "chromosome_start"),
            filters = "ensembl_gene_id",
            values = ensembl_ID,
            mart = martHS
          ) %>%
            as_tibble() %>%
            arrange(chromosome_start, desc(minor_allele_freq)) %>%
            filter(
              !is.na(minor_allele_freq),
              !duplicated(chromosome_start)
            ) %>%
            rename(
              chr_start = chromosome_start,
              PM_freq = minor_allele_freq
            ),
          error = function(e) tibble(
            chr_start = integer(),
            PM_freq = numeric()
          )
        )
        
      } else {
        PMs <- tibble(
          chr_start = integer(),
          PM_freq = numeric()
        )
      }
      
      ##If Ensembl is offline and still want to test -> load in manual test data.
      #PMs <- read.csv("~/PMs.csv")
      
      
      # ----------------------------------- milestone 13 ---------------------------
      print("milestone 13: Get human polymorfisms for RNA target")
      
      # Count Nucleobase sequence occurrences
      
      # Get the sequences
      
      if (isTRUE(use_reference_gene())){
      
      tr = target_annotation$name
      
      # Count the sequences by making it a table
      replica = table(tr)
      
      # Save it as NoRepeats
      target_annotation$NoRepeats = as.vector(replica[tr])
      
      } else {
        target_annotation$NoRepeats <- NA_integer_
      }
      
      # ----------------------------------- milestone 14 ---------------------------
      print("milestone 14: Count amount of times ASO sequence is repeated in target gene")
      
      # High-Frequency Polymorphisms analysis
      
      # Correcting end and start cord based on direction
      if (isTRUE(use_reference_gene())) {
      
        if (chr_coord$strand == 1L) {
          target_annotation$chr_start <- chr_coord$start + target_annotation$start - 1L
          target_annotation$chr_end   <- chr_coord$start + target_annotation$end - 1L
        } else {
          target_annotation$chr_start <- chr_coord$end - target_annotation$end + 1L
          target_annotation$chr_end   <- chr_coord$end - target_annotation$start + 1L
        }
      
      # ----------------------------------- milestone 15 ---------------------------
      print("milestone 15: Corrected start en end cord based on direction")

      target_annotation <- annotate_general_knockdown_region(
        target_annotation = target_annotation,
        txdb = txdb_hsa,
        ensembl_ID = ensembl_ID
      )
      
      } else {
        target_annotation$chr_start <- NA_integer_
        target_annotation$chr_end <- NA_integer_
        target_annotation$region_class <- NA_character_
        target_annotation$target_transcript <- NA_character_
      }
      
      # Keep unique names only and extract
      # Information base on chr_start from target.
      if (isTRUE(use_reference_gene()) && 
          isTRUE(input$polymorphism_input) && 
          nrow(PMs) > 0) {
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
      
      # Make table of mouse information
      if (isTRUE(use_reference_gene()) && length(RNA_target_mouse) > 0) {
        lm = width(RNA_target_mouse)
      
        oligo_lengths_current <- sort(unique(target_annotation$length))
        
       MM_tab = lapply(oligo_lengths_current, function(i) {
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
      
      target_annotation$conserved_in_mmusculus = target_annotation$name %in% RNAsitesMM
      
      } else {
        target_annotation$conserved_in_mmusculus <- NA
      }
      
      incProgress(0.25, detail = "Calculating secondary structure characteristics")
      
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
      
      nseq_pmfreq <- if (isTRUE(use_reference_gene()) && isTRUE(input$polymorphism_input)) {
        pm_rng <- pm_range()
        
        nseq_prefilter - nrow(
          target_annotation %>%
            dplyr::mutate(PM_tot_freq = tidyr::replace_na(PM_tot_freq, 0)) %>%
            dplyr::filter(PM_tot_freq >= pm_rng[1], PM_tot_freq <= pm_rng[2])
        )
      } else {
        NA_integer_
      }
      
      # remove accessibiility as filter but keep as output
      # nseq_accessible <- if (isTRUE(input$linux_input)) {
      #   acc_rng <- acc_range()
      #   nseq_prefilter - nrow(
      #     target_annotation %>%
      #       dplyr::filter(
      #         !is.na(accessibility),
      #         accessibility >= acc_rng[1],
      #         accessibility <= acc_rng[2]
      #       )
      #   )
      # } else {
      #   NA_integer_
      # }
      # 
      
     
      
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
      # remove accessibiltiy as filter but keep as output
      # if (isTRUE(input$linux_input) && isTRUE(input$Accessibility_input)) {
      #   ta_prev <- ta
      #   if (isTRUE(input$linux_input) && isTRUE(input$Accessibility_input)) {
      #     ta_prev <- ta
      #     rng <- acc_range()
      #     ta <- apply_range(ta, "accessibility", rng)
      #     
      #     if (nrow(ta) == 0) {
      #       ta <- ta_prev
      #       showNotification("Filter 'accessibility' removed all rows; reverting.", type = "warning")
      #     }
      #   }
      #   
      #   if (nrow(ta) == 0) {
      #     ta <- ta_prev
      #     message("Filter 'accessibility' removed all rows; reverting to previous dataset.")
      #     showNotification("Filter 'accessibility' removed all rows; reverting to previous dataset.", type = "warning")
      #   }
      # }
      # 
      # 4) polymorphism NA->0 (geen filter; alleen transform)
      if (isTRUE(use_reference_gene()) && isTRUE(input$polymorphism_input)) {
        ta <- ta %>%
          mutate(across(c(PM_max_freq, PM_tot_freq, PM_count), ~ tidyr::replace_na(., 0.0)))
      }
      
      # 5) polymorphism frequency filter
      if (isTRUE(use_reference_gene()) && isTRUE(input$polymorphism_input) && isTRUE(input$Poly_input)) {
        ta_prev <- ta
        if (isTRUE(use_reference_gene()) && isTRUE(input$polymorphism_input)) {
          ta <- ta %>%
            mutate(across(c(PM_max_freq, PM_tot_freq, PM_count), ~ tidyr::replace_na(., 0.0)))
        }
        
        if (isTRUE(use_reference_gene()) && isTRUE(input$polymorphism_input) && isTRUE(input$Poly_input)) {
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
      if (isTRUE(use_reference_gene()) && isTRUE(input$Conserved_input)) {
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
      if (!is.null(HS) && nrow(target_annotation) > 0) {
      
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
      
      uni_tar <-bind_rows(uni_tar_list)
      
      # ----------------------------------- milestone 21 ---------------------------
      print("milestone 21: Matched ASO sequences to potential off-targets (perfect match, one mismatch")
      
      target_annotation <- left_join(target_annotation, uni_tar, by = c("name", "length"))
      
      } else {
        target_annotation$gene_hits_pm <- NA_integer_
        target_annotation$gene_hits_1mm <- NA_integer_
      }
      
      incProgress(0.4, detail = "Searching for off-targets")
      
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
      
      target_annotation <- target_annotation[order(target_annotation$gene_hits_pm, target_annotation$gene_hits_1mm), ]
      
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
              dplyr::mutate(
                name = toupper(trimws(as.character(name))),
                distance = mismatches + deletions + insertions
              ) %>%
              dplyr::distinct(
                name,
                gene_name,
                match_string,
                query_seq,
                .keep_all = TRUE
              )
            
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
     
      off_targets_total <- tibble() 
      if (!is.null(summary_server)) {
        tmp <- tempfile(fileext = ".bgz")
        
        # Download GnomAD lof metrics by gene
        download.file(
          url <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", 
          destfile = tmp,
          mode = "wb"
        )
        
        gnomad_df <- read_tsv(
          tmp,
          show_col_types = FALSE,
          col_select = c(gene, transcript, pLI, oe_lof_upper)
        )
        
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
        #####&*#####
        # oe_lof_min <- off_targets_total %>%
        #   filter(distance <= 1) %>%
        #   group_by(name) %>%
        #   summarise(min_oe_lof = min(oe_lof, na.rm = TRUE), .groups = "drop")
        #####&*#####
        
        minloeuf_maxpli <- off_targets_total %>%
          filter(distance <= 1) %>%
          group_by(name) %>%
          summarise(
            max_pLI = if (all(is.na(pLI))) NA_real_ else max(pLI, na.rm = TRUE),
            min_LOEUF = if (all(is.na(oe_lof_upper))) NA_real_ else min(oe_lof_upper, na.rm = TRUE),
            .groups = "drop"
          )
        
        # 3. Deze twee samenvattingen samenvoegen tot één tabel per name
        off_summary <- dist_counts %>%
          left_join(minloeuf_maxpli, by = "name")
        
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
      
      target_annotation <- add_rnaseh_summary_columns(
        target_annotation,
        mod_5prime = 0,
        mod_3prime = 0
      )
      
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
      
      
      # fitlering table
      # number of ASOs that actually end up in the main results table
      n_final_selection <- nrow(target_annotation)
      
      filter_numbers <- tibble(
        Filter = c(
          "Prefiltered",
          "ASO ending with G",
          "Tox score",
          "GC content",
          "PM frequency",
          "Final selection"
        ),
        `Number of ASOs erased` = c(
          nseq_prefilter,
          nseq_ending_G,
          nseq_toxscore,
          nseq_gc,
          nseq_pmfreq,
          n_final_selection
        ),
        Percentage = c(
          "100%",
          paste0(round(100 * nseq_ending_G / nseq_prefilter), "%"),
          paste0(round(100 * nseq_toxscore / nseq_prefilter), "%"),
          paste0(round(100 * nseq_gc / nseq_prefilter), "%"),
          paste0(round(100 * nseq_pmfreq / nseq_prefilter), "%"),
          paste0(formatC(100 * n_final_selection / nseq_prefilter, format = "f", digits = 2), "%")
        )
      )
      
      output$unfiltered_results_table <- renderDT({
        
        datatable(
          filter_numbers,
          rownames = FALSE,
          options = list(
            dom = "t",
            ordering = FALSE
          ),
          class = "compact stripe"
        ) %>%
          formatStyle(
            names(filter_numbers),
            textAlign = "center"
          ) %>%
          formatStyle(
            "Filter",
            target = "row",
            fontWeight = styleEqual("Final selection", "bold")
          )
      })
      
      # Render the tables.
      ################ make download table identical to viewed table
      results1_lookup <- reactive({
        req(target_annotation)
        
        df <- target_annotation
        
        if ("off_target_score" %in% names(df)) {
          df <- df %>% dplyr::arrange(dplyr::coalesce(off_target_score, Inf))
        }
        
        df
      })
      
      results1_data <- reactive({
        df <- results1_lookup()
        
        column_order <- c(
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
          "region_class",
          "target_transcript",
          "motif_cor_score",
          "max_pLI",
          "min_LOEUF",
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
          "sec_energy",
          "duplex_energy",
          "accessibility",
          "rnaseh_max_score",
          "rnaseh_mean_score"
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
          region_class           = "Target region",
          target_transcript      = "Target transcript(s)",
          sec_energy             = "ASO self-folding energy",
          duplex_energy          = col_with_tooltip(
            "ASO duplex energy",
            "Predicted energy needed to break the binding between two ASOs. More negative values generally indicate stronger binding, which is unfavorable for ASO design."),
          max_pLI                = "Max. off-target pLI",
          min_LOEUF              = "Min. off-target LOEUF",
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
          accessibility          = "Accessibility",
          rnaseh_max_score       = "Max RNase H score",
          rnaseh_mean_score      = "Average RNase H score"
        )
        
        to_rename <- intersect(names(df), names(column_names))
        names(df)[match(to_rename, names(df))] <- column_names[to_rename]
        
        if (!isTRUE(use_reference_gene())) {
          df <- df %>% dplyr::select(-dplyr::any_of(reference_required_display_cols()))
        }
        
        if ("Off-target score" %in% names(df)) {
          df <- df %>% dplyr::arrange(dplyr::coalesce(`Off-target score`, Inf))
        }
        
        df
      })
      
      results1_export_data <- reactive({
        df <- results1_data()
        names(df) <- strip_header_html(names(df))
        df
      })
      
      observe({
        df <- results1_data()
        
        display_names <- names(df)
        plain_names <- strip_header_html(display_names)
        
        choice_map <- stats::setNames(display_names, plain_names)
        
        selected_value <- if ("Off-target score" %in% plain_names) {
          display_names[match("Off-target score", plain_names)]
        } else {
          display_names[1]
        }
        
        updateSelectInput(
          session,
          "main_table_sort_col",
          choices = choice_map,
          selected = selected_value
        )
      })
      
      results1_sorted <- reactive({
        df <- results1_data()
        
        sort_col <- input$main_table_sort_col
        sort_dir <- input$main_table_sort_dir
        
        if (is.null(sort_col) || !(sort_col %in% names(df))) {
          return(df)
        }
        
        if (identical(sort_dir, "desc")) {
          df <- df %>%
            dplyr::arrange(dplyr::desc(is.na(.data[[sort_col]])),
                           dplyr::desc(.data[[sort_col]]))
        } else {
          df <- df %>%
            dplyr::arrange(is.na(.data[[sort_col]]),
                           .data[[sort_col]])
        }
        
        df
      })
      #############
      
      output$results1 <- DT::renderDT({
        df <- results1_sorted()
        
        main_cols <- c(
          "ASO sequence",
          "Target (DNA)",
          "ASO length (nt)",
          if (isTRUE(use_reference_gene())) c("Start position in gene", "End position in gene"),
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
          escape    = FALSE,
          options   = list(
            dom = "tip",
            columnDefs = c(
              column_defs,
              list(list(className = "dt-center", targets = "_all"))
            )
          ),
          class = "compact stripe cell-border"
        )
        
        DT::formatRound(
          dt,
          columns = intersect(c("GC content (%)", "PM total freq.", "PM max freq."), names(df)),
          digits  = 0
        )
      })
      
      observe({
        req(results1_sorted())
        shinyjs::enable("toggle_cols")
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
              "pLI",
              "oe_lof_upper"
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
          deletions               = "Insertions", # these are flipped bc it is the easiest way to swicth from query insertions to subject deletion
          insertions              = "Deletions",
          offtarget_accessibility = "Accessibility (off-target)",
          pLI                     = "GnomAD pLI",
          oe_lof_upper            = "GnomAD LOEUF"
        )
        
        nm <- names(df_view)
        nm2 <- ifelse(nm %in% names(display_map), display_map[nm], nm)
        names(df_view) <- make.unique(nm2, sep = " (dup) ")
        
        dist_idx <- which(names(df_view) == "Number of Mismatches/Indels")
        if (length(dist_idx) == 0) {
          dist_idx <- which(grepl("^Number of Mismatches/Indels", names(df_view)))[1]
        }
        distance_col <- dist_idx - 1
        
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
          DT::datatable(
            rnaseh_view,
            selection = list(mode = "single", selected = 1),
            options = list(
              paging = FALSE,
              searching = FALSE,
              info = FALSE,
              dom = "t"
            )
          )
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
              "<h5 style='margin-top:0;'>Oligo sequence (ASO): </h5>",
              "<div style='font-family: monospace; white-space: pre; font-size: 20px; line-height: 1.7;'>",
              oligo_visual_fw,
              "</div>",
              "<h5 style='margin-bottom:6px;'>RNA sequence: </h5>",
              "<div style='font-family: monospace; white-space: pre; font-size: 20px; line-height: 1.7;'>",
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
          
          rnaseh_stored(rnaseh_data)
          
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
            DT::datatable(
              rnaseh_view,
              selection = list(mode = "single", selected = 1),
              options = list(
                paging = FALSE,
                searching = FALSE,
                info = FALSE,
                dom = "t"
              )
            )
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
        table_data_reactive = results1_lookup,
        off_targets_total = off_targets_total,
        selected_target = selected_target,
        oligo_sequence = oligo_sequence,
        rnaseh_stored = rnaseh_stored,
        current_seq = current_seq,
        current_mismatch = current_mismatch,
        current_offtargets = current_offtargets,
        cached_results = cached_results
      )
      
      #leftover????
      # # Table two handler call
      # handle_table_events(
      #   input = input,
      #   output = output,
      #   session = session,
      #   table_id = "results2",
      #   table_data = nucleobase_select,
      #   off_targets_total = off_targets_total,
      #   selected_target = selected_target,
      #   oligo_sequence = oligo_sequence,
      #   rnaseh_stored = rnaseh_stored,
      #   current_seq = current_seq,
      #   current_mismatch = current_mismatch,
      #   current_offtargets = current_offtargets,
      #   cached_results = cached_results
      # )
      
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
            results1_export_data(),
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
  
  #$$$#
  patient_summary_data <- reactiveVal(NULL)
  patient_variants_data <- reactiveVal(NULL)
  patient_reference_results_raw <- reactiveVal(tibble())
  patient_snv_results_raw <- reactiveVal(tibble())
  patient_offtargets_total <- reactiveVal(tibble())
  
  selected_target_patient <- reactiveVal(NULL)
  oligo_sequence_patient <- reactiveVal(NULL)
  rnaseh_stored_patient <- reactiveVal(NULL)
  
  current_seq_patient <- reactiveVal(NULL)
  current_mismatch_patient <- reactiveVal(2)
  current_offtargets_patient <- reactiveVal(NULL)
  cached_results_patient <- reactiveVal(list())
  patient_snv_only_mode <- reactive({
    identical(input$patient_analysis_mode, "snv_only")
  })
  
  # Disable Patient tab result buttons until results exist
  observe({
    shinyjs::disable("toggle_cols_patient")
    shinyjs::disable("Download_unfiltered_patient_reference")
    shinyjs::disable("Download_filtered_patient_reference")
    shinyjs::disable("toggle_cols_patient_snv")
    shinyjs::disable("Download_patient_snv_results")
    shinyjs::disable("Download_patient_ambiguous_results")
  })
  
  patient_results1_lookup <- reactive({
    df <- patient_reference_results_raw()
    
    if (is.null(df) || nrow(df) == 0) {
      return(tibble())
    }
    
    if ("off_target_score" %in% names(df)) {
      df <- df %>%
        dplyr::arrange(dplyr::coalesce(off_target_score, Inf))
    }
    
    df
  })
  
  patient_snv_results_lookup <- reactive({
    df <- patient_snv_results_raw()
    
    if (is.null(df) || nrow(df) == 0) {
      return(tibble())
    }
    
    if ("off_target_score" %in% names(df)) {
      df <- df %>%
        dplyr::arrange(
          dplyr::coalesce(off_target_score, Inf),
          variant_note,
          start
        )
    } else {
      df <- df %>%
        dplyr::arrange(
          variant_note,
          start
        )
    }
    
    df
  })
  
  patient_snv_results_data <- reactive({
    df <- patient_snv_results_lookup()
    
    if (is.null(df) || nrow(df) == 0) {
      return(tibble())
    }
    
    if (!"target_genomic_display" %in% names(df)) {
      df$target_genomic_display <- df$target_genomic_seq
    }
    if (!"target_gene_display" %in% names(df)) {
      df$target_gene_display <- df$name
    }
    
    snv_display_map <- c(
      phase_status = "Phase status",
      genotype = "Genotype",
      variant_type_patient = "Variant type",
      variant_note = "Variant note",
      oligo_seq = "ASO sequence",
      target_genomic_display = "Target genomic (forward)",
      target_gene_display = "Target gene orientation",
      length = "ASO length (nt)",
      start = "Start position in gene",
      end = "End position in gene",
      gc_content = "GC content (%)",
      tox_score = "Acute neurotox score",
      rnaseh_max_score = "Max RNase H score",
      rnaseh_mean_score = "Average RNase H score",
      region_class = "Target region",
      target_transcript = "Target transcript(s)",
      chr_start = "Chromosome start pos.",
      chr_end = "Chromosome end pos.",
      NoRepeats = "Number of repeats",
      sec_energy = "ASO self-folding energy",
      duplex_energy = "ASO duplex energy",
      off_target_score = "Off-target score",
      n_distance_0 = "Perfect matches (GGGenome)",
      n_distance_1 = "Off-targets 1 mismatch (GGGenome)",
      n_distance_2 = "Off-targets 2 mismatches (GGGenome)",
      gene_hits_pm = "Perfect matches (Pedersen)",
      gene_hits_1mm = "Off-targets 1 mismatch (Pedersen)",
      max_pLI = "Max. off-target pLI",
      min_LOEUF = "Min. off-target LOEUF"
    )
    
    for (col in names(snv_display_map)) {
      if (!col %in% names(df)) {
        df[[col]] <- NA
      }
    }
    
    df <- df %>%
      dplyr::select(dplyr::all_of(names(snv_display_map)))
    
    names(df) <- unname(snv_display_map[names(df)])
    
    df
  })
  
  observe({
    df <- patient_snv_results_data()
    
    if (is.null(df) || nrow(df) == 0) {
      return()
    }
    
    display_names <- names(df)
    plain_names <- strip_header_html(display_names)
    
    choice_map <- stats::setNames(display_names, plain_names)
    
    selected_value <- if ("Off-target score" %in% plain_names) {
      display_names[match("Off-target score", plain_names)]
    } else if ("Variant note" %in% plain_names) {
      display_names[match("Variant note", plain_names)]
    } else {
      display_names[1]
    }
    
    updateSelectInput(
      session,
      "patient_snv_table_sort_col",
      choices = choice_map,
      selected = selected_value
    )
  })
  
  patient_snv_results_sorted <- reactive({
    df <- patient_snv_results_data()
    
    if (is.null(df) || nrow(df) == 0) {
      return(tibble())
    }
    
    sort_col <- input$patient_snv_table_sort_col
    sort_dir <- input$patient_snv_table_sort_dir
    
    if (is.null(sort_col) || !(sort_col %in% names(df))) {
      return(df)
    }
    
    if (identical(sort_dir, "desc")) {
      df <- df %>%
        dplyr::arrange(
          dplyr::desc(is.na(.data[[sort_col]])),
          dplyr::desc(.data[[sort_col]])
        )
    } else {
      df <- df %>%
        dplyr::arrange(
          is.na(.data[[sort_col]]),
          .data[[sort_col]]
        )
    }
    
    df
  })
  
  handle_patient_table_events <- function(table_id, table_data_reactive) {
    observeEvent({
      list(
        input[[paste0(table_id, "_cell_clicked")]],
        input[[paste0(table_id, "_rows_selected")]]
      )
    }, {
      row <- input[[paste0(table_id, "_rows_selected")]]
      
      if (is.null(row) || length(row) == 0) {
        return()
      }
      
      df_current <- table_data_reactive()
      
      if (is.null(df_current) || nrow(df_current) == 0) {
        return()
      }
      
      row_data <- df_current[row, , drop = FALSE]
      
      if (!"name" %in% names(row_data) || !"oligo_seq" %in% names(row_data)) {
        return()
      }
      
      seq <- toupper(trimws(as.character(row_data$name[[1]])))
      
      if (!grepl("^[ACGT]+$", seq)) {
        return()
      }
      
      selected_target_patient(row_data)
      oligo_sequence_patient(row_data$oligo_seq[[1]])
      
      rnaseh_data <- rnaseh_results(
        selected_row_name = row_data$name[[1]],
        oligo_seq = row_data$oligo_seq[[1]],
        mod_5prime = input$patient_mod_5prime,
        mod_3prime = input$patient_mod_3prime
      )
      
      rnaseh_stored_patient(rnaseh_data)
      
      output$patient_rnaseh_title <- renderText({
        paste0("RNase H results for: ", row_data$name[[1]])
      })
      
      output$patient_rnaseh_info <- renderUI({
        HTML(
          paste0(
            "Length of sequence: ",
            row_data$length[[1]],
            "<br>",
            "Oligo sequence (ASO): ",
            row_data$oligo_seq[[1]]
          )
        )
      })
      
      output$patient_rnaseh_results <- DT::renderDT({
        rnaseh_view <- rename_rnaseh_cols(rnaseh_data)
        
        DT::datatable(
          rnaseh_view,
          selection = list(mode = "single", selected = 1),
          options = list(
            paging = FALSE,
            searching = FALSE,
            info = FALSE,
            dom = "t"
          )
        )
      })
      
      output$patient_cleavage_visual <- renderUI(div())
      
      ot_total <- patient_offtargets_total()
      
      current_seq_patient(seq)
      current_mismatch_patient(2)
      
      updateSelectInput(
        session,
        "patient_user_mismatch",
        selected = 2
      )
      
      output$patient_offtarget_title <- renderText({
        paste0("Off target results for: ", seq)
      })
      
      output$patient_aso_seq <- renderText({
        paste0("ASO sequence: ", row_data$oligo_seq[[1]])
      })
      
      if (!is.null(ot_total) && nrow(ot_total) > 0) {
        default_subset <- ot_total %>%
          dplyr::mutate(name_lookup = toupper(trimws(as.character(name)))) %>%
          dplyr::filter(name_lookup == seq, distance <= 2) %>%
          dplyr::select(-name_lookup)
        
        current_offtargets_patient(default_subset)
        
        output$patient_numb_offtargets <- renderText({
          paste0("# off targets: ", nrow(default_subset))
        })
        
      } else {
        current_offtargets_patient(tibble())
        
        output$patient_numb_offtargets <- renderText({
          "No off-target data available for this ASO."
        })
      }
      
      updateTabsetPanel(
        session,
        "tabs_main_patient",
        selected = "Off target results"
      )
    })
  }
  
  handle_patient_table_events(
    table_id = "patient_results1",
    table_data_reactive = patient_results1_lookup
  )
  
  handle_patient_table_events(
    table_id = "patient_snv_results",
    table_data_reactive = patient_snv_results_sorted
  )
  
  observeEvent(input$patient_rnaseh_results_rows_selected, {
    row_number <- input$patient_rnaseh_results_rows_selected
    
    if (length(row_number) == 0) {
      return()
    }
    
    row_data <- selected_target_patient()
    
    if (is.null(row_data)) {
      return()
    }
    
    rnaseh_data <- rnaseh_stored_patient()
    
    if (is.null(rnaseh_data)) {
      return()
    }
    
    selected_row <- rnaseh_data[row_number, ]
    
    reverse_string <- function(x) {
      paste0(rev(strsplit(x, "")[[1]]), collapse = "")
    }
    
    oligo_seq <- row_data$oligo_seq[[1]]
    
    mod5 <- input$patient_mod_5prime
    mod3 <- input$patient_mod_3prime
    
    oligo_len <- nchar(oligo_seq)
    
    mod5_region <- substr(oligo_seq, 1, mod5)
    mod3_region <- substr(oligo_seq, oligo_len - mod3 + 1, oligo_len)
    mid_region <- substr(oligo_seq, mod5 + 1, oligo_len - mod3)
    
    position_string <- stringr::str_split_1(selected_row$position, " - ")
    start_pos <- as.numeric(position_string[1])
    
    target_seq_fw <- row_data$name[[1]]
    target_seq_rv <- reverse_string(target_seq_fw)
    
    rna_len <- nchar(target_seq_fw)
    
    cleavage_start_fw <- start_pos
    cleavage_pos_fw <- cleavage_start_fw + 6
    cleavage_pos_rv <- rna_len - cleavage_pos_fw
    cleavage_site_up <- cleavage_pos_rv - 1
    cleavage_site_down <- cleavage_pos_rv + 7
    
    upstream <- substr(target_seq_rv, 1, cleavage_site_up - 1)
    site_start <- substr(target_seq_rv, cleavage_site_up, cleavage_pos_rv - 1)
    cut_site <- substr(target_seq_rv, cleavage_pos_rv, cleavage_pos_rv)
    site_down <- substr(target_seq_rv, cleavage_pos_rv + 1, cleavage_site_down)
    downstream <- substr(target_seq_rv, cleavage_site_down + 1, rna_len)
    
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
    
    output$patient_cleavage_visual <- renderUI({
      HTML(
        paste0(
          "<h5 style='margin-top:0;'>Oligo sequence (ASO): </h5>",
          "<div style='font-family: monospace; white-space: pre; font-size: 20px; line-height: 1.7;'>",
          oligo_visual_fw,
          "</div>",
          "<h5 style='margin-bottom:6px;'>RNA sequence: </h5>",
          "<div style='font-family: monospace; white-space: pre; font-size: 20px; line-height: 1.7;'>",
          rna_visual_rv,
          "</div>"
        )
      )
    })
  })
  
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
    
    shinyjs::disable("toggle_cols_patient")
    shinyjs::disable("Download_unfiltered_patient_reference")
    shinyjs::disable("Download_filtered_patient_reference")
    shinyjs::disable("toggle_cols_patient_snv")
    shinyjs::disable("Download_patient_snv_results")
    shinyjs::disable("Download_patient_ambiguous_results")
    
    req(input$ensemble_id_input_patient)
    req(input$patient_variant_file)
    
    run_patient_summary_only <- identical(input$patient_analysis_mode, "summary_only")
    run_patient_snv_only <- identical(input$patient_analysis_mode, "snv_only")
    
    run_reference_patient_asos <- !run_patient_summary_only && !run_patient_snv_only
    run_patient_snv_asos <- !run_patient_summary_only
    
    withProgress(message = "Preparing patient-specific input...", value = 0, {
      
      incProgress(0.05, detail = "Loading gene annotation")
      
      txdb_hsa <- tryCatch({
        loadDb("/opt/ERASOR/txdb_hsa_biomart.db")
      }, error = function(e) {
        loadDb("../txdb_hsa_biomart.db")
      })
      
      gdb_hsa <- genes(txdb_hsa)
      seqlevelsStyle(gdb_hsa) <- seqlevelsStyle(Hsapiens)
      
      common_chrs <- intersect(seqlevels(gdb_hsa), seqlevels(Hsapiens))
      gdb_hsa <- keepSeqlevels(gdb_hsa, common_chrs, pruning.mode = "coarse")
      
      print("patient milestone 1: loaded human database object")
      
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
      
      
      incProgress(0.12, detail = "Checking uploaded files")
      
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
      
      print("patient milestone 2: staged uploaded variant and index files")
      
      flank <- as.integer(input$patient_flank)
      region_string <- trimws(if (is.null(input$patient_region_input)) "" else input$patient_region_input)
      
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
        region_string <- paste0(
          gene_chr, ":",
          max(1L, gene_start - flank), "-",
          gene_end + flank
        )
      }
      
      region_info <- parse_region_string(region_string)
      
      incProgress(0.20, detail = "Reading patient variants")
      
      variants_raw <- tryCatch(
        read_patient_variants_resilient(
          variant_path  = staged$variant_path,
          region_string = region_string,
          is_indexed    = staged$is_indexed
        ),
        error = function(e) {
          showNotification(
            paste("Variant reading failed:", conditionMessage(e)),
            type = "error",
            duration = NULL
          )
          return(NULL)
        }
      )
      
      if (is.null(variants_raw)) {
        return()
      }
      
      print("patient milestone 3: read patient variants from VCF/BCF")
      
      query_mode_used <- attr(variants_raw, "query_mode_used")
      query_region_used <- attr(variants_raw, "query_region_used")
      
      target_chr_chr <- normalize_chr_style(region_info$chr)
      target_chr_nochr <- normalize_chr_style_nochr(region_info$chr)
      
      variants_region <- variants_raw %>%
        dplyr::mutate(
          chr_original = chr,
          chr_norm_chr = normalize_chr_style(chr),
          chr_norm_nochr = normalize_chr_style_nochr(chr)
        ) %>%
        dplyr::filter(
          (chr_norm_chr == target_chr_chr | chr_norm_nochr == target_chr_nochr),
          pos >= region_info$start,
          pos <= region_info$end
        ) %>%
        dplyr::arrange(pos) %>%
        dplyr::select(-chr_norm_chr, -chr_norm_nochr)
      
      variants_region <- variants_region %>%
        dplyr::filter(
          !is.na(gt),
          !gt %in% c("0|0", "0/0", "./.", ".|.")
        )
      
      variants_region_expanded <- expand_multiallelic_patient_variants(variants_region)
      
      variants_region_expanded <- validate_variants_against_reference(variants_region_expanded)
      
      bad_ref <- variants_region_expanded %>%
        dplyr::filter(!ref_matches_genome)
      
      if (nrow(bad_ref) > 0) {
        print("REF mismatches against genome used by app:")
        print(bad_ref %>% dplyr::select(chr, pos, ref, alt, ref_observed_genome) %>% head(20))
        return()
      }
      
      ##&*###
      print("Original patient variants:")
      print(variants_region %>% dplyr::select(chr, pos, ref, alt, gt, variant_type))
      
      print("Expanded patient variants:")
      print(variants_region_expanded %>% dplyr::select(chr, pos, ref, alt, gt, variant_type, variant_on))
      ###&*###
      
      incProgress(0.25, detail = "Getting reference sequence")
      
      region_gr <- GenomicRanges::GRanges(
        seqnames = region_info$chr,
        ranges   = IRanges::IRanges(
          start = region_info$start,
          end   = region_info$end
        ),
        strand = "+"
      )
      
      GenomeInfoDb::seqlevelsStyle(region_gr) <- GenomeInfoDb::seqlevelsStyle(Hsapiens)
      ref_seq <- getSeq(Hsapiens, region_gr)[[1]]
      
      print("patient milestone 4: extracted reference sequence")
      
      incProgress(0.30, detail = "Rendering summary")
      
      
      patient_summary_df <- tibble(
        Metric = c(
          "Selected gene",
          "Gene chromosome",
          "Gene start",
          "Gene end",
          "Gene strand",
          "Reference sequence length",
          "Requested region",
          "Region query actually used",
          "Variant read mode",
          "Flank",
          "Variants in selected region"
        ),
        Value = c(
          ensembl_ID,
          gene_chr,
          gene_start,
          gene_end,
          gene_strand,
          nchar(as.character(ref_seq)),
          region_string,
          ifelse(is.null(query_region_used), NA_character_, query_region_used),
          ifelse(is.null(query_mode_used), NA_character_, query_mode_used),
          flank,
          nrow(variants_region)
        )
      )
      
      patient_summary_data(patient_summary_df)
      patient_variants_data(variants_region)
      
      print("patient milestone 5: created patient summary and SNV overview")
      
      output$patient_summary <- renderTable({
        req(patient_summary_data())
        patient_summary_data()
      }, striped = TRUE, bordered = TRUE, spacing = "s")
      
      output$patient_snvs_table <- renderDT({
        req(patient_variants_data())
        
        df_raw <- patient_variants_data()
        
        chr_col <- if ("chr_original" %in% names(df_raw)) "chr_original" else "chr"
        
        df_show <- dplyr::transmute(
          df_raw,
          chr = .data[[chr_col]],
          pos = if ("pos" %in% names(df_raw)) pos else NA_integer_,
          ref = if ("ref" %in% names(df_raw)) ref else NA_character_,
          alt = if ("alt" %in% names(df_raw)) alt else NA_character_,
          variant_type = if ("variant_type" %in% names(df_raw)) variant_type else NA_character_,
          gt = if ("gt" %in% names(df_raw)) gt else NA_character_,
          hap1_allele = if ("hap1_allele" %in% names(df_raw)) hap1_allele else NA_character_,
          hap2_allele = if ("hap2_allele" %in% names(df_raw)) hap2_allele else NA_character_,
          variant_on = if ("variant_on" %in% names(df_raw)) variant_on else NA_character_,
          allele_note = if ("allele_note" %in% names(df_raw)) allele_note else NA_character_
        ) %>%
          dplyr::rename(
            "Chromosome" = chr,
            "Position" = pos,
            "Reference" = ref,
            "Patient" = alt,
            "Variant Type" = variant_type,
            "Genotype" = gt,
            "Allele 1" = hap1_allele,
            "Allele 2" = hap2_allele,
            "Variant Location" = variant_on,
            "Phase Info" = allele_note
          )
        
        DT::datatable(
          df_show,
          rownames = FALSE,
          class = "nowrap",
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            columnDefs = list(
              list(className = "dt-center", targets = "_all")
            )
          )
        )
      })
      
      output$Download_unfiltered_patient <- downloadHandler(
        filename = function() {
          paste0("patient_input_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
        },
        content = function(file) {
          req(patient_summary_data())
          write.csv(patient_summary_data(), file, row.names = FALSE)
        }
      )
      
      # -------------------------------------------------------------
      # Stop here if user only wants patient input summary
      # -------------------------------------------------------------
      if (identical(input$patient_analysis_mode, "summary_only")) {
        
        output$patient_unfiltered_results_table <- renderDT({
          datatable(
            data.frame(
              Message = "ASO generation was disabled. This run only created the patient input summary."
            ),
            rownames = FALSE,
            options = list(
              dom = "t",
              ordering = FALSE,
              searching = FALSE,
              paging = FALSE
            ),
            class = "compact stripe"
          )
        })
        
        output$patient_results1 <- renderDT({
          datatable(
            data.frame(
              Message = "No reference-gene ASO results. Switch patient analysis mode to 'Input summary + ASO generation'."
            ),
            rownames = FALSE,
            options = list(
              dom = "t",
              ordering = FALSE,
              searching = FALSE,
              paging = FALSE
            ),
            class = "compact stripe"
          )
        })
        
        output$patient_snv_results <- renderDT({
          datatable(
            data.frame(
              Message = "No patient-derived SNV ASO results. Switch patient analysis mode to 'Input summary + ASO generation'."
            ),
            rownames = FALSE,
            options = list(
              dom = "t",
              ordering = FALSE,
              searching = FALSE,
              paging = FALSE
            ),
            class = "compact stripe"
          )
        })
        
        shinyjs::disable("toggle_cols_patient_snv")
        
        output$Download_unfiltered_patient_reference <- downloadHandler(
          filename = function() {
            paste0(
              "patient_reference_unfiltered_results_",
              format(Sys.time(), "%Y%m%d_%H%M%S"),
              ".csv"
            )
          },
          content = function(file) {
            write.csv(
              data.frame(Message = "ASO generation disabled for this run."),
              file,
              row.names = FALSE
            )
          }
        )
        
        output$Download_filtered_patient_reference <- downloadHandler(
          filename = function() {
            paste0(
              "patient_reference_filtered_results_",
              format(Sys.time(), "%Y%m%d_%H%M%S"),
              ".csv"
            )
          },
          content = function(file) {
            write.csv(
              data.frame(Message = "ASO generation disabled for this run."),
              file,
              row.names = FALSE
            )
          }
        )
        
        output$Download_patient_snv_results <- downloadHandler(
          filename = function() {
            paste0(
              "patient_snv_aso_results_",
              format(Sys.time(), "%Y%m%d_%H%M%S"),
              ".csv"
            )
          },
          content = function(file) {
            write.csv(
              data.frame(Message = "ASO generation disabled for this run."),
              file,
              row.names = FALSE
            )
          }
        )
        
        shinyjs::disable("toggle_cols_patient")
        
        showNotification(
          "Patient input summary was generated. ASO generation was skipped.",
          type = "message",
          duration = 8
        )
        
        patient_reference_results_raw(tibble())
        patient_snv_results_raw(tibble())
        patient_offtargets_total(tibble())
        current_offtargets_patient(NULL)
        selected_target_patient(NULL)
        rnaseh_stored_patient(NULL)
        
        output$patient_rnaseh_title <- renderText("No patient ASO selected")
        output$patient_rnaseh_info <- renderUI(NULL)
        output$patient_rnaseh_results <- renderDT({
          datatable(
            data.frame(Message = "ASO generation was disabled for this run."),
            rownames = FALSE,
            options = list(dom = "t"),
            class = "compact stripe"
          )
        })
        output$patient_cleavage_visual <- renderUI(NULL)
        
        output$patient_offtarget_title <- renderText("No patient ASO selected")
        output$patient_aso_seq <- renderText("")
        output$patient_numb_offtargets <- renderText("ASO generation was disabled for this run.")
        output$patient_offtarget_results <- renderDT({
          datatable(
            data.frame(Message = "ASO generation was disabled for this run."),
            rownames = FALSE,
            options = list(dom = "t"),
            class = "compact stripe"
          )
        })
        
        return(NULL)
      }
      
      # -------------------------------------------------------------
      # DUPLICATED GENERAL WORKFLOW FOR PATIENT TAB
      # reference gene only, no patient-specific ASOs yet
      # -------------------------------------------------------------
      
      incProgress(0.35, detail = "Loading reference-gene ASO workflow")
      
      print("patient milestone 5.5: workflow started")
      
      if (isTRUE(input$patient_linux_input)) {
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
        
        RNAduplex_R = function(seqs) {
          sys_cmd = system('RNAduplex',
                           input = c(seqs, seqs),
                           intern = TRUE)
          as.numeric(regmatches(sys_cmd, regexpr("-?\\d+\\.\\d+", sys_cmd)))
        }
        
        RNAselffold_R = function (seqs) {
          output = system('RNAfold --noPS',
                          input = c(seqs),
                          intern = TRUE)
          output = unlist(strsplit(output[grepl('[0-9]', output)], '[(]'))
          as.double(gsub(' |[)]', '', output[grepl('[0-9]', output)]))
        }
      }
      
      incProgress(0.40, detail = "Loading sequences")
      
      HS <- getSeq(Hsapiens, gdb_hsa)
      
      RNA_target <- HS[names(HS) == ensembl_ID]
      
      if (length(RNA_target) != 1) {
        showNotification("Selected reference gene was not found.", type = "error", duration = 8)
        return(NULL)
      }
      
      chr_coord <- list(
        chr = as.character(seqnames(target_ranges)),
        start = as.integer(start(target_ranges)),
        end = as.integer(end(target_ranges)),
        strand = ifelse(as.character(strand(target_ranges)) == "+", 1L, -1L)
      )
      
      print("patient milestone 6: loaded human gene sequences. RNA target sequence and coordinates")
      
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
      
      target_annotation <- tibble()
      unfiltered_total_data_patient <- tibble()
      
      target_annotation_patient_snv <- tibble()
      target_annotation_patient_ambiguous <- tibble()
      
      if (run_reference_patient_asos) {
      
      l <- width(RNA_target)
      oligo_lengths <- input$patient_oligo_length_range[1]:input$patient_oligo_length_range[2]
      
      target_annotation <- lapply(oligo_lengths, function(i) {
        tibble(start = 1:(l - i + 1), length = i)
      }) %>%
        bind_rows() %>%
        mutate(end = start + length - 1)
      
      target_regions <- DNAStringSet(
        RNA_target[[1]],
        start = target_annotation$start,
        width = target_annotation$length
      )
      names(target_regions) = as.character(target_regions)
      target_annotation$name = names(target_regions)
      
      target_annotation$PM_tot_freq <- NA_real_
      target_annotation$PM_max_freq <- NA_real_
      target_annotation$PM_count    <- NA_real_
      
      target_annotation <- add_empty_reference_columns(target_annotation)
      
      incProgress(0.45, detail = "Calculating accessibility and core scores")
      
      if (isTRUE(input$patient_linux_input)) {
        accessibility <- RNAplfold_R(RNA_target, u.in = max(oligo_lengths)) %>%
          as_tibble() %>%
          mutate(end = 1:width(RNA_target)) %>%
          gather(length, accessibility, -end) %>%
          mutate(length = as.double(length))
        
        target_annotation <- left_join(
          target_annotation,
          accessibility,
          by = c("length", "end")
        )
      } else {
        target_annotation$accessibility <- NA_real_
      }
      
      nucleobase_seq <- reverseComplement(target_regions)
      target_annotation$oligo_seq <- as.character(nucleobase_seq[target_annotation$name])
      
      target_annotation$tox_score = calculate_acute_neurotox(target_annotation$oligo_seq)
      
      oligo_dna <- DNAStringSet(target_annotation$oligo_seq)
      gc_counts <- letterFrequency(oligo_dna, c("G", "C"))
      oligo_len <- width(oligo_dna)
      target_annotation$gc_content <- rowSums(gc_counts) / oligo_len * 100

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
      
      target_annotation$CGs = (target_annotation$length - nchar(gsub('CG', '', target_annotation$name))) / 2
      
      print("patient milestone 7: calculated accessibility, toxicity, GC content, motif score, and CpGs")
      
      incProgress(0.50, detail = "Loading BioMart data")
      
      options(timeout = 60)
      
      if (isTRUE(biomart_available())) {
        
        martHS <- tryCatch(
          useEnsembl(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "hsapiens_gene_ensembl"
          ),
          error = function(e) useEnsembl(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "hsapiens_gene_ensembl",
            mirror  = "www"
          )
        )
        
        martMM <- tryCatch(
          useEnsembl(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "mmusculus_gene_ensembl"
          ),
          error = function(e) useEnsembl(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = "mmusculus_gene_ensembl",
            mirror  = "www"
          )
        )
        
        ortho_ENS <- tryCatch(
          getBM(
            attributes = "mmusculus_homolog_ensembl_gene",
            filters = "ensembl_gene_id",
            values = ensembl_ID,
            mart = martHS,
            bmHeader = FALSE
          ),
          error = function(e) data.frame(mmusculus_homolog_ensembl_gene = character())
        )
        
        RNA_target_mouse <- tryCatch(
          DNAStringSet(
            getBM(
              attributes = c("gene_exon_intron","ensembl_gene_id"),
              filters = "ensembl_gene_id",
              values = ortho_ENS$mmusculus_homolog_ensembl_gene,
              mart = martMM
            )$gene_exon_intron
          ),
          error = function(e) DNAStringSet()
        )
        
      } else {
        martHS <- NULL
        martMM <- NULL
        ortho_ENS <- NULL
        RNA_target_mouse <- DNAStringSet()
        target_annotation$conserved_in_mmusculus <- NA
      }
      
      if (isTRUE(input$patient_polymorphism_input) && !is.null(martHS)) {
        PMs <- tryCatch(
          getBM(
            attributes = c("minor_allele_freq", "chromosome_start"),
            filters = "ensembl_gene_id",
            values = ensembl_ID,
            mart = martHS
          ) %>%
            as_tibble() %>%
            arrange(chromosome_start, desc(minor_allele_freq)) %>%
            filter(!is.na(minor_allele_freq),
                   !duplicated(chromosome_start)) %>%
            rename(chr_start = chromosome_start, PM_freq = minor_allele_freq),
          error = function(e) tibble(
            chr_start = integer(),
            PM_freq = numeric()
          )
        )
      } else {
        PMs <- tibble(
          chr_start = integer(),
          PM_freq = numeric()
        )
      }
      
      tr = target_annotation$name
      replica = table(tr)
      target_annotation$NoRepeats = as.vector(replica[tr])
      
      if (chr_coord$strand == 1L) {
        target_annotation$chr_start <- chr_coord$start + target_annotation$start - 1L
        target_annotation$chr_end   <- chr_coord$start + target_annotation$end - 1L
      } else {
        target_annotation$chr_start <- chr_coord$end - target_annotation$end + 1L
        target_annotation$chr_end   <- chr_coord$end - target_annotation$start + 1L
      }
      
      target_annotation$chr <- gene_chr
      
      target_annotation <- annotate_general_knockdown_region(
        target_annotation = target_annotation,
        txdb = txdb_hsa,
        ensembl_ID = ensembl_ID
      )
      if (isTRUE(input$patient_polymorphism_input) && nrow(PMs) > 0) {
        PM_freq <- PMs %>%
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
            PM_count = n(),
            .groups = "drop"
          )
        
        target_annotation <- left_join(
          target_annotation,
          PM_freq,
          by = c("name" = "name", "chr_start" = "chr_start_anno")
        ) %>%
          mutate(
            PM_tot_freq = coalesce(PM_tot_freq.y, PM_tot_freq.x),
            PM_max_freq = coalesce(PM_max_freq.y, PM_max_freq.x),
            PM_count    = coalesce(PM_count.y, PM_count.x)
          ) %>%
          select(
            -PM_tot_freq.x, -PM_tot_freq.y,
            -PM_max_freq.x, -PM_max_freq.y,
            -PM_count.x, -PM_count.y
          )
      }
      
      if (length(RNA_target_mouse) > 0) {
        lm = width(RNA_target_mouse)
        
        MM_tab = lapply(sort(unique(target_annotation$length)), function(i) {
          tibble(st = 1:(lm - i + 1), w = i)
        }) %>%
          bind_rows()
        
        RNAsitesMM = DNAStringSet(
          RNA_target_mouse[[1]],
          start = MM_tab$st,
          width = MM_tab$w
        )
        
        target_annotation$conserved_in_mmusculus = target_annotation$name %in% RNAsitesMM
      } else {
        target_annotation$conserved_in_mmusculus <- NA
      }
      
      if (isTRUE(input$patient_linux_input)) {
        target_annotation$sec_energy = RNAselffold_R(target_annotation$oligo_seq)
        target_annotation$duplex_energy = RNAduplex_R(target_annotation$oligo_seq)
      } else {
        target_annotation$sec_energy <- NA_real_
        target_annotation$duplex_energy <- NA_real_
      }
      
      print("patient milestone 8: retrieved BioMart polymorphism and mouse ortholog data")
      
      # -------------------------------------------------------------
      # Build separate patient-derived SNV ASO table
      # -------------------------------------------------------------
      
      print("patient milestone 9: generated patient-specific SNV ASO candidates")
  
      
      output$patient_ambiguous_results <- DT::renderDT({
        df <- target_annotation_patient_ambiguous
        
        if (is.null(df) || nrow(df) == 0) {
          return(
            datatable(
              data.frame(Message = "No ambiguous unphased multi-variant ASO windows detected."),
              rownames = FALSE,
              options = list(
                dom = "t",
                ordering = FALSE,
                searching = FALSE,
                paging = FALSE
              ),
              class = "compact stripe"
            )
          )
        }
        
        df <- df %>%
          dplyr::select(
            dplyr::any_of(c(
              "phase_status",
              "genotype",
              "variant_note",
              "ambiguity_reason",
              "length",
              "start",
              "end",
              "region_class",
              "target_transcript",
              "chr_start",
              "chr_end"
            ))
          ) %>%
          dplyr::rename(
            "Phase status" = phase_status,
            "Genotype" = genotype,
            "Variant positions" = variant_note,
            "Reason" = ambiguity_reason,
            "ASO length (nt)" = length,
            "Start position in gene" = start,
            "End position in gene" = end,
            "Target region" = region_class,
            "Target transcript(s)" = target_transcript,
            "Chromosome start pos." = chr_start,
            "Chromosome end pos." = chr_end
          )
        
        DT::datatable(
          df,
          rownames = FALSE,
          class = "compact stripe nowrap",
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            columnDefs = list(
              list(className = "dt-center", targets = "_all")
            )
          )
        )
      })
      
      
      unfiltered_total_data_patient <- target_annotation
      
      nseq_prefilter <- nrow(target_annotation)
      
      nseq_ending_G <- nseq_prefilter - nrow(
        target_annotation %>% filter(!grepl("^C", name))
      )
      
      tox_rng <- patient_tox_range()
      nseq_toxscore <- nseq_prefilter - nrow(
        target_annotation %>%
          dplyr::filter(
            !is.na(tox_score),
            tox_score >= tox_rng[1],
            tox_score <= tox_rng[2]
          )
      )
      
      gc_rng <- patient_gc_range()
      nseq_gc <- nseq_prefilter - nrow(
        target_annotation %>%
          dplyr::filter(
            !is.na(gc_content),
            gc_content >= gc_rng[1],
            gc_content <= gc_rng[2]
          )
      )
      
      nseq_pmfreq <- if (isTRUE(input$patient_polymorphism_input)) {
        pm_rng <- patient_pm_range()
        nseq_prefilter - nrow(
          target_annotation %>%
            dplyr::mutate(PM_tot_freq = tidyr::replace_na(PM_tot_freq, 0)) %>%
            dplyr::filter(PM_tot_freq >= pm_rng[1], PM_tot_freq <= pm_rng[2])
        )
      } else {
        NA_integer_
      }
      
      apply_range <- function(df, col, rng) {
        df %>%
          dplyr::filter(
            !is.na(.data[[col]]),
            .data[[col]] >= rng[1],
            .data[[col]] <= rng[2]
          )
      }
      
      ta <- target_annotation
      
      if (isTRUE(input$patient_ASO_ending_G)) {
        ta_prev <- ta
        ta <- ta %>% filter(!grepl("^C", name))
        if (nrow(ta) == 0) {
          ta <- ta_prev
          showNotification("Filter 'ASO ending with G' removed all rows; reverting to previous dataset.", type = "warning")
        }
      }
      
      if (isTRUE(input$patient_tox_input)) {
        ta_prev <- ta
        rng <- patient_tox_range()
        ta <- ta %>%
          dplyr::filter(
            !is.na(tox_score),
            tox_score >= rng[1],
            tox_score <= rng[2]
          )
        if (nrow(ta) == 0) {
          ta <- ta_prev
          showNotification("Filter 'tox_score' removed all rows; reverting to previous dataset.", type = "warning")
        }
      }
      
      if (isTRUE(input$patient_gc_input)) {
        ta_prev <- ta
        rng <- patient_gc_range()
        ta <- ta %>%
          dplyr::filter(
            !is.na(gc_content),
            gc_content >= rng[1],
            gc_content <= rng[2]
          )
        if (nrow(ta) == 0) {
          ta <- ta_prev
          showNotification("Filter 'GC content' removed all rows; reverting to previous dataset.", type = "warning")
        }
      }
      
      if (isTRUE(input$patient_polymorphism_input)) {
        ta <- ta %>%
          mutate(across(c(PM_max_freq, PM_tot_freq, PM_count), ~ tidyr::replace_na(., 0.0)))
      }
      
      if (isTRUE(input$patient_polymorphism_input) && isTRUE(input$patient_Poly_input)) {
        ta_prev <- ta
        rng <- patient_pm_range()
        ta <- apply_range(ta, "PM_tot_freq", rng)
        if (nrow(ta) == 0) {
          ta <- ta_prev
          showNotification("Filter 'PM_tot_freq' removed all rows; reverting.", type = "warning")
        }
      }
      
      if (isTRUE(input$patient_Conserved_input)) {
        ta_prev <- ta
        ta <- ta %>% filter(conserved_in_mmusculus == TRUE)
        if (nrow(ta) == 0) {
          ta <- ta_prev
          showNotification("Filter 'conserved_in_mmusculus' removed all rows; reverting to previous dataset.", type = "warning")
        }
      }
      
      target_annotation <- ta
      
      print("patient milestone 10: applied patient reference-gene ASO filters")
      
      } else {
        target_annotation <- tibble()
        unfiltered_total_data_patient <- tibble()
        
        nseq_prefilter <- 0L
        nseq_ending_G <- 0L
        nseq_toxscore <- 0L
        nseq_gc <- 0L
        nseq_pmfreq <- 0L
        
        filter_numbers_patient <- tibble(
          Filter = c(
            "Prefiltered",
            "ASO ending with G",
            "Tox score",
            "GC content",
            "PM frequency",
            "Final selection"
          ),
          `Number of ASOs erased` = c(0, 0, 0, 0, 0, 0),
          Percentage = rep("Reference ASOs disabled in SNV-only mode", 6)
        )
      }
      
      
      if (run_patient_snv_asos) {
        print("patient milestone 9: generated patient-specific SNV ASO candidates")
        
        patient_variant_rows <- build_patient_variant_aso_rows(
          target_annotation_ref = target_annotation,
          variants_region = variants_region_expanded,
          gene_strand = gene_strand,
          gene_chr = gene_chr,
          gene_start = gene_start,
          gene_end = gene_end,
          aso_lengths = input$patient_oligo_length_range[1]:input$patient_oligo_length_range[2]
        )
        
        target_annotation_patient_snv <- patient_variant_rows$resolved
        target_annotation_patient_ambiguous <- patient_variant_rows$ambiguous
      }
      
      
      # -------------------------------------------------------------
      # Annotate, score, and optionally filter patient-derived SNV ASOs
      # This must happen AFTER build_patient_variant_aso_rows()
      # -------------------------------------------------------------
      
      if (run_patient_snv_asos && nrow(target_annotation_patient_snv) > 0) {
        
        target_annotation_patient_snv <- add_empty_reference_columns(
          target_annotation_patient_snv
        )
        
        target_annotation_patient_snv <- annotate_general_knockdown_region(
          target_annotation = target_annotation_patient_snv,
          txdb = txdb_hsa,
          ensembl_ID = ensembl_ID
        )
        
        target_annotation_patient_snv$tox_score <- calculate_acute_neurotox(
          target_annotation_patient_snv$oligo_seq
        )
        
        oligo_dna_patient <- DNAStringSet(
          target_annotation_patient_snv$oligo_seq
        )
        
        gc_counts_patient <- letterFrequency(
          oligo_dna_patient,
          c("G", "C")
        )
        
        oligo_len_patient <- width(oligo_dna_patient)
        
        target_annotation_patient_snv$gc_content <-
          rowSums(gc_counts_patient) / oligo_len_patient * 100
        
        tr_patient <- target_annotation_patient_snv$name
        replica_patient <- table(tr_patient)
        
        target_annotation_patient_snv$NoRepeats <-
          as.vector(replica_patient[tr_patient])
        
        if (!"conserved_in_mmusculus" %in% names(target_annotation_patient_snv)) {
          target_annotation_patient_snv$conserved_in_mmusculus <- NA
        }
        
        target_annotation_patient_snv$target_genomic_display <- purrr::pmap_chr(
          list(
            target_annotation_patient_snv$target_genomic_seq,
            strsplit(
              target_annotation_patient_snv$target_genomic_source,
              ";",
              fixed = TRUE
            ),
            target_annotation_patient_snv$deletion_marker_after_genomic
          ),
          color_patient_alt_sequence_html
        )
        
        target_annotation_patient_snv$target_gene_display <- purrr::pmap_chr(
          list(
            target_annotation_patient_snv$target_gene_seq,
            strsplit(
              target_annotation_patient_snv$target_gene_source,
              ";",
              fixed = TRUE
            ),
            target_annotation_patient_snv$deletion_marker_after_gene
          ),
          color_patient_alt_sequence_html
        )
        
        if (isTRUE(input$patient_linux_input)) {
          target_annotation_patient_snv$sec_energy <- RNAselffold_R(
            target_annotation_patient_snv$oligo_seq
          )
          
          target_annotation_patient_snv$duplex_energy <- RNAduplex_R(
            target_annotation_patient_snv$oligo_seq
          )
        } else {
          target_annotation_patient_snv$sec_energy <- NA_real_
          target_annotation_patient_snv$duplex_energy <- NA_real_
        }
      }
      
      # Keep an empty table with the expected columns when no SNV ASOs exist
      if (!run_patient_snv_asos || nrow(target_annotation_patient_snv) == 0) {
        target_annotation_patient_snv <- tibble(
          patient_specific_aso = character(),
          allele_phasing_status = character(),
          phase_status = character(),
          genotype = character(),
          variant_overlap_count = integer(),
          variant_type_patient = character(),
          variant_note = character(),
          variant_pos = integer(),
          variant_ref = character(),
          variant_alt = character(),
          oligo_seq = character(),
          name = character(),
          target_genomic_seq = character(),
          target_gene_seq = character(),
          target_genomic_display = character(),
          target_gene_display = character(),
          target_genomic_source = character(),
          target_gene_source = character(),
          deletion_marker_after_genomic = integer(),
          deletion_marker_after_gene = integer(),
          length = integer(),
          start = integer(),
          end = integer(),
          gc_content = numeric(),
          tox_score = numeric(),
          region_class = character(),
          target_transcript = character(),
          chr_start = integer(),
          chr_end = integer(),
          NoRepeats = numeric(),
          conserved_in_mmusculus = logical(),
          sec_energy = numeric(),
          duplex_energy = numeric()
        )
      }
      
      # Apply sidebar filters to patient-derived SNV ASOs
      # This must happen AFTER gc_content, tox_score, and region columns exist
      if (isTRUE(input$patient_filter_snv_asos) &&
          nrow(target_annotation_patient_snv) > 0) {
        
        ta_snv <- target_annotation_patient_snv
        
        if (isTRUE(input$patient_ASO_ending_G)) {
          ta_prev <- ta_snv
          
          ta_snv <- ta_snv %>%
            dplyr::filter(!grepl("^C", name))
          
          if (nrow(ta_snv) == 0) {
            ta_snv <- ta_prev
            
            showNotification(
              "SNV-ASO filter 'ASO ending with G' removed all rows; reverting.",
              type = "warning"
            )
          }
        }
        
        if (isTRUE(input$patient_tox_input)) {
          ta_prev <- ta_snv
          rng <- patient_tox_range()
          
          ta_snv <- ta_snv %>%
            dplyr::filter(
              !is.na(tox_score),
              tox_score >= rng[1],
              tox_score <= rng[2]
            )
          
          if (nrow(ta_snv) == 0) {
            ta_snv <- ta_prev
            
            showNotification(
              "SNV-ASO filter 'tox score' removed all rows; reverting.",
              type = "warning"
            )
          }
        }
        
        if (isTRUE(input$patient_gc_input)) {
          ta_prev <- ta_snv
          rng <- patient_gc_range()
          
          ta_snv <- ta_snv %>%
            dplyr::filter(
              !is.na(gc_content),
              gc_content >= rng[1],
              gc_content <= rng[2]
            )
          
          if (nrow(ta_snv) == 0) {
            ta_snv <- ta_prev
            
            showNotification(
              "SNV-ASO filter 'GC content' removed all rows; reverting.",
              type = "warning"
            )
          }
        }
        
        if (isTRUE(input$patient_Conserved_input) &&
            "conserved_in_mmusculus" %in% names(ta_snv)) {
          
          ta_prev <- ta_snv
          
          ta_snv <- ta_snv %>%
            dplyr::filter(conserved_in_mmusculus == TRUE)
          
          if (nrow(ta_snv) == 0) {
            ta_snv <- ta_prev
            
            showNotification(
              "SNV-ASO filter 'conserved in mouse' removed all rows; reverting.",
              type = "warning"
            )
          }
        }
        
        target_annotation_patient_snv <- ta_snv
      }
      
      
      incProgress(0.65, detail = "Calculating  off-target ")
      
      add_pedersen_counts_patient <- function(df, HS) {
        if (is.null(df) || nrow(df) == 0) return(df)
        
        uni_tar <- df %>%
          dplyr::select(name, length) %>%
          dplyr::distinct() %>%
          split(., .$length)
        
        uni_tar_list <- future_lapply(
          X = uni_tar,
          FUN = function(X, HS) {
            dict0 <- PDict(X$name, max.mismatch = 0)
            dict1 <- PDict(X$name, max.mismatch = 1)
            
            pm <- vwhichPDict(
              pdict = dict0,
              subject = HS,
              max.mismatch = 0,
              min.mismatch = 0
            )
            X$gene_hits_pm <- tabulate(unlist(pm), nbins = nrow(X))
            
            mm1 <- vwhichPDict(
              pdict = dict1,
              subject = HS,
              max.mismatch = 1,
              min.mismatch = 1
            )
            X$gene_hits_1mm <- tabulate(unlist(mm1), nbins = nrow(X))
            
            X
          },
          HS = HS,
          future.seed = TRUE
        )
        
        uni_tar <- dplyr::bind_rows(uni_tar_list)
        
        df %>%
          dplyr::left_join(uni_tar, by = c("name", "length"))
      }
      
      if (run_reference_patient_asos && nrow(target_annotation) > 0) {
        target_annotation <- add_pedersen_counts_patient(target_annotation, HS)
      }
      
      if (run_patient_snv_asos && nrow(target_annotation_patient_snv) > 0) {
        target_annotation_patient_snv <- add_pedersen_counts_patient(target_annotation_patient_snv, HS)
      }
      
      print("patient milestone 11: calculated perfect and 1-mismatch off-target counts")
      
      
      ta_off_reference <- if (run_reference_patient_asos && nrow(target_annotation) > 0) {
        target_annotation %>%
          dplyr::select(name, length) %>%
          dplyr::distinct()
      } else {
        tibble(name = character(), length = integer())
      }
      
      ta_off_snv <- if (run_patient_snv_asos && nrow(target_annotation_patient_snv) > 0) {
        target_annotation_patient_snv %>%
          dplyr::select(name, length) %>%
          dplyr::distinct()
      } else {
        tibble(name = character(), length = integer())
      }
      
      ta_off <- dplyr::bind_rows(
        ta_off_reference,
        ta_off_snv
      ) %>%
        dplyr::mutate(
          name = toupper(trimws(as.character(name))),
          length = as.integer(length)
        ) %>%
        dplyr::filter(
          !is.na(name),
          grepl("^[ACGT]+$", name),
          !is.na(length)
        ) %>%
        dplyr::distinct(name, length)
      
      perform_offt <- nrow(ta_off) > 0 &&
        (isTRUE(input$patient_perfect_input) || isTRUE(input$patient_mismatch_input))
      
      summary_server_patient <- NULL
      off_targets_total_patient <- tibble()
      
      if (isTRUE(perform_offt)) {
        summary_server_patient <- tryCatch(
          {
            res_list <- future_lapply(
              X = seq_len(nrow(ta_off)),
              FUN = function(i) {
                seq_i <- toupper(trimws(ta_off$name[[i]]))
                len_i <- ta_off$length[[i]]
                
                df <- all_offt(seq_i, mismatches_allowed = 2)
                
                if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
                  return(NULL)
                }
                
                df$name <- seq_i
                df$length <- len_i
                df
              },
              future.seed = TRUE
            )
            
            dplyr::bind_rows(Filter(Negate(is.null), res_list)) %>%
              dplyr::mutate(
                name = toupper(trimws(as.character(name))),
                distance = mismatches + deletions + insertions
              ) %>%
              dplyr::distinct(
                name,
                gene_name,
                match_string,
                query_seq,
                .keep_all = TRUE
              )
          },
          error = function(e) {
            message("GGGenome is currently unavailable. Off-target features are disabled. ", e$message)
            NULL
          }
        )
      }
      
      if (!is.null(summary_server_patient) && nrow(summary_server_patient) > 0) {
        tmp <- tempfile(fileext = ".bgz")
        
        download.file(
          "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
          destfile = tmp,
          mode = "wb"
        )
        
        gnomad_df <- read_tsv(
          tmp,
          show_col_types = FALSE,
          col_select = c(gene, transcript, pLI, oe_lof_upper)
        )
        
        unlink(tmp)
        
        off_targets_total_patient <- dplyr::left_join(
          summary_server_patient,
          gnomad_df,
          by = c("gene_name" = "gene")
        )
        
        dist_counts <- off_targets_total_patient %>%
          dplyr::group_by(name, distance) %>%
          dplyr::summarise(n_distance = dplyr::n(), .groups = "drop") %>%
          tidyr::complete(name, distance = 0:2, fill = list(n_distance = 0)) %>%
          dplyr::mutate(distance = paste0("n_distance_", distance)) %>%
          tidyr::pivot_wider(
            names_from = distance,
            values_from = n_distance,
            values_fill = 0
          )
        
        minloeuf_maxpli <- off_targets_total_patient %>%
          dplyr::filter(distance <= 1) %>%
          dplyr::group_by(name) %>%
          dplyr::summarise(
            max_pLI = if (all(is.na(pLI))) NA_real_ else max(pLI, na.rm = TRUE),
            min_LOEUF = if (all(is.na(oe_lof_upper))) NA_real_ else min(oe_lof_upper, na.rm = TRUE),
            .groups = "drop"
          )
        
        off_summary <- dist_counts %>%
          dplyr::left_join(minloeuf_maxpli, by = "name")
        
        if (run_reference_patient_asos && nrow(target_annotation) > 0) {
          target_annotation <- target_annotation %>%
            dplyr::left_join(off_summary, by = "name") %>%
            dplyr::mutate(
              n_distance_0 = tidyr::replace_na(n_distance_0, 0L),
              n_distance_1 = tidyr::replace_na(n_distance_1, 0L),
              n_distance_2 = tidyr::replace_na(n_distance_2, 0L),
              off_target_score =
                1    * n_distance_0 +
                0.3  * n_distance_1 +
                0.02 * n_distance_2
            )
        }
        
        if (run_patient_snv_asos && nrow(target_annotation_patient_snv) > 0) {
          target_annotation_patient_snv <- target_annotation_patient_snv %>%
            dplyr::left_join(off_summary, by = "name") %>%
            dplyr::mutate(
              n_distance_0 = tidyr::replace_na(n_distance_0, 0L),
              n_distance_1 = tidyr::replace_na(n_distance_1, 0L),
              n_distance_2 = tidyr::replace_na(n_distance_2, 0L),
              off_target_score =
                1    * n_distance_0 +
                0.3  * n_distance_1 +
                0.02 * n_distance_2
            )
        }
      }
      
      if (nrow(target_annotation) > 0 && "off_target_score" %in% names(target_annotation)) {
        target_annotation <- target_annotation %>%
          dplyr::arrange(dplyr::coalesce(off_target_score, Inf))
      }
      
      if (nrow(target_annotation_patient_snv) > 0) {
        if ("off_target_score" %in% names(target_annotation_patient_snv)) {
          target_annotation_patient_snv <- target_annotation_patient_snv %>%
            dplyr::arrange(
              dplyr::coalesce(off_target_score, Inf),
              variant_note,
              start
            )
        } else {
          target_annotation_patient_snv <- target_annotation_patient_snv %>%
            dplyr::arrange(
              variant_note,
              start
            )
        }
      }
      
      patient_offtargets_total(off_targets_total_patient)
      if (run_reference_patient_asos && nrow(target_annotation) > 0) {
        target_annotation <- add_rnaseh_summary_columns(
          target_annotation,
          mod_5prime = 0,
          mod_3prime = 0
        )
      }
      patient_reference_results_raw(target_annotation)
      if (run_patient_snv_asos && nrow(target_annotation_patient_snv) > 0) {
        target_annotation_patient_snv <- add_rnaseh_summary_columns(
          target_annotation_patient_snv,
          mod_5prime = 0,
          mod_3prime = 0
        )
      }
      patient_snv_results_raw(target_annotation_patient_snv)
      
      motif_cols <- names(motif_weights)
      
      if (run_reference_patient_asos && nrow(target_annotation) == 0) {
        target_annotation <- unfiltered_total_data_patient
        showNotification(
          "No reference ASO results after filtering, exporting unfiltered reference list",
          type = "default",
          duration = NULL,
          closeButton = TRUE
        )
      }
      
      if (run_reference_patient_asos) {
        n_final_selection <- nrow(target_annotation)
        
        filter_numbers_patient <- tibble(
          Filter = c(
            "Prefiltered",
            "ASO ending with G",
            "Tox score",
            "GC content",
            "PM frequency",
            "Final selection"
          ),
          `Number of ASOs erased` = c(
            nseq_prefilter,
            nseq_ending_G,
            nseq_toxscore,
            nseq_gc,
            nseq_pmfreq,
            n_final_selection
          ),
          Percentage = c(
            "100%",
            paste0(round(100 * nseq_ending_G / nseq_prefilter), "%"),
            paste0(round(100 * nseq_toxscore / nseq_prefilter), "%"),
            paste0(round(100 * nseq_gc / nseq_prefilter), "%"),
            paste0(round(100 * nseq_pmfreq / nseq_prefilter), "%"),
            paste0(formatC(100 * n_final_selection / nseq_prefilter, format = "f", digits = 2), "%")
          )
        )
      }
      
      output$patient_unfiltered_results_table <- renderDT({
        datatable(
          filter_numbers_patient,
          rownames = FALSE,
          options = list(
            dom = "t",
            ordering = FALSE
          ),
          class = "compact stripe"
        ) %>%
          formatStyle(
            names(filter_numbers_patient),
            textAlign = "center"
          ) %>%
          formatStyle(
            "Filter",
            target = "row",
            fontWeight = styleEqual("Final selection", "bold")
          )
      })
      
      
      patient_results1_data <- reactive({
        df <- patient_results1_lookup()
        
        column_order <- c(
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
          "region_class",
          "target_transcript",
          "motif_cor_score",
          "max_pLI",
          "min_LOEUF",
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
          "sec_energy",
          "duplex_energy",
          "accessibility",
          "rnaseh_max_score",
          "rnaseh_mean_score"
        )
        
        df <- df %>%
          dplyr::select(dplyr::any_of(column_order), dplyr::everything())
        
        column_names <- c(
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
          region_class           = "Target region",
          target_transcript      = "Target transcript(s)",
          sec_energy             = "ASO self-folding energy",
          duplex_energy          = col_with_tooltip(
            "ASO duplex energy",
            "Predicted energy needed to break the binding between two ASOs. More negative values generally indicate stronger binding, which is unfavorable for ASO design."
          ),
          max_pLI                = "Max. off-target pLI",
          min_LOEUF              = "Min. off-target LOEUF",
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
          accessibility          = "Accessibility",
          rnaseh_max_score  = "Max RNase H score",
          rnaseh_mean_score = "Average RNase H score"
        )
        
        to_rename <- intersect(names(df), names(column_names))
        names(df)[match(to_rename, names(df))] <- column_names[to_rename]
        
        if ("Off-target score" %in% names(df)) {
          df <- df %>% dplyr::arrange(dplyr::coalesce(`Off-target score`, Inf))
        }
        
        df
      })
      
      patient_results1_export_data <- reactive({
        df <- patient_results1_data()
        names(df) <- strip_header_html(names(df))
        df
      })
      
      observe({
        df <- patient_results1_data()
        
        display_names <- names(df)
        plain_names <- strip_header_html(display_names)
        choice_map <- stats::setNames(display_names, plain_names)
        
        selected_value <- if ("Off-target score" %in% plain_names) {
          display_names[match("Off-target score", plain_names)]
        } else {
          display_names[1]
        }
        
        updateSelectInput(
          session,
          "patient_main_table_sort_col",
          choices = choice_map,
          selected = selected_value
        )
      })
      
      patient_results1_sorted <- reactive({
        df <- patient_results1_data()
        
        sort_col <- input$patient_main_table_sort_col
        sort_dir <- input$patient_main_table_sort_dir
        
        if (is.null(sort_col) || !(sort_col %in% names(df))) {
          return(df)
        }
        
        if (identical(sort_dir, "desc")) {
          df <- df %>%
            dplyr::arrange(dplyr::desc(is.na(.data[[sort_col]])),
                           dplyr::desc(.data[[sort_col]]))
        } else {
          df <- df %>%
            dplyr::arrange(is.na(.data[[sort_col]]),
                           .data[[sort_col]])
        }
        
        df
      })
      
      output$patient_results1 <- DT::renderDT({
        if (isTRUE(patient_snv_only_mode())) {
          return(
            datatable(
              data.frame(
                Message = "Reference-gene ASO results were disabled for this run. SNV-ASO generation only was selected."
              ),
              rownames = FALSE,
              options = list(
                dom = "t",
                ordering = FALSE,
                searching = FALSE,
                paging = FALSE
              ),
              class = "compact stripe"
            )
          )
        }
        
        df <- patient_results1_sorted()
        
        main_cols <- c(
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
        if (!show_all_cols_patient() && length(extra_idx) > 0) {
          column_defs[[1]] <- list(
            visible = FALSE,
            targets = extra_idx - 1L
          )
        }
        
        dt <- DT::datatable(
          df,
          rownames  = FALSE,
          selection = "single",
          escape    = FALSE,
          options   = list(
            dom = "tip",
            scrollX = TRUE,
            columnDefs = c(
              column_defs,
              list(list(className = "dt-center", targets = "_all"))
            )
          ),
          class = "compact stripe cell-border nowrap"
        )
        
        DT::formatRound(
          dt,
          columns = intersect(c("GC content (%)", "PM total freq.", "PM max freq."), names(df)),
          digits  = 0
        )
      })
      
      # Enable Patient reference-result buttons only when reference ASO results exist
      observe({
        if (isTRUE(patient_snv_only_mode())) {
          shinyjs::disable("toggle_cols_patient")
          shinyjs::disable("Download_unfiltered_patient_reference")
          shinyjs::disable("Download_filtered_patient_reference")
          return()
        }
        
        df <- patient_results1_sorted()
        
        if (!is.null(df) && nrow(df) > 0) {
          shinyjs::enable("toggle_cols_patient")
          shinyjs::enable("Download_unfiltered_patient_reference")
          shinyjs::enable("Download_filtered_patient_reference")
        } else {
          shinyjs::disable("toggle_cols_patient")
          shinyjs::disable("Download_unfiltered_patient_reference")
          shinyjs::disable("Download_filtered_patient_reference")
        }
      })
      
      output$patient_snv_results <- DT::renderDT({
        df <- patient_snv_results_sorted()
        
        if (is.null(df) || nrow(df) == 0) {
          return(
            datatable(
              data.frame(Message = "No patient-derived SNV ASOs detected for the current selection."),
              rownames = FALSE,
              options = list(
                dom = "t",
                ordering = FALSE,
                searching = FALSE,
                paging = FALSE
              ),
              class = "compact stripe"
            )
          )
        }
        
        main_cols <- c(
          "Phase status",
          "Genotype",
          "Variant note",
          "ASO sequence",
          "Target genomic (forward)",
          "Target gene orientation",
          "ASO length (nt)",
          "GC content (%)",
          "Acute neurotox score"
        )
        
        main_cols <- intersect(main_cols, names(df))
        main_idx  <- match(main_cols, names(df))
        extra_idx <- setdiff(seq_along(df), main_idx)
        
        column_defs <- list()
        if (!show_all_cols_patient_snv() && length(extra_idx) > 0) {
          column_defs[[1]] <- list(
            visible = FALSE,
            targets = extra_idx - 1L
          )
        }
        
        dt <- DT::datatable(
          df,
          rownames = FALSE,
          selection = "single",
          escape = FALSE,
          class = "compact stripe nowrap",
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            columnDefs = c(
              column_defs,
              list(list(className = "dt-center", targets = "_all"))
            )
          )
        )
        
        DT::formatRound(
          dt,
          columns = intersect(c("GC content (%)"), names(df)),
          digits = 0
        )
      })
      
      # Enable Patient SNV buttons only when SNV ASO results exist
      observe({
        df <- patient_snv_results_lookup()
        
        if (!is.null(df) && nrow(df) > 0) {
          shinyjs::enable("toggle_cols_patient_snv")
          shinyjs::enable("Download_patient_snv_results")
        } else {
          shinyjs::disable("toggle_cols_patient_snv")
          shinyjs::disable("Download_patient_snv_results")
        }
      })
      
      # Enable ambiguous-window download only when ambiguous results exist
      observe({
        df <- target_annotation_patient_ambiguous
        
        if (!is.null(df) && nrow(df) > 0) {
          shinyjs::enable("Download_patient_ambiguous_results")
        } else {
          shinyjs::disable("Download_patient_ambiguous_results")
        }
      })
      
      
      output$patient_offtarget_results <- DT::renderDT({
        df <- current_offtargets_patient()
        
        if (is.null(df) || nrow(df) == 0) {
          return(
            DT::datatable(
              data.frame(Message = "No off-targets available for the selected patient ASO."),
              rownames = FALSE,
              options = list(dom = "t"),
              class = "compact stripe"
            )
          )
        }
        
        df <- df %>%
          dplyr::mutate(
            match_string_display = chartr("ID", "DI", match_string),
            subject_seq_display = purrr::map2_chr(
              subject_seq,
              match_string_display,
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
              "pLI",
              "oe_lof_upper"
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
          deletions               = "Insertions",
          insertions              = "Deletions",
          offtarget_accessibility = "Accessibility (off-target)",
          pLI                     = "GnomAD pLI",
          oe_lof_upper            = "GnomAD LOEUF"
        )
        
        nm <- names(df_view)
        nm2 <- ifelse(nm %in% names(display_map), display_map[nm], nm)
        names(df_view) <- make.unique(nm2, sep = " (dup) ")
        
        DT::datatable(
          df_view,
          rownames = FALSE,
          escape = FALSE,
          options = list(scrollX = TRUE),
          class = "compact stripe nowrap"
        )
      })
      
      observeEvent(input$patient_apply_mismatch, {
        req(current_seq_patient())
        
        ot_total <- patient_offtargets_total()
        
        mm <- as.numeric(input$patient_user_mismatch)
        seq <- toupper(current_seq_patient())
        
        current_mismatch_patient(mm)
        
        if (mm %in% c(0, 1, 2)) {
          if (is.null(ot_total) || nrow(ot_total) == 0) {
            subset_df <- tibble()
          } else {
            subset_df <- ot_total %>%
              dplyr::mutate(name_lookup = toupper(trimws(as.character(name)))) %>%
              dplyr::filter(name_lookup == seq, distance <= mm) %>%
              dplyr::select(-name_lookup)
          }
          
        } else if (mm == 3) {
          subset_df <- tryCatch({
            new_res <- all_offt(seq, mismatches_allowed = 3)
            
            if (is.null(new_res) || !is.data.frame(new_res) || nrow(new_res) == 0) {
              tibble()
            } else {
              new_res$name <- seq
              new_res$length <- nchar(seq)
              
              new_res <- new_res %>%
                dplyr::mutate(
                  distance = mismatches + deletions + insertions,
                  gene_name = stringr::str_extract(line, "(?<=\\|)[^;]+")
                ) %>%
                dplyr::distinct(gene_name, match_string, query_seq, .keep_all = TRUE)
              
              tmp <- tempfile(fileext = ".bgz")
              
              download.file(
                "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
                destfile = tmp,
                mode = "wb"
              )
              
              gnomad_df <- readr::read_tsv(
                tmp,
                show_col_types = FALSE,
                col_select = c(gene, transcript, pLI, oe_lof_upper)
              )
              
              unlink(tmp)
              
              new_res %>%
                dplyr::left_join(
                  gnomad_df,
                  by = c("gene_name" = "gene")
                )
            }
          }, error = function(e) {
            showNotification(
              paste("3-mismatch GGGenome lookup failed:", conditionMessage(e)),
              type = "error",
              duration = 8
            )
            tibble()
          })
          
        } else {
          subset_df <- tibble()
        }
        
        current_offtargets_patient(subset_df)
        
        output$patient_numb_offtargets <- renderText({
          paste0("# off targets: ", nrow(subset_df))
        })
      })
      
      observeEvent(input$patient_add_mods, {
        row_data <- selected_target_patient()
        
        if (is.null(row_data)) {
          showNotification(
            "Select a patient ASO row first.",
            type = "warning",
            duration = 6
          )
          return()
        }
        
        rnaseh_data <- rnaseh_results(
          selected_row_name = row_data$name[[1]],
          oligo_seq = row_data$oligo_seq[[1]],
          mod_5prime = input$patient_mod_5prime,
          mod_3prime = input$patient_mod_3prime
        )
        
        rnaseh_stored_patient(rnaseh_data)
        
        output$patient_rnaseh_results <- DT::renderDT({
          rnaseh_view <- rename_rnaseh_cols(rnaseh_data)
          DT::datatable(
            rnaseh_view,
            selection = list(mode = "single", selected = 1),
            options = list(
              paging = FALSE,
              searching = FALSE,
              info = FALSE,
              dom = "t"
            )
          )
        })
        
        output$patient_cleavage_visual <- renderUI(div())
      })
      
      output$patient_download_offtarget <- downloadHandler(
        filename = function() {
          paste0(
            "patient_offtargets_",
            current_seq_patient(),
            "_mismatches_",
            current_mismatch_patient(),
            "_",
            Sys.Date(),
            ".csv"
          )
        },
        content = function(con) {
          data_offtarget <- current_offtargets_patient()
          req(data_offtarget)
          write.csv2(data_offtarget, con, row.names = FALSE)
        }
      )
      
      output$patient_download_rnaseh <- downloadHandler(
        filename = function() {
          row_data <- selected_target_patient()
          
          if (is.null(row_data)) {
            "patient_RNaseH_results.xlsx"
          } else {
            paste0("patient_RNaseH_results_", row_data$name[[1]], ".xlsx")
          }
        },
        content = function(file) {
          data <- rnaseh_stored_patient()
          
          if (is.null(data)) {
            data <- data.frame(Message = "No available data")
          }
          
          data[] <- lapply(data, as.character)
          
          wb <- createWorkbook()
          addWorksheet(wb, "RNaseH_results")
          writeData(wb, "RNaseH_results", data)
          saveWorkbook(wb, file, overwrite = TRUE)
        }
      )
      
      output$Download_unfiltered_patient_reference <- downloadHandler(
        filename = function() {
          paste0("patient_reference_unfiltered_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
        },
        content = function(file) {
          if (isTRUE(patient_snv_only_mode())) {
            write.csv(
              data.frame(Message = "Reference-gene ASO generation was disabled for this run."),
              file,
              row.names = FALSE
            )
          } else {
            write.csv(unfiltered_total_data_patient, file, row.names = FALSE)
          }
        }
      )
      
      output$Download_patient_ambiguous_results <- downloadHandler(
        filename = function() {
          paste0(
            "patient_ambiguous_aso_windows_",
            format(Sys.time(), "%Y%m%d_%H%M%S"),
            ".csv"
          )
        },
        content = function(file) {
          if (is.null(target_annotation_patient_ambiguous) || nrow(target_annotation_patient_ambiguous) == 0) {
            write.csv(
              data.frame(Message = "No ambiguous unphased multi-variant ASO windows detected."),
              file,
              row.names = FALSE
            )
          } else {
            write.csv(target_annotation_patient_ambiguous, file, row.names = FALSE)
          }
        }
      )
      
      
      output$Download_filtered_patient_reference <- downloadHandler(
        filename = function() {
          paste0("patient_reference_filtered_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
        },
        content = function(file) {
          if (isTRUE(patient_snv_only_mode())) {
            write.csv(
              data.frame(Message = "Reference-gene ASO generation was disabled for this run."),
              file,
              row.names = FALSE
            )
          } else {
            write.csv(patient_results1_export_data(), file, row.names = FALSE)
          }
        }
      )
      
      output$Download_patient_snv_results <- downloadHandler(
        filename = function() {
          paste0(
            "patient_snv_aso_results_",
            format(Sys.time(), "%Y%m%d_%H%M%S"),
            ".csv"
          )
        },
        content = function(file) {
          if (is.null(target_annotation_patient_snv) || nrow(target_annotation_patient_snv) == 0) {
            write.csv(
              data.frame(Message = "No patient-derived SNV ASOs detected."),
              file,
              row.names = FALSE
            )
          } else {
            write.csv(target_annotation_patient_snv, file, row.names = FALSE)
          }
        }
      )
      
      
      
      incProgress(1, detail = "Done") 
      
      print("patient milestone 13: rendered patient summary, reference ASO table, and SNV ASO table")

      showNotification(
        "Patient summary and reference-gene ASO results were generated successfully.",
        type = "message",
        duration = 8
      )
    })
  })

  
#######$$$############
  
  
#@@@#
  ###############################################################
  ###############################################################
  ###########           Specific ASO input tab     ##############
  ###############################################################
  ###############################################################
  
  single_aso_results_raw <- reactiveVal(tibble())
  single_aso_offtargets_total <- reactiveVal(tibble())
  
  selected_target_single_aso <- reactiveVal(NULL)
  oligo_sequence_single_aso <- reactiveVal(NULL)
  rnaseh_stored_single_aso <- reactiveVal(NULL)
  
  current_seq_single_aso <- reactiveVal(NULL)
  current_mismatch_single_aso <- reactiveVal(2)
  current_offtargets_single_aso <- reactiveVal(NULL)
  cached_results_single_aso <- reactiveVal(list())
  
  show_all_cols_single_aso <- reactiveVal(FALSE)
  
  # Disable Specific ASO buttons until results exist/selection exists
  observe({
    shinyjs::disable("toggle_cols_single_aso")
    shinyjs::disable("Download_single_aso_results")
    shinyjs::disable("single_aso_download_rnaseh")
    shinyjs::disable("single_aso_download_offtarget")
  })
  
  observeEvent(input$toggle_cols_single_aso, {
    current <- show_all_cols_single_aso()
    show_all_cols_single_aso(!current)
    
    if (current) {
      updateActionButton(
        session,
        inputId = "toggle_cols_single_aso",
        label = "Extended data"
      )
    } else {
      updateActionButton(
        session,
        inputId = "toggle_cols_single_aso",
        label = "Hide extended data"
      )
    }
  })
  
  observeEvent(input$run_button_single_aso, {
    
    shinyjs::disable("toggle_cols_single_aso")
    shinyjs::disable("Download_single_aso_results")
    shinyjs::disable("single_aso_download_rnaseh")
    shinyjs::disable("single_aso_download_offtarget")
    
    selected_target_single_aso(NULL)
    rnaseh_stored_single_aso(NULL)
    current_offtargets_single_aso(NULL)
    
    aso_matches <- parsed_asos()
    
    if (length(aso_matches) == 0) {
      showNotification(
        "No valid ASO sequences detected. Please enter sequences containing only A, C, G, and T.",
        type = "error",
        duration = 8
      )
      return(NULL)
    }
    
    single_mode <- input$single_aso_analysis_mode
    
    use_single_reference <- single_aso_uses_reference(single_mode)
    use_single_full_reference <- single_aso_full_reference_mode(single_mode)
    use_single_consensus_cf <- single_aso_consensus_cf_mode(single_mode)
    
    if (isTRUE(use_single_reference)) {
      req(input$ensemble_id_input_single_aso)
    }
    
    if (isTRUE(use_single_consensus_cf)) {
      req(input$single_aso_variant_file)
    }
    
    start_time <- Sys.time()
    
    showNotification(
      "Specific ASO analysis started",
      type = "default",
      duration = 10,
      closeButton = TRUE
    )
    
    withProgress(message = "Running specific ASO analysis...", value = 0, {
      
      incProgress(0.05, detail = "Preparing ASO input")
      
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
      
      target_annotation <- build_target_annotation_from_asos(aso_matches)
      
      if (nrow(target_annotation) == 0) {
        showNotification(
          "No valid ASO sequences detected after parsing.",
          type = "error",
          duration = 8
        )
        return(NULL)
      }
      
      target_annotation <- add_empty_reference_columns(target_annotation)
      
      target_annotation$accessibility <- NA_real_
      target_annotation$sec_energy <- NA_real_
      target_annotation$duplex_energy <- NA_real_
      target_annotation$conserved_in_mmusculus <- NA
      target_annotation$NoRepeats <- NA_integer_
      
      incProgress(0.15, detail = "Calculating ASO scores")
      
      target_annotation$tox_score <- calculate_acute_neurotox(
        target_annotation$oligo_seq
      )
      
      oligo_dna <- DNAStringSet(target_annotation$oligo_seq)
      
      gc_counts <- letterFrequency(
        oligo_dna,
        c("G", "C")
      )
      
      oligo_len <- width(oligo_dna)
      
      target_annotation$gc_content <- rowSums(gc_counts) / oligo_len * 100
      
      seqs <- DNAStringSet(target_annotation$oligo_seq)
      
      motif_counts_mat <- sapply(
        names(motif_weights),
        function(m) vcountPattern(m, seqs),
        simplify = "matrix"
      )
      
      if (nrow(target_annotation) == 1) {
        motif_counts_mat <- t(motif_counts_mat)
      }
      
      motif_counts_df <- as.data.frame(
        motif_counts_mat,
        stringsAsFactors = FALSE
      )
      
      motif_scores <- as.matrix(motif_counts_mat) %*% as.numeric(motif_weights)
      motif_scores <- round(motif_scores[, 1], 6)
      
      target_annotation <- bind_cols(
        target_annotation,
        motif_counts_df,
        motif_cor_score = motif_scores
      )
      
      target_annotation$CGs <-
        (target_annotation$length - nchar(gsub("CG", "", target_annotation$name))) / 2
      
      incProgress(0.30, detail = "Loading human transcriptome")
      
      txdb_hsa_single <- tryCatch({
        loadDb("/opt/ERASOR/txdb_hsa_biomart.db")
      }, error = function(e) {
        loadDb("../txdb_hsa_biomart.db")
      })
      
      gdb_hsa_single <- genes(txdb_hsa_single)
      seqlevelsStyle(gdb_hsa_single) <- seqlevelsStyle(Hsapiens)
      
      common_chrs_single <- intersect(
        seqlevels(gdb_hsa_single),
        seqlevels(Hsapiens)
      )
      
      gdb_hsa_single <- keepSeqlevels(
        gdb_hsa_single,
        common_chrs_single,
        pruning.mode = "coarse"
      )
      
      HS_single <- getSeq(Hsapiens, gdb_hsa_single)
      
      # -------------------------------------------------------------
      # Optional reference-gene or patient-consensus annotation
      # -------------------------------------------------------------
      
      if (isTRUE(use_single_reference)) {
        
        ensembl_ID_single <- input$ensemble_id_input_single_aso
        
        target_ranges_single <- gdb_hsa_single[names(gdb_hsa_single) == ensembl_ID_single]
        
        if (length(target_ranges_single) != 1) {
          showNotification(
            "Selected Specific ASO reference gene could not be uniquely identified.",
            type = "error",
            duration = NULL
          )
          return(NULL)
        }
        
        gene_chr_single <- normalize_chr_style(as.character(seqnames(target_ranges_single)))
        gene_start_single <- as.integer(start(target_ranges_single))
        gene_end_single <- as.integer(end(target_ranges_single))
        gene_strand_single <- as.character(strand(target_ranges_single))
        
        RNA_target_single <- HS_single[names(HS_single) == ensembl_ID_single]
        
        if (length(RNA_target_single) != 1) {
          showNotification(
            "Selected Specific ASO reference gene sequence was not found.",
            type = "error",
            duration = NULL
          )
          return(NULL)
        }
        
        chr_coord_single <- list(
          chr = as.character(seqnames(target_ranges_single)),
          start = gene_start_single,
          end = gene_end_single,
          strand = ifelse(gene_strand_single == "+", 1L, -1L)
        )
        
        martHS_single <- NULL
        martMM_single <- NULL
        
        if (isTRUE(biomart_available())) {
          martHS_single <- tryCatch(
            useEnsembl(
              biomart = "ENSEMBL_MART_ENSEMBL",
              dataset = "hsapiens_gene_ensembl"
            ),
            error = function(e) useEnsembl(
              biomart = "ENSEMBL_MART_ENSEMBL",
              dataset = "hsapiens_gene_ensembl",
              mirror  = "www"
            )
          )
          
          martMM_single <- tryCatch(
            useEnsembl(
              biomart = "ENSEMBL_MART_ENSEMBL",
              dataset = "mmusculus_gene_ensembl"
            ),
            error = function(e) useEnsembl(
              biomart = "ENSEMBL_MART_ENSEMBL",
              dataset = "mmusculus_gene_ensembl",
              mirror  = "www"
            )
          )
        }
        
        if (isTRUE(use_single_full_reference)) {
          
          incProgress(0.36, detail = "Mapping ASOs to selected reference gene")
          
          target_annotation <- locate_specific_asos_in_reference_gene(
            target_annotation = target_annotation,
            RNA_target = RNA_target_single,
            chr_coord = chr_coord_single
          )
          
          target_annotation$chr <- gene_chr_single
          
          target_annotation <- add_specific_aso_reference_cf_annotations(
            target_annotation = target_annotation,
            txdb_hsa = txdb_hsa_single,
            ensembl_ID = ensembl_ID_single,
            martHS = martHS_single,
            martMM = martMM_single,
            input_polymorphism = isTRUE(input$single_aso_polymorphism_input)
          )
          
          # If polymorphism lookup was requested and BioMart was available,
          # missing PM values mean no overlapping polymorphism was found.
          if (
            isTRUE(input$single_aso_polymorphism_input) &&
            !is.null(martHS_single)
          ) {
            target_annotation <- target_annotation %>%
              dplyr::mutate(
                PM_tot_freq = tidyr::replace_na(PM_tot_freq, 0),
                PM_max_freq = tidyr::replace_na(PM_max_freq, 0),
                PM_count    = tidyr::replace_na(PM_count, 0)
              )
          }
          
          # Number of repeats inside selected reference gene
          if ("reference_gene_match_count" %in% names(target_annotation)) {
            target_annotation$NoRepeats <- target_annotation$reference_gene_match_count
          }
          
        } else if (isTRUE(use_single_consensus_cf)) {
          
          incProgress(0.36, detail = "Building patient consensus sequences")
          
          region_string_single <- paste0(
            gene_chr_single,
            ":",
            gene_start_single,
            "-",
            gene_end_single
          )
          
          staged_single <- tryCatch(
            stage_patient_variant_files(
              variant_input = input$single_aso_variant_file,
              index_input = input$single_aso_variant_index_file
            ),
            error = function(e) {
              showNotification(
                paste("Specific ASO patient file preparation failed:", conditionMessage(e)),
                type = "error",
                duration = NULL
              )
              return(NULL)
            }
          )
          
          if (is.null(staged_single)) {
            return(NULL)
          }
          
          variants_raw_single <- tryCatch(
            read_patient_variants_resilient(
              variant_path = staged_single$variant_path,
              region_string = region_string_single,
              is_indexed = staged_single$is_indexed
            ),
            error = function(e) {
              showNotification(
                paste("Specific ASO variant reading failed:", conditionMessage(e)),
                type = "error",
                duration = NULL
              )
              return(NULL)
            }
          )
          
          if (is.null(variants_raw_single)) {
            return(NULL)
          }
          
          variants_region_single <- variants_raw_single %>%
            dplyr::mutate(
              chr_original = chr,
              chr_norm_chr = normalize_chr_style(chr),
              chr_norm_nochr = normalize_chr_style_nochr(chr)
            ) %>%
            dplyr::filter(
              chr_norm_chr == normalize_chr_style(gene_chr_single) |
                chr_norm_nochr == normalize_chr_style_nochr(gene_chr_single),
              pos >= gene_start_single,
              pos <= gene_end_single,
              !is.na(gt),
              !gt %in% c("0|0", "0/0", "./.", ".|.")
            ) %>%
            dplyr::arrange(pos) %>%
            dplyr::select(-chr_norm_chr, -chr_norm_nochr)
          
          variants_region_single_expanded <- expand_multiallelic_patient_variants(
            variants_region_single
          )
          
          variants_region_single_expanded <- validate_variants_against_reference(
            variants_region_single_expanded
          )
          
          bad_ref_single <- variants_region_single_expanded %>%
            dplyr::filter(!ref_matches_genome)
          
          if (nrow(bad_ref_single) > 0) {
            showNotification(
              "Specific ASO consensus build stopped because some VCF REF alleles do not match the genome used by the app.",
              type = "error",
              duration = NULL
            )
            print(
              bad_ref_single %>%
                dplyr::select(chr, pos, ref, alt, ref_observed_genome) %>%
                head(20)
            )
            return(NULL)
          }
          
          ref_seq_single_genomic_forward <- get_seq_one_based(
            chr = gene_chr_single,
            start_pos = gene_start_single,
            end_pos = gene_end_single
          )
          
          consensus_list_single <- build_patient_consensus_set(
            ref_seq = ref_seq_single_genomic_forward,
            region_chr = gene_chr_single,
            region_start = gene_start_single,
            variants_df = variants_region_single_expanded,
            gene_strand = gene_strand_single,
            max_ambiguous_combinations = as.integer(input$single_aso_consensus_max_ambiguous)
          )
          
          target_annotation <- map_specific_asos_to_consensus_set(
            target_annotation = target_annotation,
            consensus_list = consensus_list_single,
            gene_start = gene_start_single,
            gene_end = gene_end_single,
            gene_strand = gene_strand_single,
            gene_chr = gene_chr_single
          )
          
          target_annotation$chr <- gene_chr_single
          
          target_annotation <- add_specific_aso_reference_cf_annotations(
            target_annotation = target_annotation,
            txdb_hsa = txdb_hsa_single,
            ensembl_ID = ensembl_ID_single,
            martHS = martHS_single,
            martMM = martMM_single,
            input_polymorphism = isTRUE(input$single_aso_polymorphism_input)
          )
          
          # If polymorphism lookup was requested and BioMart was available,
          # missing PM values mean no overlapping polymorphism was found.
          if (
            isTRUE(input$single_aso_polymorphism_input) &&
            !is.null(martHS_single)
          ) {
            target_annotation <- target_annotation %>%
              dplyr::mutate(
                PM_tot_freq = tidyr::replace_na(PM_tot_freq, 0),
                PM_max_freq = tidyr::replace_na(PM_max_freq, 0),
                PM_count    = tidyr::replace_na(PM_count, 0)
              )
          }
          
          # In consensus-CF mode, do not display/use full-reference-only features.
          target_annotation$accessibility <- NA_real_
          target_annotation$NoRepeats <- NA_real_
        }
      }
      
      incProgress(0.45, detail = "Calculating Pedersen off-target counts")
      
      uni_tar <- target_annotation %>%
        dplyr::select(name, length) %>%
        dplyr::distinct() %>%
        split(., .$length)
      
      uni_tar_list <- future_lapply(
        X = uni_tar,
        FUN = function(X, HS) {
          dict0 <- PDict(X$name, max.mismatch = 0)
          dict1 <- PDict(X$name, max.mismatch = 1)
          
          pm <- vwhichPDict(
            pdict = dict0,
            subject = HS,
            max.mismatch = 0,
            min.mismatch = 0
          )
          
          X$gene_hits_pm <- tabulate(unlist(pm), nbins = nrow(X))
          
          mm1 <- vwhichPDict(
            pdict = dict1,
            subject = HS,
            max.mismatch = 1,
            min.mismatch = 1
          )
          
          X$gene_hits_1mm <- tabulate(unlist(mm1), nbins = nrow(X))
          
          X
        },
        HS = HS_single,
        future.seed = TRUE
      )
      
      uni_tar <- dplyr::bind_rows(uni_tar_list)
      
      target_annotation <- target_annotation %>%
        dplyr::left_join(uni_tar, by = c("name", "length"))
      
      incProgress(0.65, detail = "Searching GGGenome off-targets")
      
      ta_off <- target_annotation %>%
        dplyr::select(name, length) %>%
        dplyr::distinct()
      
      summary_server_single <- tryCatch(
        {
          res_list <- future_lapply(
            X = seq_len(nrow(ta_off)),
            FUN = function(i) {
              seq_i <- toupper(trimws(ta_off$name[[i]]))
              len_i <- ta_off$length[[i]]
              
              df <- all_offt(seq_i, mismatches_allowed = 2)
              
              if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
                return(NULL)
              }
              
              df$name <- seq_i
              df$length <- len_i
              df
            },
            future.seed = TRUE
          )
          
          dplyr::bind_rows(Filter(Negate(is.null), res_list)) %>%
            dplyr::mutate(
              name = toupper(trimws(as.character(name))),
              distance = mismatches + deletions + insertions
            ) %>%
            dplyr::distinct(
              name,
              gene_name,
              match_string,
              query_seq,
              .keep_all = TRUE
            )
        },
        error = function(e) {
          message("GGGenome failed for Specific ASO input: ", e$message)
          NULL
        }
      )
      
      off_targets_total_single <- tibble()
      
      if (!is.null(summary_server_single) && nrow(summary_server_single) > 0) {
        
        incProgress(0.80, detail = "Adding gnomAD annotation")
        
        tmp <- tempfile(fileext = ".bgz")
        
        download.file(
          "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
          destfile = tmp,
          mode = "wb"
        )
        
        gnomad_df <- read_tsv(
          tmp,
          show_col_types = FALSE,
          col_select = c(gene, transcript, pLI, oe_lof_upper)
        )
        
        unlink(tmp)
        
        off_targets_total_single <- dplyr::left_join(
          summary_server_single,
          gnomad_df,
          by = c("gene_name" = "gene")
        )
        
        dist_counts <- off_targets_total_single %>%
          dplyr::group_by(name, distance) %>%
          dplyr::summarise(n_distance = dplyr::n(), .groups = "drop") %>%
          tidyr::complete(name, distance = 0:2, fill = list(n_distance = 0)) %>%
          dplyr::mutate(distance = paste0("n_distance_", distance)) %>%
          tidyr::pivot_wider(
            names_from = distance,
            values_from = n_distance,
            values_fill = 0
          )
        
        minloeuf_maxpli <- off_targets_total_single %>%
          dplyr::filter(distance <= 1) %>%
          dplyr::group_by(name) %>%
          dplyr::summarise(
            max_pLI = if (all(is.na(pLI))) NA_real_ else max(pLI, na.rm = TRUE),
            min_LOEUF = if (all(is.na(oe_lof_upper))) NA_real_ else min(oe_lof_upper, na.rm = TRUE),
            .groups = "drop"
          )
        
        off_summary <- dist_counts %>%
          dplyr::left_join(minloeuf_maxpli, by = "name")
        
        target_annotation <- target_annotation %>%
          dplyr::left_join(off_summary, by = "name") %>%
          dplyr::mutate(
            n_distance_0 = tidyr::replace_na(n_distance_0, 0L),
            n_distance_1 = tidyr::replace_na(n_distance_1, 0L),
            n_distance_2 = tidyr::replace_na(n_distance_2, 0L),
            off_target_score =
              1    * n_distance_0 +
              0.3  * n_distance_1 +
              0.02 * n_distance_2
          )
      }
      
      target_annotation <- target_annotation %>%
        dplyr::arrange(input_order)
      
      # -------------------------------------------------------------
      # Specific ASO: ViennaRNA self-folding and duplex energies
      # -------------------------------------------------------------
      if (isTRUE(input$single_aso_linux_input)) {
        
        RNAduplex_R <- function(seqs) {
          sys_cmd <- system(
            "RNAduplex",
            input = c(seqs, seqs),
            intern = TRUE
          )
          
          as.numeric(
            regmatches(
              sys_cmd,
              regexpr("-?\\d+\\.\\d+", sys_cmd)
            )
          )
        }
        
        RNAselffold_R <- function(seqs) {
          output <- system(
            "RNAfold --noPS",
            input = c(seqs),
            intern = TRUE
          )
          
          output <- unlist(
            strsplit(
              output[grepl("[0-9]", output)],
              "[(]"
            )
          )
          
          as.double(
            gsub(
              " |[)]",
              "",
              output[grepl("[0-9]", output)]
            )
          )
        }
        
        target_annotation$sec_energy <- RNAselffold_R(
          target_annotation$oligo_seq
        )
        
        target_annotation$duplex_energy <- RNAduplex_R(
          target_annotation$oligo_seq
        )
        
      } else {
        target_annotation$sec_energy <- NA_real_
        target_annotation$duplex_energy <- NA_real_
      }
      
      motif_cols <- names(motif_weights)
      target_annotation <- target_annotation %>%
        dplyr::select(-dplyr::any_of(motif_cols))
      
      target_annotation <- add_rnaseh_summary_columns(
        target_annotation,
        mod_5prime = 0,
        mod_3prime = 0
      )
      
      single_aso_results_raw(target_annotation)
      single_aso_offtargets_total(off_targets_total_single)
      
      incProgress(1, detail = "Done")
    })
    
    end_time <- Sys.time()
    elapsed_str <- format_elapsed(start_time, end_time)
    
    showNotification(
      paste("Specific ASO analysis completed in", elapsed_str),
      type = "message",
      duration = NULL,
      closeButton = TRUE
    )
  })
  
  
  single_aso_results_lookup <- reactive({
    df <- single_aso_results_raw()
    
    if (is.null(df) || nrow(df) == 0) {
      return(tibble())
    }
    
    if ("input_order" %in% names(df)) {
      df <- df %>% dplyr::arrange(input_order)
    }
    
    df
  })
  
  single_aso_results_data <- reactive({
    df <- single_aso_results_lookup()
    
    if (is.null(df) || nrow(df) == 0) {
      return(tibble())
    }
    
    column_order <- c(
      "input_order",
      "oligo_seq",
      "name",
      "length",
      
      "reference_gene_match_status",
      "reference_gene_match_count",
      
      "consensus_match_status",
      "consensus_source",
      "consensus_match_count",
      "consensus_variant_overlap_count",
      "consensus_variant_specific_match",
      "consensus_variant_notes",
      "reference_projected_target",
      "reference_projected_match",
      
      "start",
      "end",
      "region_class",
      "target_transcript",
      "chr_start",
      "chr_end",
      
      "gc_content",
      "tox_score",
      "off_target_score",
      "n_distance_0",
      "n_distance_1",
      "n_distance_2",
      "motif_cor_score",
      "max_pLI",
      "min_LOEUF",
      "conserved_in_mmusculus",
      "mouse_conservation_status",
      "mouse_ortholog_ids",
      "CGs",
      "PM_tot_freq",
      "PM_max_freq",
      "PM_count",
      "gene_hits_pm",
      "gene_hits_1mm",
      "NoRepeats",
      "sec_energy",
      "duplex_energy",
      "accessibility",
      "rnaseh_max_score",
      "rnaseh_mean_score"
    )
    
    df <- df %>%
      dplyr::select(dplyr::any_of(column_order), dplyr::everything())
    
    column_names <- c(
      input_order                         = "Input order",
      oligo_seq                           = "ASO sequence",
      name                                = "Target (DNA)",
      length                              = "ASO length (nt)",
      
      reference_gene_match_status         = "Reference gene match status",
      reference_gene_match_count          = "Reference gene match count",
      
      consensus_match_status              = "Consensus match status",
      consensus_source                    = "Consensus source",
      consensus_match_count               = "Consensus match count",
      consensus_variant_overlap_count     = "Consensus variant-overlap count",
      consensus_variant_specific_match    = "Variant-specific consensus match",
      consensus_variant_notes             = "Consensus variant notes",
      reference_projected_target              = "Reference sequence at projected locus",
      reference_projected_match               = "Reference sequence equals target",
      
      start                               = "Start position in gene",
      end                                 = "End position in gene",
      region_class                        = "Target region",
      target_transcript                   = "Target transcript(s)",
      chr_start                           = "Chromosome start pos.",
      chr_end                             = "Chromosome end pos.",
      
      gc_content                          = "GC content (%)",
      tox_score                           = "Acute neurotox score",
      off_target_score                    = "Off-target score",
      n_distance_0                        = "Perfect matches (GGGenome)",
      n_distance_1                        = "Off-targets 1 mismatch (GGGenome)",
      n_distance_2                        = "Off-targets 2 mismatches (GGGenome)",
      motif_cor_score                     = "Motif correlation score",
      max_pLI                             = "Max. off-target pLI",
      min_LOEUF                           = "Min. off-target LOEUF",
      conserved_in_mmusculus              = "Conserved in mouse",
      mouse_conservation_status           = "Mouse conservation status",
      mouse_ortholog_ids                  = "Mouse ortholog ID(s)",
      CGs                                 = "Number of CpGs",
      PM_tot_freq                         = "PM total freq.",
      PM_max_freq                         = "PM max freq.",
      PM_count                            = "PM count",
      gene_hits_pm                        = "Perfect matches (Pedersen)",
      gene_hits_1mm                       = "Off-targets 1 mismatch (Pedersen)",
      NoRepeats                           = "Number of repeats",
      sec_energy                          = "ASO self-folding energy",
      duplex_energy                       = "ASO duplex energy",
      accessibility                       = "Accessibility",
      rnaseh_max_score                    = "Max RNase H score",
      rnaseh_mean_score                   = "Average RNase H score"
    )
    
    to_rename <- intersect(names(df), names(column_names))
    names(df)[match(to_rename, names(df))] <- column_names[to_rename]
    
    single_mode <- input$single_aso_analysis_mode
    
    if (identical(single_mode, "no_reference")) {
      df <- df %>%
        dplyr::select(-dplyr::any_of(c(
          "Reference gene match status",
          "Reference gene match count",
          "Consensus match status",
          "Consensus source",
          "Consensus match count",
          "Consensus variant-overlap count",
          "Consensus variant notes",
          "Start position in gene",
          "End position in gene",
          "Target region",
          "Target transcript(s)",
          "Chromosome start pos.",
          "Chromosome end pos.",
          "PM total freq.",
          "PM max freq.",
          "PM count",
          "Number of repeats",
          "Conserved in mouse",
          "Mouse conservation status",
          "Accessibility",
          "Max RNase H score",
          "Average RNase H score"
        )))
    }
    
    if (identical(single_mode, "with_reference_full")) {
      df <- df %>%
        dplyr::select(-dplyr::any_of(c(
          "Consensus match status",
          "Consensus source",
          "Consensus match count",
          "Consensus variant-overlap count",
          "Consensus variant notes"
        )))
    }
    
    if (identical(single_mode, "with_reference_consensus_cf")) {
      df <- df %>%
        dplyr::select(-dplyr::any_of(c(
          "Reference gene match status",
          "Reference gene match count",
          "Target region",
          "Target transcript(s)",
          "Number of repeats",
          "Accessibility"
        )))
    }
    
    df
  })
  
  single_aso_results_export_data <- reactive({
    df <- single_aso_results_data()
    names(df) <- strip_header_html(names(df))
    df
  })
  
  observe({
    df <- single_aso_results_data()
    
    if (is.null(df) || nrow(df) == 0) {
      return()
    }
    
    display_names <- names(df)
    plain_names <- strip_header_html(display_names)
    
    choice_map <- stats::setNames(display_names, plain_names)
    
    selected_value <- if ("Input order" %in% plain_names) {
      display_names[match("Input order", plain_names)]
    } else {
      display_names[1]
    }
    
    updateSelectInput(
      session,
      "single_aso_table_sort_col",
      choices = choice_map,
      selected = selected_value
    )
  })
  
  single_aso_results_sorted <- reactive({
    df <- single_aso_results_data()
    
    if (is.null(df) || nrow(df) == 0) {
      return(tibble())
    }
    
    sort_col <- input$single_aso_table_sort_col
    sort_dir <- input$single_aso_table_sort_dir
    
    if (is.null(sort_col) || !(sort_col %in% names(df))) {
      return(df)
    }
    
    if (identical(sort_dir, "desc")) {
      df <- df %>%
        dplyr::arrange(
          dplyr::desc(is.na(.data[[sort_col]])),
          dplyr::desc(.data[[sort_col]])
        )
    } else {
      df <- df %>%
        dplyr::arrange(
          is.na(.data[[sort_col]]),
          .data[[sort_col]]
        )
    }
    
    df
  })
  
  # Enable Specific ASO result buttons only when results exist
  observe({
    df <- single_aso_results_sorted()
    
    if (!is.null(df) && nrow(df) > 0) {
      shinyjs::enable("toggle_cols_single_aso")
      shinyjs::enable("Download_single_aso_results")
    } else {
      shinyjs::disable("toggle_cols_single_aso")
      shinyjs::disable("Download_single_aso_results")
    }
  })
  
  output$single_aso_results <- DT::renderDT({
    df <- single_aso_results_sorted()
    
    if (is.null(df) || nrow(df) == 0) {
      return(
        DT::datatable(
          data.frame(Message = "No Specific ASO results yet."),
          rownames = FALSE,
          options = list(dom = "t"),
          class = "compact stripe"
        )
      )
    }
    
    main_cols <- c(
      "Input order",
      "ASO sequence",
      "Target (DNA)",
      "ASO length (nt)",
      
      if (identical(input$single_aso_analysis_mode, "with_reference_full")) {
        c(
          "Reference gene match status",
          "Reference gene match count",
          "Start position in gene",
          "End position in gene"
        )
      },
      
      if (identical(input$single_aso_analysis_mode, "with_reference_consensus_cf")) {
        c(
          "Consensus match status",
          "Consensus source",
          "Consensus match count",
          "Consensus variant-overlap count",
          "Start position in gene",
          "End position in gene"
        )
      },
      
      "GC content (%)",
      "Acute neurotox score",
      
      if (!identical(input$single_aso_analysis_mode, "no_reference")) {
        c(
          "Conserved in mouse",
          "PM total freq.",
          "PM max freq.",
          "PM count"
        )
      },
      
      "Off-target score",
      "Perfect matches (GGGenome)",
      "Off-targets 1 mismatch (GGGenome)",
      "Off-targets 2 mismatches (GGGenome)"
    )
    
    main_cols <- intersect(main_cols, names(df))
    main_idx <- match(main_cols, names(df))
    extra_idx <- setdiff(seq_along(df), main_idx)
    
    column_defs <- list()
    
    if (!show_all_cols_single_aso() && length(extra_idx) > 0) {
      column_defs[[1]] <- list(
        visible = FALSE,
        targets = extra_idx - 1L
      )
    }
    
    dt <- DT::datatable(
      df,
      rownames = FALSE,
      selection = "single",
      escape = FALSE,
      options = list(
        dom = "tip",
        scrollX = TRUE,
        columnDefs = c(
          column_defs,
          list(list(className = "dt-center", targets = "_all"))
        )
      ),
      class = "compact stripe cell-border nowrap"
    )
    
    DT::formatRound(
      dt,
      columns = intersect(c("GC content (%)"), names(df)),
      digits = 0
    )
  })
  
  observeEvent(input$single_aso_results_rows_selected, {
    row <- input$single_aso_results_rows_selected
    
    if (is.null(row) || length(row) == 0) {
      return()
    }
    
    raw_df <- single_aso_results_lookup()
    display_df <- single_aso_results_sorted()
    
    if (is.null(raw_df) || nrow(raw_df) == 0) {
      return()
    }
    
    display_row <- display_df[row, , drop = FALSE]
    
    if (!"Input order" %in% names(display_row)) {
      return()
    }
    
    row_data <- raw_df
    
    if ("Input order" %in% names(display_row) &&
        "input_order" %in% names(row_data)) {
      row_data <- row_data %>%
        dplyr::filter(input_order == display_row$`Input order`[[1]])
    }
    
    if ("Start position in gene" %in% names(display_row) &&
        "start" %in% names(row_data)) {
      selected_start <- display_row$`Start position in gene`[[1]]
      
      if (!is.na(selected_start)) {
        row_data <- row_data %>%
          dplyr::filter(start == selected_start)
      }
    }
    
    if ("Chromosome start pos." %in% names(display_row) &&
        "chr_start" %in% names(row_data)) {
      selected_chr_start <- display_row$`Chromosome start pos.`[[1]]
      
      if (!is.na(selected_chr_start)) {
        row_data <- row_data %>%
          dplyr::filter(chr_start == selected_chr_start)
      }
    }
    
    if ("Consensus source" %in% names(display_row) &&
        "consensus_source" %in% names(row_data)) {
      selected_consensus <- display_row$`Consensus source`[[1]]
      
      if (!is.na(selected_consensus)) {
        row_data <- row_data %>%
          dplyr::filter(consensus_source == selected_consensus)
      }
    }
    
    row_data <- row_data %>%
      dplyr::slice(1)
    
    if (nrow(row_data) == 0) {
      return()
    }
    
    if (nrow(row_data) == 0) {
      return()
    }
    
    seq <- toupper(trimws(row_data$name[[1]]))
    
    selected_target_single_aso(row_data)
    oligo_sequence_single_aso(row_data$oligo_seq[[1]])
    
    shinyjs::enable("single_aso_download_rnaseh")
    
    rnaseh_data <- rnaseh_results(
      selected_row_name = row_data$name[[1]],
      oligo_seq = row_data$oligo_seq[[1]],
      mod_5prime = input$single_aso_mod_5prime,
      mod_3prime = input$single_aso_mod_3prime
    )
    
    rnaseh_stored_single_aso(rnaseh_data)
    
    output$single_aso_rnaseh_title <- renderText({
      paste0("RNase H results for: ", row_data$name[[1]])
    })
    
    output$single_aso_rnaseh_info <- renderUI({
      HTML(
        paste0(
          "Length of sequence: ",
          row_data$length[[1]],
          "<br>",
          "Oligo sequence (ASO): ",
          row_data$oligo_seq[[1]]
        )
      )
    })
    
    output$single_aso_rnaseh_results <- DT::renderDT({
      rnaseh_view <- rename_rnaseh_cols(rnaseh_data)
      
      DT::datatable(
        rnaseh_view,
        selection = list(mode = "single", selected = 1),
        options = list(
          paging = FALSE,
          searching = FALSE,
          info = FALSE,
          dom = "t"
        )
      )
    })
    
    output$single_aso_cleavage_visual <- renderUI(div())
    
    current_seq_single_aso(seq)
    current_mismatch_single_aso(2)
    
    updateSelectInput(
      session,
      "single_aso_user_mismatch",
      selected = 2
    )
    
    output$single_aso_offtarget_title <- renderText({
      paste0("Off target results for: ", seq)
    })
    
    output$single_aso_seq <- renderText({
      paste0("ASO sequence: ", row_data$oligo_seq[[1]])
    })
    
    ot_total <- single_aso_offtargets_total()
    
    if (!is.null(ot_total) && nrow(ot_total) > 0) {
      default_subset <- ot_total %>%
        dplyr::mutate(name_lookup = toupper(trimws(as.character(name)))) %>%
        dplyr::filter(name_lookup == seq, distance <= 2) %>%
        dplyr::select(-name_lookup)
    } else {
      default_subset <- tibble()
    }
    
    current_offtargets_single_aso(default_subset)
    
    if (!is.null(default_subset) && nrow(default_subset) > 0) {
      shinyjs::enable("single_aso_download_offtarget")
    } else {
      shinyjs::disable("single_aso_download_offtarget")
    }
    
    output$single_aso_numb_offtargets <- renderText({
      paste0("# off targets: ", nrow(default_subset))
    })
    
    updateTabsetPanel(
      session,
      "tabs_main_single_aso",
      selected = "Off target results"
    )
  })
  
  output$single_aso_offtarget_results <- DT::renderDT({
    df <- current_offtargets_single_aso()
    
    if (is.null(df) || nrow(df) == 0) {
      return(
        DT::datatable(
          data.frame(Message = "No off-targets available for the selected Specific ASO."),
          rownames = FALSE,
          options = list(dom = "t"),
          class = "compact stripe"
        )
      )
    }
    
    df <- df %>%
      dplyr::mutate(
        match_string_display = chartr("ID", "DI", match_string),
        subject_seq_display = purrr::map2_chr(
          subject_seq,
          match_string_display,
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
          "pLI",
          "oe_lof_upper"
        )),
        dplyr::everything()
      )
    
    display_map <- c(
      gene_name            = "Gene symbol",
      transcript           = "Transcript (Ensembl)",
      match_string_display = "Match (alignment)",
      subject_seq_display  = "Off-target sequence",
      off_target_nt_length = "Length (nt)",
      distance             = "Number of Mismatches/Indels",
      mismatches           = "Mismatches",
      deletions            = "Insertions",
      insertions           = "Deletions",
      pLI                  = "GnomAD pLI",
      oe_lof_upper         = "GnomAD LOEUF"
    )
    
    nm <- names(df_view)
    nm2 <- ifelse(nm %in% names(display_map), display_map[nm], nm)
    names(df_view) <- make.unique(nm2, sep = " (dup) ")
    
    DT::datatable(
      df_view,
      rownames = FALSE,
      escape = FALSE,
      options = list(scrollX = TRUE),
      class = "compact stripe nowrap"
    )
  })
  
  observeEvent(input$single_aso_apply_mismatch, {
    req(current_seq_single_aso())
    
    ot_total <- single_aso_offtargets_total()
    
    mm <- as.numeric(input$single_aso_user_mismatch)
    seq <- toupper(current_seq_single_aso())
    
    current_mismatch_single_aso(mm)
    
    if (mm %in% c(0, 1, 2)) {
      
      if (is.null(ot_total) || nrow(ot_total) == 0) {
        subset_df <- tibble()
      } else {
        subset_df <- ot_total %>%
          dplyr::mutate(name_lookup = toupper(trimws(as.character(name)))) %>%
          dplyr::filter(name_lookup == seq, distance <= mm) %>%
          dplyr::select(-name_lookup)
      }
      
    } else if (mm == 3) {
      
      subset_df <- tryCatch({
        new_res <- all_offt(seq, mismatches_allowed = 3)
        
        if (is.null(new_res) || !is.data.frame(new_res) || nrow(new_res) == 0) {
          tibble()
        } else {
          new_res$name <- seq
          new_res$length <- nchar(seq)
          
          new_res <- new_res %>%
            dplyr::mutate(
              distance = mismatches + deletions + insertions,
              gene_name = stringr::str_extract(line, "(?<=\\|)[^;]+")
            ) %>%
            dplyr::distinct(gene_name, match_string, query_seq, .keep_all = TRUE)
          
          tmp <- tempfile(fileext = ".bgz")
          
          download.file(
            "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz",
            destfile = tmp,
            mode = "wb"
          )
          
          gnomad_df <- readr::read_tsv(
            tmp,
            show_col_types = FALSE,
            col_select = c(gene, transcript, pLI, oe_lof_upper)
          )
          
          unlink(tmp)
          
          new_res %>%
            dplyr::left_join(
              gnomad_df,
              by = c("gene_name" = "gene")
            )
        }
      }, error = function(e) {
        showNotification(
          paste("3-mismatch GGGenome lookup failed:", conditionMessage(e)),
          type = "error",
          duration = 8
        )
        tibble()
      })
      
    } else {
      subset_df <- tibble()
    }
    
    current_offtargets_single_aso(subset_df)
    
    if (!is.null(subset_df) && nrow(subset_df) > 0) {
      shinyjs::enable("single_aso_download_offtarget")
    } else {
      shinyjs::disable("single_aso_download_offtarget")
    }
    
    output$single_aso_numb_offtargets <- renderText({
      paste0("# off targets: ", nrow(subset_df))
    })
  })
  
  output$single_aso_download_offtarget <- downloadHandler(
    filename = function() {
      paste0(
        "specific_aso_offtargets_",
        current_seq_single_aso(),
        "_mismatches_",
        current_mismatch_single_aso(),
        "_",
        Sys.Date(),
        ".csv"
      )
    },
    content = function(con) {
      data_offtarget <- current_offtargets_single_aso()
      req(data_offtarget)
      write.csv2(data_offtarget, con, row.names = FALSE)
    }
  )
  
  output$Download_single_aso_results <- downloadHandler(
    filename = function() {
      mode_label <- dplyr::case_when(
        identical(input$single_aso_analysis_mode, "no_reference") ~ "no_reference",
        identical(input$single_aso_analysis_mode, "with_reference_full") ~ "reference_full",
        identical(input$single_aso_analysis_mode, "with_reference_consensus_cf") ~ "consensus_cf",
        TRUE ~ "specific_aso"
      )
      
      paste0(
        "specific_aso_results_",
        mode_label,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".csv"
      )
    },
    content = function(file) {
      write.csv(
        single_aso_results_export_data(),
        file,
        row.names = FALSE
      )
    }
  )
  
  observeEvent(input$single_aso_add_mods, {
    row_data <- selected_target_single_aso()
    
    if (is.null(row_data)) {
      showNotification(
        "Select a Specific ASO row first.",
        type = "warning",
        duration = 6
      )
      return()
    }
    
    rnaseh_data <- rnaseh_results(
      selected_row_name = row_data$name[[1]],
      oligo_seq = row_data$oligo_seq[[1]],
      mod_5prime = input$single_aso_mod_5prime,
      mod_3prime = input$single_aso_mod_3prime
    )
    
    rnaseh_stored_single_aso(rnaseh_data)
    
    output$single_aso_rnaseh_results <- DT::renderDT({
      rnaseh_view <- rename_rnaseh_cols(rnaseh_data)
      
      DT::datatable(
        rnaseh_view,
        selection = list(mode = "single", selected = 1),
        options = list(
          paging = FALSE,
          searching = FALSE,
          info = FALSE,
          dom = "t"
        )
      )
    })
    
    output$single_aso_cleavage_visual <- renderUI(div())
  })
  
  observeEvent(input$single_aso_rnaseh_results_rows_selected, {
    row_number <- input$single_aso_rnaseh_results_rows_selected
    
    if (length(row_number) == 0) {
      return()
    }
    
    row_data <- selected_target_single_aso()
    
    if (is.null(row_data)) {
      return()
    }
    
    rnaseh_data <- rnaseh_stored_single_aso()
    
    if (is.null(rnaseh_data)) {
      return()
    }
    
    selected_row <- rnaseh_data[row_number, ]
    
    reverse_string <- function(x) {
      paste0(rev(strsplit(x, "")[[1]]), collapse = "")
    }
    
    oligo_seq <- row_data$oligo_seq[[1]]
    
    mod5 <- input$single_aso_mod_5prime
    mod3 <- input$single_aso_mod_3prime
    
    oligo_len <- nchar(oligo_seq)
    
    mod5_region <- substr(oligo_seq, 1, mod5)
    mod3_region <- substr(oligo_seq, oligo_len - mod3 + 1, oligo_len)
    mid_region <- substr(oligo_seq, mod5 + 1, oligo_len - mod3)
    
    position_string <- stringr::str_split_1(selected_row$position, " - ")
    start_pos <- as.numeric(position_string[1])
    
    target_seq_fw <- row_data$name[[1]]
    target_seq_rv <- reverse_string(target_seq_fw)
    
    rna_len <- nchar(target_seq_fw)
    
    cleavage_start_fw <- start_pos
    cleavage_pos_fw <- cleavage_start_fw + 6
    cleavage_pos_rv <- rna_len - cleavage_pos_fw
    cleavage_site_up <- cleavage_pos_rv - 1
    cleavage_site_down <- cleavage_pos_rv + 7
    
    upstream <- substr(target_seq_rv, 1, cleavage_site_up - 1)
    site_start <- substr(target_seq_rv, cleavage_site_up, cleavage_pos_rv - 1)
    cut_site <- substr(target_seq_rv, cleavage_pos_rv, cleavage_pos_rv)
    site_down <- substr(target_seq_rv, cleavage_pos_rv + 1, cleavage_site_down)
    downstream <- substr(target_seq_rv, cleavage_site_down + 1, rna_len)
    
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
    
    output$single_aso_cleavage_visual <- renderUI({
      HTML(
        paste0(
          "<h5 style='margin-top:0;'>Oligo sequence (ASO): </h5>",
          "<div style='font-family: monospace; white-space: pre; font-size: 20px; line-height: 1.7;'>",
          oligo_visual_fw,
          "</div>",
          "<h5 style='margin-bottom:6px;'>RNA sequence: </h5>",
          "<div style='font-family: monospace; white-space: pre; font-size: 20px; line-height: 1.7;'>",
          rna_visual_rv,
          "</div>"
        )
      )
    })
  })
  
  output$single_aso_download_rnaseh <- downloadHandler(
    filename = function() {
      row_data <- selected_target_single_aso()
      
      if (is.null(row_data)) {
        "specific_ASO_RNaseH_results.xlsx"
      } else {
        paste0("specific_ASO_RNaseH_results_", row_data$name[[1]], ".xlsx")
      }
    },
    content = function(file) {
      data <- rnaseh_stored_single_aso()
      
      if (is.null(data)) {
        data <- data.frame(Message = "No available data")
      }
      
      data[] <- lapply(data, as.character)
      
      wb <- createWorkbook()
      addWorksheet(wb, "RNaseH_results")
      writeData(wb, "RNaseH_results", data)
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
#@@@#
  
  
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



