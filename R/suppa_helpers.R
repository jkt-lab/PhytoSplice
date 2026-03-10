# R/suppa_helpers.R

#' Load SUPPA Results
#'
#' Reads .dpsi and .psivec files.
#'
#' @param dir_path Path to the directory containing output files
#' @return A list with $dpsi and $psivec dataframes, or NULL if files missing
load_suppa_data <- function(dir_path) {
  # Try "diff_results" prefix first (used by current module)
  dpsi_file <- file.path(dir_path, "diff_results.dpsi")
  psivec_file <- file.path(dir_path, "diff_results.psivec")
  
  # Fallback to "suppa" prefix
  if (!file.exists(dpsi_file)) dpsi_file <- file.path(dir_path, "suppa.dpsi")
  if (!file.exists(psivec_file)) psivec_file <- file.path(dir_path, "suppa.psivec")
  
  if (!file.exists(dpsi_file) || !file.exists(psivec_file)) {
    return(NULL)
  }
  
  tryCatch({
    # SUPPA output often has N columns of data but N-1 headers.
    # read.table automatically handles this by setting the first column as row.names.
    
    # Load dPSI
    dpsi <- read.table(dpsi_file, header=TRUE, stringsAsFactors=FALSE, sep="\t", check.names=FALSE, na.strings = c("NA", "nan"))
    if (!"Event_id" %in% names(dpsi)) {
       # If Event_id is not a column, it's likely in the rownames
       dpsi$Event_id <- rownames(dpsi)
       rownames(dpsi) <- NULL
    }
    
    # Load PSIvec
    psivec <- read.table(psivec_file, header=TRUE, stringsAsFactors=FALSE, sep="\t", check.names=FALSE, na.strings = c("NA", "nan"))
    if (!"Event_id" %in% names(psivec)) {
       psivec$Event_id <- rownames(psivec)
       rownames(psivec) <- NULL
    }
    
    # Standardize column names for dPSI (The header usually contains "control-case_dPSI" etc.)
    # We want generic names "mean_dPSI" and "p_val" for downstream plotting
    # Look for column ending in _dPSI and _p-val, or exact matches
    dpsi_col <- grep("(_dPSI$|^mean_dPSI$)", names(dpsi), value=TRUE)
    pval_col <- grep("(_p-val$|^p_val$)", names(dpsi), value=TRUE)
    
    if (length(dpsi_col) > 0) names(dpsi)[names(dpsi) == dpsi_col[1]] <- "mean_dPSI"
    if (length(pval_col) > 0) names(dpsi)[names(dpsi) == pval_col[1]] <- "p_val"
    
    # Validate final structure
    if (!"Event_id" %in% names(dpsi) || !"mean_dPSI" %in% names(dpsi) || !"p_val" %in% names(dpsi)) {
       warning("SUPPA output columns could not be parsed. Expected Event_id, *_dPSI, *_p-val.")
       return(NULL)
    }
    
    return(list(dpsi = dpsi, psivec = psivec))
  }, error = function(e) {
    print(paste("Error loading SUPPA data:", e$message))
    return(NULL)
  })
}

#' Filter SUPPA Results
#'
#' @param df Data frame from load_suppa_data
#' @param pval_cut P-value cutoff
#' @param dpsi_cut Delta PSI cutoff
#' @return Filtered data frame
filter_suppa_data <- function(df, pval_cut = 0.05, dpsi_cut = 0.1) {
  if (is.null(df)) return(NULL)
  if (!all(c("p_val", "mean_dPSI") %in% names(df))) return(df)
  
  df[df$p_val < pval_cut & abs(df$mean_dPSI) > dpsi_cut, ]
}

#' Merge Salmon Quantifications (TPM)
#' 
#' Parses multiple quant.sf files and merges TPM columns into a single matrix for SUPPA.
#'
#' @param quant_files A named character vector where names are sample IDs and values are paths to quant.sf files.
#' @param output_file Optional path to write the merged matrix (tab-separated).
#'
#' @return A dataframe of merged TPMs.
merge_salmon_quants <- function(quant_files, output_file = NULL) {
  if (length(quant_files) == 0) return(NULL)
  
  merged_df <- NULL
  
  for (sample_id in names(quant_files)) {
    fpath <- quant_files[[sample_id]]
    if (!file.exists(fpath)) {
      warning(paste("File not found:", fpath))
      next
    }
    
    # Sanitize sample_id (replace spaces with underscores) to prevent SUPPA header issues
    clean_sample_id <- gsub("\\s+", "_", sample_id)
    
    # Read quant.sf (Name, Length, EffectiveLength, TPM, NumReads)
    # Force sep="\t" to handle transcript names with spaces correctly
    df <- tryCatch({
      read.table(fpath, header = TRUE, stringsAsFactors = FALSE, sep = "\t", quote = "")
    }, error = function(e) NULL)
    
    if (is.null(df)) next
    
    # Extract ID and TPM
    sub_df <- df[, c("Name", "TPM")]
    colnames(sub_df) <- c("target_id", clean_sample_id)
    
    # Remove empty IDs or NAs
    sub_df <- sub_df[!is.na(sub_df$target_id) & sub_df$target_id != "", ]
    # Sanitize target_id (no spaces)
    sub_df$target_id <- gsub("\\s+", "_", sub_df$target_id)
    
    if (is.null(merged_df)) {
      merged_df <- sub_df
    } else {
      merged_df <- merge(merged_df, sub_df, by = "target_id", all = TRUE)
    }
  }
  
  if (!is.null(merged_df)) {
    # Remove any empty/NA target_id that might have slipped in (e.g. from merge)
    merged_df <- merged_df[!is.na(merged_df$target_id) & merged_df$target_id != "", ]
    
    # Replace NAs with 0 (if transcripts missing in some samples)
    merged_df[is.na(merged_df)] <- 0
    
    if (!is.null(output_file)) {
      write.table(merged_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
  
  return(merged_df)
}
