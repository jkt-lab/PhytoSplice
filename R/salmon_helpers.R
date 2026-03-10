# R/salmon_helpers.R

#' Build Salmon Index (Decoy-Aware)
#'
#' @param genome_fa Path to Genome FASTA
#' @param gtf_path Path to Annotation GTF
#' @param output_dir Directory where the index and intermediate files will be stored
#' @param threads Number of threads to use
#' @param callbacks List of callback functions (log, update_progress)
#'
#' @return Boolean indicating success
build_salmon_index <- function(genome_fa, gtf_path, output_dir, threads = 4, callbacks = list()) {
  
  # Helper for logging
  log_msg <- function(msg) {
    if (!is.null(callbacks$log)) callbacks$log(paste0("[ ", format(Sys.time(), "%H:%M:%S"), "] ", msg))
    cat(msg, "\n")
  }
  
  # Validate inputs
  if (!file.exists(genome_fa)) stop("Genome FASTA not found.")
  if (!file.exists(gtf_path)) stop("GTF file not found.")
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Resolve tool paths (prioritize local env)
  env_bin <- file.path(getwd(), "env", "bin")
  gffread_exe <- if (file.exists(file.path(env_bin, "gffread"))) file.path(env_bin, "gffread") else "gffread"
  salmon_exe <- if (file.exists(file.path(env_bin, "salmon"))) file.path(env_bin, "salmon") else "salmon"
  
  # Define intermediate file paths
  transcripts_fa <- file.path(output_dir, "transcripts.fa")
  decoys_txt <- file.path(output_dir, "decoys.txt")
  gentrome_fa <- file.path(output_dir, "gentrome.fa")
  index_dir <- file.path(output_dir, "salmon_index")
  
  # --- Step 2: Extract Transcriptome FASTA ---
  log_msg("Step 2: Extracting transcriptome FASTA using gffread...")
  cmd_gffread <- sprintf("'%s' '%s' -g '%s' -w '%s'", gffread_exe, gtf_path, genome_fa, transcripts_fa)
  log_msg(paste("CMD:", cmd_gffread))
  
  ret <- system(cmd_gffread)
  if (ret != 0) {
    log_msg("Error: gffread failed.")
    return(FALSE)
  }
  
  # --- Step 3: Create Decoy List ---
  log_msg("Step 3: Creating decoy list from genome...")
  # Extract only the first word (ID) from the header
  cmd_decoys <- sprintf("grep '^>' '%s' | awk '{print $1}' | sed 's/^>//' > '%s'", genome_fa, decoys_txt)
  log_msg(paste("CMD:", cmd_decoys))
  
  ret <- system(cmd_decoys)
  if (ret != 0) {
    log_msg("Error: Failed to create decoys list.")
    return(FALSE)
  }
  
  # --- Step 4: Combine Transcriptome + Genome (Gentrome) ---
  log_msg("Step 4: Combining transcriptome and genome (Gentrome)...")
  # cat transcripts.fa genome.fa > gentrome.fa
  cmd_cat <- sprintf("cat '%s' '%s' > '%s'", transcripts_fa, genome_fa, gentrome_fa)
  log_msg(paste("CMD:", cmd_cat))
  
  ret <- system(cmd_cat)
  if (ret != 0) {
    log_msg("Error: Failed to create gentrome FASTA.")
    return(FALSE)
  }
  
  # --- Step 5: Build Salmon Index ---
  log_msg("Step 5: Building Salmon Index (This may take a while)...")
  # salmon index -t gentrome.fa -d decoys.txt -i salmon_index --threads 8
  cmd_index <- sprintf("'%s' index -t '%s' -d '%s' -i '%s' --threads %d", 
                       salmon_exe, gentrome_fa, decoys_txt, index_dir, threads)
  log_msg(paste("CMD:", cmd_index))
  
  # Run Salmon (capture output to log if possible, or just system)
  ret <- system(cmd_index)
  
  if (ret == 0) {
    log_msg("✅ Salmon Index successfully created!")
    log_msg(paste("Location:", index_dir))
    return(TRUE)
  } else {
    log_msg("❌ Error: Salmon indexing failed.")
    return(FALSE)
  }
}

#' Merge Salmon Quantifications into a Single TPM Matrix
#'
#' @param quant_files Named vector of paths to quant.sf files. Names will be used as sample IDs.
#' @param output_file Path to write the merged TPM matrix.
#'
#' @return The merged data frame or NULL on error.
merge_salmon_quants <- function(quant_files, output_file) {
  if (length(quant_files) == 0) return(NULL)
  
  message("Merging ", length(quant_files), " quantification files...")
  
  # Use R-based merge for robustness against header issues and file differences
  tryCatch({
    # Initialize with the first file
    first_file <- quant_files[1]
    if (!file.exists(first_file)) {
      message("Error: First file not found: ", first_file)
      return(NULL)
    }
    
    # Read first file
    # Ensure sep="\t" and handle potential quoting issues
    base_df <- read.table(first_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
    
    if (!"Name" %in% names(base_df) || !"TPM" %in% names(base_df)) {
      message("Error: quant.sf format invalid (Missing Name or TPM columns).")
      return(NULL)
    }
    
    # Initialize merged dataframe with Target ID
    # Sanitize Target ID (remove spaces)
    base_df$Name <- gsub("\\s+", "_", base_df$Name)
    merged_df <- data.frame(target_id = base_df$Name, stringsAsFactors = FALSE)
    
    # Add first sample (Sanitize Sample Name)
    sname <- gsub("\\s+", "_", names(quant_files)[1])
    merged_df[[sname]] <- base_df$TPM
    
    # Loop through remaining files
    if (length(quant_files) > 1) {
      for (i in 2:length(quant_files)) {
        f <- quant_files[i]
        sname_raw <- names(quant_files)[i]
        sname <- gsub("\\s+", "_", sname_raw) # Sanitize
        
        if (file.exists(f)) {
          df <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
          
          # Sanitize ID
          df$Name <- gsub("\\s+", "_", df$Name)
          
          # Optimization: If IDs match exactly, just cbind (much faster)
          if (nrow(df) == nrow(merged_df) && all(df$Name == merged_df$target_id)) {
            merged_df[[sname]] <- df$TPM
          } else {
            # Robust merge
            tmp <- data.frame(target_id = df$Name, TPM = df$TPM, stringsAsFactors = FALSE)
            colnames(tmp)[2] <- sname
            merged_df <- merge(merged_df, tmp, by = "target_id", all = TRUE)
          }
        } else {
          message("Warning: File missing: ", f)
        }
      }
    }
    
    # Replace NAs with 0
    merged_df[is.na(merged_df)] <- 0
    
    # Write output
    # quote=FALSE is important for tools that don't expect quotes
    # sep="\t"
    write.table(merged_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    
    message("Merge successful.")
    return(merged_df)
    
  }, error = function(e) {
    message("Merge failed: ", e$message)
    return(NULL)
  })
}