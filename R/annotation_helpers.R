# Helper functions for annotating rMATS results using GTF or MyGene.info API

annotate_rmats_output <- function(output_dir, gtf_path) {
  if (!dir.exists(output_dir)) return(FALSE)
  
  # Initialize lookup map
  symbol_lookup <- NULL
  
  # --- Strategy 1: GTF Parsing (Local) ---
  if (file.exists(gtf_path)) {
    message("Attempting annotation using GTF: ", gtf_path)
    tryCatch({
      map_file <- file.path(output_dir, "id_symbol_map.tmp")
      
      # 1. Try gene_name (Standard)
      # Use awk for robust attribute parsing if simple grep fails, but strict grep is faster.
      # Let's stick to the pipeline that worked for gene_name before, but add robustness.
      
      cmd_gn <- paste0(
        "grep -v '^#' '", gtf_path, "' | ",
        "grep -o 'gene_id \"[^\"]*\".*gene_name \"[^\"]*\"' | ",
        "sed -e 's/gene_id \"//' -e 's/\";.*gene_name \"/\t/' -e 's/\".*//' | ",
        "sort -u > '", map_file, "'"
      )
      system(cmd_gn)
      
      # 2. Try gene_symbol
      if (file.size(map_file) == 0) {
         cmd_gs <- paste0(
          "grep -v '^#' '", gtf_path, "' | ",
          "grep -o 'gene_id \"[^\"]*\".*gene_symbol \"[^\"]*\"' | ",
          "sed -e 's/gene_id \"//' -e 's/\";.*gene_symbol \"/\t/' -e 's/\".*//' | ",
          "sort -u > '", map_file, "'"
        )
        system(cmd_gs)
      }
      
      # 3. Try description (NCBI Gnomon) - using AWK for feature type safety
      if (file.size(map_file) == 0) {
         # awk to find rows where 3rd column is 'gene', then extract attributes
         # We assume standard GTF formatting (tab delimited). 
         # We extract gene_id and description.
         cmd_desc <- paste0(
          "awk -F '\\t' '$3 == \"gene\" { print $9 }' '", gtf_path, "' | ",
          "grep -o 'gene_id \"[^\"]*\".*description \"[^\"]*\"' | ",
          "sed -e 's/gene_id \"//' -e 's/\";.*description \"/\t/' -e 's/\".*//' | ",
          "sort -u > '", map_file, "'"
        )
        system(cmd_desc)
      }
      
      # 4. Try product (fallback)
      if (file.size(map_file) == 0) {
         cmd_prod <- paste0(
          "grep -v '^#' '", gtf_path, "' | ",
          "grep -o 'gene_id \"[^\"]*\".*product \"[^\"]*\"' | ",
          "sed -e 's/gene_id \"//' -e 's/\";.*product \"/\t/' -e 's/\".*//' | ",
          "sort -u > '", map_file, "'"
        )
        system(cmd_prod)
      }
      
      if (file.size(map_file) > 0) {
        gene_map <- read.table(map_file, sep="\t", header=FALSE, col.names=c("GeneID", "Symbol"), stringsAsFactors=FALSE, quote="")
        symbol_lookup <- setNames(gene_map$Symbol, gene_map$GeneID)
        unlink(map_file)
        message("GTF mapping successful. Found ", length(symbol_lookup), " symbols.")
      } else {
        message("GTF does not contain standard name/symbol/description attributes.")
      }
    }, error = function(e) message("GTF parsing error: ", e$message))
  }
  
  # --- Strategy 2: MyGene.info API (Universal Fallback) ---
  if (is.null(symbol_lookup)) {
    message("Falling back to MyGene.info API...")
    
    # 1. Collect all unique GeneIDs from files
    files <- list.files(output_dir, pattern = "\\.MATS\\.(JC|JCEC)\\.txt$", full.names = TRUE)
    all_ids <- c()
    
    for (f in files) {
      try({
        df <- read.table(f, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE, check.names=FALSE)
        if ("GeneID" %in% names(df)) {
          ids <- gsub('"', '', as.character(df$GeneID))
          all_ids <- unique(c(all_ids, ids))
        }
      }, silent=TRUE)
    }
    
    all_ids <- unique(all_ids)
    
    if (length(all_ids) > 0) {
      message("Fetching annotation for ", length(all_ids), " unique IDs...")
      symbol_lookup <- fetch_mygene_symbols(all_ids)
    }
  }
  
  if (is.null(symbol_lookup) || length(symbol_lookup) == 0) {
    message("Could not resolve any gene symbols via GTF or Online API.")
    return(FALSE)
  }
  
  # --- Apply Annotation ---
  files <- list.files(output_dir, pattern = "\\.MATS\\.(JC|JCEC)\\.txt$", full.names = TRUE)
  
  for (f in files) {
    message("Annotating: ", basename(f))
    tryCatch({
      df <- read.table(f, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE, check.names=FALSE)
      
      if ("GeneID" %in% names(df)) {
        clean_ids <- gsub('"', '', as.character(df$GeneID))
        
        # Map Symbols
        mapped_syms <- symbol_lookup[clean_ids]
        
        # Keep original ID if no match or empty
        final_syms <- ifelse(is.na(mapped_syms) | mapped_syms == "", clean_ids, mapped_syms)
        
        # Update/Create Column
        if ("geneSymbol" %in% names(df)) {
          df$geneSymbol <- final_syms
        } else {
          gid_idx <- which(names(df) == "GeneID")
          df <- data.frame(df[,1:gid_idx, drop=FALSE], geneSymbol=final_syms, df[,(gid_idx+1):ncol(df), drop=FALSE], check.names=FALSE)
        }
        
        write.table(df, f, sep="\t", quote=FALSE, row.names=FALSE)
      }
    }, error = function(e) message("Error updating file ", f, ": ", e$message))
  }
  
  message("Annotation complete.")
  return(TRUE)
}

# Batch fetch using MyGene.info v3 API
fetch_mygene_symbols <- function(id_list) {
  # Clean input list: remove NA, empty strings, and ensure they are characters
  id_list <- unique(as.character(id_list))
  id_list <- id_list[nzchar(id_list) & !is.na(id_list)]
  
  if (length(id_list) == 0) {
    message("No valid IDs provided for annotation.")
    return(NULL)
  }
  
  lookup_map <- list()
  
  # Internal recursive function to handle batches and sub-batches on error
      query_batch <- function(chunk, batch_id, level = 0) {
      server <- "https://mygene.info/v3/query"
      
      body <- list(
        q = chunk,
        scopes = "entrezgene,ensembl.gene,symbol,refseq,accession,alias,uniprot,locus_tag",
        fields = "symbol,name", 
        species = "all"
      )    
    res_map <- list()
    
    tryCatch({
      json_body <- jsonlite::toJSON(body, auto_unbox = TRUE)
      
      r <- httr::POST(server, 
                body = json_body, 
                httr::content_type_json(),
                httr::accept_json())
      
      status <- httr::status_code(r)
      
      if (status == 200) {
        res <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"))
        
        if (!is.null(res) && (is.data.frame(res) && nrow(res) > 0)) {
          if (!"symbol" %in% names(res)) res$symbol <- NA
          if (!"name" %in% names(res)) res$name <- NA
          if (!"query" %in% names(res)) return(res_map)
          
          res_u <- res[!duplicated(res$query), ]
          labels <- ifelse(!is.na(res_u$symbol), res_u$symbol, res_u$query)
          
          is_loc <- grepl("^LOC", labels)
          has_name <- !is.na(res_u$name) & res_u$name != ""
          labels[is_loc & has_name] <- res_u$name[is_loc & has_name]
          
          valid_idx <- which(!is.na(labels))
          if (length(valid_idx) > 0) {
             res_map <- setNames(as.character(labels[valid_idx]), as.character(res_u$query[valid_idx]))
          }
        }
      } else if (status == 400 && length(chunk) > 1 && level < 5) {
        # If 400 error and we have more than one ID, split and retry to isolate problematic ID
        message(paste0("  Batch ", batch_id, " (size ", length(chunk), ") failed with 400. Splitting..."))
        mid <- ceiling(length(chunk) / 2)
        c1 <- chunk[1:mid]
        c2 <- chunk[(mid+1):length(chunk)]
        
        m1 <- query_batch(c1, paste0(batch_id, ".1"), level + 1)
        m2 <- query_batch(c2, paste0(batch_id, ".2"), level + 1)
        
        res_map <- c(m1, m2)
      } else {
        err_msg <- httr::content(r, "text", encoding = "UTF-8")
        message(paste0("MyGene.info API Error (Batch ", batch_id, "): HTTP ", status, " - ", err_msg))
        if (length(chunk) == 1) {
           message(paste0("  Skipping problematic ID: ", chunk))
        }
      }
    }, error = function(e) {
      message("MyGene.info API Exception (Batch ", batch_id, "): ", e$message)
    })
    
    return(res_map)
  }

  # Batch limit is 1000
  chunks <- split(id_list, ceiling(seq_along(id_list)/1000))
  
  withProgress(message = "Querying MyGene.info...", value = 0, {
    n <- length(chunks)
    for (i in seq_along(chunks)) {
      incProgress(1/n, detail = paste("Batch", i, "of", n))
      
      chunk <- chunks[[i]]
      if(length(chunk) == 0) next
      
      batch_results <- query_batch(chunk, as.character(i))
      if (length(batch_results) > 0) {
        # Efficiently merge named lists/vectors
        for (name in names(batch_results)) {
          lookup_map[[name]] <- batch_results[[name]]
        }
      }
      
      Sys.sleep(0.1)
    }
  })
  
  # Convert list to named vector
  if (length(lookup_map) == 0) return(NULL)
  return(unlist(lookup_map))
}