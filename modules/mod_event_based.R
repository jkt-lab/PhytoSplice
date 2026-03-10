mod_event_based_ui <- function(id) {
  ns <- NS(id)
  tagList(
    shinyjs::useShinyjs(),
    uiOutput(ns("main_ui"))
  )
}

mod_event_based_server <- function(id, manager) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    EVENT_TYPES <- c("SE", "MXE", "A3SS", "A5SS", "RI")
    EVENT_NAMES <- c("SE"="Skipped Exon", "MXE"="Mutually Exclusive", "A3SS"="Alt 3' SS", "A5SS"="Alt 5' SS", "RI"="Retained Intron")
    
    rv <- reactiveValues(
      data = list(),
      status = "no_data",
      source_dir = NULL,
      sashimi_image = NULL,
      # Popup Reactives
      popup_event_id = NULL,
      show_more_genes = FALSE
    )
    
    filters <- reactiveValues(fdr = 0.05, dpsi = 0.1, gene_list = "")
    
    observeEvent(input$apply_filters, {
      filters$fdr <- input$fdr_cut
      filters$dpsi <- input$dpsi_cut
      filters$gene_list <- input$gene_filter
      showNotification("Filters applied.", type = "message", duration = 2)
    })
    
    # --- Initialization ---
    observeEvent(manager$events$workspace_loaded, {
      req(manager$events$workspace_loaded > 0)
      req(manager$workspace_dir) 
      
      # 1. Primary Check: Standard Workspace Artifacts Path
      ws_rmats_output <- file.path(manager$workspace_dir, "artifacts", "rmats", "output")
      
      dir_to_load <- NULL
      
      # Check if standard directory exists and contains expected files
      if (dir.exists(ws_rmats_output) && length(list.files(ws_rmats_output, pattern=".MATS.JC.txt")) > 0) {
        dir_to_load <- ws_rmats_output
      } else {
        # 2. Fallback: Check Manager State (if ran in this session)
        res <- manager$get_result("rmats")
        if (!is.null(res) && !is.null(res$output_dir)) {
          candidate <- res$output_dir
          # Resolve relative path if necessary
          if (!startsWith(candidate, "/") && !grepl(":", candidate)) {
             candidate <- file.path(manager$workspace_dir, candidate)
          }
          if (dir.exists(candidate)) dir_to_load <- candidate
        }
      }
      
      if (!is.null(dir_to_load) && dir.exists(dir_to_load)) {
        load_data(dir_to_load) 
      } else {
        rv$status <- "no_file"
      }
      
      saved_params <- manager$get_params(id)
      if (!is.null(saved_params)) {
        updateNumericInput(session, "fdr_cut", value = saved_params$fdr)
        updateNumericInput(session, "dpsi_cut", value = saved_params$dpsi)
        updateTextAreaInput(session, "gene_filter", value = saved_params$gene_list)
        updateSelectInput(session, "selected_type", selected = saved_params$type)
        if (!is.null(saved_params$heatmap_top_n)) {
          updateNumericInput(session, "heatmap_top_n", value = saved_params$heatmap_top_n)
        }
        if (!is.null(saved_params$heatmap_mode)) {
          updateSelectInput(session, "heatmap_mode", selected = saved_params$heatmap_mode)
        }
        if (!is.null(saved_params$heatmap_order)) {
          updateSelectInput(session, "heatmap_order", selected = saved_params$heatmap_order)
        }
        if (!is.null(saved_params$heatmap_cluster_rows)) {
          updateCheckboxInput(session, "heatmap_cluster_rows", value = isTRUE(saved_params$heatmap_cluster_rows))
        }
        if (!is.null(saved_params$heatmap_cluster_cols)) {
          updateCheckboxInput(session, "heatmap_cluster_cols", value = isTRUE(saved_params$heatmap_cluster_cols))
        }
        
        filters$fdr <- saved_params$fdr
        filters$dpsi <- saved_params$dpsi
        filters$gene_list <- saved_params$gene_list
      }
    })
    
    observe({
      req(filters$fdr, filters$dpsi, input$selected_type)
      # We allow gene_list to be empty/null, so don't req() it
      
      # Use isolate() for heatmap inputs to prevent excessive firing if needed, 
      # but here we are in an observe, so we just collect values.
      # Note: 'req' ensures inputs are ready.
      
      manager$state$params[[id]] <- list(
        fdr = filters$fdr,
        dpsi = filters$dpsi,
        gene_list = filters$gene_list,
        type = input$selected_type,
        heatmap_top_n = input$heatmap_top_n,
        heatmap_mode = input$heatmap_mode,
        heatmap_order = input$heatmap_order,
        heatmap_cluster_rows = isTRUE(input$heatmap_cluster_rows),
        heatmap_cluster_cols = isTRUE(input$heatmap_cluster_cols)
      )
    })
    
    load_data <- function(dirpath) {
      rv$data <- list()
      found <- FALSE
      
      # Debug logging
      all_files <- list.files(dirpath)
      cat(paste("Scanning directory:", dirpath, "\nFiles found:", paste(all_files, collapse=", "), "\n"))
      
      for (t in EVENT_TYPES) {
        f <- file.path(dirpath, paste0(t, ".MATS.JC.txt"))
        if (file.exists(f)) {
          tryCatch({
            df <- read.table(f, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE, check.names=TRUE)
            
            # --- Load Novel IDs ---
            novel_pattern <- paste0("^fromGTF\\.novel.*\\.", t, "\\.txt$")
            novel_files <- list.files(dirpath, pattern=novel_pattern, full.names=TRUE)
            novel_ids <- integer(0)
            if (length(novel_files) > 0) {
               for(nf in novel_files) {
                  tryCatch({
                     tmp_df <- read.table(nf, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)
                     if ("ID" %in% names(tmp_df)) novel_ids <- c(novel_ids, tmp_df$ID)
                  }, error = function(e) {})
               }
               novel_ids <- unique(novel_ids)
            }
            df$IsNovel <- df$ID %in% novel_ids
            
            cat(paste("Loaded", t, "- Columns:", paste(names(df)[1:min(5, ncol(df))], collapse=", "), "...\n"))
                        
                        if ("FDR" %in% names(df)) df$FDR <- suppressWarnings(as.numeric(as.character(df$FDR)))
                        if ("IncLevelDifference" %in% names(df)) df$IncLevelDifference <- suppressWarnings(as.numeric(as.character(df$IncLevelDifference)))
                        
                        df$FDR[is.na(df$FDR)] <- 1
            df$IncLevelDifference[is.na(df$IncLevelDifference)] <- 0
            rv$data[[t]] <- df
            found <- TRUE
          }, error = function(e) cat(paste("Error loading", t, ":", e$message, "\n")))
        }
      }
      if (found) {
        rv$status <- "loaded"
        rv$source_dir <- dirpath
        manager$set_result(id, list(source_dir = dirpath))
      }
    }
    
    observeEvent(input$load_manual, { req(input$manual_path); load_data(input$manual_path) })
    
    # --- Annotation Handler ---
    observeEvent(input$annotate_btn, {
      req(rv$source_dir)
      
      # Try to find GTF from rMATS history
      gtf_path <- NULL
      rmats_params <- manager$get_params("rmats")
      if (!is.null(rmats_params) && !is.null(rmats_params$gtf)) {
        gtf_path <- rmats_params$gtf
      }
      
      if (is.null(gtf_path) || !file.exists(gtf_path)) {
        showNotification("GTF file path not found in analysis history. Cannot annotate.", type = "error")
        # Potential future enhancement: Ask user for GTF via modal
      } else {
        withProgress(message = "Annotating Gene Symbols...", {
          success <- annotate_rmats_output(rv$source_dir, gtf_path)
          if (success) {
            showNotification("Annotation complete. Reloading data...", type = "message")
            load_data(rv$source_dir)
          } else {
            showNotification("Annotation failed. Check logs/GTF.", type = "error")
          }
        })
      }
    })
    
    # --- Missing Genes Alert ---
    
    # State for expanding the missing genes list
    rv$show_more_genes <- FALSE
    
    observeEvent(input$toggle_missing_genes, {
      rv$show_more_genes <- !rv$show_more_genes
    })
    
    output$filter_alert <- renderUI({
      req(rv$status == "loaded")
      
      filter_str <- filters$gene_list
      if (is.null(filter_str) || !nzchar(trimws(filter_str))) return(NULL)
      
      clean_str <- gsub("[\"'`]", "", filter_str)
      clean_str <- gsub("[\n\r]", ",", clean_str)
      clean_str <- gsub("\\s+", ",", clean_str)
      targets <- trimws(unlist(strsplit(clean_str, ",")))
      targets <- targets[nzchar(targets)]
      if (length(targets) == 0) return(NULL)
      
      targets_lower <- tolower(targets)
      
      df <- get_dataset()
      if (is.null(df)) return(NULL)
      
      found_sym <- if("geneSymbol" %in% names(df)) tolower(as.character(df$geneSymbol)) else character(0)
      found_gid <- if("GeneID" %in% names(df)) tolower(as.character(df$GeneID)) else character(0)
      
      present_mask <- (targets_lower %in% found_sym) | (targets_lower %in% found_gid)
      missing_genes <- targets[!present_mask]
      
      if (length(missing_genes) > 0) {
        
        # Display Logic
        content <- NULL
        trunc_limit <- 20
        
        if (length(missing_genes) <= trunc_limit) {
           # Simple case
           content <- tags$code(paste(missing_genes, collapse = ", "))
        } else {
           # Expandable case
           if (rv$show_more_genes) {
             # Expanded
             content <- tagList(
               tags$code(paste(missing_genes, collapse = ", ")),
               br(),
               actionLink(ns("toggle_missing_genes"), "Show Less", icon = icon("chevron-up"))
             )
           } else {
             # Truncated
             shown <- head(missing_genes, trunc_limit)
             content <- tagList(
               tags$code(paste(shown, collapse = ", ")),
               span(paste("... (+", length(missing_genes) - trunc_limit, "more)")),
               br(),
               actionLink(ns("toggle_missing_genes"), "Show More", icon = icon("chevron-down"))
             )
           }
        }
        
        div(class = "alert alert-warning alert-dismissible", style="margin-top: 10px;",
            tags$button(type="button", class="close", `data-dismiss`="alert", "x"),
            icon("exclamation-triangle"), 
            strong("Warning: "), 
            "The following genes/IDs were not found in the current Event Type data:",
            br(),
            content
        )
      } else {
        NULL
      }
    })

    # --- Helpers ---
    get_dataset <- reactive({
      req(rv$status == "loaded", input$selected_type)
      df <- rv$data[[input$selected_type]]
      apply_event_filter(df)
    })
    
    apply_event_filter <- function(df) {
      if (!is.null(df) && "IsNovel" %in% names(df)) {
         novel_on <- if (is.null(input$filter_novel)) TRUE else input$filter_novel
         known_on <- if (is.null(input$filter_known)) TRUE else input$filter_known
         
         if (!novel_on) df <- df[!df$IsNovel, ]
         if (!known_on) df <- df[df$IsNovel, ]
      }
      return(df)
    }
    
    apply_gene_filter <- function(df, filter_str) {
      if (is.null(filter_str) || !nzchar(trimws(filter_str))) return(df)
      
      # Clean input: replace quotes, newlines, spaces with commas, then split
      clean_str <- gsub("[\"'`]", "", filter_str)     # Remove quotes
      clean_str <- gsub("[\n\r]", ",", clean_str)     # Replace newlines with commas
      clean_str <- gsub("\\s+", ",", clean_str)       # Replace spaces with commas
      
      # Split by comma
      targets <- trimws(unlist(strsplit(clean_str, ",")))
      targets <- targets[nzchar(targets)]
      if (length(targets) == 0) return(df)
      
      # Normalize targets to lowercase for case-insensitive matching
      targets <- tolower(targets)
      
      # Check columns
      has_sym <- "geneSymbol" %in% names(df)
      has_gid <- "GeneID" %in% names(df)
      
      if (!has_sym && !has_gid) return(df)
      
      # Filter with robust quote stripping for rMATS IDs
      df %>% dplyr::filter({
        match_sym <- if (has_sym) tolower(gsub('"', '', as.character(geneSymbol))) %in% targets else FALSE
        match_gid <- if (has_gid) tolower(gsub('"', '', as.character(GeneID))) %in% targets else FALSE
        match_sym | match_gid
      })
    }

    global_summary <- reactive({
      req(rv$status == "loaded")
      summary_df <- data.frame()
      for (t in names(rv$data)) {
        # Apply Event Filter
        df <- apply_event_filter(rv$data[[t]])
        
        # Global summary usually ignores gene filter to show overall stats, 
        # BUT if user filters, maybe they want stats for that gene set?
        # Let's keep global summary GLOBAL (unfiltered by gene) for context,
        # unless requested otherwise. The prompt says "include feature in heatmap, forest plot, volcano & scatter".
        # It doesn't explicitly mention the global bar chart. I will leave global summary UNFILTERED for now.
        
        sig <- df %>%
          dplyr::filter(FDR < filters$fdr, abs(IncLevelDifference) > filters$dpsi)
        if (nrow(sig) > 0) {
          counts <- sig %>%
            mutate(Direction = ifelse(IncLevelDifference > 0, "Included", "Skipped")) %>%
            count(Direction) %>%
            mutate(Type = t)
          summary_df <- rbind(summary_df, counts)
        }
      }
      summary_df
    })
    
    get_sample_names <- reactive({
      b1_files <- NULL
      b2_files <- NULL
      
      rmats_params <- manager$get_params("rmats")
      if (!is.null(rmats_params)) {
        if (!is.null(rmats_params$b1) && nzchar(rmats_params$b1)) {
          b1_files <- trimws(unlist(strsplit(rmats_params$b1, ",")))
        }
        if (!is.null(rmats_params$b2) && nzchar(rmats_params$b2)) {
          b2_files <- trimws(unlist(strsplit(rmats_params$b2, ",")))
        }
      }
      
      if ((length(b1_files) == 0 || length(b2_files) == 0) && !is.null(manager$workspace_dir)) {
        case_txt <- file.path(manager$workspace_dir, "artifacts", "rmats", "input", "case.txt")
        ctrl_txt <- file.path(manager$workspace_dir, "artifacts", "rmats", "input", "control.txt")
        if (file.exists(case_txt)) {
          b1_raw <- paste(readLines(case_txt, warn = FALSE), collapse = "")
          b1_files <- trimws(unlist(strsplit(b1_raw, ",")))
        }
        if (file.exists(ctrl_txt)) {
          b2_raw <- paste(readLines(ctrl_txt, warn = FALSE), collapse = "")
          b2_files <- trimws(unlist(strsplit(b2_raw, ",")))
        }
      }
      
      clean_names <- function(files) {
        if (length(files) == 0) return(character(0))
        names <- basename(files)
        # Remove file extension (.bam, .fastq, etc)
        names <- gsub("\\.[^.]+$", "", names)
        # Remove common pipeline suffixes
        names <- gsub("(\\.sorted|_sorted|\\.trimmed|_trimmed)$", "", names, ignore.case = TRUE)
        make.unique(names)
      }
      
      case_names <- clean_names(b1_files)
      ctrl_names <- clean_names(b2_files)
      
      list(case = case_names, ctrl = ctrl_names)
    })

    # --- Sashimi Generation Helper ---
    generate_sashimi_internal <- function(target_event_id, is_popup = FALSE) {
      req(target_event_id)
      
      # 1. Get BAM Paths
      b1_files <- NULL
      b2_files <- NULL
      
      rmats_params <- manager$get_params("rmats")
      
      if (!is.null(rmats_params)) {
        b1_files <- trimws(unlist(strsplit(rmats_params$b1, ",")))
        b2_files <- trimws(unlist(strsplit(rmats_params$b2, ",")))
      } else {
        # Fallback: Check Workspace Artifacts
        if (!is.null(manager$workspace_dir)) {
           case_txt <- file.path(manager$workspace_dir, "artifacts", "rmats", "input", "case.txt")
           ctrl_txt <- file.path(manager$workspace_dir, "artifacts", "rmats", "input", "control.txt")
           
           if (file.exists(case_txt) && file.exists(ctrl_txt)) {
             tryCatch({
               b1_raw <- paste(readLines(case_txt, warn=FALSE), collapse="")
               b2_raw <- paste(readLines(ctrl_txt, warn=FALSE), collapse="")
               b1_files <- trimws(unlist(strsplit(b1_raw, ",")))
               b2_files <- trimws(unlist(strsplit(b2_raw, ",")))
             }, error = function(e) warning("Failed to read BAM paths from artifacts: ", e$message))
           }
        }
      }
      
      if (is.null(b1_files) || is.null(b2_files) || length(b1_files)==0) {
        showNotification("BAM files not found. Ensure rMATS was run or workspace contains input/case.txt config.", type = "error")
        return(NULL)
      }
      
      # 2. Get Event Coordinates
      df <- get_dataset()
      
      # Input should be the ID directly (numeric)
      event_id <- suppressWarnings(as.numeric(target_event_id))
      
      if (is.na(event_id)) {
        showNotification(paste("Invalid Event ID selected:", target_event_id), type="error")
        return(NULL)
      }
      
      event_row <- df[df$ID == event_id, ][1,]
      
      # SANITIZE: Ensure coordinates are valid integers, not NA
      for(j in 1:ncol(event_row)) {
        if(is.numeric(event_row[[j]]) && is.na(event_row[[j]])) event_row[[j]] <- 0
      }
      
      # SANITIZE: Gene Symbol for Filenames/CLI
      if ("geneSymbol" %in% names(event_row)) {
         safe_sym <- gsub("[^A-Za-z0-9]", "_", as.character(event_row$geneSymbol))
         safe_sym <- gsub("_+", "_", safe_sym)
         safe_sym <- gsub("^_|_$", "", safe_sym)
         event_row$geneSymbol <- safe_sym
      }
      
      cat("Sashimi Target Event:", event_id, "\n")
      
      withProgress(message = "Generating Sashimi Plot...", detail = paste("Event:", event_id), {
        # Prep directories
        sashimi_dir <- file.path(manager$workspace_dir, "artifacts", "sashimi")
        dir.create(sashimi_dir, showWarnings=FALSE, recursive=TRUE)
        
        # Grouping file
        group_file <- file.path(sashimi_dir, "grouping.gf")
        writeLines(c(
          paste0(paste(b1_files, collapse=","), ":Case"),
          paste0(paste(b2_files, collapse=","), ":Control")
        ), group_file)
        
        # Target Event File
        temp_event_file <- file.path(sashimi_dir, "target_event.txt")
        write.table(event_row, temp_event_file, sep="\t", quote=FALSE, row.names=FALSE)
        
        # Output subfolder
        out_sub <- file.path(sashimi_dir, paste0(input$selected_type, "_", event_id))
        if(dir.exists(out_sub)) unlink(out_sub, recursive=TRUE)
        
        # Resolve executable path
        exe_path <- file.path(getwd(), "env", "bin", "rmats2sashimiplot")
        if (!file.exists(exe_path)) exe_path <- "rmats2sashimiplot" # Fallback to system PATH
        
        env_bin <- file.path(getwd(), "env", "bin")
        
        cmd_str <- paste(
          "export MPLBACKEND=Agg;",
          paste0("export PATH=", env_bin, ":$PATH;"),
          exe_path,
          "--b1", paste(b1_files, collapse=","),
          "--b2", paste(b2_files, collapse=","),
          "--event-type", input$selected_type,
          "-e", temp_event_file,
          "--l1", "Case", "--l2", "Control",
          "--exon_s", 1, "--intron_s", 5,
          "-o", out_sub
        )
        
        # Execute
        res <- system(paste("sh -c '", cmd_str, "' 2>&1"), intern=TRUE)
        
        # Find PDF result
        pdfs <- list.files(file.path(out_sub, "Sashimi_plot"), pattern="\\.pdf$", full.names=TRUE)
        if (length(pdfs) > 0) {
          # Copy PDF to www/tmp for display
          www_tmp <- "www/tmp_sashimi"
          if(!dir.exists(www_tmp)) dir.create(www_tmp, recursive=TRUE)
          target_pdf <- file.path(www_tmp, paste0("sashimi_", as.numeric(Sys.time()), ".pdf"))
          file.copy(pdfs[1], target_pdf, overwrite=TRUE)
          return(target_pdf)
        } else {
          err_msg <- paste(tail(res, 10), collapse = "\n")
          showNotification(paste("Sashimi generation failed:", err_msg), type = "error", duration = 15)
          return(NULL)
        }
      })
    }

    # --- Sashimi Events ---
    observeEvent(input$run_sashimi, {
      req(input$sashimi_event)
      pdf_path <- generate_sashimi_internal(input$sashimi_event)
      if(!is.null(pdf_path)) rv$sashimi_image <- pdf_path
    })
    
    output$popup_sashimi_ui <- renderUI({
       req(rv$popup_event_id)
       
       # Blocking generation (UI will gray out/show built-in progress)
       # We pass is_popup=TRUE to generate_sashimi_internal if needed, 
       # or just reuse the function. It generates PDF and returns path.
       pdf_path <- generate_sashimi_internal(rv$popup_event_id, is_popup=TRUE)
       
       if (is.null(pdf_path)) {
          return(h4("Error: No plot generated."))
       }
       tags$iframe(src = gsub("www/", "", pdf_path), width = "100%", height = "600px")
    })
    
    # Popup Logic
    show_sashimi_modal <- function(event_id) {
       rv$popup_event_id <- event_id
       
       # Lookup Gene Info
       df <- get_dataset()
       title_str <- paste("Sashimi Plot - Event", event_id)
       
       if (!is.null(df)) {
         # Find row. Be careful with types (event_id is numeric)
         row <- df[df$ID == as.numeric(event_id), ]
         if (nrow(row) > 0) {
           sym <- if("geneSymbol" %in% names(row)) row$geneSymbol[1] else "N/A"
           gid <- if("GeneID" %in% names(row)) row$GeneID[1] else "N/A"
           title_str <- paste0("Sashimi Plot - ", sym, " (", gid, ") : Event ", event_id)
         }
       }
       
       showModal(modalDialog(
          title = title_str,
          size = "l",
          uiOutput(ns("popup_sashimi_ui")),
          footer = modalButton("Close")
       ))
    }

    observeEvent(event_data("plotly_click", source = "volcano_source"), {
       d <- event_data("plotly_click", source = "volcano_source")
       if(!is.null(d$key)) {
         show_sashimi_modal(d$key)
         shinyjs::runjs("Shiny.setInputValue('.clientValue-plotly_click-volcano_source', null);")
       }
    })
    
    observeEvent(event_data("plotly_click", source = "scatter_source"), {
       d <- event_data("plotly_click", source = "scatter_source")
       if(!is.null(d$key)) {
         show_sashimi_modal(d$key)
         shinyjs::runjs("Shiny.setInputValue('.clientValue-plotly_click-scatter_source', null);")
       }
    })

    
    # --- Update Choices ---
    observe({
      df <- get_dataset()
      req(df)
      
      # 1. Apply Gene Filter
      df <- apply_gene_filter(df, filters$gene_list)
      
      # Filter significant
      sig <- df %>%
        dplyr::filter(FDR < filters$fdr, abs(IncLevelDifference) > filters$dpsi)
      if(nrow(sig) > 0) {
        # Values are the numeric IDs
        choices <- sig$ID
        
        # Robust Label Construction
        labels <- as.character(sig$ID)
        
        has_sym <- "geneSymbol" %in% names(sig)
        has_gid <- "GeneID" %in% names(sig)
        
        if (has_sym && has_gid) {
          labels <- paste0(sig$geneSymbol, " (", sig$GeneID, ") : ", sig$ID)
        } else if (has_gid) {
          labels <- paste0(sig$GeneID, " : ", sig$ID)
        } else if (has_sym) {
          labels <- paste0(sig$geneSymbol, " : ", sig$ID)
        }
        
        names(choices) <- labels
        updateSelectInput(session, "sashimi_event", choices = choices)
      }
    })
    
    # --- UI Components ---
    output$main_ui <- renderUI({
      if (rv$status != "loaded") {
        return(box(title="No Data", width=12, status="warning",
                   "No rMATS results loaded.",
                   textInput(ns("manual_path"), "Enter Directory Path manually:"),
                   actionButton(ns("load_manual"), "Load")
        ))
      }
      
      tagList(
        fluidRow(
          box(title = tagList("Controls", info_icon("Filter and annotate detected splicing events based on significance and gene lists.")), width=12, status="primary", collapsible=TRUE,
              column(2, 
                selectInput(ns("selected_type"), label = tagList("Event Type", info_icon("Choose the type of splicing event to analyze (e.g., SE: Skipped Exon).")), choices=names(rv$data)),
                div(style="display: flex; gap: 10px; margin-top: -10px;",
                  checkboxInput(ns("filter_novel"), label = tagList("Novel", info_icon("Include splicing events not present in the GTF annotation.")), value = TRUE),
                  checkboxInput(ns("filter_known"), label = tagList("Known", info_icon("Include splicing events matching the GTF annotation.")), value = TRUE)
                )
              ),
              column(2, numericInput(ns("fdr_cut"), label = tagList("FDR <", info_icon("What it does: Filters events based on statistical significance.\nHow it works: Uses the Benjamini-Hochberg adjusted p-value (FDR).\nWhen to use: Use 0.05 (standard) or 0.01 (strict) to minimize false positives.")), 0.05, 0, 1, 0.01)),
              column(2, numericInput(ns("dpsi_cut"), label = tagList("DeltaPSI >", info_icon("What it does: Filters events by the magnitude of splicing change.\nHow it works: Calculates the absolute difference in mean Percent Spliced In (PSI) between Case and Control.\nWhen to use: Use 0.1 for most cases, or 0.2 to focus on large biological effects.")), 0.1, 0, 1, 0.05)),
              column(3, textAreaInput(ns("gene_filter"), label = tagList("Filter Genes", info_icon("What it does: Restricts plots and tables to specific genes of interest.\nHow it works: Matches your input against 'geneSymbol' and 'GeneID' columns.")), placeholder="e.g. GAPDH, BRCA1...", rows=1, resize="none")),
              column(3, style="margin-top: 25px; display:flex; gap:5px;", 
                     actionButton(ns("apply_filters"), "Apply", class = "btn-primary btn-sm", icon=icon("refresh")),
                     actionButton(ns("annotate_btn"), "Annotate", class = "btn-info btn-sm", icon=icon("tag"))
              ),
              column(12, h5("Status: Data Loaded"))
          ),
          uiOutput(ns("filter_alert")),
        ),
        fluidRow(box(title=tagList("Global Landscape", info_icon("What it does: Shows the overall distribution of significant splicing events.\nHow it works: Aggregates counts of 'Included' (DeltaPSI > 0) vs 'Skipped' (DeltaPSI < 0) events for all types.")), width=12, 
          div(style="display:flex; justify-content:flex-end; gap:6px; margin-bottom: 5px;",
              actionButton(ns("save_global"), "Save Plots", icon=icon("save"), class="btn-xs"),
              actionButton(ns("save_global_data"), "Save Data", icon=icon("table"), class="btn-xs btn-default")
          ),
          plotly::plotlyOutput(ns("global_plot"), height="250px")
        ),
        box(title = tagList("Gene-Level Event Heatmap", info_icon("What it does: Identifies genes with complex splicing regulation (multiple events).\nHow it works: Counts significant events per gene across different event types and displays them as a density heatmap.")), width = 12, collapsible = TRUE, collapsed = FALSE,
          div(style="display:flex; justify-content:space-between; align-items:center; margin-bottom: 5px;",
              p("Heatmap showing the number of significant differential splicing events per gene across different event types.", style="margin:0;"),
              div(style="display:flex; gap:6px;",
                  actionButton(ns("save_gene_heatmap"), "Save Plots", icon=icon("save"), class="btn-xs"),
                  actionButton(ns("save_gene_heatmap_data"), "Save Data", icon=icon("table"), class="btn-xs btn-default")
              )
          ),
          plotly::plotlyOutput(ns("gene_event_heatmap"), height = "500px")
        )),
        tabsetPanel(
          tabPanel("Volcano & Scatter", 
            fluidRow(
              column(6, 
                h4(tagList("Volcano Plot", info_icon("What it does: Visualizes the trade-off between effect size (DeltaPSI) and significance (FDR).\nHow it works: Plots DeltaPSI (x) against -log10(FDR) (y). Significant events appear in the top-left and top-right corners."))),
                div(style="display:flex; justify-content:flex-end; gap:6px; margin-bottom: 5px;",
                    actionButton(ns("save_volcano"), "Save Plots", icon=icon("save"), class="btn-xs"),
                    actionButton(ns("save_volcano_data"), "Save Data", icon=icon("table"), class="btn-xs btn-default")
                ),
                plotly::plotlyOutput(ns("volcano"), height="600px")
              ),
              column(6, 
                h4(tagList("Scatter Plot", info_icon("What it does: Compares splicing inclusion levels between groups.\nHow it works: Plots Case Mean PSI (y) vs Control Mean PSI (x). Events far from the diagonal are differentially spliced."))),
                div(style="display:flex; justify-content:flex-end; gap:6px; margin-bottom: 5px;",
                    actionButton(ns("save_scatter"), "Save Plots", icon=icon("save"), class="btn-xs"),
                    actionButton(ns("save_scatter_data"), "Save Data", icon=icon("table"), class="btn-xs btn-default")
                ),
                plotly::plotlyOutput(ns("scatter"), height="600px")
              )
            )
          ),
          tabPanel("Heatmap", 
            fluidRow(
              box(title = tagList("Heatmap Controls", info_icon("What it does: Visualizes PSI variation across all samples for the top events.\nHow it works: Displays individual sample PSI values in a grid, with optional normalization and clustering.")), width=12, status="info", collapsible=TRUE,
                  column(2, numericInput(ns("heatmap_top_n"), label = tagList("Top N Events", info_icon("What it does: Limits the number of rows in the heatmap.\nHow it works: Picks the top N events based on the 'Event Order' setting (FDR or DeltaPSI).")), value = 50, min = 10, max = 500, step = 10)),
                  column(2, selectInput(ns("heatmap_mode"), label = tagList("Color Mode", info_icon("What it does: Changes how PSI values are colored.\nHow it works: 'Raw PSI' (0-1) shows absolute levels. 'Row Z-score' shows relative changes by centering each row to 0 and scaling to variance.")), choices = c("Raw PSI" = "raw", "Row Z-score" = "zscore"), selected = "raw")),
                  column(3, selectInput(ns("heatmap_order"), label = tagList("Event Order", info_icon("What it does: Determines the sorting of events (rows).\nHow it works: Ranks events by p-value (FDR) or by the magnitude of PSI change (DeltaPSI).")), choices = c("FDR (ascending)" = "fdr", "|DeltaPSI| (descending)" = "dpsi"), selected = "fdr")),
                  column(2, checkboxInput(ns("heatmap_cluster_rows"), label = tagList("Cluster Rows", info_icon("What it does: Groups similar events together.\nHow it works: Performs hierarchical clustering based on PSI similarity across samples.")), value = TRUE)),
                  column(2, checkboxInput(ns("heatmap_cluster_cols"), label = tagList("Cluster Columns", info_icon("What it does: Groups similar samples together.\nHow it works: Performs hierarchical clustering based on overall splicing profiles.")), value = FALSE)),
                  column(1, tags$div(style = "margin-top: 25px;", tags$span(class = "text-muted", "NA=gray")))
              )
            ),
            div(style="display:flex; justify-content:flex-end; gap:6px; margin-bottom: 5px;",
                actionButton(ns("save_heatmap"), "Save Plots", icon=icon("save"), class="btn-xs"),
                actionButton(ns("save_heatmap_data"), "Save Data", icon=icon("table"), class="btn-xs btn-default")
            ),
            plotly::plotlyOutput(ns("heatmap"), height="600px")
          ),
          tabPanel("PCA", 
            fluidRow(
              box(title=tagList("PCA Controls", info_icon("What it does: Projects high-dimensional splicing data into a 2D/3D space.\nHow it works: Identifies 'Principal Components' that explain the most variance in PSI values across all samples.")), width=12, status="info", collapsible=TRUE,
                  fluidRow(
                    column(3, numericInput(ns("pca_top_n"), label = tagList("Top N Features", info_icon("What it does: Selects the most informative events for PCA.\nHow it works: Ranks events by their PSI variance across all samples.")), value = 500, min = 50, max = 5000, step = 50)),
                    column(3, numericInput(ns("pca_max_na"), label = tagList("Max Missing (%)", info_icon("What it does: Filters out events with too many missing values.\nHow it works: Calculates the percentage of samples with NA (low coverage) for each event.")), value = 30, min = 0, max = 100, step = 10)),
                    column(3, numericInput(ns("pca_min_dpsi"), label = tagList("Min PSI Range", info_icon("What it does: Filters for events that actually change across samples.\nHow it works: Subtracts the minimum PSI from the maximum PSI for each event.")), value = 0.1, min = 0, max = 1, step = 0.05)),
                    column(3, numericInput(ns("pca_min_reads"), label = tagList("Min Avg Reads", info_icon("What it does: Filters for events with reliable sequencing depth.\nHow it works: Calculates average (IJC + SJC) counts across all samples.")), value = 10, min = 0, max = 500, step = 5))
                  ),
                  fluidRow(
                    column(3, selectInput(ns("pca_impute"), label = tagList("Imputation", info_icon("What it does: Handles remaining missing values (NA).\nHow it works: 'Row Mean' replaces NA with the average PSI of that event. 'Zero' replaces NA with 0.")), choices=c("Row Mean"="mean", "Zero (0.0)"="zero"))),
                    column(2, checkboxInput(ns("pca_scale"), label = tagList("Scale Data", info_icon("What it does: Standardizes values (Z-score).\nHow it works: Ensures each event has equal weight in PCA, regardless of its absolute PSI levels.")), value = FALSE)),
                    column(2, checkboxInput(ns("pca_3d"), label = tagList("3D Plot", info_icon("What it does: Adds a third dimension (PC3).\nHow it works: Renders an interactive 3D scatter plot.")), value = FALSE)),
                    column(3, checkboxInput(ns("pca_ellipse"), label = tagList("Confidence Ellipse", info_icon("What it does: Draws 95% confidence intervals around groups.\nHow it works: Assumes multivariate normal distribution for each group.")), value = TRUE)),
                    column(2, actionButton(ns("apply_pca"), "Apply Filters", class = "btn-primary btn-sm", style = "margin-top: 25px;"))
                  )
              )
            ),
            div(style="display:flex; justify-content:flex-end; gap:6px; margin-bottom: 5px;",
                actionButton(ns("save_pca"), "Save Plots", icon=icon("save"), class="btn-xs"),
                actionButton(ns("save_pca_data"), "Save Data", icon=icon("table"), class="btn-xs btn-default")
            ),
            verbatimTextOutput(ns("pca_msg")),
            tabsetPanel(id=ns("pca_tabs"),
               tabPanel("PCA Plot", plotly::plotlyOutput(ns("pca"), height="500px")),
               tabPanel("Scree Plot", plotly::plotlyOutput(ns("scree_plot"), height="500px"))
            )
          ),
          tabPanel("Forest Plot", 
             h4(tagList("Forest Plot", info_icon("What it does: Shows the point estimate and variance for the top 20 significant events.\nHow it works: Plots Case and Control mean PSI with error bars for each event side-by-side."))),
             div(style="display:flex; justify-content:space-between; align-items:center; margin-bottom: 5px;",
                 div(style="width: 250px;",
                     selectInput(ns("forest_error"), label = tagList("Error Bars", info_icon("What it does: Controls the whiskers on the plot.\nHow it works: 'Confidence Interval' shows the 95% CI around the mean. 'PSI Range' shows the actual min and max values observed.")), choices=c("Confidence Interval (95%)"="ci", "PSI Range (Min-Max)"="range"), selected="ci", width="100%")
                 ),
                 div(style="display:flex; gap:6px;",
                     actionButton(ns("save_forest"), "Save Plots", icon=icon("save"), class="btn-xs"),
                     actionButton(ns("save_forest_data"), "Save Data", icon=icon("table"), class="btn-xs btn-default")
                 )
             ),
             plotly::plotlyOutput(ns("forest_plot"), height="600px")
          ),
          tabPanel("Sashimi Plot (Read-level)",
            sidebarLayout(position = "right",
              sidebarPanel(width = 4,
                h4(tagList("Sashimi Settings", info_icon("What it does: Visualizes raw sequencing read coverage and junction counts for a specific event.\nHow it works: Uses 'rmats2sashimiplot' to extract read data from BAM files and render genomic arcs."))),
                selectInput(ns("sashimi_event"), label = tagList("Select Significant Event", info_icon("What it does: Choose which event to visualize.\nHow it works: Lists significant events (FDR < cut, DeltaPSI > cut) for the current type.")), choices = NULL),
                actionButton(ns("run_sashimi"), "Generate Plot", class = "btn-success btn-block", icon = icon("image")),
                tags$div(style = "display:flex; gap:6px; margin-top: 8px;",
                  actionButton(ns("save_sashimi"), "Save Plot", icon = icon("save"), class = "btn-info btn-sm"),
                  actionButton(ns("save_sashimi_data"), "Save Data", icon = icon("table"), class = "btn-default btn-sm")
                ),
                hr(),
                p(class="text-muted", "Generates plot from BAM files.")
              ),
              mainPanel(width = 8,
                uiOutput(ns("sashimi_display"))
              )
            )
          ),
          tabPanel("Data Table", 
             h4(tagList("Results Table", info_icon("What it does: Provides a searchable, sortable list of all significant events.\nHow it works: Direct output from rMATS filtered by your criteria."))),
             div(style="display:flex; justify-content:flex-end; margin-bottom: 5px;",
                 actionButton(ns("save_main_table"), "Save Table", icon=icon("table"), class="btn-sm btn-default")
             ),
             DTOutput(ns("table"))
          ),
          tabPanel("Debug Info", verbatimTextOutput(ns("debug_txt")))
        ),
        
        # --- Context Menu UI ---
        tags$div(id = ns("ctx_menu"), style = "display:none; position:absolute; z-index:1000; background:white; border:1px solid #ccc; padding:5px; box-shadow: 2px 2px 5px rgba(0,0,0,0.2); cursor:pointer;",
          tags$div(id = ns("ctx_menu_item"), style="padding: 5px 10px;", icon("image"), " Plot Sashimi")
        ),
        
        # --- Context Menu Script ---
        tags$script(HTML(sprintf("
          $(document).on('contextmenu', '#%s .dataTable tbody tr', function(e) {
            e.preventDefault();
            // Try to find the first cell text. 
            // Note: DataTables creates a complex structure, but usually td:first works.
            var id = $(this).find('td:first').text();
            if(id) {
              $('#%s').css({top: e.pageY, left: e.pageX}).data('evt', id).show();
            }
          });
          
          $(document).on('click', function() { 
             $('#%s').hide(); 
          });
          
          $('#%s').on('click', function(e) {
            e.stopPropagation(); // Prevent document click from hiding immediately? 
            var id = $('#%s').data('evt');
            // Send input to Shiny. We must NOT include the namespace prefix for the input name 
            // in setInputValue if we want to read it as input$name inside the module server?
            // actually no, module server input$ expects 'name' but the message must be 'ns-name'.
            // setInputValue('ns-name', val).
            Shiny.setInputValue('%s', id, {priority: 'event'});
            $('#%s').hide();
          });
        ", ns("table"), ns("ctx_menu"), ns("ctx_menu"), ns("ctx_menu_item"), ns("ctx_menu"), ns("sashimi_context_id"), ns("ctx_menu"))))
      )
    })
    
    # --- Context Menu Handler ---
    observeEvent(input$sashimi_context_id, {
       req(input$sashimi_context_id)
       show_sashimi_modal(input$sashimi_context_id)
    })
    
    # --- Plot Objects (Reactive for Display & Download) ---

    global_plot_data <- reactive({
      df <- global_summary()
      if (nrow(df) == 0) return(NULL)
      df
    })

    global_plot_obj <- reactive({
      df <- global_plot_data()
      if (is.null(df)) return(NULL)
      ggplot(df, aes(x = Type, y = n, fill = Direction)) + geom_bar(stat = "identity") + theme_minimal()
    })

    volcano_plot_data <- reactive({
      df <- get_dataset()
      req(df)
      
      df <- apply_gene_filter(df, filters$gene_list)
      
      if (nrow(df) == 0) return(NULL)

      if (!"geneSymbol" %in% names(df)) df$geneSymbol <- "N/A"
      if (!"GeneID" %in% names(df)) df$GeneID <- "N/A"
      
      df <- df %>%
        mutate(
          Sig = ifelse(FDR < filters$fdr & abs(IncLevelDifference) > filters$dpsi, "Sig", "NS")
        )

      # Smart Downsampling: Trigger only if > 50,000 points
      if (nrow(df) > 50000) {
        # Keep ALL significant points and boundary points (FDR < 0.2)
        df_keep <- df %>% dplyr::filter(Sig == "Sig" | FDR < 0.2)
        df_background <- df %>% dplyr::filter(Sig == "NS" & FDR >= 0.2)
        
        # Calculate how many background points we need to hit exactly 48,000 total 
        # ( safely in the 45k-50k range requested )
        target_total <- 48000
        needed_bg <- target_total - nrow(df_keep)
        
        # If we somehow need more background than exists (rare), take all. 
        # Otherwise, take a random sample to hit exactly 48,000 total.
        if (needed_bg > 0 && needed_bg < nrow(df_background)) {
           set.seed(42) # For visual consistency on redraws
           df_background <- df_background[sample(nrow(df_background), needed_bg), ]
        }
        
        df <- bind_rows(df_keep, df_background)
      }

      # Dynamic Ceiling Logic for FDR = 0
      # 1. Find the smallest non-zero FDR
      min_real_fdr <- min(df$FDR[df$FDR > 0], na.rm = TRUE)
      
      # 2. If min_real_fdr is NA or Inf (meaning ALL FDRs are 0 or no valid FDRs exist),
      #    fallback to a reasonable small number (e.g., 1e-10) to avoid breaking.
      if (is.na(min_real_fdr) || is.infinite(min_real_fdr)) min_real_fdr <- 1e-10
      
      # 3. Create a pseudo-zero value that is slightly smaller than the smallest real FDR.
      #    Multiplying by 0.1 shifts it up nicely by 1 log10 unit visually above the rest.
      pseudo_zero_fdr <- min_real_fdr * 0.1

      df %>%
        mutate(
          # Apply the pseudo-zero instead of relying on Inf or 1e-300
          negLog10FDR = -log10(ifelse(FDR == 0, pseudo_zero_fdr, FDR)),
          Hover = paste0("ID: ", ID, "<br>",
                         "Gene: ", geneSymbol, " (", GeneID, ")<br>",
                         "dPSI: ", round(IncLevelDifference, 4), "<br>",
                         "FDR: ", ifelse(FDR == 0, "0.0 (Ceiling)", formatC(FDR, format = "e", digits = 3)))
        )
    })

    volcano_plot_obj <- reactive({
      df <- volcano_plot_data()
      if (is.null(df)) return(NULL)

      y_cut <- -log10(filters$fdr)
      x_cut <- filters$dpsi

      ggplot(df, aes(x = IncLevelDifference, y = negLog10FDR, color = Sig, text = Hover, key = ID)) +
        geom_point(alpha = 0.6) +
        geom_hline(yintercept = y_cut, linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(-x_cut, x_cut), linetype = "dashed", color = "red") +
        theme_minimal()
    })

    scatter_plot_data <- reactive({
      df <- get_dataset()
      req(df)
      
      df <- apply_gene_filter(df, filters$gene_list)
      
      if (nrow(df) == 0) return(NULL)

      if (nrow(df) > 10000) {
        df <- df %>% arrange(FDR) %>% head(10000)
      }

      calc_mean <- function(col_str) {
        sapply(strsplit(as.character(col_str), ","), function(x) {
          x <- x[x != "NA"]
          if (length(x) == 0) return(NA)
          mean(as.numeric(x), na.rm = TRUE)
        })
      }

      col_inc1 <- grep("IncLevel1", names(df), value = TRUE, ignore.case = TRUE)[1]
      col_inc2 <- grep("IncLevel2", names(df), value = TRUE, ignore.case = TRUE)[1]
      if (is.na(col_inc1) || is.na(col_inc2)) return(NULL)

      df$PSI_Case <- calc_mean(df[[col_inc1]])
      df$PSI_Ctrl <- calc_mean(df[[col_inc2]])
      
      # Ensure bounds 0-1 (rare artifacts)
      df$PSI_Case <- pmin(pmax(df$PSI_Case, 0), 1)
      df$PSI_Ctrl <- pmin(pmax(df$PSI_Ctrl, 0), 1)
      
      # Remove NAs which cause lag in plotly WebGL k-d tree
      df <- df %>% dplyr::filter(!is.na(PSI_Case) & !is.na(PSI_Ctrl))

      if (!"geneSymbol" %in% names(df)) df$geneSymbol <- "N/A"
      if (!"GeneID" %in% names(df)) df$GeneID <- "N/A"

      df %>%
        mutate(
          Sig = ifelse(FDR < filters$fdr & abs(IncLevelDifference) > filters$dpsi, "Sig", "NS"),
          Hover = paste0("ID: ", ID, "<br>",
                         "Gene: ", geneSymbol, " (", GeneID, ")<br>",
                         "dPSI: ", round(IncLevelDifference, 4), "<br>",
                         "FDR: ", formatC(FDR, format = "e", digits = 3))
        )
    })

    scatter_plot_obj <- reactive({
      df <- scatter_plot_data()
      if (is.null(df)) return(NULL)

      d <- min(max(filters$dpsi, 0), 1)

      ggplot(df, aes(x = PSI_Ctrl, y = PSI_Case, color = Sig, text = Hover, key = ID)) +
        geom_point(alpha = 0.6) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
        geom_segment(aes(x = 0, y = d, xend = 1 - d, yend = 1), linetype = "dashed", color = "red", inherit.aes = FALSE) +
        geom_segment(aes(x = d, y = 0, xend = 1, yend = 1 - d), linetype = "dashed", color = "red", inherit.aes = FALSE) +
        labs(x = "Control Mean PSI", y = "Case Mean PSI", title = "PSI Correlation") +
        theme_minimal() +
        coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))
    })
    
    get_heatmap_data <- reactive({
      df <- get_dataset()
      req(df)
      
      df <- apply_gene_filter(df, filters$gene_list)
      
      if (nrow(df) == 0) return(NULL)

      col_inc1 <- grep("IncLevel1", names(df), value = TRUE, ignore.case = TRUE)[1]
      col_inc2 <- grep("IncLevel2", names(df), value = TRUE, ignore.case = TRUE)[1]
      if (is.na(col_inc1) || is.na(col_inc2)) return(NULL)

      parse_psi <- function(col_str) {
        vals <- trimws(unlist(strsplit(as.character(col_str), ",")))
        vals[vals %in% c("", "NA", "NaN", "nan")] <- NA
        suppressWarnings(as.numeric(vals))
      }

      df_sig <- df %>%
        dplyr::filter(FDR < filters$fdr, abs(IncLevelDifference) > filters$dpsi)
      if (nrow(df_sig) == 0) return(NULL)

      order_mode <- input$heatmap_order
      if (identical(order_mode, "dpsi")) {
        df_sig <- df_sig %>% dplyr::arrange(dplyr::desc(abs(IncLevelDifference)), FDR)
      } else {
        df_sig <- df_sig %>% dplyr::arrange(FDR, dplyr::desc(abs(IncLevelDifference)))
      }

      n_take <- min(as.integer(input$heatmap_top_n), nrow(df_sig))
      if (is.na(n_take) || n_take < 1) n_take <- min(50, nrow(df_sig))
      df_sig <- head(df_sig, n_take)

      has_sym <- "geneSymbol" %in% names(df_sig)
      has_gid <- "GeneID" %in% names(df_sig)
      if (has_sym && has_gid) {
        event_labels <- paste0(df_sig$geneSymbol, " (", df_sig$GeneID, ") : ", df_sig$ID)
      } else if (has_gid) {
        event_labels <- paste0(df_sig$GeneID, " : ", df_sig$ID)
      } else if (has_sym) {
        event_labels <- paste0(df_sig$geneSymbol, " : ", df_sig$ID)
      } else {
        event_labels <- as.character(df_sig$ID)
      }
      event_labels <- make.unique(event_labels)

      case_list <- lapply(df_sig[[col_inc1]], parse_psi)
      ctrl_list <- lapply(df_sig[[col_inc2]], parse_psi)
      n_case <- max(lengths(case_list), 0)
      n_ctrl <- max(lengths(ctrl_list), 0)
      if ((n_case + n_ctrl) == 0) return(NULL)

      sample_meta <- get_sample_names()
      c_names <- if(length(sample_meta$case) == n_case) sample_meta$case else paste0("Case_", seq_len(n_case))
      ct_names <- if(length(sample_meta$ctrl) == n_ctrl) sample_meta$ctrl else paste0("Ctrl_", seq_len(n_ctrl))
      sample_names <- c(c_names, ct_names)

      mat_raw <- matrix(NA_real_, nrow = nrow(df_sig), ncol = length(sample_names),
                        dimnames = list(event_labels, sample_names))

      for (i in seq_len(nrow(df_sig))) {
        if (length(case_list[[i]]) > 0) {
          mat_raw[i, seq_len(min(length(case_list[[i]]), n_case))] <- case_list[[i]][seq_len(min(length(case_list[[i]]), n_case))]
        }
        if (length(ctrl_list[[i]]) > 0) {
          ctrl_idx <- n_case + seq_len(min(length(ctrl_list[[i]]), n_ctrl))
          mat_raw[i, ctrl_idx] <- ctrl_list[[i]][seq_len(min(length(ctrl_list[[i]]), n_ctrl))]
        }
      }

      mode <- input$heatmap_mode
      mat_display <- mat_raw
      midpoint <- 0.5
      limits <- c(0, 1)
      legend_title <- "PSI"

      if (identical(mode, "zscore")) {
        mat_display <- t(scale(t(mat_raw)))
        mat_display[is.nan(mat_display) | is.infinite(mat_display)] <- NA
        max_abs <- max(abs(mat_display), na.rm = TRUE)
        if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
        midpoint <- 0
        limits <- c(-max_abs, max_abs)
        legend_title <- "Z-score"
      }

      row_order <- rownames(mat_display)
      col_order <- colnames(mat_display)

      make_imputed <- function(m) {
        m2 <- m
        for (r in seq_len(nrow(m2))) {
          row_mean <- mean(m2[r, ], na.rm = TRUE)
          if (is.nan(row_mean)) row_mean <- 0
          m2[r, is.na(m2[r, ])] <- row_mean
        }
        m2
      }

      if (isTRUE(input$heatmap_cluster_rows) && nrow(mat_display) > 1) {
        mat_row <- make_imputed(mat_display)
        row_order <- rownames(mat_row)[hclust(dist(mat_row))$order]
      }

      if (isTRUE(input$heatmap_cluster_cols) && ncol(mat_display) > 1) {
        mat_col <- make_imputed(mat_display)
        col_order <- colnames(mat_col)[hclust(dist(t(mat_col)))$order]
      }

      mat_display <- mat_display[row_order, col_order, drop = FALSE]
      mat_raw <- mat_raw[row_order, col_order, drop = FALSE]

      long_display <- as.data.frame(as.table(mat_display), stringsAsFactors = FALSE)
      names(long_display) <- c("EventID", "Sample", "Value")

      long_raw <- as.data.frame(as.table(mat_raw), stringsAsFactors = FALSE)
      names(long_raw) <- c("EventID", "Sample", "PSI")

      long <- dplyr::left_join(long_display, long_raw, by = c("EventID", "Sample"))
      long$Group <- ifelse(long$Sample %in% c_names, "Case", "Control")

      event_meta <- data.frame(
        EventID = event_labels,
        EventNumericID = as.character(df_sig$ID),
        Gene = if ("geneSymbol" %in% names(df_sig)) as.character(df_sig$geneSymbol) else "N/A",
        GeneID = if ("GeneID" %in% names(df_sig)) as.character(df_sig$GeneID) else "N/A",
        dPSI = as.numeric(df_sig$IncLevelDifference),
        FDR = as.numeric(df_sig$FDR),
        stringsAsFactors = FALSE
      )
      long <- dplyr::left_join(long, event_meta, by = "EventID")

      value_label <- if (identical(mode, "zscore")) "Z-score" else "PSI"
      long$Hover <- paste0(
        "Event: ", long$EventNumericID,
        "<br>Gene: ", long$Gene, " (", long$GeneID, ")",
        "<br>Sample: ", long$Sample, " [", long$Group, "]",
        "<br>PSI: ", ifelse(is.na(long$PSI), "NA", sprintf("%.3f", long$PSI)),
        "<br>", value_label, ": ", ifelse(is.na(long$Value), "NA", sprintf("%.3f", long$Value)),
        "<br>dPSI: ", ifelse(is.na(long$dPSI), "NA", sprintf("%.3f", long$dPSI)),
        "<br>FDR: ", ifelse(is.na(long$FDR), "NA", formatC(long$FDR, format = "e", digits = 3))
      )

      split_idx <- if (n_case > 0 && n_ctrl > 0) {
        case_positions <- which(colnames(mat_display) %in% c_names)
        if (length(case_positions) > 0) max(case_positions) else NA_integer_
      } else {
        NA_integer_
      }

      list(
        long = long,
        split_idx = split_idx,
        midpoint = midpoint,
        limits = limits,
        legend_title = legend_title,
        mode = mode,
        n_events = nrow(mat_display),
        n_samples = ncol(mat_display)
      )
    })

    heatmap_plot_data <- reactive({
      hm <- get_heatmap_data()
      if (is.null(hm)) return(NULL)
      out <- hm$long
      out$HeatmapMode <- ifelse(hm$mode == "zscore", "Row Z-score", "Raw PSI")
      out
    })

    heatmap_plot_obj <- reactive({
      hm <- get_heatmap_data()
      if (is.null(hm)) return(NULL)

      p <- ggplot(hm$long, aes(x = Sample, y = EventID, fill = Value, text = Hover)) +
        geom_tile(color = NA, linewidth = 0.1) +
        scale_fill_gradient2(
          low = "#2166AC",
          mid = "#F7F7F7",
          high = "#B2182B",
          midpoint = hm$midpoint,
          limits = hm$limits,
          na.value = "#BDBDBD",
          name = hm$legend_title
        ) +
        labs(
          title = paste0("Top ", hm$n_events, " Significant Events (", ifelse(hm$mode == "zscore", "Row Z-score", "Raw PSI"), ")"),
          x = "Samples",
          y = "Events"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 8)
        )

      if (!is.na(hm$split_idx)) {
        p <- p + geom_vline(xintercept = hm$split_idx + 0.5, color = "black", linewidth = 0.4)
      }

      p
    })

    pca_plot_data <- eventReactive(list(input$apply_pca, input$selected_type, manager$events$workspace_loaded), {
      df <- get_dataset()
      if (is.null(df)) {
        output$pca_msg <- renderText("No PSI data found.")
        return(NULL)
      }
      
      df <- apply_gene_filter(df, filters$gene_list)
      if (nrow(df) == 0) return(NULL)

      col_inc1 <- grep("IncLevel1", names(df), value=TRUE, ignore.case=TRUE)[1]
      col_inc2 <- grep("IncLevel2", names(df), value=TRUE, ignore.case=TRUE)[1]
      if (is.na(col_inc1) || is.na(col_inc2)) return(NULL)

      # Fast parsing
      inc1_list <- strsplit(as.character(df[[col_inc1]]), ",")
      inc2_list <- strsplit(as.character(df[[col_inc2]]), ",")
      
      if (length(inc1_list) == 0 || length(inc2_list) == 0) return(NULL)
      
      n_reps1 <- length(inc1_list[[1]])
      n_reps2 <- length(inc2_list[[1]])
      
      if (n_reps1 + n_reps2 < 3) {
         output$pca_msg <- renderText("Insufficient samples for 2D PCA. Need at least 3 samples.")
         return(NULL)
      }
      
      pad_na <- function(lst, target_len) {
          lapply(lst, function(x) {
             if (length(x) == target_len) return(x)
             if (length(x) > target_len) return(x[1:target_len])
             return(c(x, rep("NA", target_len - length(x))))
          })
      }
      
      # FILTER A: Min Average Read Depth
      min_reads <- as.numeric(input$pca_min_reads)
      if (is.na(min_reads)) min_reads <- 0
      
      col_ijc1 <- "IJC_SAMPLE_1"; col_sjc1 <- "SJC_SAMPLE_1"
      col_ijc2 <- "IJC_SAMPLE_2"; col_sjc2 <- "SJC_SAMPLE_2"
      
      if (min_reads > 0 && all(c(col_ijc1, col_sjc1, col_ijc2, col_sjc2) %in% names(df))) {
         sum_counts <- function(str_vec) {
             sapply(strsplit(as.character(str_vec), ","), function(x) {
                nums <- suppressWarnings(as.numeric(x))
                sum(nums[!is.na(nums)])
             })
         }
         
         total_reads <- sum_counts(df[[col_ijc1]]) + sum_counts(df[[col_sjc1]]) + 
                        sum_counts(df[[col_ijc2]]) + sum_counts(df[[col_sjc2]])
         avg_reads <- total_reads / (n_reps1 + n_reps2)
         
         keep_idx <- avg_reads >= min_reads
         df <- df[keep_idx, , drop=FALSE]
         
         if (nrow(df) == 0) {
            output$pca_msg <- renderText("No events passed the Minimum Average Read Count filter.")
            return(NULL)
         }
         # Re-parse lists after subsetting
         inc1_list <- strsplit(as.character(df[[col_inc1]]), ",")
         inc2_list <- strsplit(as.character(df[[col_inc2]]), ",")
      }

      inc1_list <- pad_na(inc1_list, n_reps1)
      inc2_list <- pad_na(inc2_list, n_reps2)

      mat1 <- suppressWarnings(matrix(as.numeric(unlist(inc1_list)), nrow = nrow(df), byrow = TRUE))
      mat2 <- suppressWarnings(matrix(as.numeric(unlist(inc2_list)), nrow = nrow(df), byrow = TRUE))
      mat <- cbind(mat1, mat2)

      sample_meta <- get_sample_names()
      c_names <- if(length(sample_meta$case) == n_reps1) sample_meta$case else paste0("Case_", seq_len(n_reps1))
      ct_names <- if(length(sample_meta$ctrl) == n_reps2) sample_meta$ctrl else paste0("Ctrl_", seq_len(n_reps2))
      colnames(mat) <- c(c_names, ct_names)
      
      eid <- if("GeneID" %in% names(df)) paste0(df$GeneID, "_", df$ID) else as.character(df$ID)
      rownames(mat) <- make.unique(eid)

      # FILTER B: Max NA percentage
      max_na_frac <- as.numeric(input$pca_max_na) / 100
      if (is.na(max_na_frac)) max_na_frac <- 0.3
      
      na_frac <- rowMeans(is.na(mat))
      mat <- mat[na_frac <= max_na_frac, , drop=FALSE]
      
      if (nrow(mat) == 0) {
         output$pca_msg <- renderText("No valid events after filtering NA rows.")
         return(NULL)
      }

      # IMPUTE missing values
      impute_method <- if(is.null(input$pca_impute)) "mean" else input$pca_impute
      for (i in seq_len(nrow(mat))) {
        if (any(is.na(mat[i, ]))) {
          if (impute_method == "zero") {
             mat[i, is.na(mat[i, ])] <- 0
          } else {
             m <- mean(mat[i, ], na.rm = TRUE)
             mat[i, is.na(mat[i, ])] <- if (is.nan(m)) 0 else m
          }
        }
      }

      # FILTER C: Minimum PSI Range
      min_dpsi <- as.numeric(input$pca_min_dpsi)
      if (is.na(min_dpsi)) min_dpsi <- 0.0
      
      if (min_dpsi > 0) {
         rng <- apply(mat, 1, function(x) max(x) - min(x))
         mat <- mat[rng >= min_dpsi, , drop=FALSE]
      }
      if (nrow(mat) == 0) {
         output$pca_msg <- renderText("No events passed the Minimum PSI Range filter.")
         return(NULL)
      }

      # FILTER D: Constant variance
      vars <- apply(mat, 1, var)
      mat <- mat[vars > 0, , drop=FALSE]
      vars <- vars[vars > 0]
      
      if (nrow(mat) == 0) {
        output$pca_msg <- renderText("All events have zero variance. Cannot compute PCA.")
        return(NULL)
      }
      
      # FILTER E: Top N Most Variable
      top_n <- as.numeric(input$pca_top_n)
      if (is.na(top_n) || top_n < 10) top_n <- 500
      
      if (nrow(mat) > top_n) {
         top_idx <- order(vars, decreasing = TRUE)[1:top_n]
         mat <- mat[top_idx, , drop=FALSE]
      }

      output$pca_msg <- renderText(paste0("PCA computed on Top ", nrow(mat), " most variable features."))

      do_scale <- isTRUE(input$pca_scale)
      pca <- prcomp(t(mat), scale. = do_scale)
      pcs <- as.data.frame(pca$x)

      if (ncol(pcs) < 2) {
        output$pca_msg <- renderText("PCA yielded fewer than 2 dimensions.")
        return(NULL)
      }

      pcs$Sample <- rownames(pcs)
      pcs$Group <- ifelse(pcs$Sample %in% c_names, "Case", "Control")
      
      var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
      pcs$Variance_PC1 <- if (length(var_exp) >= 1) var_exp[1] else NA_real_
      pcs$Variance_PC2 <- if (length(var_exp) >= 2) var_exp[2] else NA_real_
      pcs$Variance_PC3 <- if (length(var_exp) >= 3) var_exp[3] else NA_real_

      attr(pcs, "var_exp") <- var_exp
      pcs
    })

    pca_plot_obj <- reactive({
      pcs <- pca_plot_data()
      if (is.null(pcs)) return(NULL)
      
      var1 <- round(pcs$Variance_PC1[1] * 100, 1)
      var2 <- round(pcs$Variance_PC2[1] * 100, 1)
      
      p <- ggplot(pcs, aes(x = PC1, y = PC2, color = Group)) + 
        geom_point(aes(label = Sample), size = 3) + 
        labs(x = paste0("PC1 (", var1, "%)"), y = paste0("PC2 (", var2, "%)"), title = "PCA Plot") +
        theme_minimal()
        
      if (isTRUE(input$pca_ellipse)) {
         counts <- table(pcs$Group)
         if (any(counts >= 3)) {
            p <- p + stat_ellipse(level = 0.95, alpha=0.5, linetype="dashed")
         }
      }
      p
    })
    
    scree_plot_obj <- reactive({
       pcs <- pca_plot_data()
       if(is.null(pcs)) return(NULL)
       var_exp <- attr(pcs, "var_exp")
       if(is.null(var_exp)) return(NULL)
       
       scree_df <- data.frame(
          PC = factor(paste0("PC", seq_along(var_exp)), levels = paste0("PC", seq_along(var_exp))),
          Variance = var_exp * 100
       )
       
       if (nrow(scree_df) > 10) scree_df <- head(scree_df, 10)
       
       ggplot(scree_df, aes(x = PC, y = Variance)) +
          geom_bar(stat="identity", fill="#4682B4", alpha=0.8) +
          geom_point(color="#B22222", size=2) +
          geom_line(aes(group=1), color="#B22222", linewidth=1) +
          labs(x = "Principal Component", y = "Variance Explained (%)", title = "Scree Plot (Top 10)") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle=45, hjust=1))
    })

    forest_plot_data <- reactive({
      df <- get_dataset()
      req(df)
      
      df <- apply_gene_filter(df, filters$gene_list)

      sig_df <- df %>%
        dplyr::filter(FDR < filters$fdr, abs(IncLevelDifference) > filters$dpsi) %>%
        dplyr::arrange(FDR) %>%
        head(20)

      if (nrow(sig_df) == 0) return(NULL)

      # Ensure dependency on error bar type
      error_type <- input$forest_error
      if (is.null(error_type)) error_type <- "ci"

      get_stats <- function(inc_str, ijc_str, sjc_str) {
         inc <- as.numeric(unlist(strsplit(as.character(inc_str), ",")))
         ijc <- as.numeric(unlist(strsplit(as.character(ijc_str), ",")))
         sjc <- as.numeric(unlist(strsplit(as.character(sjc_str), ",")))

         valid <- !is.na(inc)
         inc <- inc[valid]
         if (length(ijc) >= length(inc)) ijc <- ijc[valid]
         if (length(sjc) >= length(inc)) sjc <- sjc[valid]

         if (length(inc) < 1) return(c(0, 0, 0, 0))

         m <- mean(inc)
         count <- sum(ijc, na.rm = TRUE) + sum(sjc, na.rm = TRUE)
         
         if (error_type == "range") {
           ci_min <- min(inc)
           ci_max <- max(inc)
         } else {
           # Standard 95% CI
           n <- length(inc)
           sd_val <- if (n > 1) sd(inc) else 0
           se <- if (n > 0) sd_val / sqrt(n) else 0
           t_val <- if (n > 1) qt(0.975, df = n - 1) else 0
           
           ci_min <- max(0, m - t_val * se)
           ci_max <- min(1, m + t_val * se)
         }

         return(c(m, ci_min, ci_max, count))
      }

      plot_data <- data.frame()
      col_inc1 <- grep("IncLevel1", names(sig_df), value=TRUE, ignore.case=TRUE)[1]
      col_inc2 <- grep("IncLevel2", names(sig_df), value=TRUE, ignore.case=TRUE)[1]
      col_ijc1 <- "IJC_SAMPLE_1"; col_sjc1 <- "SJC_SAMPLE_1"
      col_ijc2 <- "IJC_SAMPLE_2"; col_sjc2 <- "SJC_SAMPLE_2"

      if (is.na(col_inc1) || is.na(col_inc2)) return(NULL)

      for (i in 1:nrow(sig_df)) {
         nm <- as.character(sig_df$ID[i])
         if ("geneSymbol" %in% names(sig_df)) nm <- paste0(sig_df$geneSymbol[i], " (", nm, ")")

         # Process Case
         ijc1 <- if (col_ijc1 %in% names(sig_df)) sig_df[[col_ijc1]][i] else "0"
         sjc1 <- if (col_sjc1 %in% names(sig_df)) sig_df[[col_sjc1]][i] else "0"
         stats1 <- get_stats(sig_df[[col_inc1]][i], ijc1, sjc1)

         # Stats
         dpsi_val <- sig_df$IncLevelDifference[i]
         fdr_val <- sig_df$FDR[i]

         plot_data <- rbind(plot_data, data.frame(
            event_name = nm, Group = "Case",
            fraction_value = stats1[1], ci_min = stats1[2], ci_max = stats1[3], count_size = stats1[4],
            dpsi = dpsi_val, fdr = fdr_val,
            is_significant = TRUE
          ))

         # Process Control
         ijc2 <- if (col_ijc2 %in% names(sig_df)) sig_df[[col_ijc2]][i] else "0"
         sjc2 <- if (col_sjc2 %in% names(sig_df)) sig_df[[col_sjc2]][i] else "0"
         stats2 <- get_stats(sig_df[[col_inc2]][i], ijc2, sjc2)

         plot_data <- rbind(plot_data, data.frame(
            event_name = nm, Group = "Control",
            fraction_value = stats2[1], ci_min = stats2[2], ci_max = stats2[3], count_size = stats2[4],
            dpsi = dpsi_val, fdr = fdr_val,
            is_significant = TRUE
          ))
      }

      plot_data$event_name <- factor(plot_data$event_name, levels = rev(unique(plot_data$event_name)))
      plot_data$Group <- factor(plot_data$Group, levels = c("Case", "Control"))

      plot_data
    })

    forest_plot_obj <- reactive({
      plot_data <- forest_plot_data()
      if (is.null(plot_data)) return(NULL)
      
      n_events <- length(unique(plot_data$event_name))
      title_str <- paste0("Top ", n_events, " Significant Events (PSI Comparison)")

      plot_forest_dot(plot_data, title = title_str, xlab = "Mean PSI")
    })

    gene_event_heatmap_data <- reactive({
      req(rv$status == "loaded")
      
      all_sig <- data.frame()
      
      # Iterate all loaded event types
      for (t in names(rv$data)) {
        df <- rv$data[[t]]
        
        # Apply Novel/Known filter
        df <- apply_event_filter(df)
        
        # 1. Apply Gene Filter
        df <- apply_gene_filter(df, filters$gene_list)
        
        # 2. Filter Significant
        sig <- df %>% 
          dplyr::filter(FDR < filters$fdr, abs(IncLevelDifference) > filters$dpsi)
        
        if(nrow(sig) > 0) {
          # Ensure Gene Label
          if(!"geneSymbol" %in% names(sig)) sig$geneSymbol <- "N/A"
          if(!"GeneID" %in% names(sig)) sig$GeneID <- "N/A"
          
          # Prioritize Symbol, fallback to ID
          sig$GeneLabel <- ifelse(sig$geneSymbol != "N/A" & sig$geneSymbol != "", 
                                  sig$geneSymbol, 
                                  as.character(sig$GeneID))
          
          # Count events per gene
          tmp <- sig %>% 
            dplyr::count(GeneLabel, name="Count") %>%
            dplyr::mutate(Type = t)
          
          all_sig <- rbind(all_sig, tmp)
        }
      }
      
      if(nrow(all_sig) == 0) return(NULL)
      
      # Logic to limit rows if too many genes (and no filter active)
      is_filtered <- (filters$gene_list != "" && !is.null(filters$gene_list))
      
      if (!is_filtered && length(unique(all_sig$GeneLabel)) > 50) {
        # Take Top 50 Genes by Total Events
        top_genes <- all_sig %>% 
          dplyr::group_by(GeneLabel) %>% 
          dplyr::summarise(Total = sum(Count)) %>% 
          dplyr::arrange(dplyr::desc(Total)) %>% 
          head(50) %>% 
          dplyr::pull(GeneLabel)
        
        all_sig <- all_sig %>% dplyr::filter(GeneLabel %in% top_genes)
      }
      
      # Fill 0s for missing (Gene, Type) combinations
      all_sig <- all_sig %>% 
        tidyr::complete(GeneLabel, Type, fill = list(Count = 0))
        
      # Attach meta-data for plot title logic
      attr(all_sig, "is_filtered") <- is_filtered
      
      all_sig
    })

    gene_event_heatmap_obj <- reactive({
      all_sig <- gene_event_heatmap_data()
      if(is.null(all_sig)) return(NULL)
      
      is_filtered <- isTRUE(attr(all_sig, "is_filtered"))
      
      ggplot(all_sig, aes(x = Type, y = GeneLabel, fill = Count, 
                          text = paste("Gene:", GeneLabel, "<br>Type:", Type, "<br>Events:", Count))) + 
        geom_tile(color = NA, linewidth = 0) +
        scale_fill_viridis_c(name = "# Events") +
        labs(x = "Event Type", y = "Gene", title = if(!is_filtered) "Top 50 Genes with Diff. Splicing" else "Gene-Level Event Counts") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()
        )
    })
    
    output$gene_event_heatmap <- renderPlotly({
      p <- gene_event_heatmap_obj()
      if(is.null(p)) return(plotly::plot_ly() %>% plotly::layout(title = "No significant events found for current filters."))
      ggplotly(p, tooltip = "text")
    })

    # --- Plot Rendering (Display) ---
    output$global_plot <- renderPlotly({
      p <- global_plot_obj()
      if(is.null(p)) return(NULL)
      ggplotly(p)
    })
    
    output$volcano <- renderPlotly({
      p <- volcano_plot_obj()
      if(is.null(p)) return(NULL)
      toWebGL(ggplotly(p, tooltip="text", source = "volcano_source"))
    })
    
    output$scatter <- renderPlotly({
      p <- scatter_plot_obj()
      if(is.null(p)) return(NULL)
      toWebGL(ggplotly(p, tooltip="text", source = "scatter_source"))
    })
    
    output$heatmap <- renderPlotly({
      p <- heatmap_plot_obj()
      if(is.null(p)) { output$heatmap_msg <- renderText("No PSI data found."); return(NULL) }
      output$heatmap_msg <- renderText("")
      ggplotly(p, tooltip = "text")
    })
    
    output$pca <- renderPlotly({
      pcs <- pca_plot_data()
      if(is.null(pcs)) return(plotly::plot_ly() %>% plotly::layout(title="No PCA Data"))
      
      if (isTRUE(input$pca_3d) && "PC3" %in% names(pcs)) {
          var1 <- round(pcs$Variance_PC1[1] * 100, 1)
          var2 <- round(pcs$Variance_PC2[1] * 100, 1)
          var3 <- round(pcs$Variance_PC3[1] * 100, 1)
          
          # Native Plotly 3D scatter
          plotly::plot_ly(pcs, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Group, text = ~Sample, 
                  type = "scatter3d", mode = "markers", marker = list(size = 6)) %>%
            plotly::layout(scene = list(
              xaxis = list(title = paste0("PC1 (", var1, "%)")),
              yaxis = list(title = paste0("PC2 (", var2, "%)")),
              zaxis = list(title = paste0("PC3 (", var3, "%)"))
            ))
      } else {
          # Use ggplotly for 2D (supports stat_ellipse)
          p <- pca_plot_obj()
          if(is.null(p)) return(NULL)
          ggplotly(p, tooltip = c("label", "x", "y"))
      }
    })
    
    output$scree_plot <- renderPlotly({
      p <- scree_plot_obj()
      if(is.null(p)) return(NULL)
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    output$forest_plot <- renderPlotly({
      p <- forest_plot_obj()
      if(is.null(p)) return(NULL)
      
      gp <- ggplotly(p, tooltip="text")
      
      # Fix clipping for annotations (geom_text) outside the plot area
      # 1. Add right margin
      gp <- plotly::layout(gp, margin = list(r = 150))
      
      # 2. Iterate traces to find text traces and disable clipping
      gp$x$data <- lapply(gp$x$data, function(tr) {
        # Check if it's a text trace (geom_text conversion)
        # ggplotly usually converts geom_text to scatter with mode 'text'
        if (!is.null(tr$mode) && grepl("text", tr$mode)) {
          tr$cliponaxis <- FALSE
        }
        tr
      })
      
      gp
    })
    
    # --- Server-Side Plot Saving ---
    save_plot_to_workspace <- function(plot_reactive, plot_base_name) {
      req(plot_reactive())
      req(manager$workspace_dir)
      
      # 1. Prepare Directory
      plot_dir <- file.path(manager$workspace_dir, "artifacts", "plots", "rmats")
      if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
      
      # 2. Construct Filename
      # Include Event Type (SE, MXE...) if available from input
      type_tag <- if(!is.null(input$selected_type)) input$selected_type else "All"
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      base_filename <- paste0(plot_base_name, "_", type_tag, "_", timestamp)
      
      path_svg <- file.path(plot_dir, paste0(base_filename, ".svg"))
      path_png <- file.path(plot_dir, paste0(base_filename, ".png"))
      
      # 3. Save
      p <- plot_reactive()
      
      # Force white background for export (SVG/PNG transparency fix)
      p_export <- p + theme(plot.background = element_rect(fill = "white", color = NA))
      
      withProgress(message = "Saving plots...", detail = base_filename, {
        tryCatch({
          ggsave(path_svg, plot = p_export, device = svg, width = 12, height = 8)
          ggsave(path_png, plot = p_export, device = "png", width = 12, height = 8, dpi = 300)
          
          showNotification(paste("Saved:", base_filename), type = "message")
        }, error = function(e) {
          showNotification(paste("Save failed:", e$message), type = "error")
        })
      })
    }

    save_data_to_workspace <- function(data_reactive, data_base_name) {
      req(manager$workspace_dir)

      data_dir <- file.path(manager$workspace_dir, "artifacts", "data_exports", "rmats")
      if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

      type_tag <- if (!is.null(input$selected_type)) input$selected_type else "All"
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      filename <- paste0(data_base_name, "_", type_tag, "_", timestamp, ".csv")
      target_path <- file.path(data_dir, filename)

      withProgress(message = "Saving data...", detail = filename, {
        tryCatch({
          df <- data_reactive()
          if (is.null(df) || nrow(df) == 0) stop("No data available to save.")
          utils::write.csv(df, target_path, row.names = FALSE, na = "")
          showNotification(paste("Saved data:", filename), type = "message")
        }, error = function(e) {
          showNotification(paste("Data save failed:", e$message), type = "error")
        })
      })
    }

    get_sashimi_export_data <- reactive({
      req(input$sashimi_event)
      df <- get_dataset()
      req(df)

      event_id <- suppressWarnings(as.numeric(input$sashimi_event))
      if (is.na(event_id)) return(NULL)

      event_row <- df[df$ID == event_id, ][1, , drop = FALSE]
      if (nrow(event_row) == 0) return(NULL)

      b1_files <- character(0)
      b2_files <- character(0)
      rmats_params <- manager$get_params("rmats")

      if (!is.null(rmats_params)) {
        if (!is.null(rmats_params$b1) && nzchar(rmats_params$b1)) {
          b1_files <- trimws(unlist(strsplit(rmats_params$b1, ",")))
        }
        if (!is.null(rmats_params$b2) && nzchar(rmats_params$b2)) {
          b2_files <- trimws(unlist(strsplit(rmats_params$b2, ",")))
        }
      }

      if ((length(b1_files) == 0 || length(b2_files) == 0) && !is.null(manager$workspace_dir)) {
        case_txt <- file.path(manager$workspace_dir, "artifacts", "rmats", "input", "case.txt")
        ctrl_txt <- file.path(manager$workspace_dir, "artifacts", "rmats", "input", "control.txt")
        if (file.exists(case_txt)) {
          b1_raw <- paste(readLines(case_txt, warn = FALSE), collapse = "")
          b1_files <- trimws(unlist(strsplit(b1_raw, ",")))
        }
        if (file.exists(ctrl_txt)) {
          b2_raw <- paste(readLines(ctrl_txt, warn = FALSE), collapse = "")
          b2_files <- trimws(unlist(strsplit(b2_raw, ",")))
        }
      }

      event_row$SelectedEvent <- input$sashimi_event
      event_row$EventType <- if (!is.null(input$selected_type)) input$selected_type else "N/A"
      event_row$CaseBamCount <- length(b1_files)
      event_row$ControlBamCount <- length(b2_files)
      event_row$CaseBams <- paste(b1_files, collapse = ",")
      event_row$ControlBams <- paste(b2_files, collapse = ",")
      event_row$Workspace <- if (!is.null(manager$workspace_dir)) manager$workspace_dir else ""
      event_row$SavedSashimiPDF <- if (!is.null(rv$sashimi_image)) rv$sashimi_image else ""

      event_row
    })
    
    observeEvent(input$save_global, { save_plot_to_workspace(global_plot_obj, "global_landscape") })
    observeEvent(input$save_gene_heatmap, { save_plot_to_workspace(gene_event_heatmap_obj, "gene_event_heatmap") })
    observeEvent(input$save_volcano, { save_plot_to_workspace(volcano_plot_obj, "volcano_plot") })
    observeEvent(input$save_scatter, { save_plot_to_workspace(scatter_plot_obj, "scatter_plot") })
    observeEvent(input$save_heatmap, { save_plot_to_workspace(heatmap_plot_obj, "heatmap") })
    observeEvent(input$save_pca, { save_plot_to_workspace(pca_plot_obj, "pca_plot") })
    observeEvent(input$save_forest, { save_plot_to_workspace(forest_plot_obj, "forest_plot") })

    observeEvent(input$save_global_data, { save_data_to_workspace(global_plot_data, "global_landscape_data") })
    observeEvent(input$save_gene_heatmap_data, { save_data_to_workspace(gene_event_heatmap_data, "gene_event_heatmap_data") })
    observeEvent(input$save_volcano_data, { save_data_to_workspace(volcano_plot_data, "volcano_plot_data") })
    observeEvent(input$save_scatter_data, { save_data_to_workspace(scatter_plot_data, "scatter_plot_data") })
    observeEvent(input$save_heatmap_data, { save_data_to_workspace(heatmap_plot_data, "heatmap_data") })
    observeEvent(input$save_main_table, { save_data_to_workspace(filtered_table_data, "filtered_table") })
    observeEvent(input$save_pca_data, { save_data_to_workspace(pca_plot_data, "pca_plot_data") })
    observeEvent(input$save_forest_data, { save_data_to_workspace(forest_plot_data, "forest_plot_data") })

    # --- Sashimi Save Handler ---
    observeEvent(input$save_sashimi, {
      req(rv$sashimi_image)
      req(input$sashimi_event)
      req(manager$workspace_dir)
      
      # Prepare Directory
      plot_dir <- file.path(manager$workspace_dir, "artifacts", "plots", "rmats")
      if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
      
      # Construct Filename
      event_id <- input$sashimi_event
      type_tag <- if(!is.null(input$selected_type)) input$selected_type else "Event"
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      filename <- paste0("Sashimi_Event_", event_id, "_", type_tag, "_", timestamp, ".pdf")
      target_path <- file.path(plot_dir, filename)
      
      # Copy
      if(file.exists(rv$sashimi_image)) {
        file.copy(rv$sashimi_image, target_path, overwrite=TRUE)
        showNotification(paste("Saved Vector PDF:", filename), type = "message")
      } else {
        showNotification("No plot generated to save.", type="warning")
      }
    })

    observeEvent(input$save_sashimi_data, {
      save_data_to_workspace(get_sashimi_export_data, "sashimi_plot_data")
    })
    
    output$sashimi_display <- renderUI({
      if (is.null(rv$sashimi_image)) return(h4("Select an event and click Generate."))
      tags$iframe(src = gsub("www/", "", rv$sashimi_image), width = "100%", height = "600px")
    })
    
    filtered_table_data <- reactive({
      df <- get_dataset()
      req(df, filters$fdr, filters$dpsi)
      
      # Apply Gene Filter
      df <- apply_gene_filter(df, filters$gene_list)
      
      df %>% dplyr::filter(FDR < filters$fdr, abs(IncLevelDifference) > filters$dpsi)
    })
    
    output$table <- renderDT({
      tryCatch({
        filtered_table_data()
      }, error = function(e) {
        cat("Error rendering table:", e$message, "\n")
        data.frame(Error = paste("Error:", e$message))
      })
    }, rownames = FALSE, options=list(scrollX=TRUE))
    
    output$debug_txt <- renderText({
      df <- get_dataset()
      if(is.null(df)) return("No data.")
      paste("Rows:", nrow(df), "\nCols:", paste(names(df), collapse=", "))
    })
  })
}
