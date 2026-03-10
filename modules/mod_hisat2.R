library(shiny)
library(bslib)
library(processx)
library(shinyjs)

mod_hisat2_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(
        title = tagList("Genome Index Management", info_icon("Build or select existing HISAT2 index for alignment. Supports prebuilt downloads for common plant species.")),
        status = "warning", solidHeader = TRUE, width = 12, collapsible = TRUE,
        tabsetPanel(
          tabPanel("Select Existing Index",
            br(),
            p("Select the prefix of your existing HISAT2 index."),
            
            # --- Option 1: Prebuilt ---
            div(style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; border: 1px solid #e9ecef;",
              h5(icon("cloud-download-alt"), "Option 1: Prebuilt Indexes", 
                 info_icon("Download prebuilt indexes for standard species from Google Drive.")),
              p("Download and select standard indexes.", style="font-size: 0.9em; color: #666;"),
              uiOutput(ns("prebuilt_list_ui"))
            ),
            
            # --- Option 2: Manual ---
            h5(icon("folder"), "Option 2: Use Local Index",
               info_icon("Manually select a prefix for an existing HISAT2 index on your system.")),
            fluidRow(
              column(8, 
                textInput(ns("index_prefix"), label = tagList("Index Prefix", info_icon("Path prefix (e.g., /path/to/idx). HISAT2 looks for .ht2 files starting with this name.")), placeholder = "/path/to/index/genome_prefix", width = "100%")
              ),
              column(4, style="margin-top: 25px;",
                actionButton(ns("browse_idx"), "Select .ht2 File", icon = icon("folder-open"), class = "btn-info")
              )
            )
          ),
          tabPanel("Build New Index",
             br(),
             p("Create a splice-aware HISAT2 index from FASTA and GTF."),
             fluidRow(
               column(6,
                 tags$label(tagList("Reference Genome (FASTA)", info_icon("Genome FASTA file containing genomic sequences (e.g., genome.fa)."))),
                 div(style="display:flex",
                   textInput(ns("ref_fasta"), label = NULL, placeholder = "/path/to/genome.fa", width = "100%"),
                   div(style="margin-left:5px", actionButton(ns("browse_fasta"), "Browse", class = "btn-info"))
                 )
               ),
               column(6,
                 tags$label(tagList("GTF Annotation", info_icon("GTF annotation file containing gene structures. Used to build a splice-aware index."))),
                 div(style="display:flex",
                   textInput(ns("ref_gtf"), label = NULL, placeholder = "/path/to/annotation.gtf", width = "100%"),
                   div(style="margin-left:5px", actionButton(ns("browse_gtf"), "Browse", class = "btn-info"))
                 )
               )
             ),
             fluidRow(
               column(4, textInput(ns("build_out_name"), label = tagList("Output Index Name", info_icon("Custom name for the generated HISAT2 index files.")), value = "genome_index")),
               column(4, numericInput(ns("build_threads"), label = tagList("Threads", info_icon("Number of CPU cores to use for indexing. Speeds up the build process.")), value = 4, min = 1)),
               column(4, style="margin-top:25px;", 
                      actionButton(ns("build_index_btn"), "Build Splice-Aware Index", class = "btn-primary", icon = icon("cogs")))
             ),
             verbatimTextOutput(ns("build_log"))
          )
        )
      )
    ),
    fluidRow(
      box(
        title = tagList("Alignment Execution", info_icon("Run HISAT2 alignment on FASTQ files. Supports both Single-End and Paired-End reads.")), 
        status = "primary", solidHeader = TRUE, width = 12,
        h4("Input Reads (FASTQ)", info_icon("Sequencing reads in FASTQ format (.fq, .fastq, .gz). Can process multiple samples in batch.")),
        p("Enter comma-separated paths for multiple samples. Order must match between R1 and R2."),
        fluidRow(
          column(6,
            tags$label(tagList("Read 1 (or Single End)", info_icon("Primary FASTQ file(s). For single-end reads, only this field is required."))),
            div(style="display:flex; align-items: flex-start;",
              textAreaInput(ns("r1_files"), label = NULL, placeholder = "/path/to/sample1_1.fq.gz, ...", height = "100px", width="100%"),
              div(style="margin-left:5px", actionButton(ns("browse_r1"), "Browse", icon = icon("folder-open"), class = "btn-info"))
            )
          ),
          column(6,
            tags$label(tagList("Read 2 (Optional)", info_icon("Companion FASTQ file(s) for paired-end sequencing. Leave empty for single-end."))),
            div(style="display:flex; align-items: flex-start;",
              textAreaInput(ns("r2_files"), label = NULL, placeholder = "/path/to/sample1_2.fq.gz, ...", height = "100px", width="100%"),
              div(style="margin-left:5px", actionButton(ns("browse_r2"), "Browse", icon = icon("folder-open"), class = "btn-info"))
            )
          )
        ),
        fluidRow(
          column(4, numericInput(ns("threads"), label = tagList("Threads", info_icon("Number of CPU cores to use for alignment. More threads reduce runtime.")), value = 4, min = 1, max = 128)),
          column(4, textInput(ns("out_subdir"), label = tagList("Output Sub-directory", info_icon("Folder name inside 'artifacts' where resulting BAM files will be stored.")), value = "hisat2_bam"))
        ),
        hr(),
        h4("3. Execution"),
        verbatimTextOutput(ns("cmd_preview")),
        br(),
        actionButton(ns("run_btn"), "Run Alignment", class = "btn-success btn-lg", icon = icon("play")),
        span(textOutput(ns("status_msg")), style = "margin-left: 10px; font-weight: bold;")
      )
    ),
    fluidRow(
      box(
        title = tagList("Terminal", info_icon("Real-time output from the HISAT2 command execution and progress logs.")), 
        status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE,
        tags$script("
          // Auto-scroll log to bottom
          $(document).on('shiny:value', function(event) {
            if (event.name.endsWith('live_log')) {
              var el = document.getElementById(event.name);
              if (el) {
                 // The output is inside a <pre> tag usually, or the div wrapper
                 // We need to scroll the container div
                 var container = $(el).closest('.log-container');
                 if (container.length) {
                   container.scrollTop(container[0].scrollHeight);
                 }
              }
            }
          });
        "),
        div(class = "log-container",
          style = "background-color: #000; color: #0f0; padding: 10px; font-family: monospace; height: 300px; overflow-y: scroll;",
          verbatimTextOutput(ns("live_log"))
        )
      )
    )
  )
}

mod_hisat2_server <- function(id, manager) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      process = NULL,
      log_file = NULL,
      is_running = FALSE,
      log_content = "",
      build_log_content = "",
      build_process_steps = "",
      current_sample = NULL,
      output_dir = NULL
    )

    output$build_log <- renderText({ rv$build_log_content })
    
    # --- Initialization / Restoration ---
    observeEvent(manager$events$workspace_loaded, {
      req(manager$events$workspace_loaded > 0)
      
      # 1. Restore Parameters (Base)
      pars <- manager$get_params(id)
      
      # Initialize determined prefix
      final_prefix <- ""
      
      # --- PRIORITY 1: SAVED INPUT (User Preference) ---
      if (!is.null(pars) && !is.null(pars$index_prefix) && pars$index_prefix != "") {
        final_prefix <- pars$index_prefix
        message("[HISAT2] Priority 1 (Saved): Using saved input: ", final_prefix)
      }
      
      # --- PRIORITY 2: AUTO-DETECTION (Physical File Check) ---
      # Only if no saved input
      if (final_prefix == "" && !is.null(manager$workspace_dir)) {
        idx_dir <- file.path(manager$workspace_dir, "artifacts", "hisat2_index")
        if (dir.exists(idx_dir)) {
          candidates <- list.files(idx_dir, pattern = "\\.1\\.ht2l?$", full.names = TRUE)
          if (length(candidates) > 0) {
            # Pick the newest file
            best_cand <- candidates[1]
            if (length(candidates) > 1) {
               info <- file.info(candidates)
               best_cand <- rownames(info)[which.max(info$mtime)]
            }
            final_prefix <- sub("\\.1\\.ht2l?$", "", best_cand)
            message("[HISAT2] Priority 2 (Auto-Detect): Found index at ", final_prefix)
          }
        }
      }

      if (!is.null(pars)) {
        # Restore logs to state (always)
        if (!is.null(pars$build_log)) rv$build_log_content <- pars$build_log
        if (!is.null(pars$build_process_steps)) rv$build_process_steps <- pars$build_process_steps
        
        # --- PRIORITY 3: LOG EXTRACTION ---
        if (final_prefix == "" && !is.null(rv$build_log_content)) {
             match <- regexpr("Done\\. Index set to:\\s+(.+?)\\s*$", rv$build_log_content, perl=TRUE)
             if (match[1] != -1) {
               start <- attr(match, "capture.start")[1]
               len <- attr(match, "capture.length")[1]
               extracted_prefix <- substr(rv$build_log_content, start, start + len - 1)
               extracted_prefix <- trimws(extracted_prefix)
               
               if (nzchar(extracted_prefix)) {
                 final_prefix <- extracted_prefix
                 message("[HISAT2] Priority 3 (Log): Extracted index prefix: ", final_prefix)
               }
             }
        }
        
        # Apply Updates
        if (final_prefix != "") {
           message("[HISAT2] Updating index_prefix UI with: ", final_prefix)
           updateTextInput(session, "index_prefix", value = final_prefix)
           showNotification(paste("Index set to:", basename(final_prefix)), type = "message")
        }
        
        # Restore other inputs
        updateTextInput(session, "ref_fasta", value = pars$ref_fasta)
        updateTextInput(session, "ref_gtf", value = pars$ref_gtf)
        updateTextAreaInput(session, "r1_files", value = pars$r1_files)
        updateTextAreaInput(session, "r2_files", value = pars$r2_files)
        updateNumericInput(session, "threads", value = pars$threads)
        updateTextInput(session, "out_subdir", value = pars$out_subdir)
      }
    })
    
    # --- Save State ---
    observe({
      # Save params whenever they change
      manager$state$params[[id]] <- list(
        index_prefix = input$index_prefix,
        ref_fasta = input$ref_fasta,
        ref_gtf = input$ref_gtf,
        r1_files = input$r1_files,
        r2_files = input$r2_files,
        threads = input$threads,
        out_subdir = input$out_subdir,
        build_log = rv$build_log_content,
        build_process_steps = rv$build_process_steps
      )
    })
    
    # --- Helpers ---
    get_native_file <- function(title = "Select File", filter_pattern = NULL, multi = FALSE) {
      if (nzchar(Sys.which("zenity"))) {
        tryCatch({
          filter_arg <- ""
          if (!is.null(filter_pattern)) filter_arg <- paste0('--file-filter="', filter_pattern, '"')
          multi_arg <- if(multi) "--multiple --separator=," else ""
          cmd <- paste0('zenity --file-selection --title="', title, '" ', filter_arg, ' ', multi_arg)
          path <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
          if (length(path) > 0 && nzchar(path)) return(path)
        }, error = function(e) NA)
      } else if (capabilities("tcltk") && !is.null(Sys.getenv("DISPLAY")) && Sys.getenv("DISPLAY") != "") {
         return(tcltk::tk_choose.files(caption = title, multi = multi))
      }
      return(NA)
    }
    
    # --- Handlers ---
    observe({
      # Auto-detect threads and set max
      total_threads <- parallel::detectCores(logical = TRUE)
      if (is.na(total_threads)) total_threads <- 4 # Fallback
      
      updateNumericInput(session, "threads", max = total_threads)
      
      # Warning logic
      req(input$threads)
      if (input$threads > max(1, total_threads - 4)) {
        showNotification(
          paste("Warning: Using", input$threads, "threads leaves few resources for the system (Total:", total_threads, "). System may become unresponsive."),
          type = "warning", duration = 5, id = "thread_warning"
        )
      } else {
        removeNotification("thread_warning")
      }
    })
    
    observe({
      # Logic for Build Threads
      total_threads <- parallel::detectCores(logical = TRUE)
      if (is.na(total_threads)) total_threads <- 4
      
      updateNumericInput(session, "build_threads", max = total_threads)
      
      req(input$build_threads)
      if (input$build_threads > max(1, total_threads - 4)) {
        showNotification(
          paste("Warning: Using", input$build_threads, "threads for Indexing may leave few resources."),
          type = "warning", duration = 5, id = "build_thread_warn"
        )
      } else {
        removeNotification("build_thread_warn")
      }
    })

    # --- Prebuilt Index Manager ---
    rv$refresh_prebuilt <- 0
    
    prebuilt_genomes <- list(
      "arabidopsis_thaliana" = list(label = "Arabidopsis thaliana TAIR10", url = "https://drive.google.com/drive/folders/16WiBNHNV6SR_gn77xogqQMHtP4iUVoVH?usp=drive_link"),
      "oryza_sativa_indica" = list(label = "Oryza sativa indica ASM465v1", url = "https://drive.google.com/drive/folders/1GpEIfGtWVoKCT2XHYrK4bW1og9PbxWhr?usp=sharing"),
      "oryza_sativa_japonica" = list(label = "Oryza sativa japonica IRGSP-1.0", url = "https://drive.google.com/drive/folders/1k44eK1a0BCqqCtT1kyVH-aeM3ePUhJWd?usp=sharing"),
      "arabidopsis_halleri" = list(label = "Arabidopsis halleri Ahal2.2", url = "https://drive.google.com/drive/folders/141ypdew6z4TdTdo17NpQWfCbnNN9wzjg?usp=drive_link"),
      "malus_domestica " = list(label = "Malus domestica GDDH13", url = "https://drive.google.com/drive/folders/1YPUtAGFDQ3z9I0dnC7bzkhDKgMIh0jNA?usp=sharing"),
      "slycopersicum" = list(label = "Slycopersicum ITAG5.0", url = "https://drive.google.com/drive/folders/15PxrFMqdZONzGxvvq1cPhlJkM_Alg9kW?usp=sharing"),
      "vitis_vinifera" = list(label = "Vitis vinifera PN40024_5.1_T2T", url = "https://drive.google.com/drive/folders/1RKdawzmE8YYRrHHyw-2IWmjr5ruRphl5?usp=sharing"),
      "zea_mays" = list(label = "Zea mays Zm-B73-5.0", url = "https://drive.google.com/drive/folders/1UU0Ng4PcZVsRlEbH0qaPeis0uryjqk6i?usp=sharing"),
      "solanum_tuberosum" = list(label = "Solanum tuberosum v6.1", url = "https://drive.google.com/drive/folders/15yINHZkCulNa5wNcbL2JpnWWTsxJpoNk?usp=sharing"),
      "mangifera_indica" = list(label = "Mangifera indica CATAS_Mindica_2.1", url = "https://drive.google.com/drive/folders/15yINHZkCulNa5wNcbL2JpnWWTsxJpoNk?usp=sharing")
    )
    
    # Helper: Check status
    check_index_status <- function(code) {
       if (code == "") return(NULL)
       store_dir <- file.path(getwd(), "reference_indexes")
       target_dir <- file.path(store_dir, code)
       if (!dir.exists(target_dir)) return(NULL)
       
       # Find prefix
       cands <- list.files(target_dir, pattern = "\\.1\\.ht2l?$", full.names = TRUE, recursive = TRUE)
       if (length(cands) > 0) return(sub("\\.1\\.ht2l?$", "", cands[1]))
       # If folder exists but no index found, still return something so we can remove it? 
       # For now, stick to original logic but maybe we should allow removing empty folders too.
       # Let's assume if dir exists, it's downloaded.
       return(target_dir) 
    }
    
    # Helper: Get folder size in human readable format
    get_folder_size <- function(path) {
      if (!dir.exists(path)) return("0 B")
      files <- list.files(path, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
      if (length(files) == 0) return("0 B")
      size <- sum(file.info(files)$size, na.rm = TRUE)
      class(size) <- "object_size"
      format(size, units = "auto")
    }

    # Render Dropdown & Button
    output$prebuilt_list_ui <- renderUI({
      trigger <- rv$refresh_prebuilt
      
      # Select Input
      choices <- setNames(names(prebuilt_genomes), sapply(prebuilt_genomes, function(x) x$label))
      
      store_dir <- file.path(getwd(), "reference_indexes")
      total_size <- get_folder_size(store_dir)
      
      fluidRow(
        column(8,
           selectInput(ns("prebuilt_select"), label = tagList("Select Species", info_icon("Choose a plant species to download its pre-indexed HISAT2 genome.")), choices = c("Select..."="", choices), width="100%"),
           p(paste("Total Reference Storage:", total_size), style="font-size: 0.85em; color: #666; margin-top: -10px;")
        ),
        column(4, style="margin-top:25px;",
           uiOutput(ns("prebuilt_action_ui"))
        )
      )
    })
    
    # Dynamic Button based on selection
    output$prebuilt_action_ui <- renderUI({
       req(input$prebuilt_select)
       code <- input$prebuilt_select
       if (code == "") return(NULL)
       
       # Check if dir exists
       store_dir <- file.path(getwd(), "reference_indexes")
       target_dir <- file.path(store_dir, code)
       is_downloaded <- dir.exists(target_dir)
       
       if (is_downloaded) {
          # Attempt to find prefix to auto-populate
          cands <- list.files(target_dir, pattern = "\\.1\\.ht2l?$", full.names = TRUE, recursive = TRUE)
          if(length(cands) > 0) {
             prefix <- sub("\\.1\\.ht2l?$", "", cands[1])
             updateTextInput(session, "index_prefix", value = prefix)
          }
          
          actionButton(ns("do_remove"), "Remove Index", class="btn-danger", icon=icon("trash"), style="width:100%;")
       } else {
          actionButton(ns("do_download"), "Download Index", class="btn-primary", icon=icon("download"), style="width:100%;")
       }
    })
    
    # Remove Handler
    observeEvent(input$do_remove, {
        req(input$prebuilt_select)
        code <- input$prebuilt_select
        g <- prebuilt_genomes[[code]]
        
        store_dir <- file.path(getwd(), "reference_indexes")
        target_dir <- file.path(store_dir, code)
        
        showModal(modalDialog(
          title = "Confirm Removal",
          paste("Are you sure you want to delete the index for", g$label, "? This cannot be undone."),
          footer = tagList(
            modalButton("Cancel"),
            actionButton(ns("confirm_remove"), "Delete", class = "btn-danger")
          )
        ))
    })
    
    observeEvent(input$confirm_remove, {
        removeModal()
        req(input$prebuilt_select)
        code <- input$prebuilt_select
        
        store_dir <- file.path(getwd(), "reference_indexes")
        target_dir <- file.path(store_dir, code)
        
        if (dir.exists(target_dir)) {
           unlink(target_dir, recursive = TRUE)
           rv$refresh_prebuilt <- rv$refresh_prebuilt + 1
           showNotification("Index removed successfully.", type="message")
           updateTextInput(session, "index_prefix", value = "") # Clear input
        }
    })
    
    # Download Handler
    observeEvent(input$do_download, {
        req(input$prebuilt_select)
        code <- input$prebuilt_select
        g <- prebuilt_genomes[[code]]
        
        store_dir <- file.path(getwd(), "reference_indexes")
        if (!dir.exists(store_dir)) dir.create(store_dir, recursive = TRUE)
        
        target_dir <- file.path(store_dir, code)
        gdown_exe <- file.path(getwd(), "env", "bin", "gdown")
        
        withProgress(message = paste("Downloading", g$label, "index..."), value = 0, {
          tryCatch({
             if (!file.exists(gdown_exe)) stop("gdown executable not found in env/bin/. Please install it via conda.")
             
             incProgress(0.1, detail = "Starting gdown process...")
             
             # Use gdown to download folder
             # --folder flag is required for folder links
             args <- c(g$url, "-O", target_dir, "--folder")
             
             # Run gdown
             # capturing output to log if needed
             exit_code <- system2(gdown_exe, args = args, stdout = FALSE, stderr = FALSE)
             
             if (exit_code != 0) {
                 stop("gdown process exited with error code.")
             }
             
             if (!dir.exists(target_dir)) {
                 stop("Target directory was not created by gdown.")
             }
             
             incProgress(0.9, detail = "Locating index files...")
             
             rv$refresh_prebuilt <- rv$refresh_prebuilt + 1
             showNotification(paste(g$label, "Download Complete"), type="message")
             
             # The UI will auto-refresh and trigger the "Downloaded" state + auto-populate path
             
          }, error = function(e) {
             showNotification(paste("Download failed:", e$message), type="error")
          })
        })
    })
    
    observeEvent(input$browse_idx, {
      # Select one .ht2 file, then strip suffix
      path <- get_native_file("Select HISAT2 Index File", "*.ht2 *.ht2l", multi=FALSE)
      if (!is.na(path) && length(path) > 0) {
        clean <- sub("\\.[0-9]+\\.ht2l?$", "", path)
        clean <- sub("\\.ht2l?$", "", clean)
        updateTextInput(session, "index_prefix", value = clean)
      }
    })
    
    observeEvent(input$browse_fasta, {
      path <- get_native_file("Select Reference FASTA", "*.fa *.fasta *.fna", multi=FALSE)
      if (!is.na(path) && length(path) > 0) updateTextInput(session, "ref_fasta", value = path)
    })
    
    observeEvent(input$browse_gtf, {
      path <- get_native_file("Select GTF Annotation", "*.gtf", multi=FALSE)
      if (!is.na(path) && length(path) > 0) updateTextInput(session, "ref_gtf", value = path)
    })
    
    observeEvent(input$browse_r1, {
      paths <- get_native_file("Select Read 1 Files", "*.fq *.fastq *.gz", multi=TRUE)
      if (!is.na(paths) && length(paths) > 0) {
        # Zenity returns comma sep string if multiple, tcltk returns vector
        # Normalize to string
        val <- if(length(paths) > 1) paste(paths, collapse=", ") else paths
        current <- input$r1_files
        new_val <- if(current == "") val else paste(current, val, sep=", ")
        updateTextAreaInput(session, "r1_files", value = new_val)
      }
    })
    
    observeEvent(input$browse_r2, {
      paths <- get_native_file("Select Read 2 Files", "*.fq *.fastq *.gz", multi=TRUE)
      if (!is.na(paths) && length(paths) > 0) {
        val <- if(length(paths) > 1) paste(paths, collapse=", ") else paths
        current <- input$r2_files
        new_val <- if(current == "") val else paste(current, val, sep=", ")
        updateTextAreaInput(session, "r2_files", value = new_val)
      }
    })
    
    # --- Index Building Logic ---
    observeEvent(input$build_index_btn, {
      req(input$ref_fasta, input$ref_gtf, input$build_out_name)
      
      ws_dir <- if (is.null(manager$workspace_dir)) getwd() else manager$workspace_dir
      idx_dir <- file.path(ws_dir, "artifacts", "hisat2_index")
      if (!dir.exists(idx_dir)) dir.create(idx_dir, recursive = TRUE)
      
      prefix <- file.path(idx_dir, input$build_out_name)
      ss_file <- file.path(idx_dir, "splice_sites.txt")
      exon_file <- file.path(idx_dir, "exons.txt")
      log_file <- file.path(idx_dir, "build.log")
      
      # Environment Setup
      env_bin <- file.path(getwd(), "env", "bin")
      python_exe <- file.path(env_bin, "python3")
      script_ss <- file.path(env_bin, "hisat2_extract_splice_sites.py")
      script_ex <- file.path(env_bin, "hisat2_extract_exons.py")
      
      rv$build_process_steps <- "Starting Index Build...\n"
      rv$build_log_content <- rv$build_process_steps
      
      # 1. Extract Splice Sites & Exons
      withProgress(message = "Preparing Annotation Features...", value = 0.2, {
        tryCatch({
          # Splice Sites (Use Env Python)
          cmd_ss <- paste(python_exe, script_ss, paste0("\"", input$ref_gtf, "\""), ">", paste0("\"", ss_file, "\""))
          if(system(cmd_ss) != 0) stop("Splice site extraction failed.")
          rv$build_process_steps <- paste0(rv$build_process_steps, "Extracted splice sites.\n")
          rv$build_log_content <- rv$build_process_steps
          
          # Exons (Use Env Python)
          cmd_ex <- paste(python_exe, script_ex, paste0("\"", input$ref_gtf, "\""), ">", paste0("\"", exon_file, "\""))
          if(system(cmd_ex) != 0) stop("Exon extraction failed.")
          rv$build_process_steps <- paste0(rv$build_process_steps, "Extracted exons.\n")
          rv$build_log_content <- rv$build_process_steps
          
        }, error = function(e) {
          rv$build_process_steps <- paste0(rv$build_process_steps, "Error extracting features: ", e$message, "\n")
          rv$build_log_content <- rv$build_process_steps
        })
      })
      
      # 2. Build Index Script
      script_file <- file.path(idx_dir, "run_build.sh")
      
      cmd_build <- paste(
        "hisat2-build",
        "-p", input$build_threads,
        "--ss", paste0("\"", ss_file, "\""),
        "--exon", paste0("\"", exon_file, "\""),
        paste0("\"", input$ref_fasta, "\""),
        paste0("\"", prefix, "\"")
      )
      
      # Write script with environment activation (PATH)
      script_content <- c(
        "#!/bin/bash",
        "set -e",
        paste0("export PATH=\"", env_bin, ":$PATH\""),
        "",
        "echo \"Starting hisat2-build...\"",
        cmd_build
      )
      
      writeLines(script_content, script_file)
      Sys.chmod(script_file, "755")
      
      rv$build_process_steps <- paste0(rv$build_process_steps, "Running hisat2-build via script: ", script_file, "\n")
      rv$build_log_content <- rv$build_process_steps
      
      p <- processx::process$new("bash", args = c(script_file), stdout = log_file, stderr = "2>&1")
      
      # Polling for Build
      observe({
        if (p$is_alive()) {
           invalidateLater(1000, session)
           # Update log
           if (file.exists(log_file)) {
             lines <- readLines(log_file, warn=FALSE)
             # Show last 10 lines
             tail_lines <- if(length(lines)>10) lines[(length(lines)-10):length(lines)] else lines
             rv$build_log_content <- paste(rv$build_process_steps, paste(tail_lines, collapse="\n"), sep="\n")
           }
        } else {
           # Finished
           exit_code <- p$get_exit_status()
           if (exit_code == 0) {
             showNotification("Index Built Successfully!", type = "message")
             updateTextInput(session, "index_prefix", value = prefix)
             rv$build_log_content <- paste0(rv$build_process_steps, "\nDone. Index set to: ", prefix)
           } else {
             showNotification("Index Build Failed.", type = "error")
             rv$build_log_content <- paste0(rv$build_process_steps, "\nFAILED. See artifacts/hisat2_index/build.log for details.")
           }
        }
      })
    })
    
    generated_cmd <- reactive({
      prefix <- if (!is.null(input$index_prefix) && nzchar(input$index_prefix)) input$index_prefix else "[INDEX_PREFIX]"
      threads <- if (!is.null(input$threads)) input$threads else 4
      
      sprintf(
        "hisat2 -x %s -1 sample_1.fq -2 sample_2.fq -p %s | samtools view -bS - | samtools sort -o sample.sorted.bam",
        prefix, threads
      )
    })
    
    output$cmd_preview <- renderText({ generated_cmd() })
    
    # --- Execution Logic ---
    observeEvent(input$run_btn, {
      # 0. Validation & UI Feedback
      valid <- TRUE
      msg <- c()
      
      # Clear previous errors
      shinyjs::runjs(sprintf("$('#%s').closest('.form-group').removeClass('has-error');", ns("index_prefix")))
      shinyjs::runjs(sprintf("$('#%s').closest('.form-group').removeClass('has-error');", ns("r1_files")))
      shinyjs::runjs(sprintf("$('#%s').closest('.form-group').removeClass('has-error');", ns("r2_files")))
      
      if (is.null(input$index_prefix) || input$index_prefix == "") {
        shinyjs::runjs(sprintf("$('#%s').closest('.form-group').addClass('has-error');", ns("index_prefix")))
        valid <- FALSE
        msg <- c(msg, "Index Prefix is missing.")
      }
      
      if (is.null(input$r1_files) || input$r1_files == "") {
        shinyjs::runjs(sprintf("$('#%s').closest('.form-group').addClass('has-error');", ns("r1_files")))
        valid <- FALSE
        msg <- c(msg, "Read 1 files are missing.")
      }
      
      if (!valid) {
        showNotification(paste("Caution: Please fix the following issues:", paste(msg, collapse=" ")), type = "error")
        return()
      }
      
      # UI Feedback: Running
      shinyjs::disable("run_btn")
      updateActionButton(session, "run_btn", label = "Running...", icon = icon("spinner", class = "fa-spin"))
      
      # 1. Parse Inputs (Handle comma OR newline)
      r1s <- unlist(strsplit(input$r1_files, "[,\n]"))
      r1s <- trimws(r1s)
      r1s <- r1s[nzchar(r1s)]
      
      r2s <- unlist(strsplit(input$r2_files, "[,\n]"))
      r2s <- trimws(r2s)
      r2s <- r2s[nzchar(r2s)]
      
      # Validate Paired-End
      is_paired <- (length(r2s) > 0)
      if (is_paired && length(r1s) != length(r2s)) {
        shinyjs::runjs(sprintf("$('#%s').closest('.form-group').addClass('has-error');", ns("r1_files")))
        shinyjs::runjs(sprintf("$('#%s').closest('.form-group').addClass('has-error');", ns("r2_files")))
        
        showNotification("Caution: Number of R1 and R2 files must match for Paired-End analysis.", type = "error")
        
        # Reset Button on Error
        shinyjs::enable("run_btn")
        updateActionButton(session, "run_btn", label = "Run Alignment", icon = icon("play"))
        return()
      }
      
      # 2. Setup Output
      ws_dir <- if (is.null(manager$workspace_dir)) getwd() else manager$workspace_dir
      out_dir <- file.path(ws_dir, "artifacts", input$out_subdir)
      log_dir <- file.path(out_dir, "hisat2_logs")
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
      
      rv$output_dir <- out_dir
      rv$log_file <- file.path(out_dir, "run_alignment.log")
      file.create(rv$log_file)
      
      # 3. Generate Master Shell Script
      script_file <- file.path(out_dir, "run_alignment.sh")
      
      # Resolve App Environment Path
      app_root <- getwd()
      env_bin_path <- file.path(app_root, "env", "bin")
      
      # Header & Functions
      script_content <- c(
        "#!/usr/bin/env bash",
        "",
        "# Generated by PhytoSplice",
        "# Based on user-provided 'align.sh' and 'cutadapt_trim.sh'",
        "",
        "set -euo pipefail",
        "",
        "# Source Local Environment",
        paste0("export PATH=\"", env_bin_path, ":$PATH\""),
        "",
        paste0("HISAT2_INDEX=\"", input$index_prefix, "\""),
        paste0("HISAT2_THREADS=", input$threads),
        "SAMTOOLS_THREADS=4",  # hardcoded or derived? Using 4 to be safe.
        paste0("OUTDIR=\"", out_dir, "\""),
        paste0("LOGDIR=\"", log_dir, "\""),
        paste0("DO_TRIM=\"false\""),
        "",
        "# Check executables",
        "command -v hisat2 >/dev/null 2>&1 || { echo \"Error: hisat2 not found\"; exit 1; }",
        "command -v samtools >/dev/null 2>&1 || { echo \"Error: samtools not found\"; exit 1; }",
        "",
        "process_sample() {",
        "    local R1=$1",
        "    local R2=$2",
        "    local SAMPLE=$3",
        "    ",
        "    local LOG=\"$LOGDIR/${SAMPLE}.log\"",
        "    # Output sorted BAM directly",
        "    local OUTBAM=\"$OUTDIR/${SAMPLE}.sorted.bam\"",
        "    local NOVEL_SS=\"$OUTDIR/${SAMPLE}.novel_splicesites.txt\"",
        "    ",
        "    echo \"[$(date +'%Y-%m-%d %H:%M:%S')] Processing: $SAMPLE\"",
        "    echo \"  Log: $LOG\"",
        "    ",
        "    # --- Adapter Trimming Logic ---",
        "    if [ \"$DO_TRIM\" == \"true\" ]; then",
        "        local TRIM_DIR=\"$OUTDIR/trimmed\"",
        "        mkdir -p \"$TRIM_DIR\"",
        "        echo \"  Running Trim Galore...\"",
        "        command -v trim_galore >/dev/null 2>&1 || { echo \"Error: trim_galore not found\"; exit 1; }",
        "        ",
        "        if [ -z \"$R2\" ]; then",
        "            # Single End",
        "            trim_galore -j 4 \"$R1\" -o \"$TRIM_DIR\" >> \"$LOG\" 2>&1",
        "            # Detect Output (Trim Galore appends _trimmed.fq.gz for SE)",
        "            local base=$(basename \"$R1\" | sed 's/\\.fastq.*$//;s/\\.fq.*$//')",
        "            local trimmed=$(find \"$TRIM_DIR\" -name \"${base}*_trimmed.fq*\" | head -n 1)",
        "            if [ -n \"$trimmed\" ]; then R1=\"$trimmed\"; else echo \"Error: Trim output not found\"; exit 1; fi",
        "        else",
        "            # Paired End",
        "            trim_galore -j 4 --paired \"$R1\" \"$R2\" -o \"$TRIM_DIR\" >> \"$LOG\" 2>&1",
        "            # Detect Output (Trim Galore appends _val_1.fq.gz and _val_2.fq.gz for PE)",
        "            local base1=$(basename \"$R1\" | sed 's/\\.fastq.*$//;s/\\.fq.*$//')",
        "            local base2=$(basename \"$R2\" | sed 's/\\.fastq.*$//;s/\\.fq.*$//')",
        "            local t1=$(find \"$TRIM_DIR\" -name \"${base1}*_val_1.fq*\" | head -n 1)",
        "            local t2=$(find \"$TRIM_DIR\" -name \"${base2}*_val_2.fq*\" | head -n 1)",
        "            if [ -n \"$t1\" ] && [ -n \"$t2\" ]; then R1=\"$t1\"; R2=\"$t2\"; else echo \"Error: Trim output not found\"; exit 1; fi",
        "        fi",
        "        echo \"  Trimming complete. Using trimmed files for alignment.\"",
        "    fi",
        "    ",
        "    echo \"  Input R1: $R1\"",
        "    if [ -n \"$R2\" ]; then echo \"  Input R2: $R2\"; fi",
        "    ",
        "    # Log Command",
        "    {",
        "        echo \"----- Command start: $(date) -----\"",
        "        echo \"hisat2 --threads $HISAT2_THREADS -x $HISAT2_INDEX ... | samtools view -bS - | samtools sort -o $OUTBAM\"",
        "        echo",
        "    } >> \"$LOG\"",
        "    ",
        "    # Execution",
        "    if [ -z \"$R2\" ]; then",
        "        # Single End",
        "        hisat2 --threads \"$HISAT2_THREADS\" -x \"$HISAT2_INDEX\" -U \"$R1\" \\",
        "            --novel-splicesite-outfile \"$NOVEL_SS\" \\",
        "            2>>\"$LOG\" | samtools view -bS - | samtools sort -@ \"$SAMTOOLS_THREADS\" -o \"$OUTBAM\"",
        "    else",
        "        # Paired End",
        "        hisat2 --threads \"$HISAT2_THREADS\" -x \"$HISAT2_INDEX\" -1 \"$R1\" -2 \"$R2\" \\",
        "            --novel-splicesite-outfile \"$NOVEL_SS\" \\",
        "            2>>\"$LOG\" | samtools view -bS - | samtools sort -@ \"$SAMTOOLS_THREADS\" -o \"$OUTBAM\"",
        "    fi",
        "    ",
        "    # Index",
        "    samtools index \"$OUTBAM\"",
        "    ",
        "    if [[ $? -eq 0 ]]; then",
        "        echo \"[$(date +'%Y-%m-%d %H:%M:%S')] Finished: $SAMPLE\" | tee -a \"$LOG\"",
        "    else",
        "        echo \"[$(date +'%Y-%m-%d %H:%M:%S')] ERROR during alignment: $SAMPLE\" | tee -a \"$LOG\"",
        "        exit 1",
        "    fi",
        "    echo \"\"",
        "}",
        "",
        "echo \"[$(date +'%Y-%m-%d %H:%M:%S')] Starting HISAT2 batch alignment\""
      )
      
      # Add calls for each sample
      for (i in seq_along(r1s)) {
        sample_name <- tools::file_path_sans_ext(basename(r1s[i]))
        sample_name <- sub("(_R1|_1)$", "", sample_name)
        
        r1_path <- r1s[i]
        r2_path <- if(is_paired) r2s[i] else ""
        
        call_line <- paste0("process_sample \"", r1_path, "\" \"", r2_path, "\" \"", sample_name, "\"")
        script_content <- c(script_content, call_line)
      }
      
      script_content <- c(script_content, "echo \"[$(date +'%Y-%m-%d %H:%M:%S')] All done.\"")
      
      writeLines(script_content, script_file)
      Sys.chmod(script_file, "755")
      
      rv$log_content <- paste0("Generated Alignment Script: ", script_file, "\nRunning...\n")
      
      # 4. Run the Script
      rv$is_running <- TRUE
      rv$process <- processx::process$new("bash", args = c(script_file), stdout = rv$log_file, stderr = "2>&1")
    })
    
    # --- Polling (Simplified) ---
    observe({
      req(rv$is_running, rv$process)
      invalidateLater(500, session)
      
      if(file.exists(rv$log_file)) {
        lines <- tryCatch(readLines(rv$log_file, warn=FALSE), error=function(e) character(0))
        rv$log_content <- paste(lines, collapse="\n")
      }
      
      if (!rv$process$is_alive()) {
         rv$is_running <- FALSE
         
         # Reset Button Logic
         shinyjs::enable("run_btn")
         updateActionButton(session, "run_btn", label = "Run Alignment", icon = icon("play"))
         
         exit_code <- rv$process$get_exit_status()
         
         if (exit_code == 0) {
           msg <- "\nAlignment Batch Completed Successfully.\n"
           cat(msg, file = rv$log_file, append = TRUE)
           showNotification("Alignment Finished!", type = "message")
           
           # Register results with manager?
           # manager$set_result("hisat2", list(output_dir = rv$output_dir))
         } else {
           msg <- "\nAlignment Script Failed. See log for details.\n"
           cat(msg, file = rv$log_file, append = TRUE)
           showNotification("Alignment Failed.", type = "error")
         }
      }
    })
    
    
    output$live_log <- renderText({ rv$log_content })
    
  })
}
