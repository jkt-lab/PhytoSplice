# modules/mod_salmon.R
source("R/salmon_helpers.R")

mod_salmon_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(id = ns("tabs"),
      tabPanel("Index Generation",
        fluidRow(
          box(
            title = tagList("Salmon Index Builder (Decoy-Aware)", info_icon("What it does: Creates a transcriptome index for Salmon quantification.\nHow it works: Uses the provided FASTA and GTF to build a decoy-aware index, which improves quantification accuracy by preventing reads that originate from unannotated genomic regions from falsely mapping to transcripts.")), status = "primary", solidHeader = TRUE, width = 12,
            p("This module implements the recommended 'Decoy-Aware' indexing strategy for Salmon. It requires a Genome FASTA and an Annotation GTF."),
            
            # Input Files
            fluidRow(
              column(6,
                div(class="form-group",
                  tags$label(tagList("Genome FASTA (.fa/.fasta)", info_icon("What it does: The reference genome sequence.\nHow it works: Used to extract transcript sequences and generate the decoy sequence set."))),
                  div(style="display:flex",
                    textInput(ns("genome_fa"), label = NULL, placeholder = "/path/to/genome.fa", width = "100%"),
                    div(style="margin-left:5px", actionButton(ns("browse_genome"), "Browse", icon = icon("folder-open"), class = "btn-info"))
                  )
                )
              ),
              column(6,
                div(class="form-group",
                  tags$label(tagList("Annotation GTF (.gtf)", info_icon("What it does: The gene annotation file.\nHow it works: Defines the boundaries of transcripts to be extracted from the genome."))),
                  div(style="display:flex",
                    textInput(ns("gtf_file"), label = NULL, placeholder = "/path/to/annotation.gtf", width = "100%"),
                    div(style="margin-left:5px", actionButton(ns("browse_gtf"), "Browse", icon = icon("folder-open"), class = "btn-info"))
                  )
                )
              )
            ),
            
            # Parameters
            fluidRow(
              column(6,
                numericInput(ns("nthreads"), label = tagList("Threads", info_icon("What it does: Number of CPU cores to use.\nHow it works: Speeds up index building. Use all available cores if possible.")), value = 4, min = 1, max = 32)
              ),
              column(6,
                textInput(ns("index_name"), label = tagList("Index Output Name", info_icon("What it does: The folder name for the generated index within the 'artifacts/salmon/' directory.")), value = "salmon_index")
              )
            ),
            
            hr(),
            actionButton(ns("run_index"), "Build Salmon Index", class = "btn-danger btn-lg", icon = icon("cogs")),
            
            hr(),
            h4(tagList("Execution Log", info_icon("What it does: Shows real-time output from the Salmon index build process."))),
            verbatimTextOutput(ns("build_log"))
          )
        )
      ),
      tabPanel("Quantification",
        box(
          title = tagList("Quantification", info_icon("What it does: Quantifies transcript abundances from RNA-seq reads.\nHow it works: Maps reads against the Salmon index to estimate TPM (Transcripts Per Million). Output is automatically merged into a single matrix for SUPPA.")), status = "info", solidHeader = TRUE, width = 12,
          p("Run transcript quantification using Salmon. Output will be automatically formatted for SUPPA."),
          fluidRow(
            column(6,
              div(class="form-group",
                tags$label(tagList("Salmon Index Directory", info_icon("What it does: The path to the index generated in the 'Index Generation' tab."))),
                div(style="display:flex",
                  textInput(ns("salmon_index"), label = NULL, placeholder = "/path/to/salmon_index", width = "100%"),
                  div(style="margin-left:5px", actionButton(ns("browse_index"), "Browse", class = "btn-info"))
                )
              )
            ),
            column(6, numericInput(ns("threads_quant"), label = tagList("Threads", info_icon("What it does: Number of CPU cores to allocate per sample quantification.")), value = 4, min = 1, max = 64))
          ),
          fluidRow(
            column(6,
              div(class="form-group",
                tags$label(tagList("Read 1 Files (Comma Separated)", info_icon("What it does: FASTQ files containing the first reads (or single-end reads).\nHow it works: Supports batch processing. Separate multiple files with commas."))),
                div(style="display:flex",
                  textAreaInput(ns("r1_files"), label = NULL, placeholder = "/path/to/s1_R1.fq.gz, /path/to/s2_R1.fq.gz", height = "100px", width = "100%"),
                  div(style="margin-left:5px", actionButton(ns("browse_r1"), "Browse", icon = icon("folder-open"), class = "btn-info"))
                )
              )
            ),
            column(6,
              div(class="form-group",
                tags$label(tagList("Read 2 Files (Optional)", info_icon("What it does: FASTQ files containing the second reads for paired-end data.\nHow it works: Must match the order of Read 1 files. Leave empty for single-end."))),
                div(style="display:flex",
                  textAreaInput(ns("r2_files"), label = NULL, placeholder = "/path/to/s1_R2.fq.gz, /path/to/s2_R2.fq.gz", height = "100px", width = "100%"),
                  div(style="margin-left:5px", actionButton(ns("browse_r2"), "Browse", icon = icon("folder-open"), class = "btn-info"))
                )
              )
            )
          ),
          fluidRow(
            column(12,
              h5(tagList("Execution Command", info_icon("What it does: Preview of the bash command that will be executed."))),
              verbatimTextOutput(ns("salmon_cmd_preview")),
              br(),
              uiOutput(ns("run_quant_ui")),
              span(textOutput(ns("salmon_status")), style = "margin-left: 10px; font-weight: bold;"),
              hr(),
              h5(tagList("Quantification Log", info_icon("What it does: Shows progress and messages from the Salmon quantification run."))),
              verbatimTextOutput(ns("quant_log"))
            )
          )
        )
      )
    )
  )
}

mod_salmon_server <- function(id, manager) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      log_text = "",
      quant_log_text = "",
      is_running = FALSE,
      process = NULL
    )
    
    # Helper for logging
    append_log <- function(msg) {
      rv$log_text <- paste0(rv$log_text, msg, "\n")
    }
    
    append_quant_log <- function(msg) {
      rv$quant_log_text <- paste0(rv$quant_log_text, msg, "\n")
    }
    
    # --- Initialization / Restoration ---
    observe({
      # Restore Results (Logs)
      res <- manager$get_result(id)
      if (!is.null(res)) {
        if (!is.null(res$quant_log)) rv$quant_log_text <- res$quant_log
        if (!is.null(res$index_log)) rv$log_text <- res$index_log
      }
      
      # Restore Parameters (Inputs)
      pars <- manager$get_params(id)
      if (!is.null(pars)) {
        if (!is.null(pars$genome_fa)) updateTextInput(session, "genome_fa", value = pars$genome_fa)
        if (!is.null(pars$gtf_file)) updateTextInput(session, "gtf_file", value = pars$gtf_file)
        if (!is.null(pars$nthreads)) updateNumericInput(session, "nthreads", value = pars$nthreads)
        if (!is.null(pars$index_name)) updateTextInput(session, "index_name", value = pars$index_name)
        
        if (!is.null(pars$salmon_index)) updateTextInput(session, "salmon_index", value = pars$salmon_index)
        if (!is.null(pars$threads_quant)) updateNumericInput(session, "threads_quant", value = pars$threads_quant)
        if (!is.null(pars$r1_files)) updateTextAreaInput(session, "r1_files", value = pars$r1_files)
        if (!is.null(pars$r2_files)) updateTextAreaInput(session, "r2_files", value = pars$r2_files)
      }
    })
    
    # --- Persistence: Auto-save inputs to manager state ---
    observe({
      manager$state$params[[id]] <- list(
        genome_fa = input$genome_fa,
        gtf_file = input$gtf_file,
        nthreads = input$nthreads,
        index_name = input$index_name,
        salmon_index = input$salmon_index,
        threads_quant = input$threads_quant,
        r1_files = input$r1_files,
        r2_files = input$r2_files
      )
    })
    
    output$run_quant_ui <- renderUI({
      if (rv$is_running) {
        actionButton(ns("run_salmon_disabled"), "Running...", class = "btn-success btn-lg disabled", icon = icon("spinner", class="fa-spin"))
      } else {
        actionButton(ns("run_salmon"), "Run Salmon Quantification", class = "btn-success btn-lg", icon = icon("dna"))
      }
    })
    
    # --- Native Browser Logic ---
    get_native_file <- function(title = "Select File", filter_pattern = NULL, multi = FALSE, dir = FALSE) {
      if (nzchar(Sys.which("zenity"))) {
        tryCatch({
          if (dir) {
             cmd <- paste0('zenity --file-selection --directory --title="', title, '"')
          } else {
             filter_arg <- ""
             if (!is.null(filter_pattern)) filter_arg <- paste0('--file-filter="', filter_pattern, '"')
             multi_arg <- if(multi) "--multiple --separator=," else ""
             cmd <- paste0('zenity --file-selection --title="', title, '" ', filter_arg, ' ', multi_arg)
          }
          path <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
          if (length(path) > 0 && nzchar(path)) return(path)
        }, error = function(e) NA)
      } else if (capabilities("tcltk") && !is.null(Sys.getenv("DISPLAY")) && Sys.getenv("DISPLAY") != "") {
         if (dir) {
            return(tcltk::tk_choose.dir(caption = title))
         } else {
            return(tcltk::tk_choose.files(caption = title, multi = multi))
         }
      }
      return(NA)
    }

    # --- File Browsers (Native) ---
    observeEvent(input$browse_genome, {
      f <- get_native_file("Select Genome FASTA", filter_pattern = "*.fa *.fasta")
      if (!is.na(f)) updateTextInput(session, "genome_fa", value = f)
    })
    
    observeEvent(input$browse_gtf, {
      f <- get_native_file("Select Annotation GTF", filter_pattern = "*.gtf")
      if (!is.na(f)) updateTextInput(session, "gtf_file", value = f)
    })
    
    observeEvent(input$browse_index, {
      path <- get_native_file("Select Salmon Index Directory", dir=TRUE)
      if (!is.na(path)) updateTextInput(session, "salmon_index", value = path)
    })
    
    observeEvent(input$browse_r1, {
      paths <- get_native_file("Select Read 1 Files", "*.fq *.fastq *.gz", multi=TRUE)
      if (!is.na(paths) && length(paths) > 0) {
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
    
    # --- Build Index Logic ---
    observeEvent(input$run_index, {
      req(input$genome_fa, input$gtf_file)
      
      # 1. Check Manager (Workspace)
      if (is.null(manager$workspace_dir)) {
        showNotification("Please create/load a workspace first!", type = "error")
        return()
      }
      
      # 2. Setup Paths
      genome_path <- input$genome_fa
      gtf_path <- input$gtf_file
      
      # Output inside workspace artifacts
      out_dir <- file.path(manager$workspace_dir, "artifacts", "salmon", input$index_name)
      
      if (dir.exists(out_dir)) {
        showNotification("Index directory exists. Overwriting...", type = "warning")
        unlink(out_dir, recursive = TRUE)
      }
      
      # 3. Execution (Background via Future or simple blocking for now - Keeping blocking for simplicity as per MVP)
      # Ideally, use processx/callr in background, but for now we run directly
      
      rv$log_text <- paste0("Starting Salmon Index Build...\nTarget: ", out_dir, "\n")
      
      withProgress(message = 'Building Salmon Index...', value = 0, {
        
        # We wrap the log update in a callback
        log_callback <- function(msg) {
          append_log(msg)
          # Force UI update (hacky in blocking R)
          # In a blocking call, this won't update UI real-time unless we use processx/future.
          # For MVP, we'll see the log at the end or use a timer.
        }
        
        # Using a specialized observer for async would be better, but let's try direct call
        success <- tryCatch({
          build_salmon_index(
            genome_fa = genome_path,
            gtf_path = gtf_path,
            output_dir = out_dir,
            threads = input$nthreads,
            callbacks = list(log = log_callback)
          )
        }, error = function(e) {
          append_log(paste("Exception:", e$message))
          FALSE
        })
        
        if (success) {
          showNotification("Salmon Index Built Successfully!", type = "message")
          append_log("--- DONE ---")
        } else {
          showNotification("Salmon Index Build Failed.", type = "error")
          append_log("--- FAILED ---")
        }
        
      })
    })
    
    # --- Quantification Logic ---
    
    salmon_generated_cmd <- reactive({
      idx <- if(!is.null(input$salmon_index) && nzchar(input$salmon_index)) input$salmon_index else "[INDEX]"
      r1 <- if(!is.null(input$r1_files) && nzchar(input$r1_files)) "sample_R1.fq.gz" else "[R1]"
      r2 <- if(!is.null(input$r2_files) && nzchar(input$r2_files)) "sample_R2.fq.gz" else ""
      r2_flag <- if(nzchar(r2)) paste("-2", r2) else ""
      
      paste(
        "salmon quant -i", idx, "-l A", 
        "-1", r1, r2_flag,
        "-p", if(!is.null(input$threads_quant)) input$threads_quant else 4, 
        "--validateMappings -o output_dir/sample_name"
      )
    })
    
    output$salmon_cmd_preview <- renderText({ salmon_generated_cmd() })

    observeEvent(input$run_salmon, {
      req(input$salmon_index, input$r1_files)
      
      ws_dir <- if (is.null(manager$workspace_dir)) getwd() else manager$workspace_dir
      out_dir <- file.path(ws_dir, "artifacts", "salmon")
      log_dir <- file.path(out_dir, "logs")
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
      
      rv$quant_log_text <- paste0("Starting Salmon Quantification...\n")
      
      # Parse Files
      r1s <- trimws(unlist(strsplit(input$r1_files, "[,\n]")))
      r1s <- r1s[nzchar(r1s)]
      r2s <- trimws(unlist(strsplit(input$r2_files, "[,\n]")))
      r2s <- r2s[nzchar(r2s)]
      is_paired <- (length(r2s) > 0)
      
      # Generate Script
      script_file <- file.path(out_dir, "run_salmon.sh")
      env_bin <- file.path(getwd(), "env", "bin")
      threads <- if(!is.null(input$threads_quant)) input$threads_quant else 4
      
      script_lines <- c(
        "#!/bin/bash",
        "set -e",
        paste0("export PATH=\"", env_bin, ":$PATH\""),
        paste0("INDEX=\"", input$salmon_index, "\""),
        paste0("OUT_BASE=\"", out_dir, "\""),
        paste0("THREADS=", threads),
        "",
        "command -v salmon >/dev/null 2>&1 || { echo 'Error: salmon not found'; exit 1; }",
        "",
        "run_sample() {",
        "  R1=$1",
        "  R2=$2",
        "  NAME=$3",
        "  echo \"Processing $NAME...\"",
        "  mkdir -p \"$OUT_BASE/$NAME\"",
        "  if [ -z \"$R2\" ]; then",
        "     salmon quant -i \"$INDEX\" -l A -r \"$R1\" -p $THREADS --validateMappings -o \"$OUT_BASE/$NAME\"",
        "  else",
        "     salmon quant -i \"$INDEX\" -l A -1 \"$R1\" -2 \"$R2\" -p $THREADS --validateMappings -o \"$OUT_BASE/$NAME\"",
        "  fi",
        "}"
      )
      
      for (i in seq_along(r1s)) {
        sname <- tools::file_path_sans_ext(basename(r1s[i]))
        sname <- sub("(_R1|_1)$", "", sname)
        r2_arg <- if(is_paired) r2s[i] else ""
        script_lines <- c(script_lines, paste0("run_sample \"", r1s[i], "\" \"", r2_arg, "\" \"", sname, "\""))
      }
      
      writeLines(script_lines, script_file)
      Sys.chmod(script_file, "755")
      
      append_quant_log(paste0("Running Script: ", script_file))
      
      # Using processx for background execution
      # We need to store the process in rv to poll it
      rv$is_running <- TRUE
      rv$process <- processx::process$new("bash", args = c(script_file), stdout = file.path(log_dir, "latest.log"), stderr = "2>&1")
      
      # Polling observer
      observe({
        req(rv$process)
        invalidateLater(500, session)
        
        # Read latest log from file
        if (file.exists(file.path(log_dir, "latest.log"))) {
           lines <- readLines(file.path(log_dir, "latest.log"), warn=FALSE)
           # Simple way: just update the whole log content from file plus our header
           rv$quant_log_text <- paste("Starting Salmon Quantification...", paste(lines, collapse="\n"), sep="\n")
        }
        
        if (!rv$process$is_alive()) {
          rv$is_running <- FALSE
          
          if (rv$process$get_exit_status() == 0) {
            append_quant_log("\nSalmon finished. Merging output...")
            
            # Find quant.sf files
            quant_files <- list.files(out_dir, pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
            names(quant_files) <- basename(dirname(quant_files))
            
            merge_out <- file.path(out_dir, "all_tpm.txt")
            merged <- merge_salmon_quants(quant_files, merge_out)
            
            if (!is.null(merged)) {
               append_quant_log(paste0("Merge successful: ", merge_out))
               showNotification("Salmon Quantification & Merge Complete!", type="message")
            } else {
               append_quant_log("Merge failed.")
            }
          } else {
            showNotification("Salmon Quantification Failed", type="error")
            append_quant_log("FAILED.")
          }
          
          # Save state
          manager$set_result(id, list(
             quant_log = rv$quant_log_text,
             index_log = rv$log_text
          ))
          
          # Stop polling
          rv$process <- NULL 
        }
      })
    })
    
    output$quant_log <- renderText({
      rv$quant_log_text
    })
    
    output$build_log <- renderText({
      rv$log_text
    })
  })
}
