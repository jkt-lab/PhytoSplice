# Libraries lazy-loaded in app.R server or global.R
library(shinyjs)
source("modules/mod_event_based.R")

mod_rmats_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(id = ns("tabs"),
      # --- Tab 1: Execution ---
      tabPanel("Execution",
        fluidRow(
          box(
            title = tagList("rMATS Execution", info_icon("Run multivariate analysis of transcript splicing (rMATS) on BAM files. Supports differential splicing detection between two groups.")), 
            status = "primary", solidHeader = TRUE, width = 12,
            h4("1. Input BAM Files", info_icon("Sorted and indexed BAM files from HISAT2 or other aligners. Provide at least one file for each group.")),
            p("Select BAM files using the system dialog (Browse) or select from Workspace."),
            fluidRow(
              column(6, 
                div(class="form-group",
                  tags$label(tagList("Case (b1) BAMs", info_icon("BAM files representing the experimental or 'case' group. Comma-separated paths."))),
                  div(style="display:flex; align-items:flex-start;",
                     textAreaInput(ns("b1"), label = NULL, placeholder = "/path/to/case1.bam, ...", width = "100%", height = "100px"),
                     div(style="margin-left:5px; display:flex; flex-direction:column;", 
                         actionButton(ns("browse_b1"), "Browse", icon = icon("folder-open"), class = "btn-info", style="margin-bottom:5px;"),
                         actionButton(ns("ws_b1"), "Worksp", icon = icon("box-open"), class = "btn-warning")
                     )
                  )
                )
              ),
              column(6, 
                div(class="form-group",
                  tags$label(tagList("Control (b2) BAMs", info_icon("BAM files representing the control or 'reference' group. Comma-separated paths."))),
                  div(style="display:flex; align-items:flex-start;",
                     textAreaInput(ns("b2"), label = NULL, placeholder = "/path/to/ctrl1.bam, ...", width = "100%", height = "100px"),
                     div(style="margin-left:5px; display:flex; flex-direction:column;", 
                         actionButton(ns("browse_b2"), "Browse", icon = icon("folder-open"), class = "btn-info", style="margin-bottom:5px;"),
                         actionButton(ns("ws_b2"), "Worksp", icon = icon("box-open"), class = "btn-warning")
                     )
                  )
                )
              )
            ),
            h4("2. Parameters", info_icon("Configure rMATS run settings including read properties and discovery options.")),
            fluidRow(
              column(4,
                div(class="form-group",
                  tags$label(tagList("GTF Annotation", info_icon("GTF annotation file containing gene models used to identify splicing events. User can use provided GTF files with prebuilt indexes (present in tool_directory/reference_indexes"))),
                  div(style="display:flex",
                    textInput(ns("gtf"), label = NULL, placeholder = "/path/to/genomic.gtf", width = "100%"),
                    div(style="margin-left:5px", actionButton(ns("browse_gtf"), "Browse", icon = icon("folder-open"), class = "btn-info"))
                  )
                )
              ),
              column(4, numericInput(ns("readLength"), label = tagList("Read Length", info_icon("The length of the sequencing reads. Required for accurate splicing analysis.")), value = 150, min = 1)),
              column(4, numericInput(ns("nthread"), label = tagList("Threads", info_icon("Number of CPU cores to allocate for the rMATS analysis.")), value = 4, min = 1))
            ),
            fluidRow(
              column(4, checkboxInput(ns("libType"), label = tagList("Paired-End Mode", info_icon("Check if the input reads are paired-end. If unchecked, single-end mode is used.")), value = TRUE)),
              column(4, checkboxInput(ns("novelSS"), label = tagList("Detect Novel Splice Sites", info_icon("Enable discovery of splicing events not present in the GTF annotation.")), value = FALSE)),
              column(4, checkboxInput(ns("variableReadLength"), label = tagList("Variable Read Length", info_icon("Allow reads with different lengths in the input BAM files.")), value = FALSE))
            ),
            hr(),
            h4("3. Command Preview"),
            verbatimTextOutput(ns("cmd_preview")),
            hr(),
               uiOutput(ns("run_button_ui")),
            span(textOutput(ns("status_msg")), style = "margin-left: 10px; font-weight: bold;")
          )
        ),
        fluidRow(
          box(
            title = tagList("Console Output", info_icon("Real-time logs from the rMATS process.")), 
            status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE,
            tags$div(
              style = "background-color: #000; color: #0f0; padding: 10px; font-family: monospace; height: 300px; overflow-y: scroll;",
              verbatimTextOutput(ns("live_log"))
            )
          )
        )
      ),
      
      # --- Tab 2: Results (Merged from Event-based Analysis) ---
      tabPanel("Results",
        mod_event_based_ui(ns("results"))
      )
    )
  )
}

mod_rmats_server <- function(id, manager) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Event Constants
    EVENT_TYPES <- c("SE", "MXE", "A3SS", "A5SS", "RI")
    
    # Reactive values
    rv <- reactiveValues(
      process = NULL,
      log_file = NULL,
      is_running = FALSE,
      log_content = "",
      output_dir = NULL,
      params = NULL,
      simulated = FALSE,
      
      # Result Data
      data = list(),
      status = "no_data",
      sashimi_image = NULL
    )
    
    # --- Initialization / Restoration ---
    observe({
      res <- manager$get_result(id)
      pars <- manager$get_params(id)
      
      if (!is.null(res)) {
        rv$output_dir <- res$output_dir
        rv$simulated <- isTRUE(res$simulated)
        if (!is.null(res$log_content)) rv$log_content <- res$log_content
        
        # Load results if available
        if (!is.null(rv$output_dir) && dir.exists(rv$output_dir)) {
           load_data(rv$output_dir)
        }
      }
      
      if (!is.null(pars)) {
        rv$params <- pars
        if (manager$get_status(id) == "pending") {
          updateTextAreaInput(session, "b1", value = pars$b1)
          updateTextAreaInput(session, "b2", value = pars$b2)
          updateTextInput(session, "gtf", value = pars$gtf)
          updateNumericInput(session, "readLength", value = pars$readLength)
          updateNumericInput(session, "nthread", value = pars$nthread)
          updateCheckboxInput(session, "libType", value = pars$libType)
          updateCheckboxInput(session, "novelSS", value = pars$novelSS)
          if (!is.null(pars$variableReadLength)) updateCheckboxInput(session, "variableReadLength", value = pars$variableReadLength)
        }
      }
    })
    
    # --- Thread Management ---
    observe({
      total_threads <- parallel::detectCores(logical = TRUE)
      if (is.na(total_threads)) total_threads <- 4
      
      updateNumericInput(session, "nthread", max = total_threads)
      
      req(input$nthread)
      if (input$nthread > max(1, total_threads - 4)) {
        showNotification(
          paste("Warning: Using", input$nthread, "threads leaves few resources for the system (Total:", total_threads, "). System may become unresponsive."),
          type = "warning", duration = 5, id = "rmats_thread_warn"
        )
      } else {
        removeNotification("rmats_thread_warn")
      }
    })
    
    # --- Native Browser Logic ---
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
    
    observeEvent(input$browse_b1, {
      paths <- get_native_file("Select Case (b1) BAM Files", "*.bam", multi = TRUE)
      if (!is.na(paths) && length(paths) > 0) {
        val <- if(length(paths) > 1) paste(paths, collapse=", ") else paths
        current <- input$b1
        new_val <- if(current == "") val else paste(current, val, sep=", ")
        updateTextAreaInput(session, "b1", value = new_val)
      }
    })
    
    observeEvent(input$browse_b2, {
      paths <- get_native_file("Select Control (b2) BAM Files", "*.bam", multi = TRUE)
      if (!is.na(paths) && length(paths) > 0) {
        val <- if(length(paths) > 1) paste(paths, collapse=", ") else paths
        current <- input$b2
        new_val <- if(current == "") val else paste(current, val, sep=", ")
        updateTextAreaInput(session, "b2", value = new_val)
      }
    })
    
    observeEvent(input$browse_gtf, {
      path <- get_native_file("Select GTF Annotation File", "*.gtf", multi = FALSE)
      if (!is.na(path) && length(path) > 0) updateTextInput(session, "gtf", value = path)
    })
    
    # --- Workspace Scan Logic ---
    find_bams <- function() {
       if (is.null(manager$workspace_dir)) return(character(0))
       # Scan Recursively
       list.files(manager$workspace_dir, pattern = "\\.bam$", recursive = TRUE, full.names = TRUE)
    }
    
    # Handler for B1 Workspace
    observeEvent(input$ws_b1, {
       bams <- find_bams()
       if (length(bams) == 0) { showNotification("No BAM files found in workspace.", type="warning"); return() }
       
       # Label Generation: Name = Display (Simplified), Value = Path (Full)
       labels <- basename(bams)
       labels <- gsub("\\.bam$", "", labels, ignore.case=TRUE)
       labels <- gsub("\\.fastq.*", "", labels, ignore.case=TRUE)
       
       choices_vec <- bams
       names(choices_vec) <- labels
       
       showModal(modalDialog(
         title = "Select Case (b1) BAMs from Workspace",
         checkboxGroupInput(ns("sel_bams_b1"), "Files:", choices = choices_vec, width="100%"),
         footer = tagList(
           modalButton("Cancel"),
           actionButton(ns("add_b1"), "Add Selected", class="btn-primary")
         )
       ))
    })
    
    observeEvent(input$add_b1, {
       req(input$sel_bams_b1)
       current <- input$b1
       to_add <- paste(input$sel_bams_b1, collapse=", ")
       new_val <- if(current == "") to_add else paste(current, to_add, sep=", ")
       updateTextAreaInput(session, "b1", value = new_val)
       removeModal()
    })
    
    # Handler for B2 Workspace
    observeEvent(input$ws_b2, {
       bams <- find_bams()
       if (length(bams) == 0) { showNotification("No BAM files found in workspace.", type="warning"); return() }
       
       # Label Generation
       labels <- basename(bams)
       labels <- gsub("\\.bam$", "", labels, ignore.case=TRUE)
       labels <- gsub("\\.fastq.*", "", labels, ignore.case=TRUE)
       
       choices_vec <- bams
       names(choices_vec) <- labels
       
       showModal(modalDialog(
         title = "Select Control (b2) BAMs from Workspace",
         checkboxGroupInput(ns("sel_bams_b2"), "Files:", choices = choices_vec, width="100%"),
         footer = tagList(
           modalButton("Cancel"),
           actionButton(ns("add_b2"), "Add Selected", class="btn-primary")
         )
       ))
    })
    
    observeEvent(input$add_b2, {
       req(input$sel_bams_b2)
       current <- input$b2
       to_add <- paste(input$sel_bams_b2, collapse=", ")
       new_val <- if(current == "") to_add else paste(current, to_add, sep=", ")
       updateTextAreaInput(session, "b2", value = new_val)
       removeModal()
    })
    
    # --- Execution Logic ---
    
    output$run_button_ui <- renderUI({
      if (rv$is_running) {
        actionButton(ns("run_btn_disabled"), "Running...", class = "btn-success btn-lg disabled", icon = icon("spinner", class="fa-spin"))
      } else {
        actionButton(ns("run_btn"), "RUN rMATS", class = "btn-success btn-lg", icon = icon("play") )
      }
    })
    
    generated_cmd <- reactive({
      b1_val <- if (is.null(input$b1)) "" else input$b1
      b2_val <- if (is.null(input$b2)) "" else input$b2
      gtf_val <- if (is.null(input$gtf)) "" else input$gtf
      
      ws_dir <- if (is.null(manager$workspace_dir)) "./workspace" else manager$workspace_dir
      rmats_dir <- file.path(ws_dir, "rmats")
      
      cmd_parts <- c(
        "rmats.py",
        paste0("--b1 ", file.path(rmats_dir, "input/case.txt")),
        paste0("--b2 ", file.path(rmats_dir, "input/control.txt")),
        paste0("--gtf ", gtf_val),
        paste0("-t ", ifelse(input$libType, "paired", "single")),
        if(input$novelSS) "--novelSS" else NULL,
        if(input$variableReadLength) "--variable-read-length" else NULL,
        paste0("--readLength ", input$readLength),
        paste0("--od ", file.path(rmats_dir, "output")),
        paste0("--nthread ", input$nthread),
        paste0("--tmp ", file.path(rmats_dir, "tmp"))
      )
      paste(cmd_parts[!sapply(cmd_parts, is.null)], collapse = " ")
    })
    
    output$cmd_preview <- renderText({ generated_cmd() })
    
    observeEvent(input$run_btn, {
      # Validation
      valid <- TRUE
      msg <- c()
      
      shinyjs::runjs(sprintf("$('#%s').closest('.form-group').removeClass('has-error');", ns("b1")))
      shinyjs::runjs(sprintf("$('#%s').closest('.form-group').removeClass('has-error');", ns("b2")))
      shinyjs::runjs(sprintf("$('#%s').closest('.form-group').removeClass('has-error');", ns("gtf")))
      
      if (is.null(input$b1) || input$b1 == "") {
        shinyjs::runjs(sprintf("$('#%s').closest('.form-group').addClass('has-error');", ns("b1")))
        valid <- FALSE
        msg <- c(msg, "Case (b1) files are missing.")
      }
      if (is.null(input$b2) || input$b2 == "") {
        shinyjs::runjs(sprintf("$('#%s').closest('.form-group').addClass('has-error');", ns("b2")))
        valid <- FALSE
        msg <- c(msg, "Control (b2) files are missing.")
      }
      if (is.null(input$gtf) || input$gtf == "") {
        shinyjs::runjs(sprintf("$('#%s').closest('.form-group').addClass('has-error');", ns("gtf")))
        valid <- FALSE
        msg <- c(msg, "GTF file is missing.")
      }
      
      if (!valid) {
        showNotification(paste("Caution: Please fix the following issues:", paste(msg, collapse=" ")), type = "error")
        return()
      }
      
      rv$is_running <- TRUE
      
      # Enforce Workspace Structure
      ws_dir <- manager$workspace_dir
      if (is.null(ws_dir)) ws_dir <- getwd()
      
      artifacts_dir <- file.path(ws_dir, "artifacts")
      if (!dir.exists(artifacts_dir)) dir.create(artifacts_dir, recursive = TRUE)
      
      rmats_root <- file.path(artifacts_dir, "rmats")
      input_dir <- file.path(rmats_root, "input")
      output_dir <- file.path(rmats_root, "output")
      tmp_dir <- file.path(rmats_root, "tmp")
      log_dir <- file.path(rmats_root, "logs")
      
      dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
      
      rv$output_dir <- output_dir
      
      # Process paths
      b1_paths <- trimws(unlist(strsplit(input$b1, ",")))
      b2_paths <- trimws(unlist(strsplit(input$b2, ",")))
      
      writeLines(paste(b1_paths, collapse = ","), file.path(input_dir, "case.txt"))
      writeLines(paste(b2_paths, collapse = ","), file.path(input_dir, "control.txt"))
      
      # Reconstruct command for record using actual paths
      final_cmd_str <- paste(
        "rmats.py",
        paste0("--b1 ", file.path(input_dir, "case.txt")),
        paste0("--b2 ", file.path(input_dir, "control.txt")),
        paste0("--gtf ", input$gtf),
        paste0("-t ", ifelse(input$libType, "paired", "single")),
        if(input$novelSS) "--novelSS" else NULL,
        if(input$variableReadLength) "--variable-read-length" else NULL,
        paste0("--readLength ", input$readLength),
        paste0("--od ", output_dir),
        paste0("--nthread ", input$nthread),
        paste0("--tmp ", tmp_dir),
        collapse = " "
      )
      
      params <- list(
        b1 = input$b1, b2 = input$b2, gtf = input$gtf,
        readLength = input$readLength, nthread = input$nthread,
        libType = input$libType, novelSS = input$novelSS, 
        variableReadLength = input$variableReadLength,
        cmd = final_cmd_str
      )
      saveRDS(params, file.path(input_dir, "parameters.rds"))
      rv$params <- params
      
      rv$log_file <- file.path(log_dir, paste0("rmats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
      file.create(rv$log_file)
      
      # Determine rMATS executable
      local_rmats <- file.path(getwd(), "env", "bin", "rmats.py")
      rmats_cmd <- NULL
      
      if (file.exists(local_rmats)) {
        rmats_cmd <- local_rmats
      } else {
        # Fallback to system path
        sys_rmats <- tryCatch({ system("which rmats.py", intern = TRUE) }, warning = function(w) NULL, error = function(e) NULL)
        if (!is.null(sys_rmats) && length(sys_rmats) > 0 && nzchar(sys_rmats)) {
           rmats_cmd <- sys_rmats
        }
      }
      
      if (!is.null(rmats_cmd)) {
         clean_args <- c(
             "--b1", file.path(input_dir, "case.txt"),
             "--b2", file.path(input_dir, "control.txt"),
             "--gtf", input$gtf,
             "-t", ifelse(input$libType, "paired", "single"),
             if(input$novelSS) "--novelSS" else NULL,
             if(input$variableReadLength) "--variable-read-length" else NULL,
             "--readLength", as.character(input$readLength),
             "--od", output_dir,
             "--nthread", as.character(input$nthread),
             "--tmp", tmp_dir
         )
         clean_args <- clean_args[!sapply(clean_args, is.null)]
         
         # Note: processx needs the executable as the first arg. 
         # Explicitly use local python to avoid shebang issues (/usr/bin/env python not found)
         local_python <- file.path(getwd(), "env", "bin", "python")
         
         if (file.exists(local_python)) {
            rv$process <- processx::process$new(local_python, args = c(rmats_cmd, clean_args), stdout = rv$log_file, stderr = "2>&1")
         } else {
            rv$process <- processx::process$new(rmats_cmd, args = clean_args, stdout = rv$log_file, stderr = "2>&1")
         }
      } else {
         # Simulation
         sim_script <- file.path(input_dir, "sim_rmats.sh")
         writeLines(c("#!/bin/bash", "echo 'rMATS v4.1.2'", "echo 'Processing...'", "sleep 2", "echo 'Done'" ), sim_script)
         Sys.chmod(sim_script, "755")
         rv$process <- processx::process$new("bash", args = c(sim_script), stdout = rv$log_file, stderr = "2>&1")
         rv$simulated <- TRUE
         
         # Generate Mock Results for Simulation
         for(t in c("SE", "MXE", "A3SS", "A5SS", "RI")) {
             dummy_df <- data.frame(
               ID = 1:20,
               GeneID = paste0("Gene", 1:20),
               geneSymbol = paste0("Sym", 1:20),
               FDR = runif(20, 0, 0.05),
               IncLevelDifference = runif(20, 0.1, 0.9) * sample(c(-1, 1), 20, replace=TRUE),
               IncLevel1 = paste(round(runif(20),2), round(runif(20),2), sep=","),
               IncLevel2 = paste(round(runif(20),2), round(runif(20),2), sep=",")
             )
             write.table(dummy_df, file = file.path(output_dir, paste0(t, ".MATS.JC.txt") ), sep = "\t", row.names = FALSE, quote = FALSE)
         }
      }
    })
    
    # Polling
    observe({
      req(rv$is_running, rv$process)
      invalidateLater(500, session)
      if (file.exists(rv$log_file)) rv$log_content <- paste(readLines(rv$log_file, warn=FALSE), collapse="\n")
      
      if (!rv$process$is_alive()) {
        rv$is_running <- FALSE
        if (rv$process$get_exit_status() == 0) {
          
          rv$log_content <- paste0(rv$log_content, "\n[Post-Process] rMATS finished successfully.\n")
          rv$log_content <- paste0(rv$log_content, "Note: You can manually annotating gene symbols in the Results tab.\n")
          
          showNotification("Completed", type = "message")
          
          manager$set_result(id, list(
            output_dir = rv$output_dir,
            simulated = rv$simulated,
            log_content = rv$log_content
          ), rv$params)
          
          # Trigger reload for result submodule
          manager$events$workspace_loaded <- manager$events$workspace_loaded + 1
          
        } else {
          showNotification("Failed", type = "error")
        }
      }
    })
    
    output$live_log <- renderText({ rv$log_content })
    
    # --- Results Logic (Merged) ---
    # Delegated to mod_event_based module
    mod_event_based_server("results", manager)
    
  })
}
