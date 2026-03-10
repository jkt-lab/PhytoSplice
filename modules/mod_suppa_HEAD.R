library(shiny)
library(bslib)
# library(plotly)
library(DT)
library(processx)
library(shinyjs)

mod_suppa_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(id = ns("tabs"),
      # --- Tab 1: Workflow Execution ---
      tabPanel("Execution",
            br(),
            box(
              title = "SUPPA Configuration", status = "primary", solidHeader = TRUE, width = 12,
              fluidRow(
                column(6,
                  div(class="form-group",
                    tags$label("GTF Annotation"),
                    div(style="display:flex",
                      textInput(ns("gtf"), label = NULL, placeholder = "/path/to/annotation.gtf", width = "100%"),
                      div(style="margin-left:5px", actionButton(ns("browse_gtf"), "Browse", icon = icon("folder-open"), class = "btn-info"))
                    )
                  )
                ),
                column(6,
                   selectInput(ns("event_type"), "Event Types to Generate", 
                               choices = c("SE", "SS", "MX", "RI", "FL"), 
                               selected = c("SE", "MX", "RI"), multiple = TRUE)
                )
              ),
              hr(),
              h5("Experimental Design (Define Groups)"),
              p("Select samples from the Salmon output (merged TPM) to assign to conditions."),
              fluidRow(
                column(6,
                  selectInput(ns("ctrl_samples"), "Control Samples", choices = NULL, multiple = TRUE)
                ),
                column(6,
                  selectInput(ns("case_samples"), "Case Samples", choices = NULL, multiple = TRUE)
                )
              ),
              div(style="margin-bottom: 10px;",
                 actionButton(ns("refresh_samples"), "Refresh Sample List", icon = icon("sync"), class="btn-xs")
              ),
              hr(),
              h4("Execution"),
              verbatimTextOutput(ns("cmd_preview")),
              br(),
              uiOutput(ns("run_suppa_ui")),
              span(textOutput(ns("suppa_status")), style = "margin-left: 10px; font-weight: bold;")
            ),
            fluidRow(
              box(
                title = "Console Output", status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE,
                div(class="log-container",
                  style = "background-color: #000; color: #0f0; padding: 10px; font-family: monospace; height: 250px; overflow-y: scroll;",
                  verbatimTextOutput(ns("live_log"))
                )
              )
            )
      ),
      
      # --- Tab 2: Results ---
      tabPanel("Results",
        fluidRow(
          box(title = "Result Controls", status = "warning", solidHeader = TRUE, width = 12,
              column(4, numericInput(ns("pval_cut"), "P-value Threshold", value = 0.05, min = 0, max = 1, step = 0.01)),
              column(4, numericInput(ns("dpsi_cut"), "Delta PSI Threshold", value = 0.1, min = 0, max = 1, step = 0.05)),
              column(4, actionButton(ns("refresh_results"), "Reload Results", icon = icon("sync"), class="btn-primary", style="margin-top:25px;"))
          )
        ),
        fluidRow(
          tabBox(width = 12,
            tabPanel("Volcano Plot", 
                     plotlyOutput(ns("volcano_plot"), height = "600px")
            ),
            tabPanel("PSI Distribution", 
                     plotlyOutput(ns("psi_dist_plot"), height = "500px"),
                     p("Distribution of PSI values for significant events across conditions.")
            ),
            tabPanel("Significant Events Table",
                     DTOutput(ns("result_table"))
            )
          )
        )
      )
    )
  )
}

mod_suppa_server <- function(id, manager) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rv <- reactiveValues(
      process = NULL,
      log_file = NULL,
      is_running = FALSE,
      log_content = "",
      output_dir = NULL,
      results = NULL,
      available_samples = character(0)
    )
    
    # --- Helpers ---
    get_native_file <- function(title = "Select File", filter_pattern = NULL, dir=FALSE) {
      if (nzchar(Sys.which("zenity"))) {
        tryCatch({
          if (dir) {
            cmd <- paste0('zenity --file-selection --directory --title="', title, '"')
          } else {
            filter_arg <- if (!is.null(filter_pattern)) paste0('--file-filter="', filter_pattern, '"') else ""
            cmd <- paste0('zenity --file-selection --title="', title, '" ', filter_arg)
          }
          path <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
          if (length(path) > 0 && nzchar(path)) return(path)
        }, error = function(e) NA)
      }
      return(NA)
    }
    
    # --- Browsers ---
    observeEvent(input$browse_gtf, {
      path <- get_native_file("Select GTF Annotation", "*.gtf")
      if (!is.na(path)) updateTextInput(session, "gtf", value = path)
    })
    
    # --- Initialization ---
    observeEvent(manager$events$workspace_loaded, {
      # Try to load existing merged TPMs to populate sample lists
      ws_dir <- manager$workspace_dir
      if (!is.null(ws_dir)) {
        tpm_file <- file.path(ws_dir, "artifacts", "salmon", "all_tpm.txt")
        if (file.exists(tpm_file)) {
          df <- read.table(tpm_file, header=TRUE, nrows=1)
          rv$available_samples <- setdiff(names(df), "target_id")
        }
      }
      
      # Restore Params and Logs
      pars <- manager$get_params(id)
      res <- manager$get_result(id)
      
      if (!is.null(res) && !is.null(res$log_content)) {
         rv$log_content <- res$log_content
      }
      
      if (!is.null(pars)) {
        if (!is.null(pars$gtf)) updateTextInput(session, "gtf", value = pars$gtf)
        if (!is.null(pars$salmon_index)) updateTextInput(session, "salmon_index", value = pars$salmon_index)
        if (!is.null(pars$ctrl_samples)) updateSelectInput(session, "ctrl_samples", selected = pars$ctrl_samples)
        if (!is.null(pars$case_samples)) updateSelectInput(session, "case_samples", selected = pars$case_samples)
      } else {
        # Autofill GTF from rMATS
        rmats_pars <- manager$get_params("rmats")
        if (!is.null(rmats_pars) && !is.null(rmats_pars$gtf)) updateTextInput(session, "gtf", value = rmats_pars$gtf)
      }
    })
    
    observe({
      updateSelectInput(session, "ctrl_samples", choices = rv$available_samples, selected = input$ctrl_samples)
      updateSelectInput(session, "case_samples", choices = rv$available_samples, selected = input$case_samples)
    })
    
    observeEvent(input$refresh_samples, {
      ws_dir <- if (is.null(manager$workspace_dir)) getwd() else manager$workspace_dir
      tpm_file <- file.path(ws_dir, "artifacts", "salmon", "all_tpm.txt")
      if (file.exists(tpm_file)) {
        df <- tryCatch(read.table(tpm_file, header=TRUE, nrows=1), error=function(e) NULL)
        if(!is.null(df)) {
          rv$available_samples <- setdiff(names(df), "target_id")
          showNotification("Sample list refreshed from all_tpm.txt", type="message")
        }
      } else {
        showNotification("all_tpm.txt not found. Run Salmon first.", type="warning")
      }
    })

    output$run_suppa_ui <- renderUI({
      if (rv$is_running) {
        actionButton(ns("run_suppa_disabled"), "Running...", class = "btn-primary btn-lg disabled", icon = icon("spinner", class="fa-spin"))
      } else {
        actionButton(ns("run_suppa"), "Run SUPPA Analysis", class = "btn-primary btn-lg", icon = icon("play"))
      }
    })

    # --- SUPPA Logic ---
    generated_cmd <- reactive({
      paste(
        "# 1. Generate Events",
        paste("suppa.py generateEvents -i", input$gtf, "-o events -e", paste(input$event_type, collapse=","), "-f ioe"),
        "# 2. Calculate PSI",
        paste("suppa.py psiPerEvent -i events.ioe -e all_tpm.txt -o project_PSI"),
        "# 3. DiffSplice",
        paste("suppa.py diffSplice -m project_PSI.psi -i events.ioe -p all_tpm.txt -gc Control -pc Case -o diff_results"),
        sep="\n"
      )
    })
    output$cmd_preview <- renderText({ generated_cmd() })
    
    observeEvent(input$run_suppa, {
      req(input$gtf, input$ctrl_samples, input$case_samples)
      
      ws_dir <- if (is.null(manager$workspace_dir)) getwd() else manager$workspace_dir
      suppa_dir <- file.path(ws_dir, "artifacts", "suppa")
      dir.create(suppa_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Check for Merged TPM
      salmon_tpm <- file.path(ws_dir, "artifacts", "salmon", "all_tpm.txt")
      if (!file.exists(salmon_tpm)) {
        showNotification("Merged TPM file (all_tpm.txt) not found in artifacts/salmon. Run Salmon first.", type="error")
        return()
      }
      
      # Auto-Repair: Sanitize Header (Remove spaces to prevent SUPPA parsing errors)
      tryCatch({
        # Attempt to read as tab-separated
        tpm_df <- read.table(salmon_tpm, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE, quote="")
        
        orig_cols <- colnames(tpm_df)
        # Replace spaces with underscores
        new_cols <- gsub("\\s+", "_", orig_cols)
        
        if (any(orig_cols != new_cols)) {
           colnames(tpm_df) <- new_cols
           write.table(tpm_df, salmon_tpm, sep="\t", quote=FALSE, row.names=FALSE)
           showNotification("Repaired all_tpm.txt header (removed spaces).", type="warning")
        }
      }, error = function(e) {
         # If read fails, we can't fix it safely. Warn user.
         showNotification(paste("Warning: Could not verify TPM file integrity:", e$message), type="warning")
      })
      
      rv$output_dir <- suppa_dir
      rv$log_file <- file.path(suppa_dir, "suppa_run.log")
      rv$is_running <- TRUE
      
      # Save Params
      manager$state$params[[id]] <- list(
        gtf = input$gtf, ctrl_samples = input$ctrl_samples, case_samples = input$case_samples
      )
      
      # Generate Script
      script_file <- file.path(suppa_dir, "run_suppa.sh")
      env_bin <- file.path(getwd(), "env", "bin")
      
      events_prefix <- file.path(suppa_dir, "events")
      psi_prefix <- file.path(suppa_dir, "project")
      diff_prefix <- file.path(suppa_dir, "diff_results")
      
      # Comma-sep strings for samples (SUPPA expects sample names matching header)
      ctrl_str <- paste(input$ctrl_samples, collapse=",")
      case_str <- paste(input$case_samples, collapse=",")
      types_str <- paste(input$event_type, collapse=" ") # space separated for suppa? NO, check syntax.
      # suppa.py generateEvents -e SE SS ... (space separated)
      
      script_lines <- c(
        "#!/bin/bash",
        "set -e",
        paste0("export PATH=\"", env_bin, ":$PATH\""),
        paste0("GTF=\"", input$gtf, "\""),
        paste0("TPM=\"", salmon_tpm, "\""),
        "",
        "echo \"Step 1: Generate Events\"",
        paste0("suppa.py generateEvents -i \"$GTF\" -o \"", events_prefix, "\" -e ", paste(input$event_type, collapse=" "), " -f ioe"),
        "",
        "echo \"Step 2: PSI Calculation\"",
        # Merge IOE files? generateEvents produces one per type (events_SE.ioe). 
        # psiPerEvent can take multiple or we loop. 
        # Easier to concat them for single run? 
        # "awk 'FNR==1 && NR!=1{next;}{print}' events_*.ioe > events.ioe"
        paste0("cat \"", events_prefix, "_\"*.ioe | awk '!seen[$0]++' > \"", suppa_dir, "/all_events.ioe\""),
        paste0("suppa.py psiPerEvent -i \"", suppa_dir, "/all_events.ioe\" -e \"$TPM\" -o \"", psi_prefix, "\""),
        "",
        "echo \"Step 3: Differential Splicing\"",
        paste0("suppa.py diffSplice -m \"", psi_prefix, ".psi\" -i \"", suppa_dir, "/all_events.ioe\" -p \"$TPM\" -gc ", ctrl_str, " -pc ", case_str, " -o \"", diff_prefix, "\"")
      )
      
      writeLines(script_lines, script_file)
      Sys.chmod(script_file, "755")
      
      rv$log_content <- paste0("Starting SUPPA Workflow: ", script_file, "\n")
      rv$process <- processx::process$new("bash", args = c(script_file), stdout = rv$log_file, stderr = "2>&1")
    })
    
    # Polling Observer (Outside observeEvent)
    observe({
      req(rv$is_running, rv$process)
      invalidateLater(500, session)
      
      # Update Log
      if(file.exists(rv$log_file)) {
        lines <- tryCatch(readLines(rv$log_file, warn=FALSE), error=function(e) character(0))
        rv$log_content <- paste(lines, collapse="\n")
      }
      
      if (!rv$process$is_alive()) {
        rv$is_running <- FALSE
        
        if (rv$process$get_exit_status() == 0) {
           showNotification("SUPPA Analysis Complete!", type="message")
           rv$log_content <- paste0(rv$log_content, "\nDone.\n")
           load_results(rv$output_dir)
           manager$set_result(id, list(output_dir = rv$output_dir, log_content = rv$log_content))
        } else {
           showNotification("SUPPA Analysis Failed.", type="error")
           rv$log_content <- paste0(rv$log_content, "\nFAILED.\n")
           manager$set_result(id, list(output_dir = rv$output_dir, log_content = rv$log_content))
        }
        
        # Cleanup
        rv$process <- NULL
      }
    })

    # --- Results (Existing Logic Preserved/Adapted) ---
    load_results <- function(dir) {
      rv$results <- load_suppa_data(dir)
    }
    
    output$live_log <- renderText({ rv$log_content })
    
    filtered_data <- reactive({
      req(rv$results)
      filter_suppa_data(rv$results, input$pval_cut, input$dpsi_cut)
    })
    
    output$result_table <- renderDT({
      req(filtered_data())
      filtered_data() %>% dplyr::select(Event_id, mean_dPSI, p_val) %>% datatable(selection = 'single')
    })
    
    output$volcano_plot <- renderPlotly({
      req(rv$results)
      df <- rv$results
      df$Sig <- ifelse(df$p_val < input$pval_cut & abs(df$mean_dPSI) > input$dpsi_cut, "Significant", "NS")
      p <- ggplot(df, aes(x=mean_dPSI, y=-log10(p_val), color=Sig, text=Event_id)) +
        geom_point(alpha=0.6) + theme_minimal() +
        geom_hline(yintercept = -log10(input$pval_cut), linetype="dashed", color="gray") +
        geom_vline(xintercept = c(-input$dpsi_cut, input$dpsi_cut), linetype="dashed", color="gray")
      ggplotly(p, tooltip="text")
    })
    
    output$psi_dist_plot <- renderPlotly({
      # Basic violin plot of significant events
      req(filtered_data())
      # (Note: This requires access to the PSI values which are in .psivec, loaded by helper)
      NULL 
    })
  })
}
