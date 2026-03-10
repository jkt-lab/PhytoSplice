library(shiny)
library(bslib)
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
              title = tagList("SUPPA Configuration", info_icon("What it does: Sets up parameters to extract splicing events and calculate differential splicing.\nHow it works: SUPPA first generates an IOE (Isoform Overlapping Event) file from the GTF, then uses Salmon TPMs to calculate PSI (Percent Spliced In) and diffSplice.")), status = "primary", solidHeader = TRUE, width = 12,
              fluidRow(
                column(6,
                  div(class="form-group",
                    tags$label(tagList("GTF Annotation", info_icon("What it does: Defines the gene and transcript structures.\nHow it works: Used to generate the reference list of all possible alternative splicing events."))),
                    div(style="display:flex",
                      textInput(ns("gtf"), label = NULL, placeholder = "/path/to/annotation.gtf", width = "100%"),
                      div(style="margin-left:5px", actionButton(ns("browse_gtf"), "Browse", icon = icon("folder-open"), class = "btn-info"))
                    )
                  )
                ),
                column(6,
                   selectInput(ns("event_type"), label = tagList("Event Types to Generate", info_icon("What it does: Selects which types of alternative splicing events to look for.\nHow it works: SE: Skipping Exon, SS: Alternative Splice Site (A5SS/A3SS), MX: Mutually Exclusive, RI: Retained Intron, FL: Alternative First/Last exon.")), 
                               choices = c("SE", "SS", "MX", "RI", "FL"), 
                               selected = c("SE", "MX", "RI"), multiple = TRUE)
                )
              ),
              hr(),
              h5(tagList("Experimental Design (Define Groups)", info_icon("What it does: Assigns biological samples to experimental conditions for statistical testing.\nHow it works: Reads available samples from the Salmon merged TPM file."))),
              p("Select samples from the Salmon output (merged TPM) to assign to conditions."),
              fluidRow(
                column(6,
                  selectInput(ns("ctrl_samples"), label = tagList("Control Samples", info_icon("What it does: The baseline or reference group.")), choices = NULL, multiple = TRUE)
                ),
                column(6,
                  selectInput(ns("case_samples"), label = tagList("Case Samples", info_icon("What it does: The treatment or experimental group.")), choices = NULL, multiple = TRUE)
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
                title = tagList("Console Output", info_icon("What it does: Shows execution progress from the SUPPA python scripts.")), status = "info", solidHeader = TRUE, width = 12, collapsible = TRUE,
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
          box(title = tagList("Result Controls", info_icon("What it does: Filters the SUPPA output for visualizations.")), status = "warning", solidHeader = TRUE, width = 12,
              column(4, numericInput(ns("pval_cut"), label = tagList("P-value Threshold", info_icon("What it does: Filters out events with a p-value strictly greater than this number.")), value = 0.05, min = 0, max = 1, step = 0.01)),
              column(4, numericInput(ns("dpsi_cut"), label = tagList("Delta PSI Threshold", info_icon("What it does: Filters out events with an absolute difference in PSI strictly smaller than this number.")), value = 0.1, min = 0, max = 1, step = 0.05)),
              column(4, actionButton(ns("refresh_results"), "Reload Results", icon = icon("sync"), class="btn-primary", style="margin-top:25px;"))
          )
        ),
        fluidRow(
          tabBox(width = 12,
            tabPanel("dPSI Distribution",
                     plotly::plotlyOutput(ns("dpsi_hist_plot"), height = "500px"),
                     p("Global distribution of splicing changes (Delta PSI).")
            ),
            tabPanel("Volcano Plot", 
                     plotly::plotlyOutput(ns("volcano_plot"), height = "600px")
            ),
            tabPanel("Heatmap",
                     plotly::plotlyOutput(ns("heatmap_plot"), height = "600px"),
                     p("PSI values for top 50 significant events.")
            ),
            tabPanel("Event Types",
                     plotly::plotlyOutput(ns("event_type_plot"), height = "600px"),
                     p("Distribution of splicing changes by event type.")
            ),
            tabPanel("PSI Distribution (Single Event)", 
                     plotly::plotlyOutput(ns("psi_dist_plot"), height = "500px"),
                     p("Distribution of PSI values for selected event across conditions.")
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
    
    get_event_type <- function(event_ids) {
      # Extract type between first ; and : (e.g. gene;SE:...)
      sub(".*;([A-Z0-9]+):.*", "\\1", event_ids)
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
      
      # Determine output directory (priority: result state > guessed from workspace)
      target_dir <- NULL
      if (!is.null(res) && !is.null(res$output_dir)) {
         target_dir <- res$output_dir
      } else if (!is.null(manager$workspace_dir)) {
         guess_dir <- file.path(manager$workspace_dir, "artifacts", "suppa")
         if (dir.exists(guess_dir)) target_dir <- guess_dir
      }
      
      if (!is.null(target_dir)) {
         rv$output_dir <- target_dir
         load_results(rv$output_dir)
      }
      
      if (!is.null(pars)) {
        if (!is.null(pars$gtf)) updateTextInput(session, "gtf", value = pars$gtf)
        if (!is.null(pars$event_type)) updateSelectInput(session, "event_type", selected = pars$event_type)
        if (!is.null(pars$ctrl_samples)) updateSelectInput(session, "ctrl_samples", selected = pars$ctrl_samples)
        if (!is.null(pars$case_samples)) updateSelectInput(session, "case_samples", selected = pars$case_samples)
        
        # Result controls
        if (!is.null(pars$pval_cut)) updateNumericInput(session, "pval_cut", value = pars$pval_cut)
        if (!is.null(pars$dpsi_cut)) updateNumericInput(session, "dpsi_cut", value = pars$dpsi_cut)
      } else {
        # Autofill GTF from rMATS
        rmats_pars <- manager$get_params("rmats")
        if (!is.null(rmats_pars) && !is.null(rmats_pars$gtf)) updateTextInput(session, "gtf", value = rmats_pars$gtf)
      }
    })
    
    # --- Persistence: Auto-save inputs to manager state ---
    observe({
      manager$state$params[[id]] <- list(
        gtf = input$gtf,
        event_type = input$event_type,
        ctrl_samples = input$ctrl_samples,
        case_samples = input$case_samples,
        pval_cut = input$pval_cut,
        dpsi_cut = input$dpsi_cut
      )
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

    observeEvent(input$refresh_results, {
      if (is.null(rv$output_dir) && !is.null(manager$workspace_dir)) {
        guess_dir <- file.path(manager$workspace_dir, "artifacts", "suppa")
        if (dir.exists(guess_dir)) rv$output_dir <- guess_dir
      }
      
      req(rv$output_dir)
      load_results(rv$output_dir)
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
        showNotification("Merged TPM file (all_tpm.txt) not found. Run Salmon first.", type="error")
        return()
      }
      
      # --- Pre-processing: Create Condition-Specific TPM Files ---
      # We use R to split the TPM file because 'suppa.py joinFiles' is not reliable for raw quant.sf
      tryCatch({
        tpm_df <- read.table(salmon_tpm, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE, quote="")
        
        # --- ID Harmonization ---
        # Clean Transcript IDs to match GTF (Remove prefixes, extra info in headers)
        tpm_df$target_id <- gsub("^transcript:", "", tpm_df$target_id, ignore.case = TRUE)
        tpm_df$target_id <- gsub("\\|.*", "", tpm_df$target_id)
        tpm_df$target_id <- trimws(tpm_df$target_id)
        
        # Validate samples exist
        missing_ctrl <- setdiff(input$ctrl_samples, names(tpm_df))
        missing_case <- setdiff(input$case_samples, names(tpm_df))
        
        if (length(missing_ctrl) > 0 || length(missing_case) > 0) {
          showNotification("Selected samples not found in TPM file.", type="error")
          return()
        }
        
        # Write Control TPM
        ctrl_df <- tpm_df[, c("target_id", input$ctrl_samples)]
        # SUPPA requires header with only sample names (no column name for the ID column)
        # Format:
        # sample1\tsample2
        # tx1\t10\t20
        
        ctrl_file <- file.path(suppa_dir, "control.tpm")
        cat(paste(input$ctrl_samples, collapse="\t"), "\n", file=ctrl_file)
        write.table(ctrl_df, ctrl_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
        
        # Write Case TPM
        case_df <- tpm_df[, c("target_id", input$case_samples)]
        
        case_file <- file.path(suppa_dir, "case.tpm")
        cat(paste(input$case_samples, collapse="\t"), "\n", file=case_file)
        write.table(case_df, case_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
        
      }, error = function(e) {
        showNotification(paste("Error preparing TPM files:", e$message), type="error")
        return()
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
      # Output prefixes for PSI calculation
      psi_ctrl_prefix <- file.path(suppa_dir, "control")
      psi_case_prefix <- file.path(suppa_dir, "case")
      diff_prefix <- file.path(suppa_dir, "diff_results")
      
      script_lines <- c(
        "#!/bin/bash",
        "set -e",
        paste0("export PATH=\"", env_bin, ":$PATH\""),
        paste0("GTF=\"", input$gtf, "\""),
        paste0("CTRL_TPM=\"", suppa_dir, "/control.tpm\""),
        paste0("CASE_TPM=\"", suppa_dir, "/case.tpm\""),
        "",
        "echo \"Step 1: Generate Events\"",
        paste0("suppa.py generateEvents -i \"$GTF\" -o \"", events_prefix, "\" -e ", paste(input$event_type, collapse=" "), " -f ioe"),
        "",
        "echo \"Step 2: Consolidate Events\"",
        # Combine all generated IOE files into one
        paste0("cat \"", events_prefix, "_\"*.ioe | awk '!seen[$0]++' > \"", suppa_dir, "/all_events.ioe\""),
        "",
        "echo \"Step 3: PSI Calculation (Control)\"",
        paste0("suppa.py psiPerEvent -i \"", suppa_dir, "/all_events.ioe\" -e \"$CTRL_TPM\" -o \"", psi_ctrl_prefix, "\""),
        "",
        "echo \"Step 4: PSI Calculation (Case)\"",
        paste0("suppa.py psiPerEvent -i \"", suppa_dir, "/all_events.ioe\" -e \"$CASE_TPM\" -o \"", psi_case_prefix, "\""),
        "",
        "echo \"Step 5: Differential Splicing\"",
        # Note: psiPerEvent adds .psi extension automatically
        paste0("suppa.py diffSplice -m empirical -i \"", suppa_dir, "/all_events.ioe\" --psi \"", psi_ctrl_prefix, ".psi\" \"", psi_case_prefix, ".psi\" --tpm \"$CTRL_TPM\" \"$CASE_TPM\" -o \"", diff_prefix, "\"")
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
      if (is.null(dir) || !dir.exists(dir)) {
        showNotification("SUPPA output directory not found.", type = "error")
        return(NULL)
      }
      
      res <- load_suppa_data(dir)
      if (is.null(res)) {
        showNotification("Could not parse SUPPA results in the artifacts directory. Check if files exist.", type = "warning")
      } else {
        # Try to reload log content if available
        log_path <- file.path(dir, "suppa_run.log")
        if (file.exists(log_path)) {
           log_lines <- tryCatch(readLines(log_path, warn=FALSE), error=function(e) NULL)
           if (!is.null(log_lines)) rv$log_content <- paste(log_lines, collapse="\n")
        }
        
        # Update available samples from psivec headers
        cols <- names(res$psivec)
        rv$available_samples <- setdiff(cols, "Event_id")
        showNotification("SUPPA results loaded successfully.", type = "message")
        
        # Sync with AnalysisManager to update progress indicator
        if (manager$get_status(id) != "completed") {
           manager$set_result(id, list(output_dir = dir, log_content = rv$log_content))
        }
      }
      rv$results <- res
    }
    
    output$live_log <- renderText({ rv$log_content })
    
    filtered_data <- reactive({
      req(rv$results, rv$results$dpsi)
      filter_suppa_data(rv$results$dpsi, input$pval_cut, input$dpsi_cut)
    })
    
    output$result_table <- renderDT({
      shiny::validate(shiny::need(rv$results, "No results loaded."))
      df <- filtered_data()
      if (is.null(df)) return(NULL)
      df %>% dplyr::select(Event_id, mean_dPSI, p_val) %>% datatable(selection = 'single')
    })
    
    output$dpsi_hist_plot <- renderPlotly({
      shiny::validate(shiny::need(rv$results, "No results loaded."))
      
      # Filter by p-value only, so we can show dPSI distribution for all significant events
      # highlighting those that pass/fail the dPSI cutoff
      df <- rv$results$dpsi
      if (is.null(df)) return(NULL)
      
      # Filter by P-value
      if ("p_val" %in% names(df)) {
        df <- df[df$p_val < input$pval_cut, ]
      }
      
      if (nrow(df) == 0) {
        return(plotly_empty() %>% layout(title = "No events pass the p-value filter"))
      }
      
      df <- df[!is.na(df$mean_dPSI), ]
      
      # Define Status based on dPSI cutoff
      df$Status <- ifelse(abs(df$mean_dPSI) >= input$dpsi_cut, "Pass", "Below Cutoff")
      
      plot_ly(data = df, x = ~mean_dPSI, type = "histogram", 
              color = ~Status,
              colors = c("Pass" = "#377EB8", "Below Cutoff" = "orange"),
              marker = list(line = list(color = "white", width = 0.5))) %>%
        layout(
          title = "Significant Delta PSI Distribution",
          xaxis = list(title = "Mean Delta PSI"),
          yaxis = list(title = "Count"),
          barmode = "stack"
        )
    })
    
    output$heatmap_plot <- renderPlotly({
      shiny::validate(shiny::need(rv$results, "No results loaded."))
      
      # Top 50 significant events
      sig <- filtered_data()
      if (nrow(sig) == 0) return(plotly_empty() %>% layout(title = "No significant events found"))
      
      top50 <- head(sig[order(sig$p_val), ], 50)
      
      # Merge with PSI data
      psivec <- rv$results$psivec
      merged <- merge(top50[, "Event_id", drop=FALSE], psivec, by="Event_id")
      
      if (nrow(merged) == 0) return(NULL)
      
      # Matrix for heatmap (exclude Event_id)
      # Identify sample columns: intersect all columns with ctrl/case samples
      samples <- c(input$ctrl_samples, input$case_samples)
      valid_samples <- intersect(names(merged), samples)
      
      if (length(valid_samples) == 0) {
         # Fallback: Use all numeric columns
         valid_samples <- names(merged)[sapply(merged, is.numeric)]
      }
      
      if (length(valid_samples) == 0) return(NULL)
      
      mat <- as.matrix(merged[, valid_samples])
      rownames(mat) <- merged$Event_id
      
      # Simple Z-score scaling per row for better visualization
      mat_scaled <- t(scale(t(mat)))
      
      plot_ly(z = mat_scaled, x = valid_samples, y = rownames(mat), type = "heatmap", colors = "RdBu") %>%
        layout(
          title = "Top 50 Significant Events (Z-score of PSI)",
          xaxis = list(title = "Samples"),
          yaxis = list(title = "Events", showticklabels = FALSE)
        )
    })
    
    output$event_type_plot <- renderPlotly({
      shiny::validate(shiny::need(rv$results, "No results loaded."))
      df <- filtered_data()
      if (is.null(df) || nrow(df) == 0) {
        return(plotly_empty() %>% layout(title = "No events pass the current filters"))
      }
      
      df <- df[!is.na(df$mean_dPSI), ]
      
      df$Type <- get_event_type(df$Event_id)
      
      # 1. Bar Plot (Counts)
      counts <- as.data.frame(table(df$Type))
      
      p1 <- plot_ly(data = counts, x = ~Var1, y = ~Freq, type = "bar", marker = list(color = "#4DAF4A"), name = "Count") %>%
        layout(yaxis = list(title = "Count"))
      
      # 2. Box Plot (dPSI)
      p2 <- plot_ly(data = df, x = ~Type, y = ~mean_dPSI, type = "box", color = ~Type, showlegend = FALSE) %>%
        layout(yaxis = list(title = "Delta PSI"))
      
      subplot(p1, p2, nrows = 2, shareX = TRUE, titleY = TRUE) %>%
        layout(title = "Splicing Changes by Event Type (Filtered)")
    })
    
    output$volcano_plot <- renderPlotly({
      shiny::validate(shiny::need(rv$results, "No results loaded."))
      df <- rv$results$dpsi
      if (is.null(df) || nrow(df) == 0) return(NULL)
      
      # Basic cleaning
      df$mean_dPSI <- as.numeric(df$mean_dPSI)
      df$p_val <- as.numeric(df$p_val)
      # Only remove rows with invalid stats, preserve rows with missing sample PSIs
      df <- df[!is.na(df$mean_dPSI) & !is.na(df$p_val), ]
      
      if (nrow(df) == 0) return(NULL)
      
      # Significance Grouping
      df$Status <- "Not Sig"
      df$Status[df$p_val < input$pval_cut & df$mean_dPSI >= input$dpsi_cut] <- "Up"
      df$Status[df$p_val < input$pval_cut & df$mean_dPSI <= -input$dpsi_cut] <- "Down"
      
      # Colors
      pal <- c("Not Sig" = "grey", "Up" = "#E41A1C", "Down" = "#377EB8")
      
      # Threshold Lines logic
      max_y <- max(-log10(df$p_val)[is.finite(-log10(df$p_val))], na.rm=TRUE)
      min_x <- min(df$mean_dPSI, na.rm=TRUE)
      max_x <- max(df$mean_dPSI, na.rm=TRUE)
      
      plot_ly(data = df, x = ~mean_dPSI, y = ~-log10(p_val), 
              type = 'scatter', mode = 'markers',
              color = ~Status, colors = pal,
              text = ~paste("<b>Event:</b>", Event_id, 
                            "<br><b>dPSI:</b>", round(mean_dPSI, 3), 
                            "<br><b>P-val:</b>", formatC(p_val, format = "e", digits = 2)),
              hoverinfo = "text",
              marker = list(opacity = 0.7, size = 6)) %>%
        layout(
          title = "Volcano Plot",
          xaxis = list(title = "Mean Delta PSI"),
          yaxis = list(title = "-Log10 P-value")
        )
    })
    
    output$psi_dist_plot <- renderPlotly({
      shiny::validate(shiny::need(rv$results, "No results loaded."))
      
      # Check selection
      s <- input$result_table_rows_selected
      if (length(s) == 0) {
        return(plotly_empty() %>% layout(title = "Select an event from the table to view PSI distribution"))
      }
      
      # Get Event ID (map from filtered data row index to ID)
      dat <- filtered_data()
      event_id <- dat$Event_id[s]
      
      # Get Full Data row from PSIvec
      # The dataframe is in rv$results$psivec
      psivec <- rv$results$psivec
      row_data <- psivec[psivec$Event_id == event_id, ]
      
      if(nrow(row_data) == 0) return(NULL)
      
      # Extract PSI values
      # We rely on column names matching the input samples
      ctrl_cols <- intersect(names(row_data), input$ctrl_samples)
      case_cols <- intersect(names(row_data), input$case_samples)
      
      if(length(ctrl_cols) == 0 || length(case_cols) == 0) {
         # Fallback search for control/case in names
         ctrl_cols <- grep("control", names(row_data), ignore.case = TRUE, value = TRUE)
         case_cols <- grep("case", names(row_data), ignore.case = TRUE, value = TRUE)
      }
      
      if(length(ctrl_cols) == 0 || length(case_cols) == 0) {
         return(plotly_empty() %>% layout(title = "Sample columns not found in results"))
      }
      
      ctrl_vals <- as.numeric(row_data[1, ctrl_cols])
      case_vals <- as.numeric(row_data[1, case_cols])
      
      plot_df <- data.frame(
        Sample = c(rep("Control", length(ctrl_vals)), rep("Case", length(case_vals))),
        PSI = c(ctrl_vals, case_vals)
      )
      plot_df <- na.omit(plot_df)
      
      plot_ly(plot_df, x = ~Sample, y = ~PSI, color = ~Sample, colors = c("Control"="#377EB8", "Case"="#E41A1C"),
              type = "box", boxpoints = "all", jitter = 0.3, pointpos = -1.8,
              marker = list(size = 8)) %>%
        layout(
          title = paste(event_id),
          yaxis = list(title = "PSI Value", range = c(-0.05, 1.05))
        )
    })
  })
}
