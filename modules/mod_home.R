mod_home_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Hero Section
    fluidRow(
      column(12,
        div(style = "background-color: #ecf0f5; padding: 20px; border-radius: 5px; text-align: center; margin-bottom: 20px; border-left: 5px solid #3c8dbc; position: relative;",
          h2("PhytoSplice", style = "font-weight: 700; color: #2c3e50; margin-top: 10px;"),
          p("A comprehensive, modular platform for Alternative Splicing detection, quantification, and functional enrichment.", 
            style = "font-size: 16px; color: #7f8c8d; margin-bottom: 10px;")
        )
      )
    ),
    
    # Dynamic Workspace Status
    uiOutput(ns("workspace_status_ui")),
    
    br(),
    
    # Workspace Management
    fluidRow(
      box(
        title = "Create New Workspace", status = "success", solidHeader = TRUE, width = 4, icon = icon("plus"),
        p("Start a fresh analysis project with empty state."),
        textInput(ns("create_name"), "Project Name", value = "NewProject", width = "100%"),
        br(),
        actionButton(ns("create_btn"), "Create New...", class = "btn-success btn-block", icon = icon("plus"))
      ),
      box(
        title = "Save Copy As...", status = "primary", solidHeader = TRUE, width = 4, icon = icon("copy"),
        p("Save a copy of the current session to a new project folder."),
        textInput(ns("save_name"), "New Project Name", value = "MyAnalysis_Copy", width = "100%"),
        br(),
        actionButton(ns("save_btn"), "Save Copy...", class = "btn-primary btn-block", icon = icon("save"))
      ),
      box(
        title = "Switch Workspace", status = "warning", solidHeader = TRUE, width = 4, icon = icon("folder-open"),
        p("Close the current session and load an existing workspace folder."),
        br(), br(), 
        actionButton(ns("load_btn"), "Open Workspace...", class = "btn-warning btn-block", icon = icon("folder-open"))
      )
    ),

    br(),
    hr(),
    
    # Footer
    fluidRow(
      column(12,
        div(style = "text-align: center; color: #777; padding: 20px;",
          p(strong("Developed By:"), " JKT Lab"),
          p(strong("Authors:"), " Sourabh Kumawat, Yogeshwar Dhar, Jitendra K Thakur"),
          p("© 2026 | Licensed under GPL 3.0"),
          p(a(href = "https://www.jktlab.org/", "https://www.jktlab.org/", target = "_blank"))
        )
      )
    )
  )
}

mod_home_server <- function(id, manager) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # --- Dynamic Status UI ---
    output$workspace_status_ui <- renderUI({
      # Trigger on workspace load or status update
      manager$events$workspace_loaded
      manager$events$status_changed
      
      ws_name <- if (!is.null(manager$workspace_dir)) basename(manager$workspace_dir) else "No Workspace Loaded"
      last_saved <- if (!is.null(manager$last_saved)) manager$last_saved else "Never"
      
      # Check Analysis Status
      # We check specific keys in manager$state$status
      get_badge <- function(key, label) {
        status <- manager$state$status[[key]]
        if (!is.null(status) && status == "completed") {
          tags$span(class = "badge bg-green", style="margin-right: 5px; font-size: 1em;", icon("check"), label)
        } else {
          tags$span(class = "badge bg-gray", style="margin-right: 5px; font-size: 1em;", icon("clock"), label)
        }
      }
      
      tagList(
        h3("Current Workspace Dashboard"),
        fluidRow(
          valueBox(
            value = ws_name, 
            subtitle = paste0("Active Project (Last Saved: ", format(as.POSIXct(last_saved), "%H:%M"), ")"), 
            icon = icon("folder"), 
            color = "blue", 
            width = 6
          ),
          box(
            title = "Analysis Progress", status = "success", solidHeader = TRUE, width = 6,
            div(style="padding: 10px;",
              h5("Completed Modules:"),
              get_badge("rmats", "rMATS"),
              get_badge("suppa", "SUPPA")
            )
          )
        )
      )
    })
    
    get_native_dir <- function(title = "Select Directory", default_path = getwd()) {
      if (nzchar(Sys.which("zenity"))) {
        tryCatch({
          path <- system(paste0('zenity --file-selection --directory --title="', title, '"'), intern = TRUE, ignore.stderr = TRUE)
          if (length(path) > 0 && nzchar(path)) return(path)
        }, error = function(e) NA)
      } else if (capabilities("tcltk") && !is.null(Sys.getenv("DISPLAY")) && Sys.getenv("DISPLAY") != "") {
        return(tcltk::tk_choose.dir(default = default_path, caption = title))
      }
      return(NA)
    }
    
    # Create New
    observeEvent(input$create_btn, {
      req(input$create_name)
      
      parent_dir <- get_native_dir("Select Parent Directory for New Project", default_path = getwd())
      
      if (!is.na(parent_dir)) {
        withProgress(message = "Creating New Workspace...", {
          tryCatch({
            # Reset State for Fresh Start
            manager$state <- list(
              data = list(),
              params = list(),
              status = list()
            )
            manager$is_loaded <- FALSE
            
            manager$create_workspace(parent_dir, input$create_name)
            
            showNotification(paste("New Workspace created at:", manager$workspace_dir), type = "message")
          }, error = function(e) {
            showNotification(paste("Creation failed:", e$message), type = "error")
          })
        })
      }
    })

    # Save
    observeEvent(input$save_btn, {
      req(input$save_name)
      
      parent_dir <- get_native_dir("Select Parent Directory for New Project Copy", default_path = manager$workspace_dir)
      
      if (!is.na(parent_dir)) {
        withProgress(message = "Saving Copy...", {
          tryCatch({
            # For "Save As", we might need a specific method in manager if we want to switch to it,
            # or just save the state to a new location.
            # Currently manager$save_workspace saves to manager$workspace_dir.
            # To implement "Save Copy", we should ideally copy the current state to a new dir.
            # But manager$save_workspace(target_path) isn't fully implemented in the snippet I saw 
            # (the snippet I modified earlier only took no args for save_workspace).
            # Wait, I need to check AnalysisManager.R again.
            
            # Implementation for Save Copy:
            # 1. Create new dir
            # 2. Save current state there
            # 3. (Optional) Switch to it? The user prompt implies "Save As" usually switches.
            # Let's assume we want to switch to the new one.
            
            target_path <- file.path(parent_dir, input$save_name)
            
            # We can reuse create_workspace to setup dir, then save.
            # Or manually.
            if (dir.exists(target_path)) stop("Directory already exists.")
            dir.create(target_path, recursive = TRUE)
            dir.create(file.path(target_path, "artifacts"))
            
            # Update manager to point to new one?
            # Or just save a snapshot? 
            # Usually "Save As" switches context.
            manager$workspace_dir <- target_path
            manager$save_workspace() # Now saves to new path
            
            manager$events$workspace_loaded <- manager$events$workspace_loaded + 1 # Refresh UI
            
            showNotification(paste("Workspace saved to:", target_path), type = "message", duration = 5)
          }, error = function(e) {
            showNotification(paste("Save failed:", e$message), type = "error")
          })
        })
      }
    })
    
    # Load
    observeEvent(input$load_btn, {
      target_dir <- get_native_dir("Select Workspace Folder", default_path = manager$workspace_dir)
      
      if (!is.na(target_dir)) {
        withProgress(message = "Loading Workspace...", {
          tryCatch({
            manager$load_workspace(target_dir)
            showNotification("Workspace loaded successfully!", type = "message")
          }, error = function(e) {
            showNotification(paste("Load failed:", e$message), type = "error")
          })
        })
      }
    })
  })
}