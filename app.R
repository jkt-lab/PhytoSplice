source("global.R")
library(tcltk)

ui <- dashboardPage(
  dashboardHeader(
    title = "PhytoSplice",
    tags$li(class = "dropdown", 
      actionButton("header_save", "Save Workspace", icon = icon("save"), class = "btn-success", style = "margin-top: 8px; margin-right: 10px;")
    )
  ),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("HISAT2 Alignment", tabName = "hisat2", icon = icon("align-center")),
      menuItem("Salmon Quantification", tabName = "salmon", icon = icon("fish")),
      menuItem("rMATS Analysis", tabName = "rmats", icon = icon("terminal")),
      menuItem("SUPPA Analysis", tabName = "suppa", icon = icon("dna"))
    )
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
        /* Loading Screen */
        #loading_screen {
          position: fixed;
          top: 0;
          left: 0;
          width: 100%;
          height: 100%;
          background: #222d32;
          z-index: 10000;
          display: flex;
          flex-direction: column;
          justify-content: center;
          align-items: center;
          color: white;
          font-family: 'Source Sans Pro', sans-serif;
        }
        .spinner {
          width: 50px;
          height: 50px;
          border: 5px solid #ecf0f5;
          border-top: 5px solid #3c8dbc;
          border-radius: 50%;
          animation: spin 1s linear infinite;
          margin-bottom: 20px;
        }
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }

        .main-header {
          position: fixed;
          width: 100%;
          top: 0;
          z-index: 1000;
        }
        .main-sidebar {
          position: fixed !important;
          height: calc(100vh - 50px) !important;
          top: 50px !important;
          overflow-y: auto;
        }
        .content-wrapper, .right-side {
          margin-top: 50px !important;
        }
        /* Match header to sidebar background */
        .skin-blue .main-header .logo {
          background-color: #222d32;
        }
        .skin-blue .main-header .logo:hover {
          background-color: #1e282c;
        }
        .skin-blue .main-header .navbar {
          background-color: #222d32;
        }
        /* Reposition Notifications to Bottom-Left (Inside Sidebar) */
        #shiny-notification-panel {
          position: fixed !important;
          bottom: 10px !important;
          left: 5px !important;
          right: auto !important;
          top: auto !important;
          width: 220px !important; /* Fits within 230px sidebar */
          max-height: 50vh !important; /* Max height 50% of screen */
          overflow-y: auto !important; /* Scroll if too many */
          z-index: 99999 !important;
        }
        
        /* Dark Grey Theme for Notifications */
        .shiny-notification {
          background-color: #2c3b41 !important; /* Dark grey matching sidebar */
          color: #ffffff !important;
          border: 1px solid #4b646f !important;
          border-radius: 3px !important;
          opacity: 0.95 !important;
          box-shadow: 0 1px 1px rgba(0,0,0,0.1) !important;
        }
        
        /* White Close Button */
        .shiny-notification-close {
          color: #ffffff !important;
          opacity: 0.6 !important;
        }
        .shiny-notification-close:hover {
          opacity: 1 !important;
        }
      ")),
      tags$script(HTML("
        $(document).ready(function(){
          $('body').tooltip({
            selector: '[data-toggle=\"tooltip\"]',
            container: 'body',
            trigger: 'hover',
            html: true
          });
        });
      "))
    ),
    useShinyjs(),
    div(id = "loading_screen", 
      div(class = "spinner"),
      h2("PhytoSplice"),
      p("Initializing environment...")
    ),
    uiOutput("body_content") 
  )
)

server <- function(input, output, session) {
  
  # --- Lazy Library Loading ---
  # Load heavy libraries here so they don't block the initial UI flush
  library(jsonlite)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(stringr)
  library(httr)
  library(plotly)
  
  try({ library(GenomicRanges) }, silent=TRUE)

  manager <- AnalysisManager$new()
  rv <- reactiveValues(workspace_loaded = FALSE)
  
  # Register Modules
  manager$register_module("home", "Home", mod_home_ui, mod_home_server, "home")
  manager$register_module("hisat2", "HISAT2 Alignment", mod_hisat2_ui, mod_hisat2_server, "align-center")
  manager$register_module("salmon", "Salmon Quantification", mod_salmon_ui, mod_salmon_server, "fish")
  manager$register_module("rmats", "rMATS Analysis", mod_rmats_ui, mod_rmats_server, "terminal")
  manager$register_module("suppa", "SUPPA Analysis", mod_suppa_ui, mod_suppa_server, "dna")
  
  # Initialize Servers
  for (mod_id in names(manager$modules)) {
    mod <- manager$modules[[mod_id]]
    mod$server(mod$id, manager)
  }

  # --- Remove Loading Screen & Show Modal ---
  observe({
    shinyjs::delay(100, {
      shinyjs::hide("loading_screen", anim = TRUE, animType = "fade", time = 0.5)
      
      showModal(modalDialog(
        title = "Welcome to PhytoSplice",
        size = "m",
        footer = NULL,
        easyClose = FALSE,
        fade = FALSE,
        div(style = "text-align: center;",
          h4("Workspace Setup"),
          p("A workspace is required to organize analysis files and results."),
          br(),
          actionButton("start_create", "Create New Workspace", class = "btn-primary btn-lg", icon = icon("plus")),
          br(), br(),
          actionButton("start_load", "Load Existing Workspace", class = "btn-warning btn-lg", icon = icon("folder-open"))
        )
      ))
    })
  })
  
  # --- UI Rendering ---
  output$body_content <- renderUI({
    # If not loaded, show nothing (Modal is active) or a placeholder
    if (!rv$workspace_loaded) {
      return(div(style="padding: 20px; text-align: center; color: #777;", h2("Please Create or Load a Workspace to begin.")))
    }
    
    req(input$tabs)
    mod <- manager$modules[[input$tabs]]
    if (!is.null(mod)) {
      mod$ui(mod$id)
    } else {
      h3("Module not found")
    }
  })
  
  # --- Startup Modal ---
  
  # Native Ubuntu File Picker Helper
  get_native_dir <- function(title = "Select Directory") {
    # check for zenity
    if (nzchar(Sys.which("zenity"))) {
      tryCatch({
        path <- system(paste0('zenity --file-selection --directory --title="', title, '"'), intern = TRUE, ignore.stderr = TRUE)
        if (length(path) > 0 && nzchar(path)) return(path)
      }, error = function(e) NA)
    } else if (capabilities("tcltk") && !is.null(Sys.getenv("DISPLAY")) && Sys.getenv("DISPLAY") != "") {
      # Fallback
      return(tcltk::tk_choose.dir(default = getwd(), caption = title))
    }
    return(NA)
  }
  
  # Old Modal Call Removed (Moved to Loading Sequence)
  
  # --- Workspace Handlers ---
  
  observeEvent(input$start_create, {
    removeModal()
    showModal(modalDialog(
      title = "New Workspace",
      textInput("new_ws_name", "Project Name", value = "NewProject"),
      footer = tagList(
        actionButton("confirm_create", "Create", class = "btn-success"),
        actionButton("cancel_start", "Cancel")
      )
    ))
  })
  
  observeEvent(input$confirm_create, {
    req(input$new_ws_name)
    parent_dir <- getwd()
    
    # Use native picker
    sel <- get_native_dir("Select Parent Directory")
    if (!is.na(sel)) parent_dir <- sel
    
    tryCatch({
      manager$create_workspace(parent_dir, input$new_ws_name)
      removeModal()
      rv$workspace_loaded <- TRUE
      showNotification(paste("Workspace created:", manager$workspace_dir), type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  observeEvent(input$start_load, {
    target_dir <- get_native_dir("Select Workspace Folder")
    if (!is.na(target_dir)) {
      tryCatch({
        manager$load_workspace(target_dir)
        removeModal()
        rv$workspace_loaded <- TRUE
        showNotification("Workspace loaded!", type = "message")
      }, error = function(e) {
        showNotification(paste("Load failed:", e$message), type = "error")
      })
    }
  })
  
  observeEvent(input$cancel_start, { session$reload() })
  
  # --- Global Save ---
  observeEvent(input$header_save, {
    if (manager$is_loaded) {
      manager$save_workspace()
      showNotification("Workspace Saved.", type = "message")
    } else {
      showNotification("No active workspace.", type = "warning")
    }
  })
}

shinyApp(ui, server)