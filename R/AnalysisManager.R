library(R6)

# AnalysisManager Class
# Manages the global state of the application with robust persistence.
#
# Workspace Structure (Directory based):
# [ProjectName]/
#   ├── manifest.json       (Metadata, Versioning, Integrity Info)
#   ├── state.rds           (Serialized R state: Params, UI settings)
#   └── artifacts/          (Heavy data, Logs, Results)
#       ├── rmats/
#       └── logs/

AnalysisManager <- R6Class("AnalysisManager",
  public = list(
    # Fields
    SCHEMA_VERSION = "1.0.0",
    modules = list(),
    
    state = list(
      data = list(),
      params = list(),
      status = list()
    ),
    
    workspace_dir = NULL,
    is_loaded = FALSE, 
    events = NULL,
    last_saved = NULL,
    
    # Constructor
    initialize = function() {
      self$modules <- list()
      self$state$status <- list()
      self$is_loaded <- FALSE
      self$events <- shiny::reactiveValues(workspace_loaded = 0, status_changed = 0)
      self$last_saved <- "Never"
    },
    
    # --- Module Management ---
    register_module = function(id, name, ui_func, server_func, icon_name) {
      self$modules[[id]] <- list(id=id, name=name, ui=ui_func, server=server_func, icon=icon_name)
      if (is.null(self$state$status[[id]])) self$state$status[[id]] <- "pending"
    },
    
    # --- State Management ---
    set_result = function(module_id, result_data, parameters = NULL) {
      self$state$data[[module_id]] <- result_data
      if (!is.null(parameters)) self$state$params[[module_id]] <- parameters
      self$state$status[[module_id]] <- "completed"
      self$events$status_changed <- self$events$status_changed + 1
    },
    
    get_result = function(module_id) self$state$data[[module_id]],
    get_params = function(module_id) self$state$params[[module_id]],
    get_status = function(module_id) self$state$status[[module_id]],
    
    # --- Workspace Operations ---
    
    create_workspace = function(parent_dir, name) {
      target_dir <- file.path(parent_dir, name)
      if (dir.exists(target_dir)) stop("Workspace directory already exists.")
      
      dir.create(target_dir, recursive = TRUE)
      dir.create(file.path(target_dir, "artifacts"))
      
      # Initialize empty state
      self$workspace_dir <- target_dir
      self$is_loaded <- TRUE
      
      # Initial Save to establish structure
      self$save_workspace()
      
      self$events$workspace_loaded <- self$events$workspace_loaded + 1
      return(target_dir)
    },
    
    save_workspace = function() {
      if (!self$is_loaded || is.null(self$workspace_dir)) stop("No active workspace.")
      
      # Manifest
      self$last_saved <- as.character(Sys.time())
      manifest <- list(
        version = self$SCHEMA_VERSION,
        updated_at = self$last_saved
      )
      jsonlite::write_json(manifest, file.path(self$workspace_dir, "manifest.json"), auto_unbox = TRUE)
      
      # State
      saveRDS(self$state, file.path(self$workspace_dir, "state.rds"))
      
      message(paste("Workspace saved to:", self$workspace_dir))
      return(self$workspace_dir)
    },
    
    load_workspace = function(target_dir) {
      if (!dir.exists(target_dir)) stop("Directory not found.")
      if (!file.exists(file.path(target_dir, "manifest.json"))) stop("Not a valid workspace.")
      
      manifest <- jsonlite::read_json(file.path(target_dir, "manifest.json"))
      if (!is.null(manifest$updated_at)) self$last_saved <- manifest$updated_at
      
      self$state <- readRDS(file.path(target_dir, "state.rds"))
      
      # Restore paths (Basic rehydration)
      if (!is.null(self$state$data$rmats$output_dir)) {
         self$state$data$rmats$output_dir <- file.path(target_dir, "artifacts/rmats/output")
      }
      if (!is.null(self$state$data$suppa$output_dir)) {
         self$state$data$suppa$output_dir <- file.path(target_dir, "artifacts/suppa")
      }
      
      self$workspace_dir <- target_dir
      self$is_loaded <- TRUE
      
      self$events$workspace_loaded <- self$events$workspace_loaded + 1
      return(TRUE)
    }
  )
)
