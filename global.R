# Lazy-loaded in app.R server
library(DT)
library(shinyjs)
library(shinyWidgets)
library(bslib)
library(colourpicker)
library(R6)
library(shiny)
library(shinydashboard)

# Source core logic
source("R/AnalysisManager.R")
source("R/sashimi_helpers.R")
source("R/suppa_helpers.R")
source("R/annotation_helpers.R")
source("R/visualization_helpers.R")

# Initialize the global manager instance
# In a real multi-user Shiny Server/Connect env, this should be inside server function 
# or handled carefully. For single-process/local use, global is fine, 
# BUT to support multiple sessions with independent state, we should instantiate it inside server.
# However, for module registration (which is static), we can use a global registry pattern 
# or just register them in app.R.
# Let's keep the *Definition* global, but instantiation inside server for session isolation.

# Source all modules
module_files <- list.files('modules', full.names = TRUE, pattern = "\\.R$")
if (length(module_files) > 0) {
  sapply(module_files, source)
}