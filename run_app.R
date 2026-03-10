library(shiny)

# Set host to 0.0.0.0 to allow remote access
options(shiny.host = '0.0.0.0')
options(shiny.port = 3838)

# Run the app from the current directory
runApp('.', launch.browser = FALSE)