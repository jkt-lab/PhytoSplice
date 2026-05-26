PhytoSplice
A comprehensive, modular Shiny application for Alternative Splicing (AS) analysis. This tool provides a user-friendly interface for the entire workflow: from read alignment (HISAT2) and quantification (Salmon), to differential splicing analysis (rMATS, SUPPA2), and interactive visualization.
Features
rMATS Execution Module:
Configure and run rMATS for differential splicing analysis directly from the UI.
Supports both Paired-End and Single-End reads.
Real-time execution logging.
rMATS Results (Event-Based Analysis):
Integrated Results Section: Event-based analysis is provided inside the rMATS analysis workflow as the results view.
Interactive Visualizations: Volcano Plots, Heatmaps, and PCA plots based on PSI values.
Sashimi Plots: Integrated visualization (wrapping `rmats2sashimiplot`) for inspecting read coverage across splice junctions.
Data Exploration: Interactive tables with FDR and DeltaPSI filtering.
SUPPA Analysis Module:
Workflow for SUPPA2 (Generate Events -> PSI Calculation -> Differential Splicing).
Visualizations: Volcano plots, Heatmaps, and PSI distribution (Violin/Box) plots for significant events.
HISAT2 Alignment Module:
Prebuilt genome index downloader for common plant species.
Build custom splice-aware indexes from FASTA and GTF.
Batch alignment execution.
Salmon Quantification Module:
Decoy-aware index building.
Batch quantification and automated merging of TPM values for downstream SUPPA analysis.
Workspace System:
Project Management: Create, Load, and Save workspaces to organize specific analysis runs.
Persistence: Analysis state, parameters, and paths are saved, allowing you to resume work later.
Prerequisites
The application is designed to run within a Conda environment.
Conda Environment Setup
We provide an `environment.yml` file to set up all necessary dependencies (R packages and command-line tools) in an isolated environment.
```bash
conda env create -f environment.yml
conda activate ./env
```
If you prefer to manage R packages manually, the core dependencies are:
```r
install.packages(c("shiny", "shinydashboard", "shinyjs", "shinyWidgets", "colourpicker", "bslib", "plotly", "tidyverse", "DT", "R6", "processx"))
```
External Tools
The following tools are managed by the conda environment, but if installing manually, they must be accessible in your system `PATH`:
rMATS (v4.3.0+) - `rmats.py`
SUPPA2 - `suppa.py`
rmats2sashimiplot
HISAT2
Salmon
Samtools
gffread
gdown
Usage
Clone the Repository
    ```bash
    git clone <repository-url>
    cd SplicingShinySuite
    ```
Run the Application
You can start the app using the provided shell script or directly via R:
    ```bash
    # Option 1: Desktop Application Mode (Linux)
    ./install_desktop_shortcut.sh
    # You can then launch 'PhytoSplice' from your application menu.
    
    # Option 2: Shell script (CLI)
    ./run.sh
    ```
Access the Interface
Open your web browser and navigate to the address shown in the console (default: `http://0.0.0.0:3838`).
Start Analysis
Click Create New Workspace to initialize a project folder.
Navigate through the tabs (HISAT2, Salmon, rMATS, SUPPA) to perform your analysis.
Project Structure
Core Application
`app.R`: The main entry point for the Shiny application. It defines the `ui` (dashboard layout) and `server` logic, integrating all modules.
`global.R`: Responsible for the global environment setup. It loads necessary R libraries, sources helper scripts from `R/`, and initializes the `AnalysisManager` R6 class.
`run_app.R`: A dedicated R script to launch the application on host `0.0.0.0` and port `3838`, commonly used by the shell launchers.
Modules (`modules/`)
These files contain the modularized logic for each tab in the application.
`mod_home.R`: The landing page module. Manages workspace creation, loading, and status display.
`mod_rmats.R`: User interface and server logic for configuring and executing rMATS.
`mod_event_based.R`: Handles the event-based visualization of rMATS results.
`mod_suppa.R`: Implements the workflow for SUPPA2.
`mod_hisat2.R`: Module for HISAT2 alignment and index management.
`mod_salmon.R`: Module for Salmon quantification and index generation.
Logic & Helpers (`R/`)
`AnalysisManager.R`: A central R6 class that acts as the "State Manager". It handles workspace persistence (saving/loading `state.rds`), logging, and data sharing between modules.
`sashimi_helpers.R`: Functions specifically designed to format data for and execute `rmats2sashimiplot`.
`suppa_helpers.R`: Data parsing, formatting, and filtering functions specific to SUPPA.
`salmon_helpers.R`: Functions for building Salmon indexes and merging quantification outputs.
`annotation_helpers.R`: Functions for mapping gene IDs to gene symbols using GTF files or APIs.
`visualization_helpers.R`: General-purpose plotting utilities used across multiple modules.
Execution
`run.sh`: The primary shell script to start the application. It automatically detects and activates the Conda environment before launching R.
`splash.py`: A Python/Tkinter script that displays a graphical splash screen while the R environment loads.
`install_desktop_shortcut.sh`: Helper script to create a Linux desktop shortcut (`.desktop` file) and launch wrapper.
