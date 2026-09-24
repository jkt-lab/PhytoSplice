# PhytoSplice

  A Hybrid Analytical and Machine Learning Platform for Plant Alternative Splicing & NMD Prediction
  PhytoSplice is a comprehensive, modular Shiny application engineered to streamline Alternative Splicing (AS) analysis. Going beyond standard AS detection, PhytoSplice bridges the gap between transcript quantification and functional transcript
  fate by introducing a rigorous, leakage-audited Machine Learning engine to predict Nonsense-Mediated Decay (NMD) susceptibility.
  Designed with an interactive interface, it unifies the entire computational workflow: from short/long-read alignment and assembly, to differential splicing analysis, ML-based NMD prediction, and interactive visual analytics.

  ## Key Features & Modules

  PhytoSplice integrates industry-standard bioinformatics tools and bespoke machine learning classifiers into a single, seamless graphical interface anchored to canonical gene locus identifiers.
  
   Module                              | Core Tool(s)                        | Capabilities
  -------------------------------------|-------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
   Workspace Manager                   | PhytoSplice                         | Create, load, and save projects with full state persistence. Save analysis parameters and directory paths to easily resume work later.
   Read Alignment                      | HISAT2 / Minimap2                   | Perform splice-aware short-read alignment (HISAT2) or long-read Iso-Seq alignment (minimap2 -ax splice).
   Assembly                            | StringTie2                          | Reconstruct transcriptomes using reference-guided or long-read modes (-L, -G, --merge).
   Quantification                      | Salmon                              | Build decoy-aware indexes and perform batch alignment-free quantification. Automatically merges TPM values for downstream SUPPA2 analysis.
   Differential Splicing               | rMATS / SUPPA2                      | Execute robust statistical AS event detection for multiple conditions. Filter results via interactive tables using FDR and DeltaPSI.
   NMD Prediction (v2.0)               | XGBoost / Python                    | NEW: Classify transcripts as NMD Targets or Escapes. Features dual-track prediction (ML model + Rule-based structural heuristic) and real-time SHAP interpretability.
   Visualization                       | Plotly / ggsashimi                  | Interactive Volcano Plots, Heatmaps, PCA, and PSI distribution curves. Integrated in-browser Sashimi plot generator.

  ## The NMD Prediction Engine
  
  The newest addition to PhytoSplice is a highly calibrated, sequence-based NMD prediction module designed with strict biological rigor:

  • Dual-Track Inference: Choose between our pre-trained XGBoost models (validated on Arabidopsis and Rice) or a species-agnostic rule-based heuristic (PTC–EJC distance, 3′ UTR length, downstream EJCs, uORFs) for highly divergent genomes.
  • Leakage-Audited: Built to avoid the "tautological leakage" common in genomic ML. Label-defining structural features are rigorously pruned, relying strictly on an independent matrix of 19 sequence-composition features (k-mers, GC content,
  uORFs, intron length).
  • Actionable Precision: Thresholds are widened (≥0.75 for Targets, ≤0.25 for Escapes) to prioritize high-confidence bench-actionable discoveries rather than inflated headline accuracy.
  • SHAP Interpretability: "Black-box" predictions are a thing of the past. The dashboard instantly generates SHAP feature attributions, explaining exactly which sequence motifs drove the NMD prediction for every single isoform.
  ──────
  ## 🛠️ Prerequisites & Installation

  PhytoSplice relies on both R and Python environments to orchestrate the UI and the deep learning backend. It is designed to run within a dedicated Conda environment to ensure seamless dependency management.

  ### 1. Conda Environment Setup (Recommended)

  An environment.yml file is provided to automatically install all necessary R packages, Python ML libraries (xgboost, shap, scikit-learn), and external command-line tools.

    conda env create -f environment.yml
    conda activate ./env

  │ Note on External Tools: The Conda environment automatically manages dependencies like rmats.py, suppa.py, HISAT2, minimap2, StringTie2, Salmon, and Samtools. If you choose to install these manually, they must be accessible in your system's
  │ PATH.

  ### 2. Manual R & Python Installation (Alternative)

  If you prefer to manage packages manually outside of Conda:

  • R Dependencies:
    install.packages(c("shiny", "shinydashboard", "shinyjs", "shinyWidgets", "colourpicker", "bslib", "plotly", "tidyverse", "DT", "R6", "processx"))
    BiocManager::install(c("Biostrings", "GenomicRanges"))

  • Python Dependencies (Python 3.12+):
    pip install xgboost shap scikit-learn numpy pandas matplotlib

  
  ## Usage

  Step 1: Clone the Repository

    git clone <repository-url>
    cd shinysplicingsuite

  Step 2: Launch the Application
  You can start the app via the command line or set it up as a desktop application:

    # Option 1: Command Line Interface (CLI)
    ./run.sh

    # Option 2: Install Desktop Shortcut (Linux only)
    ./install_desktop_shortcut.sh
    # Launch 'PhytoSplice' directly from your application menu.

  Step 3: Access the Interface
  Open your web browser and navigate to the address shown in your console (default is http://0.0.0.0:3838).

  Step 4: Begin Your Analysis
  Click Create New Workspace to initialize a project folder, then navigate sequentially through the top tabs (HISAT2, StringTie, Salmon, AS Detection, NMD Prediction) to progress through your workflow.
  
  ## Project Architecture

  The application is structured to ensure modularity and maintainability:

  • Core Application
      • app.R: Main Shiny application entry point defining the UI and server logic.
      • global.R: Loads libraries, sources helper scripts, and initializes the central AnalysisManager.
      • run_app.R: Dedicated script to launch the app on host 0.0.0.0 and port 3838.
  • Modules (modules/)
      • Standalone UI and server logic for each dashboard tab: mod_home.R, mod_rmats.R, mod_suppa.R, mod_hisat2.R, mod_salmon.R, mod_stringtie.R, and the new mod_nmd.R.
  • Logic & Helpers (R/)
      • AnalysisManager.R: Central R6 class managing state persistence, logging, and data sharing between modules.
      • Helper scripts for specific tools: sashimi_helpers.R, suppa_helpers.R, annotation_helpers.R, etc.
  • Python Backend (python/)
      • predict_nmd.py & extract_features.py: Handles sequence parsing, XGBoost execution, and SHAP value calculation, communicating securely with the R frontend via JSON.
  • Execution Scripts
      • run.sh: Primary launch script that activates the Conda environment before starting R.
      • splash.py: Python/Tkinter graphical splash screen displayed while R loads.
      • install_desktop_shortcut.sh: Utility to create a Linux .desktop file.
