# PhytoSplice

PhytoSplice is a comprehensive, modular Shiny application engineered to streamline Alternative Splicing (AS) analysis. Designed with a user-friendly interface, it unifies the entire computational workflow from read alignment and transcript quantification to differential splicing analysis and interactive visualization.

## 🚀 Key Features

PhytoSplice integrates industry-standard bioinformatics tools into a seamless graphical interface:

| Module | Core Tool | Capabilities |
| --- | --- | --- |
| **Workspace Manager** | PhytoSplice | Create, load, and save projects with full state persistence. Save analysis parameters and paths to easily resume work later. |
| **Read Alignment** | HISAT2 | Download prebuilt genome indexes for common plant species or build custom splice-aware indexes from FASTA and GTF files. Execute alignments in batch. |
| **Quantification** | Salmon | Build decoy-aware indexes and perform batch quantification. Automatically merge TPM values for downstream SUPPA analysis. |
| **Differential Splicing** | rMATS (v4.3.0+) | Configure and run analysis for Paired-End/Single-End reads with real-time execution logging. Filter results via interactive tables using FDR and DeltaPSI. |
| **Alternative Splicing** | SUPPA2 | Complete SUPPA2 workflow: Generate Events, Calculate PSI, and perform Differential Splicing analysis. |
| **Visualization** | Integrated | Interactive Volcano Plots, Heatmaps, and PCA plots based on PSI values. Integrated Sashimi plots (wrapping `rmats2sashimiplot`) and PSI distribution plots (Violin/Box). |

---

## 🛠️ Prerequisites & Installation

PhytoSplice is designed to run within a dedicated Conda environment to ensure seamless dependency management.

### 1. Conda Environment Setup (Recommended)

An `environment.yml` file is provided to automatically install all necessary R packages and external command-line tools.

```bash
conda env create -f environment.yml
conda activate ./env

```

> **Note on External Tools:** The Conda environment manages dependencies like `rmats.py` (v4.3.0+), `suppa.py`, `rmats2sashimiplot`, `HISAT2`, `Salmon`, `Samtools`, `gffread`, and `gdown`. If you choose to install these manually, they must be accessible in your system's `PATH`.

### 2. Manual R Package Installation (Alternative)

If you prefer to manage R packages manually outside of Conda, install the core dependencies:

```r
install.packages(c("shiny", "shinydashboard", "shinyjs", "shinyWidgets", "colourpicker", "bslib", "plotly", "tidyverse", "DT", "R6", "processx"))

```

*(Note: You must still ensure the external tools listed above are installed and in your `PATH`)*

---

## 💻 Usage

* **Step 1: Clone the Repository**
```bash
git clone <repository-url>
cd SplicingShinySuite

```


*(Note: Ensure your directory matches the repository name)*
* **Step 2: Launch the Application**
You can start the app via the command line or set it up as a desktop application:
```bash
# Option 1: Command Line Interface (CLI)
./run.sh

# Option 2: Install Desktop Shortcut (Linux only)
./install_desktop_shortcut.sh
# Launch 'PhytoSplice' directly from your application menu.

```


* **Step 3: Access the Interface**
Open your web browser and navigate to the address shown in your console (default is `http://0.0.0.0:3838`).
* **Step 4: Begin Your Analysis**
Click **Create New Workspace** to initialize a project folder, then navigate through the top tabs (HISAT2, Salmon, rMATS, SUPPA) to progress through your workflow.

---

## 📂 Project Architecture

The application is structured to ensure modularity and maintainability:

* **Core Application**
* `app.R`: Main Shiny application entry point defining the UI and server logic.
* `global.R`: Loads libraries, sources helper scripts, and initializes the central `AnalysisManager`.
* `run_app.R`: Dedicated script to launch the app on host `0.0.0.0` and port `3838`.


* **Modules (`modules/`)**
* Contains standalone UI and server logic for each dashboard tab: `mod_home.R`, `mod_rmats.R`, `mod_event_based.R`, `mod_suppa.R`, `mod_hisat2.R`, and `mod_salmon.R`.


* **Logic & Helpers (`R/`)**
* `AnalysisManager.R`: Central R6 class managing state persistence, logging, and data sharing between modules.
* Helper scripts for specific tools and functions: `sashimi_helpers.R`, `suppa_helpers.R`, `salmon_helpers.R`, `annotation_helpers.R`, and `visualization_helpers.R`.


* **Execution Scripts**
* `run.sh`: Primary launch script that activates the Conda environment before starting R.
* `splash.py`: Python/Tkinter graphical splash screen displayed while R loads.
* `install_desktop_shortcut.sh`: Utility to create a Linux `.desktop` file.
