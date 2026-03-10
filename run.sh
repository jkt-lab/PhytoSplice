#!/bin/bash

# Initialize Conda for script execution
# We attempt to source the conda.sh script to enable 'conda activate'
if command -v conda >/dev/null 2>&1; then
    CONDA_BASE=$(conda info --base)
    if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        conda activate ./env
    else
        # Fallback: manually add to PATH if conda.sh isn't found
        export PATH="$PWD/env/bin:$PATH"
    fi
else
    # Fallback: manually add to PATH if conda command isn't found
    export PATH="$PWD/env/bin:$PATH"
fi

# Define the R executable (now just Rscript since it's in PATH)
R_EXEC="Rscript"

# Check if R exists in PATH
if ! command -v "$R_EXEC" >/dev/null 2>&1; then
    echo "Error: Rscript not found."
    echo "Please ensure the conda environment is in 'env/' directory."
    exit 1
fi

echo "Starting PhytoSplice..."
echo "Using R at: $R_EXEC"

# Run the app
# We use run_app.R to set specific host/port
"$R_EXEC" run_app.R
