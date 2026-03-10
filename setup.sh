#!/bin/bash
set -e

echo "Cloning repository into 'phytosplice'..."
# Git clone will automatically create the folder if it doesn't exist 
# and put the repo contents directly into it.
git clone --depth 1 https://github.com/jkt-lab/PhytoSplice.git phytosplice

echo "Changing directory to 'phytosplice'..."
cd phytosplice

echo "Creating local conda environment in './env'..."
# Creating the environment in the ./env directory as per project guidelines
conda env create -p ./env -f ../environment.yml 

echo "Making install_desktop_shortcut.sh executable and running it..."
chmod +x install_desktop_shortcut.sh
./install_desktop_shortcut.sh

echo "Setup complete!"
