#!/bin/bash

# --- CONFIGURATION ---
# Set to 'true' to show a terminal window with logs and keep it open.
# Set to 'false' to run in the background (headless), logging to ~/.PhytoSplice.log.
SHOW_TERMINAL=false
# ---------------------

# Get the absolute path of the current project directory
PROJECT_DIR="$(pwd)"
ICON_NAME="$PROJECT_DIR/logo.png"

echo "Setting up PhytoSplice desktop shortcut..."
echo "Terminal Mode: $SHOW_TERMINAL"

# 1. Start creating run_desktop.sh (Common Header)
cat > run_desktop.sh << 'EOF'
#!/bin/bash

# Define the path to the Conda R executable
R_EXEC="./env/bin/Rscript"
PY_EXEC="./env/bin/python"
# We need to know where we are if run from absolute path
PROJECT_DIR="$(dirname "$(realpath "$0")")"
cd "$PROJECT_DIR"

# Check if R exists
if [ ! -f "$R_EXEC" ]; then
    echo "Error: R executable not found at $R_EXEC"
    echo "Please ensure the conda environment is in 'env/' directory."
    echo "Press Enter to exit..."
    read
    exit 1
fi

echo "Starting PhytoSplice..."
echo "Using R at: $R_EXEC"

# Cleanup: Kill any existing instance on port 3838
if command -v fuser >/dev/null 2>&1; then
    fuser -k 3838/tcp >/dev/null 2>&1
fi

# Set default browser for R (fallback)
export R_BROWSER=xdg-open

# Helper function to launch browser
launch_browser() {
    local APP_URL="$1"
    echo "Launching browser for $APP_URL..."

    # Create a temporary user data dir to ensure the browser process blocks/isolates
    local BROWSER_DATA_DIR=$(mktemp -d)

    # Use --class or --name to match the StartupWMClass in .desktop file
    if command -v microsoft-edge >/dev/null 2>&1; then
        microsoft-edge --app="$APP_URL" --user-data-dir="$BROWSER_DATA_DIR" --class="PhytoSplice" --name="PhytoSplice"
    elif command -v google-chrome >/dev/null 2>&1; then
        google-chrome --app="$APP_URL" --user-data-dir="$BROWSER_DATA_DIR" --class="PhytoSplice" --name="PhytoSplice"
    elif command -v chromium-browser >/dev/null 2>&1; then
        chromium-browser --app="$APP_URL" --user-data-dir="$BROWSER_DATA_DIR" --class="PhytoSplice" --name="PhytoSplice"
    elif command -v chromium >/dev/null 2>&1; then
        chromium --app="$APP_URL" --user-data-dir="$BROWSER_DATA_DIR" --class="PhytoSplice" --name="PhytoSplice"
    elif command -v firefox >/dev/null 2>&1; then
        # Firefox doesn't have a strict --app mode, but --new-window is close
        # --no-remote ensures a separate process so we can wait for it
        # Firefox ignores --class often, but we try anyway.
        firefox --new-window "$APP_URL" --profile "$BROWSER_DATA_DIR" --no-remote --class "PhytoSplice" --name "PhytoSplice"
    else
        xdg-open "$APP_URL"
    fi
    
    # Cleanup profile
    rm -rf "$BROWSER_DATA_DIR"
}

# --- START SPLASH SCREEN ---
echo "Launching Splash Screen..."
"$PY_EXEC" splash.py > "$HOME/.splash_debug.log" 2>&1 &
SPLASH_PID=$!
# Ensure splash is killed if script exits early
trap "kill $SPLASH_PID 2>/dev/null" EXIT
# ---------------------------
EOF

# 2. Append Mode-Specific Logic to run_desktop.sh
if [ "$SHOW_TERMINAL" = "true" ]; then
    # --- TERMINAL MODE ---
    # Output to stdout, Wait for user Enter to close
    cat >> run_desktop.sh << 'EOF'

echo "Starting R Server (Foreground Logs)..."
"$R_EXEC" -e "library(shiny); options(shiny.host='0.0.0.0', shiny.port=3838); runApp('.', launch.browser=FALSE)" &
SERVER_PID=$!

echo "Waiting for server to initialize..."
# Dynamic check: Wait for port 3838 to be open (TCP connect)
MAX_RETRIES=30
SERVER_READY=false
for i in $(seq 1 $MAX_RETRIES); do
    if "$PY_EXEC" -c "import socket; s = socket.socket(); s.settimeout(1); print('OPEN' if s.connect_ex(('127.0.0.1', 3838)) == 0 else 'CLOSED'); s.close()" 2>/dev/null | grep -q "OPEN"; then
        echo "Server port is open!"
        SERVER_READY=true
        break
    fi
    sleep 1
done

if [ "$SERVER_READY" = "false" ]; then
    echo "Warning: Server port did not open within $MAX_RETRIES seconds."
fi

# --- KILL SPLASH SCREEN ---
kill $SPLASH_PID 2>/dev/null
trap - EXIT
# --------------------------

# Launch Browser in BACKGROUND so the terminal script doesn't block on it
APP_URL="http://127.0.0.1:3838"
launch_browser "$APP_URL" &

echo ""
echo "=========================================================="
echo " PhytoSplice is running."
echo " Server PID: $SERVER_PID"
echo " Logs are being displayed above."
echo ""
echo " -> Press [ENTER] to stop the server and close this window."
echo "=========================================================="
read dummy

echo "Stopping server..."
kill $SERVER_PID
EOF

else
    # --- HEADLESS MODE ---
    # Redirect logs, Wait for browser to close
    cat >> run_desktop.sh << 'EOF'

echo "Starting R Server (Background Logs -> ~/.PhytoSplice.log)..."
"$R_EXEC" -e "library(shiny); options(shiny.host='0.0.0.0', shiny.port=3838); runApp('.', launch.browser=FALSE)" > "$HOME/.PhytoSplice.log" 2>&1 &
SERVER_PID=$!

# Dynamic check: Wait for port 3838 to be open (TCP connect)
MAX_RETRIES=30
SERVER_READY=false
for i in $(seq 1 $MAX_RETRIES); do
    if "$PY_EXEC" -c "import socket; s = socket.socket(); s.settimeout(1); print('OPEN' if s.connect_ex(('127.0.0.1', 3838)) == 0 else 'CLOSED'); s.close()" 2>/dev/null | grep -q "OPEN"; then
        echo "Server port is open!"
        SERVER_READY=true
        break
    fi
    sleep 1
done

if [ "$SERVER_READY" = "false" ]; then
    echo "Warning: Server port did not open within $MAX_RETRIES seconds."
fi

# --- KILL SPLASH SCREEN ---
kill $SPLASH_PID 2>/dev/null
trap - EXIT
# --------------------------

# Launch Browser in FOREGROUND (Blocking)
APP_URL="http://127.0.0.1:3838"
launch_browser "$APP_URL"

echo "Browser closed. Stopping server..."
kill $SERVER_PID
EOF

fi

# Make the wrapper executable
chmod +x run_desktop.sh

# 3. Determine Exec Command and Terminal Flag for .desktop file
if [ "$SHOW_TERMINAL" = "true" ]; then
    # Try to find a known terminal emulator to launch explicitly
    if command -v gnome-terminal >/dev/null 2>&1; then
        # Explicitly launch gnome-terminal running our script
        # Note: We assume the script is executable
        EXEC_CMD="gnome-terminal -- $PROJECT_DIR/run_desktop.sh"
        USE_TERMINAL_FLAG="false" 
    elif command -v x-terminal-emulator >/dev/null 2>&1; then
        EXEC_CMD="x-terminal-emulator -e $PROJECT_DIR/run_desktop.sh"
        USE_TERMINAL_FLAG="false"
    elif command -v konsole >/dev/null 2>&1; then
        EXEC_CMD="konsole -e $PROJECT_DIR/run_desktop.sh"
        USE_TERMINAL_FLAG="false"
    elif command -v xfce4-terminal >/dev/null 2>&1; then
        EXEC_CMD="xfce4-terminal -e $PROJECT_DIR/run_desktop.sh"
        USE_TERMINAL_FLAG="false"
    elif command -v xterm >/dev/null 2>&1; then
        EXEC_CMD="xterm -e $PROJECT_DIR/run_desktop.sh"
        USE_TERMINAL_FLAG="false"
    else
        # Fallback: Let the desktop environment try
        EXEC_CMD="/usr/bin/env bash $PROJECT_DIR/run_desktop.sh"
        USE_TERMINAL_FLAG="true"
    fi
else
    # Headless mode
    EXEC_CMD="/usr/bin/env bash $PROJECT_DIR/run_desktop.sh"
    USE_TERMINAL_FLAG="false"
fi

# 4. Create the .desktop entry
DESKTOP_FILE="$HOME/.local/share/applications/PhytoSplice.desktop"

cat > "$DESKTOP_FILE" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=PhytoSplice
Comment=Run PhytoSplice Analysis Suite
Exec=$EXEC_CMD
Icon=$ICON_NAME
Path=$PROJECT_DIR
Terminal=$USE_TERMINAL_FLAG
StartupWMClass=PhytoSplice
Categories=Science;Education;
EOF

# Make the .desktop file executable (trusted)
chmod +x "$DESKTOP_FILE"

echo "-------------------------------------------------------"
echo "Success! The shortcut 'PhytoSplice' has been added."
echo "You should now find it in your Applications/Dash."
if [ "$SHOW_TERMINAL" = "true" ]; then
    echo "Mode: Terminal (Window should open)"
    echo "Exec Command: $EXEC_CMD"
else
    echo "Mode: Headless (Logs in ~/.PhytoSplice.log)"
fi
echo "-------------------------------------------------------"