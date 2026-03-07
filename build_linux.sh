#!/bin/bash
# =============================================================================
# LOREON — Linux / Ubuntu Build Script
# Produces: dist/LOREON/          (executable folder)
#           dist/LOREON_linux.tar.gz  (distributable archive)
#
# Requirements (on the build machine):
#   - Python 3.8+  with pip
#   - PyQt5 system libraries:
#       sudo apt install python3-pyqt5 libxcb-xinerama0 libxcb-cursor0
#   - All project Python dependencies:
#       pip install -r requirements.txt
#
# Usage:
#   chmod +x build_linux.sh
#   ./build_linux.sh
#
# End-user installation:
#   tar -xzf LOREON_linux.tar.gz
#   cd LOREON
#   ./LOREON                      # launch directly
#   # OR for desktop integration:
#   ./install_desktop.sh
# =============================================================================

set -e
cd "$(dirname "$0")"

# ── Helpers ──────────────────────────────────────────────────────────────────
info()  { echo "[INFO]  $*"; }
error() { echo "[ERROR] $*" >&2; exit 1; }

# ── Checks ───────────────────────────────────────────────────────────────────
[ -f "pipeline_gui.py" ]        || error "pipeline_gui.py not found."
[ -f "loreon_app_icon_2.png" ]  || error "loreon_app_icon_2.png not found."
[ -f "loreon.jpeg" ]            || error "loreon.jpeg not found."
[ -f "credits.md" ]             || error "credits.md not found."

# ── Install build tools ───────────────────────────────────────────────────────
info "Installing PyInstaller..."
python3 -m pip install pyinstaller --quiet

# ── PyInstaller build ─────────────────────────────────────────────────────────
info "Running PyInstaller..."
python3 -m PyInstaller \
    --onedir \
    --windowed \
    --name LOREON \
    --hidden-import=_socket \
    --hidden-import=select \
    --hidden-import=PyQt5.sip \
    --hidden-import=encodings.utf_8 \
    --hidden-import=encodings.ascii \
    --collect-all=PyQt5 \
    --clean \
    -y \
    pipeline_gui.py

# ── Copy data files ───────────────────────────────────────────────────────────
DEST="dist/LOREON"
info "Copying data files into $DEST ..."

cp Dockerfile              "$DEST/"
cp requirements.txt        "$DEST/"
cp loreon_app_icon_2.png   "$DEST/"
cp loreon.jpeg             "$DEST/"
cp credits.md              "$DEST/"
cp metaGenomics_new.py     "$DEST/"
cp OtuUtils.py             "$DEST/"
cp ResultsReader.py        "$DEST/"
cp PipelineLogger.py       "$DEST/"
cp report_generator.py     "$DEST/"
cp pipeline_profiler.py    "$DEST/"
cp template.html           "$DEST/"
[ -f .dockerignore ] && cp .dockerignore "$DEST/"

# ── Launcher shell script ─────────────────────────────────────────────────────
info "Creating launcher script..."
cat > "$DEST/loreon.sh" << 'LAUNCHER'
#!/bin/bash
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$DIR/LOREON" "$@"
LAUNCHER
chmod +x "$DEST/loreon.sh"

# ── Desktop integration helper ────────────────────────────────────────────────
info "Creating desktop integration script..."
cat > "$DEST/install_desktop.sh" << 'DESKTOP_INSTALLER'
#!/bin/bash
# Installs a .desktop entry so LOREON appears in the application menu.
INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DESKTOP_FILE="$HOME/.local/share/applications/loreon.desktop"
mkdir -p "$(dirname "$DESKTOP_FILE")"
cat > "$DESKTOP_FILE" << EOF
[Desktop Entry]
Name=LOREON
Comment=Long-Read ONT Metagenomic Pipeline
Exec=$INSTALL_DIR/LOREON
Icon=$INSTALL_DIR/loreon_app_icon_2.png
Terminal=false
Type=Application
Categories=Science;Biology;Education;
StartupWMClass=LOREON
EOF
chmod +x "$DESKTOP_FILE"
update-desktop-database "$HOME/.local/share/applications/" 2>/dev/null || true
echo "Desktop entry installed: $DESKTOP_FILE"
echo "LOREON should now appear in your application menu."
DESKTOP_INSTALLER
chmod +x "$DEST/install_desktop.sh"

# ── Create tar.gz archive ─────────────────────────────────────────────────────
info "Creating LOREON_linux.tar.gz ..."
rm -f dist/LOREON_linux.tar.gz
tar -czf dist/LOREON_linux.tar.gz -C dist LOREON

echo ""
echo "══════════════════════════════════════════"
echo "  BUILD COMPLETE"
echo "  Executable folder : dist/LOREON/"
echo "  Archive           : dist/LOREON_linux.tar.gz"
echo "══════════════════════════════════════════"
echo ""
echo "Distribution instructions:"
echo "  1. Copy LOREON_linux.tar.gz to the target machine"
echo "  2. tar -xzf LOREON_linux.tar.gz"
echo "  3. cd LOREON && ./LOREON"
echo "  4. For app-menu entry: ./install_desktop.sh"
echo ""
echo "If the app fails to start, install Qt system libs:"
echo "  sudo apt install libxcb-xinerama0 libxcb-cursor0 libgl1"
