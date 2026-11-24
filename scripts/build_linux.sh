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
#       pip install -r sources/requirements.txt
#
# Usage:
#   chmod +x scripts/build_linux.sh
#   scripts/build_linux.sh
#
# End-user installation:
#   tar -xzf LOREON_linux.tar.gz
#   cd LOREON
#   ./LOREON                      # launch directly
#   # OR for desktop integration:
#   ./install_desktop.sh
# =============================================================================

set -e
cd "$(dirname "$0")/.."

# Folder of the application files, relative to the root of the repository
SRC="sources"

# ── Helpers ──────────────────────────────────────────────────────────────────
info()  { echo "[INFO]  $*"; }
error() { echo "[ERROR] $*" >&2; exit 1; }

# ── Checks ───────────────────────────────────────────────────────────────────
[ -f "$SRC/pipeline_gui.py" ]        || error "$SRC/pipeline_gui.py not found."
[ -f "$SRC/loreon_app_icon_2.png" ]  || error "$SRC/loreon_app_icon_2.png not found."
[ -f "$SRC/loreon.jpeg" ]            || error "$SRC/loreon.jpeg not found."
[ -f "$SRC/credits.md" ]             || error "$SRC/credits.md not found."

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
    --paths "$SRC" \
    --clean \
    -y \
    "$SRC/pipeline_gui.py"

# ── Copy data files ───────────────────────────────────────────────────────────
DEST="dist/LOREON"
info "Copying data files into $DEST ..."

cp "$SRC/Dockerfile"              "$DEST/"
cp "$SRC/requirements.txt"        "$DEST/"
cp "$SRC/loreon_app_icon_2.png"   "$DEST/"
cp "$SRC/loreon.jpeg"             "$DEST/"
cp "$SRC/credits.md"              "$DEST/"
cp LICENSE                 "$DEST/"
cp "$SRC/metaGenomics_new.py"     "$DEST/"
cp "$SRC/OtuUtils.py"             "$DEST/"
cp "$SRC/ResultsReader.py"        "$DEST/"
cp "$SRC/PipelineLogger.py"       "$DEST/"
cp "$SRC/report_generator.py"     "$DEST/"
cp "$SRC/pipeline_profiler.py"    "$DEST/"
cp "$SRC/template.html"           "$DEST/"
cp -R "$SRC/example_data"       "$DEST/"
[ -f "$SRC/.dockerignore" ] && cp "$SRC/.dockerignore" "$DEST/"

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
