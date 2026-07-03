#!/bin/bash
# =============================================================================
# LOREON — macOS DMG Installer Build Script
# Produces: dist/LOREON-1.0-macOS.dmg
#
# Creates a professional DMG with:
#   - LOREON.app bundle (PyInstaller --onedir)
#   - Symbolic link to /Applications for drag-and-drop install
#   - Custom volume name and compressed format
#
# Requirements (on the build machine):
#   - Python 3.8+  with pip
#   - PyQt5 and all project dependencies (pip install -r requirements.txt)
#   - Xcode Command Line Tools  (for sips / iconutil / hdiutil)
#
# Runtime dependencies (user must install separately):
#   - minimap2  (brew install minimap2  OR  conda install -c bioconda minimap2)
#   - samtools  (brew install samtools  OR  conda install -c bioconda samtools)
#
# Usage:
#   chmod +x build_dmg.sh
#   ./build_dmg.sh
#
# Gatekeeper note:
#   The resulting .app is NOT code-signed. First launch:
#   right-click → Open, or run:
#       xattr -d com.apple.quarantine /Applications/LOREON.app
# =============================================================================

set -e
cd "$(dirname "$0")"

# ── Configuration ────────────────────────────────────────────────────────────
APP_NAME="LOREON"
APP_VERSION="1.0"
BUNDLE_ID="com.loreon.pipeline"
DMG_NAME="${APP_NAME}-${APP_VERSION}-macOS.dmg"
DMG_VOLNAME="${APP_NAME} ${APP_VERSION}"

# Files to bundle inside the .app
DATA_FILES=(
    Dockerfile
    requirements.txt
    loreon_app_icon_2.png
    loreon.jpeg
    credits.md
    metaGenomics_new.py
    OtuUtils.py
    ResultsReader.py
    PipelineLogger.py
    report_generator.py
    pipeline_profiler.py
    pipeline_worker.py
    template.html
    fileorganizer.py
    venv_manager.py
)

# ── Helpers ──────────────────────────────────────────────────────────────────
info()    { printf "\033[1;34m[INFO]\033[0m  %s\n" "$*"; }
success() { printf "\033[1;32m[OK]\033[0m    %s\n" "$*"; }
warn()    { printf "\033[1;33m[WARN]\033[0m  %s\n" "$*"; }
error()   { printf "\033[1;31m[ERROR]\033[0m %s\n" "$*" >&2; exit 1; }

# ── Pre-flight checks ───────────────────────────────────────────────────────
info "Checking required files..."
[ -f "pipeline_gui.py" ]        || error "pipeline_gui.py not found."
[ -f "loreon_app_icon_2.png" ]  || error "loreon_app_icon_2.png not found."
[ -f "loreon.jpeg" ]            || error "loreon.jpeg not found."
[ -f "credits.md" ]             || error "credits.md not found."
[ -f "template.html" ]          || error "template.html not found."

# Check build tools
command -v python3   >/dev/null || error "python3 not found in PATH."
command -v sips      >/dev/null || error "sips not found. Install Xcode Command Line Tools."
command -v iconutil  >/dev/null || error "iconutil not found. Install Xcode Command Line Tools."
command -v hdiutil   >/dev/null || error "hdiutil not found (macOS only)."

# ── Clean previous build ────────────────────────────────────────────────────
info "Cleaning previous build artifacts..."
rm -rf build dist "${APP_NAME}.spec" 2>/dev/null || true

# ── Install build tools ─────────────────────────────────────────────────────
info "Installing/upgrading PyInstaller..."
python3 -m pip install pyinstaller --quiet --upgrade

# ── Build .icns icon ────────────────────────────────────────────────────────
info "Building .icns icon from loreon_app_icon_2.png..."
ICONSET="${APP_NAME}_build.iconset"
rm -rf "$ICONSET"
mkdir "$ICONSET"

for SIZE in 16 32 64 128 256 512; do
    sips -z $SIZE $SIZE loreon_app_icon_2.png --out "$ICONSET/icon_${SIZE}x${SIZE}.png" >/dev/null 2>&1
done
# Retina variants
sips -z 32   32   loreon_app_icon_2.png --out "$ICONSET/icon_16x16@2x.png"    >/dev/null 2>&1
sips -z 64   64   loreon_app_icon_2.png --out "$ICONSET/icon_32x32@2x.png"    >/dev/null 2>&1
sips -z 256  256  loreon_app_icon_2.png --out "$ICONSET/icon_128x128@2x.png"  >/dev/null 2>&1
sips -z 512  512  loreon_app_icon_2.png --out "$ICONSET/icon_256x256@2x.png"  >/dev/null 2>&1
sips -z 1024 1024 loreon_app_icon_2.png --out "$ICONSET/icon_512x512@2x.png"  >/dev/null 2>&1

ICNS_FILE="${APP_NAME}_build.icns"
iconutil -c icns "$ICONSET" -o "$ICNS_FILE"
rm -rf "$ICONSET"
success "Icon built: $ICNS_FILE"

# ── PyInstaller build ───────────────────────────────────────────────────────
info "Running PyInstaller (this may take a few minutes)..."
python3 -m PyInstaller \
    --onedir \
    --windowed \
    --name "$APP_NAME" \
    --icon "$ICNS_FILE" \
    --hidden-import=_socket \
    --hidden-import=select \
    --hidden-import=PyQt5.sip \
    --hidden-import=encodings.utf_8 \
    --hidden-import=encodings.ascii \
    --hidden-import=encodings.latin_1 \
    --collect-all=PyQt5 \
    --hidden-import=venv_manager \
    --clean \
    -y \
    pipeline_gui.py

success "PyInstaller build complete."

# ── Set bundle identifier in Info.plist ──────────────────────────────────────
PLIST="dist/${APP_NAME}.app/Contents/Info.plist"
if [ -f "$PLIST" ]; then
    info "Updating Info.plist with bundle identifier and version..."
    /usr/libexec/PlistBuddy -c "Set :CFBundleIdentifier ${BUNDLE_ID}" "$PLIST" 2>/dev/null || \
    /usr/libexec/PlistBuddy -c "Add :CFBundleIdentifier string ${BUNDLE_ID}" "$PLIST"

    /usr/libexec/PlistBuddy -c "Set :CFBundleShortVersionString ${APP_VERSION}" "$PLIST" 2>/dev/null || \
    /usr/libexec/PlistBuddy -c "Add :CFBundleShortVersionString string ${APP_VERSION}" "$PLIST"

    /usr/libexec/PlistBuddy -c "Set :CFBundleVersion ${APP_VERSION}" "$PLIST" 2>/dev/null || \
    /usr/libexec/PlistBuddy -c "Add :CFBundleVersion string ${APP_VERSION}" "$PLIST"
fi

# ── Copy data files into the .app bundle ────────────────────────────────────
BUNDLE_DIR="dist/${APP_NAME}.app/Contents/MacOS"
info "Copying data files into ${APP_NAME}.app bundle..."

for f in "${DATA_FILES[@]}"; do
    if [ -f "$f" ]; then
        cp "$f" "$BUNDLE_DIR/"
    else
        warn "File not found, skipping: $f"
    fi
done
[ -f .dockerignore ] && cp .dockerignore "$BUNDLE_DIR/"

success "Data files copied."

# ── Create DMG with Applications symlink ────────────────────────────────────
info "Creating DMG installer..."

# Stage directory for DMG contents
DMG_STAGE="dist/dmg_stage"
rm -rf "$DMG_STAGE"
mkdir -p "$DMG_STAGE"

# Copy .app into stage
cp -R "dist/${APP_NAME}.app" "$DMG_STAGE/"

# Create symbolic link to /Applications for drag-and-drop install
ln -s /Applications "$DMG_STAGE/Applications"

# Create a README for the DMG
cat > "$DMG_STAGE/README.txt" << 'READMETXT'
LOREON - Long-Read ONT Metagenomic Pipeline
============================================

PREREQUISITES:
  Python 3.8+ must be installed on your system.
  Install via: https://www.python.org/downloads/
  Or via Homebrew: brew install python@3

INSTALLATION:
  Drag LOREON.app into the Applications folder.

FIRST LAUNCH:
  The app is not code-signed. On first launch:
  - Right-click LOREON.app → Open → Open
  - Or run in Terminal:
    xattr -d com.apple.quarantine /Applications/LOREON.app

  On first launch, LOREON will automatically create a Python
  virtual environment and install all required dependencies.
  This may take 1-2 minutes (internet connection required).
  Subsequent launches will be instant.

RUNTIME DEPENDENCIES:
  minimap2 and samtools must be installed and available in PATH.
  Install via Homebrew:
    brew install minimap2 samtools
  Or via Conda:
    conda install -c bioconda minimap2 samtools

For more info: https://github.com/iscka/loreon
READMETXT

# Create the DMG
rm -f "dist/$DMG_NAME"
hdiutil create \
    -volname "$DMG_VOLNAME" \
    -srcfolder "$DMG_STAGE" \
    -ov \
    -format UDZO \
    -imagekey zlib-level=9 \
    "dist/$DMG_NAME"

# Cleanup stage
rm -rf "$DMG_STAGE"
rm -f "$ICNS_FILE"

# ── Summary ──────────────────────────────────────────────────────────────────
DMG_SIZE=$(du -sh "dist/$DMG_NAME" | cut -f1)

echo ""
echo "══════════════════════════════════════════════════════════"
echo "  BUILD COMPLETE"
echo "══════════════════════════════════════════════════════════"
echo "  App bundle  : dist/${APP_NAME}.app"
echo "  DMG image   : dist/${DMG_NAME}  (${DMG_SIZE})"
echo "══════════════════════════════════════════════════════════"
echo ""
echo "  Distribution:"
echo "    1. Share dist/${DMG_NAME}"
echo "    2. User opens DMG, drags LOREON to Applications"
echo "    3. First launch: right-click → Open"
echo ""
echo "  Runtime dependencies (user must install):"
echo "    brew install minimap2 samtools"
echo ""
