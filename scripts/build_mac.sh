#!/bin/bash
# =============================================================================
# LOREON — macOS Build Script
# Produces: dist/LOREON.app  +  dist/LOREON_mac.dmg
#
# Requirements (on the build machine):
#   - Python 3.8+  with pip
#   - PyQt5 and all project dependencies (pip install -r sources/requirements.txt)
#   - Xcode Command Line Tools  (for sips / iconutil / hdiutil)
#
# Usage:
#   chmod +x scripts/build_mac.sh
#   scripts/build_mac.sh
#
# Gatekeeper note:
#   The resulting .app is NOT code-signed. The first time you open it,
#   right-click → Open, or run:
#       xattr -d com.apple.quarantine dist/LOREON.app
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
info "Installing PyInstaller and Pillow..."
python3 -m pip install pyinstaller pillow --quiet

# ── Build ICNS icon ───────────────────────────────────────────────────────────
info "Building .icns icon from loreon_app_icon_2.png..."
ICONSET="LOREON_build.iconset"
rm -rf "$ICONSET"
mkdir "$ICONSET"

sips -z 16   16   "$SRC/loreon_app_icon_2.png" --out "$ICONSET/icon_16x16.png"      >/dev/null
sips -z 32   32   "$SRC/loreon_app_icon_2.png" --out "$ICONSET/icon_16x16@2x.png"   >/dev/null
sips -z 32   32   "$SRC/loreon_app_icon_2.png" --out "$ICONSET/icon_32x32.png"      >/dev/null
sips -z 64   64   "$SRC/loreon_app_icon_2.png" --out "$ICONSET/icon_32x32@2x.png"   >/dev/null
sips -z 128  128  "$SRC/loreon_app_icon_2.png" --out "$ICONSET/icon_128x128.png"    >/dev/null
sips -z 256  256  "$SRC/loreon_app_icon_2.png" --out "$ICONSET/icon_128x128@2x.png" >/dev/null
sips -z 256  256  "$SRC/loreon_app_icon_2.png" --out "$ICONSET/icon_256x256.png"    >/dev/null
sips -z 512  512  "$SRC/loreon_app_icon_2.png" --out "$ICONSET/icon_256x256@2x.png" >/dev/null
sips -z 512  512  "$SRC/loreon_app_icon_2.png" --out "$ICONSET/icon_512x512.png"    >/dev/null
iconutil -c icns "$ICONSET" -o loreon_build.icns
rm -rf "$ICONSET"

# ── PyInstaller build ─────────────────────────────────────────────────────────
info "Running PyInstaller..."
python3 -m PyInstaller \
    --onedir \
    --windowed \
    --name LOREON \
    --icon loreon_build.icns \
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

# ── Copy data files into the .app bundle ─────────────────────────────────────
BUNDLE="dist/LOREON.app/Contents/MacOS"
info "Copying data files into $BUNDLE ..."

cp "$SRC/Dockerfile"              "$BUNDLE/"
cp "$SRC/requirements.txt"        "$BUNDLE/"
cp "$SRC/loreon_app_icon_2.png"   "$BUNDLE/"
cp "$SRC/loreon.jpeg"             "$BUNDLE/"
cp "$SRC/credits.md"              "$BUNDLE/"
cp LICENSE                 "$BUNDLE/"
cp "$SRC/metaGenomics_new.py"     "$BUNDLE/"
cp "$SRC/OtuUtils.py"             "$BUNDLE/"
cp "$SRC/ResultsReader.py"        "$BUNDLE/"
cp "$SRC/PipelineLogger.py"       "$BUNDLE/"
cp "$SRC/report_generator.py"     "$BUNDLE/"
cp "$SRC/pipeline_profiler.py"    "$BUNDLE/"
cp "$SRC/template.html"           "$BUNDLE/"
cp -R "$SRC/example_data"       "$BUNDLE/"
[ -f "$SRC/.dockerignore" ] && cp "$SRC/.dockerignore" "$BUNDLE/"

# ── Create .dmg ───────────────────────────────────────────────────────────────
info "Creating LOREON_mac.dmg ..."
rm -f dist/LOREON_mac.dmg
hdiutil create \
    -volname  "LOREON" \
    -srcfolder "dist/LOREON.app" \
    -ov \
    -format UDZO \
    "dist/LOREON_mac.dmg"

# ── Cleanup ───────────────────────────────────────────────────────────────────
rm -f loreon_build.icns

echo ""
echo "══════════════════════════════════════════"
echo "  BUILD COMPLETE"
echo "  App bundle : dist/LOREON.app"
echo "  Disk image : dist/LOREON_mac.dmg"
echo "══════════════════════════════════════════"
echo ""
echo "NOTE: The app is not code-signed."
echo "  First launch: right-click → Open"
echo "  Or: xattr -d com.apple.quarantine dist/LOREON.app"
