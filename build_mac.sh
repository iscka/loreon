#!/bin/bash
# =============================================================================
# LOREON — macOS Build Script
# Produces: dist/LOREON.app  +  dist/LOREON_mac.dmg
#
# Requirements (on the build machine):
#   - Python 3.8+  with pip
#   - PyQt5 and all project dependencies (pip install -r requirements.txt)
#   - Xcode Command Line Tools  (for sips / iconutil / hdiutil)
#
# Usage:
#   chmod +x build_mac.sh
#   ./build_mac.sh
#
# Gatekeeper note:
#   The resulting .app is NOT code-signed. The first time you open it,
#   right-click → Open, or run:
#       xattr -d com.apple.quarantine dist/LOREON.app
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
info "Installing PyInstaller and Pillow..."
python3 -m pip install pyinstaller pillow --quiet

# ── Build ICNS icon ───────────────────────────────────────────────────────────
info "Building .icns icon from loreon_app_icon_2.png..."
ICONSET="LOREON_build.iconset"
rm -rf "$ICONSET"
mkdir "$ICONSET"

sips -z 16   16   loreon_app_icon_2.png --out "$ICONSET/icon_16x16.png"      >/dev/null
sips -z 32   32   loreon_app_icon_2.png --out "$ICONSET/icon_16x16@2x.png"   >/dev/null
sips -z 32   32   loreon_app_icon_2.png --out "$ICONSET/icon_32x32.png"      >/dev/null
sips -z 64   64   loreon_app_icon_2.png --out "$ICONSET/icon_32x32@2x.png"   >/dev/null
sips -z 128  128  loreon_app_icon_2.png --out "$ICONSET/icon_128x128.png"    >/dev/null
sips -z 256  256  loreon_app_icon_2.png --out "$ICONSET/icon_128x128@2x.png" >/dev/null
sips -z 256  256  loreon_app_icon_2.png --out "$ICONSET/icon_256x256.png"    >/dev/null
sips -z 512  512  loreon_app_icon_2.png --out "$ICONSET/icon_256x256@2x.png" >/dev/null
sips -z 512  512  loreon_app_icon_2.png --out "$ICONSET/icon_512x512.png"    >/dev/null
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
    --clean \
    -y \
    pipeline_gui.py

# ── Copy data files into the .app bundle ─────────────────────────────────────
BUNDLE="dist/LOREON.app/Contents/MacOS"
info "Copying data files into $BUNDLE ..."

cp Dockerfile              "$BUNDLE/"
cp requirements.txt        "$BUNDLE/"
cp loreon_app_icon_2.png   "$BUNDLE/"
cp loreon.jpeg             "$BUNDLE/"
cp credits.md              "$BUNDLE/"
cp metaGenomics_new.py     "$BUNDLE/"
cp OtuUtils.py             "$BUNDLE/"
cp ResultsReader.py        "$BUNDLE/"
cp PipelineLogger.py       "$BUNDLE/"
cp report_generator.py     "$BUNDLE/"
cp pipeline_profiler.py    "$BUNDLE/"
cp template.html           "$BUNDLE/"
[ -f .dockerignore ] && cp .dockerignore "$BUNDLE/"

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
