#!/bin/bash
# =============================================================================
# LOREON — Ubuntu/Debian .deb Package Build Script
# Produces: dist/loreon_2.6_amd64.deb
#
# Creates a proper Debian package with:
#   - PyInstaller-bundled executable in /opt/loreon/
#   - Desktop entry for application menu integration
#   - Symbolic link /usr/local/bin/loreon for CLI access
#   - Dependency declarations for system libraries
#   - Pre/post install scripts for proper integration
#
# Requirements (on the build machine):
#   - Ubuntu 20.04+ / Debian 11+
#   - Python 3.8+  with pip
#   - dpkg-deb  (usually pre-installed)
#   - fakeroot  (sudo apt install fakeroot)
#   - PyQt5 system libraries:
#       sudo apt install python3-pyqt5 libxcb-xinerama0 libxcb-cursor0 libgl1
#   - All project Python dependencies:
#       pip install -r requirements.txt
#
# Runtime dependencies (automatically declared in .deb, installed separately):
#   - minimap2  (conda install -c bioconda minimap2)
#   - samtools  (sudo apt install samtools  OR  conda install -c bioconda samtools)
#
# Usage:
#   chmod +x build_deb.sh
#   ./build_deb.sh
#
# Install the .deb:
#   sudo dpkg -i dist/loreon_2.6_amd64.deb
#   sudo apt-get install -f    # resolve any missing dependencies
#
# Uninstall:
#   sudo dpkg -r loreon
# =============================================================================

set -e
cd "$(dirname "$0")"

# ── Configuration ────────────────────────────────────────────────────────────
APP_NAME="loreon"
APP_NAME_DISPLAY="LOREON"
APP_VERSION="2.6"
ARCH="amd64"
DEB_NAME="${APP_NAME}_${APP_VERSION}_${ARCH}.deb"
INSTALL_PREFIX="/opt/loreon"

MAINTAINER="Roberto Scarponi <roberto.scarponi@dottorandi.unipg.it>"
DESCRIPTION="Long-Read ONT Metagenomic Pipeline"
HOMEPAGE="https://github.com/iscka/loreon"

# Files to bundle inside the package
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
command -v dpkg-deb  >/dev/null || error "dpkg-deb not found. Install: sudo apt install dpkg"
command -v fakeroot  >/dev/null || error "fakeroot not found. Install: sudo apt install fakeroot"

# Detect architecture
MACHINE_ARCH=$(dpkg --print-architecture 2>/dev/null || echo "amd64")
if [ "$MACHINE_ARCH" != "$ARCH" ]; then
    warn "Build machine is $MACHINE_ARCH but package targets $ARCH."
    warn "Set ARCH variable if building for different architecture."
    ARCH="$MACHINE_ARCH"
    DEB_NAME="${APP_NAME}_${APP_VERSION}_${ARCH}.deb"
fi

# ── Clean previous build ────────────────────────────────────────────────────
info "Cleaning previous build artifacts..."
rm -rf build dist "${APP_NAME_DISPLAY}.spec" deb_build 2>/dev/null || true

# ── Install build tools ─────────────────────────────────────────────────────
info "Installing/upgrading PyInstaller..."
pipx install pyinstaller

# ── PyInstaller build ───────────────────────────────────────────────────────
info "Running PyInstaller (this may take a few minutes)..."
python3 -m PyInstaller \
    --onedir \
    --windowed \
    --name "$APP_NAME_DISPLAY" \
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

# ── Copy data files into PyInstaller output ─────────────────────────────────
PYINST_DIR="dist/${APP_NAME_DISPLAY}"
info "Copying data files into PyInstaller output..."

for f in "${DATA_FILES[@]}"; do
    if [ -f "$f" ]; then
        cp "$f" "$PYINST_DIR/"
    else
        warn "File not found, skipping: $f"
    fi
done
[ -f .dockerignore ] && cp .dockerignore "$PYINST_DIR/"

success "Data files copied."

# ── Build Debian package structure ──────────────────────────────────────────
info "Building Debian package structure..."

DEB_ROOT="deb_build"
rm -rf "$DEB_ROOT"

# Application files → /opt/loreon/
mkdir -p "$DEB_ROOT${INSTALL_PREFIX}"
cp -R "$PYINST_DIR/"* "$DEB_ROOT${INSTALL_PREFIX}/"

# Symbolic link → /usr/local/bin/loreon
mkdir -p "$DEB_ROOT/usr/local/bin"
# Will be created by postinst script (symlinks in deb can be tricky)

# Desktop entry → /usr/share/applications/
mkdir -p "$DEB_ROOT/usr/share/applications"
cat > "$DEB_ROOT/usr/share/applications/loreon.desktop" << EOF
[Desktop Entry]
Name=${APP_NAME_DISPLAY}
GenericName=Metagenomic Pipeline
Comment=${DESCRIPTION}
Exec=${INSTALL_PREFIX}/${APP_NAME_DISPLAY}
Path=${INSTALL_PREFIX}
Icon=${INSTALL_PREFIX}/loreon_app_icon_2.png
Terminal=false
Type=Application
Categories=Science;Biology;Education;
StartupWMClass=${APP_NAME_DISPLAY}
Keywords=metagenomics;bioinformatics;ONT;nanopore;fungi;
EOF

# Icon for desktop integration (multiple sizes via hicolor)
for SIZE in 48 128 256; do
    ICON_DIR="$DEB_ROOT/usr/share/icons/hicolor/${SIZE}x${SIZE}/apps"
    mkdir -p "$ICON_DIR"
    if command -v convert >/dev/null 2>&1; then
        # ImageMagick available: resize icon
        convert loreon_app_icon_2.png -resize ${SIZE}x${SIZE} "$ICON_DIR/loreon.png"
    else
        # Fallback: copy original
        cp loreon_app_icon_2.png "$ICON_DIR/loreon.png"
    fi
done

# Pixmaps fallback (legacy)
mkdir -p "$DEB_ROOT/usr/share/pixmaps"
cp loreon_app_icon_2.png "$DEB_ROOT/usr/share/pixmaps/loreon.png"

# Man page (optional, minimal)
mkdir -p "$DEB_ROOT/usr/share/man/man1"
cat > "$DEB_ROOT/usr/share/man/man1/loreon.1" << MANPAGE
.TH LOREON 1 "$(date +"%B %Y")" "${APP_VERSION}" "User Commands"
.SH NAME
loreon \- Long-Read ONT Metagenomic Pipeline
.SH SYNOPSIS
.B loreon
.br
.B loreon-cli
\-i INPUT \-o OUTPUT \-d DATABASE [OPTIONS]
.SH DESCRIPTION
LOREON is a bioinformatics pipeline for processing Oxford Nanopore
Technologies (ONT) long-read metagenomic sequencing data. It integrates
minimap2 alignment, samtools filtering, DuckDB-based OTU aggregation,
and an interactive HTML report generator into a PyQt5 graphical interface.
.SH GUI MODE
Launch the graphical interface:
.PP
.B loreon
.SH CLI MODE
Run the pipeline from command line:
.PP
.B loreon-cli
\-i /path/to/fastq \-o /path/to/output \-d /path/to/database.fasta
.SH DEPENDENCIES
.B minimap2
and
.B samtools
must be installed and available in PATH.
.SH AUTHOR
Roberto Scarponi (University of Perugia)
.SH SEE ALSO
.UR https://github.com/iscka/loreon
GitHub Repository
.UE
MANPAGE
gzip -f "$DEB_ROOT/usr/share/man/man1/loreon.1"

# Documentation
mkdir -p "$DEB_ROOT/usr/share/doc/loreon"
cp credits.md "$DEB_ROOT/usr/share/doc/loreon/"
cat > "$DEB_ROOT/usr/share/doc/loreon/copyright" << COPYR
Format: https://www.debian.org/doc/packaging-manuals/copyright-format/1.0/
Upstream-Name: LOREON
Upstream-Contact: Roberto Scarponi <roberto.scarponi@unipg.it>
Source: https://github.com/iscka/loreon

Files: *
Copyright: $(date +%Y) Roberto Scarponi
License: MIT
 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 .
 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.
 .
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
COPYR

# ── DEBIAN control files ────────────────────────────────────────────────────
mkdir -p "$DEB_ROOT/DEBIAN"

# Calculate installed size (in KB)
INSTALLED_SIZE=$(du -sk "$DEB_ROOT" | cut -f1)

# control
cat > "$DEB_ROOT/DEBIAN/control" << CTRL
Package: ${APP_NAME}
Version: ${APP_VERSION}
Section: science
Priority: optional
Architecture: ${ARCH}
Installed-Size: ${INSTALLED_SIZE}
Depends: python3 (>= 3.8), python3-venv, python3-pip, libxcb-xinerama0, libxcb-cursor0, libgl1, libglib2.0-0, libfontconfig1, libxkbcommon0
Recommends: samtools
Suggests: minimap2
Maintainer: ${MAINTAINER}
Homepage: ${HOMEPAGE}
Description: ${DESCRIPTION}
 LOREON (Long-Read ONT Metagenomic Pipeline) is a comprehensive
 bioinformatics pipeline for processing Oxford Nanopore Technologies
 (ONT) long-read metagenomic sequencing data.
 .
 Features:
  - PyQt5 graphical interface for easy operation
  - minimap2-based long-read alignment
  - DuckDB-powered OTU table aggregation
  - Interactive HTML reports with Plotly visualizations
  - Support for UNITE, SILVA, Eukariome, CBS databases
 .
 Note: minimap2 and samtools must be installed separately.
 Install via: conda install -c bioconda minimap2 samtools
CTRL

# CLI launcher wrapper (more robust than a bare symlink)
mkdir -p "$DEB_ROOT/usr/local/bin"
cat > "$DEB_ROOT/usr/local/bin/loreon" << 'LAUNCHER'
#!/bin/bash
# LOREON GUI launcher — ensures correct working directory and library paths
APP_DIR="/opt/loreon"
export LD_LIBRARY_PATH="${APP_DIR}/_internal${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd "$APP_DIR"
exec "$APP_DIR/LOREON" "$@"
LAUNCHER
chmod 755 "$DEB_ROOT/usr/local/bin/loreon"

# postinst — runs after installation
cat > "$DEB_ROOT/DEBIAN/postinst" << 'POSTINST'
#!/bin/bash
set -e

# CLI script symlink
ln -sf /opt/loreon/metaGenomics_new.py /usr/local/bin/loreon-cli

# Make CLI script executable
chmod +x /opt/loreon/metaGenomics_new.py 2>/dev/null || true

# Update desktop database
if command -v update-desktop-database >/dev/null 2>&1; then
    update-desktop-database /usr/share/applications/ 2>/dev/null || true
fi

# Update icon cache
if command -v gtk-update-icon-cache >/dev/null 2>&1; then
    gtk-update-icon-cache -f -t /usr/share/icons/hicolor/ 2>/dev/null || true
fi

# Update man database
if command -v mandb >/dev/null 2>&1; then
    mandb -q 2>/dev/null || true
fi

echo ""
echo "══════════════════════════════════════════════════════════"
echo "  LOREON installed successfully!"
echo "══════════════════════════════════════════════════════════"
echo ""
echo "  Launch GUI:  loreon"
echo "  CLI mode:    loreon-cli -i INPUT -o OUTPUT -d DATABASE"
echo ""
echo "  On first launch, LOREON will automatically create a"
echo "  Python virtual environment and install dependencies."
echo "  This requires an internet connection (1-2 minutes)."
echo ""
echo "  Runtime dependencies (install if not present):"
echo "    conda install -c bioconda minimap2 samtools"
echo "    # or: sudo apt install samtools"
echo ""

POSTINST
chmod 755 "$DEB_ROOT/DEBIAN/postinst"

# prerm — runs before removal
cat > "$DEB_ROOT/DEBIAN/prerm" << 'PRERM'
#!/bin/bash
set -e

# Remove launcher wrapper and CLI symlink
rm -f /usr/local/bin/loreon
rm -f /usr/local/bin/loreon-cli

PRERM
chmod 755 "$DEB_ROOT/DEBIAN/prerm"

# postrm — runs after removal
cat > "$DEB_ROOT/DEBIAN/postrm" << 'POSTRM'
#!/bin/bash
set -e

# Update desktop database
if command -v update-desktop-database >/dev/null 2>&1; then
    update-desktop-database /usr/share/applications/ 2>/dev/null || true
fi

# Update icon cache
if command -v gtk-update-icon-cache >/dev/null 2>&1; then
    gtk-update-icon-cache -f -t /usr/share/icons/hicolor/ 2>/dev/null || true
fi

# On purge, remove /opt/loreon completely (config/logs may remain)
if [ "$1" = "purge" ]; then
    rm -rf /opt/loreon
fi

POSTRM
chmod 755 "$DEB_ROOT/DEBIAN/postrm"

# conffiles (empty — no config files in /etc)
touch "$DEB_ROOT/DEBIAN/conffiles"

# ── Set permissions ─────────────────────────────────────────────────────────
info "Setting file permissions..."

# All files owned by root:root, readable by all
find "$DEB_ROOT" -type d -exec chmod 755 {} \;
find "$DEB_ROOT" -type f -exec chmod 644 {} \;

# Executables
chmod 755 "$DEB_ROOT${INSTALL_PREFIX}/${APP_NAME_DISPLAY}"
chmod 755 "$DEB_ROOT${INSTALL_PREFIX}/metaGenomics_new.py"
find "$DEB_ROOT${INSTALL_PREFIX}" -name "*.so" -exec chmod 755 {} \; 2>/dev/null || true
find "$DEB_ROOT${INSTALL_PREFIX}" -name "*.so.*" -exec chmod 755 {} \; 2>/dev/null || true
# PyInstaller helper binaries inside _internal/ (e.g. python3, QtWebEngineProcess)
find "$DEB_ROOT${INSTALL_PREFIX}/_internal" -type f -executable -exec chmod 755 {} \; 2>/dev/null || true
find "$DEB_ROOT${INSTALL_PREFIX}/_internal" -type f -name "python*" -exec chmod 755 {} \; 2>/dev/null || true
find "$DEB_ROOT${INSTALL_PREFIX}/_internal" -type f -name "Qt*Process" -exec chmod 755 {} \; 2>/dev/null || true

# DEBIAN scripts
chmod 755 "$DEB_ROOT/DEBIAN/postinst"
chmod 755 "$DEB_ROOT/DEBIAN/prerm"
chmod 755 "$DEB_ROOT/DEBIAN/postrm"

# ── Build .deb package ──────────────────────────────────────────────────────
info "Building .deb package (this may take a moment)..."
mkdir -p dist
fakeroot dpkg-deb --build "$DEB_ROOT" "dist/$DEB_NAME"

success ".deb package built."

# ── Verify package ──────────────────────────────────────────────────────────
info "Verifying package..."
dpkg-deb --info "dist/$DEB_NAME"
echo ""
dpkg-deb --contents "dist/$DEB_NAME" | head -30
echo "  ... (truncated)"

# ── Lint (optional) ─────────────────────────────────────────────────────────
if command -v lintian >/dev/null 2>&1; then
    info "Running lintian checks..."
    lintian "dist/$DEB_NAME" --no-tag-display-limit 2>/dev/null || true
fi

# ── Cleanup ─────────────────────────────────────────────────────────────────
rm -rf "$DEB_ROOT"

# ── Summary ─────────────────────────────────────────────────────────────────
DEB_SIZE=$(du -sh "dist/$DEB_NAME" | cut -f1)

echo ""
echo "══════════════════════════════════════════════════════════"
echo "  BUILD COMPLETE"
echo "══════════════════════════════════════════════════════════"
echo "  Package     : dist/${DEB_NAME}  (${DEB_SIZE})"
echo "  Architecture: ${ARCH}"
echo "══════════════════════════════════════════════════════════"
echo ""
echo "  Install:"
echo "    sudo dpkg -i dist/${DEB_NAME}"
echo "    sudo apt-get install -f    # resolve dependencies"
echo ""
echo "  Or with apt directly:"
echo "    sudo apt install ./dist/${DEB_NAME}"
echo ""
echo "  Uninstall:"
echo "    sudo dpkg -r loreon"
echo ""
echo "  After install:"
echo "    loreon              # launch GUI"
echo "    loreon-cli --help   # CLI mode"
echo ""
echo "  Runtime dependencies (install separately):"
echo "    conda install -c bioconda minimap2 samtools"
echo "    # or: sudo apt install samtools"
echo ""
