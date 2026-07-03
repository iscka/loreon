#!/bin/bash
# run_loreon.sh - Setup directories and launch LOREON Docker container

DOCS_DIR="$HOME/Documents"
LOREON_DIR="$DOCS_DIR/loreon"
DATA_DIR="$LOREON_DIR/data"
RESULTS_DIR="$LOREON_DIR/results"
DB_DIR="$LOREON_DIR/db"

# Create directories if they don't exist
mkdir -p "$DATA_DIR" "$RESULTS_DIR" "$DB_DIR"

echo "=== LOREON Pipeline ==="
echo "  Input:   $DATA_DIR"
echo "  Results: $RESULTS_DIR"
echo "  DB:      $DB_DIR"
echo ""

# Require at least -d (database file) argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 -d <db_filename.fasta> [options]"
    echo ""
    echo "Required:"
    echo "  -d <file>     FASTA database filename (must be inside $DB_DIR)"
    echo ""
    echo "Optional:"
    echo "  -f <format>   Taxonomy format: unite|silva|eukariome|cbs|none (default: unite)"
    echo "  -min <int>    Minimum read length (default: 200)"
    echo "  -max <int>    Maximum read length"
    echo "  -t <int>      Threads (default: all available)"
    echo ""
    echo "Input FASTQ files must be placed in: $DATA_DIR"
    echo "Results will be written to:          $RESULTS_DIR"
    exit 1
fi

MSYS_NO_PATHCONV=1 docker run --rm \
    -v "$DATA_DIR":/data/input \
    -v "$RESULTS_DIR":/data/output \
    -v "$DB_DIR":/data/db \
    loreon:latest \
    python3 /app/metaGenomics_new.py \
        -i /data/input \
        -o /data/output \
        "$@"
