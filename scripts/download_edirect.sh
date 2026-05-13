#!/usr/bin/env bash
# Download GenBank records by taxon using NCBI EDirect tools.
#
# Usage:
#   bash scripts/download_edirect.sh -t 6270 -e you@example.com -k YOUR_KEY
#   bash scripts/download_edirect.sh -t 6270 -e you@example.com -r  # resume
#
# Requirements: edirect (esearch, efetch)
# Output: batch_001.gb, batch_002.gb, ... in output directory

set -euo pipefail

# Defaults
TAXON=""
OUTPUT="./gb_download"
BATCH_SIZE=500
EMAIL=""
API_KEY=""
RESUME=false

usage() {
    cat <<EOF
Usage: $0 -t TAXON_ID -e EMAIL [-o OUTPUT] [-b BATCH_SIZE] [-k API_KEY] [-r]

Options:
  -t TAXON_ID    NCBI Taxon ID (e.g. 6270)
  -e EMAIL       Email for NCBI
  -o OUTPUT      Output directory (default: ./gb_download)
  -b BATCH_SIZE  Records per batch (default: 500)
  -k API_KEY     NCBI API key
  -r             Resume from previous download
EOF
    exit 1
}

while getopts "t:e:o:b:k:rh" opt; do
    case $opt in
        t) TAXON="$OPTARG" ;;
        e) EMAIL="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        b) BATCH_SIZE="$OPTARG" ;;
        k) API_KEY="$OPTARG" ;;
        r) RESUME=true ;;
        h|*) usage ;;
    esac
done

if [ -z "$TAXON" ] || [ -z "$EMAIL" ]; then
    echo "Error: --taxon (-t) and --email (-e) are required"
    usage
fi

# Check edirect tools
for cmd in esearch efetch; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "Error: '$cmd' not found. Install edirect first:"
        echo "  sh -c \"\$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)\""
        exit 1
    fi
done

mkdir -p "$OUTPUT"

# Rate limit delay
if [ -n "$API_KEY" ]; then
    DELAY=0.1
    API_PARAM="-api_key $API_KEY"
else
    DELAY=0.34
    API_PARAM=""
fi

PROGRESS_FILE="$OUTPUT/download_progress.txt"

# Get total record count
COUNT=$(esearch -db nuccore -query "txid${TAXON}[Organism]" $API_PARAM -email "$EMAIL" 2>/dev/null \
    | xtract -pattern Count -element Count)

if [ -z "$COUNT" ] || [ "$COUNT" -eq 0 ]; then
    echo "No records found for txid${TAXON}"
    exit 1
fi

NUM_BATCHES=$(( (COUNT + BATCH_SIZE - 1) / BATCH_SIZE ))

# Resume: read last completed batch
START_BATCH=0
if $RESUME && [ -f "$PROGRESS_FILE" ]; then
    LAST=$(cat "$PROGRESS_FILE")
    START_BATCH=$((LAST + 1))
    echo "Resuming from batch $((START_BATCH)) ($LAST completed)"
fi

echo "Taxon: txid${TAXON}, ${COUNT} records, ${NUM_BATCHES} batches of ${BATCH_SIZE}"

for BATCH_IDX in $(seq $START_BATCH $((NUM_BATCHES - 1))); do
    RETSTART=$((BATCH_IDX * BATCH_SIZE))
    RETMAX=$BATCH_SIZE
    BATCH_NUM=$((BATCH_IDX + 1))
    BATCH_FILE="$OUTPUT/batch_$(printf '%04d' $BATCH_NUM).gb"

    END=$((RETSTART + RETMAX))
    if [ "$END" -gt "$COUNT" ]; then
        END=$COUNT
    fi

    printf "  Batch %d/%d (%d-%d/%d)... " "$BATCH_NUM" "$NUM_BATCHES" "$((RETSTART + 1))" "$END" "$COUNT"

    if esearch -db nuccore -query "txid${TAXON}[Organism]" $API_PARAM -email "$EMAIL" 2>/dev/null \
        | efetch -format gb -retstart "$RETSTART" -retmax "$RETMAX" \
        > "$BATCH_FILE" 2>/dev/null; then

        SIZE=$(stat -f%z "$BATCH_FILE" 2>/dev/null || stat -c%s "$BATCH_FILE" 2>/dev/null || echo "0")
        SIZE_KB=$((SIZE / 1024))
        echo "OK (${SIZE_KB} KB)"

        # Save progress
        echo "$BATCH_IDX" > "$PROGRESS_FILE"
    else
        echo "FAILED"
        echo "Run with -r to resume from batch $BATCH_NUM"
        exit 1
    fi

    sleep "$DELAY"
done

# Summary
NUM_FILES=$(ls "$OUTPUT"/batch_*.gb 2>/dev/null | wc -l)
echo ""
echo "Done! $NUM_FILES batch files in $OUTPUT"
