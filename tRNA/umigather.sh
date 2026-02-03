#!/usr/bin/env bash

# Usage:
#   ./extract_all_umi_logs.sh <input_folder> <output_csv>
# Example:
#   ./extract_all_umi_logs.sh /path/to/logs summary.csv

INPUT_FOLDER="$1"
OUTPUT_CSV="$2"

echo "sample,input_reads,reads_out" > "$OUTPUT_CSV"

for LOGFILE in "$INPUT_FOLDER"/*.umi_log; do
    BASENAME=$(basename "$LOGFILE")
    SAMPLE=${BASENAME%.umi_log}

    # Extract numbers using robust AWK
    INPUT_READS=$(awk '/Reads: Input Reads:/ {print $NF}' "$LOGFILE" | head -n 1)
    READS_OUT=$(awk '/Number of reads out:/ {print $NF}' "$LOGFILE" | head -n 1)

    # Handle missing values
    if [[ -z "$INPUT_READS" ]]; then
        INPUT_READS=NA
    fi
    if [[ -z "$READS_OUT" ]]; then
        READS_OUT=NA
    fi

    echo "${SAMPLE},${INPUT_READS},${READS_OUT}" >> "$OUTPUT_CSV"
done
