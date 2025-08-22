#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
OUTPUT_FILE="${SCRIPT_DIR}/kegg_compounds.txt"

URL="https://rest.kegg.jp/list/compound"

curl -# -o "$OUTPUT_FILE" "$URL"

if [ $? -eq 0 ]; then
    echo "Download completed"
else
    echo "Download failed"
    exit 1
fi
