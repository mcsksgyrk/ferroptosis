#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
OUTPUT_FILE="${SCRIPT_DIR}/omnipath_interactions.txt"
URL="https://omnipathdb.org/interactions/?fields=sources&fields=references"

curl -# -o "$OUTPUT_FILE" "$URL"

if [ $? -eq 0 ]; then
    echo "Download completed"
else
    echo "Download failed"
    exit 1
fi
