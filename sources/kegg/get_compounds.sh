#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
COMPOUND_FILE="${SCRIPT_DIR}/kegg_compounds.txt"
DRUG_FILE="${SCRIPT_DIR}/kegg_drugs.txt"

COMPOUND_URL="https://rest.kegg.jp/list/compound"
DRUG_URL="https://rest.kegg.jp/list/drug"

curl -# -o "$COMPOUND_FILE" "$COMPOUND_URL"
if [ $? -eq 0 ]; then
    echo "Compunds download completed"
else
    echo "Compounds download failed"
    exit 1
fi

curl -# -o "$DRUG_FILE" "$DRUG_URL"
if [ $? -eq 0 ]; then
    echo "Drugs download completed"
else
    echo "Drugs download failed"
    exit 1
fi
