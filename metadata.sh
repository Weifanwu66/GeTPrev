#!/bin/bash
FULL_METADATA="database/assembly_summary_bacteria.txt"
GENOME_DIR="database/complete_genomes"
METADATA_DIR="database/metadata"
mkdir -p "$METADATA_DIR"
# Creata metadata directory structure mirroring complete_genomes
rsync -a --exclude='*' "$GENOME_DIR/" "$METADATA_DIR/"
extract_metadata() {
local search_dir="$1"
local output_file="$2"
if [[ ! -d "$search_dir" ]]; then
return
fi
genome_accessions=$(find "$search_dir" -type f -name "*_genomic.fna" -exec basename {} \; | awk -F'_' '{print $1"_"$2}')
grep -F -w -f <(echo "$genome_accessions") "$FULL_METADATA" > "$output_file"
echo "metadata saved: $output_file"
}
find "$GENOME_DIR" -mindepth 2 -type d | while read -r dir; do
relative_path="${dir#$GENOME_DIR/}"
metadata_path="$METADATA_DIR/$relative_path"
mkdir -p "$metadata_path"
metadata_file="${metadata_path}/metadata.csv"
extract_metadata "$dir" "$metadata_file"
done
