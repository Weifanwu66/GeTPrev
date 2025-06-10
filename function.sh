#!/bin/bash
WORKDIR="$(pwd)"
MONOPHASIC_TYPHIMURIUM_LIST="$WORKDIR/database/monophasic_Typhimurium_list.txt"
DUPLICATE_SEROTYPE_LIST="$WORKDIR/database/duplicate_sal_serotypes.txt"
METADATA_FILE="$WORKDIR/database/assembly_summary_bacteria.txt"
BLAST_DB_DIR="$WORKDIR/database/complete_blast_db"
ASSEMBLY_LEVEL="complete"

# check if dependencies are all installed
function check_dependencies() {
local missing=0
local deps=(blastn makeblastdb ncbi-genome-download esearch efetch xtract)
for cmd in "${deps[@]}"; do
if ! command -v "$cmd" >/dev/null 2>&1; then
echo "error: Required command '$cmd' not found in PATH" >&2
missing=1
fi
done
if (( missing )); then
echo "One or more dependencies are missing. Please activate the correct environment." >&2
exit 1
else
echo "All dependencies are installed."
fi
}

# get CPU count
function get_cpus() {
if [[ -n "${SLURM_CPUS_ON_NODE:-}" ]]; then
echo "$SLURM_CPUS_ON_NODE"
else
nproc
fi
}

# protection mechanism with download retry
function download_with_retry() {
local attempts="${RETRY_ATTEMPTS:-3}"
local delay="${RETRY_DELAY:-5}"
local count=1
while (( count <= attempts )); do
"$@" && return 0
echo "Attempt $count/$attempts failed for: $*" >&2
sleep "$delay"
((count++))
done
echo "Command failed after $attempts attempts: $*" >> "$FAILED_FLAG"
return 1
}

# download genus-level genomes
function download_single_genus() {
local genus="$1"
local output_dir="$2"
local genus_dir="${output_dir}/${genus}"
if [[ -n "$(find "$genus_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for $genus, skipping"
return
fi
mkdir -p "$genus_dir"
echo "Downloading genus: $genus"
download_with_retry ncbi-genome-download bacteria --genera "$genus" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$genus_dir" --verbose --flat-output
find "$genus_dir" -type f -name "*_genomic.fna.gz" -exec gzip -d {} \;
echo "Downloaded and organized genomes for $genus"
}

function download_genus() {
local input="$1"
local output_dir="$2"
local total_cpus=$(get_cpus)
local cores_half=$(( total_cpus / 2 ))
local max_parallel_genus=$(( cores_half < 1 ? 1 : (cores_half > 6 ? 6 : cores_half) ))
if [[ -f "$input" ]]; then
while read -r genus; do
[[ -z "$genus" ]] && continue
(
download_single_genus "$genus" "$output_dir"
) &
while (( $(jobs -r | wc -l) >= max_parallel_genus )); do sleep 1; done
done < "$input"
wait
else
download_single_genus "$input" "$output_dir"
fi
}

# get all species names under each target genus
function get_species_list() {
local target_genus="$1"
local output_file="$2"
local max_attempts=5
local attempt=1
local success=false
mkdir -p "$(dirname "$output_file")"
while [[ $attempt -le $max_attempts ]]; do
echo "Fetching species list for $target_genus (Attempt $attempt/$max_attempts)..."
esearch -db taxonomy -query "${target_genus}[Subtree]" | \
efetch -format xml | \
xtract -pattern Taxon -element ScientificName | \
awk 'NF == 2 && $1 ~ /^[A-Z][a-z]+$/ && $2 ~ /^[a-z]+$/ {print $0}' > "$output_file.tmp"
if [[ -s "$output_file.tmp" ]]; then
mv "$output_file.tmp" "$output_file"
echo "Saved all species names of $target_genus to $output_file"
success=true
break
else
echo "Attempt $attempt failed or returned empty. Retrying after $((attempt * 3))s..."
sleep $((attempt * 3))
((attempt++))
fi
done
if [[ "$success" != true ]]; then
echo "ERROR: Failed to retrieve species list for $target_genus after $max_attempts attempts." >&2
rm -f "$output_file.tmp"
touch "$output_file"
fi
}

# download species-level genomes
function download_species() {
local target_input="$1"
local output_dir="$2"
local total_cpus=$(get_cpus)
local cores_half=$(( total_cpus / 2 ))
local max_parallel_species=$(( cores_half < 1 ? 1 : (cores_half > 4 ? 4 : cores_half) ))
download_single_species() {
local species="$1"
local genus=$(echo "$species" | awk '{print $1}')
local clean_species="${species// /_}"
species_dir="${output_dir}/${clean_species}"
mkdir -p "$species_dir"
if [[ -n "$(find "$species_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for $species, skipping downloading"
return
fi
echo "Downloading genomes for $species"
download_with_retry ncbi-genome-download bacteria --genera "$species" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$species_dir" --flat-output --verbose
find "$species_dir" -type f -name "*_genomic.fna.gz" -exec gzip -d {} \;
if [[ -z "$(find "$species_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "No genomes found for $species. Remove empty directory."
rm -rf "$species_dir"
fi
}
if [[ -f "$target_input" ]]; then
while read -r species; do
[[ -z "$species" ]] && continue
( download_single_species "$species" ) &
while (( $(jobs -r | wc -l) >= max_parallel_species )); do sleep 1; done
done < "$target_input"
wait
else
download_single_species "$target_input"
fi
echo "Downloaded and organized species genomes"
}

# get all subspecies names under Salmonella enterica
function get_salmonella_subsp_list() {
local output_file="$1"
local output_dir
output_dir=$(dirname "$output_file")
mkdir -p "$output_dir"
esearch -db taxonomy -query "Salmonella enterica[Subtree]" | efetch -format xml | \
xtract -pattern Taxon -element ScientificName | grep -v -E "serovar|str\." | sort -u | grep -v "Salmonella enterica$" | sed 's/Salmonella enterica subsp. //' > "$output_file"
echo "Saved all Salmonella enterica subsp. names in $output_file"
}

function download_salmonella_subsp() {
local subspecies_list_file="$1"
echo "Downloading Salmonella enterica subspecies"
while read -r subspecies; do
[[ -z "$subspecies" ]] && continue
safe_subsp=$(echo "$subspecies" | tr -d "'\"" | tr -cs '[:alnum:]_:.-' '_' | sed 's/^_//;s/_$//')
local subspecies_dir="${GENOME_DIR}/Salmonella_enterica_subsp_${safe_subsp}"
mkdir -p "$subspecies_dir"
if [[ -n "$(find "$subspecies_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for Salmonella enterica subsp. $subspecies"
continue
fi
echo "Downloading genomes for Salmonella enterica subsp. $subspecies"
download_with_retry ncbi-genome-download bacteria --genera "Salmonella enterica subsp. $subspecies" \
 --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$subspecies_dir" --verbose --flat-output
find "$subspecies_dir" -type f -name "*_genomic.fna.gz" -exec gzip -d {} \;
if [[ -z "$(find "$subspecies_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "No genomes found for $subspecies. Removing empty directory."
rm -rf "$subspecies_dir"
fi
done < "$subspecies_list_file"
}

# get all Salmonella serotype names
function get_salmonella_serotype_list() {
local output_file="$1"
local output_dir
output_dir=$(dirname "$output_file")
mkdir -p "$output_dir"
esearch -db taxonomy -query "Salmonella enterica subsp. enterica[Subtree]" | efetch -format xml | xtract -pattern Taxon -element ScientificName | grep -v -E "str\.|var\." | sort -u | grep -v "Salmonella enterica subsp. enterica$" > "$output_file"
echo "Saved all serotype names in $output_file"
# remove all Salmonella enterica subsp. enterica serovar prefix
sed -i 's/Salmonella enterica subsp. enterica serovar //g' "$output_file"
# delete known serotype names names by name match that are already stored in monophasic_Typhimurium_list.txt and Typhimurium_list.txt
grep -vFf "$DUPLICATE_SEROTYPE_LIST" "$output_file" > tmp && mv tmp "$output_file"
} 

# download Salmonella serotypes
function download_salmonella_serotype() {
local serotype_list_file="$1"
local total_cpus=$(get_cpus)
local max_parallel_serotypes=$(( total_cpus / 4 ))
(( max_parallel_serotypes < 1 )) && max_parallel_serotypes=1
(( max_parallel_serotypes > 10 )) && max_parallel_serotypes=10

download_single_serotype() {
local serotype="$1"
local clean_serotype="${serotype// /_}"
local serotype_dir="${GENOME_DIR}/Salmonella_${clean_serotype}"
mkdir -p "$serotype_dir"
if [[ -n "$(find "$serotype_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for Salmonella $serotype, skipping."
return
fi
if [[ "$serotype" == "monophasic Typhimurium" ]]; then
echo "Downloading all monophasic Typhimurium genomes"
while read -r mono; do
[[ -z "$mono" ]] && continue
echo "Downloading: $mono"
download_with_retry ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $mono" \
 --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$serotype_dir" --flat-output
done < "$MONOPHASIC_TYPHIMURIUM_LIST"
find "$serotype_dir" -type f -name "*_genomic.fna.gz" -exec gunzip -f {} +
find "$serotype_dir" -type f -name "*.fna" -exec grep -Eil "serovar 43:a:1,7" {} \; -exec rm -f {} \;
return
fi
echo "Downloading genomes for: Salmonella enterica subsp. enterica serovar $serotype"
download_with_retry ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $serotype" \
 --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$serotype_dir" --flat-output
find "$serotype_dir" -type f -name "*_genomic.fna.gz" -exec gunzip -f {} +
if [[ "$serotype" == "Typhi" ]]; then
for fna in "$serotype_dir"/*.fna; do
if ! grep -Eiq "serovar Typhi([^a-zA-Z]|$)" "$fna" 2>/dev/null; then
echo "Deleting non-Typhi genome: $(basename "$fna")"
rm -f "$fna"
fi
done
fi
if [[ "$serotype" == "Typhimurium" ]]; then
for fna in "$serotype_dir"/*.fna; do
if ! grep -Eiq "serovar Typhimurium([^[:alpha:]]|$)" "$fna" 2>/dev/null || grep -Eiq "serovar Typhimurium.*var\." "$fna" 2>/dev/null; then
echo "Deleting non-Typhimurium genome: $(basename "$fna")"
rm -f "$fna"
fi
done
fi
if [[ -z "$(find "$serotype_dir" -type f -name "*.fna" 2>/dev/null)" ]]; then
echo "No usable genomes for $serotype. Removing directory."
rm -rf "$serotype_dir"
fi
}

while read -r serotype; do
[[ -z "$serotype" ]] && continue
( download_single_serotype "$serotype" ) &
while (( $(jobs -r | wc -l) >= max_parallel_serotypes )); do sleep 1; done
done < "$serotype_list_file"
wait
echo "Serotype downloads complete."
}

# classify Salmonella genus into sub-levels based on metadata
function classify_salmonella_by_metadata() {
local subspecies_list_file="$1"
local serotype_list_file="$2"
local genome_dir="$(pwd)/database/complete_genomes/Salmonella"
echo "Classifying Salmonella genomes using metadata"
mkdir -p "$genome_dir/Salmonella_enterica"
mkdir -p "$genome_dir/Salmonella_bongori"
mkdir -p "$genome_dir/unclassified"
mapfile -t subspecies_list < "$subspecies_list_file"
mapfile -t serotype_list < "$serotype_list_file"
mapfile -t mono_aliases < "$MONOPHASIC_TYPHIMURIUM_LIST"
mapfile -t deduped_seros < "$DUPLICATE_SEROTYPE_LIST"
for subsp in "${subspecies_list[@]}"; do
[[ -z "$subsp" ]] && continue
subsp_dir="$genome_dir/Salmonella_enterica/$subsp"
mkdir -p "$subsp_dir"
if [[ "$subsp" == "enterica" ]]; then
mkdir -p "$subsp_dir/monophasic_Typhimurium"
mkdir -p "$subsp_dir/Typhi"
mkdir -p "$subsp_dir/unclassified"
for sero in "${serotype_list[@]}"; do
[[ -z "$sero" ]] && continue
if printf '%s\n' "${deduped_seros[@]}" | grep -Fxq -- "$sero" 2>/dev/null; then
continue
fi
mkdir -p "$subsp_dir/$sero"
done
fi
done
declare -A accession_to_organism
while IFS=$'\t' read -r col1 _ _ _ _ _ _ col8 _; do
[[ "$col1" == \#* ]] && continue
accession_to_organism["$col1"]="$col8"
done < "$METADATA_FILE"
find "$genome_dir" -maxdepth 1 -type f -name "GCA_*_genomic.fna" | while read -r fasta_file; do
base_filename=$(basename "$fasta_file")
accession=$(echo "$base_filename" | sed -E 's/^(GCA_[0-9]+\.[0-9]+).*/\1/')
organism_name="${accession_to_organism[$accession]}"
target_dir="$genome_dir/unclassified"
if [[ "$organism_name" == Salmonella\ enterica* ]]; then
base_path="$genome_dir/Salmonella_enterica"
matched_subsp=false
for subsp in "${subspecies_list[@]}"; do
if [[ "$organism_name" == *subsp.*"$subsp"* ]]; then
matched_subsp=true
base_path="$base_path/$subsp"
matched_sero=false
for alias in "${mono_aliases[@]}"; do
if [[ "$organism_name" == *"$alias"* ]]; then
target_dir="$base_path/monophasic_Typhimurium"
matched_sero=true
break
fi
done
if [[ "$matched_sero" == false && "$organism_name" =~ serovar[[:space:]]*Typhi([^a-zA-Z]|$) ]]; then
target_dir="$base_path/Typhi"
matched_sero=true
fi
if [[ "$matched_sero" == false ]]; then
for sero in "${serotype_list[@]}"; do
if [[ "$organism_name" =~ serovar[[:space:]]*$sero([^a-zA-Z]|$) ]]; then
if printf '%s\n' "${deduped_seros[@]}" | grep -Fxq -- "$sero" 2>/dev/null; then
continue
fi
safe_sero_name=$(echo "$sero" | tr -d "'\"" | tr -cs '[:alnum:]_:.-' '_' | sed 's/^_//;s/_$//')
target_dir="$base_path/$safe_sero_name"
matched_sero=true
break
fi
done
fi
if [[ "$matched_sero" == false && "$subsp" == "enterica" ]]; then
target_dir="$base_path/unclassified"
fi
break
fi
done
if [[ "$matched_subsp" == false ]]; then
target_dir="$genome_dir/Salmonella_enterica/unclassified"
fi
elif [[ "$organism_name" == Salmonella\ bongori* ]]; then
target_dir="$genome_dir/Salmonella_bongori"
fi
mkdir -p "$target_dir"
mv "$fasta_file" "$target_dir/"
done
for subsp in "${subspecies_list[@]}"; do
subsp_dir="$genome_dir/Salmonella_enterica/$subsp"
if [[ "$subsp" == "enterica" ]]; then
for sero in "${serotype_list[@]}"; do
[[ -z "$sero" ]] && continue
if printf '%s\n' "${deduped_seros[@]}" | grep -Fxq -- "$sero" 2>/dev/null; then
continue
fi
sero_name=$(echo "$sero")
safe_sero_name=$(echo "$sero" | tr -d "'\"" | tr -cs '[:alnum:]_:.-' '_' | sed 's/^_//;s/_$//')
sero_dir="$subsp_dir/$safe_sero_name"
sero_old_dir="$subsp_dir/$sero_name"
[[ -d "$sero_dir" && -z "$(ls -A "$sero_dir")" ]] && rmdir "$sero_dir"
[[ -d "$sero_old_dir" && -z "$(ls -A "$sero_old_dir")" ]] && rmdir "$sero_old_dir"
done
[[ -d "$subsp_dir/unclassified" && -z "$(ls -A "$subsp_dir/unclassified")" ]] && rmdir "$subsp_dir/unclassified"
fi
[[ -d "$subsp_dir" && -z "$(ls -A "$subsp_dir")" ]] && rmdir "$subsp_dir"
done
echo "Finished classifying Salmonella genomes."
}

# classify other genus into species level based on metadata
function classify_genus_by_metadata() {
local genus="$1"
local genome_dir="$(pwd)/database/complete_genomes/$genus"
species_list_file="$genome_dir/species_list.txt"
mkdir -p "$genome_dir/unclassified"
declare -A accession_to_organism
while IFS=$'\t' read -r col1 _ _ _ _ _ _ col8 _; do
[[ "$col1" == \#* ]] && continue
accession_to_organism["$col1"]="$col8"
done < "$METADATA_FILE"
mapfile -t species_list < "$species_list_file"
for species in "${species_list[@]}"; do
[[ -z "$species" ]] && continue
clean_species=$(echo "$species" | sed 's/ /_/g')
mkdir -p "$genome_dir/$clean_species"
done
find "$genome_dir" -maxdepth 1 -type f -name "GCA_*_genomic.fna" | while read -r fasta_file; do
base_filename=$(basename "$fasta_file")
accession=$(echo "$base_filename" | sed -E 's/^(GCA_[0-9]+\.[0-9]+).*/\1/')
organism_name="${accession_to_organism[$accession]}"
matched=false
for species in "${species_list[@]}"; do
[[ -z "$species" ]] && continue
[[ "$organism_name" == "$species" ]] || [[ "$organism_name" == "$species "* ]] && {
clean_species=$(echo "$species" | sed 's/ /_/g')
target_dir="$genome_dir/$clean_species"
mv "$fasta_file" "$target_dir/"
matched=true
break
}
done
if [[ "$matched" == false ]]; then
mv "$fasta_file" "$genome_dir/unclassified/"
fi
done
[[ -d "$genome_dir/unclassified" && -z "$(ls -A "$genome_dir/unclassified")" ]] && rmdir "$genome_dir/unclassified"
for species in "${species_list[@]}"; do
clean_species=$(echo "$species" | sed 's/ /_/g')
[[ -d "$genome_dir/$clean_species" && -z "$(ls -A "$genome_dir/$clean_species")" ]] && rmdir "$genome_dir/$clean_species"
done
echo "Finished classifying $genus genomes."
}

# generate genome list for each directory
function generate_directory_csv() {
local base_dir="$1"
echo "[$(date)] Generating genomes_list.csv under $base_dir"
find "$base_dir" -type d | while read -r dir; do
genome_ids=()
while IFS= read -r file; do
accession=$(basename "$file" | sed -E 's/^(GCA_[0-9]+\.[0-9]+).*/\1/')
genome_ids+=("$accession")
done < <(find "$dir" -maxdepth 1 -type f \( -name "*_genomic.fna" -o -name "*_genomic.fna.gz" \))
if [[ ${#genome_ids[@]} -eq 0 ]]; then
continue
fi
output_csv="$dir/genomes_list.csv"
{
grep '^#' "$METADATA_FILE" | head -n 1 | sed 's/^#//'
for acc in "${genome_ids[@]}"; do
grep "^${acc}[[:space:]]" "$METADATA_FILE"
done
} > "$output_csv"
echo "[$(date)] Saved: $output_csv"
done
}

# check for duplicates
function sanity_check_unique_genomes() {
local genome_dir="$1"
local temp_all_files="$WORKDIR/all_fasta_files.txt"
local temp_all_accessions="$WORKDIR/all_accessions.txt"
> "$temp_all_files"
> "$temp_all_accessions"
find "$genome_dir" -type f -name "*_genomic.fna" >> "$temp_all_files"
while IFS= read -r f; do
accession=$(basename "$f" | sed -E 's/^(GCA_[0-9]+\.[0-9]+).*/\1/')
echo "$accession" >> "$temp_all_accessions"
done < "$temp_all_files"
total=$(wc -l < "$temp_all_accessions")
unique=$(sort "$temp_all_accessions" | uniq | wc -l)
echo "Sanity Check: Total genomes = $total, Unique = $unique"
if [[ "$total" -ne "$unique" ]]; then
echo "Warning: Duplicate genomes detected!" >> "$FAILED_FLAG"
fi
rm -f "$temp_all_files" "$temp_all_accessions"
}

# create blastdb
function create_blastdb() {
local INPUT_FASTA="$1"
local DB_PATH="$2"
makeblastdb -in "$INPUT_FASTA" -dbtype nucl -out "${DB_PATH}"
echo "BLAST database created for ${DB_PATH}"
}

function build_blastdb_for_EB_default() {
input_dir="$1"
output_dir="$2"
mkdir -p "$output_dir"
echo "Building BLAST database from genomes in $input_dir"
find "$input_dir" -mindepth 1 -maxdepth 1 -type d | while read -r TAXON_DIR; do
taxon_name=$(basename "$TAXON_DIR")
[[ "$taxon_name" == "unclassified" ]] && continue
subdir_count=$(find "$TAXON_DIR" -mindepth 1 -maxdepth 1 -type d ! -name "unclassified" | wc -l)
if [[ "$subdir_count" -eq 0 ]]; then
echo "Detected species-level directory: $taxon_name"
species_fasta="$TAXON_DIR/${taxon_name}_all_genomes.fna"
find "$TAXON_DIR" -type f -name "*_genomic.fna" | while read -r file; do accession=$(basename "$file" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+'); awk -v acc="$accession" '/^>/{print ">" acc "|" substr($0,2); next} 1' "$file"; done > "$species_fasta"
create_blastdb "$species_fasta" "$output_dir/$taxon_name"
else
echo "Detected genus-level directory: $taxon_name"
num_children=$(find "$TAXON_DIR" -mindepth 1 -maxdepth 1 -type d ! -name "unclassified" | wc -l)
if [[ "$num_children" -gt 1 ]]; then
genus_fasta="$TAXON_DIR/${taxon_name}_all_genomes.fna"
find "$TAXON_DIR" -type f -name "*_genomic.fna" | while read -r file; do accession=$(basename "$file" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+'); awk -v acc="$accession" '/^>/{print ">" acc "|" substr($0,2); next} 1' "$file"; done > "$genus_fasta"
create_blastdb "$genus_fasta" "$output_dir/$taxon_name"
else
echo "Skipping genus-level BLAST DB for $taxon_name (only one species/subdir)"
fi
find "$TAXON_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r species_dir; do
species_name=$(basename "$species_dir")
[[ "$species_name" == "unclassified" ]] && continue
species_fasta="$species_dir/${species_name}_all_genomes.fna"
find "$species_dir" -type f -name "*_genomic.fna" | while read -r file; do accession=$(basename "$file" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+'); awk -v acc="$accession" '/^>/{print ">" acc "|" substr($0,2); next} 1' "$file"; done > "$species_fasta"
create_blastdb "$species_fasta" "$output_dir/$species_name"
if [[ "$species_name" == "Salmonella_enterica" ]]; then
find "$species_dir" -mindepth 1 -maxdepth 1 -type d | while read -r subspecies_dir; do
subspecies_name=$(basename "$subspecies_dir")
[[ "$subspecies_name" == "unclassified" ]] && continue
echo "Detected Salmonella subspecies: $subspecies_name"
subspecies_fasta="$subspecies_dir/Salmonella_enterica_subsp_${subspecies_name}_all_genomes.fna"
find "$subspecies_dir" -type f -name "*_genomic.fna" | while read -r file; do accession=$(basename "$file" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+'); awk -v acc="$accession" '/^>/{print ">" acc "|" substr($0,2); next} 1' "$file"; done > "$subspecies_fasta"
create_blastdb "$subspecies_fasta" "$output_dir/Salmonella_enterica_subsp_${subspecies_name}"
if [[ "$subspecies_name" == "enterica" ]]; then
find "$subspecies_dir" -mindepth 1 -maxdepth 1 -type d | while read -r serotype_dir; do
serotype_name=$(basename "$serotype_dir")
[[ "$serotype_name" == "unclassified" ]] && continue
echo "Detected Salmonella enterica serotype: $serotype_name"
serotype_fasta="$serotype_dir/${serotype_name}_genomes.fna"
find "$serotype_dir" -type f -name "*_genomic.fna" | while read -r file; do accession=$(basename "$file" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+'); \
awk -v acc="$accession" '/^>/{print ">" acc "|" substr($0,2); next} 1' "$file"; done > "$serotype_fasta"
create_blastdb "$serotype_fasta" "$output_dir/Salmonella_${serotype_name}"
done
fi
done
fi
done
fi
done
find "$input_dir" -type f \( -name "*_genomic.fna" -o -name "*_all_genomes.fna" \) | while read -r f; do gzip -f "$f" & done
wait
echo "Finished building BLAST databases for default EB genomes"
}

function build_custom_blastdb() {
input_dir="$1"
output_dir="$2"
mkdir -p "$output_dir"
echo "Building BLAST databases for custom panel from: $input_dir"
find "$input_dir" -mindepth 1 -maxdepth 1 -type d | while read -r taxon_dir; do
taxon_name=$(basename "$taxon_dir")
[[ "$taxon_name" == "unclassified" ]] && continue
echo "Processing: $taxon_name"
combined_fasta="$taxon_dir/${taxon_name}_all_genomes.fna"
find "$taxon_dir" -type f -name "*_genomic.fna" | while read -r file; do
accession=$(basename "$file" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+')
awk -v acc="$accession" '/^>/{print ">" acc "|" substr($0,2); next} 1' "$file"
done > "$combined_fasta"
if [[ "$taxon_name" =~ ^Salmonella_enterica_subsp_([a-zA-Z0-9_]+)$ ]]; then
db_name="Salmonella_enterica_subsp_${BASH_REMATCH[1]}"
elif [[ "$taxon_name" =~ ^Salmonella_([A-Z][a-zA-Z0-9_]+)$ ]]; then
db_name="Salmonella_${BASH_REMATCH[1]}"
else
db_name="$taxon_name"
fi
create_blastdb "$combined_fasta" "$output_dir/$db_name"
find "$taxon_dir" -mindepth 1 -maxdepth 1 -type d | while read -r species_dir; do
species_name=$(basename "$species_dir")
species_fasta="$species_dir/${species_name}_all_genomes.fna"
find "$species_dir" -type f -name "*_genomic.fna" | while read -r file; do
accession=$(basename "$file" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+')
awk -v acc="$accession" '/^>/{print ">" acc "|" substr($0,2); next} 1' "$file"
done > "$species_fasta"
create_blastdb "$species_fasta" "$output_dir/$species_name"
done
done
find "$input_dir" -type f \( -name "*_genomic.fna" -o -name "*_all_genomes.fna" \) | while read -r f; do gzip -f "$f" & done
wait
echo "Finished building custom BLAST databases"
}

# clean up Salmonella serotype names
function extract_taxon_info() {
local input="$1"
local taxon_name=""
read -ra words <<< "$input"
if [[ "${words[1]}" == "monophasic" ]]; then
taxon_name="Salmonella enterica subsp. enterica serovar monophasic Typhimurium"
elif [[ "${words[1]}" =~ ^[A-Z] || "${words[1]}" =~ ":" ]]; then
taxon_name="Salmonella enterica subsp. enterica serovar ${words[1]}"
else
taxon_name="$input"
fi
echo "$taxon_name"
}

# delete genomes don't belong to that taxonomic group falsely downloaded by ncbi-genome-download
function get_strict_accessions() {
local target_name="$1" 
local accession_column=1
local name_column=8
local serovar=$(echo "$target_name" | awk '{for (i=6; i<=NF; i++) printf $i" "; print ""}' | sed 's/ *$//')
awk -F '\t' -v acc_col="$accession_column" -v name_col="$name_column" -v sv="$serovar" '
tolower($name_col) ~ ("^salmonella enterica subsp\\. enterica serovar " tolower(sv) "( |$)") &&
tolower($name_col) !~ /var\./ {print $acc_col
  }' "$METADATA_FILE"
}

# count the total downloaded (screened) genomes of a taxonomic group
function get_total_genomes_count() {
local input="$1"
local assembly_level="$2"
local total_genomes=0
local query="$(extract_taxon_info "$input")"
local clean_query="${query// /_}"
clean_query="${clean_query//./}"
local db_name=""
if [[ "$query" == "Salmonella enterica subsp. enterica serovar monophasic Typhimurium" ]]; then
db_name="Salmonella_monophasic_Typhimurium"
elif [[ "$clean_query" =~ ^Salmonella_enterica_subsp_enterica_serovar_([A-Z][a-zA-Z0-9_]+)$ ]]; then
db_name="Salmonella_${BASH_REMATCH[1]}"
elif [[ "$clean_query" =~ ^Salmonella_enterica_subsp_([a-zA-Z0-9_]+)$ ]]; then
db_name="Salmonella_enterica_subsp_${BASH_REMATCH[1]}"
else
db_name="$clean_query"
fi
local db_path="$BLAST_DB_DIR/$db_name"
if [[ "$assembly_level" == "complete" ]]; then
if [[ -f "${db_path}.nhr" || -f "${db_path}.00.nhr" ]]; then
total_genomes=$(blastdbcmd -db "$db_path" -entry all -outfmt "%t" | cut -d'|' -f1 | sort -u | wc -l)
else
echo "No genomes found for $input" >&2
total_genomes=0
fi
else
if [[ "$query" == "Salmonella enterica subsp. enterica serovar Typhi" ]]; then
all_accessions=$(ncbi-genome-download bacteria --genera "$query" --assembly-level "$assembly_level" --section genbank --dry-run | tail -n +2 | awk -F '/' '{print $NF}' | awk '{print $1}')
strict_accessions=$(get_strict_accessions "$query")
total_genomes=$(comm -12 <(echo "$all_accessions" | sort) <(echo "$strict_accessions" | sort) | wc -l)
elif [[ "$query" == "Salmonella enterica subsp. enterica serovar Typhimurium" ]]; then
all_accessions=$(ncbi-genome-download bacteria --genera "$query" --assembly-level "$assembly_level" --section genbank --dry-run | tail -n +2 | awk -F '/' '{print $NF}' | awk '{print $1}')
strict_accessions=$(get_strict_accessions "$query") 
total_genomes=$(comm -12 <(echo "$all_accessions" | sort) <(echo "$strict_accessions" | sort) | wc -l)
elif [[ "$query" == "Salmonella enterica subsp. enterica serovar monophasic Typhimurium" ]]; then
total_genomes=0
while read -r name; do
count=$(ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $name" \
 --assembly-level "$assembly_level" --section genbank --dry-run | tail -n +2 | grep -vi 'serovar 43:a:1,7' | wc -l)
((total_genomes += count))
done < "$MONOPHASIC_TYPHIMURIUM_LIST"
else
total_genomes=$(ncbi-genome-download bacteria --genera "$query" --assembly-level "$assembly_level" --section genbank --dry-run | tail -n +2 | wc -l)
fi
fi
echo "$total_genomes"
}

# decide the sample size per iteration and the number of iterations using Cochran's formula for draft genomes
function calculate_sample_size_and_iterations(){
local total_genomes="$1"
local max_iterations=20
if [[ "$total_genomes" -eq 0 ]]; then
echo "0 0"
return
fi
if [[ "$total_genomes" -le 100 ]]; then
echo "$total_genomes 1"
return
fi
# Cochran's formula constants
local Z=1.96
local p=0.5
local e=0.05
local n0=$(echo "scale=6; ($Z^2 * $p * (1 - $p)) / ($e^2)" | bc -l)
# finite population correction (FPC)
local sample_size=$(echo "scale=6; ($n0 * $total_genomes) / ($total_genomes + $n0 -1)" | bc -l)
sample_size=$(echo "$sample_size" | awk '{printf "%.0f", $1}')
# compute number of iterations using square-root scaling
if [[ "$sample_size" -gt 0 ]]; then
local iterations=$(echo "scale=6; sqrt($total_genomes / (2 * $sample_size))" | bc -l)
iterations=$(echo "$iterations" | awk '{printf "%.0f", $1}')
else
local iterations=1
fi
if [[ "$iterations" -lt 1 ]]; then
iterations=1
elif [[ "$iterations" -gt "$max_iterations" ]]; then
iterations="$max_iterations"
fi
echo "$sample_size $iterations"
}

# download random draft genomes
function download_random_draft_genomes() {
local input="$1"
local sample_size="$2"
local output_dir="$3"
local iteration="$4"
local query="$(extract_taxon_info "$input")"
local taxon_dir_name="${query// /_}"
taxon_dir_name="${taxon_dir_name//./}"
local iteration_dir="${output_dir}/${taxon_dir_name}/genomes_${iteration}"
mkdir -p "$iteration_dir"
local accessions=""
if [[ "$query" == "Salmonella enterica subsp. enterica serovar monophasic Typhimurium" ]]; then
while read -r actual_name; do
accessions+=$(ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $actual_name" \
 --assembly-level contig --section genbank --dry-run | tail -n +2 | \
 grep -vi 'serovar 43:1:1,7' | awk -F '/' '{print $NF}' | awk '{print $1}')$'\n'
done < "$MONOPHASIC_TYPHIMURIUM_LIST"
elif [[ "$query" == "Salmonella enterica subsp. enterica serovar Typhimurium" ]]; then
all_accessions=$(ncbi-genome-download bacteria --genera "$query" --assembly-level contig --section genbank --dry-run | tail -n +2 | awk -F '/' '{print $NF}' | awk '{print $1}')
strict_accessions=$(get_strict_accessions "$query") 
accessions=$(comm -12 <(echo "$all_accessions" | sort) <(echo "$strict_accessions" | sort))
elif [[ "$query" == "Salmonella enterica subsp. enterica serovar Typhi" ]]; then
all_accessions=$(ncbi-genome-download bacteria --genera "$query" --assembly-level contig --section genbank --dry-run | tail -n +2 | awk -F '/' '{print $NF}' | awk '{print $1}')
strict_accessions=$(get_strict_accessions "$query")
accessions=$(comm -12 <(echo "$all_accessions" | sort) <(echo "$strict_accessions" | sort))
else
accessions=$(ncbi-genome-download bacteria --genera "$query" --assembly-level contig --section genbank --dry-run | tail -n +2 | awk -F '/' '{print $NF}' | awk '{print $1}')
fi
local valid_accessions=$(echo "$accessions" | grep -E '^GCA_[0-9]+\.[0-9]+$' | shuf -n "$sample_size")
echo "$valid_accessions" > "${iteration_dir}/selected_accessions.txt"
download_with_retry ncbi-genome-download bacteria --assembly-accessions "${iteration_dir}/selected_accessions.txt" --formats fasta --assembly-level contig --section genbank --output-folder "$iteration_dir" --flat-output
find "$iteration_dir" -type f -name "*_genomic.fna.gz" -exec gunzip -f {} +
}

function perform_blast(){
local query_gene="$1"
local perc_identity="$2"
local output_dir="$3"
local iteration="$4" # only used for draft genomes
local taxon="$5"
if [[ -n "$iteration" ]]; then
echo "Processing draft genomes for iteration $iteration: $taxon"
local standard_taxon="$(extract_taxon_info "$taxon")"
local safe_taxon="${standard_taxon// /_}"
safe_taxon="${safe_taxon//./}"
local genome_dir="${DRAFT_GENOMES_DIR}/${safe_taxon}/genomes_${iteration}"
local blast_db="${DRAFT_BLAST_DB_DIR}/${safe_taxon}/iteration_${iteration}"
mkdir -p "$(dirname "$blast_db")"
if [ ! -d "$genome_dir" ] || [ -z "$(find "$genome_dir" -type f -name "*_genomic.fna" 2>/dev/null)" ]; then
echo "Error: No genomic.fna files found in $genome_dir. Check download process." >&2
exit 1
fi
local concatenated_genome="${genome_dir}/combined.fna"
echo "Concatenating genome FASTAs for $taxon..."
find "$genome_dir" -name "*_genomic.fna" -exec cat {} + > "$concatenated_genome"
echo "Building BLAST DB for $taxon..."
makeblastdb -in "$concatenated_genome" -dbtype nucl -out "$blast_db"
find "$genome_dir" -type f -name "*_genomic.fna" -exec gzip -f {} +
gzip -f "$concatenated_genome"
local blast_output="${output_dir}/${safe_taxon}/iteration_${iteration}_draft_blast_results.txt"
mkdir -p "$(dirname "$blast_output")"
echo "Running BLAST for $taxon..."
blastn -query "$query_gene" -db "$blast_db" -out "$blast_output" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
-perc_identity "$perc_identity" -max_target_seqs 10000
echo "BLAST results saved to: $blast_output"
else
# skip empty lines
[[ -z "$taxon" ]] && return
local blast_db_name="${taxon// /_}"
blast_db_name="${blast_db_name//./}"
local blast_db="${BLAST_DB_DIR}/${blast_db_name}"
local clean_taxon="$(extract_taxon_info "$taxon")"
local blast_output_name="${clean_taxon// /_}"
blast_output_name="${blast_output_name//./}"
local blast_output="${output_dir}/${blast_output_name}_complete_blast_results.txt"
blastn -query "$query_gene" -db "$blast_db" -out "$blast_output" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
-perc_identity "$perc_identity" -max_target_seqs 10000
echo "BLAST results saved to $blast_output"
fi
echo "BLAST analysis completed. Results saved in: $output_dir"
}

function filter_blast_results() {
local blast_result_file="$1"  # full path containing blast results
local output_dir="$2"
local coverage_threshold="$3"
local genome_type="$4"
local base_name="$(basename "$blast_result_file")"
if [[ "$genome_type" == "draft" ]]; then
local taxon="$(basename "$(dirname "$blast_result_file")")"
local iteration="$(grep -o '[0-9]\+' <<< "$base_name")"
local filtered_result="${output_dir}/${taxon}/filtered_iteration_${iteration}_draft_blast_results.txt"
mkdir -p "$(dirname "$filtered_result")"
awk -v cov="$coverage_threshold" '($13 > 0) && (($4 / $13 * 100) >= cov) {print $0}' "$blast_result_file" > "$filtered_result"
echo "Filtered BLAST results for draft genomes are saved to $filtered_result"
else
local filtered_result="$output_dir/filtered_${base_name}"
awk -v cov="$coverage_threshold" '($13 > 0) && (($4 / $13 * 100) >= cov) {print $0}' "$blast_result_file" > "$filtered_result"
echo "Filtered BLAST results for complete genomes are saved to $filtered_result"
fi
}
