#!/bin/bash
MONOPHASIC_TYPHIMURIUM_LIST="$(pwd)/database/monophasic_Typhimurium_list.txt"
TYPHIMURIUM_LIST="$(pwd)/database/Typhimurium_list.txt"
DUPLICATE_SEROTYPE_LIST="$(pwd)/database/duplicate_sal_serotypes.txt"
ASSEMBLY_LEVEL="complete"
# Get CPU count (SLURM-aware)
function get_cpus() {
if [[ -n "${SLURM_CPUS_ON_NODES:-}" ]]; then
echo "$SLURM_CPUS_ON_NODE"
else
nproc
fi
}
# Build a function to download genus-level genomes using ncbi-genome-download
# Download a single genus
function download_single_genus() {
local genus="$1"
local output_dir="$2"
local genus_dir="${output_dir}/${genus}/aggregated"
if [[ -n "$(find "$genus_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for $genus, skipping"
return
fi
mkdir -p "$genus_dir"
echo "Downloading genus: $genus"
ncbi-genome-download bacteria --genera "$genus" \
 --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$genus_dir" --verbose --flat-output
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

# Build a function to get all species names under each target genus
function get_species_list() {
local target_genus="$1"
local output_file="$2"
mkdir -p "$(dirname "$output_file")"
esearch -db taxonomy -query "${target_genus}[Subtree]" | efetch -format xml | xtract -pattern Taxon -element ScientificName | \
awk 'NF == 2 && $1 ~ /^[A-Z][a-z]+$/ && $2 ~ /^[a-z]+$/ {print $0}' > "$output_file"
echo "Saved all species names of $target_genus to $output_file"
}

# Build a function to download species-level genomes using ncbi-genome-download
function download_species() {
local target_input="$1"
local output_dir="$2"
local total_cpus=$(get_cpus)
local cores_half=$(( total_cpus / 2 ))
local max_parallel_species=$(( cores_half < 1 ? 1 : (cores_half > 4 ? 4 : cores_half) ))
download_single_species() {
local species="$1"
local genus=$(echo "$species" | awk '{print $1}')
local clean_species=$(echo "$species" | sed 's/ /_/g')
species_dir="${output_dir}/${genus}/${clean_species}"
[[ "$species" == "Salmonella enterica" ]] && species_dir="$species_dir/aggregated"
mkdir -p "$species_dir"
if [[ -n "$(find "$species_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for $species, skipping donwloading"
return
fi
echo "Downloading genomes for $species"
ncbi-genome-download bacteria --genera "$species" \
 --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$species_dir" --verbose --flat-output
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

# Build a function to get all subspecies names under Salmonella enterica
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
local subspecies_dir="$GENOME_DIR/Salmonella/Salmonella_enterica/$subspecies"
[[ "$subspecies" == "enterica" ]] && subspecies_dir="$subspecies_dir/aggregated"
mkdir -p "$subspecies_dir"
if [[ -n "$(find "$subspecies_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for *Salmonella enterica* subsp. $subspecies"
continue
fi
echo "Downloading genomes for Salmonella enterica subsp. $subspecies"
ncbi-genome-download bacteria --genera "Salmonella enterica subsp. $subspecies" \
 --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$subspecies_dir" --verbose --flat-output
find "$subspecies_dir" -type f -name "*_genomic.fna.gz" -exec gzip -d {} \;
if [[ -z "$(find "$subspecies_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "No genomes found for $subspecies. Remove empty directory."
rm -rf "$subspecies_dir"
fi
done < "$subspecies_list_file"
}

# Build a function to download Salmonella serotype
function get_salmonella_serotype_list() {
local output_file="$1"
local output_dir
output_dir=$(dirname "$output_file")
mkdir -p "$output_dir"
esearch -db taxonomy -query "Salmonella enterica subsp. enterica[Subtree]" | efetch -format xml | xtract -pattern Taxon -element ScientificName | grep -v -E "str\.|var\." | sort -u | grep -v "Salmonella enterica subsp. enterica$" > "$output_file"
echo "Saved all serotype names in $output_file"
# Remove all Salmonella enterica subsp. enterica serovar prefix
sed -i 's/Salmonella enterica subsp. enterica serovar //g' "$output_file"
# Delete known serotype names names by name match that are already stored in monophasic_Typhimurium_list.txt and Typhimurium_list.txt
grep -vFf "$DUPLICATE_SEROTYPE_LIST" "$output_file" > tmp && \
mv tmp "$output_file"
} 

function download_salmonella_serotype() {
local serotype_list_file="$1"
local total_cpus=$(get_cpus)
local max_parallel_serotypes=$(( total_cpus / 4 ))
(( max_parallel_serotypes < 1 )) && max_parallel_serotypes=1
(( max_parallel_serotypes > 10 )) && max_parallel_serotypes=10
echo "Downloading Salmonella enterica subsp. enterica serotype complete genomes"
local monophasic_dir="$GENOME_DIR/Salmonella/Salmonella_enterica/enterica/monophasic_Typhimurium"
mkdir -p "$monophasic_dir"
local typhimurium_dir="$GENOME_DIR/Salmonella/Salmonella_enterica/enterica/Typhimurium"
mkdir -p "$typhimurium_dir"
download_single_serotype() {
local serotype="$1"
local serotype_dir="${GENOME_DIR}/Salmonella/Salmonella_enterica/enterica/${serotype}"
mkdir -p "$serotype_dir"
if [[ -n "$(find "$serotype_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for Salmonella enterica $serotype, skipping."
return
fi
echo "Downloading genomes for Salmonella $serotype"
ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $serotype" \
 --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$serotype_dir" --verbose --flat-output
find "$serotype_dir" -type f -name "*_genomic.fna.gz" -exec gzip -d {} \;
if [[ -z "$(find "$serotype_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "No genomes found for $serotype. Remove empty directory."
rm -rf "$serotype_dir"
fi
}
while read -r serotype; do
[[ -z "$serotype" ]] && continue
(download_single_serotype "$serotype" ) &
while (( $(jobs -r | wc -l) >= max_parallel_serotypes )); do sleep 1; done
done < "$serotype_list_file"
wait
echo "Downloading all monophasic Typhimurium genomes"
while read -r monophasic; do
[[ -z "$monophasic" ]] && continue
ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $monophasic" \
 --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$monophasic_dir" --verbose --flat-output
find "$monophasic_dir" -type f -name "*_genomic.fna.gz" -exec gzip -d {} \;
done < "$MONOPHASIC_TYPHIMURIUM_LIST"
echo "Downloading all monophasic Typhimurium genomes"
while read -r typhimurium; do
[[ -z "$typhimurium" ]] && continue
ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $typhimurium" \                                                           --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$monophasic_dir" --verbose --flat-output
find "$typhimurium_dir" -type f -name "*_genomic.fna.gz" -exec gzip -d {} \;
done < "$TYPHIMURIUM_LIST"
}

function move_unclassified_genomes() {
local output_dir="$1"
echo "Scanning for aggregated directories under: $output_dir"
mapfile -t aggregated_dirs < <(find "$output_dir" -type d -name "aggregated")
for aggregated_dir in "${aggregated_dirs[@]}"; do
echo "Found aggregated dir: $aggregated_dir"
local level_dir
level_dir=$(dirname "$aggregated_dir")
local unclassified_dir="${level_dir}/unclassified"
mkdir -p "$unclassified_dir"
declare -A classified_files
while IFS= read -r -d '' genome_path; do
genome_basename=$(basename "$genome_path" | tr -d '\r\n ')
classified_files["$genome_basename"]=1
done < <(find "$level_dir" -mindepth 1 -maxdepth 1 -type d \
  ! \( -name "aggregated" -o -name "unclassified" \) \
  -exec find {} -type f -name "*_genomic.fna" -print0 \;)
while IFS= read -r -d '' genome_file; do
genome_basename=$(basename "$genome_file" | tr -d '\r\n ')
if [[ -z "${classified_files["$genome_basename"]+exists}" ]]; then
echo "Moving unclassified genome: $genome_basename"
mv -f "$genome_file" "$unclassified_dir/" || echo "Failed to move $genome_file"
fi
done < <(find "$aggregated_dir" -type f -name "*_genomic.fna" -print0)
rm -rf "$aggregated_dir"
echo "Removed: $aggregated_dir"
done
echo "Finished organizing unclassified genomes in: $output_dir"
}

function create_blastdb() {
local INPUT_FASTA="$1"
local DB_PATH="$2"
makeblastdb -in "$INPUT_FASTA" -dbtype nucl -out "${DB_PATH}"
echo "BLAST database created for ${DB_PATH}"
}

function build_blastdb() {
local input_dir="$1"
local output_dir="$2"
mkdir -p "$output_dir"
echo "Building BLAST database from genomes in $input_dir"
# Create BLAST database for all genus
find "$input_dir" -mindepth 1 -maxdepth 1 -type d | while read -r GENUS_DIR; do
GENUS=$(basename "$GENUS_DIR")
echo "Processing genus: $GENUS"
GENUS_FASTA="${GENUS_DIR}/${GENUS}_all_genomes.fna"
find "$GENUS_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$GENUS_FASTA"
create_blastdb "$GENUS_FASTA" "${output_dir}/${GENUS}"
# Species-level databases
find "$GENUS_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SPECIES_DIR; do
SPECIES=$(basename "$SPECIES_DIR")
[[ "$SPECIES" == "unclassified" ]] && continue
SPECIES_FASTA="${SPECIES_DIR}/${SPECIES}_all_genomes.fna"
find "$SPECIES_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SPECIES_FASTA"
create_blastdb "$SPECIES_FASTA" "${output_dir}/${SPECIES}"
if [[ "$SPECIES" == "Salmonella_enterica" ]]; then
find "$SPECIES_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SUBSPECIES_DIR; do
SUBSPECIES=$(basename "$SUBSPECIES_DIR")
[[ "$SUBSPECIES" == "unclassified" ]] && continue
SUBSPECIES_FASTA="${SUBSPECIES_DIR}/Salmonella_enterica_subsp_${SUBSPECIES}_all_genomes.fna"
find "$SUBSPECIES_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SUBSPECIES_FASTA"
create_blastdb "$SUBSPECIES_FASTA" "${output_dir}/Salmonella_enterica_subsp_${SUBSPECIES}"
if [[ "$SUBSPECIES" == "enterica" ]]; then
find "$SUBSPECIES_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SEROTYPE_DIR; do
SEROTYPE=$(basename "$SEROTYPE_DIR")
[[ "$SEROTYPE" == "unclassified" ]] && continue
SEROTYPE_FASTA="${SEROTYPE_DIR}/${SEROTYPE}_genomes.fna"
find "$SEROTYPE_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SEROTYPE_FASTA"
create_blastdb "$SEROTYPE_FASTA" "${output_dir}/Salmonella_${SEROTYPE}"
done
fi
done
fi
done
done
echo "Finished building BLAST databases"
}

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

function get_total_genomes_count() {
local input="$1"
local assembly_level="$2"
local total_genomes=0
local query="$(extract_taxon_info "$input")"
if [[ "$query" == "Salmonella enterica subsp. enterica serovar monophasic Typhimurium" ]]; then
while read -r actual_name; do
count=$(( $(ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $actual_name" --assembly-level "$assembly_level" --section genbank --dry-run  | wc -l) - 1 ))
((total_genomes += count))
done < "$MONOPHASIC_TYPHIMURIUM_LIST"
else
total_genomes=$(( $(ncbi-genome-download bacteria --genera "$query" --assembly-level "$assembly_level" --section genbank --dry-run --verbose | wc -l) - 1 ))
total_genomes=$(( total_genomes < 0 ? 0 : total_genomes ))
fi
echo "$total_genomes"
}

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
# Finite population correction (FPC)
local sample_size=$(echo "scale=6; ($n0 * $total_genomes) / ($total_genomes + $n0 -1)" | bc -l)
sample_size=$(echo "$sample_size" | awk '{printf "%.0f", $1}')
# Compute number of iterations using square-root scaling
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
accessions+=$(ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $actual_name" --assembly-level contig --section genbank --dry-run | tail -n +2 | awk -F '/' '{print $NF}' | awk '{print $1}')
done < "$MONOPHASIC_TYPHIMURIUM_LIST"
else
accessions=$(ncbi-genome-download bacteria --genera "$query" --assembly-level contig --section genbank --dry-run | tail -n +2 | awk -F '/' '{print $NF}' | awk '{print $1}')
fi
valid_accessions=$(echo "$accessions" | grep -E '^GCA_[0-9]+\.[0-9]+$' | shuf -n "$sample_size")
echo "$valid_accessions" > "${iteration_dir}/selected_accessions.txt"
ncbi-genome-download bacteria --assembly-accessions "$iteration_dir/selected_accessions.txt" --formats fasta --assembly-level contig --section genbank  --output-folder "$iteration_dir"
find "$iteration_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$0" && mv "${0%.gz}" "'"$iteration_dir"'"' {} \;
rm -rf "$iteration_dir/genbank"
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
local blast_output="${output_dir}/${taxon}/iteration_${iteration}_draft_blast_results.txt"
mkdir -p "$(dirname "$blast_output")"
echo "Running BLAST for $taxon..."
blastn -query "$query_gene" -db "$blast_db" -out "$blast_output" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
-perc_identity "$perc_identity" -max_target_seqs 7000
echo "BLAST results saved to: $blast_output"
else
# Skip empty lines
[[ -z "$taxon" ]] && return
local blast_db_name="${taxon// /_}"
blast_db_name="${blast_db_name//./}"
local blast_db="${BLAST_DB_DIR}/${blast_db_name}"
local clean_taxon="$(extract_taxon_info "$taxon")"
local blast_output_name="${clean_taxon// /_}"
blast_output_name="${blast_output_name//./}"
local blast_output="${output_dir}/${blast_output_name}_complete_blast_results.txt"
blastn -query "$query_gene" -db "$blast_db" -out "$blast_output" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
-perc_identity "$perc_identity" -max_target_seqs 7000
echo "BLAST results saved to $blast_output"
fi
echo "BLAST analysis completed. Results saved in: $output_dir"
}

function filter_blast_results() {
local blast_result_file="$1"  # Full path containing blast results
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
