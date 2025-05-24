#!/bin/bash
#SBATCH --job-name=build_EB_complete_genomes_db
#SBATCH --output=slurm_build_EB_db_%j.out
#SBATCH --error=slurm_build_EB_db_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
set -euo pipefail
# ============================
# Set working directories and paths
# ============================
WORKDIR=$(pwd)
DATABASE_DIR="$WORKDIR/database"
GENOME_DIR="$DATABASE_DIR/complete_genomes"
BLAST_DB_DIR="$DATABASE_DIR/complete_blast_db"
FAILED_FLAG="$WORKDIR/build_EB_db_failed.flag"
MONOPHASIC_TYPHIMURIUM_LIST="$DATABASE_DIR/monophasic_Typhimurium_list.txt"
mkdir -p "$GENOME_DIR"
mkdir -p "$BLAST_DB_DIR"
> "$FAILED_FLAG"

source "${WORKDIR}/function.sh"
TOTAL_CPUS=$(get_cpus)
# ============================
# Step 1: Download genus and species in parallel
# ============================
echo "Starting genus and species download in parallelization"
printf "%s\n" "Escherichia" "Salmonella" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter" > "$GENOME_DIR/genus_list.txt"
max_parallel_genus=6
max_parallel_species=4
running_genus_jobs=0
while read -r GENUS; do
(
if ! download_single_genus "$GENUS" "$GENOME_DIR"; then
echo "Failed to download genus: $GENUS" >> "$FAILED_FLAG"
exit 1
fi
temp_species_list="$GENOME_DIR/species_list_${GENUS}.txt"
if ! get_species_list "$GENUS" "$temp_species_list"; then
echo "Failed to get species list for: $GENUS" >> "$FAILED_FLAG"
exit 1
fi
sort "$temp_species_list" | uniq > "${temp_species_list}.dedup"
mv "${temp_species_list}.dedup" "$temp_species_list"
species_jobs=0
while IFS= read -r species; do
[[ -z "$species" ]] && continue
(
if ! download_species "$species" "$GENOME_DIR"; then
echo "Failed to download $species under genus: $GENUS" >> "$FAILED_FLAG"
fi
) &
((species_jobs++))
if (( species_jobs >= max_parallel_species )); then
wait
species_jobs=0
fi
done < "$temp_species_list"
wait
) &
((running_genus_jobs++))
if (( running_genus_jobs >= max_parallel_genus )); then
wait
running_genus_jobs=0
fi
done < "$GENOME_DIR/genus_list.txt"
wait
cat "$GENOME_DIR"/species_list_*.txt | sort | uniq > "$GENOME_DIR/species_list.txt"
rm -f "$GENOME_DIR"/species_list_*.txt "$GENOME_DIR/genus_list.txt"
# ============================
# Step 2: Download Salmonella subspecies (sequential) and serotypes in parallel
# ============================
echo "Starting Salmonella subspecies and serotype download"
if ! get_salmonella_subsp_list "$GENOME_DIR"; then
echo "Failed to get Salmonella subspecies list" >> "$FAILED_FLAG"
exit 1
fi

if ! download_salmonella_subsp "$GENOME_DIR"; then
echo "Failed to download Salmonella subspecies" >> "$FAILED_FLAG"
exit 1
fi

if ! get_salmonella_serotype_list "$GENOME_DIR"; then
echo "Failed to get Salmonella serotype list" >> "$FAILED_FLAG"
exit 1
fi

if ! download_salmonella_serotype "$GENOME_DIR"; then
echo "Failed to download Salmonella serotypes" >> "$FAILED_FLAG"
exit 1
fi

# ============================
# Step 3: Fail-safe check before continuing
# ============================
if [[ -s "$FAILED_FLAG" ]]; then
echo "One or more genome downloads failed. See $FAILED_FLAG for details."
echo "BLAST database will NOT be built."
exit 1
fi

# ============================
# Step 4: Organize unclassified genomes and build BLAST database
# ============================
echo "Organizing unclassified genomes"
move_unclassified_genomes "$GENOME_DIR"

echo "Building BLAST databases"
build_blastdb "$GENOME_DIR" "$BLAST_DB_DIR"

echo "Genome download and BLAST database build completed successfully."
