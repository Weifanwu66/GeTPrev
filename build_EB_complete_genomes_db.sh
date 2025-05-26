#!/bin/bash
#SBATCH --job-name=build_EB_complete_genomes_db
#SBATCH --output=slurm_build_EB_db_%j.out
#SBATCH --error=slurm_build_EB_db_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=72

set -euo pipefail

WORKDIR=$(pwd)
DATABASE_DIR="$WORKDIR/database"
GENOME_DIR="$DATABASE_DIR/complete_genomes"
BLAST_DB_DIR="$DATABASE_DIR/complete_blast_db"
FAILED_FLAG="$WORKDIR/build_EB_db_failed.flag"
MONOPHASIC_TYPHIMURIUM_LIST="$DATABASE_DIR/monophasic_Typhimurium_list.txt"
export FORCE_UPDATE="true"

mkdir -p "$GENOME_DIR"
mkdir -p "$BLAST_DB_DIR"
> "$FAILED_FLAG"

if [[ "$FORCE_UPDATE" == "true" ]]; then
echo "Cleaning up old genome and BLAST database directories..."
rm -rf "$GENOME_DIR"/*
rm -rf "$BLAST_DB_DIR"/*
fi

source "${WORKDIR}/function.sh"
TOTAL_CPUS=$(get_cpus)
MAX_PARALLEL_JOBS=$(( TOTAL_CPUS * 2 / 3 ))
max_parallel_genus=$(( MAX_PARALLEL_JOBS / 6 ))
max_parallel_species=$(( MAX_PARALLEL_JOBS / max_parallel_genus ))
max_parallel_serotypes=$(( MAX_PARALLEL_JOBS / 3 ))

(
echo "Processing Salmonella"
if ! download_single_genus "Salmonella" "$GENOME_DIR"; then
echo "Failed to download genus: Salmonella" >> "$FAILED_FLAG"; exit 1
fi
SALMONELLA_DIR="$GENOME_DIR/Salmonella"
SALMONELLA_SPECIES_LIST="$SALMONELLA_DIR/species_list.txt"
if ! get_species_list "Salmonella" "$SALMONELLA_SPECIES_LIST"; then
echo "Failed to get species list for Salmonella" >> "$FAILED_FLAG"; exit 1
fi

sort "$SALMONELLA_SPECIES_LIST" | uniq > "${SALMONELLA_SPECIES_LIST}.dedup"
mv "${SALMONELLA_SPECIES_LIST}.dedup" "$SALMONELLA_SPECIES_LIST"

species_pids=()
while IFS= read -r species; do
[[ -z "$species" ]] && continue
(
echo "Salmonella species: $species"
if ! download_species "$species" "$GENOME_DIR"; then
echo "Failed to download $species under Salmonella" >> "$FAILED_FLAG"
fi
) &
species_pids+=($!)
if (( ${#species_pids[@]} >= max_parallel_species )); then
wait "${species_pids[@]}"
species_pids=()
fi
done < "$SALMONELLA_SPECIES_LIST"
wait "${species_pids[@]}"

if ! get_salmonella_subsp_list "$SALMONELLA_DIR"; then
echo "Failed to get Salmonella subspecies list" >> "$FAILED_FLAG"; exit 1
fi
if ! download_salmonella_subsp "$GENOME_DIR"; then
echo "Failed to download Salmonella subspecies" >> "$FAILED_FLAG"; exit 1
fi

if ! get_salmonella_serotype_list "$SALMONELLA_DIR"; then
echo "Failed to get Salmonella serotype list" >> "$FAILED_FLAG"; exit 1
fi

serotype_pids=()
while IFS= read -r serotype; do
[[ -z "$serotype" ]] && continue
(
echo "Downloading Salmonella serotype: $serotype"
if ! download_salmonella_serotype "$GENOME_DIR" "$serotype"; then
echo "Failed to download serotype: $serotype" >> "$FAILED_FLAG"
fi
) &
serotype_pids+=($!)
if (( ${#serotype_pids[@]} >= max_parallel_serotypes )); then
wait "${serotype_pids[@]}"
serotype_pids=()
fi
done < "$SALMONELLA_DIR/serotype_list.txt"
wait "${serotype_pids[@]}"

move_unclassified_genomes "$SALMONELLA_DIR"
build_blastdb "$SALMONELLA_DIR" "$BLAST_DB_DIR"
find "$SALMONELLA_DIR" -type f -name "*_genomic.fna" -exec gzip -f {} \;
) &
SALMONELLA_PID=$!

GENUS_LIST=("Escherichia" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter")
genus_pids=()

for GENUS in "${GENUS_LIST[@]}"; do
(
echo "Processing $GENUS"
if ! download_single_genus "$GENUS" "$GENOME_DIR"; then
echo "Failed to download genus: $GENUS" >> "$FAILED_FLAG"; exit 1
fi
GENUS_DIR="$GENOME_DIR/$GENUS"
SPECIES_LIST_FILE="$GENUS_DIR/species_list.txt"
if ! get_species_list "$GENUS" "$SPECIES_LIST_FILE"; then
echo "Failed to get species list for $GENUS" >> "$FAILED_FLAG"; exit 1
fi

sort "$SPECIES_LIST_FILE" | uniq > "${SPECIES_LIST_FILE}.dedup"
mv "${SPECIES_LIST_FILE}.dedup" "$SPECIES_LIST_FILE"

echo "$GENUS: downloading species"
species_pids=()
while IFS= read -r species; do
[[ -z "$species" ]] && continue
(
echo "$GENUS: $species"
if ! download_species "$species" "$GENOME_DIR"; then
echo "Failed to download $species under $GENUS" >> "$FAILED_FLAG"
fi
) &
species_pids+=($!)
if (( ${#species_pids[@]} >= max_parallel_species )); then
wait "${species_pids[@]}"
species_pids=()
fi
done < "$SPECIES_LIST_FILE"
wait "${species_pids[@]}"

echo "$GENUS: organizing unclassified genomes"
move_unclassified_genomes "$GENUS_DIR"

echo "$GENUS: building BLAST db"
build_blastdb "$GENUS_DIR" "$BLAST_DB_DIR"

echo "$GENUS: compressing fasta files"
find "$GENUS_DIR" -type f -name "*_genomic.fna" -exec gzip -f {} \;
) &
genus_pids+=($!)

done

if (( ${#genus_pids[@]} > 0 )); then
echo "Waiting on remaining genus jobs"
for pid in "${genus_pids[@]}"; do
if ! wait "$pid"; then
echo "Genus job with PID $pid failed." >> "$FAILED_FLAG"
fi
done
fi

echo "Waiting on Salmonella job"
if ! wait "$SALMONELLA_PID"; then
echo "Salmonella job failed." >> "$FAILED_FLAG"
fi

if [[ -s "$FAILED_FLAG" ]]; then
echo "One or more genome downloads failed. See $FAILED_FLAG for details."
exit 1
fi

echo "All genome downloads, database builds, and compression steps completed."
