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
TYPHIMURIUM_LIST="$DATABASE_DIR/Typhimurium_list.txt"
DUPLICATE_SEROTYPE_LIST="$DATABASE_DIR/duplicate_sal_serotypes.txt"
SEROTYPE_LIST_FILE="$GENOME_DIR/Salmonella/salmonella_serotype_list.txt"
SUBSPECIES_LIST_FILE="$GENOME_DIR/Salmonella/salmonella_subspecies_list.txt"
METADATA_FILE="$DATABASE_DIR/assembly_summary_bacteria.txt"
export FORCE_UPDATE="true"

mkdir -p "$GENOME_DIR"
mkdir -p "$BLAST_DB_DIR"
> "$FAILED_FLAG"

if [[ "$FORCE_UPDATE" == "true" ]]; then
 echo "[$(date '+%F %T')] Cleaning up old genome and BLAST database directories..."
  rm -rf "$GENOME_DIR"/*
  rm -rf "$BLAST_DB_DIR"/*
fi

# Download latest metadata
echo "[$(date '+%F %T')] Downloading latest assembly metadata"
wget -O "$METADATA_FILE" "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"

source "${WORKDIR}/function.sh"
TOTAL_CPUS=$(get_cpus)
MAX_PARALLEL_JOBS=$(( TOTAL_CPUS * 2 / 3 ))
max_parallel_genus=$(( MAX_PARALLEL_JOBS / 6 ))

(
echo "[$(date '+%F %T')] Processing Salmonella"
if ! download_single_genus "Salmonella" "$GENOME_DIR"; then
echo "Failed to download genus: Salmonella" >> "$FAILED_FLAG"
fi
if ! get_salmonella_subsp_list "$SUBSPECIES_LIST_FILE"; then
echo "Failed to get Salmonella subspecies list" >> "$FAILED_FLAG"
fi
if ! get_salmonella_serotype_list "$SEROTYPE_LIST_FILE"; then
echo "Failed to get Salmonella serotype list" >> "$FAILED_FLAG"
fi
classify_salmonella_by_metadata "$SUBSPECIES_LIST_FILE" "$SEROTYPE_LIST_FILE"
echo "[$(date '+%F %T')] Finished processing Salmonella"
) &
SALMONELLA_PID=$!

GENUS_LIST=("Escherichia" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter")
genus_pids=()
for GENUS in "${GENUS_LIST[@]}"; do
(
echo "[$(date '+%F %T')] Processing $GENUS"
if ! download_single_genus "$GENUS" "$GENOME_DIR"; then
echo "Failed to download genus: $GENUS" >> "$FAILED_FLAG"
fi
SPECIES_LIST_FILE="$(pwd)/database/complete_genomes/${GENUS}/species_list.txt"
if ! get_species_list "$GENUS" "$SPECIES_LIST_FILE"; then
echo "Failed to get species list for $GENUS" >> "$FAILED_FLAG"
fi
classify_genus_by_metadata "$GENUS" 
echo "[$(date '+%F %T')] Finished processing $GENUS"
) &
genus_pids+=($!)
done

if (( ${#genus_pids[@]} > 0 )); then
echo "[$(date '+%F %T')] Waiting on remaining genus jobs"
for pid in "${genus_pids[@]}"; do
if ! wait "$pid"; then
echo "Genus job with PID $pid failed." >> "$FAILED_FLAG"
fi
done
fi

echo "[$(date '+%F %T')] Waiting on Salmonella job"
if ! wait "$SALMONELLA_PID"; then
echo "Salmonella job failed." >> "$FAILED_FLAG"
fi

echo "[$(date '+%F %T')] Running genome duplication sanity check"
sanity_check_unique_genomes "$GENOME_DIR"

echo "[$(date '+%F %T')] Building BLAST DB for all genus"
build_blastdb_for_EB_default "$GENOME_DIR" "$BLAST_DB_DIR"

generate_directory_csv "$GENOME_DIR" 
echo "[$(date '+%F %T')] Default database construction completed."
