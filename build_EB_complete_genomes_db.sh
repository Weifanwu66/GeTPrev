#!/bin/bash
set -euo pipefail
# ==============================================================
# self-submitting header for build_EB_complete_genomes.sh
# supports SLURM, SGE, and direct local execution with RAM check
# ==============================================================
JOB_NAME=${JOB_NAME:-build_EB_complete_genomes_db}
QUEUE=${QUEUE:-}
ACCOUNT=${ACCOUNT:-}
RUNTIME=${RUNTIME:-12:00:00}
MEMORY=${RAM:-128G}
THREADS=${CPUS:-16}

WORKDIR=$(pwd)
DATABASE_DIR="$WORKDIR/database"
GENOME_DIR="$DATABASE_DIR/complete_genomes"
BLAST_DB_DIR="$DATABASE_DIR/complete_blast_db"
FAILED_FLAG="$WORKDIR/build_EB_db_failed.flag"
MONOPHASIC_TYPHIMURIUM_LIST="$DATABASE_DIR/monophasic_Typhimurium_list.txt"
DUPLICATE_SEROTYPE_LIST="$DATABASE_DIR/duplicate_sal_serotypes.txt"
SEROTYPE_LIST_FILE="$GENOME_DIR/Salmonella/salmonella_serotype_list.txt"
SUBSPECIES_LIST_FILE="$GENOME_DIR/Salmonella/salmonella_subspecies_list.txt"
METADATA_FILE="$DATABASE_DIR/assembly_summary_bacteria.txt"
export FORCE_UPDATE="true"

if [[ -n "${SLURM_JOB_ID:-}" || -n "${JOB_ID:-}" ]]; then
echo "Running inside a job environment."
source function.sh || { echo "Error sourcing function.sh" >&2; exit 1; }
mkdir -p "$GENOME_DIR"
mkdir -p "$BLAST_DB_DIR"
> "$FAILED_FLAG"

if [[ "$FORCE_UPDATE" == "true" ]]; then
echo "[$(date '+%F %T')] Cleaning up old genome and BLAST database directories"
rm -rf "$GENOME_DIR"/*
rm -rf "$BLAST_DB_DIR"/*
fi

echo "[$(date '+%F %T')] Downloading latest assembly metadata"
download_with_retry wget -O "$METADATA_FILE" "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"

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
SPECIES_LIST_FILE="$GENOME_DIR/${GENUS}/species_list.txt"
if ! get_species_list "$GENUS" "$SPECIES_LIST_FILE"; then
echo "Failed to get species list for $GENUS" >> "$FAILED_FLAG"
fi
classify_genus_by_metadata "$GENUS"
echo "[$(date '+%F %T')] Finished processing $GENUS"
) &
genus_pids+=($!)
done

if (( ${#genus_pids[@]} > 0 )); then
echo "[$(date '+%F %T')] Waiting on genus downloads"
for pid in "${genus_pids[@]}"; do
if ! wait "$pid"; then
echo "Genus job with PID $pid failed." >> "$FAILED_FLAG"
fi
done
fi

echo "[$(date '+%F %T')] Waiting on Salmonella download"
if ! wait "$SALMONELLA_PID"; then
echo "Salmonella job failed." >> "$FAILED_FLAG"
fi

echo "[$(date '+%F %T')] Checking for duplicated genomes"
sanity_check_unique_genomes "$GENOME_DIR"

echo "[$(date '+%F %T')] Building BLAST database"
build_blastdb_for_EB_default "$GENOME_DIR" "$BLAST_DB_DIR"

generate_directory_csv "$GENOME_DIR"
echo "[$(date '+%F %T')] Default Enterobacteriaceae database build complete"
exit 0
fi

if command -v sbatch &>/dev/null && sinfo &>/dev/null; then
echo "SLURM detected. Submitting job"
cp slurm.sh slurm2.sh
sed -i "s/name/$JOB_NAME/g" slurm2.sh
[[ -n "$QUEUE" ]] && sed -i "s/queue/$QUEUE/g" slurm2.sh
sed -i "s/runtime/$RUNTIME/g" slurm2.sh
sed -i "s/RAM/$MEMORY/g" slurm2.sh
sed -i "s/hpctasks/$THREADS/g" slurm2.sh
[[ -n "$ACCOUNT" ]] && sed -i "s/account/$ACCOUNT/g" slurm2.sh
sed -i "s%vars%%g" slurm2.sh
sed -i "s%command%bash $(realpath "$0")%g" slurm2.sh
sbatch slurm2.sh
rm -f slurm2.sh
exit 0
elif command -v qsub &>/dev/null && qhost &>/dev/null; then
echo "SGE detected. Submitting job"
cp sge.sh sge2.sh
sed -i "s/name/$JOB_NAME/g" sge2.sh
[[ -n "$QUEUE" ]] && sed -i "s/queue/$QUEUE/g" sge2.sh
sed -i "s/runtime/$RUNTIME/g" sge2.sh
sed -i "s/RAM/$MEMORY/g" sge2.sh
sed -i "s/hpctasks/$THREADS/g" sge2.sh
[[ -n "$ACCOUNT" ]] && sed -i "s/account/$ACCOUNT/g" sge2.sh
sed -i "s%vars%%g" sge2.sh
sed -i "s%command%bash $(realpath "$0")%g" sge2.sh
qsub sge2.sh
rm -f sge2.sh
exit 0
else
echo "No job scheduler detected, running locally"
REQUIRED_RAM_GB=64
AVAILABLE_RAM_GB=$(free -g | awk '/^Mem:/ {print $7}')
echo "Available RAM: ${AVAILABLE_RAM_GB}GB"
if (( AVAILABLE_RAM_GB < REQUIRED_RAM_GB )); then
echo "Requires at least ${REQUIRED_RAM_GB}GB RAM. Only ${AVAILABLE_RAM_GB}GB available"
exit 1
fi
echo "Resources are sufficient. Running script locally"
exec bash "$0"
fi
