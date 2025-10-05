#!/bin/bash
# Set initial variables and paths
MODE="light"
MIN_COVERAGE=80
MIN_IDENTITY=90
WORKDIR=$(pwd)
TAXON_FILE=""
GENE_FILE=""
DOWNLOAD_FILE=""
FAILED_FLAG="${WORKDIR}/build_custom_db_failed.flag"
> "$FAILED_FLAG"
OUTPUT_FILE="${WORKDIR}/result/gene_summary.csv"
DATABASE_DIR="${WORKDIR}/database"
BLAST_DB_DIR="${DATABASE_DIR}/complete_blast_db"
CUSTOM_GENOMES_DIR="${DATABASE_DIR}/complete_genomes"
BLAST_RESULT_DIR="${WORKDIR}/result/complete_blast_results"
FILTERED_BLAST_RESULT_DIR="${WORKDIR}/result/filtered_complete_blast_results"
DRAFT_GENOMES_DIR="${DATABASE_DIR}/draft_genomes"
DRAFT_BLAST_DB_DIR="${DATABASE_DIR}/draft_blast_db"
DRAFT_BLAST_RESULT_DIR="${WORKDIR}/result/draft_blast_results"
FILTERED_DRAFT_BLAST_RESULT_DIR="${WORKDIR}/result/filtered_draft_blast_results"
OVERWRITE=false
FORCE_REBUILD=false
GET_ALL_SPECIES=false
# job scheduler variables
runtime=24:00:00; hpcmem=360GB; hpcthreads=72; hpc=F; queue=NA; account=NA
# usage setup
usage() {
    echo "Usage: $0 -g GENE_FILE [-t TAXON_FILE] [-d DOWNLOAD_FILE] [-c COVERAGE] [-i IDENTITY] [-o OUTPUT_FILE] [-p HPC_CLUSTER] [-q QUEUE] [-r RUNTIME] [-m MEMORY] [-C THREADS] [-a ACCOUNT] [-H MODE] [-O OVERWRITE] [-F FORCE_REBUILD] [--get-all-species] [-h|--help]"
    echo ""
    echo "Required arguments:"
    echo "-g GENE_FILE      : FASTA file containing target gene sequences."
    echo ""
    echo "Optional arguments:"
    echo "-t TAXON_FILE     : File containing target species (one per line) or a single taxon name."
    echo "-d DOWNLOAD_FILE  : File listing target genera (one per line) to download complete genomes."
    echo "-c COVERAGE       : Minimum genome coverage threshold (default: 80%)."
    echo "-i IDENTITY       : Minimum percentage identity threshold (default: 90%)."
    echo "-o OUTPUT_FILE    : Output filename for results (default: gene_summary.tsv)."
    echo "-p HPC_CLUSTER    : HPC system name (optional; future compatibility)."
    echo "-q QUEUE          : Queue/partition name (e.g., ceres, short, long)."
    echo "-r RUNTIME        : SLURM walltime request (e.g., 04:00:00)."
    echo "-m MEMORY         : Memory request for SLURM job (e.g., 16G)."
    echo "-C THREADS        : Number of CPU cores to request for SLURM."
    echo "-a ACCOUNT        : SLURM account/project (if needed)."
    echo "-H MODE           : Analysis mode ('light' or 'heavy'). Default is light."
    echo "-O OVERWRITE      : Set to true to overwrite previous results (default: false)."
    echo "-F FORCE_REBUILD  : Set to true to rebuild previous custom complete genomes database (default: false)."
    echo "--get-all-species : Automatically expand genus-level input to species classification if genus only during custom targets database construction (default: false)." 
    echo "-h, --help        : Show this help message and exit."
}
# Parse argument (adapted)
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -g) GENE_FILE="$2"; shift 2 ;;
        -t) TAXON_FILE="$2"; shift 2 ;;
        -d) DOWNLOAD_FILE="$2"; shift 2 ;;
        -c) MIN_COVERAGE="$2"; shift 2 ;;
        -i) MIN_IDENTITY="$2"; shift 2 ;;
        -o) OUTPUT_FILE="$2"; shift 2 ;;
        -p) hpc="$2"; shift 2 ;;
        -q) queue="$2"; shift 2 ;;
        -r) runtime="$2"; shift 2 ;;
        -m) hpcmem="$2"; shift 2 ;;
        -C) hpcthreads="$2"; shift 2 ;;
        -a) account="$2"; shift 2 ;;
        -H) MODE="$2"; shift 2 ;;
        -O) OVERWRITE="$2"; shift 2 ;;
        -F) FORCE_REBUILD="$2"; shift 2 ;;
        --get-all-species) GET_ALL_SPECIES=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Invalid option: $1"; usage; exit 1 ;;
    esac
done
# convert key file paths to absolute form
GENE_FILE="$(readlink -f "$GENE_FILE")"
[[ -n "$TAXON_FILE" && -f "$TAXON_FILE" ]] && TAXON_FILE="$(readlink -f "$TAXON_FILE")"
[[ -n "$DOWNLOAD_FILE" && -f "$DOWNLOAD_FILE" ]] && DOWNLOAD_FILE="$(readlink -f "$DOWNLOAD_FILE")"
OUTPUT_FILE="$(readlink -f "$OUTPUT_FILE" || echo "$OUTPUT_FILE")"
source "${WORKDIR}/function.sh" || { echo "Error sourcing function.sh" >&2; exit 1; }
check_dependencies
#adapted from GEAbash_v1.0.0; seems to be working as expected
#while getopts ':g:t::c::i::o::p::q::r::m::C::a::H::O::h::' flag; do
#  case "${flag}" in
#    g) GENE_FILE="${OPTARG}" ;;    t) TAXON_FILE="${OPTARG}" ;;    c) MIN_COVERAGE="${OPTARG}" ;;
#    i) MIN_IDENTITY="${OPTARG}" ;;    o) OUTPUT_FILE="${OPTARG}" ;;    p) hpc="${OPTARG}" ;;
#    q) queue="${OPTARG}" ;;    r) runtime="${OPTARG}" ;;    m) hpcmem="${OPTARG}" ;;
#    C) hpcthreads="${OPTARG}" ;;    a) account="${OPTARG}" ;;    H) MODE="${OPTARG}" ;;
#    O) OVERWRITE="${OPTARG}" ;;    F) FORCE_REBUILD="${OPTARG}"
#    -h|--help) usage; exit 0 ;;
#    *) echo "Invalid option: $1"; usage
#       exit 1 ;;    esac; done


while [ $hpc == F ]; do

# detect and utilize slurm or sge job manager. tested on atlas/ceres (slurm) and hank (sge) hpcs and putatively any other systems as requested by users.
# progresses to an expected `find` error line 292
# suspected dependencies: 1) conda (unless some testing/development done independent). e.g. the batch script slurm template in GEAbash calls 
# module load apptainer and I expect this slurm template will ultimately need to call module load miniconda (ceres) or module load miniconda3 (atlas)
# 2) The batch template will need to be in the same directory as EGP.sh and both will need to be in $(pwd)
# more testing to be done after have THE database.
if [ ! $queue == "NA" ]; then 

#detect job submission system
if [ -n "$(sinfo 2>/dev/null)" ]; then submsys=slurm
#detect sge
elif [ -n "$(qhost 2>/dev/null)" ]; then submsys=sge
#exit if queue specified but no job scheduler detected
else echo -n "queue specified but hpc system "
echo "uncertain. exiting"; exit 1; fi

#if sge, remove any logs from previous runs
if [[ " ${submsys} " = " sge " ]]; then 
if [ -f sge2.sh ]; then rm sge2_getprev.* >/dev/null 2>&1; sleep 60; fi
#then go get a new sge template file to use for THIS getprev run
cp sge.sh sge2.sh
#if slurm, remove any logs from previous runs
elif [[ " ${submsys} " = " slurm " ]]; then
if [ -f slurm2.sh ]; then rm slurm2*.out >/dev/null 2>&1; sleep 60; fi
#then go get a new slurm template file to use for THIS GEAbash run
cp slurm.sh slurm2.sh
#otherwise exit. this line should be unnecessary
else echo "system uncertain. exiting"; exit 1
fi

###if sge...
if [[ " ${submsys} " = " sge " ]]
#alert the user
then echo -n "preparing to run GeTPrev in hpc cluster mode. "
echo -n "GeTPrev log outputs will be in the hpc submission system "
echo "log files for sge. e.g., getprev.e* & getprev.o*"
#edit the sge template with local variables and user arguments
sed -i "s/name/sge2_getprev/g" sge2.sh #local
sed -i "s/queue/$queue/g" sge2.sh #user
sed -i "s/runtime/$runtime/g" sge2.sh #user
sed -i "s/RAM/$hpcmem/g" sge2.sh #user
sed -i "s/hpctasks/$hpcthreads/g" sge2.sh #user
#write the getprev command to the sge template
#to use this design syntax, all user options need single dash shortcuts
#e.g --mode heavy and --overwrite become -H heavy and -O true respectively.
#lines 201-205 assume -m <RAM> -q <queue> -r <runtime> in 
#user supplied arguments to getprev or getprev defaults
#lines 214,221 assume -g -c -i -o -H -O -t in 
#user supplied arguments to getprev or getprev defaults
if [[ "$GET_ALL_SPECIES" == true ]]; then
sed -i \
's%command%bash getprev.sh -g "$GENE_FILE" -c "$MIN_COVERAGE" -i "$MIN_IDENTITY" -o "$OUTPUT_FILE" -H "$MODE" -O "$OVERWRITE" -t "$TAXON_FILE" -d "$DOWNLOAD_FILE" -F "$FORCE_REBUILD" --get-all-species -p T%g' \
sge2.sh
else
sed -i \
's%command%bash getprev.sh -g "$GENE_FILE" -c "$MIN_COVERAGE" -i "$MIN_IDENTITY" -o "$OUTPUT_FILE" -H "$MODE" -O "$OVERWRITE" -t "$TAXON_FILE" -d "$DOWNLOAD_FILE" -F "$FORCE_REBUILD" -p T%g' \
sge2.sh
fi
#write the needed variables to the sge template
#the last variable tells EGP that ${"hpc"} = T
#so that this entire while loop will be skipped
#by the EGP resubmission
sed -i \
"s%vars%OVERWRITE='$OVERWRITE'; GENE_FILE='$GENE_FILE'; MIN_COVERAGE='$MIN_COVERAGE'; MIN_IDENTITY='$MIN_IDENTITY'; OUTPUT_FILE='$OUTPUT_FILE'; MODE='$MODE'; TAXON_FILE='$TAXON_FILE'; DOWNLOAD_FILE='$DOWNLOAD_FILE'; FORCE_REBUILD='$FORCE_REBUILD'; GET_ALL_SPECIES='$GET_ALL_SPECIES'%g" \
sge2.sh

#make sure the user provided account is written to the template
if [ $account == "NA" ]
then other='##'
sed -i "s/account/$other/g" sge2.sh
sed -i "s/-P #/###/g" sge2.sh
else sed -i "s/account/$account/g" sge2.sh #an optional user supplied variable
fi
#submit the sge batch job containing the EGP
#command from line 214
qsub sge2.sh
#exit normally without error
exit 0

###if slurm...
else echo -n "preparing to run GeTPrev in hpc cluster mode. "
echo -n "GeTPrev log outputs will be in the hpc submission system "
echo "log files for slurm. e.g., slurm2-*.out"
#edit the slurm template with local variables and user arguments
sed -i "s/name/slurm2_getprev/g" slurm2.sh #local
sed -i "s/queue/$queue/g" slurm2.sh #user
sed -i "s/runtime/$runtime/g" slurm2.sh #user
sed -i "s/RAM/$hpcmem/g" slurm2.sh #user
sed -i "s/hpctasks/$hpcthreads/g" slurm2.sh #user
sed -i 's/other/$other/g' slurm2.sh #local.
#write the GEA command to the slurm template
if [[ "$GET_ALL_SPECIES" == true ]]; then
sed -i \
's%command%bash getprev.sh -g "$GENE_FILE" -c "$MIN_COVERAGE" -i "$MIN_IDENTITY" -o "$OUTPUT_FILE" -H "$MODE" -O "$OVERWRITE" -t "$TAXON_FILE" -d "$DOWNLOAD_FILE" -F "$FORCE_REBUILD" --get-all-species -p T%g' \
slurm2.sh
else
sed -i \
's%command%bash getprev.sh -g "$GENE_FILE" -c "$MIN_COVERAGE" -i "$MIN_IDENTITY" -o "$OUTPUT_FILE" -H "$MODE" -O "$OVERWRITE" -t "$TAXON_FILE" -d "$DOWNLOAD_FILE" -F "$FORCE_REBUILD" -p T%g' \
slurm2.sh
fi
#write the needed variables to the slurm template
#the last variable tells EGP that ${"hpc"} = T
#so that this entire while loop will be skipped
#by the EGP resubmission
sed -i \
"s%vars%OVERWRITE='$OVERWRITE'; GENE_FILE='$GENE_FILE'; MIN_COVERAGE='$MIN_COVERAGE'; MIN_IDENTITY='$MIN_IDENTITY'; OUTPUT_FILE='$OUTPUT_FILE'; MODE='$MODE'; TAXON_FILE='$TAXON_FILE'; DOWNLOAD_FILE='$DOWNLOAD_FILE'; FORCE_REBUILD='$FORCE_REBUILD'; GET_ALL_SPECIES='$GET_ALL_SPECIES'%g" \
slurm2.sh

#make sure the user provided account is written to the template
if [ $account == "NA" ]
then other='##'
sed -i "s/account/$other/g" slurm2.sh
sed -i "s/-A #/###/g" slurm2.sh
else sed -i "s/account/$account/g" slurm2.sh
fi
#submit the slurm batch job containing the EGP
#command from line 250
sbatch slurm2.sh
#exit normally without error
exit 0
#finish the if statement started line 195.
fi
#finish the if statement started line 170.
fi
#if the user has NOT specified a queue then prepare
#to run GeTPrev normally
if [ $queue == "NA" ]
then echo "GeTPrev continuing without job submission system"
#if a job submission system is detected,
#alert the user and exit with error.
if [ -n "$(sinfo 2>/dev/null)" ]
then echo -n "queue not specified for running "
echo "in hpc mode. exiting"; exit 1
elif [ -n "$(qhost 2>/dev/null)" ]
then echo -n "queue not specified for running "
echo "in hpc mode. exiting"; exit 1
fi
#finish the if statement started line 161
fi
#if made it this far, set hpc to T to break while loop.
hpc=T
done

# determine number of threads for blastn operations
if [[ "$hpc" == "T" ]]; then
BLAST_THREADS="$hpcthreads"
else
BLAST_THREADS="$(get_cpus)"
fi

if [[ "$OVERWRITE" == true ]]; then
rm -rf "$BLAST_RESULT_DIR" "$FILTERED_BLAST_RESULT_DIR" "$DRAFT_BLAST_RESULT_DIR" "$DRAFT_BLAST_DB_DIR" "$FILTERED_BLAST_RESULT_DIR" "$DRAFT_GENOMES_DIR" "$FILTERED_DRAFT_BLAST_RESULT_DIR" "$OUTPUT_FILE"
fi
# Ensure gene file is provided
if [[ -z "$GENE_FILE" ]]; then
echo "Error! Please provide a gene sequence file!"
usage
exit 1
fi

# Format the fasta file if the gene file ends with cds_from_genomic.fna
if [[ "$GENE_FILE" == "*cds_from_genomic.fna" ]]; then
echo "Detected NCBI CDS fasta file. Reformatting headers.."
formatted_fasta="${WORKDIR}/formatted_gene.fasta"
normalize_fasta_headers "$GENE_FILE" "$formatted_fasta"
GENE_FILE="$formatted_fasta"
fi

# Ensure if heavy mode is enabled, a taxon file must be provided
if [[ "$MODE" == "heavy" && -z "$TAXON_FILE" && -z "$DOWNLOAD_FILE" ]]; then
echo "Error: A taxon file is required in heavy mode using default database."
exit 1
fi

if [[ -z "$DOWNLOAD_FILE" ]]; then
if [[ ! -d "$BLAST_DB_DIR" || -z $(ls -A "$BLAST_DB_DIR"/*.nsq 2>/dev/null) ]]; then
echo "Error: No complete genome BLAST database detected in $BLAST_DB_DIR."
echo "Please download or build the default BLAST database first, or include -d to specify a custom database."
exit 1
fi
fi

# Output directory setup
mkdir -p "$BLAST_RESULT_DIR" "$FILTERED_BLAST_RESULT_DIR"
# Download the genbank assembly bacteria metadata
METADATA_FILE="$DATABASE_DIR/assembly_summary_bacteria.txt"
echo "Downloading latest assembly metadata"
download_with_retry wget -q -O "$METADATA_FILE" "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt" 

# if a single value is provided for -d instead of a file, create a temporary file
if [[ -n "$DOWNLOAD_FILE" && ! -f "$DOWNLOAD_FILE" ]]; then
tmp_download_file=$(mktemp)
echo "$DOWNLOAD_FILE" > "$tmp_download_file"
DOWNLOAD_FILE="$tmp_download_file"
trap "rm -f '$tmp_download_file'" EXIT
fi

# Handle custom download if provided
if [[ -n "${DOWNLOAD_FILE:-}" ]]; then
echo "Custom panel download requested"
GENOME_DIR="$CUSTOM_GENOMES_DIR"
CUSTOM_PANEL_CHECKPOINT="$GENOME_DIR/.custom_download_complete"
mkdir -p "$GENOME_DIR"
mkdir -p "$BLAST_DB_DIR"
> "$FAILED_FLAG"

if [[ -f "$CUSTOM_PANEL_CHECKPOINT" && "${FORCE_REBUILD,,}" != "true" ]]; then
echo "Custom panel database already downloaded. Skipping rebuild."
else
if [[ "${FORCE_REBUILD,,}" == "true" ]]; then
echo "FORCE_REBUILD is true. Removing previous custom genomes and BLAST DB..."
rm -rf "$GENOME_DIR"/*
rm -rf "$BLAST_DB_DIR"/*
rm -rf "$DRAFT_GENOMES_DIR"/*
rm -rf "$DRAFT_BLAST_DB_DIR"/*
fi

TOTAL_CPUS=$(get_cpus)
MAX_PARALLEL_JOBS=$(( TOTAL_CPUS * 2 / 3 ))

if [[ "$GET_ALL_SPECIES" == true ]]; then
> "$GENOME_DIR/expanded_species_list.txt"
cp "$DOWNLOAD_FILE" "$GENOME_DIR/expanded_species_list.txt"
fi

mapfile -t taxa < "$DOWNLOAD_FILE"
for taxon in "${taxa[@]}"; do
taxon=$(echo "$taxon" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | tr -d '\r')
[[ -z "$taxon" ]] && continue
(
if [[ "$taxon" =~ ^Salmonella\ enterica\ subsp\.?[[:space:]]+([a-zA-Z0-9_]+)$ ]]; then
subsp="${BASH_REMATCH[1]}"
echo "Detected Salmonella enterica subspecies: $subsp"
tmp_subsp_file=$(mktemp "${GENOME_DIR}/temp_subsp_list.XXXXXX")
echo "$subsp" > "$tmp_subsp_file"
download_salmonella_subsp "$tmp_subsp_file"
rm -f "$tmp_subsp_file"
elif [[ "$taxon" =~ ^Salmonella\ enterica\ subsp\.?\ enterica\ serovar[[:space:]]+([A-Za-z0-9_]+)$ ]]; then
serotype="${BASH_REMATCH[1]}"
echo "Detected Salmonella serotype: $serotype"
tmp_sero_file=$(mktemp "${GENOME_DIR}/temp_serotype_list.XXXXXX")
echo "$serotype" > "$tmp_sero_file"
download_salmonella_serotype "$tmp_sero_file"
rm -f "$tmp_sero_file"
# Serotype detection
elif [[ "$taxon" =~ ^Salmonella\ ([A-Z][a-zA-Z0-9_]+)$ ]]; then
serotype="${BASH_REMATCH[1]}"
echo "Detected Salmonella serotype (shorthand): $serotype"
tmp_sero_file=$(mktemp "${GENOME_DIR}/temp_serotype_list.XXXXXX")
echo "$serotype" > "$tmp_sero_file"
download_salmonella_serotype "$tmp_sero_file"
rm -f "$tmp_sero_file"
elif [[ "$taxon" =~ ^Salmonella[[:space:]]+[Mm]onophasic[[:space:]]+Typhimurium$ ]]; then
serotype="monophasic Typhimurium"
echo "Detected Salmonella serotype (monophasic shorthand): $serotype"
tmp_sero_file=$(mktemp "${GENOME_DIR}/temp_serotype_list.XXXXXX")
echo "$serotype" > "$tmp_sero_file"
download_salmonella_serotype "$tmp_sero_file"
rm -f "$tmp_sero_file"
elif [[ "$taxon" =~ ^[A-Z][a-z]+\ [a-z]+$ ]]; then
echo "Detected species: $taxon"
download_species "$taxon" "$GENOME_DIR"
elif [[ "$taxon" =~ ^[A-Z][a-z]+$ ]]; then
echo "Detected genus: $taxon"
GENUS_DIR="$GENOME_DIR/$taxon"
mkdir -p "$GENUS_DIR"
download_single_genus "$taxon" "$GENOME_DIR"
if [[ "$GET_ALL_SPECIES" == true && "$taxon" =~ ^[A-Z][a-z]+$ ]]; then
echo "Expanding genus $taxon to species level"
SPECIES_LIST_FILE="$GENUS_DIR/species_list.txt"
get_species_list "$taxon" "$SPECIES_LIST_FILE"
classify_genus_by_metadata "$taxon"
find "$GENUS_DIR" -mindepth 1 -maxdepth 1 -type d ! -name "unclassified" | while read -r species_dir; do
if compgen -G "$species_dir/*_genomic.fna" > /dev/null; then
species_name=$(basename "$species_dir" | tr '_' ' ')
echo "$species_name" >> "$GENOME_DIR/expanded_species_list.txt"
fi
done
fi
else
echo "Unrecognized or unsupported taxon format: $taxon" >> "$FAILED_FLAG"
fi
) 2>/dev/null &
while (( $(jobs -r | wc -l) >= MAX_PARALLEL_JOBS )); do sleep 1; done
done
wait
echo "Building BLAST databases for custom panel"
build_custom_blastdb "$GENOME_DIR" "$BLAST_DB_DIR"
generate_directory_csv "$GENOME_DIR"
if [[ -s "$FAILED_FLAG" ]]; then
echo "Custom panel completed with some failures. See $FAILED_FLAG"
else
echo "Custom panel build complete."
touch "$CUSTOM_PANEL_CHECKPOINT"
fi
fi
else
GENOME_DIR=""
echo "Default mode: using prebuilt BLAST database."
fi
# if a single value is provided for -t instead of a file, create a temporary file
if [[ -n "$TAXON_FILE" && ! -f "$TAXON_FILE" ]]; then
tmp_taxon_file=$(mktemp)
echo "$TAXON_FILE" > "$tmp_taxon_file"
TAXON_FILE="$tmp_taxon_file"
trap "rm -f '$tmp_taxon_file'" EXIT
fi
# Set delimiter as a space
DELIMITER=" "
if [[ -n "$TAXON_FILE" ]]; then
sed 's/[\t]\+/ /g' "$TAXON_FILE" > "${TAXON_FILE}_processed"
# Ensure if a taxon file is provided, the genus should be one of the target Enterobacteriaceae
if [[ -z "$DOWNLOAD_FILE" ]]; then
awk -v delim="$DELIMITER" '
BEGIN { FS=delim; OFS=delim }
{ 
gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
genus = $1
species_or_serotype = ($2 != "") ? $2 : ""
if (genus !~ /^(Salmonella|Escherichia|Citrobacter|Enterobacter|Klebsiella|Shigella|Cronobacter)$/) {
print "Error: Invalid genus in taxon file. Allowed: Salmonella, Escherichia, Citrobacter, Enterobacter, Klebsiella, Shigella, Cronobacter in default. To proceed with other genus, use custom panel flag -d."
print "Your line:", $0
exit 1 
}
}' "${TAXON_FILE}_processed" || exit 1
fi
rm -f "${TAXON_FILE}_processed"
fi
if [[ "$MODE" == "light" && -z "$DOWNLOAD_FILE" && -z "$TAXON_FILE" ]]; then
echo "Default database mode (light): no target taxa list provided, scanning all available taxa in the prebuilt BLAST database."
elif [[ -n "$DOWNLOAD_FILE" ]]; then
echo "Custom database mode: skipping taxon genus restriction check (user-provided genome list)."
elif [[ -n "$TAXON_FILE" ]]; then
echo "Default database mode: restricted to taxa listed in the provided file."
fi
# start with core processing functions
process_complete_genomes() {
local taxon="$1"
echo "Processing complete genomes for $taxon"
local clean_taxon="$(extract_taxon_info "$taxon")"
local blast_output_name="${clean_taxon// /_}"
blast_output_name="${blast_output_name//./}"
local blast_output="${BLAST_RESULT_DIR}/${blast_output_name}_complete_blast_results.txt"
if [[ ! -s "$blast_output" ]]; then
perform_blast "$GENE_FILE" "$MIN_IDENTITY" "$BLAST_RESULT_DIR" "" "$taxon" "$BLAST_THREADS"
echo "BLAST result saved to $blast_output"
else
echo "Skipping BLAST: results already exist for $taxon"
fi
for blast_result_file in "$blast_output"; do
filter_blast_results "$blast_result_file" "$FILTERED_BLAST_RESULT_DIR" "$MIN_COVERAGE" "complete"
done
echo "Finished processing complete genomes for $taxon"
}

process_draft_genomes() {
local taxon="$1"
local clean_taxon="$(extract_taxon_info "$taxon")"
local blast_db_name="${clean_taxon// /_}"
blast_db_name="${blast_db_name//./}"
local total_draft_genomes=$(get_total_genomes_count "$taxon" "contig")
read -r sample_size iterations <<< "$(calculate_sample_size_and_iterations "$total_draft_genomes")"
echo "Processing $taxon | Total draft genomes: $total_draft_genomes. Running $iterations iterations (max 20)."
local max_parallel_iter=10
(( iterations < max_parallel_iter )) && max_parallel_iter=$iterations
local i
for ((i=1; i<=iterations; i++)); do
(
local iter="$i"
echo "Starting iterations $iter/$iterations for $taxon"
download_random_draft_genomes "$taxon" "$sample_size" "$DRAFT_GENOMES_DIR" "$iter"
perform_blast "$GENE_FILE" "$MIN_IDENTITY" "$DRAFT_BLAST_RESULT_DIR" "$iter" "$taxon" "$BLAST_THREADS"
mkdir -p "${DRAFT_BLAST_RESULT_DIR}"
blast_result_file="${DRAFT_BLAST_RESULT_DIR}/${blast_db_name}/iteration_${iter}_draft_blast_results.txt"
filter_blast_results "$blast_result_file" "$FILTERED_DRAFT_BLAST_RESULT_DIR" "$MIN_COVERAGE" "draft"
) 2>/dev/null &
while (( $(jobs -r | wc -l) >= max_parallel_iter )); do sleep 1; done
done
wait
}

# Set up the output file headers
if [[ "$MODE" == "heavy" ]]; then
echo -e "Organism,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,Draft_genomes_sample_size_per_iteration,Number_of_iterations,Complete_genomes_with_target_genes,Draft_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes,Percentage_with_target_genes_draft_genomes" > "$OUTPUT_FILE"
else
echo -e "Organism,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,Complete_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes" > "$OUTPUT_FILE"
fi
# Initialize TAXON_LIST
if [[ -n "$TAXON_FILE" ]]; then
TAXON_LIST=$(< "$TAXON_FILE")
elif [[ "$GET_ALL_SPECIES" == true && -f "$CUSTOM_GENOMES_DIR/expanded_species_list.txt" ]]; then
TAXON_LIST=$(< "$CUSTOM_GENOMES_DIR/expanded_species_list.txt")
elif [[ -n "$DOWNLOAD_FILE" ]]; then
TAXON_LIST=$(< "$DOWNLOAD_FILE")
else
# Rmove numeric suffix and convert underscores and dots to spaces
TAXON_LIST=$(find "$BLAST_DB_DIR" -name "*.nsq" -exec basename {} .nsq \; | \
sed -E '/[._][0-9]{2}$/s/[._][0-9]{2}$//' | sort -u |  sed 's/_/ /g' | sed -E 's/subsp /subsp. /g' | sed -E 's/ ([0-9]+)$/,\1/')
fi
# Process each taxon in TAXON_LIST using parallel summarization
MAX_PARALLEL_SUMMARY=$(( $(get_cpus) * 2 / 3 ))
(( MAX_PARALLEL_SUMMARY < 1 )) && MAX_PARALLEL_SUMMARY = 1
TMP_SUMMARY_DIR=$(mktemp -d "${WORKDIR}/summary_tmp.XXXXXX")
echo "$TAXON_LIST" | while IFS= read -r taxon; do
[[ -z "$taxon" ]] && continue
(
echo "Processing $taxon in $MODE mode"
process_complete_genomes "$taxon"
[[ "$MODE" == "heavy" ]] && process_draft_genomes "$taxon"
TOTAL_COMPLETE_GENOMES=$(get_total_genomes_count "$taxon" "complete") 
TOTAL_DRAFT_GENOMES=$(get_total_genomes_count "$taxon" "contig")
clean_taxon="$(extract_taxon_info "$taxon")"
local_taxon="${clean_taxon// /_}"
local_taxon="${local_taxon//./}"
if [[ "$MODE" == "heavy" && "$TOTAL_DRAFT_GENOMES" -gt 0 ]]; then
read -r DRAFT_SAMPLE_SIZE ITERATIONS <<< "$(calculate_sample_size_and_iterations "$TOTAL_DRAFT_GENOMES")"
else
DRAFT_SAMPLE_SIZE=0
ITERATIONS=0
fi
GENE_WITH_HITS=$(awk '{print $1}' "$FILTERED_BLAST_RESULT_DIR/filtered_${local_taxon}_complete_blast_results.txt" 2>/dev/null)
[[ "$MODE" == "heavy" ]] && GENE_WITH_HITS+=$'\n'$(awk '{print $1}' "$FILTERED_DRAFT_BLAST_RESULT_DIR/$local_taxon"/* 2>/dev/null)
GENE_WITH_HITS=$(echo "$GENE_WITH_HITS" | sort -u)
# Process all genes in query gene file
mapfile -t ALL_GENES < <(grep "^>" "$GENE_FILE" | sed 's/>//' | awk '{print $1}')
TMP_OUT="${TMP_SUMMARY_DIR}/${local_taxon}.csv"
for GENE_ID in "${ALL_GENES[@]}"; do
if grep -qx -- "$GENE_ID" <<< "$GENE_WITH_HITS" 2>/dev/null || true; then
COMPLETE_GENOMES_WITH_TARGET_GENES=$(awk -v gene="$GENE_ID" '$1 == gene {split($2, a, "|"); print a[1]}' "$FILTERED_BLAST_RESULT_DIR/filtered_${local_taxon}_complete_blast_results.txt" 2>/dev/null | sort -u | wc -l)
if [[ "$MODE" == "heavy" && "$TOTAL_DRAFT_GENOMES" -gt 0 ]]; then
TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES=0
for ((i=1; i<="$ITERATIONS"; i++)); do
DRAFT_GENOMES_WITH_TARGET_GENES=$(awk -v gene="$GENE_ID" '$1 == gene {print $2}' "$FILTERED_DRAFT_BLAST_RESULT_DIR/${local_taxon}/filtered_iteration_${i}_draft_blast_results.txt" 2>/dev/null | cut -c1-8 | sort -u | wc -l)
TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES=$((TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES + DRAFT_GENOMES_WITH_TARGET_GENES))
done
AVERAGE_DRAFT_GENOMES_WITH_TARGET_GENES=$(echo "$TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES/$ITERATIONS" | bc 2>/dev/null)
PERCENT_WITH_TARGET_GENES_DRAFT_GENOMES=$(echo "scale=2; ($AVERAGE_DRAFT_GENOMES_WITH_TARGET_GENES * 100)/$DRAFT_SAMPLE_SIZE" | bc 2>/dev/null)
else
DRAFT_GENOMES_WITH_TARGET_GENES=0
PERCENT_WITH_TARGET_GENES_DRAFT_GENOMES=0
fi
if [[ "$TOTAL_COMPLETE_GENOMES" -gt 0 ]]; then
PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES=$(echo "scale=2; ($COMPLETE_GENOMES_WITH_TARGET_GENES * 100)/$TOTAL_COMPLETE_GENOMES" | bc 2>/dev/null)
else
PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES=0
fi
if [[ "$MODE" == "heavy" ]]; then
echo -e "\"$taxon\",$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$DRAFT_SAMPLE_SIZE,$ITERATIONS,$COMPLETE_GENOMES_WITH_TARGET_GENES,$AVERAGE_DRAFT_GENOMES_WITH_TARGET_GENES,${PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES}%,${PERCENT_WITH_TARGET_GENES_DRAFT_GENOMES}%" >> "$TMP_OUT"
else
echo -e "\"$taxon\",$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$COMPLETE_GENOMES_WITH_TARGET_GENES,${PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES}%"  >> "$TMP_OUT"
fi
else
# If no hits were found for the gene in a serotype, output is recorded as 0%
if [[ "$MODE" == "heavy" ]]; then
echo -e "\"$taxon\",$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$DRAFT_SAMPLE_SIZE,$ITERATIONS,0,0,0%,0%" >> "$TMP_OUT"
else
echo -e "\"$taxon\",$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,0,0%" >> "$TMP_OUT"
fi
fi
done
) &
while (( $(jobs -r | wc -l) >= MAX_PARALLEL_SUMMARY )); do sleep 1; done
done
wait
cat "$TMP_SUMMARY_DIR"/*.csv >> "$OUTPUT_FILE"
rm -rf "$TMP_SUMMARY_DIR"
echo "Analysis complete. Results saved in $OUTPUT_FILE"
