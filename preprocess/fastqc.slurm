#!/usr/bin/env bash
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition="all"
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=2G

set -euo pipefail

# =============================================================================
# FastQC + MultiQC Quality Control Pipeline
# =============================================================================
source ../../dev/sequencing/preprocess/jobinfo.sh

# --- Defaults ----------------------------------------------------------------
RUN_NAME="${RUN_NAME:-test}"
EXT="${EXT:-.fastq.gz}"
output_dir="./metrics/fastqc"
run_name="test"
fastqc="${FASTQC:-fastqc}"      # assumes fastqc is on $PATH
multiqc="${MULTIQC:-multiqc}"   # assumes multiqc is on $PATH


# --- Logs: use job name if set, fallback to RUN_NAME -------------------------
# LOG_DIR="$(pwd)/logs"
# mkdir -p "${LOG_DIR}"

# LOG_PREFIX="${SLURM_JOB_NAME:-${RUN_NAME}}_${SLURM_JOB_ID:-local}"

# exec > >(tee -a "${LOG_DIR}/${LOG_PREFIX}.out") \
#      2> >(tee -a "${LOG_DIR}/${LOG_PREFIX}.err" >&2)


# --- Usage -------------------------------------------------------------------
usage() {
    echo "
Usage: bash $(basename "$0") -r <path> [-fqc <path>] [-mqc <path>] [-o <path>] [-n <str>] [-e <str>]

Required:
    -r|--raw-dir      : directory containing input read files

Optional:
    -fqc|--fastqc     : path to FastQC executable  [default: fastqc]
    -mqc|--multiqc    : path to MultiQC executable [default: multiqc]
    -o|--output-dir   : base output directory       [default: ./metrics/fastqc]
    -n|--run-name     : name of the run             [default: test]
    -e|--ext          : input file extension        [default: .fastq.gz]
                        e.g: .fq.gz .fastq .fq
    -h|--help         : show this message and exit
"
    exit 1
}

if [[ $# -lt 1 ]]; then
    usage
fi

# --- Argument parsing --------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--raw-dir)    raw_dir="${2%/}";    shift 2 ;;
        -fqc|--fastqc)   fastqc="$2";         shift 2 ;;
        -mqc|--multiqc)  multiqc="$2";        shift 2 ;;
        -o|--output-dir) output_dir="${2%/}"; shift 2 ;;
        -n|--run-name)   run_name="$2";       shift 2 ;;
        -e|--ext)        EXT="$2";            shift 2 ;;
        -h|--help)       usage ;;
        *)
            echo "ERROR: Unknown argument: $1"
            usage
            ;;
    esac
done

# --- Validate required args --------------------------------------------------
missing=()
[[ -z "${raw_dir:-}"  ]] && missing+=("-r/--raw-dir")
[[ -z "${fastqc:-}"   ]] && missing+=("-fqc/--fastqc")
[[ -z "${multiqc:-}"  ]] && missing+=("-mqc/--multiqc")

if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: Missing required argument(s): ${missing[*]}"
    usage
fi

run_dir="${output_dir}/${run_name}"
fastqc_out="${run_dir}/fastqc"
multiqc_out="${run_dir}/multiqc"

# --- Logging helper ----------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# --- Input validation --------------------------------------------------------
log "Validating inputs..."

input_files=( "${raw_dir}"/*"${EXT}" )

if [[ ! -e "${input_files[0]}" ]]; then
    log "ERROR: No *${EXT} files found in ${raw_dir}"
    exit 1
fi

log "Run name   : ${run_name}"
log "Input dir  : ${raw_dir}"
log "Extension  : ${EXT}"
log "Run dir    : ${run_dir}"
log "FastQC out : ${fastqc_out}"
log "MultiQC out: ${multiqc_out}"
log "Samples (${#input_files[@]}):"
for f in "${input_files[@]}"; do
    log "  $(basename "${f}")"
done

mkdir -p "${fastqc_out}" "${multiqc_out}"

# --- FastQC ------------------------------------------------------------------
log "Running FastQC..."

fastqc_cmd=(
    "${fastqc}"
    --threads "${SLURM_CPUS_PER_TASK:-4}"
    --outdir  "${fastqc_out}"
    "${input_files[@]}"
)

# set JVM heap size if SLURM memory is available
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    export _JAVA_OPTIONS="-Xmx${SLURM_MEM_PER_NODE}m"
    log "JVM memory : -Xmx${SLURM_MEM_PER_NODE}m"
fi

log "FastQC command: ${fastqc_cmd[*]}"
"${fastqc_cmd[@]}"

log "FastQC complete."

# --- MultiQC -----------------------------------------------------------------
log "Running MultiQC..."

"${multiqc}" \
    "${fastqc_out}" \
    --outdir   "${multiqc_out}" \
    --filename "${run_name}" \
    --verbose

log "MultiQC complete. Report: ${multiqc_out}/${run_name}.html"

# --- Job info ----------------------------------------------------------------
log "Pipeline finished."

jobinfo