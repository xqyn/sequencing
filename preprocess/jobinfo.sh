#!/bin/bash

# jobinfo
#
# Description:
#   Prints a structured summary of SLURM job-related information, including job metadata,
#   resource allocation (CPU, memory, GPU), and wall-clock runtime with millisecond precision.
#
# Globals:
#   Uses SLURM environment variables, such as:
#     SLURM_JOB_NAME, SLURM_JOBID, SLURM_JOB_USER,
#     SLURM_SUBMIT_HOST, SLURMD_NODENAME, SLURM_SUBMIT_DIR,
#     SLURM_NNODES, SLURM_JOB_NUM_NODES, SLURM_NTASKS,
#     SLURM_NTASKS_PER_NODE, SLURM_TASKS_PER_NODE, SLURM_NPROCS,
#     SLURM_CPUS_ON_NODE, SLURM_CPUS_PER_TASK, SLURM_JOB_CPUS_PER_NODE,
#     SLURM_JOB_GPUS, SLURM_MEM_PER_NODE, SLURM_MEM_PER_CPU
#
#   Also relies on global variables set externally:
#     date_start - Start time in human-readable format (e.g., from `date`)
#     start      - Start time in seconds.nanoseconds format (`date +%s.%N`)
#
# Output:
#   Prints formatted job information and runtime to standard output.
#   Includes millisecond-level timing for precise runtime diagnostics.
#
# Usage:
#   Call this function after your SLURM job completes or at desired checkpoints.
#
# Example:
#   date_start=$(date)
#   start=$(date +%s.%N)
#   ...
#   jobinfo


date_start="${date_start:-$(date)}"
start="${start:-$(date +%s.%N)}"


jobinfo() {
    local date_end=$(date)
    local runtime=$(echo "$(date +%s.%N) - $start" | bc -l)
    local runtime_min=$(echo "scale=2; $runtime / 60" | bc)

    printf "%s\n" 
    printf "%s\n" 
    echo "# -------------------- JOB INFORMATION --------------------"

    printf "%-22s: %s\n" "Job Name"             "${SLURM_JOB_NAME:-N/A}"
    printf "%-22s: %s\n" "Job ID"               "${SLURM_JOBID:-N/A}"
    printf "%-22s: %s\n" "Job User"             "${SLURM_JOB_USER:-N/A}"

    printf "\n%-22s: %s\n" "Submit Host"        "${SLURM_SUBMIT_HOST:-N/A}"
    printf "%-22s: %s\n" "Run Node"             "${SLURMD_NODENAME:-N/A}"
    printf "%-22s: %s\n" "Submit Dir"           "${SLURM_SUBMIT_DIR:-N/A}"

    printf "\n%-22s: %s\n" "Nodes"              "${SLURM_NNODES:-N/A}"
    printf "%-22s: %s\n" "Jobs per Node"        "${SLURM_JOB_NUM_NODES:-N/A}"
    printf "%-22s: %s\n" "Tasks"                "${SLURM_NTASKS:-N/A}"
    printf "%-22s: %s\n" "Tasks per Node"       "${SLURM_NTASKS_PER_NODE:-N/A}"
    printf "%-22s: %s\n" "Tasks per Node (alt)" "${SLURM_TASKS_PER_NODE:-N/A}"
    printf "%-22s: %s\n" "Processors"           "${SLURM_NPROCS:-N/A}"

    # printf "Nodes: %s | Jobs per Node: %s | Tasks: %s | Tasks per Node: %s | Tasks per Node (alt): %s | Processors: %s\n" \
    #     "$SLURM_NNODES" "$SLURM_JOB_NUM_NODES" "$SLURM_NTASKS" "$SLURM_NTASKS_PER_NODE" "$SLURM_TASKS_PER_NODE" "$SLURM_NPROCS"

    printf "\n%-22s: %s\n" "CPUs per Node"      "${SLURM_CPUS_ON_NODE:-N/A}"
    printf "%-22s: %s\n" "CPUs per Task"        "${SLURM_CPUS_PER_TASK:-N/A}"
    printf "%-22s: %s\n" "CPUs per Node (Job)"  "${SLURM_JOB_CPUS_PER_NODE:-N/A}"
    
    printf "%s\n" 
    if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
        mem_node_gb=$((SLURM_MEM_PER_NODE / 1024))
        printf "%-22s: %s MB (%s GB)\n" "Memory per Node" "$SLURM_MEM_PER_NODE" "$mem_node_gb"
    else
        printf "%-22s: %s\n" "Memory per Node" "N/A"
    fi

    printf "%-22s: %s\n" "Memory per CPU" "${SLURM_MEM_PER_CPU:-N/A}"

    # For GPUs per Node with fallback to N/A
    if [[ -n "${SLURM_JOB_GPUS:-}" ]]; then
        printf "%-22s: %s\n" "GPUs per Node" "$SLURM_JOB_GPUS"
    else
        printf "%-22s: %s\n" "GPUs per Node" "N/A"
    fi


    echo -e "\n# -------------------- RUNTIME INFORMATION ----------------"

    printf "%-22s: %s\n" "Date Start"           "$date_start"
    printf "%-22s: %s\n" "Date End"             "$date_end"
    printf "%-22s: %.2f minutes (%.3f seconds)\n" "Run Time" "$runtime_min" "$runtime"

    # Add SLURM-reported RunTime and TimeLimit if JobID is available
    if [[ -n "${SLURM_JOBID:-}" ]]; then
        local RUNTIME=$(scontrol show job "$SLURM_JOBID" | grep -oP 'RunTime=\K[0-9:]+')
        local TIMELIMIT=$(scontrol show job "$SLURM_JOBID" | grep -oP 'TimeLimit=\K[0-9:]+')
        printf "%-22s: %s\n" "SLURM Run Time"    "$RUNTIME"
        printf "%-22s: %s\n" "SLURM Time Limit"  "$TIMELIMIT"
    fi

    echo -e "\n# ==================== JOB END ============================"
    printf "%s\n" 
    printf "%s\n" 
}