#!/bin/bash
# ============================================================================
# FRESH Stage 3c — VISTA mgCST classification
# Reuses the vista conda env + VISTA repo/data built in the audit run.
# NOTE: run_VISTA.R arg2 must be the dir CONTAINING VISTA_data/ (the repo dir),
#       NOT VISTA_data itself; and it writes outputs to CWD -> we cd into OUTDIR.
# Run AFTER 03b compile:  sbatch --dependency=afterok:<03b_jobid> 03c_vista_fresh.sh
# ============================================================================
#SBATCH --job-name=fresh_vista
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/sahlab/Megan/metagenomics_fresh/logs/vista_%j.out
set -euo pipefail

FRESH=/scratch/sahlab/Megan/metagenomics_fresh
AUDIT=/scratch/sahlab/Megan/metagenomics_fgt
VISTA_REPO=$AUDIT/resources/VISTA        # contains VISTA_data/  (this is run_VISTA.R arg2)
VISTA_ENV=$AUDIT/envs/vista
COMPILED=$FRESH/03_cervicotype/virgo2/VIRGO2_Compiled.summary.NR.txt
OUTDIR=$FRESH/03_cervicotype/vista
mkdir -p "$OUTDIR"

if [ ! -s "$COMPILED" ]; then echo "ERROR: compiled matrix not found: $COMPILED (run Stage 3b first)"; exit 1; fi

set +u
source /ref/sahlab/software/miniforge3/bin/activate
conda activate "$VISTA_ENV"
set -u

cd "$OUTDIR"    # VISTA writes outputs to getwd()
echo "[$(date "+%F %T")] VISTA on $COMPILED"
Rscript "$VISTA_REPO/run_VISTA.R" "$COMPILED" "$VISTA_REPO"
echo "[$(date "+%F %T")] VISTA done -> $OUTDIR"
ls -lh "$OUTDIR"
