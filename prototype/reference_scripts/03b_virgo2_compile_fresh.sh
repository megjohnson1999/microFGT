#!/bin/bash
# ============================================================================
# FRESH Stage 3b — VIRGO2 compile: all per-sample .out -> VISTA input matrix
# Run AFTER the 03a array finishes:
#   sbatch --dependency=afterok:<03a_arrayjobid> 03b_virgo2_compile_fresh.sh
# OUT: 03_cervicotype/virgo2/VIRGO2_Compiled.summary.NR.txt
# ============================================================================
#SBATCH --job-name=fresh_v2compile
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/sahlab/Megan/metagenomics_fresh/logs/virgo2_compile_%j.out
set -euo pipefail

WD=/scratch/sahlab/Megan/metagenomics_fresh
OUT=$WD/03_cervicotype/virgo2
VIRGO2_DIR=/scratch/sahlab/Megan/metagenomics_fgt/resources/VIRGO2

set +u
source /ref/sahlab/software/miniforge3/bin/activate
conda activate fgt_preprocessing_cervicotyping
set -u

n=$(ls "$OUT"/*.out 2>/dev/null | wc -l)
echo "[$(date "+%F %T")] VIRGO2 compile over ${n} per-sample .out files"
python3 "$VIRGO2_DIR/VIRGO2.py" compile -i "$OUT" -o "$OUT/VIRGO2_Compiled"
echo "[$(date "+%F %T")] done -> $OUT/VIRGO2_Compiled.summary.NR.txt"
