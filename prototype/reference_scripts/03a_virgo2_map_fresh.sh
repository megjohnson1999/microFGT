#!/bin/bash
# ============================================================================
# FRESH Stage 3a — VIRGO2 mapping (ARRAY: one sample per task)
# IN : 02_preproc/<id>/<id>_nonhost_R1/R2.fastq.gz
# OUT: 03_cervicotype/virgo2/<id>.out  (gene counts)
# Reuses the VIRGO2 install + bowtie2 index built in the audit run.
# Submit: sbatch --array=1-11 ...   (validation: --array=1)
# ============================================================================
#SBATCH --job-name=fresh_virgo2
#SBATCH --array=1-11
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/sahlab/Megan/metagenomics_fresh/logs/virgo2_%A_%a.out
set -euo pipefail

WD=/scratch/sahlab/Megan/metagenomics_fresh
IDS=$WD/metadata/sample_ids.txt
PREP=$WD/02_preproc
OUT=$WD/03_cervicotype/virgo2
READS_TMP=$WD/03_cervicotype/tmp_combined
VIRGO2_DIR=/scratch/sahlab/Megan/metagenomics_fgt/resources/VIRGO2   # reuse built index (Index/VIRGO2.*.bt2)
mkdir -p "$OUT" "$READS_TMP"
THREADS=${SLURM_CPUS_PER_TASK:-16}

acc=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IDS")
[ -z "$acc" ] && { echo "no sample at array index ${SLURM_ARRAY_TASK_ID}"; exit 1; }

set +u
source /ref/sahlab/software/miniforge3/bin/activate
conda activate fgt_preprocessing_cervicotyping   # bowtie2 + samtools + pandas/numpy
set -u

r1=$PREP/$acc/${acc}_nonhost_R1.fastq.gz
r2=$PREP/$acc/${acc}_nonhost_R2.fastq.gz
combined=$READS_TMP/${acc}_combined.fq.gz

if [ -s "$OUT/${acc}.out" ]; then echo "[$(date "+%F %T")] $acc .out exists, skipping"; exit 0; fi
if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then echo "ERROR: missing nonhost reads for $acc (run Stage 2 first)"; exit 1; fi

echo "[$(date "+%F %T")] VIRGO2 map $acc  (single-end: R1+R2 concatenated)"
cat "$r1" "$r2" > "$combined"
python3 "$VIRGO2_DIR/VIRGO2.py" map -r "$combined" -o "$OUT/$acc" -p "$THREADS"
rm -f "$combined"
echo "[$(date "+%F %T")] $acc done -> $OUT/${acc}.out"
