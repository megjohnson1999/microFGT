#!/bin/bash
# ============================================================================
# FRESH Stage 2 — Preprocessing (fastp trim + minimap2/samtools host removal)
# ARRAY JOB: one sample per task.  Validation: sbatch --array=1 ...
#            full run:            sbatch --array=1-11 ...
# IN : reads/<id>_1/_2.fastq.gz (symlinks -> /lts FRESH fastqs)
# OUT: 02_preproc/<id>/<id>_nonhost_R1/R2.fastq.gz
# Deep samples (12-175M pairs) -> array parallelism + big time/mem vs the ENA run.
# ============================================================================
#SBATCH --job-name=fresh_preproc
#SBATCH --array=1-11
#SBATCH --mem=72G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/sahlab/Megan/metagenomics_fresh/logs/preproc_%A_%a.out
set -euo pipefail

WD=/scratch/sahlab/Megan/metagenomics_fresh
READS=$WD/reads
IDS=$WD/metadata/sample_ids.txt
OUT=$WD/02_preproc
HOST_REF=/ref/sahlab/data/GRCh38.fna.gz
THREADS=${SLURM_CPUS_PER_TASK:-16}

acc=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$IDS")
[ -z "$acc" ] && { echo "no sample at array index ${SLURM_ARRAY_TASK_ID}"; exit 1; }

set +u
source /ref/sahlab/software/miniforge3/bin/activate
conda activate fgt_preprocess
set -u

s=$OUT/$acc; mkdir -p "$s"
r1=$READS/${acc}_1.fastq.gz
r2=$READS/${acc}_2.fastq.gz
t1=$s/${acc}_trim_R1.fastq.gz
t2=$s/${acc}_trim_R2.fastq.gz
n1=$s/${acc}_nonhost_R1.fastq.gz
n2=$s/${acc}_nonhost_R2.fastq.gz

if [ -s "$n1" ] && [ -s "$n2" ]; then
  echo "[$(date "+%F %T")] $acc already done, skipping"; exit 0
fi
echo "[$(date "+%F %T")] preprocessing $acc  (host_ref=$HOST_REF, threads=$THREADS)"

# --- 2a. fastp: adapter (auto PE) + quality trim ---
fastp --in1 "$r1" --in2 "$r2" --out1 "$t1" --out2 "$t2" \
  --thread "$THREADS" \
  --json "$s/${acc}_fastp.json" --html "$s/${acc}_fastp.html" \
  2> "$s/${acc}_fastp.log"

# --- 2b. host removal: keep pairs where BOTH mates are unmapped to GRCh38 ---
minimap2 -ax sr -t "$THREADS" "$HOST_REF" "$t1" "$t2" \
  | samtools view -@ "$THREADS" -b -f 12 -F 256 - \
  | samtools sort -n -@ "$THREADS" -m 2G -T "$s/${acc}_sorttmp" - \
  | samtools fastq -@ "$THREADS" -1 "$n1" -2 "$n2" -0 /dev/null -s /dev/null - \
  2> "$s/${acc}_hostremoval.log"

# --- tallies: raw vs clean, and host-removed fraction ---
raw=$(( $(zcat "$r1" | wc -l) / 4 ))
clean=$(( $(zcat "$n1" | wc -l) / 4 ))
pct=$(awk -v r="$raw" -v c="$clean" 'BEGIN{ if(r>0) printf "%.1f", 100*(r-c)/r; else print "NA" }')
echo "  raw_pairs=$raw  clean_pairs=$clean  host_removed=${pct}%" | tee -a "$s/${acc}_hostremoval.log"

# tidy the big trimmed intermediates once host removal succeeded
rm -f "$t1" "$t2"
echo "[$(date "+%F %T")] $acc done"
