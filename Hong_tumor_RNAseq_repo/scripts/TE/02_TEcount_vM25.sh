#!/bin/bash
#SBATCH --job-name=tumor_TEcount_vM25
#SBATCH --array=0-11%6
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --partition=HTC
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=ZIW183@pitt.edu
#SBATCH --mail-type=END,FAIL
# Hong tumor RNA-seq: TEcount (TEtranscripts, subfamily level) on the vM25 STAR BAMs
# Submit after 01 finishes:  sbatch --dependency=afterok:<01_JOBID> 03_tumor_TEcount_vM25.sh
# (can run at the same time as 02)

source ~/.bashrc
source /ix1/hzhang/SteveWang/TEtranscripts_venv/bin/activate
conda activate /vast/bnacev/ziw183/conda_envs/nacevlab
set -eo pipefail

# --- Configuration ---
BASE=/ix1/hzhang/SteveWang/Hong_data/tumor_RNA_Seq
STAR_OUT=$BASE/STAR_vM25
TECOUNT_OUT=$BASE/TE_res_vM25/TEcount
GTF_PATH=/ix1/bnacev/ziw183/ref/grcm38/gencode.vM25.annotation.gtf
TE_GTF="/ix1/bnacev/ziw183/ref/mm10/mm10_rmsk_TE.gtf"   # Hammell-format TE GTF (same source as mm10_rmsk_TE.gtf.locInd)
STRANDED="no"       # Salmon detected IU = unstranded

SAMPLES=(NB_2R NB_5L NB_6L  ND_9L ND_10R_2 ND_11L_2  S7B_24R S7B_25R S7B_28L  S7D_30L S7D_30R S7D_34L)
# --- End of Configuration ---

S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
BAM=$STAR_OUT/$S/${S}_Aligned.sortedByCoord.out.bam
[ -f "$TE_GTF" ] || { echo "ERROR: TE GTF not found: $TE_GTF"; exit 1; }
samtools quickcheck "$BAM" || { echo "ERROR: BAM missing/truncated: $BAM"; exit 1; }
mkdir -p "$TECOUNT_OUT"

echo "TEcount: $S  (--stranded $STRANDED)"
TEcount --sortByPos \
        --format BAM \
        --mode multi \
        --stranded "$STRANDED" \
        -b "$BAM" \
        --GTF "$GTF_PATH" \
        --TE "$TE_GTF" \
        --project "$TECOUNT_OUT/$S"

{
  echo "sample: $S   date: $(date)"
  echo "$(TEcount --version 2>&1 | head -1)"
  echo "TEcount: --sortByPos --format BAM --mode multi --stranded $STRANDED | GTF=$GTF_PATH | TE=$TE_GTF"
} > "$TECOUNT_OUT/${S}_run_info.txt"
echo "Finished $S -> $TECOUNT_OUT/${S}.cntTable"
