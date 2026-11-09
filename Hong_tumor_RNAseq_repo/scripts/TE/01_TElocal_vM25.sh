#!/bin/bash
#SBATCH --job-name=tumor_TElocal_vM25
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
# Hong tumor RNA-seq: TElocal (locus level) on the vM25 STAR BAMs from script 01
# Submit after 01 finishes:  sbatch --dependency=afterok:<01_JOBID> 02_tumor_TElocal_vM25.sh

source ~/.bashrc
source /ix1/hzhang/SteveWang/TEtranscripts_venv/bin/activate
conda activate /vast/bnacev/ziw183/conda_envs/nacevlab
set -eo pipefail

# --- Configuration ---
BASE=/ix1/hzhang/SteveWang/Hong_data/tumor_RNA_Seq
STAR_OUT=$BASE/STAR_vM25
TELOCAL_OUT=$BASE/TE_res_vM25/TElocal
GTF_PATH=/ix1/bnacev/ziw183/ref/grcm38/gencode.vM25.annotation.gtf
TE_LOCIND=/ix1/bnacev/ziw183/ref/mm10/mm10_rmsk_TE.gtf.locInd   # same source as TEcount TE GTF (ref/mm10/mm10_rmsk_TE.gtf)
STRANDED="no"       # Salmon detected IU = unstranded

SAMPLES=(NB_2R NB_5L NB_6L  ND_9L ND_10R_2 ND_11L_2  S7B_24R S7B_25R S7B_28L  S7D_30L S7D_30R S7D_34L)
# --- End of Configuration ---

S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
BAM=$STAR_OUT/$S/${S}_Aligned.sortedByCoord.out.bam
samtools quickcheck "$BAM" || { echo "ERROR: BAM missing/truncated: $BAM"; exit 1; }
mkdir -p "$TELOCAL_OUT"

echo "TElocal: $S  (--stranded $STRANDED)"
TElocal --sortByPos \
        -b "$BAM" \
        --GTF "$GTF_PATH" \
        --TE "$TE_LOCIND" \
        --stranded "$STRANDED" \
        --project "$TELOCAL_OUT/$S"

{
  echo "sample: $S   date: $(date)"
  echo "$(TElocal --version 2>&1 | head -1)"
  echo "TElocal: --sortByPos --stranded $STRANDED (mode: multi) | GTF=$GTF_PATH | TE=$TE_LOCIND"
  echo "BAM: $BAM (STAR vM25, --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200)"
} > "$TELOCAL_OUT/${S}_run_info.txt"
echo "Finished $S -> $TELOCAL_OUT/${S}.cntTable"
