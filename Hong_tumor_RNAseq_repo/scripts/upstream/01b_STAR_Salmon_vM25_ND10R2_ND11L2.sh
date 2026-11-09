#!/bin/bash
#SBATCH --job-name=tumor_STAR_Salmon_vM25_ND2
#SBATCH --array=4,5
#SBATCH --time=16:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --partition=HTC
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-user=ZIW183@pitt.edu
#SBATCH --mail-type=END,FAIL
# Hong tumor RNA-seq (NT/sg7 x PBS/PD1): retrim from RAW -> validate -> STAR + Salmon (GENCODE vM25)
# The old *_trimmed_R*.fq.gz files are corrupt (interleaved records), so this script
# re-trims from raw and records whether the old trimmed files were corrupt.
# Submit:  sbatch 01b_STAR_Salmon_vM25_ND10R2_ND11L2.sh   (runs ONLY tasks 4,5 = ND_10R_2, ND_11L_2)
# Rerun only some samples, e.g.:  sbatch --array=4,5 01_STAR_Salmon_vM25.sh
# Reruns are incremental: validated trimmed fastq and valid STAR BAMs are reused (FORCE=1 redoes all).

source ~/.bashrc
source /ix1/hzhang/SteveWang/TEtranscripts_venv/bin/activate
conda activate /vast/bnacev/ziw183/conda_envs/nacevlab
ulimit -n 65535
set -eo pipefail

# --- Configuration ---
BASE=/ix1/hzhang/SteveWang/Hong_data/tumor_RNA_Seq
# Raw data: one folder per sample containing *_1.fq.gz / *_2.fq.gz (Novogene layout).
# Add more folders here if some samples (e.g. ND_10R_2, ND_11L_2) live elsewhere.
RAW_DIRS=(/ix1/hzhang/HongCao/RNAseq_202604_org_tumor/tumor/01.RawData
          /ix1/hzhang/HongCao/RNAseq_202604_org_tumor/org/01.RawData)   # ND_10R_2, ND_11L_2 live here
OLD_TRIM_DIR=$BASE/Salmon_Results_grcm38           # old (corrupt?) trimmed fastq, checked only
TRIM_OUT=$BASE/trimmed_vM25
STAR_OUT=$BASE/STAR_vM25
SALMON_OUT=$BASE/Salmon_Results_vM25
REF=/ix1/bnacev/ziw183/ref/grcm38/vM25
STAR_IDX=$REF/STAR_index_vM25
SALMON_IDX=$REF/salmon_index_vM25_decoy

ADAPTER="AGATCGGAAGAGC"     # Illumina TruSeq
NEXTSEQ_Q=20                # 2-colour chemistry (NovaSeq X): poly-G aware quality trim
MIN_LEN=20

SAMPLES=(NB_2R NB_5L NB_6L  ND_9L ND_10R_2 ND_11L_2  S7B_24R S7B_25R S7B_28L  S7D_30L S7D_30R S7D_34L)
THREADS=$SLURM_CPUS_PER_TASK
FORCE=0
# --- End of Configuration ---

S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "========== $S =========="
for t in cutadapt STAR salmon samtools; do command -v $t >/dev/null || { echo "ERROR: $t not in env"; exit 1; }; done
GZ=gzip; command -v pigz >/dev/null && GZ="pigz -p $THREADS"

# one-pass FASTQ check: 4-line records, '@' header, '+' separator, seq/qual equal length
check_fq() {   # $1 = file (.gz or plain); prints record count, returns 1 if malformed
    local cat_cmd=cat; [[ "$1" == *.gz ]] && cat_cmd=zcat
    $cat_cmd "$1" | awk '
        NR%4==1 { if (substr($0,1,1)!="@") {bad=1; exit} }
        NR%4==2 { s=length($0) }
        NR%4==3 { if (substr($0,1,1)!="+") {bad=1; exit} }
        NR%4==0 { if (length($0)!=s) {bad=1; exit} }
        END { if (bad || NR%4) { print "BAD_at_line_" NR; exit 1 } print NR/4 }'
}

# ---- 0. Was the OLD trimmed pair corrupt? (for the record) ----
OLD1=$OLD_TRIM_DIR/$S/${S}_trimmed_R1.fq.gz; OLD2=$OLD_TRIM_DIR/$S/${S}_trimmed_R2.fq.gz
old_status="not_found"
if [ -f "$OLD1" ] && [ -f "$OLD2" ]; then
    o1=$(check_fq "$OLD1" || true); o2=$(check_fq "$OLD2" || true)
    if [[ "$o1" == BAD* || "$o2" == BAD* ]]; then old_status="CORRUPT (R1:$o1 R2:$o2)"
    elif [ "$o1" != "$o2" ]; then old_status="CORRUPT (unequal reads R1=$o1 R2=$o2)"
    else old_status="ok ($o1 pairs)"; fi
fi
echo "Old trimmed fastq: $old_status"

T1=$TRIM_OUT/$S/${S}_trimmed_1.fq; T2=$TRIM_OUT/$S/${S}_trimmed_2.fq
MARK=$TRIM_OUT/$S/${S}.trim_validated
mkdir -p "$TRIM_OUT/$S"
# reuse a previously validated trim (plain, or compressed by a finished run)
if [ "$FORCE" -eq 0 ] && [ -f "$MARK" ]; then
    if [ ! -f "$T1" ] && [ -f "$T1.gz" ]; then T1=$T1.gz; T2=$T2.gz; fi
fi
# trimmed files left by an earlier (e.g. OOM-killed) run without a marker: re-validate once
if [ "$FORCE" -eq 0 ] && [ ! -f "$MARK" ] && [ -s "$T1" ] && [ -s "$T2" ]; then
    echo "Found unmarked trimmed fastq from an earlier run — validating..."
    if c1=$(check_fq "$T1") && c2=$(check_fq "$T2") && [ "$c1" = "$c2" ]; then
        echo "NA $c1" > "$MARK"; echo "  valid ($c1 pairs) — will reuse"
    else
        echo "  invalid — will re-trim"
    fi
fi
if [ "$FORCE" -eq 0 ] && [ -f "$MARK" ] && [ -f "$T1" ] && [ -f "$T2" ]; then
    read -r raw_pairs n1 < "$MARK"
    RAW1=("(reused trim)"); RAW2=("(reused trim)")
    echo "Reusing validated trimmed fastq: $raw_pairs raw / $n1 trimmed pairs"
else
# ---- 1. Locate raw fastq (merge lanes if several) ----
RAW1=(); RAW2=()
for d in "${RAW_DIRS[@]}"; do
    if [ -d "$d/$S" ]; then
        mapfile -t RAW1 < <(ls "$d/$S"/*_1.fq.gz 2>/dev/null | sort)
        mapfile -t RAW2 < <(ls "$d/$S"/*_2.fq.gz 2>/dev/null | sort)
        [ ${#RAW1[@]} -gt 0 ] && break
    fi
done
if [ ${#RAW1[@]} -eq 0 ] || [ ${#RAW1[@]} -ne ${#RAW2[@]} ]; then
    echo "ERROR: raw fastq for $S not found in: ${RAW_DIRS[*]}  (add its folder to RAW_DIRS)"; exit 1
fi
echo "Raw R1: ${RAW1[*]}"
if [ ${#RAW1[@]} -gt 1 ]; then          # gzip files can be concatenated directly
    cat "${RAW1[@]}" > "$TRIM_OUT/$S/${S}_raw_merged_1.fq.gz"
    cat "${RAW2[@]}" > "$TRIM_OUT/$S/${S}_raw_merged_2.fq.gz"
    IN1=$TRIM_OUT/$S/${S}_raw_merged_1.fq.gz; IN2=$TRIM_OUT/$S/${S}_raw_merged_2.fq.gz
else
    IN1=${RAW1[0]}; IN2=${RAW2[0]}
fi
raw_pairs=$(check_fq "$IN1") || { echo "ERROR: RAW R1 malformed: $raw_pairs"; exit 1; }

# ---- 2. Trim to UNCOMPRESSED output (avoids the gzip race), then validate ----
rm -f "$T1" "$T2" "$T1.gz" "$T2.gz" "$MARK"
cutadapt -a "$ADAPTER" -A "$ADAPTER" \
         --nextseq-trim=$NEXTSEQ_Q \
         -m ${MIN_LEN}:${MIN_LEN} \
         -j 8 \
         -o "$T1" -p "$T2" \
         "$IN1" "$IN2" > "$TRIM_OUT/$S/${S}_cutadapt.log"
n1=$(check_fq "$T1") || { echo "ERROR: trimmed R1 malformed ($n1)"; exit 1; }
n2=$(check_fq "$T2") || { echo "ERROR: trimmed R2 malformed ($n2)"; exit 1; }
[ "$n1" = "$n2" ] || { echo "ERROR: R1/R2 read counts differ ($n1 vs $n2)"; exit 1; }
echo "Raw pairs: $raw_pairs | trimmed pairs: $n1  (validated)"
[ ${#RAW1[@]} -gt 1 ] && rm -f "$IN1" "$IN2"
echo "$raw_pairs $n1" > "$MARK"
fi
RFC=(); [[ "$T1" == *.gz ]] && RFC=(--readFilesCommand zcat)

# ---- 3. STAR (TE-ready multimapping; intronMotif for unstranded data) ----
mkdir -p "$STAR_OUT/$S"
P=$STAR_OUT/$S/${S}_
BAM=${P}Aligned.sortedByCoord.out.bam
if [ "$FORCE" -eq 0 ] && samtools quickcheck "$BAM" 2>/dev/null && [ -f "${P}Log.final.out" ] && [ -f "${P}ReadsPerGene.out.tab" ]; then
    echo "Valid STAR BAM exists — skipping STAR."
else
rm -rf "$BAM" "${P}_STARtmp"
STAR --runThreadN "$THREADS" \
     --genomeDir "$STAR_IDX" \
     --readFilesIn "$T1" "$T2" "${RFC[@]}" \
     --outFileNamePrefix "$P" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMstrandField intronMotif \
     --outFilterMultimapNmax 100 \
     --winAnchorMultimapNmax 200 \
     --quantMode GeneCounts \
     --limitBAMsortRAM 20000000000
samtools quickcheck "$BAM"
fi
[ -f "$BAM.bai" ] || samtools index -@ "$THREADS" "$BAM"

# ---- 4. Salmon (selective alignment, decoy-aware) ----
rm -rf "$SALMON_OUT/$S"
salmon quant -i "$SALMON_IDX" -l A \
             -1 "$T1" -2 "$T2" \
             --gcBias --seqBias \
             -p "$THREADS" \
             -o "$SALMON_OUT/$S"

# ---- 5. compress trimmed fastq (single-threaded per file = safe) ----
if [[ "$T1" != *.gz ]]; then $GZ "$T1"; $GZ "$T2"; fi

{
  echo "sample: $S   date: $(date)"
  echo "raw: ${RAW1[*]} | ${RAW2[*]}"
  echo "old_trimmed_fastq_status: $old_status"
  echo "raw_pairs: $raw_pairs   trimmed_pairs: $n1"
  echo "cutadapt $(cutadapt --version): -a $ADAPTER -A $ADAPTER --nextseq-trim=$NEXTSEQ_Q -m ${MIN_LEN}:${MIN_LEN}"
  echo "STAR $(STAR --version) | salmon $(salmon --version 2>&1) | $(samtools --version | head -1)"
  echo "STAR: --outSAMstrandField intronMotif --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 --quantMode GeneCounts"
  echo "salmon: -l A --gcBias --seqBias (selective alignment, default) | index=$SALMON_IDX"
  grep -E "Uniquely mapped reads %|% of reads mapped to multiple loci" "${P}Log.final.out"
  grep -iE '"percent_mapped"|"num_processed"|"num_mapped"|decoy' "$SALMON_OUT/$S/aux_info/meta_info.json"
  grep '"expected_format"' "$SALMON_OUT/$S/lib_format_counts.json"
} > "$STAR_OUT/$S/${S}_run_info.txt"
cat "$STAR_OUT/$S/${S}_run_info.txt"
echo "Finished $S"
