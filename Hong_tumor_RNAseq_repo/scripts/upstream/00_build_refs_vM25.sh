#!/bin/bash
#SBATCH --job-name=build_vM25_refs
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --partition=HTC
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=your_email
#SBATCH --mail-type=END,FAIL
# Build ALL references from ONE annotation (GENCODE vM25, GRCm38/mm10):
#   - Salmon decoy-aware index (vM25 transcripts + genome decoys, --gencode)
#   - STAR index (genome + vM25 GTF)
#   - tx2gene_vM25.tsv (transcript_id, gene_id, gene_name, gene_type)
# Submit:  sbatch 00_build_refs_vM25.sh
# If the compute node has no internet, run the wget line on the login node first.

source ~/.bashrc
source /ix1/hzhang/SteveWang/TEtranscripts_venv/bin/activate
conda activate /vast/bnacev/ziw183/conda_envs/nacevlab
set -eo pipefail

# --- Configuration ---
REF=/ix1/bnacev/ziw183/ref/grcm38
GENOME=$REF/GRCm38.primary_assembly.genome.fa
GTF=$REF/gencode.vM25.annotation.gtf
TX_FA=$REF/gencode.vM25.transcripts.fa.gz
TX_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz
OUT=$REF/vM25
SJDB_OVERHANG=149          # read length - 1 (PE150); 100 also works fine
THREADS=$SLURM_CPUS_PER_TASK
# --- End of Configuration ---

for t in STAR salmon samtools; do command -v $t >/dev/null || { echo "ERROR: $t not in env"; exit 1; }; done
mkdir -p "$OUT"

# 1. vM25 transcript FASTA
[ -s "$TX_FA" ] || wget -O "$TX_FA" "$TX_URL"

# 2. tx2gene (versions kept: they match Salmon --gencode names exactly)
awk -F'\t' '$3=="transcript"{
    match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16)
    match($9,/gene_id "[^"]+"/);       g=substr($9,RSTART+9, RLENGTH-10)
    match($9,/gene_name "[^"]+"/);     n=substr($9,RSTART+11,RLENGTH-12)
    match($9,/gene_type "[^"]+"/);     y=substr($9,RSTART+11,RLENGTH-12)
    print t"\t"g"\t"n"\t"y }' "$GTF" | sort -u > "$OUT/tx2gene_vM25.tsv"
sed -i '1i transcript_id\tgene_id\tgene_name\tgene_type' "$OUT/tx2gene_vM25.tsv"
echo "tx2gene rows: $(($(wc -l < "$OUT/tx2gene_vM25.tsv") - 1))"

# 3. Salmon decoy-aware index
grep '^>' "$GENOME" | cut -d ' ' -f1 | sed 's/^>//' > "$OUT/decoys_vM25.txt"
cat <(zcat "$TX_FA") "$GENOME" > "$OUT/gentrome_vM25.fa"
salmon index -t "$OUT/gentrome_vM25.fa" -d "$OUT/decoys_vM25.txt" \
             -i "$OUT/salmon_index_vM25_decoy" -k 31 --gencode -p "$THREADS"
rm -f "$OUT/gentrome_vM25.fa"

# 4. STAR index
mkdir -p "$OUT/STAR_index_vM25"
STAR --runMode genomeGenerate \
     --runThreadN "$THREADS" \
     --genomeDir "$OUT/STAR_index_vM25" \
     --genomeFastaFiles "$GENOME" \
     --sjdbGTFfile "$GTF" \
     --sjdbOverhang "$SJDB_OVERHANG" \
     --limitGenomeGenerateRAM 80000000000

# 5. Record what was built
{
  echo "built: $(date)"
  echo "genome: $GENOME"; echo "gtf: $GTF"; echo "transcripts: $TX_FA ($TX_URL)"
  echo "salmon $(salmon --version 2>&1) | STAR $(STAR --version)"
  echo "salmon index: -k 31 --gencode, decoys = all genome sequences"
  echo "STAR index: --sjdbOverhang $SJDB_OVERHANG"
  echo "transcripts in FASTA: $(zcat "$TX_FA" | grep -c '^>')"
} > "$OUT/BUILD_INFO_vM25.txt"
cat "$OUT/BUILD_INFO_vM25.txt"
