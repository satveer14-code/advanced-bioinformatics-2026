 #!/bin/bash
set -euo pipefail

#######################################
# Section 2 NGS pipeline
#
# This script reproduces the successful workflow used in Section 2.
# It performs:
#   1. Raw read QC
#   2. Read trimming
#   3. Post-trim QC
#   4. Reference indexing
#   5. Alignment with BWA-MEM
#   6. BAM conversion, sorting, indexing
#   7. Duplicate marking
#   8. BAM filtering
#   9. Alignment statistics
#  10. Variant calling with Freebayes
#  11. Variant compression and indexing
#  12. Variant filtering
#  13. ANNOVAR conversion
#  14. ANNOVAR annotation
#  15. SnpEff annotation
#  16. Basic prioritisation
#
# This script ensures:
# - required tools are checked 
# - required input files are checked
# - skips completed steps when outputs already exist
# - writes all console output to a log file

#######################################
# Logging setup
#######################################

# Create a timestamp for the run log
RUN_TIME="$(date +%Y%m%d_%H%M%S)"

# Temporary default log location until BASE_DIR is defined
TMP_LOG="$HOME/section2_ngs_pipeline_${RUN_TIME}.log"

# Send all output to both screen and log file
exec > >(tee -a "$TMP_LOG") 2>&1

echo "Starting Section 2 NGS pipeline at $(date)"
echo "Temporary log file: $TMP_LOG"

#######################################
# Define project paths
#######################################

# Base project directory
BASE_DIR="$HOME/ngs_assessment"

# Input data directories
RAW_DIR="$BASE_DIR/data/raw"
TRIM_DIR="$BASE_DIR/data/trimmed"
REF_DIR="$BASE_DIR/data/reference"
BED_DIR="$BASE_DIR/data/annotation"

# Output directories
RESULTS_DIR="$BASE_DIR/results"
FASTQC_RAW_DIR="$RESULTS_DIR/fastqc_raw"
FASTQC_TRIM_DIR="$RESULTS_DIR/fastqc_trimmed"
ALIGN_DIR="$RESULTS_DIR/aligned"
VARIANT_DIR="$RESULTS_DIR/variants"
ANNOT_OUT_DIR="$RESULTS_DIR/annotation"
LOG_DIR="$BASE_DIR/logs"

# Main log file inside project
PIPELINE_LOG="$LOG_DIR/section2_ngs_pipeline_${RUN_TIME}.log"

# Ensure log directory exists, then switch logging to project log
mkdir -p "$LOG_DIR"
exec > >(tee -a "$PIPELINE_LOG") 2>&1

echo "Project log file: $PIPELINE_LOG"

#######################################
# Define input files
#######################################

# Raw paired-end FASTQ files
R1_RAW="$RAW_DIR/NGS0001.R1.fastq.gz"
R2_RAW="$RAW_DIR/NGS0001.R2.fastq.gz"

# Reference FASTA and BED target file
REF_FASTA="$REF_DIR/hg19.fa"
TARGET_BED="$BED_DIR/annotation.bed"

#######################################
# Define output files
#######################################

# Trimmed FASTQ files
R1_PAIRED="$TRIM_DIR/NGS0001_R1_paired.fastq.gz"
R1_UNPAIRED="$TRIM_DIR/NGS0001_R1_unpaired.fastq.gz"
R2_PAIRED="$TRIM_DIR/NGS0001_R2_paired.fastq.gz"
R2_UNPAIRED="$TRIM_DIR/NGS0001_R2_unpaired.fastq.gz"

# Alignment files
SAM_FILE="$ALIGN_DIR/NGS0001.sam"
BAM_FILE="$ALIGN_DIR/NGS0001.bam"
SORTED_BAM="$ALIGN_DIR/NGS0001_sorted.bam"
MARKED_BAM="$ALIGN_DIR/NGS0001_sorted_marked.bam"
FILTERED_BAM="$ALIGN_DIR/NGS0001_sorted_filtered.bam"

echo "Using existing filtered BAM file: $FILTERED_BAM"

# Variant files
RAW_VCF="$VARIANT_DIR/NGS0001_bcftools_raw.vcf"
RAW_VCF_GZ="$VARIANT_DIR/NGS0001_bcftools_raw.vcf.gz"
RAW_VCF_TBI="$RAW_VCF_GZ.tbi"
FILTERED_VCF="$VARIANT_DIR/NGS0001_bcftools_filtered.vcf"

# ANNOVAR files
ANNOVAR_DIR="$HOME/annovar"
ANNOVAR_DB="$ANNOVAR_DIR/humandb"
ANNOVAR_INPUT="$ANNOT_OUT_DIR/NGS0001.avinput"
ANNOVAR_PREFIX="$ANNOT_OUT_DIR/NGS0001_annovar"
ANNOVAR_CSV="${ANNOVAR_PREFIX}.hg19_multianno.csv"

# SnpEff files
SNPEFF_DIR="$HOME/snpEff"
SNPEFF_JAR="$SNPEFF_DIR/snpEff.jar"
SNPEFF_DATA_DIR="$SNPEFF_DIR/data"
SNPEFF_GENOME="GRCh37.75"
SNPEFF_VCF="$ANNOT_OUT_DIR/NGS0001_snpeff.vcf"

# Prioritised output
PRIORITISED_CSV="$ANNOT_OUT_DIR/NGS0001_prioritized.csv"

# Statistics files
FLAGSTAT_TXT="$LOG_DIR/NGS0001_flagstat.txt"
IDXSTATS_TXT="$LOG_DIR/NGS0001_idxstats.txt"
DEPTH_TXT="$LOG_DIR/NGS0001_depth.txt"
MARKDUP_METRICS="$LOG_DIR/NGS0001_markdup_metrics.txt"
INSERT_METRICS="$LOG_DIR/NGS0001_insert_metrics.txt"
INSERT_HISTOGRAM="$LOG_DIR/NGS0001_insert_histogram.pdf"

# Trimmomatic adapters file that worked in your run
ADAPTERS="/home/ubuntu/anaconda3/pkgs/trimmomatic-0.40-hfd78af_0/share/trimmomatic-0.40-0/adapters/NexteraPE-PE.fa"

#######################################
# Create output directories
#######################################

mkdir -p "$FASTQC_RAW_DIR" \
         "$FASTQC_TRIM_DIR" \
         "$ALIGN_DIR" \
         "$VARIANT_DIR" \
         "$ANNOT_OUT_DIR" \
         "$LOG_DIR" \
         "$TRIM_DIR"

echo "Output directories checked/created."

#######################################
# Helper functions
#######################################

# Print section headers to make the log easier to read
print_section() {
    echo
    echo "=================================================="
    echo "$1"
    echo "=================================================="
}

# Exit with a clear error message
die() {
    echo "ERROR: $1" >&2
    exit 1
}

# Check that a command exists in PATH
check_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

# Check that a file exists
check_file() {
    [[ -f "$1" ]] || die "Required file not found: $1"
}

#######################################
# Check required commands
#######################################

print_section "Checking required tools"

check_command fastqc
check_command trimmomatic
check_command bwa
check_command samtools
check_command picard
check_command bcftools
check_command bgzip
check_command tabix
check_command vcffilter
check_command java
check_command awk

echo "All required command-line tools were found."

#######################################
# Check required input files
#######################################

print_section "Checking required input files"

# Tools/resources always needed for downstream annotation
check_file "$ANNOVAR_DIR/convert2annovar.pl"
check_file "$ANNOVAR_DIR/table_annovar.pl"
check_file "$ANNOVAR_DB/hg19_refGene.txt"
check_file "$ANNOVAR_DB/hg19_ensGene.txt"
check_file "$SNPEFF_JAR"
check_file "$SNPEFF_DATA_DIR/$SNPEFF_GENOME/snpEffectPredictor.bin"

# Choose the highest-level available checkpoint
if [[ -f "$FILTERED_VCF" ]]; then
    echo "Existing filtered VCF detected. Pipeline will resume from filtered VCF."
elif [[ -f "$FILTERED_BAM" && -f "${FILTERED_BAM}.bai" ]]; then
    echo "Existing filtered BAM detected. Pipeline will resume from BAM."
    check_file "$REF_FASTA"
    check_file "$TARGET_BED"
elif [[ -f "$R1_RAW" && -f "$R2_RAW" ]]; then
    echo "Raw FASTQ inputs detected. Pipeline will run from raw reads."
    check_file "$REF_FASTA"
    check_file "$TARGET_BED"
    check_file "$ADAPTERS"
else
    die "No usable input found. Need either filtered VCF, filtered BAM + index, or raw FASTQ files."
fi

echo "Required input files were found."

# Decide where to resume from
if [[ -f "$FILTERED_VCF" ]]; then
    RESUME_FROM="vcf"
elif [[ -f "$FILTERED_BAM" && -f "${FILTERED_BAM}.bai" ]]; then
    RESUME_FROM="bam"
else
    RESUME_FROM="raw"
fi

echo "Resume mode: $RESUME_FROM"

if [[ "$RESUME_FROM" == "raw" ]]; then

#######################################
# Step 1: FastQC on raw reads
#######################################

print_section "Step 1: FastQC on raw reads"

RAW_QC1="$FASTQC_RAW_DIR/NGS0001.R1_fastqc.zip"
RAW_QC2="$FASTQC_RAW_DIR/NGS0001.R2_fastqc.zip"

if [[ -f "$RAW_QC1" && -f "$RAW_QC2" ]]; then
    echo "Raw FastQC outputs already exist. Skipping Step 1."
else
    fastqc -t 4 "$R1_RAW" "$R2_RAW" -o "$FASTQC_RAW_DIR"
    echo "Raw read FastQC completed."
fi

#######################################
# Step 2: Trim reads
#######################################

print_section "Step 2: Trimming reads"

if [[ -f "$R1_PAIRED" && -f "$R1_UNPAIRED" && -f "$R2_PAIRED" && -f "$R2_UNPAIRED" ]]; then
    echo "Trimmed FASTQ outputs already exist. Skipping Step 2."
else
    trimmomatic PE -threads 4 -phred33 \
      "$R1_RAW" "$R2_RAW" \
      "$R1_PAIRED" "$R1_UNPAIRED" \
      "$R2_PAIRED" "$R2_UNPAIRED" \
      ILLUMINACLIP:"$ADAPTERS":2:30:10 \
      TRAILING:25 MINLEN:50
    echo "Read trimming completed."
fi

#######################################
# Step 3: FastQC on trimmed reads
#######################################

print_section "Step 3: FastQC on trimmed reads"

TRIM_QC1="$FASTQC_TRIM_DIR/NGS0001_R1_paired_fastqc.zip"
TRIM_QC2="$FASTQC_TRIM_DIR/NGS0001_R2_paired_fastqc.zip"

if [[ -f "$TRIM_QC1" && -f "$TRIM_QC2" ]]; then
    echo "Trimmed FastQC outputs already exist. Skipping Step 3."
else
    fastqc -t 4 "$R1_PAIRED" "$R2_PAIRED" -o "$FASTQC_TRIM_DIR"
    echo "Trimmed read FastQC completed."
fi

#######################################
# Step 4: Reference indexing
#######################################

print_section "Step 4: Reference indexing"

# Build BWA index if missing
if [[ -f "${REF_FASTA}.bwt" ]]; then
    echo "BWA index already exists. Skipping BWA indexing."
else
    bwa index "$REF_FASTA"
    echo "BWA index completed."
fi

# Build samtools FASTA index if missing
if [[ -f "${REF_FASTA}.fai" ]]; then
    echo "FASTA index already exists. Skipping samtools faidx."
else
    samtools faidx "$REF_FASTA"
    echo "samtools faidx completed."
fi

#######################################
# Step 5: Align reads with BWA-MEM
#######################################

print_section "Step 5: Align reads with BWA-MEM"

if [[ -f "$SAM_FILE" ]]; then
    echo "SAM alignment file already exists. Skipping Step 5."
else
    bwa mem \
      -R '@RG\tID:NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera\tPU:unit1' \
      "$REF_FASTA" \
      "$R1_PAIRED" \
      "$R2_PAIRED" \
      > "$SAM_FILE"
    echo "Alignment completed."
fi

#######################################
# Step 6: Convert SAM to BAM, sort, index
#######################################

print_section "Step 6: Convert, sort and index BAM"

if [[ -f "$SORTED_BAM" && -f "${SORTED_BAM}.bai" ]]; then
    echo "Sorted BAM and index already exist. Skipping Step 6."
else
    [[ -f "$BAM_FILE" ]] || samtools view -Sb "$SAM_FILE" > "$BAM_FILE"
    samtools sort "$BAM_FILE" -o "$SORTED_BAM"
    samtools index "$SORTED_BAM"
    echo "BAM conversion, sorting and indexing completed."
fi

#######################################
# Step 7: Mark duplicates
#######################################

print_section "Step 7: Mark duplicates"

if [[ -f "$MARKED_BAM" && -f "${MARKED_BAM}.bai" && -f "$MARKDUP_METRICS" ]]; then
    echo "Duplicate-marked BAM already exists. Skipping Step 7."
else
    picard MarkDuplicates \
      I="$SORTED_BAM" \
      O="$MARKED_BAM" \
      M="$MARKDUP_METRICS"

    samtools index "$MARKED_BAM"
    echo "Duplicate marking completed."
fi

#######################################
# Step 8: Filter BAM
#######################################

print_section "Step 8: Filter BAM"

if [[ -f "$FILTERED_BAM" && -f "${FILTERED_BAM}.bai" ]]; then
    echo "Filtered BAM already exists. Skipping Step 8."
else
    samtools view \
      -h -b \
      -q 20 \
      -F 1796 \
      "$MARKED_BAM" \
      > "$FILTERED_BAM"

    samtools index "$FILTERED_BAM"
    echo "Filtered BAM completed."
fi

fi

if [[ "$RESUME_FROM" == "raw" || "$RESUME_FROM" == "bam" ]]; then

#######################################
# Step 9: Generate alignment statistics
#######################################

print_section "Step 9: Generate alignment statistics"

[[ -f "$FLAGSTAT_TXT" ]] || samtools flagstat "$FILTERED_BAM" > "$FLAGSTAT_TXT"
[[ -f "$IDXSTATS_TXT" ]] || samtools idxstats "$FILTERED_BAM" > "$IDXSTATS_TXT"
[[ -f "$DEPTH_TXT" ]] || samtools depth "$FILTERED_BAM" > "$DEPTH_TXT"

if [[ -f "$INSERT_METRICS" && -f "$INSERT_HISTOGRAM" ]]; then
    echo "Insert size metrics already exist."
else
    picard CollectInsertSizeMetrics \
      I="$FILTERED_BAM" \
      O="$INSERT_METRICS" \
      H="$INSERT_HISTOGRAM" \
      M=0.5
fi

echo "Alignment statistics completed."

#######################################
# Step 10: Variant calling with bcftools
#######################################

print_section "Step 10: Variant calling with bcftools mpileup + bcftools call"

if [[ -f "$RAW_VCF_GZ" && -f "$RAW_VCF_TBI" ]]; then
    echo "Compressed/indexed bcftools raw VCF already exists. Skipping raw VCF generation."
elif [[ -f "$RAW_VCF" ]]; then
    echo "bcftools raw VCF already exists. Compression/indexing will be handled in Step 11."
else
    bcftools mpileup \
      -Ou \
      -f "$REF_FASTA" \
      -R "$TARGET_BED" \
      "$FILTERED_BAM" \
    | bcftools call \
      -mv \
      -Ov \
      -o "$RAW_VCF"

    echo "bcftools variant calling completed."
fi

#######################################
# Step 11: Compress and index raw VCF
#######################################

print_section "Step 11: Compress and index bcftools raw VCF"

if [[ -f "$RAW_VCF_GZ" && -f "$RAW_VCF_TBI" ]]; then
    echo "Compressed/indexed bcftools raw VCF already exists. Skipping Step 11."
else
    bgzip -f "$RAW_VCF"
    tabix -p vcf "$RAW_VCF_GZ"
    echo "bcftools raw VCF compression and indexing completed."
fi

#######################################
# Step 12: Filter variants
#######################################

print_section "Step 12: Filter variants"

if [[ -f "$FILTERED_VCF" ]]; then
    echo "Filtered VCF already exists. Skipping Step 12."
else
    vcffilter \
      -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
      "$RAW_VCF_GZ" \
      > "$FILTERED_VCF"
    echo "Variant filtering completed."
fi

fi

#######################################
# Step 13: Convert VCF to ANNOVAR input
#######################################

print_section "Step 13: Convert VCF to ANNOVAR input"

if [[ -f "$ANNOVAR_INPUT" ]]; then
    echo "ANNOVAR input already exists. Skipping Step 13."
else
    "$ANNOVAR_DIR/convert2annovar.pl" -format vcf4 \
      "$FILTERED_VCF" \
      > "$ANNOVAR_INPUT"
    echo "ANNOVAR input conversion completed."
fi

#######################################
# Step 14: ANNOVAR annotation
#######################################

print_section "Step 14: ANNOVAR annotation"

if [[ -f "$ANNOVAR_CSV" ]]; then
    echo "ANNOVAR CSV output already exists. Skipping Step 14."
else
    "$ANNOVAR_DIR/table_annovar.pl" \
      "$ANNOVAR_INPUT" \
      "$ANNOVAR_DB/" \
	      -buildver hg19 \
      -out "$ANNOVAR_PREFIX" \
      -remove \
      -protocol refGene,ensGene \
      -operation g,g \
      -otherinfo \
      -nastring . \
      -csvout
    echo "ANNOVAR annotation completed."
fi

#######################################
# Step 15: SnpEff annotation
#######################################

print_section "Step 15: SnpEff annotation"

if [[ -f "$SNPEFF_VCF" ]]; then
    echo "SnpEff VCF already exists. Skipping Step 15."
else
    java -Xmx4g -jar "$SNPEFF_JAR" \
      -dataDir "$SNPEFF_DATA_DIR" \
      "$SNPEFF_GENOME" \
      "$FILTERED_VCF" \
      > "$SNPEFF_VCF"
    echo "SnpEff annotation completed."
fi

#######################################
# Step 16: Basic prioritisation
#######################################

print_section "Step 16: Basic prioritisation"

awk -F',' 'NR==1 || ($6=="exonic" && $14==".")' \
"$ANNOVAR_CSV" > "$PRIORITISED_CSV"

echo "Basic prioritisation completed."

#######################################
# Step 17: Final summary
#######################################

print_section "Pipeline completed"

echo "Key outputs:"
echo "  Filtered BAM:        $FILTERED_BAM"
echo "  Filtered VCF:        $FILTERED_VCF"
echo "  ANNOVAR CSV:         $ANNOVAR_CSV"
echo "  SnpEff VCF:          $SNPEFF_VCF"
echo "  Prioritised CSV:     $PRIORITISED_CSV"
echo "  Pipeline log:        $PIPELINE_LOG"

echo
echo "Line count of prioritised CSV:"
wc -l "$PRIORITISED_CSV"

echo
echo "Section 2 NGS pipeline finished successfully at $(date)"
