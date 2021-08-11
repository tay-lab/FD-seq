#!/bin/bash
#SBATCH --job-name=Van_MOI1_alignment
#SBATCH --output=Van_MOI1_alignment_sbatch.out
#SBATCH --error=Van_MOI1_alignment_sbatch.err
#SBATCH --partition=broadwl
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mem-per-cpu=4000

#Directory variable of working Directory
my_dir=/project2/tays/Van/20200526

#Change directory**
mkdir -p $my_dir/MOI1
cd $my_dir/MOI1

# Load picard and STAR**
module load java
module load picard
module load STAR

# change fastq.gz to bam files, using picard
java -jar $PICARD FastqToSam \
F1=$my_dir/fastq/ST-HVP-3S-DROP-HVP-MOI-1_S1_R1_001.fastq.gz \
F2=$my_dir/fastq/ST-HVP-3S-DROP-HVP-MOI-1_S1_R2_001.fastq.gz \
O=$my_dir/MOI1/MOI1.bam \
SM=OC43_MOI1 \
TMP_DIR=$my_dir/sbatch_files

# Tag cell barcode
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar TagBamWithReadSequenceExtended \
INPUT=$my_dir/MOI1/MOI1.bam \
OUTPUT=$my_dir/MOI1/MOI1_tag_Cell.bam \
SUMMARY=$my_dir/MOI1/MOI1_tag_Cell_summary.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1

# Tag UMI (molecular barcode)
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar TagBamWithReadSequenceExtended \
INPUT=$my_dir/MOI1/MOI1_tag_Cell.bam \
OUTPUT=$my_dir/MOI1/MOI1_tag_CellUMI.bam \
SUMMARY=$my_dir/MOI1/MOI1_tag_CellUMI_summary.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=True \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1

# Quality filter
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar FilterBam \
TAG_REJECT=XQ \
INPUT=$my_dir/MOI1/MOI1_tag_CellUMI.bam \
OUTPUT=$my_dir/MOI1/MOI1_tag_filtered.bam

# Trim SMART adapter sequence at 5'
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar TrimStartingSequence \
INPUT=$my_dir/MOI1/MOI1_tag_filtered.bam \
OUTPUT=$my_dir/MOI1/MOI1_tag_filtered_trimmed_adapter.bam \
OUTPUT_SUMMARY=$my_dir/MOI1/MOI1_tag_filtered_trimmed_adapter_report.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5

# Trim poly A at 3'
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar PolyATrimmer \
INPUT=$my_dir/MOI1/MOI1_tag_filtered_trimmed_adapter.bam \
OUTPUT=$my_dir/MOI1/MOI1_tag_filtered_trimmed_adapter_pA.bam \
OUTPUT_SUMMARY=$my_dir/MOI1/MOI1_tag_filtered_trimmed_adapter_pA_report.txt \
MISMATCHES=0 \
NUM_BASES=6 \
USE_NEW_TRIMMER=true

java -jar $PICARD SamToFastq \
INPUT=$my_dir/MOI1/MOI1_tag_filtered_trimmed_adapter_pA.bam \
FASTQ=$my_dir/MOI1/MOI1_for_mapping.fastq \
TMP_DIR=$my_dir/sbatch_files

STAR \
--runThreadN 12 \
--genomeDir /project2/tays/Van/reference_genomes/human_OC43_genome/genomeDir \
--readFilesIn $my_dir/MOI1/MOI1_for_mapping.fastq \
--outFileNamePrefix $my_dir/MOI1/STAR_MOI1_ \
--outSAMtype BAM Unsorted

java -jar $PICARD SortSam \
I=$my_dir/MOI1/STAR_MOI1_Aligned.out.bam \
O=$my_dir/MOI1/STAR_MOI1_Aligned.sorted.bam \
SO=queryname \
TMP_DIR=$my_dir/sbatch_files

java -jar $PICARD MergeBamAlignment \
REFERENCE_SEQUENCE=/project2/tays/Van/reference_genomes/human_OC43_genome/GRCh38_OC43.fa \
UNMAPPED_BAM=$my_dir/MOI1/MOI1_tag_filtered_trimmed_adapter_pA.bam \
ALIGNED_BAM=$my_dir/MOI1/STAR_MOI1_Aligned.sorted.bam \
OUTPUT=$my_dir/MOI1/STAR_MOI1_merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false \
TMP_DIR=$my_dir/sbatch_files

# # TagReadWithGeneExon from v1 dropseq tool
# java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar TagReadWithGeneExon \
# I=$my_dir/MOI1/STAR_MOI1_merged.bam \
# O=$my_dir/MOI1/STAR_MOI1_exon_tagged.bam \
# ANNOTATIONS_FILE=/project2/tays/Van/reference_genomes/human_GRCh38/Human_reference_genome/gencode.v27.annotation.gtf \
# TAG=GE

# TagReadWithGeneFunction from v2 dropseq tool
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar TagReadWithGeneFunction \
I=$my_dir/MOI1/STAR_MOI1_merged.bam \
O=$my_dir/MOI1/STAR_MOI1_exon_tagged.bam \
ANNOTATIONS_FILE=/project2/tays/Van/reference_genomes/human_OC43_genome/GRCh38_OC43.gtf

# Detect bead barcode substitution error
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar DetectBeadSubstitutionErrors \
I=$my_dir/MOI1/STAR_MOI1_exon_tagged.bam \
O=$my_dir/MOI1/STAR_MOI1_clean_substitution.bam \
MIN_UMIS_PER_CELL=20 \
OUTPUT_REPORT=$my_dir/MOI1/STAR_MOI1_clean_substitution_report.txt

# Detect bead barcode synthesis error
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar DetectBeadSynthesisErrors \
I=$my_dir/MOI1/STAR_MOI1_clean_substitution.bam \
O=$my_dir/MOI1/STAR_MOI1_clean_substitution_synthesis.bam \
MIN_UMIS_PER_CELL=20 \
REPORT=$my_dir/MOI1/STAR_MOI1_synthesis_error_report.txt \
OUTPUT_STATS=$my_dir/MOI1/_synthesis_error_stats.txt \
SUMMARY=$my_dir/MOI1/_synthesis_error_summary.txt \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

# Histogram for knee plot
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar BamTagHistogram \
I=$my_dir/MOI1/STAR_MOI1_clean_substitution_synthesis.bam \
O=$my_dir/MOI1/STAR_MOI1_readcounts.txt.gz \
TAG=XC


#THIS PART IN PYTHON****************************
#Get the cell barcodes from the knee plot
#Export the barcodes as text file: SEEneg-cDNA_knee_barcodes.txt
#END PART IN PYTHON****************************

# # Export to DGE
# java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar DigitalExpression \
# I=$my_dir/MOI1/STAR_MOI1_clean_substitution_synthesis.bam \
# O=$my_dir/MOI1/MOI1_dge.txt.gz \
# SUMMARY=$my_dir/MOI1/MOI1_summary_dge.txt \
# CELL_BC_FILE=$my_dir/MOI1/STAR_MOI1_knee_barcodes.txt
