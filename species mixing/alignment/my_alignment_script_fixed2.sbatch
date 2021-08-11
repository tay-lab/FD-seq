#!/bin/bash
#SBATCH --job-name=Van_fixed2
#SBATCH --output=Van_fixed2_sbatch.out
#SBATCH --error=Van_fixed2_sbatch.err
#SBATCH --partition=bigmem2
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mem=128G

#Change directory**
my_dir=/project2/tays/Van/20190212
mkdir -p $my_dir/Alignments/fixed
cd $my_dir/sbatch_files

#Make a folder for each sample in the Alignments directory**
# mkdir fixed

#Load picard and STAR**
module load picard
module load STAR

# #change fastq to bam files, using picard**
# java -jar $PICARD FastqToSam \
# F1=$my_dir/fastq/fixed_R1.fastq.gz \
# F2=$my_dir/fastq/fixed_R2.fastq.gz \
# O=$my_dir/Alignments/fixed/fixed.bam \
# SM=fixed_mix
#
# java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar TagBamWithReadSequenceExtended \
# INPUT=$my_dir/Alignments/fixed/fixed.bam \
# OUTPUT=$my_dir/Alignments/fixed/fixed_tag_cell.bam \
# SUMMARY=$my_dir/Alignments/fixed/fixed_tag_cell.bam_summary.txt \
# BASE_RANGE=1-12 \
# BASE_QUALITY=10 \
# BARCODED_READ=1 \
# DISCARD_READ=False \
# TAG_NAME=XC \
# NUM_BASES_BELOW_QUALITY=1
#
# java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar TagBamWithReadSequenceExtended \
# INPUT=$my_dir/Alignments/fixed/fixed_tag_cell.bam \
# OUTPUT=$my_dir/Alignments/fixed/fixed_tag_cellUMI.bam \
# SUMMARY=$my_dir/Alignments/fixed/fixed_tag_cellUMI.bam_summary.txt \
# BASE_RANGE=13-20 \
# BASE_QUALITY=10 \
# BARCODED_READ=1 \
# DISCARD_READ=True \
# TAG_NAME=XM \
# NUM_BASES_BELOW_QUALITY=1
#
# java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar FilterBAM \
# TAG_REJECT=XQ \
# INPUT=$my_dir/Alignments/fixed/fixed_tag_cellUMI.bam \
# OUTPUT=$my_dir/Alignments/fixed/fixed_tag_cellUMI_filtered.bam
#
# java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar TrimStartingSequence \
# INPUT=$my_dir/Alignments/fixed/fixed_tag_cellUMI_filtered.bam \
# OUTPUT=$my_dir/Alignments/fixed/fixed_tag_cellUMI_filtered_trimmed.bam \
# OUTPUT_SUMMARY=$my_dir/Alignments/fixed/fixed_tag_cellUMI_filtered_trimmed.bam_summary.txt \
# SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
# MISMATCHES=0 \
# NUM_BASES=5
#
# java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar PolyATrimmer \
# INPUT=$my_dir/Alignments/fixed/fixed_tag_cellUMI_filtered_trimmed.bam \
# OUTPUT=$my_dir/Alignments/fixed/fixed_tag_cellUMI_filtered_trimmed_pA.bam \
# OUTPUT_SUMMARY=$my_dir/Alignments/fixed/fixed_tag_cellUMI_filtered_trimmed_pA.bam_summary.txt \
# MISMATCHES=0 \
# NUM_BASES=6
#
# java -jar $PICARD SamToFastq \
# INPUT=$my_dir/Alignments/fixed/fixed_tag_cellUMI_filtered_trimmed_pA.bam \
# FASTQ=$my_dir/Alignments/fixed/fixed_for_mapping.fastq

STAR \
--runThreadN 10 \
--genomeDir /project2/tays/Van/reference_genomes/human_mouse_mix/genome_dir \
--readFilesIn $my_dir/Alignments/fixed/fixed_for_mapping.fastq \
--outFileNamePrefix $my_dir/Alignments/fixed/STAR_fixed_ \
--outSAMtype BAM Unsorted

java -jar $PICARD SortSam \
I=$my_dir/Alignments/fixed/STAR_fixed_Aligned.out.bam \
O=$my_dir/Alignments/fixed/STAR_fixed_Aligned.sorted.bam \
SO=queryname

java -jar $PICARD MergeBamAlignment \
REFERENCE_SEQUENCE=/project2/tays/Van/reference_genomes/human_mouse_mix/hg19_mm10_transgenes.fasta \
UNMAPPED_BAM=$my_dir/Alignments/fixed/fixed_tag_cellUMI_filtered_trimmed_pA.bam \
ALIGNED_BAM=$my_dir/Alignments/fixed/STAR_fixed_Aligned.sorted.bam \
OUTPUT=$my_dir/Alignments/fixed/STAR_fixed_merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar TagReadWithGeneExon \
I=$my_dir/Alignments/fixed/STAR_fixed_merged.bam \
O=$my_dir/Alignments/fixed/STAR_fixed_merged_exon_tagged.bam \
ANNOTATIONS_FILE=/project2/tays/Van/reference_genomes/human_mouse_mix/hg19_mm10_transgenes.gtf \
TAG=GE

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar DetectBeadSynthesisErrors \
I=$my_dir/Alignments/fixed/STAR_fixed_merged_exon_tagged.bam \
O=$my_dir/Alignments/fixed/STAR_fixed_merged_exon_tagged_clean.bam \
OUTPUT_STATS=$my_dir/Alignments/fixed/my.synthesis_stats.txt \
SUMMARY=$my_dir/Alignments/fixed/my.synthesis_stats.summary.txt \
NUM_BARCODES=10000 \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar FilterBAM \
I=$my_dir/Alignments/fixed/STAR_fixed_merged_exon_tagged_clean.bam \
O=$my_dir/Alignments/fixed/STAR_fixed_merged_exon_tagged_clean_human.bam \
REF_SOFT_MATCHED_RETAINED=HUMAN

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar FilterBAM \
I=$my_dir/Alignments/fixed/STAR_fixed_merged_exon_tagged_clean.bam \
O=$my_dir/Alignments/fixed/STAR_fixed_merged_exon_tagged_clean_mouse.bam \
REF_SOFT_MATCHED_RETAINED=MOUSE

# java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar BAMTagHistogram \
# I=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_exon_tagged_clean.bam \
# O=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_readcounts.txt.gz \
# TAG=XC
