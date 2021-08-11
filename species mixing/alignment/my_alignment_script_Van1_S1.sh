#!/bin/bash
#Change directory**
cd /project2/tays/Van/20190212/Alignments

#Make a folder for each sample in the Alignments directory**
mkdir AB-HE190129-Van1_S1

#Load picard and STAR**
module load picard
module load STAR

#change fastq to bam files, using picard**

java -jar $PICARD FastqToSam \
F1=/project2/tays/Van/20190212/fastq/AB-HE190129-Van1_S1_R1_001.fastq \
F2=/project2/tays/Van/20190212/fastq/AB-HE190129-Van1_S1_R2_001.fastq \
O=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1.bam \
SM=for_tool_testing

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar TagBamWithReadSequenceExtended \
INPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1.bam \
OUTPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_Cell.bam \
SUMMARY=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_Cell.bam_summary.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=False \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar TagBamWithReadSequenceExtended \
INPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_Cell.bam \
OUTPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI.bam \
SUMMARY=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI.bam_summary.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=True \
TAG_NAME=XM \
NUM_BASES_BELOW_QUALITY=1

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar FilterBAM \
TAG_REJECT=XQ \
INPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI.bam \
OUTPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI_filtered.bam

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar TrimStartingSequence \
INPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI_filtered.bam \
OUTPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI_trimmed.bam \
OUTPUT_SUMMARY=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI_trimmed.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar PolyATrimmer \
INPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI_trimmed.bam \
OUTPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI_trimmed_pA.bam \
OUTPUT_SUMMARY=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_pA_trimming_report.txt \
MISMATCHES=0 \
NUM_BASES=6

java -jar $PICARD SamToFastq \
INPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI_trimmed_pA.bam \
FASTQ=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_for_mapping.fastq

STAR \
--runThreadN 6 \
--genomeDir /project2/tays/Van/reference_genomes/human_mouse_mix/genome_dir \
--readFilesIn /project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_for_mapping.fastq \
--outFileNamePrefix /project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_ \
--outSAMtype BAM Unsorted

java -jar $PICARD SortSam \
I=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_Aligned.out.bam \
O=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_Aligned.sorted.bam \
SO=queryname

java -jar $PICARD MergeBamAlignment \
REFERENCE_SEQUENCE=/project2/tays/Van/reference_genomes/human_mouse_mix/hg19_mm10_transgenes.fasta \
UNMAPPED_BAM=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/AB-HE190129-Van1_S1_tag_CellUMI_trimmed_pA.bam \
ALIGNED_BAM=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_Aligned.sorted.bam \
OUTPUT=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_merged.bam \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar TagReadWithGeneExon \
I=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_merged.bam \
O=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_exon_tagged.bam \
ANNOTATIONS_FILE=/project2/tays/Van/reference_genomes/human_mouse_mix/hg19_mm10_transgenes.gtf \
TAG=GE

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar DetectBeadSynthesisErrors \
I=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_exon_tagged.bam \
O=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_exon_tagged_clean.bam \
OUTPUT_STATS=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/my.synthesis_stats.txt \
SUMMARY=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/my.synthesis_stats.summary.txt \
NUM_BARCODES= 10000 \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar FilterBAM \
I=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_exon_tagged_clean.bam \
O=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_exon_tagged_clean_human.bam \
REF_SOFT_MATCHED_RETAINED=HUMAN

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar FilterBAM \
I=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_exon_tagged_clean.bam \
O=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_exon_tagged_clean_mouse.bam \
REF_SOFT_MATCHED_RETAINED=MOUSE

java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-1.13/jar/dropseq.jar BAMTagHistogram \
I=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_exon_tagged_clean.bam \
O=/project2/tays/Van/20190212/Alignments/AB-HE190129-Van1_S1/STAR_AB-HE190129-Van1_S1_readcounts.txt.gz \
TAG=XC
