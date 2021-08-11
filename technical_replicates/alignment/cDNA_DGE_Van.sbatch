#!/bin/bash
#SBATCH --job-name=Van_cDNA_DGE
#SBATCH --output=/project2/tays/Van/20210608/sbatch_filesVan_cDNA_DGE_sbatch.out
#SBATCH --error=/project2/tays/Van/20210608/sbatch_files/Van_cDNA_DGE_sbatch.err
#SBATCH --partition=broadwl
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=32G

# Directory variable of working Directory
my_dir=/project2/tays/Van/20210608

module load java

# Export to DGE
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar DigitalExpression \
I=$my_dir/Van_cDNA/STAR_cDNA_rep1_clean_substitution_synthesis.bam \
O=$my_dir/Van_cDNA/cDNA_rep1_dge.txt.gz \
SUMMARY=$my_dir/Van_cDNA/cDNA_rep1_summary_dge.txt \
CELL_BC_FILE=$my_dir/Van_cDNA/cDNA_rep1_knee_barcodes.txt

# Export to DGE
java -jar /project2/tays/Van/DropSeq_software/Drop-seq_tools-2.3.0/jar/dropseq.jar DigitalExpression \
I=$my_dir/Van_cDNA/STAR_cDNA_rep2_clean_substitution_synthesis.bam \
O=$my_dir/Van_cDNA/cDNA_rep2_dge.txt.gz \
SUMMARY=$my_dir/Van_cDNA/cDNA_rep2_summary_dge.txt \
CELL_BC_FILE=$my_dir/Van_cDNA/cDNA_rep2_knee_barcodes.txt
