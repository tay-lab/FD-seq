#!/bin/bash
#SBATCH --job-name=Van_preprocessing_20210608
#SBATCH --output=/project2/tays/Van/20210608/sbatch_files/Van_preprocessing_sbatch.out
#SBATCH --error=/project2/tays/Van/20210608/sbatch_files/Van_preprocessing_sbatch.err
#SBATCH --partition=broadwl
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G

# Directory variable of working Directory
my_dir=/project2/tays/Van/20210608/FastQ
cd $my_dir

# # Download from basespace
# module load bs-cli
# bs download project -i 270860608 -o /project2/tays/Van/20210608/FastQ/ -vvv

# rename bulk directories
for i in $my_dir/*_ds\.*
do
  mv $i/*.fastq.gz $my_dir
done

rm -rf $my_dir/*_ds\.*

# concatenate
for i1 in $my_dir/*_L001_*.fastq.gz
do
  i2="${i1/L001/L002}"
  i3="${i1/L001/L003}"
  i4="${i1/L001/L004}"
  ifinal="${i1/_L001_/_}"
  cat $i1 $i2 $i3 $i4 > $ifinal
  # echo $i1 $i2 $i3 $i4 $ifinal
done

rm -f $my_dir/*_L00[1-4]_*.fastq.gz
