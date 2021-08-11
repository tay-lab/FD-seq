# FD-seq
Fixed single-cell RNA sequencing for understanding virus infection and host response

This repository store the code used to align sequencing data, and analyze count data.

### Species-mixing data
* alignment folder: alignment code.
* mix_species_analysis.py: analysis of species mixing data.

### Technical replicates data
* alignment folder: alignment code.
* data_analysis.R and data_analysis.py: normalization and analysis of technical replicates data.

### KSHV data
* alignment folder: alignment code.
* Van_analysis_monocle_20200819.R: analysis of FD-seq single-cell data.
* viral_gene_correlation.py: calculate correlation among viral genes from FD-seq single-cell data.
* data_analysis.py: analysis of fluorescence intensity data.
* data_analysis_v3.py: analysis of qPCR data.

### OC43 data
* alignment folder: alignment code.
* 20200630.R: analysis of FD-seq single-cell data.
* plot_data.py: additional plotting of FD-seq single-cell data.

### RNA velocity analysis
* RNA velocity commands.txt: get spliced and unspliced count data using dropEst.
* positive_velocyto.loom and negative_velocyto.loom: spliced and unspliced count data.
* Van_analysis_loom.py: plotting of data processed by velocyto.
