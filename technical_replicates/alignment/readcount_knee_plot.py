# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 09:54:42 2019

@author: Van Local
"""

# Import usual packages
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

# Set directory
my_dir = "Z:/Van/20210531 - FD-seq technical replicates/rcc_alignment/"

#%% rep 1 sample
# Read txt.gz file
df = pd.read_csv(my_dir+"STAR_cDNA_rep1_readcounts.txt.gz", sep='\t')

# Calculate normalized cumulative read counts
rc = df.iloc[:,0].values
rc_cumsum = np.cumsum(rc)
rc_cumsum = rc_cumsum/max(rc_cumsum)

# Plot read counts
n_cells = 1500
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(1,1,1)
ax.plot(rc_cumsum)
ax.axvline(x=n_cells, color='red', linestyle='--')
ax.set_xlim(0,10e3)
ax.set_xlabel('Cells')
ax.set_ylabel('Normalized cumulative read counts')

# Write the first n_cells cells
df.iloc[:n_cells,1].to_csv(my_dir+"cDNA_rep1_knee_barcodes.txt", index=False, header=False)

#%% rep 2 sample
# Read txt.gz file
df = pd.read_csv(my_dir+"STAR_cDNA_rep2_readcounts.txt.gz", sep='\t')

# Calculate normalized cumulative read counts
rc = df.iloc[:,0].values
rc_cumsum = np.cumsum(rc)
rc_cumsum = rc_cumsum/max(rc_cumsum)

# Plot read counts
n_cells = 1350
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(1,1,1)
ax.plot(rc_cumsum)
ax.axvline(x=n_cells, color='red', linestyle='--')
ax.set_xlim(0,10e3)
ax.set_xlabel('Cells')
ax.set_ylabel('Normalized cumulative read counts')

# Write the first n_cells cells
df.iloc[:n_cells,1].to_csv(my_dir+"cDNA_rep2_knee_barcodes.txt", index=False, header=False)
