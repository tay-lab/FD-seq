# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:15:53 2018

@author: Van Local
"""

# Import packages here
import pandas as pd
import numpy as np
import scipy.stats as stats

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker
import seaborn as sns

#*****
mpl.rcdefaults()
# Set font to be arial
mpl.rc('font', **{'sans-serif':'Arial', 'size':12})
mpl.rcParams['mathtext.rm'] = 'sans' # to have non-italic greek letter, use r'$\mathrm{\alpha}$', does NOT work with f-string
mpl.rcParams['axes.titlesize'] = 12
mpl.rcParams['hatch.linewidth'] = 1.5  # previous svg hatch linewidth
# Set default tick size
mpl.rcParams['xtick.major.size'] = 5.5
mpl.rcParams['ytick.major.size'] = 5.5
mpl.rcParams['xtick.minor.size'] = 2.5
mpl.rcParams['ytick.minor.size'] = 2.5
# Set default seaborn palette
sns.set_palette('Dark2')
#*****

# Set directory
my_dir = 'Z:/Van/20210531 - FD-seq technical replicates/data_analysis/'


# =============================================================================
# # Continue data analysis from monocle: data_analysis.R
# =============================================================================
# Import normalized data
df_all = pd.read_csv(my_dir+'all_norm_gene_expr.txt.gz', sep='\t', index_col=0)

# Set the minimum expression level to 0.1
df_all[df_all < 0.1] = 0.1

# First 1483 cells: rep1
# Next 1338 cells: rep2
df1 = df_all.iloc[:,:1483]
df2 = df_all.iloc[:,1483:]

# Plot mean values
fig, ax = plt.subplots(figsize=(3.3,3.3), tight_layout=True)
ax.scatter(df1.mean(axis=1), df2.mean(axis=1), color='grey', edgecolor='none', s=10)
ax.plot([0.9e-1,1e2], [0.9e-1,1e2], color='red')
rho = stats.pearsonr(np.log10(df1.mean(axis=1)), np.log10(df2.mean(axis=1)))
ax.text(0.07, 0.9, r'$\mathrm{\rho}$ = '+f'{rho[0]:.3f}',
        transform=ax.transAxes)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks([1e-1,1e0,1e1,1e2])
ax.set_xlim(.7e-1, 1.5e2)
ax.set_ylim(.7e-1, 1.5e2)
ax.minorticks_off()
ax.set_xlabel("Normalized expression\n(technical replicate 1)")
ax.set_ylabel("Normalized expression\n(technical replicate 2)")
sns.despine(fig=fig)
fig.savefig(my_dir+'technical_replicate_expression.png', bbox_inches='tight', dpi=600, pad_inches=0)