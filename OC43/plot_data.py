# -*- coding: utf-8 -*-
"""
@author: Van Local

Plot data from fixed A549/OC43 drop-seq
"""
import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib as mpl

import os
my_dir = "Z:/Van/20200526 - Drop-seq fixed A549 OC43/data_analysis/"
os.chdir(my_dir)

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
# loc = mplticker.MultipleLocator(base=0.5)
# ggplot2_colors = ['#f8766dff','#7cae00ff','#00bfc4ff','#c77cffff']
ggplot2_colors = ['#f8766dff','#00ba38ff','#619cffff']
#*****

#%% 
# Import relative expression level of viral genes in infected cells
df = pd.read_csv(my_dir+"infected_relative_viral_genes.csv", index_col=0)
# Sort by percentage viral genes
df = df.loc[df["pct.viral"].sort_values().index,:]

# Import raw count of viral genes in all cells
df2 = pd.read_csv(my_dir+"all_samples_viral_genes.csv", index_col=0)
df2 = df2.T
viral_genes = df2.columns
df2["infected"] = df2.sum(axis=1) > 0
df2["sample"] = [s.split('_')[0] for s in df2.index]

# Import relative data of viral genes in moi1 cells
df3 = pd.read_csv(my_dir+"moi1_relative_viral_genes.csv", index_col=0)
# Sort by total.viral.log
df3 = df3.loc[df3["total.viral.log"].sort_values().index,:]

#%% Plot expression level of viral genes
# Heatmap
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.5,3))
sns.heatmap(df.loc[:,df.columns != 'pct.viral'].T, xticklabels=False, linewidth=0,
            cmap=sns.cubehelix_palette(light=0.95, dark=0.35, as_cmap=True),
            cbar_kws={"orientation":"horizontal", "label":"Expression level"}, ax=ax)
ax.tick_params(axis="both", length=0)
ax2 = ax.twinx()
ax2_color = 'tab:blue'
ax2.plot(range(df.shape[0]), df['pct.viral'].to_numpy(), color=ax2_color)
ax2.tick_params(axis='y', colors=ax2_color)
ax2.set_ylabel('Percentage of\nviral transcripts', color=ax2_color)
fig.tight_layout()
fig.savefig(my_dir+"heatmap_pct_viral.png", dpi=900, bbox_inches="tight")

#%% Plot % of infected cells
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(2.05,2.85))
for count, value in enumerate(['mock','MOI1']):
    temp = sum((df2['sample']==value)&(df2['infected']))/sum(df2['sample']==value)*100
    ax.bar(count, temp, width=0.5, color='grey')
    ax.text(count, temp+2, f'{temp:.1f}', horizontalalignment='center')
ax.set_xlabel("Sample")
ax.set_ylabel("Percentage of infected cells")
ax.set_ylim(0,80)
ax.set_xticks(range(2))
ax.set_xticklabels(['mock','MOI 1'])
sns.despine(fig=fig)
fig.savefig(my_dir+"pct_infected.svg", bbox_inches="tight")

#%% Plot relative expression level of viral genes in moi1
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
for count, value in enumerate(viral_genes):
    row = int(count/3)
    col = count % 3
    sns.scatterplot(x='total.viral.log', y=value, hue='cluster', data=df3,
                    palette=ggplot2_colors, s=7.5,
                    edgecolors=None, linewidth=0, alpha=0.8,
                    legend=False, ax=ax[row,col])
    # Plot rolling average of 50 cells
    ax[row,col].plot(df3['total.viral.log'].to_numpy(),
                     df3[value].rolling(window=50).mean().to_numpy(), c='k', alpha=0.75)
    ax[row,col].set_xlabel("log10 (total viral UMI + 1)")
    ax[row,col].set_ylabel("Expression level")
    ax[row,col].set_title(value)
sns.despine(fig=fig)
fig.tight_layout()
fig.savefig(my_dir+"8viral_genes_log.png", bbox_inches="tight", dpi=900)

#%% Plot KEGG pathway
# df_kegg = pd.DataFrame({"kegg":["TNF signaling","IL-17 signaling","NF-kB signaling", "Ribosome\nbiogenesis","Ribosome"],
#                         "p_val":[2.111e-7, 4.405e-7, 1.001e-6, 6.893e-15, 1.711e-38],
#                         "cluster":[0,0,0, 2, 3]})
df_kegg = pd.DataFrame({"kegg":["IL-17 signaling","TNF signaling","NF-kB signaling", "Ribosome\nbiogenesis"],
                        "p_val":[6.525e-8, 5.484e-7, 2.376e-6, 6.708e-13],
                        "cluster":[1,1,1, 2]})

# Bar plot of kegg p values
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(2.5,2.7))
for i, row in df_kegg.iterrows():
    ax.barh(i, -np.log10(row['p_val']), height=0.7, color=ggplot2_colors[row['cluster']])
ax.axvline(x=-np.log10(0.05), ls='--', color='k')
ax.set_xlabel("-log10 (adjusted p-value)")
ax.set_yticks(range(df_kegg.shape[0]))
ax.set_yticklabels(df_kegg['kegg'])
ax.invert_yaxis()
ax.set_xticks([0,5,10,15])
ax.tick_params(axis="y", length=0)
sns.despine(fig=fig)
fig.savefig(my_dir+"kegg_pval.svg", bbox_inches="tight")
