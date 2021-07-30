# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:15:53 2018

@author: Van Local
"""

# Import packages here
import pandas as pd
import numpy as np
import math

import scipy.cluster.hierarchy as sch
import scipy.spatial as sp
import scipy.stats as stats
from scipy.optimize import curve_fit

import matplotlib as mpl
import matplotlib.pyplot as plt
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
my_dir = 'Z:/Van/20190212 - Species mix live and fixed cells/'

# =============================================================================
# # Make mix species plot
# =============================================================================
# Import data as data frame
df1_human = pd.read_csv(my_dir+'STAR_AB-HE190129-Van1_S1_human.dge.txt.gz', sep='\t', index_col=0) # Live
df1_mouse = pd.read_csv(my_dir+'STAR_AB-HE190129-Van1_S1_mouse.dge.txt.gz', sep='\t', index_col=0) # Live
df2_human = pd.read_csv(my_dir+'STAR_AB-HE190129-Van2_S2_human.dge.txt.gz', sep='\t', index_col=0) # Fixed
df2_mouse = pd.read_csv(my_dir+'STAR_AB-HE190129-Van2_S2_mouse.dge.txt.gz', sep='\t', index_col=0) # Fixed

# Calculate species transcripts
df1_UMI = pd.DataFrame(data={'human':df1_human.sum(axis=0), 'mouse':df1_mouse.sum(axis=0)})
df2_UMI = pd.DataFrame(data={'human':df2_human.sum(axis=0), 'mouse':df2_mouse.sum(axis=0)})

# Label the cells by species
# More than 90% of transcripts coming from a species = species-specific
df1_UMI['label'] = 'Mixed'
df1_UMI.loc[df1_UMI['human']/(df1_UMI['human']+df1_UMI['mouse'])>0.9,'label'] = 'Human'
df1_UMI.loc[df1_UMI['mouse']/(df1_UMI['human']+df1_UMI['mouse'])>0.9,'label'] = 'Mouse'
df2_UMI['label'] = 'Mixed'
df2_UMI.loc[df2_UMI['human']/(df2_UMI['human']+df2_UMI['mouse'])>0.9,'label'] = 'Human'
df2_UMI.loc[df2_UMI['mouse']/(df2_UMI['human']+df2_UMI['mouse'])>0.9,'label'] = 'Mouse'

# Plot
#my_color = sns.color_palette()[0:3]
fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(1,2,1)
#for (counter, value) in enumerate(['Human','Mouse','Mixed']):
#    ax1.scatter(x=df1_UMI.loc[df1_UMI['label']==value,'human'], y=df1_UMI.loc[df1_UMI['label']==value,'mouse'],
#                c=my_color[counter], s=15)
#ax1.legend(labels=['Human','Mouse','Mixed'], loc='upper right').get_frame().set_facecolor('none')
sns.scatterplot(x='human',y='mouse',hue='label',hue_order=['Human','Mouse','Mixed'], data=df1_UMI,
                edgecolor='none', s=15, ax=ax1)
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles=handles[1:], labels=[f"{x} ({(df1_UMI['label']==x).sum()})" for x in labels[1:]], loc='upper right', frameon=False, handlelength=1, handleheight=1.125)
ax1.set_xlabel('Human transcripts')
ax1.set_ylabel('Mouse transcripts')
ax1.set_xticklabels(['{:,}'.format(int(x)) for x in ax1.get_xticks().tolist()])
ax1.set_yticklabels(['{:,}'.format(int(x)) for x in ax1.get_yticks().tolist()])
ax1.set_title('Live cells')
sns.despine(ax=ax1)

ax2 = fig.add_subplot(1,2,2)
#for (counter, value) in enumerate(['Human','Mouse','Mixed']):
#    ax2.scatter(x=df2_UMI.loc[df2_UMI['label']==value,'human'], y=df2_UMI.loc[df2_UMI['label']==value,'mouse'],
#                c=my_color[counter], s=15)
#ax2.legend(labels=['Human','Mouse','Mixed'], loc='upper right').get_frame().set_facecolor('none')
sns.scatterplot(x='human',y='mouse',hue='label',hue_order=['Human','Mouse','Mixed'], data=df2_UMI,
                edgecolor='none', s=15, ax=ax2)
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles=handles[1:], labels=[f"{x} ({(df2_UMI['label']==x).sum()})" for x in labels[1:]], loc='upper right', frameon=False, handlelength=1, handleheight=1.125)
ax2.set_xlabel('Human transcripts')
ax2.set_ylabel('Mouse transcripts')
ax2.set_xticklabels(['{:,}'.format(int(x)) for x in ax2.get_xticks().tolist()])
ax2.set_yticklabels(['{:,}'.format(int(x)) for x in ax2.get_yticks().tolist()])
ax2.set_title('Fixed cells')
sns.despine(ax=ax2)
fig.tight_layout()
fig.savefig(my_dir+'mix_species_plot.svg', bbox_inches='tight')

# =============================================================================
# # Plot number of UMIs and genes detected
# =============================================================================
# Filter human and mouse cells separately
df1_human_only = df1_human.loc[:,df1_UMI['label']=='Human']
df1_mouse_only = df1_mouse.loc[:,df1_UMI['label']=='Mouse']
df2_human_only = df2_human.loc[:,df2_UMI['label']=='Human']
df2_mouse_only = df2_mouse.loc[:,df2_UMI['label']=='Mouse']
# Make a list of data frames
df_list = [df1_human_only.copy(), df1_mouse_only.copy(), df2_human_only.copy(), df2_mouse_only.copy()]

# Calculate number of transcripts
df_UMI_list = []
for (counter, value) in enumerate(df_list):
    df_temp = pd.DataFrame(value.sum(axis=0), columns=['n_UMIs'])
    # Label species
    if (counter == 0) | (counter == 2):
        df_temp['Species'] = 'Human'
    else:
        df_temp['Species'] = 'Mouse'
    # Label ilve or fixed
    if counter <= 1:
        df_temp['LorF'] = 'Live'
    else:
        df_temp['LorF'] = 'Fixed'
    df_UMI_list.append(df_temp)
# Combine the number of UMIs into 1 data frame
df_UMI_all = pd.concat(df_UMI_list, ignore_index=False)

# Violin plot number of UMIs
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(1,1,1)
sns.violinplot(x='Species', y='n_UMIs', hue='LorF', data=df_UMI_all, ax=ax, saturation=1)
#sns.boxplot(x='Species', y='n_UMIs', hue='LorF', data=df_UMI_all, width=0.2, showfliers=False, ax=ax)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, loc='upper left', frameon=False)
ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
ax.set_ylabel('Number of transcripts')
ax.xaxis.label.set_visible(False)
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(my_dir+'/mix_species_num_UMIs.svg', bbox_inches='tight')

# Clean up the DGE matrices to calculate the number of genes
for (counter, value) in enumerate(df_list):
    # Keep genes detected in at least 5 cells
    df_list[counter] = value.loc[(value>0).sum(axis=1)>=5,:]

# Calculate number of genes
df_gene_list = []
for (counter, value) in enumerate(df_list):
    df_temp = pd.DataFrame((value>0).sum(axis=0), columns=['n_genes'])
    # Label species
    if (counter == 0) | (counter == 2):
        df_temp['Species'] = 'Human'
    else:
        df_temp['Species'] = 'Mouse'
    # Label ilve or fixed
    if counter <= 1:
        df_temp['LorF'] = 'Live'
    else:
        df_temp['LorF'] = 'Fixed'
    df_gene_list.append(df_temp)
# Combine the number of genes into 1 data frame
df_gene_all = pd.concat(df_gene_list, ignore_index=False)

# Violin plot number of genes
fig = plt.figure(figsize=(4,3.7))
ax = fig.add_subplot(1,1,1)
sns.violinplot(x='Species', y='n_genes', hue='LorF', data=df_gene_all, ax=ax, saturation=1)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, loc='upper left', frameon=False)
ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
ax.xaxis.label.set_visible(False)
#ax.set_xlabel(' ')
ax.set_ylabel('Number of genes')
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(my_dir+'mix_species_num_genes.svg', bbox_inches='tight')


# =============================================================================
# # Downsample UMIs, and calculate the number of genes detected
# =============================================================================
# Only keep cells with at least 1500 UMI
df1_human_only_1000 = df1_human_only.loc[:,df1_human_only.sum(axis=0)>=1500]
df1_mouse_only_1000 = df1_mouse_only.loc[:,df1_mouse_only.sum(axis=0)>=1500]
df2_human_only_1000 = df2_human_only.loc[:,df2_human_only.sum(axis=0)>=1500]
df2_mouse_only_1000 = df2_mouse_only.loc[:,df2_mouse_only.sum(axis=0)>=1500]

# Loop through each data frame and subsample 1000 UMIs
df_subsample_list = [] # list to store the data frame of number of genes
np.random.seed(1)
for counter, df_i in enumerate([df1_human_only_1000, df1_mouse_only_1000, df2_human_only_1000, df2_mouse_only_1000]):
    
    # label the species and live/fixed
    if (counter % 2) == 0:
        species_i = 'Human'
    else:
        species_i = 'Mouse'
        
    if counter <= 1:
        LorF_i = 'Live'
    else:
        LorF_i = 'Fixed'
    
    # Initialize a data frame to store the number of genes
    df_subsample_list.append(pd.DataFrame(data={'n_genes':0, 'Species':species_i, 'LorF':LorF_i}, index=df_i.columns))
    
    # Assign a number ID to each gene
    gene_num_ID = np.arange(df_i.shape[0])
    for col_i in df_i:
        my_temp = np.repeat(gene_num_ID, df_i[col_i].values)
        my_temp = np.random.choice(my_temp, size=1000, replace=False)
        df_subsample_list[counter].at[col_i,'n_genes'] = len(set(my_temp))
        
# Combined into 1 data frames
df_subsample_all = pd.concat(df_subsample_list)

# Plot number of genes
# Violin plot number of genes
fig = plt.figure(figsize=(4,3.7))
ax = fig.add_subplot(1,1,1)
sns.violinplot(x='Species', y='n_genes', hue='LorF', data=df_subsample_all, ax=ax, saturation=1)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=handles, labels=labels, loc='lower right', frameon=False)
#ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
ax.xaxis.label.set_visible(False)
#ax.set_xlabel(' ')
ax.set_ylabel('Number of genes')
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(my_dir+'mix_species_num_genes_downsample.svg', bbox_inches='tight')

# =============================================================================
# # Export the species-specific cells for analysis with R's monocle
# # Genes are not filtered yet
# =============================================================================
# df1_human_only.to_csv(my_dir+'/Van1_S1_human_only.txt.gz', sep='\t')
# df1_mouse_only.to_csv(my_dir+'/Van1_S1_mouse_only.txt.gz', sep='\t')
# df2_human_only.to_csv(my_dir+'/Van2_S2_human_only.txt.gz', sep='\t')
# df2_mouse_only.to_csv(my_dir+'/Van2_S2_mouse_only.txt.gz', sep='\t')


# =============================================================================
# # Compare the spearman coefficient
# # Genes are filtered
# =============================================================================
# =============================================================================
# # Check if there is any overlapping cell barcodes
# # set() = no overlapping
# print(set(df_list[0].columns).intersection(df_list[2].columns))
# print(set(df_list[1].columns).intersection(df_list[3].columns))
# 
# # Combine live and fixed cells into 1 dataframe by species
# # Inner join: only keep genes detected in both conditions
# df3_human = pd.concat([df_list[i] for i in (0,2)], axis=1, join='inner', sort=False)
# df3_human.fillna(value=0, inplace=True)
# df3_mouse = pd.concat([df_list[i] for i in (1,3)], axis=1, join='inner', sort=False)
# df3_mouse.fillna(value=0, inplace=True)
# 
# # Save cell barcode label: live or fixed
# temp_x = pd.DataFrame(df_list[0].columns, columns=['barcode'])
# temp_x['LorF'] = 'Live'
# temp_y = pd.DataFrame(df_list[2].columns, columns=['barcode'])
# temp_y['LorF'] = 'Fixed'
# df3_human_LorF = pd.concat([temp_x,temp_y], ignore_index=True)
# temp_x = pd.DataFrame(df_list[1].columns, columns=['barcode'])
# temp_x['LorF'] = 'Live'
# temp_y = pd.DataFrame(df_list[3].columns, columns=['barcode'])
# temp_y['LorF'] = 'Fixed'
# df3_mouse_LorF = pd.concat([temp_x,temp_y], ignore_index=True)
# 
# # Calculate pairwise spearman correlation between columns ***SLOWW***
# df3_human_corr = df3_human.corr(method='spearman')
# # Cluster
# network_lut = {'Live':sns.color_palette()[0], 'Fixed':sns.color_palette()[1]}
# network_colors = df3_human_LorF['LorF'].map(network_lut)
# network_colors = network_colors.values # convert Panda series to numpy array to avoid row/col color labels
# 
# df3_human_dist = 1 - df3_human_corr
# linkage = sch.linkage(sp.distance.squareform(df3_human_dist), method='ward') # need to convert to condensed distance matrix, because linkage requires condensed distance matrix
# 
# # Plot clustermap
# fig = plt.figure(figsize = (5,5))
# g = sns.clustermap(df3_human_corr, cmap='viridis', row_linkage=linkage, col_linkage=linkage,
#                    row_colors=network_colors, col_colors=network_colors,
#                    linewidths=0, cbar_kws={'label': 'Spearman\ncorrelation'})
# =============================================================================


# =============================================================================
# # Continue data analysis from monocle: mix_species_analysis.R
# =============================================================================
# Import normalized data
df_all_human_norm = pd.read_csv(my_dir+'all_human_norm_gene_expr.txt', sep='\t', index_col=0)
df_all_mouse_norm = pd.read_csv(my_dir+'all_mouse_norm_gene_expr.txt', sep='\t', index_col=0)

# Calculate means by grouping of live or fixed
# Transpose
df_all_human_norm_mean = df_all_human_norm.T
df_all_mouse_norm_mean = df_all_mouse_norm.T
# Add label
df_all_human_norm_mean = df_all_human_norm_mean.join(df_UMI_all['LorF'], how='left')
df_all_mouse_norm_mean = df_all_mouse_norm_mean.join(df_UMI_all['LorF'], how='left')
# Group mean
df_all_human_norm_mean = df_all_human_norm_mean.groupby(by='LorF', sort=False).mean()
df_all_mouse_norm_mean = df_all_mouse_norm_mean.groupby(by='LorF', sort=False).mean()
# Transpose
df_all_human_norm_mean = df_all_human_norm_mean.T
df_all_mouse_norm_mean = df_all_mouse_norm_mean.T

# Set the minimum expression level to 0.1
df_all_human_norm_mean[df_all_human_norm_mean <= 0.1] = 0.1
df_all_mouse_norm_mean[df_all_mouse_norm_mean <= 0.1] = 0.1

# Plot human
fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(1,2,1)
ax1.set_xscale('log')
ax1.set_yscale('log')
sns.scatterplot(x='Live', y='Fixed', data=df_all_human_norm_mean,
                color='grey', edgecolor='none', s=10, ax=ax1)
ax1.set_aspect('equal')
ax1.set_xlim(0.75e-1, 3e1)
ax1.set_ylim(0.75e-1, 3e1)
ax1.plot([0.75e-1, 3e1*0.9], [0.75e-1, 3e1*0.9], color='red')
ax1.set_xlabel('Normalized expression (live cells)')
ax1.set_ylabel('Normalized expression (fixed cells)')
#ax1.set_title('Human cells', fontdict={'fontsize':12})
# Show pearson correlation of log10(values)
ax1.text(0.07, 0.81, r'$\mathrm{\rho}$ = '+f'{np.log10(df_all_human_norm_mean).corr().iloc[0,1]:.3f}', transform=ax1.transAxes)
ax1.text(0.5, 1, 'Human cells', horizontalalignment='center', verticalalignment='top', transform=ax1.transAxes)
sns.despine(ax=ax1)

# Plot mouse
ax2 = fig.add_subplot(1,2,2)
ax2.set_xscale('log')
ax2.set_yscale('log')
sns.scatterplot(x='Live', y='Fixed', data=df_all_mouse_norm_mean,
                color='grey', edgecolor='none', s=10, ax=ax2)
ax2.set_aspect('equal')
ax2.set_xlim(0.75e-1, 4.5e1)
ax2.set_ylim(0.75e-1, 4.5e1)
ax2.plot([0.75e-1, 4.5e1*0.9], [0.75e-1, 4.5e1*0.9], color='red')
ax2.set_xlabel('Normalized expression (live cells)')
ax2.set_ylabel('Normalized expression (fixed cells)')
#ax2.set_title('Mouse cells', fontdict={'fontsize':12})
# Show pearson correlation of log10(values)
ax2.text(0.07, 0.81, r'$\mathrm{\rho}$ = '+f'{np.log10(df_all_mouse_norm_mean).corr().iloc[0,1]:.3f}', transform=ax2.transAxes)
ax2.text(0.5, 1, 'Mouse cells', horizontalalignment='center', verticalalignment='top', transform=ax2.transAxes)
sns.despine(ax=ax2)
fig.tight_layout()
#fig.savefig(my_dir+'mix_species_livevsfixed_expression.svg', bbox_inches='tight') # vector graphic file is too large
fig.savefig(my_dir+'mix_species_livevsfixed_expression.png', bbox_inches='tight', dpi=600, pad_inches=0.01)

# =============================================================================
# # Compare UMI vs read counts
# =============================================================================
# Import read count data
df1_readcounts = pd.read_csv(my_dir+'STAR_AB-HE190129-Van1_S1_knee_barcodes_readcounts.txt', sep=',', index_col=1, header=None, names=['read_counts'])
df2_readcounts = pd.read_csv(my_dir+'STAR_AB-HE190129-Van2_S2_knee_barcodes_readcounts.txt', sep=',', index_col=1, header=None, names=['read_counts'])

# Combine the two readcounts df
df_readcounts_all = pd.concat([df1_readcounts, df2_readcounts], ignore_index=False)

# Add the read counts data to df_UMI_all
df_UMI_readcounts_all = df_UMI_all.join(df_readcounts_all, how='left')

# Plot
#g = sns.FacetGrid(data=df_UMI_readcounts_all, col='Species', sharey=False, sharex=False, height=4)
#g.map(sns.scatterplot, 'read_counts', 'n_UMIs', 'LorF', s=15, edgecolor='none')
#g.set_xlabels('Read counts')
#g.set_ylabels('Number of transcripts')
## Add thousand separator
#for counter, value in enumerate(g.axes[0]):
#    g.axes[0][counter].set_xticklabels(['{:,}'.format(int(x)) for x in value.get_xticks().tolist()])
#    g.axes[0][counter].set_yticklabels(['{:,}'.format(int(x)) for x in value.get_yticks().tolist()])
    
# Plot UMIs vs read counts
g = sns.lmplot(x='read_counts', y='n_UMIs', hue='LorF', hue_order=['Live','Fixed'], data=df_UMI_readcounts_all, col='Species', sharey=False, sharex=False, height=4,
               ci=None, scatter_kws={'s':15,'alpha':0.5,'edgecolors':'None'})
g.set_xlabels('Read counts')
g.set_ylabels('Number of transcripts')
g.set_titles(template='{col_name}')
for counter, value in enumerate(g.axes[0]):
    # Add thousand separator
#    g.axes[0][counter].set_xticklabels(['{:,}'.format(int(x)) for x in value.get_xticks().tolist()])
    g.axes[0][counter].set_xticklabels([f'{int(x):,}' for x in value.get_xticks().tolist()])
    g.axes[0][counter].set_yticklabels([f'{int(x):,}' for x in value.get_yticks().tolist()])
    
    # Change axis limit
    xmin, xmax = value.get_xlim()
    g.axes[0][counter].set_xlim(-xmax/50, xmax)
    ymin, ymax = value.get_ylim()
    g.axes[0][counter].set_ylim(-ymax/50, ymax)
    
    # Set x label
#    g.axes[0][counter].set_ylabel('Number of transcripts')
g.savefig(my_dir+'UMI_vs_readcounts.svg', bbox_inches='tight')

# Violin plot number of UMIs pere 1,000 reads
#df_UMI_readcounts_all = df_UMI_readcounts_all.assign(UMI_per_readcount=df_UMI_readcounts_all['n_UMIs'].values/df_UMI_readcounts_all['read_counts'].values*1000)
#g = sns.FacetGrid(data=df_UMI_readcounts_all, col='Species', col_order=['Human','Mouse'], sharey=False, sharex=False, height=4)
#g.map(sns.violinplot, 'LorF', 'UMI_per_readcount', order=['Live','Fixed'], s=15, edgecolor='None', alpha=0.5)
#    
#g.set_xlabels('')
#g.set_ylabels('Number of transcripts\nper 1,000 reads')
#g.set_titles(template='{col_name}')
#g.axes[0][0].set_ylim(450,850)
#g.axes[0][1].set_ylim(400,900)
#g.axes[0][-1].legend({key+' cells': value for value, key in zip(*g.axes[0][1].get_legend_handles_labels()) if key!='LorF'}, loc='lower right', frameon=False, handlelength=1, handleheight=1.125)
#g.savefig(my_dir+'num_UMIs_per_readcounts.svg', bbox_inches='tight')

df_UMI_readcounts_all = df_UMI_readcounts_all.assign(UMI_per_readcount=df_UMI_readcounts_all['n_UMIs'].values/df_UMI_readcounts_all['read_counts'].values*1000)
np.random.seed(5) # jitter random seed
fig = plt.figure(figsize=(6,3.5))
ax1 = fig.add_subplot(1,2,1)
sns.violinplot(x='LorF', y='UMI_per_readcount', order=['Live','Fixed'],
               data=df_UMI_readcounts_all.loc[df_UMI_readcounts_all['Species']=='Human',:],
               color=sns.color_palette('Dark2',2)[0], inner=None, ax=ax1)
sns.stripplot(x='LorF', y='UMI_per_readcount', order=['Live','Fixed'],
              data=df_UMI_readcounts_all.loc[df_UMI_readcounts_all['Species']=='Human',:],
              jitter=0.3, color='k', size=3, ax=ax1)
ax1.set_xlabel('')
ax1.set_ylabel('Number of transcripts\nper 1,000 reads')
ax1.set_ylim(400,900)
ax1.set_title('Human cells')
ax2 = fig.add_subplot(1,2,2)
sns.violinplot(x='LorF', y='UMI_per_readcount', order=['Live','Fixed'],
               data=df_UMI_readcounts_all.loc[df_UMI_readcounts_all['Species']=='Mouse',:],
               color=sns.color_palette('Dark2',2)[1], inner=None, ax=ax2)
sns.stripplot(x='LorF', y='UMI_per_readcount', order=['Live','Fixed'],
              data=df_UMI_readcounts_all.loc[df_UMI_readcounts_all['Species']=='Mouse',:],
              jitter=0.3, color='k', size=3, ax=ax2)
ax2.set_xlabel('')
ax2.set_ylabel('Number of transcripts\nper 1,000 reads')
ax2.set_ylim(400,900)
ax2.set_title('Mouse cells')
sns.despine(fig=fig)
fig.tight_layout()
# fig.savefig(my_dir+'num_UMI_per_readcounts.svg', bbox_inches='tight')
fig.savefig(my_dir+'num_UMI_per_readcounts.png', dpi=600, bbox_inches='tight')

# =============================================================================
# # Compare genes vs read counts
# =============================================================================
# Add the read counts data to df_gene_all
df_gene_readcounts_all = df_gene_all.join(df_readcounts_all, how='left')

# Fit a saturation curve to each group: y = ax/(b+x)
def my_sat_fun(x, a, b):
    return a*x/(b+x) # use p0=[1e4,1]
#    return a*(1 - np.exp(-x/b)) # use p0=[1e4,4e4]

# Plot # genes vs read counts
g = sns.FacetGrid(data=df_gene_readcounts_all, col='Species', col_order=['Human','Mouse'], sharey=False, sharex=False, height=4)
g.map(sns.scatterplot, 'read_counts', 'n_genes', 'LorF', s=15, edgecolor='None', alpha=0.5)
# Curve fitting
for counter1, value1 in enumerate(g.axes[0]):
    xmin, xmax = value1.get_xlim()
    for value2 in ['Live','Fixed']:
        xdata = df_gene_readcounts_all.loc[(df_gene_readcounts_all['Species']==['Human','Mouse'][counter1]) & (df_gene_readcounts_all['LorF']==value2),'read_counts'].values
        ydata = df_gene_readcounts_all.loc[(df_gene_readcounts_all['Species']==['Human','Mouse'][counter1]) & (df_gene_readcounts_all['LorF']==value2),'n_genes'].values
        popt, _ = curve_fit(my_sat_fun, xdata=xdata, ydata=ydata, p0=[1e4,1])
        
        
        g.axes[0][counter1].plot(np.arange(0,xmax), my_sat_fun(np.arange(0,xmax),popt[0],popt[1]))
        
    # Add thousand separator
    g.axes[0][counter1].set_xticklabels(['{:,}'.format(int(x)) for x in value1.get_xticks().tolist()])
    g.axes[0][counter1].set_yticklabels(['{:,}'.format(int(x)) for x in value1.get_yticks().tolist()])
    
g.set_xlabels('Read counts')
g.set_ylabels('Number of genes')
g.set_titles(template='{col_name}')
g.axes[0][-1].legend({key+' cells': value for value, key in zip(*g.axes[0][1].get_legend_handles_labels()) if key!='LorF'}, loc='lower right', frameon=False, handlelength=1, handleheight=1.125)

# Violin plot number of genes pere 1,000 reads
df_gene_readcounts_all = df_gene_readcounts_all.assign(gene_per_readcount=df_gene_readcounts_all['n_genes'].values/df_gene_readcounts_all['read_counts'].values*1000)
g = sns.FacetGrid(data=df_gene_readcounts_all, col='Species', col_order=['Human','Mouse'], sharey=False, sharex=False, height=4)
g.map(sns.violinplot, 'LorF', 'gene_per_readcount', order=['Live','Fixed'], s=15, edgecolor='None', alpha=0.5)
    
g.set_xlabels('')
g.set_ylabels('Number of genes\nper 1,000 reads')
g.set_titles(template='{col_name}')
g.axes[0][0].set_ylim(100,620)
g.axes[0][-1].legend({key+' cells': value for value, key in zip(*g.axes[0][1].get_legend_handles_labels()) if key!='LorF'}, loc='lower right', frameon=False, handlelength=1, handleheight=1.125)
g.savefig(my_dir+'num_genes_per_readcounts.svg', bbox_inches='tight')

plt.show()