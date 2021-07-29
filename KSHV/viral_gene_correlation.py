# -*- coding: utf-8 -*-
"""
This script is for analyzing scaled expression levels of viral genes

@author: Van Local
"""
# Import libraries
import numpy as np
import pandas as pd

import scipy.spatial as sp
import scipy.cluster.hierarchy as sch
import scipy.stats as stats

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import seaborn as sns

# from sklearn.preprocessing import StandardScaler
# from sklearn.decomposition import PCA
# from sklearn.cross_decomposition import CCA # canonical correlation analysis

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

my_dir = 'Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/DGE files/'

# Import data
df = pd.read_csv(my_dir+'scaled_expr_pos_kshv.csv', index_col=0) # scaled epxression of viral genes of positive cells
df_pct_kshv = pd.read_csv(my_dir+'pct_viral_pos_kshv.csv', index_col=0) # % of viral transcripts

# Transpose such that columns are genes/features
df = df.T

# Only keep genes detected in >=50% of cells
df_filter = df.copy()
df_filter = df_filter.loc[:,(df_filter>0).sum(axis=0)/df_filter.shape[0]>=0.9]


# Label viral genes by stage (ref: https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1003847 )
kshv_latent = ['KSHV:' + s for s in ['K1','K2','K12','ORF71','ORF72','ORF73']]
kshv_imme = ['KSHV:'+ s for s in['ORF11','ORF70','K4','K4.1','K4.2','K5','K6','ORF16','ORF45',
                                 'ORF50','K8','ORF57']] # immediate-early
kshv_early = ['KSHV:' + s for s in ['ORF4','ORF6','K3','ORF17.5','ORF18','ORF34','ORF35','ORF36',
                                    'ORF37','ORF38','ORF39','ORF46','ORF47','ORF58','ORF59',
                                    'ORF60','ORF61','K14','ORF74']]
kshv_late = ['KSHV:' + s for s in ['ORF8','ORF9','ORF10','K3A','ORF17','ORF21','ORF22','ORF23',
                                   'ORF24','ORF25','ORF26','ORF27','ORF28','ORF29','ORF30',
                                   'ORF31','ORF32','ORF33','ORF40','ORF41','ORF42','ORF43',
                                   'ORF44','ORF45.1','K8.1','ORF52','ORF53','ORF54','ORF55',
                                   'ORF56','K9','K10','K10.5','K11','ORF62','ORF65','ORF66',
                                   'ORF67','ORF67.5','ORF68','ORF69','ORF75']]
# Combine into 1 dictionary
kshv_stage_dict = {'Latent':kshv_latent, 'Immediate-early':kshv_imme, 'Early':kshv_early, 'Late':kshv_late}

#---------------------------------#
# Do PCA
#---------------------------------#
# =============================================================================
# # Standardize each gene and do pca
# df_scaled = StandardScaler().fit_transform(np.log10(df))
# my_pca = PCA(n_components=10, random_state=2019)
# df_pca = my_pca.fit_transform(df_scaled)
# 
# # Plot PCA
# fig, ax = plt.subplots(figsize=(4,4))
# ax.scatter(df_pca[:,0], df_pca[:,1], c=df_pct_kshv.iloc[:,0], cmap='viridis')
# =============================================================================

#---------------------------------#
# Plot correlation between genes
#---------------------------------#

# Calculate pairwise correlation between viral genes
df_spearman = df.corr(method='spearman')

# Extract the unique pairwise correlation
pairwise_corr = df_spearman.to_numpy()[np.triu_indices(df_spearman.shape[0], 1)]

# # Calculate correlation between viral genes with % viral transcripts
# pct_kshv_spearman = np.array([df_pct_kshv.iloc[:,0].corr(df.loc[:,s], method='spearman') for s in df])

# Calculate pairwise correlation between viral genes of the same stage
df_spearman_bystage = {}
for key, value in kshv_stage_dict.items():
    df_temp = df.loc[:,[i in value for i in df.columns]].corr(method='spearman')
    df_spearman_bystage[key] = df_temp.to_numpy()[np.triu_indices(df_temp.shape[0],1)]

# Calculate correlation between K8.1 and late viral genes
k81_late_spearman = np.array([stats.spearmanr(df_filter['KSHV:K8.1'].values, df_filter[s])[0] for s in df_filter.columns if s in kshv_late])

# Plot histogram of spearman by stage
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(6.2,6.2), tight_layout=True)
for counter, i in enumerate(['Latent','Immediate-early','Early','Late']):
    irow = counter // 2
    icol = counter % 2
    ax[irow,icol].hist(df_spearman_bystage[i], bins=20, color='grey', edgecolor='k')
    ax[irow,icol].set_xlabel("Spearman correlation coefficient")
    ax[irow,icol].set_ylabel("Number of gene pairs")
    ax[irow,icol].set_title(f"{i} viral genes")
    if i == 'Latent':
        ax[irow,icol].set_xlim(-0.1,0.2)
    elif i == 'Immediate-early':
        ax[irow,icol].set_ylim(0,12)
    elif i == 'Early':
        ax[irow,icol].set_ylim(0,40)
    else:
        ax[irow,icol].set_xlim(-0.25,0.75)
sns.despine(fig=fig)
fig.savefig("Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/hist_viral_correlation_bystage.svg", bbox_inches='tight')

# Bin cells by percentage of viral transcripts, then calculate spearman coefficients
temp_bins = pd.cut(df_pct_kshv.iloc[:,0], bins=5, retbins=True)
temp_bins_corr = {}
for i in temp_bins[0].cat.categories:
    temp_bins_df = df.loc[temp_bins[0].index[temp_bins[0]==i],:]
    temp_df_corr = temp_bins_df.corr(method='spearman')
    temp_bins_corr[i] = temp_df_corr.to_numpy()[np.triu_indices(temp_df_corr.shape[0],1)]
temp_bins_corr = pd.DataFrame(temp_bins_corr)
temp_bins_corr2 = pd.melt(temp_bins_corr, var_name='bin', value_name='corr')

# Plot box plot
fig, ax = plt.subplots(figsize=(4.2,3.2), tight_layout=True)
sns.boxplot(y='bin', x='corr', data=temp_bins_corr2, ax=ax, fliersize=3)
for box in ax.artists:
    box.set_edgecolor('black')
    box.set_facecolor('white')
labels = [f"{s.left*100:.2f}-{s.right*100:.2f}%" for s in temp_bins[0].cat.categories]
ax.set_yticklabels(labels)
ax.set_ylabel("Percentage of viral transcripts")
ax.set_xlabel("Spearman correlation coefficient")
sns.despine(fig=fig)
fig.savefig("Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/boxplot_viral_correlation_binned.svg", bbox_inches='tight')

# # Plot histogram of spearman correlation among late viral genes, depending on how highly expressed the genes are
# df_late = df.loc[:,[s in kshv_late for s in df.columns]]
# fig, ax = plt.subplots(figsize=(4,4))
# for i, j in zip([0.1,0.5,0.9], ['black','blue','red']):
#     df_late_temp = df_late.loc[:,(df_late>0).sum(axis=0)/df_late.shape[0]>=i]
#     # Pairwise correlation
#     df_late_temp = df_late_temp.corr(method='spearman')
#     corr_temp = np.triu(df_late_temp.values,1).flatten()
#     corr_temp = corr_temp[corr_temp!=0]
#     sns.kdeplot(corr_temp, color=j, ax=ax)
# ax.set_xlabel('Spearman correlation')
# ax.set_ylabel('Density')
# #ax.set_xlim(-0.22,1)
# sns.despine(ax=ax)
# # Manually add legend
# black_line = mlines.Line2D([], [], color='black', label=r'$\geq$10% of cells')
# blue_line = mlines.Line2D([], [], color='blue', label=r'$\geq$50% of cells')
# red_line = mlines.Line2D([], [], color='red', label=r'$\geq$90% of cells')
# ax.legend(handles=[black_line,blue_line,red_line], loc='best', edgecolor='black')
# fig.tight_layout()
# fig.savefig('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/late_viral_genes_correlation.svg', bbox_inches='tight')

# =============================================================================
# # Calculate correlation between K8.1 and all other viral genes
# df_k81 = df.loc[:, 'KSHV:K8.1'].values
# df_non_k81 = df.loc[:, df.columns != 'KSHV:K8.1'].values
# k81_spearman = np.array([stats.spearmanr(df_k81,s, nan_policy='omit')[0] for s in df_non_k81.T])
# 
# 
# # Plot histogram
# fig, ax = plt.subplots(figsize=(5,4))
# #sns.distplot(pairwise_corr, bins=20, kde=False, label='All pair-wise correlations', norm_hist=True, color='blue', hist_kws={'edgecolor':'blue'})
# #sns.distplot(k81_spearman, bins=20, kde=False, label='Correlations with K8.1 gene', norm_hist=True, color='gray', hist_kws={'edgecolor':'k'})
# sns.kdeplot(pairwise_corr, color='blue', ax=ax)
# #sns.kdeplot(pct_kshv_spearman, color='red', ax=ax)
# sns.kdeplot(k81_spearman, color='k', ax=ax)
# ax.set_xlabel('Spearman correlation')
# ax.set_ylabel('Density')
# sns.despine(ax=ax)
# blue_line = mlines.Line2D([], [], color='blue', label='All pair-wise correlations')
# red_line = mlines.Line2D([], [], color='k', label='Correlations with K8.1')
# ax.legend(handles = [blue_line, red_line], loc=(0.6,0.8), edgecolor='k')
# fig.tight_layout()
# fig.savefig('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/viral_genes_correlation.svg', bbox_inches='tight')
# =============================================================================

#---------------------------------#
# Scatter plot of viral genes
#---------------------------------#
df_scatter = df.copy()
df_scatter[df_scatter < 0.1] = 0.1
# Plot two highly expressed early viral genes
fig, ax = plt.subplots(figsize=(3.5,3.5), tight_layout=True)
ax.scatter(df['KSHV:ORF58'].values, df['KSHV:K7'].values,
           c='grey', alpha=0.8, edgecolor="None", s=10)
ax.text(0.4, 0.9, r'$\mathrm{\rho}_{\mathrm{Spearman}}$ = '+f'{df[["KSHV:ORF58","KSHV:K7"]].corr(method="spearman").iloc[0,1]:.3f}', transform=ax.transAxes)
ax.set_xlabel('Expression of ORF58')
ax.set_ylabel("Expression of K7")
sns.despine(ax=ax)
fig.savefig('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/k7_vs_orf58.png',
            dpi=600, bbox_inches='tight')


fig, ax = plt.subplots(figsize=(3.5,3.5), tight_layout=True)
ax.scatter(df['KSHV:K8.1'].values, df['KSHV:ORF52'].values,
           c='grey', alpha=0.8, edgecolor="None", s=10)
ax.text(0.4, 0.05, r'$\mathrm{\rho}_{\mathrm{Spearman}}$ = '+f'{df[["KSHV:K8.1","KSHV:ORF52"]].corr(method="spearman").iloc[0,1]:.3f}', transform=ax.transAxes)
ax.set_xlabel('Expression of K8.1')
ax.set_ylabel("Expression of ORF52")
sns.despine(ax=ax)
fig.savefig('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/orf52_vs_k8.1.png',
            dpi=600, bbox_inches='tight')

# All pair-wise correlation of viral genes
fig, ax = plt.subplots(figsize=(6,6))
sns.heatmap(df.corr(method="spearman"), ax=ax)

#---------------------------------#
# Cluster viral genes
#---------------------------------#
# Annotate the stage of viral genes
kshv_stage_each_gene_dict = {s:key for key, value in kshv_stage_dict.items() for s in value}
kshv_stage_all_annotated = []
for value in df.columns:
    if value in kshv_stage_each_gene_dict.keys():
        kshv_stage_all_annotated.append(kshv_stage_each_gene_dict[value])
    else:
        kshv_stage_all_annotated.append('Unannotated')  

my_colors = sns.color_palette('Set1',5)
my_colors_stage_dict = {'Latent':0, 'Immediate-early':1, 'Early':2, 'Late':3, 'Unannotated':4}
network_colors = [my_colors[my_colors_stage_dict[s]] for s in kshv_stage_all_annotated]

# Hierarchical clustering using 1 - spearman correlation
df_dist = 1 - df_spearman
        
linkage = sch.linkage(sp.distance.squareform(df_dist), method='complete') # need to convert to condensed distance matrix, because linkage requires condensed distance matrix
g = sns.clustermap(df_spearman, cmap='RdYlBu_r', row_linkage=linkage, col_linkage=linkage, row_colors=network_colors, col_colors=network_colors,
                   vmin=-1, vmax=1, xticklabels=False, yticklabels=False, figsize=(10,10))
# g.savefig('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/clustermap_viral_genes.svg', bbox_inches='tight')
g.savefig('Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/clustermap_viral_genes.png', dpi=600, bbox_inches='tight')

fig, ax = plt.subplots(figsize=(3.2,3.1), tight_layout=True)
ax.hist(df_spearman.to_numpy()[np.triu_indices(df_spearman.shape[0],1)],
        bins=20, color='grey', edgecolor='k')
ax.set_xlabel("Pairwise Spearman coefficient")
ax.set_ylabel("Number of gene pairs")
ax.set_xticks([-0.5, 0, 0.5, 1])
ax.set_ylim(0,1200)
sns.despine(fig=fig)
fig.savefig("Z:/Van/20180927 - Drop-seq fixed BC3 reactivated/figures/viral_genes_spearman_hist.svg", bbox_inches='tight')

# =============================================================================
# # Hierarchical clustering using 1 - spearman correlation
# df_dist2 = sp.distance.squareform(sp.distance.pdist(df_scaled.T, metric='euclidean'))
# linkage = sch.linkage(sp.distance.squareform(df_dist2), method='complete')
# fig = plt.figure(figsize = (5,5))
# g = sns.clustermap(df_dist2, cmap='vlag', row_linkage=linkage, col_linkage=linkage, row_colors=network_colors, col_colors=network_colors)
# =============================================================================

#---------------------------------#
# Import host + viral genes
#---------------------------------#
#df_all = pd.read_csv(my_dir+'scaled_expr_pos.csv', index_col=0)
#
## Keep cells with more than 20% viral transcripts, and remove genes expressed in less than 5% of the cells
#df_all = df_all.loc[(df_all>0).sum(axis=1)>=0.05*df_all.shape[1],df_pct_kshv.index]