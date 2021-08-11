# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:14:20 2019

@author: Van Local

Analyze velocyto loom files
"""

import velocyto as vcy
import numpy as np
import loompy
from sklearn.manifold import TSNE

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

# Combine loom objects
#loompy.combine(["/project2/tays/Van/20180927/Alignments/negative_velocyto/negative_velocyto.loom","/project2/tays/Van/20180927/Alignments/positive_velocyto/positive_velocyto.loom"],
#               "/project2/tays/Van/20180927/Alignments/combined_velocyto.loom", key="Accession")

# Plot fractions of spliced/unspliced
vlm = vcy.VelocytoLoom("/project2/tays/Van/20180927/Alignments/combined_velocyto.loom")
# Manually add SampleID
sampleID = np.repeat("K8.1-", vlm.ca["CellID"].shape[0]).astype(object)
sampleID[["positive" in s for s in vlm.ca["CellID"]]] = "K8.1+"
vlm.ca["SampleID"] = sampleID
vlm.plot_fractions()

# Create loom object
vlm_neg = vcy.VelocytoLoom("/project2/tays/Van/20180927/Alignments/negative_velocyto/negative_velocyto.loom")
vlm_pos = vcy.VelocytoLoom("/project2/tays/Van/20180927/Alignments/positive_velocyto/positive_velocyto.loom")

# Annotate sample
vlm_neg.ca["SampleID"] = np.repeat("K8.1-", len(vlm_neg.ca["CellID"])).astype(object)
vlm_pos.ca["SampleID"] = np.repeat("K8.1+", len(vlm_pos.ca["CellID"])).astype(object)

# Set cluster
vlm_neg.set_clusters(vlm_neg.ca["SampleID"])
vlm_pos.set_clusters(vlm_pos.ca["SampleID"])

# Plot spliced vs unspliced reads
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(1,2,1)
ax.scatter(vlm_neg.initial_cell_size, vlm_neg.initial_Ucell_size, alpha=0.5, s=5, c=sns.color_palette()[1])
ax.set_xlabel("Spliced molecules")
ax.set_ylabel("Unspliced molecules")
ax.set_title("K8.1-")
sns.despine(ax=ax)
ax = fig.add_subplot(1,2,2)
ax.scatter(vlm_pos.initial_cell_size, vlm_pos.initial_Ucell_size, alpha=0.5, s=5, c=sns.color_palette()[0])
ax.set_xlabel("Spliced molecules")
ax.set_ylabel("Unspliced molecules")
ax.set_title("K8.1+")
sns.despine(ax=ax)
fig.tight_layout()
fig.savefig("/project2/tays/Van/20180927/velocyto_python_analysis/spliced_unspliced_counts.svg", bbox_inches='tight')

# =============================================================================
# # Filter cells
# vlm_neg.filter_cells(bool_array=vlm_neg.initial_Ucell_size > np.percentile(vlm_neg.initial_Ucell_size, 0.5))
# vlm_pos.filter_cells(bool_array=(vlm_pos.initial_Ucell_size > np.percentile(vlm_pos.initial_Ucell_size, 0.5)) &
#                      (vlm_pos.initial_cell_size > 0) & (vlm_pos.initial_Ucell_size > 0))
# 
# # Filter genes
# vlm_neg.score_detection_levels(min_expr_counts=40, min_cells_express=30)
# vlm_neg.filter_genes(by_detection_levels=True)
# vlm_pos.score_detection_levels(min_expr_counts=40, min_cells_express=30)
# vlm_pos.filter_genes(by_detection_levels=True)
# 
# 
# # Feature selection
# vlm_neg.score_cv_vs_mean(1000, plot=True, max_expr_avg=50)
# vlm_neg.filter_genes(by_cv_vs_mean=True)
# vlm_pos.score_cv_vs_mean(2000, plot=True, max_expr_avg=50)
# vlm_pos.filter_genes(by_cv_vs_mean=True)
# 
# # Normalize by size
# vlm_neg._normalize_S(relative_size=vlm_neg.initial_cell_size, target_size=np.mean(vlm_neg.initial_cell_size))
# vlm_neg._normalize_U(relative_size=vlm_neg.initial_Ucell_size, target_size=np.mean(vlm_neg.initial_Ucell_size))
# vlm_pos._normalize_S(relative_size=vlm_pos.initial_cell_size, target_size=np.mean(vlm_pos.initial_cell_size))
# vlm_pos._normalize_U(relative_size=vlm_pos.initial_Ucell_size, target_size=np.mean(vlm_pos.initial_Ucell_size))
# 
# # PCA
# vlm_neg.perform_PCA()
# vlm_pos.perform_PCA()
# fig, ax = plt.subplots(1,1, figsize=(4,4))
# ax.plot(np.cumsum(vlm_pos.pca.explained_variance_ratio_[:100]))
# 
# # Gamma fit
# k = 50
# vlm_neg.knn_imputation(n_pca_dims=20, k=k, balanced=True, b_sight=k*4, b_maxl=k*3, n_jobs=16)
# k = 200
# vlm_pos.knn_imputation(n_pca_dims=20, k=k, balanced=True, b_sight=k*4, b_maxl=k*3, n_jobs=16)
# #vlm_pos.knn_imputation(n_pca_dims=20, k=k, balanced=False, n_jobs=16)
# vlm_neg.fit_gammas()
# vlm_pos.fit_gammas()
# #fig = plt.figure(figsize=(8,4))
# #ax = fig.add_subplot(1,2,1)
# #vlm_neg.plot_phase_portraits(["ISCU"])
# #ax.set_title("K8.1-")
# fig, ax = plt.subplots(1,1, figsize=(4,4))
# vlm_pos.plot_phase_portraits(["ISCU"])
# ax.set_title("K8.1+")
# 
# 
# 
# # Calculate velocity
# vlm_pos.predict_U()
# vlm_pos.calculate_velocity()
# vlm_pos.calculate_shift(assumption="constant_velocity")
# vlm_pos.extrapolate_cell_at_t(delta_t=1.)
# 
# # t-SNE
# vlm_pos.ts = TSNE(random_state=2019).fit_transform(vlm_pos.pcs[:,:20])
# vlm_pos.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1,
#                                  n_neighbors=50, knn_random=True, sampled_fraction=0.9)
# vlm_pos.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)
# vlm_pos.calculate_grid_arrows(smooth=0.5, steps=(40, 40), n_neighbors=15)
# 
# # Plot velocity arrows
# fig = plt.figure(figsize=(4,4))
# vlm_pos.plot_grid_arrows(quiver_scale=0.01,
#                          scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True}, min_mass=1, angles='xy', scale_units='xy',
#                          headaxislength=2.75, headlength=5, headwidth=4.8, minlength=1,
#                          plot_random=True, scale_type="absolute")
# fig, ax = plt.subplots(1,1,figsize=(4,4))
# ax.scatter(vlm_pos.flow_embedding[:,0], vlm_pos.flow_embedding[:,1], **{"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True})
# ax.quiver(vlm_pos.embedding[:,0],vlm_pos.embedding[:,1], vlm_pos.delta_embedding[:,0], vlm_pos.delta_embedding[:,1])
# ax.set_title("Embedding")
# fig, ax = plt.subplots(1,1,figsize=(4,4))
# ax.scatter(vlm_pos.flow_embedding[:,0], vlm_pos.flow_embedding[:,1], **{"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True})
# ax.quiver(vlm_pos.flow_grid[:,0],vlm_pos.flow_grid[:,1], vlm_pos.flow[:,0], vlm_pos.flow[:,1], scale=0.2, angles='xy', scale_units='xy')
# ax.set_title("Grid")
# =============================================================================

plt.show()