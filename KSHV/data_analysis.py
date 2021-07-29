# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 15:15:53 2018

@author: Van Local
"""

# Import packages here
import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import scipy.stats as stats

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

my_dir = 'Z:/Van/20200113 - Michiel qPCR data/'

# =============================================================================
# # 570 experiment: overexpression
# =============================================================================
# nls-GFP fold induction data (sample 115)
# Columns are replicates, rows are 24, 48, 72, 96, 120 hours
df_115 = np.array([ [0.78, 1.04, 1.24],
                    [0.41, 0.28, 0.33],
                    [3.49, 4.84, 5.51],
                    [6.60, 6.13, 6.31],
                    [6.72, 6.97, 7.93] ])

# TMEM119 over expression (sample 117)
df_117 = np.array([ [0.96, 1.24, 0.84],
                    [0.25, 0.24, 0.33],
                    [5.27, 4.66, 3.82],
                    [11.18, 11.43, 8.20],
                    [9.39, 14.0, 13.93] ])

# Time points
xt = np.array([24,48,72,96,120])

# Plot error bars
fig,ax = plt.subplots(figsize=(12,3.1), ncols=4)
ax[0].errorbar(x=xt,y=np.mean(df_115, axis=1), yerr=np.std(df_115, axis=1, ddof=1)/np.sqrt(3),
               marker="o", ms=3, capsize=3, lw=2, capthick=2, label="nls-GFP")
ax[0].errorbar(x=xt,y=np.mean(df_117, axis=1), yerr=np.std(df_117, axis=1, ddof=1)/np.sqrt(3),
               marker="o", ms=3, capsize=3, lw=2, capthick=2, label="TMEM119")
ax[0].set_xticks(xt)
ax[0].set_xlabel("Hours post transfection")
ax[0].set_ylabel("ORF74 relative expression")
ax[0].set_ylim(0,15)
ax[0].legend(edgecolor="black")

# one-sided welch t-test
ttest1 = stats.ttest_ind(df_117[-1,:], df_115[-1,:], axis=None, equal_var=False)
print(ttest1)
ttest1[1]/2

#fig.tight_layout()
#fig.savefig(my_dir+'overexpression.svg', bbox_inches='tight')
#fig.savefig(my_dir+'overexpression.png', bbox_inches='tight', dpi=300)
#plt.show()

# =============================================================================
# # 583 experiment: siRNA
# =============================================================================
# ORF50 ===============
# si control
# Columns are replicates, rows are 0, 12, 24, 36 and 48 hours
df_ctrl = np.array([ [1.03, 1.07, 0.9],
                     [6.71, 7.85, 9.39],
                     [30.03, 33.81, 34.86],
                     [38.38, 32.95, 35.32],
                     [20.67, 19.10, 19.89] ])
# si TMEM119
df_tmem = np.array([ [0.72, 0.76, 0.73],
                     [4.34, 4.02, 4.22],
                     [21.79, 20.68, 21.16],
                     [25.39, 23.53, 23.54],
                     [12.31, 10.57, 13.41] ])

# Time points
xt = np.array([0, 12, 24, 36, 48])

ax[1].errorbar(x=xt,y=np.mean(df_ctrl, axis=1), yerr=np.std(df_ctrl, axis=1, ddof=1)/np.sqrt(3), marker="o", ms=3, capsize=3, lw=2, capthick=2, label="si-ctrl")
ax[1].errorbar(x=xt,y=np.mean(df_tmem, axis=1), yerr=np.std(df_tmem, axis=1, ddof=1)/np.sqrt(3), marker="o", ms=3, capsize=3, lw=2, capthick=2, label="si-TMEM119")
ax[1].set_xticks(xt)
ax[1].set_xlabel("Hours post stimulation")
ax[1].set_ylabel("ORF50 relative expression")
ax[1].set_ylim(0,40)
ax[1].legend(edgecolor="black")

# one-sided welch t-test
ttest2 = stats.ttest_ind(df_ctrl[-1,:], df_tmem[-1,:], axis=None, equal_var=False)
print(ttest2)
ttest2[1]/2

# ORF57 ===============
# si control
df_ctrl = np.array([ [1.02, 1.08, 0.91],
                     [11.82, 14.20, 17.45],
                     [372.48, 385.34, 379.25],
                     [149.60, 106.82, 122.36],
                     [98.02, 83.81, 69.55] ])
# si TMEM119
df_tmem = np.array([ [0.67, 0.71, 0.70],
                     [6.52, 5.65, 5.82],
                     [202.81, 182.53, 196.99],
                     [74.75, 77.60, 69.79],
                     [66.30, 64.27, 68.50] ])

# Time points
xt = np.array([0, 12, 24, 36, 48])

ax[2].errorbar(x=xt,y=np.mean(df_ctrl, axis=1), yerr=np.std(df_ctrl, axis=1, ddof=1)/np.sqrt(3), marker="o", ms=3, capsize=3, lw=2, capthick=2, label="si-ctrl")
ax[2].errorbar(x=xt,y=np.mean(df_tmem, axis=1), yerr=np.std(df_tmem, axis=1, ddof=1)/np.sqrt(3), marker="o", ms=3, capsize=3, lw=2, capthick=2, label="si-TMEM119")
ax[2].set_xticks(xt)
ax[2].set_xlabel("Hours post stimulation")
ax[2].set_ylabel("ORF57 relative expression")
ax[2].set_ylim(0,400)
#ax[2].legend(edgecolor="black")

# one-sided welch t-test
ttest3 = stats.ttest_ind(df_ctrl[2,:], df_tmem[2,:], axis=None, equal_var=False)
print(ttest3)
ttest3[1]/2

# ORF74 ===============
# si control
df_ctrl = np.array([ [1.02, 1.08, 0.91],
                     [44.21, 48.51, 41.74],
                     [44.15, 33.60, 34.35],
                     [34.04, 25.29, 32.10] ])
# si TMEM119
df_tmem = np.array([ [0.67, 0.71, 0.70],
                     [13.86, 13.68, 14.73],
                     [10.14, 13.32, 13.24],
                     [13.95, 12.97, 13.45] ])

# Time points
xt = np.array([0, 24, 36, 48])

ax[3].errorbar(x=xt,y=np.mean(df_ctrl, axis=1), yerr=np.std(df_ctrl, axis=1, ddof=1)/np.sqrt(3), marker="o", ms=3, capsize=3, lw=2, capthick=2, label="si-ctrl")
ax[3].errorbar(x=xt,y=np.mean(df_tmem, axis=1), yerr=np.std(df_tmem, axis=1, ddof=1)/np.sqrt(3), marker="o", ms=3, capsize=3, lw=2, capthick=2, label="si-TMEM119")
ax[3].set_xticks(np.array([0, 12, 24, 36, 48]))
ax[3].set_xlabel("Hours post stimulation")
ax[3].set_ylabel("ORF74 relative expression")
ax[3].set_ylim(0,50)

# one-sided welch t-test
ttest4 = stats.ttest_ind(df_ctrl[-1,:], df_tmem[-1,:], axis=None, equal_var=False)
print(ttest4)
ttest4[1]/2

sns.despine(fig=fig)
fig.tight_layout()
fig.savefig(my_dir+'over_si.svg', bbox_inches='tight')
# fig.savefig(my_dir+'over_si.png', bbox_inches='tight', dpi=300)
# plt.show()
