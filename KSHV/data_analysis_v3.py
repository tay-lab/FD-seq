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

my_dir = 'Z:/Van/20190501 - BC3 293 KSHV bulk RT-qPCR/'

# =============================================================================
# # Calculate relative mRNA expression level by delta delta Ct method
# =============================================================================
# Import data as data frame
df = pd.read_csv(my_dir+'Van_2019-05-01 14-45-19_CT030406 -  Quantification Cq Results.csv')
df = df[['Well','Cq']] # retain well number and Cq values only

# Only keep rows with gene
#my_genes = {'A':'GAPDH', 'B':'ORF26', 'C':'ISCU', 'D':'CORO1C', 'E':'KCTD2', 'F':'TMEM119', 'G':'CDH1'}
my_genes = {'A':'GAPDH', 'B':'ORF26', 'C':'ISCU', 'D':'CORO1C', 'F':'TMEM119', 'G':'CDH1'}
df = df.loc[[any(s1 in s2 for s1 in my_genes.keys()) for s2 in df['Well']],:]

# Fill not-detected values with 40
#df.loc[df['Cq'].isna(),'Cq'] = 40

# Add a column to describe target gene
df['Gene'] = [my_genes[s[0]] for s in df['Well']]

# Add a column to describe cell line
df['Cell line'] = '293'
df.loc[[int(s[1:]) > 12 for s in df['Well']],'Cell line'] = 'BC3'
# Add a column to describe whether the cells were sorted (the same as live or fixed)
df['Sort'] = 'Non-sorted'
df.loc[[(int(s[1:])-1) % 12 >= 6 for s in df['Well']],'Sort'] = 'Sorted'
# Add a column to describe status of infection
df['Treatment'] = '-'
df.loc[[int(s[1:]) % 2 == 0 for s in df['Well']],'Treatment'] = '+'

# Concatenate cell line and LorF
df['Group'] = df['Cell line'] + '_' + df['Sort']

# Calculate 2^(- delta delta Ct): ref https://doi.org/10.1006/meth.2001.1262
for i in df['Group'].unique().tolist():
    mean_Ct_GAPDH =  np.mean(df.loc[(df['Group']==i) & (df['Treatment']=='-') & (df['Gene']=='GAPDH'),'Cq']) # calculate mean Ct of GAPDH of untreated cells
    
    for j in my_genes.values():
        mean_Ct_gene =  np.mean(df.loc[(df['Group']==i) & (df['Treatment']=='-') & (df['Gene']==j),'Cq']) # calculate mean Ct of the gene of untreated cells
        
        # delta delta Ct
        df.loc[(df['Group']==i) & (df['Gene']==j),'ddCt'] = (df.loc[(df['Group']==i) & (df['Gene']==j),'Cq'].values - df.loc[(df['Group']==i) & (df['Gene']=='GAPDH'),'Cq'].values) - (mean_Ct_gene - mean_Ct_GAPDH)
        
        # SD of delta delta Ct

# Calculate relative expression
df['Rel_expr'] = 2**(-df['ddCt'])

# =============================================================================
# # Plot relative mRNA expression level
# =============================================================================
## Plot bar chart
#fig = plt.figure(figsize=(8,8))
#for counter, value in enumerate(df['Group'].unique().tolist()):
#    ax = fig.add_subplot(2,2,counter+1)
#    ax.set_yscale('log')
#    sns.barplot(x='Gene', y='Rel_expr', hue='Treatment', data=df.loc[(df['Group']==value) & (df['Gene']!='GAPDH'),:],
#                order=['ISCU','CDH1','CORO1C','TMEM119','ORF26'], ci='sd', ax=ax)
#    
#    xticklabel = ax.get_xticklabels()
#    ax.set_xticklabels(xticklabel, rotation=45, horizontalalignment='right')
#    ax.xaxis.label.set_visible(False)
#    ax.set_ylabel('Relative expression')
#    ax.set_title(value)
#    if (counter+1) == 2:
#        ymin, ymax = ax.get_ylim()
#        ax.set_ylim(1e-1,1e2)
#    elif (counter+1) == 4:
#        ymin, ymax = ax.get_ylim()
#        ax.set_ylim(1e-1,ymax)
#    
#    sns.despine(ax=ax)
#    
#fig.tight_layout()
#fig.savefig(my_dir+'data_analysis_v2/relative_expression.svg', bbox_inches='tight')
#fig.savefig(my_dir+'data_analysis_v2/relative_expression.png', bbox_inches='tight', dpi=300)


# =============================================================================
# # Import repeated experiments on 6/12/2019
# # Genes: GAPDH, CDH1, TMEM119, TRIM43 (new primers)
# # Only TMEM119 and GAPDH are usable
# =============================================================================
df2 = pd.read_csv('Z:/Van/20190612 - BC3 293 KSHV bulk RT-qPCR/Van_2019-06-12 15-08-08_CT030406 -  Quantification Cq Results.csv')
df2 = df2[['Well','Cq']] # retain well number and Cq values only

# Only keep rows with gene
#my_genes2 = {'A':'GAPDH', 'B':'CDH1', 'C':'TMEM119', 'D':'KCTD2', 'E':'DUX4', 'F':'TRIM43'}
my_genes2 = {'A':'GAPDH', 'B':'CDH1', 'C':'TMEM119', 'D':'TRIM43'}
df2 = df2.loc[[any(s1 in s2 for s1 in my_genes2.keys()) for s2 in df2['Well']],:]

# Add a column to describe target gene
df2['Gene'] = [my_genes2[s[0]] for s in df2['Well']]

# Fill not-detected values for TMEM119 with 40
df2.loc[df2['Cq'].isna() & (df2['Gene']=='TMEM119'),'Cq'] = 40

# Add a column to describe cell line
df2['Cell line'] = '293'
df2.loc[[int(s[1:]) > 12 for s in df2['Well']],'Cell line'] = 'BC3'
# Add a column to describe whether the cells were sorted (the same as live or fixed)
df2['Sort'] = 'Non-sorted'
df2.loc[[(int(s[1:])-1) % 12 >= 6 for s in df2['Well']],'Sort'] = 'Sorted'
# Add a column to describe status of infection
df2['Treatment'] = '-'
df2.loc[[int(s[1:]) % 2 == 0 for s in df2['Well']],'Treatment'] = '+'

# Concatenate cell line and LorF
df2['Group'] = df2['Cell line'] + '_' + df2['Sort']

# Calculate 2^(- delta delta Ct): ref https://doi.org/10.1006/meth.2001.1262
for i in df2['Group'].unique().tolist():
    mean_Ct_GAPDH =  np.mean(df2.loc[(df2['Group']==i) & (df2['Treatment']=='-') & (df2['Gene']=='GAPDH'),'Cq']) # calculate mean Ct of GAPDH of untreated cells
    
    for j in my_genes2.values():
        mean_Ct_gene =  np.mean(df2.loc[(df2['Group']==i) & (df2['Treatment']=='-') & (df2['Gene']==j),'Cq']) # calculate mean Ct of the gene of untreated cells
        
        # delta delta Ct
        df2.loc[(df2['Group']==i) & (df2['Gene']==j),'ddCt'] = (df2.loc[(df2['Group']==i) & (df2['Gene']==j),'Cq'].values - df2.loc[(df2['Group']==i) & (df2['Gene']=='GAPDH'),'Cq'].values) - (mean_Ct_gene - mean_Ct_GAPDH)

# Calculate relative expression
df2['Rel_expr'] = 2**(-df2['ddCt'])

# Merge df and df2
df_temp = df.loc[[s in ['ORF26', 'ISCU', 'CORO1C'] for s in df['Gene']],:]
df3 = pd.concat([df_temp, df2.loc[df2['Gene']=='TMEM119',:]], sort=False)
#df3 = pd.concat([df_temp, df2], sort=False)

# =============================================================================
# # Import repeated experiments on 6/14/2019
# # Genes: GAPDH, CDH1, TRIM43 (new primers, different from 6/12/2019)
# # CDH1 and TRIM43 didn't work for sorted samples
# =============================================================================
# =============================================================================
# df2 = pd.read_csv('Z:/Van/20190614 - BC3 293 KSHV bulk RT-qPCR/Van_2019-06-14 14-15-49_CT030406 -  Quantification Cq Results.csv')
# df2 = df2[['Well','Cq']] # retain well number and Cq values only
# 
# # Only keep rows with gene
# #my_genes2 = {'A':'GAPDH', 'B':'CDH1', 'C':'TMEM119', 'D':'KCTD2', 'E':'DUX4', 'F':'TRIM43'}
# my_genes2 = {'A':'GAPDH', 'B':'CDH1', 'C':'TRIM43'}
# df2 = df2.loc[[any(s1 in s2 for s1 in my_genes2.keys()) for s2 in df2['Well']],:]
# 
# # Add a column to describe target gene
# df2['Gene'] = [my_genes2[s[0]] for s in df2['Well']]
# 
# # Fill not-detected values with 40
# #df2.loc[df2['Cq'].isna(),'Cq'] = 40
# 
# # Add a column to describe cell line
# df2['Cell line'] = '293'
# df2.loc[[int(s[1:]) > 12 for s in df2['Well']],'Cell line'] = 'BC3'
# # Add a column to describe whether the cells were sorted (the same as live or fixed)
# df2['Sort'] = 'Non-sorted'
# df2.loc[[(int(s[1:])-1) % 12 >= 6 for s in df2['Well']],'Sort'] = 'Sorted'
# # Add a column to describe status of infection
# df2['Treatment'] = '-'
# df2.loc[[int(s[1:]) % 2 == 0 for s in df2['Well']],'Treatment'] = '+'
# 
# # Concatenate cell line and LorF
# df2['Group'] = df2['Cell line'] + '_' + df2['Sort']
# 
# # Calculate 2^(- delta delta Ct): ref https://doi.org/10.1006/meth.2001.1262
# for i in df2['Group'].unique().tolist():
#     mean_Ct_GAPDH =  np.mean(df2.loc[(df2['Group']==i) & (df2['Treatment']=='-') & (df2['Gene']=='GAPDH'),'Cq']) # calculate mean Ct of GAPDH of untreated cells
#     
#     for j in my_genes2.values():
#         mean_Ct_gene =  np.mean(df2.loc[(df2['Group']==i) & (df2['Treatment']=='-') & (df2['Gene']==j),'Cq']) # calculate mean Ct of the gene of untreated cells
#         
#         # delta delta Ct
#         df2.loc[(df2['Group']==i) & (df2['Gene']==j),'ddCt'] = (df2.loc[(df2['Group']==i) & (df2['Gene']==j),'Cq'].values - df2.loc[(df2['Group']==i) & (df2['Gene']=='GAPDH'),'Cq'].values) - (mean_Ct_gene - mean_Ct_GAPDH)
# 
# # Calculate relative expression
# df2['Rel_expr'] = 2**(-df2['ddCt'])
# 
# # Merge df3 and df2
# df3 = pd.concat([df3, df2.loc[df2['Gene']!='GAPDH',:]], sort=False)
# =============================================================================

# =============================================================================
# # Import repeated experiments on 7/1/2019
# # Genes: GAPDH, CDH1, TRIM43 (new primers, different from 6/14/2019)
# # TRIM43, CDH1 didn't work for sorted BC3 -> try pre-amp
# # Use TRIM43 and CDH1 for 293 and non-sorted BC3
# =============================================================================
# =============================================================================
# df2 = pd.read_csv('Z:/Van/20190701 - RT-qPCR/Van_2019-07-01 11-36-30_CT030406 -  Quantification Cq Results.csv')
# df2 = df2[['Well','Cq']] # retain well number and Cq values only
# 
# # Only keep rows with gene
# #my_genes2 = {'A':'GAPDH', 'B':'CDH1', 'C':'TMEM119', 'D':'KCTD2', 'E':'DUX4', 'F':'TRIM43'}
# my_genes2 = {'A':'GAPDH', 'B':'CDH1', 'C':'TRIM43'}
# df2 = df2.loc[[any(s1 in s2 for s1 in my_genes2.keys()) for s2 in df2['Well']],:]
# 
# # Add a column to describe target gene
# df2['Gene'] = [my_genes2[s[0]] for s in df2['Well']]
# 
# # Add a column to describe cell line
# df2['Cell line'] = '293'
# df2.loc[[int(s[1:]) > 12 for s in df2['Well']],'Cell line'] = 'BC3'
# # Add a column to describe whether the cells were sorted (the same as live or fixed)
# df2['Sort'] = 'Non-sorted'
# df2.loc[[(int(s[1:])-1) % 12 >= 6 for s in df2['Well']],'Sort'] = 'Sorted'
# # Add a column to describe status of infection
# df2['Treatment'] = '-'
# df2.loc[[int(s[1:]) % 2 == 0 for s in df2['Well']],'Treatment'] = '+'
# 
# # Fill not-detected values with 40 in 293 cells
# df2.loc[(df2['Cq'].isna()) & (df2['Cell line']=='293'),'Cq'] = 40
# 
# # Concatenate cell line and LorF
# df2['Group'] = df2['Cell line'] + '_' + df2['Sort']
# 
# # Calculate 2^(- delta delta Ct): ref https://doi.org/10.1006/meth.2001.1262
# for i in df2['Group'].unique().tolist():
#     mean_Ct_GAPDH =  np.mean(df2.loc[(df2['Group']==i) & (df2['Treatment']=='-') & (df2['Gene']=='GAPDH'),'Cq']) # calculate mean Ct of GAPDH of untreated cells
#     
#     for j in my_genes2.values():
#         mean_Ct_gene =  np.mean(df2.loc[(df2['Group']==i) & (df2['Treatment']=='-') & (df2['Gene']==j),'Cq']) # calculate mean Ct of the gene of untreated cells
#         
#         # delta delta Ct
#         df2.loc[(df2['Group']==i) & (df2['Gene']==j),'ddCt'] = (df2.loc[(df2['Group']==i) & (df2['Gene']==j),'Cq'].values - df2.loc[(df2['Group']==i) & (df2['Gene']=='GAPDH'),'Cq'].values) - (mean_Ct_gene - mean_Ct_GAPDH)
# 
# # Calculate relative expression
# df2['Rel_expr'] = 2**(-df2['ddCt'])
# 
# # Merge df3 and df2
# df3 = pd.concat([df3, df2.loc[df2['Gene']!='GAPDH',:]], sort=False)
# =============================================================================

# =============================================================================
# # Import repeated experiments on 7/3/2019
# # Genes: GAPDH, CDH1, TRIM43 (BC3 pre-amp)
# =============================================================================
# =============================================================================
# df2 = pd.read_csv('Z:/Van/20190703 - BC3 KSHV bulk RT-qPCR pre-amp/Jing_2019-07-03 14-15-53_CT030406 -  Quantification Cq Results.csv')
# df2 = df2[['Well','Cq']] # retain well number and Cq values only
# 
# # Only keep rows with gene
# #my_genes2 = {'A':'GAPDH', 'B':'CDH1', 'C':'TMEM119', 'D':'KCTD2', 'E':'DUX4', 'F':'TRIM43'}
# my_genes2 = {'A':'GAPDH', 'B':'CDH1', 'C':'TRIM43'}
# df2 = df2.loc[[any(s1 in s2 for s1 in my_genes2.keys()) & (int(s2[1:])<=6) for s2 in df2['Well']],:]
# 
# # Add a column to describe target gene
# df2['Gene'] = [my_genes2[s[0]] for s in df2['Well']]
# 
# # Add a column to describe cell line
# df2['Cell line'] = 'BC3'
# # Add a column to describe whether the cells were sorted (the same as live or fixed)
# df2['Sort'] = 'Sorted'
# # Add a column to describe status of infection
# df2['Treatment'] = '-'
# df2.loc[[int(s[1:]) % 2 == 0 for s in df2['Well']],'Treatment'] = '+'
# 
# # Fill not-detected values with 40 in 293 cells
# df2.loc[df2['Cq'].isna(),'Cq'] = 40
# 
# # Concatenate cell line and LorF
# df2['Group'] = df2['Cell line'] + '_' + df2['Sort']
# 
# # Calculate 2^(- delta delta Ct): ref https://doi.org/10.1006/meth.2001.1262
# for i in df2['Group'].unique().tolist():
#     mean_Ct_GAPDH =  np.mean(df2.loc[(df2['Group']==i) & (df2['Treatment']=='-') & (df2['Gene']=='GAPDH'),'Cq']) # calculate mean Ct of GAPDH of untreated cells
#     
#     for j in my_genes2.values():
#         mean_Ct_gene =  np.mean(df2.loc[(df2['Group']==i) & (df2['Treatment']=='-') & (df2['Gene']==j),'Cq']) # calculate mean Ct of the gene of untreated cells
#         
#         # delta delta Ct
#         df2.loc[(df2['Group']==i) & (df2['Gene']==j),'ddCt'] = (df2.loc[(df2['Group']==i) & (df2['Gene']==j),'Cq'].values - df2.loc[(df2['Group']==i) & (df2['Gene']=='GAPDH'),'Cq'].values) - (mean_Ct_gene - mean_Ct_GAPDH)
# 
# # Calculate relative expression
# df2['Rel_expr'] = 2**(-df2['ddCt'])
# 
# # Merge df3 and df2
# df3 = df3.loc[~((df3['Group']=='BC3_Sorted') & ((df3['Gene']=='CDH1') | (df3['Gene']=='TRIM43'))),:]
# df3 = pd.concat([df3, df2.loc[df2['Gene']!='GAPDH',:]], sort=False)
# =============================================================================

# =============================================================================
# # Plot relative mRNA expression level
# =============================================================================
# Plot bar chart
fig = plt.figure(figsize=(6.2,6.2))
for counter, value in enumerate(df3['Group'].unique().tolist()):
    ax = fig.add_subplot(2,2,counter+1)
    ax.set_yscale('log')
    sns.barplot(x='Gene', y='Rel_expr', hue='Treatment', data=df3.loc[(df3['Group']==value) & (df3['Gene']!='GAPDH'),:],
                palette=sns.color_palette('Dark2',2)[::-1],
                order=['ISCU','CORO1C','TMEM119','ORF26'],
                ci='sd', errcolor='black', capsize=0.075, ax=ax)
    sns.swarmplot(x='Gene', y='Rel_expr', hue='Treatment', data=df3.loc[(df3['Group']==value) & (df3['Gene']!='GAPDH'),:],
                  color='black', dodge=True, size=4.5,
                  order=['ISCU','CORO1C','TMEM119','ORF26'], ax=ax)
    
    xticklabel = ax.get_xticklabels()
    ax.set_xticklabels(xticklabel, rotation=45, horizontalalignment='right')
    ax.xaxis.label.set_visible(False)
    ax.set_ylabel('Relative expression')
    ax.set_title(value)
    if value == 'BC3_Sorted':
#        ymin, ymax = ax.get_ylim()
        ax.set_ylim(1e-1,1e4)
    elif value == '293_Non-sorted':
#        ymin, ymax = ax.get_ylim()
        ax.set_ylim(1e-1,4e1)
    
    sns.despine(ax=ax)
fig.tight_layout()
fig.savefig(my_dir+'data_analysis_v3/relative_expression.svg', bbox_inches='tight')
# fig.savefig(my_dir+'data_analysis_v3/relative_expression.png', bbox_inches='tight', dpi=300)

plt.show()