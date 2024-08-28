#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 09:54:25 2024

@author: Ferenc Kagan
Description: This script processes metabolomics data for PCA and OPLS-DA analysis,
             followed by pathway analysis using Reactome and KEGG pathways.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pypls import cross_validation, plotting
import sspa

# Change directory to the specified path
os.chdir("/Users/ferenc.kagan/Documents/Projects/Baiqun_metabolomics/")

# Load the Excel files
df1 = pd.read_excel('input/240717_Baiqun_large_small_worms_Bz_QE_Plus_TF_R_sample_amount_Norm_tables.xlsx', 
                    sheet_name="normalised_data")
df2 = pd.read_excel('input/240717_Baiqun_large_small_worms_IC_HF_TF_R_sample_amount_Norm_tables.xlsx', 
                    sheet_name="normalised_data")

# Combine the dataframes
df = pd.concat([df1, df2])

# Print the shape and columns of the DataFrame
print(f"Data shape: {df.shape}")
print(f"Columns: {df.columns}")

# Check for missing values
print(f"Missing values:\n{df.isna().sum()}")

# Filter out columns containing 'Pool'
filtered_cols = [col for col in df.columns if 'Pool' not in col]
filtered_df = df[filtered_cols].dropna()

# Display the filtered DataFrame
print(filtered_df.head())

# Extract metabolomics data matrix (assuming first column is sample names)
metab_mat = filtered_df.iloc[:, 1:].values.T  # Transpose for PCA

# Create categories based on sample names
categories = ['small' if 'small' in sample else 'large' for sample in filtered_df.columns[1:]]

####################
### PCA ANALYSIS ###
####################

# Perform PCA
n_components = 5
pca = PCA(n_components=n_components)
metab_pca = pca.fit_transform(metab_mat)

# Create plots for PCA
fig, axes = plt.subplots(2, 1, figsize=(12, 12))

# Scree plot
axes[0].plot(range(1, n_components + 1), pca.explained_variance_ratio_, marker='o')
axes[0].set_title('Scree Plot')
axes[0].set_xlabel('Principal Component')
axes[0].set_ylabel('Explained Variance Ratio')

# PCA scatter plot by category
categories_df = pd.DataFrame({
    'PC1': metab_pca[:, 0], 
    'PC2': metab_pca[:, 1], 
    'Category': categories
})
sns.scatterplot(x='PC1', y='PC2', hue='Category', data=categories_df, ax=axes[1], palette="Set1", legend='full')
axes[1].set_title('PCA of Metabolomics Measurement')
axes[1].set_xlabel(f'PC1 ({round(pca.explained_variance_ratio_[0] * 100, 2)}%)')
axes[1].set_ylabel(f'PC2 ({round(pca.explained_variance_ratio_[1] * 100, 2)}%)')

# Save and display the PCA plot
plt.savefig('output/pca.png', dpi=300, bbox_inches='tight')
plt.show()

########################
### OPLS-DA ANALYSIS ###
########################

# Prepare data for OPLS-DA
transposed_df = filtered_df.iloc[:, 1:].T
transposed_df.columns = filtered_df.iloc[:, 0]

# Initialize cross-validation for OPLS-DA
cv = cross_validation.CrossValidation(kfold=3, estimator="opls")

# Standardize the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(transposed_df.values)

# Define target labels
y = ["large"] * 3 + ["small"] * 3

# Fit and evaluate OPLS-DA model
cv.fit(scaled_data, y)
cv.permutation_test()
cv.p(metric="q2")

# Plotting OPLS-DA results
plots = plotting.Plots(cv)
plots.plot_cv_errors()
plots.jackknife_loading_plot()
plots.splot()
plots.vip_plot()
plots.permutation_plot()

# Plotting VIP scores
vip_scores = cv.vip
sorted_indices = np.argsort(vip_scores)[::-1]
sorted_vip_scores = vip_scores[sorted_indices]
sorted_feature_names = np.array(transposed_df.columns)[sorted_indices]

plt.figure(figsize=(10, 6))
plt.barh(np.arange(len(sorted_vip_scores)), sorted_vip_scores)
plt.yticks(np.arange(len(sorted_vip_scores)), sorted_feature_names, fontsize=3.5)
plt.axvline(x=1, color='red', linestyle='--', linewidth=2, label='VIP score\ncut-off')
plt.title('VIP Scores')
plt.xlabel('VIP Score')
plt.ylabel('Metabolites')
plt.legend(bbox_to_anchor=(0.8, 0.15))
plt.savefig('output/vip.png', dpi=300, bbox_inches='tight')
plt.show()

# Identify important features
important_features = np.where(vip_scores > 1)[0]
print(f"Important features: {important_features}")

# Plot clustered heatmaps
plt.figure(figsize=(12, 8))
sns.clustermap(transposed_df.iloc[:, 1:].T, cmap='magma', annot=False, linewidths=0.8, method='average', 
               metric='euclidean', z_score=0, dendrogram_ratio=(.1, .1), cbar_pos=(-0.1, 0.4, 0.03, 0.3))
plt.savefig('output/heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure(figsize=(12, 8))
sns.clustermap(transposed_df.iloc[:, important_features].T, cmap='magma', annot=False, linewidths=0.8, method='average', 
               metric='euclidean', z_score=0, dendrogram_ratio=(.1, .1), cbar_pos=(-0.1, 0.4, 0.03, 0.3))
plt.savefig('output/selected_features_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

# Export unscaled data to CSV
transposed_df.to_csv('output/metabolomics_areas_unscaled.csv')

########################
### PATHWAY ANALYSIS ###
########################

# Load and process pathway data
reactome_pathways = sspa.process_reactome(organism="Mus musculus")
kegg_human_pathways = sspa.process_kegg(organism="hsa")

# Load conversion table manually processed via Metaboanalyst
conversion_table = pd.read_csv("input/name_map.csv")
conversion_table.ChEBI = conversion_table.ChEBI.astype(str).str.split('.', expand=True)[0]

# Map identifiers to dataset and perform over-representation analysis (ORA)
significant_matrix = transposed_df.iloc[:, important_features]
processed_data_mapped = sspa.map_identifiers(conversion_table, output_id_type="ChEBI", matrix=significant_matrix)
processed_data_mapped = processed_data_mapped.drop(columns=processed_data_mapped.filter(regex='nan').columns)

ora = sspa.sspa_ora(processed_data_mapped, y, reactome_pathways, 0.05, DA_testtype='ttest')
significant_results = ora.DA_test_res.loc[ora.DA_test_res["P-adjust"] <= 0.05]
ora_res = ora.over_representation_analysis()