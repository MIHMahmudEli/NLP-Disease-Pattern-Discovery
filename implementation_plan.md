# Implementation Plan: Disease-Specific Knowledge Mapping (Unsupervised Clustering) in R

This project aims to group biomedical documents (from `pubmed_dataset.csv`) into clusters using unsupervised learning to identify disease-related topics and patterns.

## 1. Environment Setup
- Install and load necessary R packages:
  - `tidyverse`: Data manipulation and visualization.
  - `tidytext`: Text mining and TF-IDF calculation.
  - `tm`: Text preprocessing.
  - `umap`: Dimensionality reduction.
  - `cluster`: K-means and silhouette analysis.
  - `dbscan`: HDBSCAN clustering.
  - `fpc`: Cluster validation.
  - `wordcloud`: Cluster interpretation.

## 2. Data Loading & Preparation
- Read `pubmed_dataset.csv`.
- Combine `Title` and `Abstract` into a single text column.
- Handle missing values if any.

## 3. Preprocessing
- Tokenization.
- Convert to lowercase.
- Remove punctuation, numbers, and common stop words.
- Apply stemming or lemmatization to reduce word variations.

## 4. Feature Extraction (TF-IDF)
- Calculate Term Frequency-Inverse Document Frequency (TF-IDF) to convert text into a numerical matrix.
- Filter out extremely rare or common words to reduce noise.

## 5. Dimensionality Reduction (UMAP)
- Apply UMAP to the TF-IDF matrix to project high-dimensional text data into 2D or 3D space.
- UMAP is preferred over PCA for preserving local and global structures in non-linear data.

## 6. Clustering
- **K-Means**: Find an optimal number of clusters (K) using the Elbow method or Silhouette analysis.
- **HDBSCAN**: An alternative density-based method that handles noise and varying cluster densities.

## 7. Evaluation
- **Silhouette Score**: Measure how similar an object is to its own cluster compared to other clusters.
- **Davies–Bouldin Index**: Evaluate the separation and compactness of clusters.

## 8. Knowledge Mapping & Interpretation
- Identify top keywords for each cluster.
- Map clusters to disease categories (using the ground truth `Disease` column for validation if available).
- Generate word clouds for each cluster.

## 9. Visualization
- Plot the UMAP coordinates colored by cluster.
- Interactive visualization (optional, using `plotly`).
