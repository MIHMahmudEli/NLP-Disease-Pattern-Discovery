# Disease-Specific Knowledge Mapping (Unsupervised Clustering)

This project performs unsupervised clustering on biomedical literature (PubMed abstracts) to identify underlying disease patterns and topics across different medical conditions like Diabetes, Cancer, Cardiovascular, and Neurological diseases.

## 🚀 Project Overview
The goal is to map large-scale medical knowledge by grouping related documents together without prior labeling. We use advanced NLP and dimensionality reduction techniques to visualize and interpret these clusters.

## 🛠️ Technology Stack
- **Language**: R
- **Key Libraries**: `tidyverse` (Data manipulation), `tidytext` (NLP), `uwot` (UMAP), `dbscan` (HDBSCAN), `cluster` (K-means).

## 📋 Pipeline Workflow

### 1. Data Loading & Preparation
- **Input**: `pubmed_dataset.csv`
- **Action**: Load the dataset and extract the `Abstract` column as the primary source of text.

### 2. Preprocessing
- **Cleaning**: Lowercasing, expanding contractions (e.g., "don't" → "do not"), removing HTML tags (e.g., `<b>`, `<i>`), punctuation, and numbers.
- **Stemming**: Reducing words to their root form (e.g., "clustering" → "cluster") using the Snowball algorithm.

### 3. Feature Extraction (TF-IDF)
- **TF-IDF**: Calculates the importance of words relative to documents.
- **Filtering**: Words appearing in fewer than 10 documents are removed to reduce noise and speed up processing.
- **Output**: A high-dimensional sparse matrix representing document content.

### 4. Dimensionality Reduction (UMAP)
- **Algorithm**: UMAP (Uniform Manifold Approximation and Projection) using the `uwot` engine.
- **Metric**: Cosine Similarity (best for text data).
- **Goal**: Compress 8,000+ dimensions down to a 2D plane while preserving local structures.

### 5. Clustering
- **K-Means**: Partitions the 2D data into 4 distinct clusters (targeting the major disease types).
- **HDBSCAN**: A density-based method used to identify clusters of varying shapes and handle outliers.

### 6. Evaluation
- **Silhouette Score**: Measures how well-defined and separated the clusters are (Scale -1 to 1).
- **Davies-Bouldin Index**: Measures the ratio of within-cluster scatter to between-cluster separation (Lower is better).

### 7. Keyword Extraction & Mapping
- **Interpretation**: Extracts the top 10 words per cluster to identify the medical theme (e.g., Cluster 1 = "Oncology").
- **Knowledge Mapping**: Compares identified clusters with original disease labels to validate the model's accuracy.

## 📊 How to Run
1. Ensure `pubmed_dataset.csv` is in the project directory.
2. Open `disease_mapping.R` in RStudio.
3. Run the entire script. It will automatically install missing dependencies.

## 📁 Outputs
- `umap_clusters_plot.png`: A scatter plot showing document clusters in 2D space.
- `disease_mapping_plot.png`: A bar chart validating how well the clusters align with known disease categories.
- **Console Log**: Detailed step-by-step progress and data previews.

---
**Author**: Antigravity AI  
**Date**: May 2026
