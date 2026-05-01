# Disease-Specific Knowledge Mapping (Unsupervised Clustering)

This project groups biomedical documents (PubMed abstracts) into clusters to identify underlying disease-related topics and patterns across medical literature.

## Pipeline
1. **Preprocessing**: Text cleaning, tokenization, stop word removal, and lemmatization using `tm` and `textstem`.
2. **Feature Extraction**: TF-IDF weighting with sparsity reduction.
3. **Dimensionality Reduction**: UMAP (Uniform Manifold Approximation and Projection) for 2D visualization and clustering support.
4. **Clustering**: 
   - **KMeans**: Centroid-based clustering with Elbow Method for optimal $k$.
   - **HDBSCAN**: Density-based clustering for identifying natural groupings and noise.
5. **Evaluation**: Silhouette Score and Davies–Bouldin Index.
6. **Keyword Extraction**: Identifying top terms per cluster to interpret disease categories.
7. **Visualization**: Interactive-ready UMAP plots with cluster labels.

## Results
The project successfully maps clusters to major medical condition categories:
- Diabetes / Metabolic
- Oncology / Cancer
- Cardiovascular
- Neurology
- Immunology
- and more...

## How to Run
Ensure R is installed with the required packages:
```r
source("disease_mapping.R")
```

## Visualizations
Generated plots:
- `disease_knowledge_map.png`: 2D projection of medical literature clusters.
- `top_words_frequency.png`: Global dataset term distribution.
- `cluster_keywords_comparison.png`: Topic-specific keywords per cluster.
