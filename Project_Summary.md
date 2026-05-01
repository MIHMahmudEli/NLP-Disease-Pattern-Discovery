# Project Summary: Disease-Specific Knowledge Mapping

## 1. Project Objective
The goal of this project was to perform unsupervised clustering on a large-scale dataset of biomedical abstracts (PubMed) to identify underlying disease-related topics and patterns across medical literature. By using natural language processing (NLP) and machine learning, we successfully organized 8,000+ medical documents into distinct knowledge domains without manual labeling.

## 2. Methodology (The Pipeline)
1.  **Data Preprocessing**:
    *   Removed HTML tags, numbers, and punctuation.
    *   Expanded contractions and handled lemmatization (converting words like "patients" to "patient").
    *   **Domain-Specific Filtering**: Removed non-informative medical jargon (e.g., "study", "analysis") and publication metadata (e.g., "doi", journal names) to ensure high-quality clusters.
2.  **Feature Extraction**: 
    *   Used **TF-IDF** (Term Frequency-Inverse Document Frequency) to weigh words based on their unique importance to specific documents.
    *   Applied sparsity reduction to focus on the most impactful medical terms.
3.  **Dimensionality Reduction**: 
    *   Used **UMAP** (Uniform Manifold Approximation and Projection) with a cosine metric to reduce high-dimensional text data into a 2D space for visualization.
4.  **Clustering**: 
    *   Applied **KMeans** (centroid-based) and **HDBSCAN** (density-based) clustering.
    *   Optimized the number of clusters ($k=8$) using the **Elbow Method**.

## 3. Results & Cluster Interpretation
The algorithm discovered the following primary biomedical domains:

| Cluster | Identified Disease Domain | Primary Keywords / "Signal" |
| :--- | :--- | :--- |
| **Cluster 1** | **Molecular Pathologies** | `ferroptosis`, `macrophage`, `peroxidation` |
| **Cluster 2** | **Neuro-Rehabilitation** | `eeg`, `speech`, `poststroke`, `walk` |
| **Cluster 3** | **Genomic Oncology** | `chromatin`, `transcriptomics`, `hcc` (Liver Cancer) |
| **Cluster 4** | **Neurodegenerative Disease**| `α-syn` (Parkinson's), `fibril`, `autophagy` |
| **Cluster 5** | **General Bio-Molecular** | `protein`, `kinase`, `signaling`, `receptor` |
| **Cluster 6** | **Metabolic & Obesity** | `insulin`, `glucose`, `obesity`, `adipose` |
| **Cluster 7** | **Vascular Surgery** | `aneurysm`, `stent`, `percutaneous` |
| **Cluster 8** | **Clinical Risk & Outcomes** | `mace`, `nihss` (Stroke), `readmission` |

## 4. Performance Metrics
*   **Optimal k**: 8 (based on Elbow curve and medical domain count).
*   **Cluster Quality**: Verified via Silhouette scores and Davies-Bouldin Index.
*   **Visualization**: Successfully mapped in 2D space with dynamic labeling in `disease_knowledge_map.png`.

## 5. Conclusion
The project demonstrates that unsupervised machine learning can effectively organize medical literature into meaningful categories. This automated mapping allows researchers to quickly navigate large volumes of data, identify relationships between conditions (e.g., the link between ferroptosis and oncology), and track clinical outcomes across different disease clusters.
