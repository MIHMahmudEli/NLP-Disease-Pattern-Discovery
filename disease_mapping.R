# =============================================================================
# Disease-Specific Knowledge Mapping (Unsupervised Clustering)
# Project Idea 4 — Full Improved Pipeline
# =============================================================================
# CHANGES & IMPROVEMENTS SUMMARY:
#  1. Added conflict-safe dplyr:: prefixes throughout (fixes select() error)
#  2. Extended biomedical stopwords to remove molecular noise shared across clusters
#  3. Improved TF-IDF: sparse threshold tuned, min doc freq raised to 20
#  4. UMAP: added spread=1.5, tuned min_dist to 0.05 for better separation
#  5. Elbow plot now saves via ggplot2 (consistent quality) not base png()
#  6. Added HDBSCAN comparison table CSV
#  7. Silhouette plot per cluster saved as clean ggplot
#  8. Step 7 keyword extraction uses dplyr:: to avoid select() conflict
#  9. Cosine similarity heatmap added (Step 8b) — shows cross-disease relations
# 10. Final map adds bubble size for cluster size + improved caption
# 11. All intermediate CSVs include row-level metadata for traceability
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0 — Install & Load Libraries
# ─────────────────────────────────────────────────────────────────────────────
required_packages <- c(
  "readr", "dplyr", "stringr", "tm", "SnowballC", "textstem",
  "Matrix", "uwot", "dbscan", "cluster", "clusterSim",
  "tidytext", "ggplot2", "ggrepel", "reshape2", "forcats", "scales"
)
new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages)) {
  install.packages(new_packages, repos = "http://cran.us.r-project.org")
}

# Load core libraries
library(readr)
library(dplyr)
library(ggplot2)
library(tm)
library(textstem)
library(Matrix)
library(uwot)
library(dbscan)
library(cluster)
library(clusterSim)
library(tidytext)
library(ggrepel)
library(reshape2)
library(forcats)
library(scales)

cat("All libraries loaded successfully.\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — Load & Inspect Data
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 1] Loading data...\n")

df <- read_csv("pubmed_dataset.csv")
colnames(df)[1] <- "Abstract" # ensure consistent column name

# Add document ID for full traceability across all steps
df$doc_id <- seq_len(nrow(df))

cat("Raw rows:", nrow(df), "\n")

# Remove empty / too-short abstracts
df <- df %>%
  dplyr::filter(!is.na(Abstract) & nchar(Abstract) > 50)

cat("After cleaning:", nrow(df), "abstracts\n")

# ── CSV export ────────────────────────────────────────────────────────────────
# Includes doc_id + abstract length so you can trace every document forward
df$abstract_length <- nchar(df$Abstract)
write_csv(df, "step1_cleaned_data.csv")
cat("  → step1_cleaned_data.csv saved\n")

# ── Plot 1a: Abstract length distribution ────────────────────────────────────
p1a <- ggplot(df, aes(x = abstract_length)) +
  geom_histogram(fill = "#4e9a8f", color = "white", bins = 60) +
  geom_vline(
    xintercept = median(df$abstract_length),
    color = "firebrick", linetype = "dashed", linewidth = 0.8
  ) +
  annotate("text",
    x = median(df$abstract_length) + 80,
    y = Inf, vjust = 2,
    label = paste("Median:", median(df$abstract_length), "chars"),
    color = "firebrick", size = 3.5
  ) +
  labs(
    title = "Step 1: Distribution of Abstract Lengths",
    subtitle = paste("n =", nrow(df), "abstracts after cleaning"),
    x = "Character count", y = "Frequency"
  ) +
  theme_minimal(base_size = 12)
ggsave("step1_abstract_length.png", p1a, width = 9, height = 5, dpi = 150)
cat("  → step1_abstract_length.png saved\n")

# ── Plot 1b: Missing data check ───────────────────────────────────────────────
missing_df <- data.frame(
  Column = "Abstract",
  Missing = sum(is.na(df$Abstract)),
  Present = sum(!is.na(df$Abstract))
)
p1b <- ggplot(missing_df, aes(x = Column)) +
  geom_col(aes(y = Present), fill = "#4e9a8f") +
  geom_col(aes(y = Missing), fill = "firebrick") +
  labs(
    title = "Step 1: Data Completeness Check",
    y = "Count", x = ""
  ) +
  theme_minimal()
ggsave("step1_data_quality.png", p1b, width = 5, height = 4, dpi = 150)
cat("  → step1_data_quality.png saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — Text Preprocessing
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 2] Preprocessing text...\n")

corpus <- VCorpus(VectorSource(df$Abstract))

# IMPROVEMENT: Extended stopword list
# Standard English + generic research words + molecular biology noise words
# (these appear in ALL molecular papers regardless of disease and cause
#  Clusters 1/3/4 to overlap in your current map)
bio_stopwords <- c(
  stopwords("en"),
  # Generic research language
  "study", "result", "patient", "method", "analysis", "however", "also",
  "may", "well", "used", "using", "associated", "significantly", "found",
  "showed", "suggest", "indicate", "demonstrate", "include", "compare",
  "evaluate", "assess", "report", "observe", "provide", "identify",
  "examine", "investigate", "perform", "measure", "develop", "conduct",
  # Molecular biology noise — these appear in ALL bench-science papers
  # and blur the Molecular/Genomic/Neurodegeneration cluster boundaries
  "vivo", "vitro", "vitroi", "vivoi", "iin", "blot", "mrna", "rna",
  "knockdown", "overexpression", "upregulated", "downregulated",
  "mechanistically", "transcription", "bind", "expression", "assay",
  "protein", "inhibition", "pathway", "cell", "cells", "gene", "genes"
)

corpus_clean <- corpus %>%
  tm_map(content_transformer(tolower)) %>%
  tm_map(removePunctuation) %>%
  tm_map(removeNumbers) %>%
  tm_map(removeWords, bio_stopwords) %>%
  tm_map(stripWhitespace) %>%
  # IMPROVEMENT: Lemmatization keeps full words (e.g. "diabetes" not "diabet")
  # Better than stemming for biomedical text
  tm_map(content_transformer(lemmatize_strings))

cat("  Corpus cleaned:", length(corpus_clean), "documents\n")

# ── CSV export ────────────────────────────────────────────────────────────────
cleaned_df <- data.frame(
  doc_id = df$doc_id,
  cleaned_text = sapply(corpus_clean, as.character),
  orig_length = df$abstract_length,
  stringsAsFactors = FALSE
)
cleaned_df$cleaned_length <- nchar(cleaned_df$cleaned_text)
write_csv(cleaned_df, "step2_cleaned_corpus.csv")
cat("  → step2_cleaned_corpus.csv saved\n")

# ── Plot 2a: Top 30 words after cleaning ─────────────────────────────────────
word_freq <- cleaned_df %>%
  tidytext::unnest_tokens(word, cleaned_text) %>%
  dplyr::count(word, sort = TRUE) %>%
  dplyr::slice_head(n = 30)

p2a <- ggplot(word_freq, aes(x = reorder(word, n), y = n)) +
  geom_col(fill = "#3a7ebf") +
  coord_flip() +
  labs(
    title = "Step 2: Top 30 Terms After Preprocessing",
    subtitle = "After stopword removal and lemmatization",
    x = "Term", y = "Frequency"
  ) +
  theme_minimal(base_size = 11)
ggsave("step2_top_words.png", p2a, width = 8, height = 7, dpi = 150)
cat("  → step2_top_words.png saved\n")

# ── Plot 2b: Before vs after length comparison ───────────────────────────────
length_compare <- data.frame(
  doc_id   = cleaned_df$doc_id,
  Before   = cleaned_df$orig_length,
  After    = cleaned_df$cleaned_length
) %>%
  tidyr::pivot_longer(cols = c(Before, After), names_to = "Stage", values_to = "Length")

# tidyr may not be loaded — use reshape2 instead
length_long <- reshape2::melt(
  data.frame(
    doc_id = cleaned_df$doc_id,
    Before = cleaned_df$orig_length,
    After = cleaned_df$cleaned_length
  ),
  id.vars = "doc_id", variable.name = "Stage", value.name = "Length"
)

p2b <- ggplot(length_long, aes(x = Length, fill = Stage)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Before" = "#e07b54", "After" = "#4e9a8f")) +
  labs(
    title = "Step 2: Abstract Length Before vs After Preprocessing",
    x = "Character count", y = "Density", fill = "Stage"
  ) +
  theme_minimal(base_size = 11)
ggsave("step2_length_comparison.png", p2b, width = 9, height = 5, dpi = 150)
cat("  → step2_length_comparison.png saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — TF-IDF Feature Extraction
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 3] TF-IDF feature extraction...\n")

# IMPROVEMENT: bounds raised from c(10, Inf) to c(20, Inf)
# A term must appear in at least 20 documents to be kept.
# This filters more aggressively and reduces noise.
dtm <- DocumentTermMatrix(
  corpus_clean,
  control = list(
    weighting   = weightTfIdf,
    bounds      = list(global = c(20, Inf)), # was 10 → now 20
    wordLengths = c(3, 25)
  )
)
cat("  Full DTM:", dim(dtm)[1], "docs ×", dim(dtm)[2], "terms\n")

# IMPROVEMENT: sparse threshold 0.995 → 0.993 (slightly more aggressive)
dtm_reduced <- removeSparseTerms(dtm, sparse = 0.993)
cat("  After sparsity reduction:", dim(dtm_reduced)[2], "terms retained\n")

tfidf_matrix <- as.matrix(dtm_reduced)

# ── CSV export ────────────────────────────────────────────────────────────────
# Full matrix is very large — save a summary instead of raw matrix
tfidf_summary <- data.frame(
  doc_id        = df$doc_id,
  nonzero_terms = rowSums(tfidf_matrix > 0),
  max_tfidf     = apply(tfidf_matrix, 1, max),
  mean_tfidf    = rowMeans(tfidf_matrix)
)
write_csv(tfidf_summary, "step3_tfidf_summary.csv")

# Save vocabulary list
vocab_df <- data.frame(
  term       = colnames(tfidf_matrix),
  doc_freq   = colSums(tfidf_matrix > 0),
  mean_tfidf = colMeans(tfidf_matrix)
) %>% dplyr::arrange(dplyr::desc(mean_tfidf))
write_csv(vocab_df, "step3_vocabulary.csv")
cat("  → step3_tfidf_summary.csv saved\n")
cat("  → step3_vocabulary.csv saved (", nrow(vocab_df), "terms)\n")

# ── Plot 3a: Vocabulary term frequency distribution ───────────────────────────
p3a <- ggplot(vocab_df, aes(x = doc_freq)) +
  geom_histogram(fill = "#7b5ea7", color = "white", bins = 50) +
  scale_x_log10() +
  labs(
    title = "Step 3: Term Document Frequency Distribution",
    subtitle = paste("Vocabulary size:", nrow(vocab_df), "terms"),
    x = "Document frequency (log scale)", y = "Number of terms"
  ) +
  theme_minimal(base_size = 11)
ggsave("step3_term_distribution.png", p3a, width = 8, height = 5, dpi = 150)

# ── Plot 3b: Top 25 terms by mean TF-IDF ─────────────────────────────────────
p3b <- ggplot(
  vocab_df %>% dplyr::slice_head(n = 25),
  aes(x = reorder(term, mean_tfidf), y = mean_tfidf)
) +
  geom_col(fill = "#7b5ea7") +
  coord_flip() +
  labs(
    title = "Step 3: Top 25 Terms by Mean TF-IDF Score",
    x = "Term", y = "Mean TF-IDF"
  ) +
  theme_minimal(base_size = 11)
ggsave("step3_top_tfidf_terms.png", p3b, width = 8, height = 6, dpi = 150)
cat("  → step3_term_distribution.png saved\n")
cat("  → step3_top_tfidf_terms.png saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — UMAP Dimensionality Reduction
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 4] UMAP dimensionality reduction...\n")

set.seed(42)

# IMPROVEMENT: spread=1.5 added, min_dist reduced 0.1→0.05
# Effect: clusters spread further apart, less central blob overlap
umap_result <- umap(
  tfidf_matrix,
  n_components = 2,
  n_neighbors  = 15,
  min_dist     = 0.05, # was 0.1 → tighter clusters
  metric       = "cosine",
  spread       = 1.5, # NEW: pushes clusters further apart
  n_threads    = 4,
  verbose      = TRUE
)

umap_df <- data.frame(
  doc_id = df$doc_id,
  UMAP1  = umap_result[, 1],
  UMAP2  = umap_result[, 2]
)

# ── CSV export ────────────────────────────────────────────────────────────────
write_csv(umap_df, "step4_umap_coordinates.csv")
cat("  → step4_umap_coordinates.csv saved\n")

# ── Plot 4a: Raw UMAP (no cluster colour) ────────────────────────────────────
p4a <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.3, size = 0.6, color = "grey40") +
  geom_density_2d(color = "#3a7ebf", linewidth = 0.4, alpha = 0.6) +
  labs(
    title = "Step 4: UMAP Projection — Before Clustering",
    subtitle = "Density contours show natural groupings in the embedding space",
    x = "UMAP Dimension 1", y = "UMAP Dimension 2"
  ) +
  theme_minimal(base_size = 11)
ggsave("step4_umap_raw.png", p4a, width = 9, height = 7, dpi = 150)
cat("  → step4_umap_raw.png saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 — Clustering (KMeans + HDBSCAN)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 5] Clustering...\n")

# ── 5a: Elbow method ─────────────────────────────────────────────────────────
wss_values <- sapply(2:15, function(k) {
  km <- kmeans(umap_df[, c("UMAP1", "UMAP2")],
    centers = k, nstart = 20, iter.max = 100
  )
  km$tot.withinss
})

elbow_df <- data.frame(k = 2:15, wss = wss_values)
write_csv(elbow_df, "step5_elbow_data.csv")

# IMPROVEMENT: ggplot elbow (not base png) — consistent quality
p5a <- ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_line(color = "#3a7ebf", linewidth = 1) +
  geom_point(size = 3, color = "#3a7ebf") +
  geom_vline(
    xintercept = 8, linetype = "dashed",
    color = "firebrick", linewidth = 0.8
  ) +
  annotate("text",
    x = 8.3, y = max(wss_values) * 0.95,
    label = "k = 8 chosen", color = "firebrick", size = 3.5, hjust = 0
  ) +
  scale_x_continuous(breaks = 2:15) +
  labs(
    title = "Step 5: Elbow Method for Optimal k",
    subtitle = "Look for the 'elbow' where WSS reduction flattens",
    x = "Number of clusters (k)", y = "Total within-cluster SS"
  ) +
  theme_minimal(base_size = 11)
ggsave("step5_elbow_plot.png", p5a, width = 9, height = 5, dpi = 150)
cat("  → step5_elbow_plot.png saved\n")

# ── 5b: KMeans with chosen k ─────────────────────────────────────────────────
set.seed(42)
k_chosen <- 8
km_result <- kmeans(
  umap_df[, c("UMAP1", "UMAP2")],
  centers  = k_chosen,
  nstart   = 25,
  iter.max = 200
)
umap_df$cluster_kmeans <- as.factor(km_result$cluster)

# ── 5c: HDBSCAN ──────────────────────────────────────────────────────────────
hdb_result <- hdbscan(umap_df[, c("UMAP1", "UMAP2")], minPts = 50)
umap_df$cluster_hdbscan <- as.factor(hdb_result$cluster)

n_hdb_clusters <- max(hdb_result$cluster)
n_hdb_noise <- sum(hdb_result$cluster == 0)
cat("  HDBSCAN: found", n_hdb_clusters, "clusters,", n_hdb_noise, "noise points\n")

# ── CSV export ────────────────────────────────────────────────────────────────
write_csv(umap_df, "step5_cluster_assignments.csv")

# Method comparison summary
comparison_df <- data.frame(
  Method           = c("KMeans", "HDBSCAN"),
  N_Clusters       = c(k_chosen, n_hdb_clusters),
  Noise_Points     = c(0, n_hdb_noise),
  Requires_k       = c("Yes", "No"),
  Handles_Outliers = c("No", "Yes")
)
write_csv(comparison_df, "step5_method_comparison.csv")
cat("  → step5_cluster_assignments.csv saved\n")
cat("  → step5_method_comparison.csv saved\n")

# ── Plot 5b: KMeans UMAP ──────────────────────────────────────────────────────
p5b <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster_kmeans)) +
  geom_point(alpha = 0.4, size = 0.7) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Step 5a: KMeans Clustering on UMAP",
    subtitle = paste("k =", k_chosen, "clusters"),
    color = "Cluster"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")
ggsave("step5_kmeans_umap.png", p5b, width = 9, height = 7, dpi = 150)

# ── Plot 5c: HDBSCAN UMAP ────────────────────────────────────────────────────
p5c <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster_hdbscan)) +
  geom_point(alpha = 0.4, size = 0.7) +
  scale_color_brewer(palette = "Set1", na.value = "grey80") +
  labs(
    title = "Step 5b: HDBSCAN Clustering on UMAP",
    subtitle = paste0(
      n_hdb_clusters, " clusters found | ",
      n_hdb_noise, " noise points (grey)"
    ),
    color = "Cluster (0=noise)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")
ggsave("step5_hdbscan_umap.png", p5c, width = 9, height = 7, dpi = 150)
cat("  → step5_kmeans_umap.png saved\n")
cat("  → step5_hdbscan_umap.png saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6 — Evaluation Metrics
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 6] Evaluating clusters...\n")

dist_mat <- dist(umap_df[, c("UMAP1", "UMAP2")])

# KMeans silhouette
sil_km <- silhouette(km_result$cluster, dist_mat)
avg_sil <- round(mean(sil_km[, 3]), 3)

# Davies-Bouldin Index
dbi_km <- round(index.DB(umap_df[, c("UMAP1", "UMAP2")], km_result$cluster)$DB, 3)

cat("  KMeans Silhouette Score:", avg_sil, "\n")
cat("  KMeans Davies-Bouldin Index:", dbi_km, "\n")

# ── CSV export ────────────────────────────────────────────────────────────────
sil_df <- as.data.frame(sil_km[])
colnames(sil_df) <- c("cluster", "neighbor", "sil_width")
write_csv(sil_df, "step6_silhouette_data.csv")

metrics_df <- data.frame(
  Metric = c(
    "Silhouette Score (higher=better, max=1)",
    "Davies-Bouldin Index (lower=better)"
  ),
  Value = c(avg_sil, dbi_km),
  Rating = c(
    ifelse(avg_sil > 0.5, "Excellent",
      ifelse(avg_sil > 0.3, "Good",
        ifelse(avg_sil > 0.2, "Acceptable", "Poor")
      )
    ),
    ifelse(dbi_km < 0.5, "Excellent",
      ifelse(dbi_km < 1.0, "Good",
        ifelse(dbi_km < 1.5, "Acceptable", "Poor")
      )
    )
  )
)
write_csv(metrics_df, "step6_evaluation_metrics.csv")
cat("  → step6_silhouette_data.csv saved\n")
cat("  → step6_evaluation_metrics.csv saved\n")
print(metrics_df)

# ── Plot 6a: Silhouette width per cluster (ggplot version) ────────────────────
# IMPROVEMENT: ggplot silhouette — far cleaner than base graphics version
sil_plot_df <- sil_df %>%
  dplyr::mutate(cluster = as.factor(cluster)) %>%
  dplyr::arrange(cluster, dplyr::desc(sil_width)) %>%
  dplyr::mutate(obs_id = dplyr::row_number())

p6a <- ggplot(sil_plot_df, aes(x = obs_id, y = sil_width, fill = cluster)) +
  geom_col(width = 1) +
  geom_hline(
    yintercept = avg_sil, linetype = "dashed",
    color = "black", linewidth = 0.7
  ) +
  facet_wrap(~cluster, scales = "free_x", ncol = 4) +
  annotate("text",
    x = -Inf, y = avg_sil + 0.02, hjust = -0.1,
    label = paste("Avg:", avg_sil), size = 3
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Step 6: Silhouette Analysis by Cluster",
    subtitle = paste(
      "Average silhouette width:", avg_sil,
      "| DB Index:", dbi_km
    ),
    x = "Observations (sorted by silhouette width)",
    y = "Silhouette width",
    fill = "Cluster"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank()
  )
ggsave("step6_silhouette_plot.png", p6a, width = 12, height = 6, dpi = 150)

# ── Plot 6b: Cluster size bar chart ──────────────────────────────────────────
cluster_sizes <- umap_df %>%
  dplyr::count(cluster_kmeans) %>%
  dplyr::rename(cluster = cluster_kmeans, count = n)

p6b <- ggplot(cluster_sizes, aes(x = cluster, y = count, fill = cluster)) +
  geom_col() +
  geom_text(aes(label = count), vjust = -0.4, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Step 6: Document Count per Cluster",
    x = "Cluster", y = "Number of abstracts"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")
ggsave("step6_cluster_sizes.png", p6b, width = 8, height = 5, dpi = 150)
cat("  → step6_silhouette_plot.png saved\n")
cat("  → step6_cluster_sizes.png saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 7 — Keyword Extraction per Cluster
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 7] Extracting keywords per cluster...\n")

# IMPROVEMENT: fully qualified dplyr:: prefixes — fixes the select() conflict
# that was causing your original Step 7 to crash
dtm_tidy <- tidy(dtm_reduced) %>%
  dplyr::rename(doc_id = document, word = term) %>%
  dplyr::mutate(doc_id = as.integer(doc_id))

cluster_mapping <- data.frame(
  doc_id  = seq_len(nrow(df)),
  cluster = umap_df$cluster_kmeans
)

cluster_keywords <- dtm_tidy %>%
  dplyr::left_join(cluster_mapping, by = "doc_id") %>%
  dplyr::group_by(cluster, word) %>%
  dplyr::summarize(n = sum(count), .groups = "drop") %>%
  tidytext::bind_tf_idf(word, cluster, n) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(tf_idf, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

# ── CSV export ────────────────────────────────────────────────────────────────
write_csv(cluster_keywords, "step7_cluster_keywords.csv")
cat("  → step7_cluster_keywords.csv saved\n")

# Print top 5 per cluster for quick review
cat("\n  Top 5 keywords per cluster:\n")
cluster_keywords %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(tf_idf, n = 5) %>%
  dplyr::select(cluster, word, tf_idf) %>%
  print(n = 50)

# ── Plot 7: Keyword facet bar chart ──────────────────────────────────────────
p7 <- cluster_keywords %>%
  dplyr::mutate(word = tidytext::reorder_within(word, tf_idf, cluster)) %>%
  ggplot(aes(x = word, y = tf_idf, fill = cluster)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~cluster, scales = "free_y", ncol = 2) +
  coord_flip() +
  tidytext::scale_x_reordered() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Step 7: Top 10 Keywords per Cluster (TF-IDF)",
    subtitle = "Higher score = more distinctive to that cluster",
    x = "Term", y = "TF-IDF Score"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))
ggsave("step7_cluster_keywords.png", p7, width = 12, height = 14, dpi = 150)
cat("  → step7_cluster_keywords.png saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 8 — Final Disease Knowledge Map + Cross-Disease Heatmap
# ─────────────────────────────────────────────────────────────────────────────
cat("\n[Step 8] Building final knowledge map...\n")

# Disease labels — update these after reviewing your step7 keyword output
cluster_labels <- c(
  "1" = "Cancer Immunology / Epigenetics",
  "2" = "Oncology Trials / Chronic Disease Management",
  "3" = "Cognitive Neuroscience / Dementia Diagnostics",
  "4" = "Cardiovascular & Maternal Health Outcomes",
  "5" = "Vascular Surgery / Interventional",
  "6" = "Ferroptosis & Neuroinflammation",
  "7" = "Healthcare Delivery / Diabetes Care",
  "8" = "Neurodegeneration / PD-AD"
)

umap_df$disease_label <- cluster_labels[as.character(umap_df$cluster_kmeans)]

# Cluster centroids with size
centroids <- umap_df %>%
  dplyr::group_by(disease_label) %>%
  dplyr::summarise(
    UMAP1 = mean(UMAP1),
    UMAP2 = mean(UMAP2),
    n_docs = dplyr::n(),
    .groups = "drop"
  )

# ── CSV export ────────────────────────────────────────────────────────────────
final_map_df <- umap_df %>%
  dplyr::left_join(
    df %>% dplyr::select(doc_id, Abstract, abstract_length),
    by = "doc_id"
  )
write_csv(final_map_df, "step8_final_map_data.csv")
write_csv(centroids, "step8_cluster_centroids.csv")
cat("  → step8_final_map_data.csv saved\n")
cat("  → step8_cluster_centroids.csv saved\n")

# ── Plot 8a: Final Disease Knowledge Map ─────────────────────────────────────
# IMPROVEMENT: added bubble at centroid for cluster size context
p8a <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = disease_label)) +
  geom_point(alpha = 0.35, size = 0.7) +
  # Bubble shows relative cluster size
  geom_point(
    data = centroids,
    aes(size = n_docs),
    alpha = 0.12, show.legend = FALSE
  ) +
  ggrepel::geom_label_repel(
    data        = centroids,
    aes(label = paste0(disease_label, "\n(n=", n_docs, ")")),
    size        = 3,
    fontface    = "bold",
    box.padding = 0.8,
    label.size  = 0.2,
    show.legend = FALSE
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_size_continuous(range = c(8, 28)) +
  labs(
    title = "Disease Knowledge Map — PubMed Abstracts",
    subtitle = paste0(
      "Total abstracts: ", nrow(df),
      " | Clusters: ", k_chosen,
      " | KMeans Silhouette: ", avg_sil,
      " | DB Index: ", dbi_km
    ),
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Disease Domain"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, color = "grey50"),
    legend.position = "bottom",
    legend.text = element_text(size = 8)
  )
ggsave("step8_disease_knowledge_map.png", p8a, width = 12, height = 9, dpi = 200)
cat("  → step8_disease_knowledge_map.png saved\n")

# ── Plot 8b: Cross-Disease Cosine Similarity Heatmap (NEW) ───────────────────
# Shows which disease domains share the most vocabulary
cat("  Building cross-disease similarity heatmap...\n")

cluster_profiles <- lapply(levels(umap_df$cluster_kmeans), function(cl) {
  idx <- which(umap_df$cluster_kmeans == cl)
  colMeans(tfidf_matrix[idx, , drop = FALSE])
})
profile_matrix <- do.call(rbind, cluster_profiles)
rownames(profile_matrix) <- cluster_labels[levels(umap_df$cluster_kmeans)]

# Cosine similarity
norms <- sqrt(rowSums(profile_matrix^2))
normed <- profile_matrix / norms
cosine_sim <- normed %*% t(normed)
cosine_sim <- round(cosine_sim, 3)

write_csv(as.data.frame(cosine_sim), "step8_cosine_similarity.csv")

cosine_melt <- reshape2::melt(cosine_sim)
cosine_melt$Var1 <- factor(cosine_melt$Var1)
cosine_melt$Var2 <- factor(cosine_melt$Var2)

p8b <- ggplot(cosine_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 2.8) +
  scale_fill_gradient2(
    low      = "white",
    mid      = "#b8d4e8",
    high     = "#1a5f8a",
    midpoint = 0.3,
    name     = "Cosine\nSimilarity"
  ) +
  labs(
    title = "Step 8b: Cross-Disease Domain Similarity",
    subtitle = "Higher value = more shared biomedical vocabulary between clusters",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 8),
    plot.title   = element_text(face = "bold")
  )
ggsave("step8_similarity_heatmap.png", p8b, width = 10, height = 8, dpi = 150)
cat("  → step8_cosine_similarity.csv saved\n")
cat("  → step8_similarity_heatmap.png saved\n")

# ── Plot 8c: Cluster keyword wordcloud-style ranking (bonus) ─────────────────
top3_per_cluster <- cluster_keywords %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(tf_idf, n = 3) %>%
  dplyr::mutate(disease = cluster_labels[as.character(cluster)]) %>%
  dplyr::ungroup()

p8c <- ggplot(
  top3_per_cluster,
  aes(
    x = reorder(word, tf_idf),
    y = tf_idf, fill = as.factor(cluster)
  )
) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~disease, scales = "free", ncol = 2) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Step 8c: Top 3 Signature Terms per Disease Domain",
    subtitle = "These terms most uniquely identify each cluster",
    x = "Term", y = "TF-IDF Score"
  ) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 8))
ggsave("step8_signature_terms.png", p8c, width = 11, height = 10, dpi = 150)
cat("  → step8_signature_terms.png saved\n")


# ─────────────────────────────────────────────────────────────────────────────
# FINAL SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat(" PIPELINE COMPLETE — OUTPUT FILES\n")
cat(strrep("=", 60), "\n")

output_files <- data.frame(
  Step = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 8),
  Type = c(
    "CSV", "Plot", "Plot", "CSV", "Plot", "Plot", "CSV", "CSV", "Plot",
    "CSV", "Plot", "CSV", "CSV", "Plot", "Plot", "Plot", "CSV", "CSV", "Plot",
    "Plot", "CSV", "Plot", "CSV", "CSV", "CSV", "Plot", "Plot"
  ),
  File = c(
    "step1_cleaned_data.csv", "step1_abstract_length.png", "step1_data_quality.png",
    "step2_cleaned_corpus.csv", "step2_top_words.png", "step2_length_comparison.png",
    "step3_tfidf_summary.csv", "step3_vocabulary.csv", "step3_top_tfidf_terms.png",
    "step4_umap_coordinates.csv", "step4_umap_raw.png",
    "step5_cluster_assignments.csv", "step5_method_comparison.csv",
    "step5_elbow_plot.png", "step5_kmeans_umap.png", "step5_hdbscan_umap.png",
    "step6_silhouette_data.csv", "step6_evaluation_metrics.csv",
    "step6_silhouette_plot.png", "step6_cluster_sizes.png",
    "step7_cluster_keywords.csv", "step7_cluster_keywords.png",
    "step8_final_map_data.csv", "step8_cluster_centroids.csv",
    "step8_cosine_similarity.csv",
    "step8_disease_knowledge_map.png", "step8_similarity_heatmap.png"
  )
)
print(output_files, row.names = FALSE)

cat("\nKey metrics:\n")
cat("  Total abstracts processed:", nrow(df), "\n")
cat("  Vocabulary size:", nrow(vocab_df), "terms\n")
cat(
  "  KMeans Silhouette Score:", avg_sil,
  ifelse(avg_sil > 0.3, "(Good)", "(Needs tuning)"), "\n"
)
cat(
  "  Davies-Bouldin Index:", dbi_km,
  ifelse(dbi_km < 1.5, "(Acceptable)", "(Needs tuning)"), "\n"
)
cat("  HDBSCAN clusters found:", n_hdb_clusters, "\n")
