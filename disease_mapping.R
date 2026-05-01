# Disease-Specific Knowledge Mapping (Unsupervised Clustering)
# Project: Disease-Specific Knowledge Mapping
# Author: Antigravity AI
# Date: 2026-05-01

# 1. Load Libraries -------------------------------------------------------
cat("\n--- STEP 1: LOADING LIBRARIES ---\n")
required_packages <- c("tidyverse", "tidytext", "tm", "uwot", "dbscan", "cluster", "fpc", "wordcloud", "SnowballC", "Matrix")

safe_install <- function(pkgs) {
  lock_dir <- file.path(.libPaths()[1], "00LOCK")
  if (dir.exists(lock_dir)) {
    cat("Removing existing lock directory:", lock_dir, "\n")
    unlink(lock_dir, recursive = TRUE, force = TRUE)
  }
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs)) {
    cat("Installing missing packages:", paste(new_pkgs, collapse = ", "), "\n")
    install.packages(new_pkgs, dependencies = TRUE, INSTALL_opts = "--no-lock")
  }
}

safe_install(required_packages)

library(tidyverse)
library(tidytext)
library(tm)
library(uwot)
library(dbscan)
library(cluster)
library(fpc)
library(wordcloud)
library(SnowballC)
library(Matrix)
cat("Libraries loaded successfully.\n")

# 2. Load Dataset ---------------------------------------------------------
cat("\n--- STEP 2: LOADING DATASET ---\n")
data_path <- "pubmed_dataset.csv"
if (!file.exists(data_path)) stop("Dataset 'pubmed_dataset.csv' not found!")

df <- read.csv(data_path, stringsAsFactors = FALSE)
cat("Dataset loaded. Total rows:", nrow(df), "\n")

# Use only Abstract as per user requirement
df <- df %>%
  mutate(text = Abstract) %>%
  mutate(doc_id = row_number())

cat("Data structure preview (First 3 rows):\n")
print(head(df[, c("doc_id", "Disease", "text")], 3))

# 3. Preprocessing --------------------------------------------------------
cat("\n--- STEP 3: PREPROCESSING TEXT ---\n")

clean_text <- function(text) {
  text <- tolower(text)
  text <- gsub("<.*?>", "", text) # Remove HTML tags
  text <- removePunctuation(text)
  text <- removeNumbers(text)
  text <- removeWords(text, stopwords("en"))
  text <- stripWhitespace(text)
  return(text)
}

cat("Cleaning text and removing HTML tags...\n")
df$cleaned_text <- clean_text(df$text)

cat("Tokenizing and Stemming...\n")
tidy_docs <- df %>%
  unnest_tokens(word, cleaned_text) %>%
  mutate(word = wordStem(word)) %>%
  filter(!word %in% stop_words$word)

cat("Tokenized data preview:\n")
print(head(tidy_docs, 10))

# 4. TF-IDF Feature Extraction -------------------------------------------
cat("\n--- STEP 4: TF-IDF FEATURE EXTRACTION ---\n")

cat("Filtering rare words (appearing in < 10 documents)...\n")
word_counts <- tidy_docs %>%
  count(word) %>%
  filter(n >= 10)

cat("Unique words retained:", nrow(word_counts), "\n")

tf_idf_matrix <- tidy_docs %>%
  filter(word %in% word_counts$word) %>%
  count(doc_id, word) %>%
  bind_tf_idf(word, doc_id, n)

cat("Generating sparse matrix...\n")
sparse_matrix <- tf_idf_matrix %>%
  cast_sparse(doc_id, word, tf_idf)

# Remove empty documents
row_sums <- Matrix::rowSums(sparse_matrix)
valid_docs <- which(row_sums > 0)
sparse_matrix_clean <- sparse_matrix[valid_docs, ]

cat("Final matrix dimensions for UMAP:", nrow(sparse_matrix_clean), "x", ncol(sparse_matrix_clean), "\n")

# 5. Dimensionality Reduction (UMAP) -------------------------------------
cat("\n--- STEP 5: UMAP DIMENSIONALITY REDUCTION ---\n")
cat("Reducing dimensions to 2D using Cosine Similarity (this may take a moment)...\n")

umap_results <- umap(as.matrix(sparse_matrix_clean), 
                     n_neighbors = 15, 
                     min_dist = 0.1, 
                     metric = "cosine",
                     n_threads = parallel::detectCores() - 1)

umap_df <- as.data.frame(umap_results)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$doc_id <- as.integer(rownames(sparse_matrix_clean))

cat("UMAP coordinates preview:\n")
print(head(umap_df, 5))

# 6. Clustering -----------------------------------------------------------
cat("\n--- STEP 6: CLUSTERING (K-MEANS & HDBSCAN) ---\n")

set.seed(42)
k_value <- 4 # Target disease categories
cat("Running K-means with K =", k_value, "...\n")
kmeans_res <- kmeans(umap_df[, c("UMAP1", "UMAP2")], centers = k_value, nstart = 25)
umap_df$cluster_kmeans <- as.factor(kmeans_res$cluster)

cat("Running HDBSCAN (Density-based)...\n")
hdbscan_res <- hdbscan(umap_df[, c("UMAP1", "UMAP2")], minPts = 10)
umap_df$cluster_hdbscan <- as.factor(hdbscan_res$cluster)

cat("Clustered data preview:\n")
print(head(umap_df, 5))

# 7. Evaluation Metrics ---------------------------------------------------
cat("\n--- STEP 7: EVALUATION METRICS ---\n")

cat("Calculating Silhouette Score...\n")
sil <- silhouette(kmeans_res$cluster, dist(umap_df[, c("UMAP1", "UMAP2")]))
avg_sil <- mean(sil[, 3])
cat("Average Silhouette Width (K-means):", round(avg_sil, 4), "\n")

cat("Calculating Davies-Bouldin Index...\n")
db_index <- cluster.stats(dist(umap_df[, c("UMAP1", "UMAP2")]), kmeans_res$cluster)$db
cat("Davies-Bouldin Index (K-means):", round(db_index, 4), "\n")

# 8. Keyword Extraction per Cluster --------------------------------------
cat("\n--- STEP 8: KEYWORD EXTRACTION ---\n")

cluster_info <- umap_df %>% select(doc_id, cluster_kmeans)
tidy_clusters <- tidy_docs %>%
  inner_join(cluster_info, by = "doc_id")

cat("Top 10 keywords per cluster:\n")
top_keywords <- tidy_clusters %>%
  group_by(cluster_kmeans, word) %>%
  count() %>%
  arrange(cluster_kmeans, desc(n)) %>%
  group_by(cluster_kmeans) %>%
  slice_max(n, n = 10)

print(top_keywords)

# 9. Visualization --------------------------------------------------------
cat("\n--- STEP 9: GENERATING VISUALIZATIONS ---\n")

cat("Creating UMAP cluster plot...\n")
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster_kmeans)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Disease-Specific Knowledge Mapping (UMAP + K-means)",
       subtitle = paste("Clusters identified:", k_value),
       color = "Cluster")

ggsave("umap_clusters_plot.png", umap_plot, width = 10, height = 7)

if ("Disease" %in% colnames(df)) {
  cat("Creating Disease Mapping plot...\n")
  disease_mapping <- umap_df %>%
    inner_join(df %>% select(doc_id, Disease), by = "doc_id") %>%
    group_by(cluster_kmeans, Disease) %>%
    count() %>%
    arrange(cluster_kmeans, desc(n))
  
  mapping_plot <- ggplot(disease_mapping, aes(x = cluster_kmeans, y = n, fill = Disease)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = "Disease Mapping to Identified Clusters",
         x = "Cluster", y = "Document Count")
  
  ggsave("disease_mapping_plot.png", mapping_plot, width = 10, height = 7)
}

cat("\n--- PROJECT COMPLETED SUCCESSFULLY ---\n")
cat("Outputs saved:\n")
cat("1. umap_clusters_plot.png\n")
cat("2. disease_mapping_plot.png\n\n")
