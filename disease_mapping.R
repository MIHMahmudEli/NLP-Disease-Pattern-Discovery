#Step 0 — Install all packages
install.packages(c(
  "readr", "dplyr", "stringr",
  "tm", "SnowballC", "textstem",
  "Matrix", "uwot",
  "stats",           # kmeans (built-in)
  "dbscan",          # HDBSCAN
  "cluster",         # silhouette score
  "clusterSim",      # Davies-Bouldin Index
  "tidytext",        # tidy keyword extraction
  "ggplot2", "ggrepel"
))

#Step 1 — Load data
library(readr)
library(dplyr)

df <- read_csv("pubmed_dataset.csv")
cat("Rows:", nrow(df), "\n")   # Should print 8934
cat("Columns:", names(df), "\n")

# Remove empty abstracts
df <- df %>% filter(!is.na(Abstract) & nchar(Abstract) > 50)
cat("After cleaning:", nrow(df), "\n")


#Step 2 — Text preprocessing
library(tm)
library(textstem)   # for lemmatization (better than stemming for biomedical text)

# Build a corpus from the Abstract column
corpus <- VCorpus(VectorSource(df$Abstract))

# Custom biomedical stopwords (supplement the standard list)
bio_stopwords <- c(
  stopwords("en"),
  "study", "result", "patient", "method", "analysis",
  "however", "also", "may", "well", "used", "using",
  "associated", "significantly", "found", "showed"
)

corpus_clean <- corpus %>%
  tm_map(content_transformer(tolower)) %>%
  tm_map(removePunctuation) %>%
  tm_map(removeNumbers) %>%
  tm_map(removeWords, bio_stopwords) %>%
  tm_map(stripWhitespace) %>%
  tm_map(content_transformer(lemmatize_strings))  # lemmatize (not stem)

# Inspect one abstract after cleaning
inspect(corpus_clean[[1]])


#Step 3 — TF-IDF feature extraction
library(Matrix)

# Build Document-Term Matrix, then apply TF-IDF weighting
dtm <- DocumentTermMatrix(
  corpus_clean,
  control = list(
    weighting  = weightTfIdf,
    bounds     = list(global = c(10, Inf)),  # term must appear in ≥10 docs
    wordLengths = c(3, 25)                   # filter very short/long tokens
  )
)

cat("DTM dimensions:", dim(dtm), "\n")  # docs × terms

# Remove sparse terms (keeps ~top features)
dtm_reduced <- removeSparseTerms(dtm, sparse = 0.995)
cat("After sparsity reduction:", dim(dtm_reduced), "\n")

# Convert to a regular matrix for UMAP
tfidf_matrix <- as.matrix(dtm_reduced)



#Step 4 — UMAP dimensionality reduction
library(uwot)

set.seed(42)

umap_result <- umap(
  tfidf_matrix,
  n_components  = 2,      # 2D for visualization
  n_neighbors   = 15,     # local neighborhood size
  min_dist      = 0.1,    # controls cluster tightness
  metric        = "cosine",  # best for text/TF-IDF data
  n_threads     = 4,
  verbose       = TRUE
)

# umap_result is an N×2 matrix
umap_df <- data.frame(
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2]
)


#Step 5a — KMeans clustering
# First, find the optimal k using the elbow method
wss <- sapply(2:15, function(k) {
  kmeans(umap_df, centers = k, nstart = 20, iter.max = 100)$tot.withinss
})

plot(2:15, wss, type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster SS",
     main = "Elbow Method for Optimal k")

# Apply KMeans with chosen k (e.g. k=8 for 8 disease categories)
set.seed(42)
k_chosen <- 8
km_result <- kmeans(umap_df, centers = k_chosen, nstart = 25, iter.max = 200)

umap_df$cluster_kmeans <- as.factor(km_result$cluster)


#Step 5b — HDBSCAN clustering (alternative)
library(dbscan)

hdb_result <- hdbscan(
  umap_df[, c("UMAP1", "UMAP2")],
  minPts = 50   # min documents to form a cluster; tune this
)

cat("HDBSCAN clusters found:", max(hdb_result$cluster), "\n")
cat("Noise points (cluster=0):", sum(hdb_result$cluster == 0), "\n")

umap_df$cluster_hdbscan <- as.factor(hdb_result$cluster)


#Step 6 — Evaluation metrics
library(cluster)
library(clusterSim)

# --- Silhouette Score (higher is better, range -1 to 1) ---
sil_km <- silhouette(km_result$cluster, dist(umap_df[, c("UMAP1", "UMAP2")]))
cat("KMeans Silhouette Score:", mean(sil_km[, 3]), "\n")

# For HDBSCAN, exclude noise points (cluster = 0)
valid_idx <- hdb_result$cluster != 0
if (sum(valid_idx) > 1) {
  sil_hdb <- silhouette(
    hdb_result$cluster[valid_idx],
    dist(umap_df[valid_idx, c("UMAP1", "UMAP2")])
  )
  cat("HDBSCAN Silhouette Score:", mean(sil_hdb[, 3]), "\n")
}

# --- Davies-Bouldin Index (lower is better) ---
dbi_km <- index.DB(
  umap_df[, c("UMAP1", "UMAP2")],
  km_result$cluster
)$DB
cat("KMeans Davies-Bouldin Index:", dbi_km, "\n")





#Step 7 — Keyword extraction per cluster
library(tidytext)
library(dplyr)

# Attach cluster labels
df$cluster <- umap_df$cluster_kmeans

cluster_keywords <- df %>%
  dplyr::select(Abstract, cluster) %>%          # <- explicit namespace
  unnest_tokens(word, Abstract) %>%
  filter(!word %in% stopwords("en"), nchar(word) > 3) %>%
  count(cluster, word, sort = TRUE) %>%
  bind_tf_idf(word, cluster, n) %>%
  group_by(cluster) %>%
  slice_max(tf_idf, n = 10) %>%
  ungroup()


# Print top keywords for each cluster
cluster_keywords %>%
  arrange(cluster, desc(tf_idf)) %>%
  print(n = 80)




#Step 8 — UMAP visualization
library(ggplot2)
library(ggrepel)

# Disease label mapping (update after inspecting cluster keywords above)
cluster_labels <- c(
  "1" = "Diabetes / Metabolic",
  "2" = "Oncology / Cancer",
  "3" = "Cardiovascular",
  "4" = "Neurology",
  "5" = "Infectious Disease",
  "6" = "Kidney / Renal",
  "7" = "Obesity / Nutrition",
  "8" = "Immunology"
)

umap_df$disease_label <- cluster_labels[as.character(umap_df$cluster_kmeans)]

# Compute cluster centroids for label placement
centroids <- umap_df %>%
  group_by(disease_label) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = disease_label)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_label_repel(
    data      = centroids,
    aes(label = disease_label),
    size      = 3,
    fontface  = "bold",
    box.padding = 0.5,
    show.legend = FALSE
  ) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title    = "Disease knowledge map — PubMed abstracts",
    subtitle = paste0(nrow(df), " abstracts, k=", k_chosen, " clusters"),
    x        = "UMAP dimension 1",
    y        = "UMAP dimension 2",
    color    = "Disease cluster"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("disease_knowledge_map.png", width = 10, height = 8, dpi = 300)