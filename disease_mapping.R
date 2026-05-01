# Disease-Specific Knowledge Mapping (Unsupervised Clustering)
# Project Idea 4 Implementation

# Step 0 — Install/Load Libraries
required_packages <- c("readr", "dplyr", "stringr", "tm", "SnowballC", "textstem", 
                       "Matrix", "uwot", "dbscan", "cluster", "clusterSim", 
                       "tidytext", "ggplot2", "ggrepel")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos='http://cran.us.r-project.org')

library(readr)
library(dplyr)
library(ggplot2)

# Step 1 — Load data
cat("Step 1: Loading Data...\n")
df <- read_csv("pubmed_dataset.csv")
colnames(df)[1] <- "Abstract" # Ensure naming consistency

# Visualization: Abstract Length Distribution
length_plot <- ggplot(df, aes(x = nchar(Abstract))) +
  geom_histogram(fill = "#69b3a2", color = "white", bins = 50) +
  labs(title = "Step 1: Distribution of Abstract Lengths", x = "Character Count", y = "Frequency") +
  theme_minimal()
ggsave("step1_length_distribution.png", length_plot)

# Remove empty/short abstracts
df <- df %>% filter(!is.na(Abstract) & nchar(Abstract) > 50)
write_csv(df, "step1_raw_cleaned_data.csv")


# Step 2 — Text preprocessing
cat("Step 2: Preprocessing Text...\n")
library(tm)
library(textstem)

corpus <- VCorpus(VectorSource(df$Abstract))
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
  tm_map(content_transformer(lemmatize_strings))

# Export Cleaned Corpus as CSV
cleaned_df <- data.frame(doc_id = 1:length(corpus_clean), 
                         text = sapply(corpus_clean, as.character), 
                         stringsAsFactors = FALSE)
write_csv(cleaned_df, "step2_cleaned_corpus.csv")

# Visualization: Word Frequencies
library(tidytext)
word_freq <- cleaned_df %>%
  unnest_tokens(word, text) %>%
  count(word, sort = TRUE) %>%
  head(30)

freq_plot <- ggplot(word_freq, aes(x = reorder(word, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Step 2: Top 30 Words After Preprocessing", x = "Word", y = "Frequency") +
  theme_minimal()
ggsave("step2_word_frequencies.png", freq_plot)


# Step 3 — TF-IDF feature extraction
cat("Step 3: Feature Extraction (TF-IDF)...\n")
library(Matrix)

dtm <- DocumentTermMatrix(
  corpus_clean,
  control = list(
    weighting  = weightTfIdf,
    bounds     = list(global = c(10, Inf)),
    wordLengths = c(3, 25)
  )
)

dtm_reduced <- removeSparseTerms(dtm, sparse = 0.995)
tfidf_matrix <- as.matrix(dtm_reduced)
write_csv(as.data.frame(tfidf_matrix), "step3_tfidf_matrix.csv")


# Step 4 — UMAP dimensionality reduction
cat("Step 4: Dimensionality Reduction (UMAP)...\n")
library(uwot)
set.seed(42)

umap_result <- umap(
  tfidf_matrix,
  n_components  = 2,
  n_neighbors   = 15,
  min_dist      = 0.1,
  metric        = "cosine",
  n_threads     = 4
)

umap_df <- data.frame(
  UMAP1 = umap_result[, 1],
  UMAP2 = umap_result[, 2]
)
write_csv(umap_df, "step4_umap_coordinates.csv")

# Visualization: Raw UMAP
umap_raw_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.4, color = "grey50") +
  labs(title = "Step 4: UMAP Projection (Before Clustering)", subtitle = "Natural density of PubMed abstracts") +
  theme_minimal()
ggsave("step4_umap_raw.png", umap_raw_plot)


# Step 5a — KMeans clustering
cat("Step 5: Clustering...\n")
wss <- sapply(2:15, function(k) {
  kmeans(umap_df, centers = k, nstart = 20)$tot.withinss
})

# Save Elbow Plot
png("step5_elbow_plot.png")
plot(2:15, wss, type = "b", pch = 19, xlab = "Number of Clusters (k)", ylab = "WSS", main = "Elbow Method")
dev.off()

set.seed(42)
k_chosen <- 8
km_result <- kmeans(umap_df, centers = k_chosen, nstart = 25)
umap_df$cluster_kmeans <- as.factor(km_result$cluster)
write_csv(umap_df, "step5_kmeans_results.csv")

# Step 5b — HDBSCAN clustering
library(dbscan)
hdb_result <- hdbscan(umap_df[, c("UMAP1", "UMAP2")], minPts = 50)
umap_df$cluster_hdbscan <- as.factor(hdb_result$cluster)
write_csv(umap_df, "step5_hdbscan_results.csv")


# Step 6 — Evaluation metrics
cat("Step 6: Evaluation...\n")
library(cluster)
library(clusterSim)

# KMeans Metrics
sil_km <- silhouette(km_result$cluster, dist(umap_df[, c("UMAP1", "UMAP2")]))
write_csv(as.data.frame(sil_km[]), "step6_silhouette_data.csv")

# Visualization: Silhouette Plot
png("step6_silhouette_plot.png", width = 800, height = 600)
plot(sil_km, main = "Step 6: KMeans Silhouette Analysis", col = 1:k_chosen, border = NA)
dev.off()

dbi_km <- index.DB(umap_df[, c("UMAP1", "UMAP2")], km_result$cluster)$DB
write_lines(dbi_km, "step6_davies_bouldin.txt")


# Step 7 — Keyword extraction per cluster
cat("Step 7: Keyword Extraction...\n")
dtm_tidy <- tidy(dtm_reduced) %>%
  rename(doc_id = document, word = term) %>%
  mutate(doc_id = as.integer(doc_id))

cluster_mapping <- data.frame(doc_id = 1:nrow(df), cluster = umap_df$cluster_kmeans)

cluster_keywords <- dtm_tidy %>%
  left_join(cluster_mapping, by = "doc_id") %>%
  group_by(cluster, word) %>%
  summarize(n = sum(count), .groups = 'drop') %>%
  bind_tf_idf(word, cluster, n) %>%
  group_by(cluster) %>%
  slice_max(tf_idf, n = 10, with_ties = FALSE) %>%
  ungroup()

write_csv(cluster_keywords, "step7_cluster_keywords.csv")

# Visualization: Keyword Facets
keyword_viz <- ggplot(cluster_keywords, aes(x = reorder_within(word, tf_idf, cluster), y = tf_idf, fill = cluster)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~cluster, scales = "free") +
  coord_flip() +
  scale_x_reordered() +
  labs(title = "Step 7: Top Cluster Keywords (TF-IDF)", x = "Terms", y = "Score") +
  theme_minimal()
ggsave("step7_keyword_facets.png", keyword_viz, width = 12, height = 8)


# Step 8 — Final Knowledge Map
cat("Step 8: Final Visualization...\n")
library(ggrepel)

cluster_labels <- c(
  "1" = "Molecular Pathologies / Ferroptosis",
  "2" = "Neuro-Rehabilitation / Cognitive",
  "3" = "Genomic Oncology / HCC",
  "4" = "Neurodegeneration / PD-AD",
  "5" = "Bio-Molecular Signaling",
  "6" = "Metabolic / Obesity & Diabetes",
  "7" = "Vascular Surgery / Interventional",
  "8" = "Clinical Outcomes / Stroke-CVD"
)

umap_df$disease_label <- cluster_labels[as.character(umap_df$cluster_kmeans)]
centroids <- umap_df %>% group_by(disease_label) %>% summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

final_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = disease_label)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_label_repel(data = centroids, aes(label = disease_label), size = 3.5, fontface = "bold") +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Final Disease Knowledge Map — PubMed Abstracts", 
       subtitle = paste0("Total Abstracts: ", nrow(df), " | Clusters: ", k_chosen),
       caption = paste0("KMeans Silhouette Score: ", round(mean(sil_km[, 3]), 3)),
       x = "UMAP Dimension 1", y = "UMAP Dimension 2", color = "Disease Domain") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave("final_disease_knowledge_map.png", final_plot, width = 11, height = 8)

cat("\n[Pipeline Complete] All step-by-step CSVs and Graphs generated successfully.\n")