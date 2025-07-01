import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import umap.umap_ as umap
from sklearn.manifold import trustworthiness
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from math import sqrt

# --- User-Specified Parameters ---
# UMAP parameters:
umap_n_neighbors = 30   # e.g., 15
umap_min_dist = 0.0     # e.g., 0.1
umap_n_components = 2   # For visualization
umap_metric = "cosine"  # e.g. cosine

# GMM parameters:
gmm_n_components = 3    # Number of clusters/components
gmm_covariance_type = "spherical"  # Options: "full", "tied", "diag", "spherical"
gmm_n_init = 5          # Number of initializations

# File paths:
input_csv = "proteins_scaled_10years.csv"            # Replace with your input CSV filename
output_metrics_csv = "gmm_umap_metrics.csv"  # CSV file with evaluation metrics
output_plot = "umap_gmm_clustering.svg"       # Plot output (you can change this to pdf, png, svg)

# --- Create results directory if needed ---
results_dir = "results"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# --- Optional: Patch UMAP's internal Euclidean distance (for numerical safety) ---
import umap.distances as ud
def safe_euclidean(x, y):
    result = 0.0
    for i in range(len(x)):
        diff = x[i] - y[i]
        result += diff * diff
    return sqrt(max(result, 0))
ud.euclidean = safe_euclidean

# --- Define additional evaluation metric functions ---

def compute_pac(probabilities, threshold=0.9):
    """
    Proportion of Ambiguous Clustering (PAC):
    Proportion of points for which the maximum membership probability is less than the threshold.
    """
    max_probs = probabilities.max(axis=1)
    pac = np.mean(max_probs < threshold)
    return pac

def compute_co_membership(labels):
    """
    Compute a simple co-membership index:
    Average proportion of pairs of points that are assigned to the same cluster.
    """
    n = len(labels)
    same = 0
    total = 0
    for i in range(n):
        for j in range(i+1, n):
            total += 1
            if labels[i] == labels[j]:
                same += 1
    return same / total if total > 0 else 0

def compute_gap_statistic(X, labels, B=10):
    """
    Compute a basic GAP statistic.
    Calculate the within-cluster dispersion for X and compare it to B reference datasets 
    generated uniformly within the bounds of X.
    """
    unique_labels = np.unique(labels)
    Wk = 0.0
    for label in unique_labels:
        cluster_points = X[labels == label]
        if len(cluster_points) > 0:
            centroid = cluster_points.mean(axis=0)
            Wk += np.sum(np.linalg.norm(cluster_points - centroid, axis=1) ** 2)
    Wk = Wk / len(X)
    
    # Generate B reference datasets and compute dispersion
    mins = X.min(axis=0)
    maxs = X.max(axis=0)
    Wk_refs = []
    for b in range(B):
        X_ref = np.random.uniform(low=mins, high=maxs, size=X.shape)
        ref_labels = np.random.choice(unique_labels, size=len(X))
        Wk_ref = 0.0
        for label in unique_labels:
            ref_cluster = X_ref[ref_labels == label]
            if len(ref_cluster) > 0:
                centroid = ref_cluster.mean(axis=0)
                Wk_ref += np.sum(np.linalg.norm(ref_cluster - centroid, axis=1) ** 2)
        Wk_refs.append(Wk_ref / len(X_ref))
    Wk_ref_mean = np.mean(Wk_refs)
    gap = np.log(Wk_ref_mean) - np.log(Wk)
    return gap

# -------------------------
# 1. Load the Data
# -------------------------
# Assume your CSV has a header and the first column (e.g., participant ID) as index.
df = pd.read_csv(input_csv, index_col=0)
# Assume the data are numeric (or select the numeric columns)
X = df.select_dtypes(include=[np.number]).values
n_samples = X.shape[0]

# -------------------------
# 2. Run UMAP with user-specified parameters
# -------------------------
umap_reducer = umap.UMAP(n_components=umap_n_components,
                         n_neighbors=umap_n_neighbors,
                         min_dist=umap_min_dist,
                         metric=umap_metric,
                         random_state=42)
umap_embedding = umap_reducer.fit_transform(X)
umap_tw = trustworthiness(X, umap_embedding, n_neighbors=umap_n_neighbors)
print(f"UMAP done: n_neighbors={umap_n_neighbors}, min_dist={umap_min_dist}, trustworthiness={umap_tw:.4f}")

# -------------------------
# 3. Run GMM on the UMAP embedding
# -------------------------
gmm = GaussianMixture(n_components=gmm_n_components,
                      covariance_type=gmm_covariance_type,
                      random_state=42,
                      n_init=gmm_n_init)
gmm.fit(umap_embedding)
gmm_labels = gmm.predict(umap_embedding)
sil_umap = silhouette_score(umap_embedding, gmm_labels, metric="euclidean")
sil_raw = silhouette_score(X, gmm_labels, metric="euclidean")
aic = gmm.aic(umap_embedding)
bic = gmm.bic(umap_embedding)
log_likelihood = gmm.score(umap_embedding) * n_samples

# Compute additional metrics
probs = gmm.predict_proba(umap_embedding)
pac_val = compute_pac(probs, threshold=0.9)
co_mem_val = compute_co_membership(gmm_labels)
gap_val = compute_gap_statistic(umap_embedding, gmm_labels, B=10)

# -------------------------
# 4. Save evaluation metrics to CSV
# -------------------------
metrics = {
    "umap_n_neighbors": umap_n_neighbors,
    "umap_min_dist": umap_min_dist,
    "umap_trustworthiness": umap_tw,
    "gmm_n_components": gmm_n_components,
    "gmm_covariance_type": gmm_covariance_type,
    "silhouette_umap": sil_umap,
    "silhouette_raw": sil_raw,
    "AIC": aic,
    "BIC": bic,
    "log_likelihood": log_likelihood,
    "PAC": pac_val,
    "CoMembership": co_mem_val,
    "GAP": gap_val
}
metrics_df = pd.DataFrame([metrics])
metrics_df.to_csv(os.path.join(results_dir, output_metrics_csv), index=False)
print(f"Saved evaluation metrics to {os.path.join(results_dir, output_metrics_csv)}")

# -------------------------
# 5. Plot UMAP embedding with GMM clustering
# -------------------------
plt.figure(figsize=(8,6))
sns.scatterplot(x=umap_embedding[:,0], y=umap_embedding[:,1],
                hue=gmm_labels, palette="tab10", s=60, alpha=0.8)
plt.title(f"UMAP + GMM Clustering\n(n_neighbors={umap_n_neighbors}, min_dist={umap_min_dist}, "
          f"trustworthiness={umap_tw:.2f}, silhouette={sil_umap:.2f})")
plt.xlabel("UMAP Dimension 1")
plt.ylabel("UMAP Dimension 2")
plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig(os.path.join(results_dir, output_plot), dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved UMAP + GMM clustering plot to {os.path.join(results_dir, output_plot)}")

# -------------------------
# 6. Saving .csv file for further downstream analysis of clusters
# -------------------------

df_with_clusters = df.copy()
df_with_clusters["Cluster"] = gmm_labels

# Reset index so that the original index becomes the first column
df_with_clusters.reset_index(inplace=True)

# Save the resulting DataFrame to a CSV file (index is now a column)
output_with_clusters = "input_with_clusters.csv"
df_with_clusters.to_csv(os.path.join(results_dir, output_with_clusters), index=False)

print(f"Saved original data with cluster labels to {os.path.join(results_dir, output_with_clusters)}")
