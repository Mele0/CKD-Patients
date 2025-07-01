import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Import UMAP and related functions
import umap.umap_ as umap
from sklearn.manifold import trustworthiness

from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from mpl_toolkits.mplot3d import Axes3D  # for 3D plotting
from math import sqrt

# --- Create results directory ---
results_dir = "results"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# --- Patch UMAP's internal Euclidean distance ---
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
        # For reference clustering, randomly assign points to clusters with the same sizes as original.
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
df = pd.read_csv("proteins_scaled_10years.csv", index_col=0)
X = df.select_dtypes(include=[np.number]).values
n_samples = X.shape[0]

# -------------------------
# 2. Single UMAP Run with User-Specified Parameters
# -------------------------
# Set the UMAP parameters you want to use:
umap_nn = 30        # e.g., number of neighbors
umap_md = 0.0       # e.g., min_dist
umap_metric = "cosine"     # e.g., cosine

reducer = umap.UMAP(n_components=2, n_neighbors=umap_nn, min_dist=umap_md, metric=umap_metric,  random_state=42)
best_embedding = reducer.fit_transform(X)
print(f"UMAP completed with n_neighbors={umap_nn}, min_dist={umap_md}, metric={umap_metric}")

# -------------------------
# 3. GMM Grid Search on Best UMAP Embedding
# -------------------------
gmm_results = []
n_components_list = range(2, 11)
cov_types = ["full", "tied", "diag", "spherical"]

print("=== GMM Grid Search on Best UMAP Embedding ===")
for n_components in n_components_list:
    for cov_type in cov_types:
        gmm = GaussianMixture(n_components=n_components, covariance_type=cov_type, random_state=42, n_init=5)
        gmm.fit(best_embedding)
        labels = gmm.predict(best_embedding)
        sil_umap = silhouette_score(best_embedding, labels, metric="euclidean")
        sil_raw = silhouette_score(X, labels, metric="euclidean")
        aic = gmm.aic(best_embedding)
        bic = gmm.bic(best_embedding)
        log_likelihood = gmm.score(best_embedding) * best_embedding.shape[0]
        
        # Compute PAC based on membership probabilities from GMM
        probs = gmm.predict_proba(best_embedding)
        pac_val = compute_pac(probs, threshold=0.9)
        
        # Compute co-membership index
        co_mem_val = compute_co_membership(labels)
        
        # Compute GAP statistic
        gap_val = compute_gap_statistic(best_embedding, labels, B=10)
        
        gmm_results.append({
            "n_components": n_components,
            "cov_type": cov_type,
            "silhouette_umap": sil_umap,
            "silhouette_raw": sil_raw,
            "AIC": aic,
            "BIC": bic,
            "log_likelihood": log_likelihood,
            "PAC": pac_val,
            "CoMembership": co_mem_val,
            "GAP": gap_val
        })
        print(f"n_components={n_components}, cov_type={cov_type}, sil(UMAP)={sil_umap:.3f}, sil(Raw)={sil_raw:.3f}, "
              f"AIC={aic:.1f}, BIC={bic:.1f}, LL={log_likelihood:.1f}, PAC={pac_val:.3f}, CoMem={co_mem_val:.3f}, GAP={gap_val:.3f}")

gmm_results_df = pd.DataFrame(gmm_results)
gmm_results_df.sort_values(by="silhouette_umap", ascending=False, inplace=True)
gmm_results_df.to_csv(os.path.join(results_dir, "gmm_all_configs_on_UMAP.csv"), index=False)
print("Saved all GMM configurations to gmm_all_configs_on_UMAP.csv")

top10_gmm = gmm_results_df.head(10)
top10_gmm.to_csv(os.path.join(results_dir, "gmm_top10_silhouette_on_UMAP.csv"), index=False)
print("Saved top 10 GMM configurations by silhouette to gmm_top10_silhouette_on_UMAP.csv")

# -------------------------
# 4. Plot Evaluation Metrics vs. Number of Clusters for Each Covariance Type
# -------------------------
def plot_metric(metric_col, ylabel, filename):
    plt.figure(figsize=(10,6))
    for cov in cov_types:
        subset = gmm_results_df[gmm_results_df["cov_type"] == cov].sort_values("n_components")
        plt.plot(subset["n_components"], subset[metric_col], marker="o", label=cov)
    plt.xlabel("Number of Clusters (n_components)")
    plt.ylabel(ylabel)
    plt.title(f"{ylabel} vs. Number of Clusters (GMM on UMAP)")
    plt.legend(title="Covariance Type")
    plt.savefig(os.path.join(results_dir, filename), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved plot: {filename}")

plot_metric("silhouette_umap", "Silhouette Score (UMAP)", "gmm_silhouette_vs_ncomponents.png")
plot_metric("AIC", "AIC", "gmm_AIC_vs_ncomponents.png")
plot_metric("BIC", "BIC", "gmm_BIC_vs_ncomponents.png")
plot_metric("log_likelihood", "Log-Likelihood", "gmm_loglikelihood_vs_ncomponents.png")

# -------------------------
# 5. Combined 3D Scatter Plot of Evaluation Metrics (Silhouette vs AIC vs Log-Likelihood)
# -------------------------
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection="3d")
colors = {"full": "red", "tied": "blue", "diag": "green", "spherical": "purple"}
for cov in cov_types:
    subset = gmm_results_df[gmm_results_df["cov_type"] == cov].sort_values("n_components")
    ax.scatter(subset["silhouette_umap"], subset["AIC"], subset["log_likelihood"],
               color=colors[cov], label=cov, s=50, alpha=0.8)
ax.set_xlabel("Silhouette (UMAP)")
ax.set_ylabel("AIC")
ax.set_zlabel("Log-Likelihood")
plt.title("GMM Evaluation: Silhouette vs AIC vs Log-Likelihood")
plt.legend(title="Covariance Type")
plt.savefig(os.path.join(results_dir, "gmm_evaluation_3d.png"), dpi=300, bbox_inches="tight")
plt.close()
print("Saved 3D evaluation metrics plot as gmm_evaluation_3d.png")

# -------------------------
# 6. Visualize Top 5 GMM Clustering Solutions on UMAP Embedding
# -------------------------
top5_gmm = top10_gmm.head(5)
for i, row in top5_gmm.iterrows():
    n_components = int(row["n_components"])
    cov_type = row["cov_type"]
    gmm = GaussianMixture(n_components=n_components, covariance_type=cov_type, random_state=42, n_init=5)
    gmm.fit(best_embedding)
    labels = gmm.predict(best_embedding)
    
    plt.figure(figsize=(8,6))

    sns.scatterplot(x=best_embedding[:,0],
                    y=best_embedding[:,1],
                    hue=labels,
                    palette="tab10",
                    s=60,
                    alpha=0.8)
    plt.title(f"GMM on UMAP: n_components={n_components}, cov={cov_type}\n"
              f"Sil(UMAP)={row['silhouette_umap']:.3f}, Sil(Raw)={row['silhouette_raw']:.3f}, AIC={row['AIC']:.1f}, "
              f"BIC={row['BIC']:.1f}, LL={row['log_likelihood']:.1f}")
    plt.xlabel("UMAP Dim 1")
    plt.ylabel("UMAP Dim 2")
    plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    filename = f"gmm_umap_top_{n_components}_{cov_type}.png"
    plt.savefig(os.path.join(results_dir, filename), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved cluster visualization: {filename}")

print("\nUMAP + GMM analysis complete. Check the 'results' directory for output files.")
