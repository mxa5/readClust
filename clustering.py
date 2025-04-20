import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from joblib import Parallel, delayed
from tqdm import tqdm

dists = []
sample_rep = []

def cluster_chromosome_kdtree(df_chrom, all_samples, dist, rep):

    df_chrom = df_chrom.sort_values(["start_pos", "end_pos"]).reset_index(drop=True)
    n = len(df_chrom)
    assigned = np.zeros(n, dtype=bool)
    cluster_ids = np.full(n, -1, dtype=int)
    
    # Build a KD-tree on (start_pos, end_pos)
    coords = df_chrom[["start_pos", "end_pos"]].values
    tree = cKDTree(coords)
    
    current_cluster_id = 1
    required_sample_count = max(1, int(rep * len(all_samples)))
    
    for i in range(n):
        if assigned[i]:
            continue
            
        lower_start = coords[i, 0]
        lower_end = coords[i, 1]
        center = np.array([lower_start + dist, lower_end + dist])
        candidate_indices = tree.query_ball_point(center, dist, p=np.inf)
        candidate_indices = [j for j in candidate_indices if not assigned[j]]
        
        if len(candidate_indices) == 0:
            continue
        
        
        candidate_samples = set(df_chrom.loc[candidate_indices, "sample"])
        
        if len(candidate_samples) >= required_sample_count:
            
            for j in candidate_indices:
                cluster_ids[j] = current_cluster_id
                assigned[j] = True
            current_cluster_id += 1

    df_chrom["cluster_id"] = cluster_ids
    return df_chrom

def main():
    
    for dist in dists:
        for rep in sample_rep:
            
            required_cols = {"read_id", "start_pos", "end_pos", "chrom", "sample"}
            if not required_cols.issubset(df.columns):
                raise ValueError("Input data must contain columns: " + ", ".join(required_cols))
            
            
            all_samples = set(df["sample"].unique())
            
            
            chrom_groups = list(df.groupby("chrom"))
            
            results = Parallel(n_jobs=64)(
                delayed(cluster_chromosome_kdtree)(group.copy(), all_samples, dist, rep)
                for chrom, group in tqdm(chrom_groups, desc="Clustering chromosomes", total=len(chrom_groups))
            )
            
            
            clustered_df = pd.concat(results, ignore_index=True)

            
            clustered_df = clustered_df[clustered_df['cluster_id'] != -1]

            
            chrom_cluster_reads = clustered_df.groupby("chrom").size()
            
            print("Clustering complete.")

if __name__ == "__main__":
    main()
