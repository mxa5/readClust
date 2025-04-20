from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

mw_results = []

# Group the data by unique cluster_id
for cluster, df_cluster in clusters_merged.groupby('uniq_id'):
  
    control = df_cluster[df_cluster['type'] == 'control']['% methylated']
    disease = df_cluster[df_cluster['type'] == 'disease']['% methylated']
    
    
    if len(control) > 0 and len(disease) > 0:
      
        # Perform a two-sided Mann–Whitney U test
        stat, p_value = mannwhitneyu(control, disease, alternative='two-sided')
    else:
        stat, p_value = None, None
        
    mw_results.append({
        'cluster_id': cluster,
        'statistic': stat,
        'p_value': p_value,
        'control_mean': control.mean(),
        'disease_mean': disease.mean()
    })

df_mw_clean = df_mw.dropna(subset=['p_value'])

rejected, pvals_corrected, _, _ = multipletests(df_mw_clean['p_value'], method='fdr_bh')

df_mw_clean['p_adj'] = pvals_corrected
df_mw_clean['significant'] = rejected


df_mw = df_mw.merge(df_mw_clean[['cluster_id', 'p_adj', 'significant']], on='cluster_id', how='left')
