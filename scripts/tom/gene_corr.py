from multiprocessing import Pool
from itertools import chain
import pandas as pd
import numpy as np
import ateam.data.shiny as shiny
from ateam.analysis import ols
import sys

variables_long = [
    'SeuratMapping',
    'seurat_cluster',
    'L23_cell_depth',
    'max_path_length',
        'area'
]
variables = [
    'cluster',
    'cluster',
    'depth',
    'max_length',
    'membrane_area'
]

ephys_path = '/home/tom.chartrand/projects/ephys_analysis/data/human_mouse_ephys_all_0127.csv'
ephys_df = pd.read_csv(ephys_path, index_col=0)

md_path = "/home/tom.chartrand/projects/metadata/consolidated_clusters_and_metadata/human_IVSCC_excitatory_L23_consolidated_0130.csv"
md_df = pd.read_csv(md_path, index_col=0)
human_df = (md_df.join(ephys_df)
#             .drop(columns=['cluster'])
            .rename(columns=dict(zip(variables_long, variables)))
            .loc[lambda df: df.cluster.str.contains("Exc L2")]
           )

genes = pd.read_csv('/home/tom.chartrand/projects/ephys_analysis/data/gene_names_bad.csv', index_col=0).iloc[:,0].values
efeatures = list(ephys_df.columns) + ["depth"]
shiny_dir = shiny.shiny_directory('human')
join_on = 'sample_id'

def calc(gene):
    genes_df = (pd.read_feather(shiny_dir + '/data.feather', columns=[gene, join_on])
                .set_index(join_on)
                .apply(lambda x: np.log2(x+1))
               )
    data = human_df.join(genes_df, on='sample_id')
    if data.groupby('cluster')[gene].apply(lambda x: sum(x>1)).min() < 5:
        return []
    results = []
    for efeature in efeatures:
        results.append(ols.anova_all(data, efeature, gene, f2='cluster', cov_type='HC3'))
    return results

# out = map(calc, genes[:20])
# n=1

n = sys.argv[1]
step = 1000
start = int(n)*step
end = min(start+step, len(genes))
pool = Pool()
out = pool.map(calc, genes[start:end])
metrics = pd.DataFrame.from_records(chain(*out)).drop(columns=["cluster"])
metrics.to_csv(f'/home/tom.chartrand/projects/ephys_analysis/data/genes/gene_ephys_new{n}.csv')

# [3,10,1,19,35,29,30,40,36,46]