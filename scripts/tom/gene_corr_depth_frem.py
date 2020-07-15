from multiprocessing import Pool
from itertools import chain
import pandas as pd
import numpy as np
import ateam.data.shiny as shiny
from ateam.analysis import ols
import sys
from scipy.stats import pearsonr
from pingouin import partial_corr

import os.path
n = sys.argv[1]
out_path = f'/home/tom.chartrand/projects/data/genes/frem_gene_depth_ephys_{n}.csv'
if os.path.exists(out_path):
    quit()

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

ephys_path = '/home/tom.chartrand/projects/data/human_mouse_ephys_all_0127.csv'
ephys_df = pd.read_csv(ephys_path, index_col=0)

md_path = "/home/tom.chartrand/projects/metadata/consolidated_clusters_and_metadata/human_IVSCC_excitatory_L23_consolidated_0130.csv"
md_df = pd.read_csv(md_path, index_col=0)
human_df = (md_df.join(ephys_df)
#             .drop(columns=['cluster'])
            .rename(columns=dict(zip(variables_long, variables)))
            .loc[lambda df: df.cluster.str.contains("FREM")]
           )

genes = pd.read_csv('/home/tom.chartrand/projects/data/gene_names.csv', index_col=0).iloc[:,0].values
efeatures = list(ephys_df.columns)
shiny_dir = shiny.shiny_directory('human')
join_on = 'sample_id'

def calc(genes):
    genes_df = (pd.read_feather(shiny_dir + '/data.feather', columns=np.append(genes, [join_on]))
                .set_index(join_on)
                .apply(lambda x: np.log2(x+1))
               )
    data = human_df.join(genes_df, on='sample_id')

    results = []
    for gene in genes:
        if data[gene].pipe(lambda x: sum(x>1)) < 5:
            continue
        for efeature in efeatures:
            df = data.dropna(subset=[gene, efeature])
            try:
                pcorr = partial_corr(data=df, x=gene, y=efeature, covar=['depth'])
                pcorr = {
                    "r_pcorr": pcorr['r'][0],
                    "rsquared_pcorr": pcorr['r2'][0],
                    "p_pcorr": pcorr['p-val'][0],
                }
            except np.linalg.LinAlgError:
                pcorr = {
                    "r_pcorr": np.nan,
                    "rsquared_pcorr": np.nan,
                    "p_pcorr": np.nan,
                }
            r, p = pearsonr(df[gene], df[efeature])
            out = {
            "gene":gene,
            "feature":efeature,
            'r_corr':r,
            'rsquared_corr':r**2,
            'p_corr':p,
            }
            out.update(pcorr)
            results.append(out)
    return results
# n=1
# out = list(map(calc, np.array_split(genes[:50], 10)))

# n = sys.argv[1]
step = 1000
nproc = 24
start = int(n)*step
end = min(start+step, len(genes))
pool = Pool(processes=nproc)
out = pool.map(calc, np.array_split(genes[start:end], nproc))

metrics = pd.DataFrame.from_records(chain(*out))
metrics.to_csv(out_path)
