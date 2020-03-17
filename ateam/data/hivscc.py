import pandas as pd
import numpy as np
from ateam.data import shiny

depth = "L23_cell_depth"
cluster = "SeuratMapping"
order=[ 'LTK', 'GLP2R', 'FREM3', 'CARM1P1', 'COL22A1',
      'Adamts2', 'Rrad', 'Agmat', ]

def load_data():
    ephys_df = pd.read_csv('../data/human_mouse_ephys_all_0127.csv', index_col=0).dropna(how='all')
    morph_df = pd.concat([
        pd.read_csv("../data/current/All_Mouse_Cells_Lockdown_All_raw_features.csv", index_col=0),
        pd.read_csv("../data/current/All_L23_Lockdown_all_raw_features.csv", index_col=0),
    ]).assign(
        total_area=lambda data: data['apical_dendrite_total_surface_none'] + data['basal_dendrite_total_surface_none'])

    md_path = "../metadata/consolidated_clusters_and_metadata/human_IVSCC_excitatory_L23_consolidated_0131.csv"
    md_df = pd.read_csv(md_path, index_col=0)
    human_df = (md_df.join(ephys_df).join(morph_df)
                .assign(species='human')
                .pipe(fix_df)
            )

    print(f"{len(human_df)} human records")

    md_path = "../metadata/consolidated_clusters_and_metadata/mouse_IVSCC_excitatory_L23_consolidated_0129.csv"
    md_df = pd.read_csv(md_path, index_col=0)
    mouse_df = (md_df.join(ephys_df).join(morph_df)
                .assign(species='mouse')
                .pipe(fix_df)
            )
    print(f"{len(mouse_df)} mouse records")
    return human_df, mouse_df, ephys_df, morph_df.drop('total_area', axis=1)

def fix_df(df):
    return (df.assign(cluster=lambda df: df[cluster]
                        .apply(lambda name: name.split(' ')[-1])
                        .astype('category', categories=order, ordered=True),
                    depth=lambda df: df[depth])
    #             .sort_values('cluster')
                .sample(frac=1, random_state=42))

def join_gene_data(df, genes):
    shiny_dir = shiny.shiny_directory('human')
    join_on = 'sample_id'
    genes_df = (pd.read_feather(shiny_dir + '/data.feather', columns=genes+[join_on])
                .set_index(join_on)
                .apply(lambda x: np.log2(x+1))
               )
    return df.join(genes_df, on='sample_id')