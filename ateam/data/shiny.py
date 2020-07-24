import pandas as pd
import os.path

def load_shiny_data(species, csv_path=None, drop_offpipeline=True, nms_pass=True):
    shiny_df = _load_shiny_data(species, csv_path)
    shiny_df = filter_shiny_data(shiny_df, drop_offpipeline=drop_offpipeline, nms_pass=nms_pass)
    return shiny_df

def _load_shiny_data(species=None, directory=None, csv_path=None):
    if csv_path:
        shiny_df = pd.read_csv(csv_path)
    else:
        directory = directory or shiny_directory(species)
        path = os.path.join(directory, 'anno.feather')
        shiny_df = pd.read_feather(path)
    shiny_df.drop(columns=[col for col in shiny_df.columns 
        if col.endswith('_color') or (col.endswith('_id') and not col in ['spec_id','sample_id'])], 
        inplace=True)
    shiny_df.rename(axis=1, mapper=lambda col: col.replace('_label',''), inplace=True)
    shiny_df.replace('ZZ_Missing', float('nan'), inplace=True)
    if 'spec_id' in shiny_df.columns:
        shiny_df = shiny_df[shiny_df['spec_id'] != 'ZZ_Missing'] # may be nan now
        shiny_df = shiny_df.dropna(subset=['spec_id'])
        assert shiny_df['spec_id'].is_unique
        shiny_df.index = shiny_df['spec_id'].astype(int)

    # add a few helpful columns for the human data
    if species=='human':
        shiny_df["leaf_matched_seurat"] = shiny_df.seurat_cluster == shiny_df.topLeaf
        shiny_df["leaf_mapped"] = shiny_df.topLeaf == shiny_df.cluster
        shiny_df["L23_depth_normalized"] = shiny_df.L23_cell_depth / shiny_df.L23_total_thickness

    return shiny_df

def shiny_directory(species):
    feather_path = '/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/{}'
    # Maybe hard-code a date since columns may change with version?
    if species=='human':
        feather_path = feather_path.format('human/human_patchseq_MTG_current')
    elif species=='mouse':
        feather_path = feather_path.format('mouse_patchseq_VISp_current')
    return feather_path

def filter_shiny_data(shiny_df, drop_offpipeline=True, nms_pass=True):
    project_col = 'cell_specimen_project'
    if drop_offpipeline:
        # ends with MET, not METx, METc etc
        shiny_df = shiny_df[~shiny_df[project_col].isna() & shiny_df[project_col].str.endswith("MET")]
    if nms_pass:
        shiny_df = shiny_df[shiny_df['Norm_Marker_Sum.0.4']=='TRUE']
    return shiny_df

def load_genes_shiny(genes, species=None, directory=None, csv_path=None, drop_offpipeline=False, nms_pass=False):
    shiny_df = _load_shiny_data(species, directory, csv_path)

    join_on = 'sample_id'
    genes_df = pd.read_feather(os.path.join(directory, 'data.feather'), columns=genes+[join_on]).set_index(join_on)
    shiny_df = shiny_df.join(genes_df, on=join_on)
    shiny_df = filter_shiny_data(shiny_df, drop_offpipeline=drop_offpipeline, nms_pass=nms_pass)
    return shiny_df