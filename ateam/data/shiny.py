import pandas as pd

def load_shiny_data(species, drop_offpipeline=True, nms_pass=True):
    shiny_df = _load_shiny_data(species)
    shiny_df = filter_shiny_data(shiny_df, drop_offpipeline=drop_offpipeline, nms_pass=nms_pass)
    return shiny_df

def _load_shiny_data(species):
    shiny_dir = shiny_directory(species)
    shiny_df = pd.read_feather(shiny_dir + '/anno.feather')
    shiny_df.drop(columns=[col for col in shiny_df.columns if not col.endswith('_label')], inplace=True)
    shiny_df.rename(axis=1, mapper=lambda col: col.replace('_label',''), inplace=True)

    # add a few helpful columns for the human data
    if species=='human':
        shiny_df["leaf_matched_seurat"] = shiny_df.seurat_cluster == shiny_df.cluster
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
    shiny_df = shiny_df[shiny_df['spec_id'] != 'ZZ_Missing']
    assert shiny_df['spec_id'].is_unique
    shiny_df.index = shiny_df['spec_id'].astype(int)
    if drop_offpipeline:
        # ends with MET, not METx, METc etc
        shiny_df = shiny_df[shiny_df[project_col].str.endswith("MET")]
    if nms_pass:
        shiny_df = shiny_df[shiny_df['Norm_Marker_Sum.0.4']=='TRUE']
    return shiny_df

def load_genes_shiny(species, genes, drop_offpipeline=True, nms_pass=True):
    shiny_dir = shiny_directory(species)
    genes_df = pd.read_feather(shiny_dir + '/data.feather', columns=genes)
    shiny_df = _load_shiny_data(species).join(genes_df)
    shiny_df = filter_shiny_data(shiny_df, drop_offpipeline=drop_offpipeline, nms_pass=nms_pass)
    return shiny_df