"""Tools for accessing data in web products using AllenSDK
"""

from allensdk.api.queries.biophysical_api import BiophysicalApi
from allensdk.api.queries.cell_types_api import CellTypesApi
from bmtk.simulator.utils.config import ConfigDict
import os.path

_model_types_dict = {"active": 491455321, "peri": 329230710}
bp = BiophysicalApi()

def download_morph_files(specimen_ids, folder_path, add_cre_info=False):
    """Download morphology for a list of cells in Cell Types database
    
    Arguments:
        specimen_ids {iterable} -- IDs as list of int or string
        folder_path {str} -- folder to save to
    """
    ct = CellTypesApi()
    for specimen_id in specimen_ids:
        if add_cre_info:
            info = ct.get_cell(specimen_id)
            filename = info['line_name'].split('-')[0]+'_'+'{}.swc'.format(specimen_id)
        else:
            filename = '{}.swc'.format(specimen_id)
        save_path = os.path.join(folder_path, filename)
        if not os.path.exists(save_path):
            ct.save_reconstruction(specimen_id, save_path)

def download_models_to_config(specimen_ids, config_path, model_type='peri'):
    """Download model files for a list of cells in Cell Types database.
    Saves to the paths specified for a simulation by its config.json
    
    Arguments:
        specimen_ids {iterable} -- IDs as list of int or string
        config_path {str} -- path to BMTK config file
    
    Keyword Arguments:
        model_type {str} -- 'active' or 'peri' (default: {'peri'})
    """
    config_folder = os.path.dirname(config_path)
    conf = ConfigDict.load(config_path)
    fit_path = os.path.normpath(os.path.join(config_folder, conf.biophysical_neuron_models_dir))
    morph_path = os.path.normpath(os.path.join(config_folder, conf.morphologies_dir))
    download_cell_model_files(specimen_ids, fit_path, morph_path, model_type)

def download_cell_model_files(specimen_ids, fit_path, morph_path, model_type='peri'):
    """Download model files for a list of cells in Cell Types database
    
    Arguments:
        specimen_ids {iterable} -- IDs as list of int or string
        fit_path {str} -- path to save fit params
        morph_path {str} -- path to save morphology SWCs
    
    Keyword Arguments:
        model_type {str} -- 'active' or 'peri' (default: {'peri'})
    """
    models = bp.get_neuronal_models(specimen_ids, model_type_ids=[_model_types_dict[model_type]])
    for model in models:
        file_dict = bp.get_well_known_file_ids(model['id'])
        specimen = model['specimen_id']
        _download_file(file_dict['fit'], fit_path, rename=specimen)
        _download_file(file_dict['morphology'], morph_path, rename=specimen)

def _download_file(id_dict, working_directory, rename=None):
    for well_known_id, filename in id_dict.items():
        well_known_file_url = bp.construct_well_known_file_download_url(well_known_id)
        if rename:
            ext = os.path.splitext(filename)[1]
            filename = str(rename) + ext
        save_path = os.path.join(working_directory, filename)
        bp.retrieve_file_over_http(well_known_file_url, save_path)
