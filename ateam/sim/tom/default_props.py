
def fit_file_active(cell_id):
    return 'optim_param_%s.json' % cell_id

def fit_file_peri(cell_id):
    return '%s.json' % cell_id

def cell_name(cell_id):
    return 'Human_%s' % cell_id

def morph_file(cell_id):
    return '%s.swc' % cell_id

def cellprops_active(cell_id):
    biophys_props = {
        'cell_name': cell_name(cell_id),
        'morphology': morph_file(cell_id),
        'dynamics_params': fit_file_active(cell_id),
        'model_type': 'biophysical',
        'model_template': 'ctdb:Biophys1.hoc',
        'model_processing': 'aibs_allactive'
    }
    return biophys_props

def cellprops_peri(cell_id):
    biophys_props = {
        'cell_name': cell_name(cell_id),
        'morphology': morph_file(cell_id),
        'dynamics_params': fit_file_peri(cell_id),
        'model_type': 'biophysical',
        'model_template': 'ctdb:Biophys1.hoc',
        'model_processing': 'aibs_perisomatic'
    }
    return biophys_props

cellprops_virtual = {
    'model_type': 'virtual'
}