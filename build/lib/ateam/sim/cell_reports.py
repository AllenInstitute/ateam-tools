
from bmtk.simulator.bionet.modules import MembraneReport
from bmtk.simulator.bionet.modules.sim_module import SimulatorMod
import pandas as pd

def save_morph_single(sim, filename):
    gids = sim.biophysical_gids
    props = morph_props(sim, gids[0])
    df = pd.DataFrame.from_dict(props)
    df.to_csv(filename)

def morph_props(sim, gid):
    cell = sim.net.get_cell_gid(gid)
    seg_props = cell.morphology.seg_prop
    nsegs = len(seg_props['x'])
    seg_coords = cell._seg_coords
    for name in ['p0', 'p1']:
        seg_props.update({name+'_x': seg_coords[name][0], 
                          name+'_y': seg_coords[name][1], 
                          name+'_z': seg_coords[name][2]})
    seg_props.update(seg_index=range(nsegs), gid=nsegs*[gid])
    return seg_props

# class SaveMorphology(SimulatorMod):
#     def initialize(self, sim):

# class EcpCurrent(MembraneReport):
#     """Special case for when only needing to save the soma variable"""
#     def __init__(self, tmp_dir, file_name, variable_name, cells, sections='all', buffer_data=True):
#         super(EcpCurrent, self).__init__(tmp_dir=tmp_dir, file_name=file_name, variable_name=variable_name, cells=cells,
#                                          sections=sections, buffer_data=buffer_data)