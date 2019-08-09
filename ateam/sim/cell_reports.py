
# from bmtk.simulator.bionet.modules import MembraneReport
# from bmtk.simulator.bionet.modules.sim_module import SimulatorMod
from ateam.sim.morph import morph_props
import pandas as pd

def save_morph_single(sim, filename):
    gids = sim.biophysical_gids
    props = morph_props(sim, gids[0])
    df = pd.DataFrame.from_dict(props)
    df.to_csv(filename)

# class SaveMorphology(SimulatorMod):
#     def initialize(self, sim):

# class EcpCurrent(MembraneReport):
#     """Special case for when only needing to save the soma variable"""
#     def __init__(self, tmp_dir, file_name, variable_name, cells, sections='all', buffer_data=True):
#         super(EcpCurrent, self).__init__(tmp_dir=tmp_dir, file_name=file_name, variable_name=variable_name, cells=cells,
#                                          sections=sections, buffer_data=buffer_data)