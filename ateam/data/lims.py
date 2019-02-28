"""Tools for accessing internal data directly through the LIMS database
"""
import warnings
import os.path
import pandas as pd
import sqlalchemy as sql

URI_LIMS = 'postgresql://limsreader:limsro@limsdb2/lims2'
TRIPLE = 'Short Square - Triple'

class LimsReader(object):
    """Class for reading cell and sweep info from LIMS ephys records"""

    def __init__(self, base_uri=URI_LIMS):
        super(LimsReader, self).__init__()
        self.engine = sql.create_engine(base_uri)

    def get_cells_df(self, project_id=None, project_code=None, cells_list=None):
        """Get info for cell specimens in LIMS with ephys data
        Select by project (code or id) or list of cell specimen IDs
        
        Keyword Arguments:
            project_id -- numerical project ID
            project_code -- project code (e.g. hIVSCC-MET)
            cells_list {list} -- list of cell specimen IDs (as int or str)
        
        Returns:
            DataFrame -- cell specimen info, indexed by specimen ID
        """
        with open(os.path.join(os.path.dirname(__file__), 'cells.sql'), 'r') as sqlfile:
            sql = sqlfile.read()
        if project_id:
            sql += " AND sp.project_id = {}".format(project_id)
        if project_code:
            sql += " AND projects.code = '{}'".format(project_code)
        if cells_list is not None:
            sql += " AND sp.id IN ({})".format(", ".join([str(cell) for cell in cells_list]))
        cells_df = pd.read_sql(sql, self.engine, index_col='id')
        return cells_df

    def get_sweeps(self, cell_id, sweep_type):
        """Get a list of sweeps for a single cell specimen, by sweep type name
        """
        sql = """SELECT sw.sweep_number
            FROM ephys_sweeps sw
            JOIN ephys_stimuli stim ON stim.id = sw.ephys_stimulus_id
            JOIN ephys_stimulus_types stype ON stype.id = stim.ephys_stimulus_type_id
            WHERE sw.specimen_id = %s AND stype.name LIKE '%s'
            """ % (cell_id, sweep_type)
        test = self.engine.execute(sql)
        sweeps = [s[0] for s in test.fetchall()]
        return sweeps

    def get_sweep_info(self, cell_id, sweep_type=None):
        """Get a table of information for all sweeps of a given cell specimen
        """
        sql = """SELECT sw.sweep_number, sw.workflow_state,
            stype.name, stim.description
            FROM ephys_sweeps sw
            JOIN ephys_stimuli stim ON stim.id = sw.ephys_stimulus_id
            JOIN ephys_stimulus_types stype ON stype.id = stim.ephys_stimulus_type_id
            WHERE sw.specimen_id = %s
            """ % (cell_id)
        if sweep_type:
            sql += " AND stype.name LIKE '{}'".format(sweep_type)
        return pd.read_sql(sql, self.engine, index_col='sweep_number')
  
    def get_nwb_path_from_lims(self, cell_id):
        sql = """
            SELECT nwb.storage_directory || nwb.filename AS nwb_path
            FROM specimens sp
            JOIN ephys_roi_results err ON sp.ephys_roi_result_id = err.id
            JOIN well_known_files nwb ON nwb.attachable_id = err.id
            JOIN well_known_file_types ftype ON nwb.well_known_file_type_id = ftype.id
            WHERE sp.id = {id}
            AND nwb.attachable_type = 'EphysRoiResult'
            AND ftype.name = 'NWB'
            """
        return self.single_result_query(sql.format(id=cell_id))
    
    def get_swc_path_from_lims(self, cell_id, manual_only=True):
        sql = """
            SELECT f.storage_directory || f.filename AS nwb_path FROM 
            neuron_reconstructions n JOIN well_known_files f ON n.id = f.attachable_id 
            AND n.specimen_id = {id} AND f.well_known_file_type_id = 303941301
            AND NOT n.superseded
            """
        # Not sure whether this manual flag is needed, found it in Nathan's code...
        if manual_only:
            sql += "AND n.manual"
        return self.single_result_query(sql.format(id=cell_id))

    def single_result_query(self, sql):
        test = self.engine.execute(sql)
        results = [s[0] for s in test.fetchall()]
        if len(results)>1:
            warnings.warn("Multiple results found in LIMS, expected single.")
        if len(results)==0:
            warnings.warn("No results found in LIMS.")
            return None
        return results[0]



from allensdk.core.nwb_data_set import NwbDataSet
import matplotlib.pyplot as plt
import numpy as np

def plot_sweep_lims(cell_id, sweep_num):
    v, i, t = get_sweep_v_i_t_lims(cell_id, sweep_num)
    plot_sweep(v, i, t)

def plot_sweep(v, i, t):
    fig, axes = plt.subplots(2, 1, sharex=True)
    axes[0].plot(t, v, color='black')
    axes[1].plot(t, i, color='gray')

    axes[0].set_ylabel("mV")
    axes[0].set_title("Voltage")
    axes[1].set_title("Stimulus")
    axes[1].set_ylabel("pA")
    axes[1].set_xlabel("ms")

def get_sweep_v_i_t_from_set(data_set, sweep_num):
    sweep_data = data_set.get_sweep(sweep_num)
    i = sweep_data["stimulus"] # in A
    v = sweep_data["response"] # in V
    i *= 1e12 # to pA
    v *= 1e3 # to mV
    sampling_rate = sweep_data["sampling_rate"] # in Hz
    t = np.arange(0, len(v)) * (1.0 / sampling_rate)
    return v, i, t

def get_sweep_v_i_t_lims(cell_id, sweep_num):
    lr = LimsReader()
    nwb_path = lr.get_nwb_path_from_lims(cell_id)
    dataset = NwbDataSet(nwb_path)
    v, i, t = get_sweep_v_i_t_from_set(dataset, sweep_num)
    return v, i, t