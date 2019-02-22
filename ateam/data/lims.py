"""Tools for accessing internal data directly through the LIMS database
"""
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

    def get_cells(self, project_id=None, project_code=None):
        """Get info for cell specimens in LIMS with ephys data, by project        
        """
        with open(os.path.join(os.path.dirname(__file__), 'cells.sql'), 'r') as sqlfile:
            sql = sqlfile.read()
        if project_id:
            sql += " AND sp.project_id = {}".format(project_id)
        if project_code:
            sql += " AND projects.code = '{}'".format(project_code)
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
        test = self.engine.execute(sql.format(id=cell_id))
        results = [s[0] for s in test.fetchall()]
        return results[0]

