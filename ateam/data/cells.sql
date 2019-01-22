SELECT sp.id, sp.name, sp.cell_depth,
st.name as structure, donors.name AS donor_name, 
layer.name AS layer,
projects.code as project_code,
nwb.storage_directory || nwb.filename AS nwb_path, nwb.filename
FROM specimens sp
JOIN ephys_roi_results err ON sp.ephys_roi_result_id = err.id
LEFT JOIN projects ON projects.id = sp.project_id
LEFT JOIN donors ON donors.id = sp.donor_id
LEFT JOIN structures st ON sp.structure_id = st.id
LEFT JOIN structures layer ON sp.cortex_layer_id = layer.id
JOIN well_known_files nwb ON nwb.attachable_id = err.id
JOIN well_known_file_types ftype ON nwb.well_known_file_type_id = ftype.id
WHERE nwb.attachable_type = 'EphysRoiResult'
AND ftype.name = 'NWB'