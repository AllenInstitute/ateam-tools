with dendrite_type as
    (
    select sts.specimen_id, st.name
    from specimen_tags_specimens sts
    join specimen_tags st on sts.specimen_tag_id = st.id
    where st.name like 'dendrite type%%'
    )
SELECT sp.id, sp.name, sp.created_at,
csl.normalized_depth, csl.soma_depth_um,
struct.name AS structure, 
donors.name AS donor_name, 
dt.name AS dendrite_type,
layer.name AS layer,
projects.code AS project_code,
nwb.storage_directory || nwb.filename AS nwb_path
FROM specimens sp
JOIN ephys_roi_results err ON sp.ephys_roi_result_id = err.id
LEFT JOIN projects ON projects.id = sp.project_id
LEFT JOIN donors ON donors.id = sp.donor_id
LEFT JOIN structures struct ON sp.structure_id = struct.id
LEFT JOIN structures layer ON sp.cortex_layer_id = layer.id
LEFT JOIN cell_soma_locations csl ON csl.specimen_id = sp.id
LEFT JOIN dendrite_type dt ON sp.id = dt.specimen_id
LEFT JOIN neuronal_models nm ON sp.id = nm.specimen_id
LEFT JOIN
    (
        SELECT nr.specimen_id, f.filename FROM 
        neuron_reconstructions nr JOIN well_known_files f ON nr.id = f.attachable_id 
        WHERE f.well_known_file_type_id = 303941301
        AND NOT nr.superseded
    ) AS nr
ON sp.id = nr.specimen_id
JOIN well_known_files nwb ON nwb.attachable_id = err.id
JOIN well_known_file_types ftype ON nwb.well_known_file_type_id = ftype.id
WHERE nwb.attachable_type = 'EphysRoiResult'
AND ftype.name = 'NWB'