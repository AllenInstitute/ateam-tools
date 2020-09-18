from ipfx.bin.pipeline_from_specimen_id import run_pipeline_from_id
from ateam.data.convert_ipfx_output_full import extract_fx_output
import allensdk.core.json_utilities as ju
import ateam.data.shiny as shiny
import os
from functools import partial

run_name = 'pipeline_inh_l1'
output_dir = f"/home/tom.chartrand/projects/ephys_analysis/data/{run_name}"

import psutil
from multiprocessing import Pool
pool = Pool(psutil.cpu_count(logical=False))

shiny_df = shiny.load_shiny_data('human', drop_offpipeline=True, nms_pass=False)
inh_df = shiny_df[shiny_df.topLeaf.str.contains('Inh L1') & shiny_df.seurat_cluster.str.contains('Inh L1')].copy()
cells = inh_df.index.astype(str).values
cells = [cell for cell in cells if not os.path.isfile(os.path.join(output_dir, cell, 'pipeline_output.json'))]
pool.map(partial(run_pipeline_from_id, output_dir=output_dir), cells[:])

cells = os.listdir(output_dir)
records = []
for cell in cells:
    path = os.path.join(output_dir, cell, 'pipeline_output.json')
    record = extract_fx_output(ju.read(path), v2=False, extra=[])
    record["specimen_id"] = cell
    records.append(record)
ephys_df = pd.DataFrame.from_records(records, index="specimen_id")
ephys_df.to_csv(f"/home/tom.chartrand/projects/ephys_analysis/data/{run_name}.csv")

print(F"{len(records)} cells processed successfully")

# run_hpc 'python /home/tom.chartrand/work/ateam-tools/test/test_ipfx_rerun.py' --conda ipfx3_home -t 04:00:00 --ppn 32 --mem 64g
