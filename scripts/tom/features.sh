#!/bin/sh
source activate ipfx3
cd /home/tom.chartrand/projects/ephys_analysis/data/ipfx
python -m ipfx.bin.run_nathan_feature_collection --input_file mouse_tea_cells.txt --output_file mouse_tea_features.csv --log_level INFO
python -m ipfx.bin.run_nathan_feature_collection --input_file mouse_v1_cells.txt --output_file mouse_v1_features.csv --log_level INFO
python -m ipfx.bin.run_nathan_feature_collection --input_file human_l23_cells.txt --output_file human_l23_features.csv --log_level INFO