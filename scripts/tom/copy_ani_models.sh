#!/bin/sh
cd "/allen/programs/celltypes/workgroups/humancolumn_ephysmodeling/anin/Optimizations_HPC/Human/Prod_1/"

EXP="./[0-9]+/Stage2/fitted_params/optim_param_[0-9]+.json"
find . "$BASE" -type f -regex "$EXP" -exec cp {} /home/tom.chartrand/network/bmtk_networks/biophys_components_shared/biophys_params/ani_new_met/ \;


EXP="./[0-9]+/cell_types/[0-9]+.swc"
find . "$BASE" -type f -regex "$EXP" -exec cp {} /home/tom.chartrand/network/bmtk_networks/biophys_components_shared/morphologies/ \;