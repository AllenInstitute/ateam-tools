# Granule cell ephys feature extraction and data analysis code

This code contains two notebooks describing the feature extraction. 

The first notebook 'extract_features_and_save.ipynb' extracts features from standard nwb files and saves them as a feature matrix in csv format 'ALL_EPHYS_FEATURES_EXP_hipp_GC.csv'. File names are taken as cell ids and saved in 'ALL_EPHYS_EXP_CELL_CELL_NAMES_hipp_GC.csv'. The order of features is fixed and is the same for all the cells and is saved into 'ALL_EPHYS_EXP_FEATURE_NAMES_hipp_GC.csv'. The original data was taken from 8 patients suffering from hippocampal epilepsy. All granule cells and corresponding data points are labeled based on the patient, the order is saved in 'ALL_EPHYS_EXP_PATIENT_LABELS_hipp_GC.csv'.

The second notebook 'analyze_features.ipynb' contains the feature analysis code. The first part corresponds to the feature download from .csv files. The feature extraction procedure takes a lot of time, therefore it is better to save the results in csv files and start the analysis from there. All results in various analysis types will be saved in eps and svg formats in the corresponding folder.

The original nwb files take too much space and could not be uploaded on github directly. If you would like to run feature extraction yourself, you could find the files in nwb subfolder of each directory:

H16.06.008
https://drive.google.com/open?id=15XlUQFNo_km1eel08zPQ3Y4EBB1S30tV

H16.06.013
https://drive.google.com/open?id=16_jsZMSUQG00zH3Vh3cjIBXjOX0W2Btl

H17.06.012
https://drive.google.com/open?id=1qhMmbuqw8ElAamPMmIg6qOLj-AQtCjsC

H17.06.014
https://drive.google.com/open?id=1Ot2huVWtI2SP1Rg6jCM7w5SaftvGespf

H17.06.015
https://drive.google.com/open?id=1XgonSxgF4fMDnOpcnEQ-ywK9v7UiaMjx

H18.06.366
https://drive.google.com/open?id=1QE0Oldc051JC478LZblARIhJ4ypflw28

H18.06.368
https://drive.google.com/open?id=1S7CnxNzZdWrg6shldMLRbazmi1xbjtvi

H18.06.371
https://drive.google.com/open?id=19gD6CrjI444DRPNlZfFU3crDa9hgFz1f


If files are not converted from IGOR pro to standard NWB, one could use the converter 'convert_igor_nwb.py'. Put the script into the directory with nwb files and run it in python. It would scan the directory for nwb files and try convert them based on the aibs_stimulus_description value in every sweep:

stimulus_path=str('acquisition/timeseries/Sweep_')+str(num)+str('/aibs_stimulus_description/')

Then the results will be saved in -converted.nwb files. The stimulus classification is done based on the following description:

http://confluence.corp.alleninstitute.org/display/IT/IVSCC+Ephys+NWB+stimulus+summary


## Getting Started

To start data analysis using 'analyze_features.ipynb' it is enougth to have the libraries installed.


### Prerequisites

Before you start it is important to have the following Python libraries:

matplotlib
numpy
scipy
pandas
csv
math
h5py
allensdk
sklearn


## Authors

* **Anatoly Buchin** - *The ephys data analysis workflow for hippocampal granule cells* - [abuchin](https://github.com/abuchin)


## License

Copyright 2017-2019. Allen Institute. All rights reserved

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
following conditions are met:

1 Redistributions of source code must retain the above copyright notice, this list of conditions and the following
disclaimer.

2 Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided with the distribution.

3 Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Acknowledgments

This code has been developed with help of Costas Anastassiou, Nathan Gouwens, David Feng and other contributors of AllenSDK

https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/costas-anastassiou/
https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/nathan-gouwens/
https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/david-feng/

https://allensdk.readthedocs.io/en/latest/
