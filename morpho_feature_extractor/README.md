# Morpho feature extraction and data analysis code

This code contains two notebooks describing the feature extraction for morphological data, swc files. The 'analysis_morpho.ipyb' notebook contains the analysis code that loads the feature matrix from 'WG1TS_WG1_WG4_extended.csv'.

The 'feature_extractor_CK' folder contains the code to use in feature extractor. The main program code is located in "bin" subfolder. The input files need to be copied into "swc" subfolder, while the results of the analysis will be saved in "OUT_swc" subfolder.


## Getting Started the analysis

To start data analysis using 'analysis_morpho.ipyb' it is enougth to have the libraries installed.


## Getting Started the feature extraction

1) Go to the 'feature_extractor_CK' folder.
2) Create the folder with the name 'swc'. Copy your swc files into this 'swc' subfolder. The swc files should contain the cell ID name in their name, this way they will be connected to LIMS.
3) Run the following command, which will upload the swc files into LIMS, extract features and save them.

python ./bin/upload.py ./swc ./OUT_swc >& log.txt

4) After feature extraction, check out log.txt, where one could find possible errors. The results will be saved into OUT_swc/morphology_data/. One could access html folder that would contain the visualisations and feature folder, that would contain the feautures. All morphological (independent from translation and rotation) features will be saved in 'non_upright_features.csv' file in the matrix form, where row correspond to cells and columns correspond to features.


### Prerequisites

Before you start it is important to have the following Python libraries:

matplotlib
numpy
scipy
pandas
csv
math
sklearn


## Authors

* **Anatoly Buchin** - *The morpho data analysis workflow* - [abuchin](https://github.com/abuchin)
* **Keith Godfrey** - *The feature extraction code* - (https://www.linkedin.com/in/keith-godfrey-phd-3b96295/)
* **Changkyu CK Lee** - *The feature extraction code* - (https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/changkyu-ck-lee/)


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

This code has been developed with help of Costas Anastassiou, Changkyu Ck Lee, Nathan Gouwens and Staci Sorensen and other contributors of AllenSDK.

https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/costas-anastassiou/
https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/changkyu-ck-lee/
https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/nathan-gouwens/
https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/staci-sorensen/

https://allensdk.readthedocs.io/en/latest/
