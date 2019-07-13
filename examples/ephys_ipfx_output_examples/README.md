# Analysis of ipfx output

This code analysis the folder and parses the information from the json files inside the folders, where the numbers correspond to the cell ids.


## Getting Started

This code analyses the output ipfx files and saves them as the csv file. To run the analysis, copy the files from ateam/analysis/ephys_ipfx_output folder and run save_jsons_to_csv.py script. Do not forget to uncompress the json files from zip. The output will be saved in the csv file, where columns correspond to cell specimen IDs and rows correspond to cells.


## Remarks

Currently working on getting more additional features, like bursting index. This section is commented in parse_json_ipfx function in utils_ipfx_linear.py. Feel free to uncomment it if you would like to try to make it work.


### Prerequisites

Before you start it is important to have the following Python libraries:

numpy
pandas
lims
allensdk


## Authors

* **Anatoly Buchin** - *The ephys data analysis workflow for hippocampal granule cells* - [abuchin](https://github.com/abuchin)


## License

Copyright 2017-2019. Allen Institute. All rights reserved

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1 Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2 Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3 Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## Acknowledgments

This code has been developed with help of Tom Chartrand, Sergey Gratiy, Costas Anastassiou, Nathan Gouwens, David Feng and other contributors of AllenSDK

https://alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/sergey-gratiy/
https://alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/tom-chartrand/
https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/costas-anastassiou/
https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/nathan-gouwens/
https://www.alleninstitute.org/what-we-do/brain-science/about/team/staff-profiles/david-feng/

https://allensdk.readthedocs.io/en/latest/
