

import os
from pprint import pprint as pp

from ice.analysis import single_sanger_analysis, multiple_sanger_analysis



definition_file = os.path.abspath('/Users/nicholas.rossi/Documents/Timeline/2019/05/ICE_V2_pull_request/ying_data/mznqj4xcyudg3d29/truncated.xlsx')
data_directory = os.path.abspath('/Users/nicholas.rossi/Documents/Timeline/2019/05/ICE_V2_pull_request/ying_data/mznqj4xcyudg3d29/ab1s')
output_dir = '/Users/nicholas.rossi/Documents/Timeline/2019/05/ICE_V2_pull_request/ying_data/mznqj4xcyudg3d29/output'


job_args = (definition_file, output_dir)
job_kwargs = {
    'verbose': True,
    'data_dir': data_directory
}

multiple_sanger_analysis(*job_args, **job_kwargs)