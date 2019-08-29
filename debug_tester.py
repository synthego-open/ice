

import os
from pprint import pprint as pp

from ice.analysis import single_sanger_analysis, multiple_sanger_analysis



definition_file = os.path.abspath('/Users/nicholas.rossi/Documents/Timeline/2019/08/ying_clonal_regression/ICE temp.xlsx')
data_directory = os.path.abspath('/Users/nicholas.rossi/Documents/Timeline/2019/08/ying_clonal_regression/30-274247844_ab1/')
output_dir = '/Users/nicholas.rossi/Documents/Timeline/2019/08/ying_clonal_regression/lasso_output'


job_args = (definition_file, output_dir)
job_kwargs = {
    'verbose': True,
    'data_dir': data_directory
}

multiple_sanger_analysis(*job_args, **job_kwargs)