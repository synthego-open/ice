

import os
from pprint import pprint as pp

from ice.analysis import single_sanger_analysis, multiple_sanger_analysis



upper_dir='/Users/nicholas.rossi/Documents/Repos/2020/3/allspice_reviews/'
definition_file = os.path.abspath(upper_dir+'debug_jansen.xlsx')
data_directory = os.path.abspath(upper_dir+'Janssen_Ratio_AB1_20200319/')
output_dir = upper_dir+'allspice'


job_args = (definition_file, output_dir)
job_kwargs = {
    'verbose': True,
    'data_dir': data_directory
}

multiple_sanger_analysis(*job_args, **job_kwargs)