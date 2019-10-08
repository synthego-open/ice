"""
Copyright 2018 Synthego Corporation All Rights Reserved

The Synthego ICE software was developed at Synthego Corporation.

Permission to use, copy, modify and distribute any part of Synthego ICE for
educational, research and non-profit purposes, without fee, and without a
written agreement is hereby granted, provided that the above copyright notice,
this paragraph and the following paragraphs appear in all copies.

Those desiring to incorporate this Synthego ICE software into commercial
products or use for commercial purposes should contact Synthego support at Ph:
(888) 611-6883 ext:1, E-MAIL: support@synthego.com.

In no event shall Synthego Corporation be liable to any party for direct,
indirect, special, incidental, or consequential damages, including lost
profits, arising out of the use of Synthego ICE, even if Synthego Corporation
has been advised of the possibility of such damage.

The Synthego ICE tool provided herein is on an "as is" basis, and the Synthego
Corporation has no obligation to provide maintenance, support, updates,
enhancements, or modifications. The Synthego Corporation makes no
representations and extends no warranties of any kind, either implied or
express, including, but not limited to, the implied warranties of
merchantability or fitness for a particular purpose, or that the use of
Synthego ICE will not infringe any patent, trademark or other rights.
"""


import os
from pprint import pprint as pp

from ice.analysis import single_sanger_analysis, multiple_sanger_analysis

# Running a single analysis

control_path = os.path.abspath('./ice/tests/test_data/good_example_control.ab1')
sample_path = os.path.abspath('./ice/tests/test_data/good_example_edited.ab1')
guide = 'AACCAGTTGCAGGCGCCCCA'
base_outputname = './results/good_example'
donor = None
verbose = True

results = single_sanger_analysis(control_path=control_path,
                                 sample_path=sample_path,
                                 base_outputname=base_outputname,
                                 guide=guide,
                                 donor=donor,
                                 verbose=verbose)

pp(results)


# Running a batch analysis

definition_file = os.path.abspath('./ice/tests/test_data/batch_example.xlsx')
data_directory = os.path.abspath('./ice/tests/test_data/')
output_dir = './results/batch_example'

job_args = (definition_file, output_dir)
job_kwargs = {
    'verbose': True,
    'data_dir': data_directory
}

multiple_sanger_analysis(*job_args, **job_kwargs)
