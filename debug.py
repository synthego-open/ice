import os
from pprint import pprint as pp

from ice.analysis import single_sanger_analysis, multiple_sanger_analysis
#
# # debug sample
# upper_dir='/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/multiguide_titration'
# # # upper_dir='/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/alt_pool_sanger/ab1s'
# #
# #
# # definition_file = '/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/multiguide_titration/file_mapping.xlsx'
# definition_file = '/Users/nicholas.rossi/Documents/Repos/2021/05/allspice_refactor/titration_debug.xlsx'
# #definition_file ='/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/alt_pool_sanger_submission.xlsx'
# # data_directory = '/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/alt_pool_sanger/ab1s'
# data_directory='/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/multiguide_titration/ab1s'
# output_dir = '/Users/nicholas.rossi/Documents/Repos/2021/05/allspice_refactor/debug'
#
#
# job_args = (definition_file, output_dir)
#
# job_kwargs = {
#     'verbose': True,
#     'data_dir': data_directory
# }
#
# results=multiple_sanger_analysis(*job_args, **job_kwargs)



# working on cas12

def titraion_debug():


    definition_file = '/Users/nicholas.rossi/Documents/Repos/2021/05/allspice_refactor/titration_debug.xlsx'
    data_directory='/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/multiguide_titration/ab1s'
    output_dir = '/Users/nicholas.rossi/Documents/Repos/2021/05/allspice_refactor/titration_debug'


    job_args = (definition_file, output_dir)

    job_kwargs = {
        'verbose': True,
        'data_dir': data_directory
    }

    results=multiple_sanger_analysis(*job_args, **job_kwargs)


def cas12():
    input_sheet='/Users/nicholas.rossi/Documents/Repos/2021/05/allspice_refactor/cas12a_single_debug.csv'
    ab1s='/Users/nicholas.rossi/Documents/Repos/2020/3/allspice_reviews/Janssen_Ratio_AB1_20200319/'
    output='/Users/nicholas.rossi/Documents/Repos/2021/05/allspice_refactor/debug'
    job_args = (input_sheet, output)

    job_kwargs = {
        'verbose': True,
        'data_dir': ab1s
    }

    results=multiple_sanger_analysis(*job_args, **job_kwargs)

def sample_titration_4(upper_out):
    definition_file = '/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/multiguide_titration/4_sample_titration.csv'
    data_directory='/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/multiguide_titration/ab1s'
    output_dir = os.path.join(upper_out,'4_sample_titration_new')


    job_args = (definition_file, output_dir)

    job_kwargs = {
        'verbose': False,
        'data_dir': data_directory
    }

    results=multiple_sanger_analysis(*job_args, **job_kwargs)


if __name__ == '__main__':
    sample_titration_4('/Users/nicholas.rossi/Documents/Repos/2021/05/systematic_allspice_reckoning/')
#runnning all single plex

# input_sheet='/Volumes/GoogleDrive/Shared drives/Bioinformatics/data/ICEvNGS/ICE/no_donor/ice_input.xlsx'
# ab1_dir='/Volumes/GoogleDrive/Shared drives/Bioinformatics/data/ICEvNGS/ICE/no_donor/ab1s'
# output='/Users/nicholas.rossi/Documents/Repos/2021/05/allspice_refactor/singleplex'
# job_args = (input_sheet, output)
#
# job_kwargs = {
#     'verbose': True,
#     'data_dir': ab1_dir
# }
#
# results=multiple_sanger_analysis(*job_args, **job_kwargs)