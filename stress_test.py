
from ice.analysis import single_sanger_analysis, multiple_sanger_analysis
import os


### SINGLE PLEX

def run_single_plex():

    input_sheet='/Volumes/GoogleDrive/Shared drives/Bioinformatics/data/ICEvNGS/ICE/no_donor/ice_input.xlsx'
    ab1_dir='/Volumes/GoogleDrive/Shared drives/Bioinformatics/data/ICEvNGS/ICE/no_donor/ab1s'
    output='/Users/nicholas.rossi/Documents/Repos/2021/05/allspice_refactor/singleplex'
    job_args = (input_sheet, output)

    job_kwargs = {
        'verbose': False,
        'data_dir': ab1_dir
    }

    results=multiple_sanger_analysis(*job_args, **job_kwargs)


def run_mg_titration(upper_out):
    definition_file = '/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/multiguide_titration/file_mapping.xlsx'
    data_directory='/Users/nicholas.rossi/Documents/Repos/2021/04/allspice_thunderdome/multiguide_titration/ab1s'
    output_dir = os.path.join(upper_out,'mg_titration')


    job_args = (definition_file, output_dir)

    job_kwargs = {
        'verbose': False,
        'data_dir': data_directory
    }

    results=multiple_sanger_analysis(*job_args, **job_kwargs)

def run_multi_plex(upper_out):
    definition_file='/Volumes/GoogleDrive/Shared\ drives/Bioinformatics/data/multiplex/all_multiplex_truncated.xlsx'
    ab1s='/Users/nicholas.rossi/Documents/Repos/2020/3/multiguide_data_roundup/abs1/'
    output_dir = os.path.join(upper_out,'multiplex')
    job_args = (definition_file, output_dir)

    job_kwargs = {
        'verbose': False,
        'data_dir': ab1s
    }

    results=multiple_sanger_analysis(*job_args, **job_kwargs)


def run_cas12a(upper_out):
    definition_file='/Volumes/GoogleDrive/Shared\ drives/Bioinformatics/data/multiplex/all_multiplex_truncated.xlsx'
    ab1s='/Users/nicholas.rossi/Documents/Repos/2020/3/multiguide_data_roundup/abs1/'
    output_dir = os.path.join(upper_out,'multiplex')
    job_args = (definition_file, output_dir)

    job_kwargs = {
        'verbose': False,
        'data_dir': ab1s
    }

    results=multiple_sanger_analysis(*job_args, **job_kwargs)


def run_single_clone_deletion(upper_out):
    definition_file='/Users/nicholas.rossi/Documents/Repos/2021/05/deletion_mapping/clone_deletions/file_mapping.csv'
    ab1s='/Users/nicholas.rossi/Documents/Repos/2021/05/deletion_mapping/clone_deletions/'
    output_dir = os.path.join(upper_out,'sg_clone')
    job_args = (definition_file, output_dir)

    job_kwargs = {
        'verbose': False,
        'data_dir': ab1s
    }

    results=multiple_sanger_analysis(*job_args, **job_kwargs)

### titration experiments
if __name__ == '__main__':
    # upper_out='/Users/nicholas.rossi/Documents/Repos/2021/05/inference_window_experiments/if_adjusted_no_wt/'
    # run_mg_titration(upper_out)

    upper_out='/Users/nicholas.rossi/Documents/Repos/2021/05/systematic_allspice_reckoning/sg_clone/'
    run_single_clone_deletion(upper_out)
