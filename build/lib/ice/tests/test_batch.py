import os

from pprint import pprint as pp
import pytest

from ice.analysis import multiple_sanger_analysis
from ice.tests.fixtures import data_dir, temp_dir
from ice.classes.edit_proposal_creator import EditProposalCreator


def test_batch_analysis(temp_dir, data_dir):
    """ Execute a batch analysis"""

    definition_file = os.path.join(data_dir, 'batch_example.xlsx')
    data_directory = data_dir
    output_dir = os.path.join(temp_dir, 'batch_analysis')

    job_args = (definition_file, output_dir)
    job_kwargs = {
        'verbose': True,
        'data_dir': data_directory
    }
    results = multiple_sanger_analysis(*job_args, **job_kwargs)
    pp(results)

    for idx, result in enumerate(results):
        if idx == 1:
            assert result['sample_name'] == 'multiplex_example1'
            guides = []
            expected_guides = ['UGCCAGGAUCACCUCCGAGA'.replace('U', 'T'),
                               'CGAUAGGGGGGCCUUCUCGG'.replace('U', 'T'),
                               'GCGUCCUCUUAUCUUCUGCC'.replace('U', 'T')]
            for g in result['guides']:
                guides.append(g['sequence'])
            assert set(expected_guides) == set(guides)
        if idx == 6:
            assert result['sample_name'] == 'test_revcomp'
            assert result['guides'][0]['sequence'] == 'CCAGAGGCTGATGCTCACCA'
        if idx == 8:
            msg = "Could not analyze homologous recombination case: "
            msg += "Homology arms of length {} not found in control sequence".format(
                EditProposalCreator.MIN_HOMOLOGY_ARM)
            assert result['notes'] == msg


def test_bad_batch(temp_dir, data_dir):
    definition_file = os.path.join(data_dir, 'bad_batch_definition.xlsx')
    data_directory = data_dir
    output_dir = os.path.join(temp_dir, 'bad_batch_analysis')

    job_args = (definition_file, output_dir)
    job_kwargs = {
        'verbose': True,
        'data_dir': data_directory
    }
    results = multiple_sanger_analysis(*job_args, **job_kwargs)
    assert not results
