
import os
from pprint import pprint as pp

import pytest

from ice.analysis import single_sanger_analysis
from ice.classes.sanger_analysis import SangerAnalysis
from ice.tests.fixtures import data_dir, temp_dir


def test_high_quality_sample(temp_dir):
    """ Running a good sample"""

    control = os.path.join(data_dir(), "good_example_control.ab1")
    sample = os.path.join(data_dir(), "good_example_edited.ab1")

    output_path = os.path.join(temp_dir, 'high_quality')
    guide = "AACCAGTTGCAGGCGCCCCA"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True, 'allprops':True}


    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice'] == 78
    allproposals = os.path.join(temp_dir, 'high_quality.allproposals.json')
    assert os.path.exists(allproposals)


def test_low_quality_control(temp_dir):
    #Low quality control but analysis still possible
    control = os.path.join(data_dir(), "low_quality_control2.ab1")
    sample = os.path.join(data_dir(), "low_quality_edited2.ab1")
    guide = 'CUGGUCGCGCACGAUCAGGG'
    output_path = os.path.join(temp_dir, 'low_quality')

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)
    assert 'Low quality control trace' in results['notes']

def test_low_quality_control_fail(temp_dir):
    #Low quality control trace --> no analysis possible
    bad_quality = os.path.join(data_dir(), "low_quality_control.ab1")
    sample = os.path.join(data_dir(), "good_example_edited.ab1")
    guide = 'CCCGCCTCGGCCTCCCTAA'
    output_path = os.path.join(temp_dir, 'low_quality')

    job_args = (bad_quality, sample , output_path, guide)
    job_kwargs = {'verbose': True}
    results = single_sanger_analysis(*job_args, **job_kwargs)

    assert results['status'] == 'succeeded'
    assert 'Control ab1 trace quality scores too low' in results['notes']


def test_sanger_analysis_bad_path(temp_dir):
    sa = SangerAnalysis()
    bad_file = os.path.join(data_dir(), "no_file_here.xyz")
    control = os.path.join(data_dir(), "good_example_control.ab1")
    sample = os.path.join(data_dir(), "good_example_edited.ab1")

    guide = "CGGCUCAAAAGGACCAGACU"
    with pytest.raises(Exception):
        sa.initialize_with(bad_file, sample, guide)
    with pytest.raises(Exception):
        sa.initialize_with(control, bad_file, guide)
    with pytest.raises(Exception):
        sa.initialize_with(control, sample, 123)
    with pytest.raises(Exception):
        sa.initialize_with(control, sample, "1234")

def test_low_quality_sample(temp_dir):
    """ Running a low quality sample"""

    control = os.path.join(data_dir(), "low_quality_control.ab1")
    sample = os.path.join(data_dir(), "low_quality_edited.ab1")

    output_path = os.path.join(temp_dir, 'low_quality')
    guide = "CGGCUCAAAAGGACCAGACU"

    # guide = 'AAAGUCAUCCAAGCCAAGUC'
    # guide = 'CAAAAGGACCAGACUUGGCU'

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    #low quality sample
    results = single_sanger_analysis(*job_args, **job_kwargs)
    assert results['status'] == 'succeeded'
    assert 'Sample ab1 low_quality_edited.ab1 quality scores too low' in results['notes']



def test_bad_paths(temp_dir):
    control = os.path.join(data_dir(), "does_not_exist.abcxyz")
    sample = os.path.join(data_dir(), "alignment_fail_edited.ab1")
    output_path = os.path.join(temp_dir, 'bad_paths')

    guide = "AACCAGTTGCAGGCGCCCCA"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}
    with pytest.raises(Exception):
        results = single_sanger_analysis(*job_args, **job_kwargs)

    good_control = os.path.join(data_dir(), "low_quality_control.ab1")
    bad_sample = os.path.join(data_dir(), "does_not_exist.abcxyz")
    job_args = (good_control, bad_sample, output_path, guide)
    job_kwargs = {'verbose': True}
    with pytest.raises(Exception):
        results = single_sanger_analysis(*job_args, **job_kwargs)

def test_alignment_fail(temp_dir):

    control = os.path.join(data_dir(), "alignment_fail_control.ab1")
    sample = os.path.join(data_dir(), "alignment_fail_edited.ab1")

    output_path = os.path.join(temp_dir, 'alignment_fail')
    guide = "AACCAGTTGCAGGCGCCCCA"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    try:
        results = single_sanger_analysis(*job_args, **job_kwargs)
        assert False

    except Exception as e:
        pp(results)
        assert results['status'] == 'succeeded'
        assert "No alignment found between control and edited sample" in results['notes']


def test_sequence_not_found(temp_dir):

    control = os.path.join(data_dir(), "bad_guide_sequence_control.ab1")
    sample = os.path.join(data_dir(), "bad_guide_sequence_edited.ab1")

    output_path = os.path.join(temp_dir, 'bad_sequence')
    guide = "AACCAGTTGCAGGCGCCCCA"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    try:
        results = single_sanger_analysis(*job_args, **job_kwargs)
        assert False

    except Exception as e:
        pp(results)
        assert results['status'] == 'succeeded'
        assert 'guide AACCAGTTGCAGGCGCCCCA not found in control sequence' in results['notes']


def test_donor_example(temp_dir):
    """Running a good sample with HDR"""

    control = os.path.join(data_dir(), "donor_example_control.ab1")
    sample = os.path.join(data_dir(), "donor_example_knockin.ab1")

    output_path = os.path.join(temp_dir, 'donor_example')
    guide = "AAGTGCAGCTCGTCCGGCGT"
    donor = 'ATCCTCCCGGGAACGTCTCCACCAGCTTCCCTTCCAGCCGACGAGATTGATCTCCGACCCGACGAGCTGCACTTCCTGTCCAAGCACTTCCGCAGCTCAGAGAA'

    job_args = (control, sample, output_path, guide, donor)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['hdr_pct'] > 20



def test_multiplex_0_0(temp_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir(), "multiplex0_control.ab1")
    sample = os.path.join(data_dir(), "multiplex0_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex0')
    guide = "UGCCAGGAUCACCUCCGAGA"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice_d'] == 66



def test_multiplex_0_1(temp_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir(), "multiplex0_control.ab1")
    sample = os.path.join(data_dir(), "multiplex0_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex0')
    guide = "CGAUAGGGGGGCCUUCUCGG"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice_d'] == 66

def test_multiplex_0_2(temp_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir(), "multiplex0_control.ab1")
    sample = os.path.join(data_dir(), "multiplex0_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex0')
    guide = "GCGUCCUCUUAUCUUCUGCC"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}


    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice_d'] == 66


def test_multiplex_1_0(temp_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir(), "multiplex1_control.ab1")
    sample = os.path.join(data_dir(), "multiplex1_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex1')
    guide = "UUCCCCUUAUUUAGAUAACU"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice_d'] == 67


def test_multiplex_1_1(temp_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir(), "multiplex1_control.ab1")
    sample = os.path.join(data_dir(), "multiplex1_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex1')
    guide = "UCCCCUUAUUUAGAUAACUC"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}


    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice_d'] == 67


def test_multiplex_three_guides(temp_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir(), "multiplex1_control.ab1")
    sample = os.path.join(data_dir(), "multiplex1_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex1')
    guide = "UUCCCCUUAUUUAGAUAACU,UCCCCUUAUUUAGAUAACUC,UGUUAACCAAAUUAUGAUGA"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}


    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice'] == 69
    assert len(results['guides']) == 3
