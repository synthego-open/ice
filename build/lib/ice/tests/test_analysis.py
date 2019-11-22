import os
from pprint import pprint as pp

import pytest
from unittest.mock import MagicMock

from ice.analysis import single_sanger_analysis
from ice.classes.sanger_analysis import SangerAnalysis
from ice.tests.fixtures import data_dir, temp_dir


def test_high_quality_sample(temp_dir, data_dir):
    """ Running a good sample"""

    control = os.path.join(data_dir, "good_example_control.ab1")
    sample = os.path.join(data_dir, "good_example_edited.ab1")

    output_path = os.path.join(temp_dir, 'high_quality')
    guide = "AACCAGTTGCAGGCGCCCCA"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True, 'allprops': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice'] == 77
    allproposals = os.path.join(temp_dir, 'high_quality.allproposals.json')
    assert os.path.exists(allproposals)


def test_lowercase_guide(temp_dir, data_dir):
    """ Running a good sample"""

    control = os.path.join(data_dir, "good_example_control.ab1")
    sample = os.path.join(data_dir, "good_example_edited.ab1")

    output_path = os.path.join(temp_dir, 'high_quality_lowercase')
    guide = "aaccagttgcaggcgcccca"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True, 'allprops': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice'] == 77
    allproposals = os.path.join(temp_dir, 'high_quality_lowercase.allproposals.json')
    assert os.path.exists(allproposals)


def test_low_quality_control(temp_dir, data_dir):
    # Low quality control but analysis still possible
    control = os.path.join(data_dir, "low_quality_control2.ab1")
    sample = os.path.join(data_dir, "low_quality_edited2.ab1")
    guide = 'CUGGUCGCGCACGAUCAGGG'
    output_path = os.path.join(temp_dir, 'low_quality')

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)
    assert 'Low quality control trace' in results['notes']


def test_low_quality_control_fail(temp_dir, data_dir):
    # Low quality control trace --> no analysis possible
    bad_quality = os.path.join(data_dir, "low_quality_control.ab1")
    sample = os.path.join(data_dir, "good_example_edited.ab1")
    guide = 'CCCGCCTCGGCCTCCCTAA'
    output_path = os.path.join(temp_dir, 'low_quality')

    job_args = (bad_quality, sample, output_path, guide)
    job_kwargs = {'verbose': True}
    results = single_sanger_analysis(*job_args, **job_kwargs)

    assert results['status'] == 'succeeded'
    assert 'Control ab1 trace quality scores too low' in results['notes']


def test_sanger_analysis_bad_path(temp_dir, data_dir):
    sa = SangerAnalysis()
    bad_file = os.path.join(data_dir, "no_file_here.xyz")
    control = os.path.join(data_dir, "good_example_control.ab1")
    sample = os.path.join(data_dir, "good_example_edited.ab1")

    guide = "CGGCUCAAAAGGACCAGACU"
    with pytest.raises(Exception):
        sa.initialize_with(bad_file, sample, guide)
    with pytest.raises(Exception):
        sa.initialize_with(control, bad_file, guide)
    with pytest.raises(Exception):
        sa.initialize_with(control, sample, 123)
    with pytest.raises(Exception):
        sa.initialize_with(control, sample, "1234")


def test_low_quality_sample(temp_dir, data_dir):
    """ Running a low quality sample"""

    control = os.path.join(data_dir, "low_quality_control.ab1")
    sample = os.path.join(data_dir, "low_quality_edited.ab1")

    output_path = os.path.join(temp_dir, 'low_quality')
    guide = "CGGCUCAAAAGGACCAGACU"

    # guide = 'AAAGUCAUCCAAGCCAAGUC'
    # guide = 'CAAAAGGACCAGACUUGGCU'

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    # low quality sample
    results = single_sanger_analysis(*job_args, **job_kwargs)
    assert results['status'] == 'succeeded'
    assert 'Sample ab1 low_quality_edited.ab1 quality scores too low' in results['notes']


def test_bad_paths(temp_dir, data_dir):
    control = os.path.join(data_dir, "does_not_exist.abcxyz")
    sample = os.path.join(data_dir, "alignment_fail_edited.ab1")
    output_path = os.path.join(temp_dir, 'bad_paths')

    guide = "AACCAGTTGCAGGCGCCCCA"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}
    with pytest.raises(Exception):
        results = single_sanger_analysis(*job_args, **job_kwargs)

    good_control = os.path.join(data_dir, "low_quality_control.ab1")
    bad_sample = os.path.join(data_dir, "does_not_exist.abcxyz")
    job_args = (good_control, bad_sample, output_path, guide)
    job_kwargs = {'verbose': True}
    with pytest.raises(Exception):
        results = single_sanger_analysis(*job_args, **job_kwargs)


def test_alignment_fail(temp_dir, data_dir):

    control = os.path.join(data_dir, "alignment_fail_control.ab1")
    sample = os.path.join(data_dir, "alignment_fail_edited.ab1")

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


def test_sequence_not_found(temp_dir, data_dir):

    control = os.path.join(data_dir, "bad_guide_sequence_control.ab1")
    sample = os.path.join(data_dir, "bad_guide_sequence_edited.ab1")

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


def test_donor_example(temp_dir, data_dir):
    """ Running a good sample with HDR """

    control = os.path.join(data_dir, 'donor_example_control.ab1')
    sample = os.path.join(data_dir, 'donor_example_knockin.ab1')

    output_path = os.path.join(temp_dir, 'donor_example')

    # this should result in a +14 insertion
    guide = 'AAGTGCAGCTCGTCCGGCGT'
    donor = 'ATCCTCCCGGGAACGTCTCCACCAGCTTCCCTTCCAGCCGACGAGATTGATCTCCGACCCGACGAGCTGCACTTCCTGTCCAAGCACTTCCGCAGCTCAGAGAA'

    job_args = (control, sample, output_path, guide, donor)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert pytest.approx(results['hdr_pct']) == 33


def test_donor_substitution_example(temp_dir, data_dir):
    """ Running a good sample with donor containing only a 2bp substitution """
    control = os.path.join(data_dir, 'donor_sub_example_control.ab1')
    sample = os.path.join(data_dir, 'donor_sub_example_edited.ab1')

    output_path = os.path.join(temp_dir, 'donor_example')

    # this should result in a 2bp substitution
    guide = 'GCTGCTTCCTGGGGGCGCCT'
    donor = 'AGCCTCAGGGCCTCACACCAGCCCATGTGGATGACCTGAGGGTCCTGTTTCCCATCCCACttCAGGCGCCCCCAGGAAGCAGCGGCGGGAGCGCACCACCTTCACCCGGAGCCAACTGGAGG'

    job_args = (control, sample, output_path, guide, donor)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert pytest.approx(results['hdr_pct']) == 44


def test_multiplex_0_0(temp_dir, data_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir, "multiplex0_control.ab1")
    sample = os.path.join(data_dir, "multiplex0_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex0')
    guide = "UGCCAGGAUCACCUCCGAGA"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice_d'] == 66


def test_multiplex_0_1(temp_dir, data_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir, "multiplex0_control.ab1")
    sample = os.path.join(data_dir, "multiplex0_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex0')
    guide = "CGAUAGGGGGGCCUUCUCGG"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice_d'] == 66


def test_multiplex_0_2(temp_dir, data_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir, "multiplex0_control.ab1")
    sample = os.path.join(data_dir, "multiplex0_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex0')
    guide = "GCGUCCUCUUAUCUUCUGCC"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice_d'] == 66


def test_multiplex_1_0(temp_dir, data_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir, "multiplex1_control.ab1")
    sample = os.path.join(data_dir, "multiplex1_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex1')
    guide = "UUCCCCUUAUUUAGAUAACU"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice_d'] == 67


def test_multiplex_1_1(temp_dir, data_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir, "multiplex1_control.ab1")
    sample = os.path.join(data_dir, "multiplex1_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex1')
    guide = "UCCCCUUAUUUAGAUAACUC"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice_d'] == 67


def test_multiplex_three_guides(temp_dir, data_dir):
    """ Running a good multiplex sample"""

    control = os.path.join(data_dir, "multiplex1_control.ab1")
    sample = os.path.join(data_dir, "multiplex1_edited.ab1")

    output_path = os.path.join(temp_dir, 'multiplex1')
    guide = "UUCCCCUUAUUUAGAUAACU,UCCCCUUAUUUAGAUAACUC,UGUUAACCAAAUUAUGAUGA"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice'] == 63
    assert len(results['guides']) == 3


def test_multiplex_three_guides_no_editing(temp_dir, data_dir):
    """
    Running a zero editing multiplex sample, with guides that are far apart
    """

    control = os.path.join(data_dir, "multiplex1_control.ab1")
    sample = os.path.join(data_dir, "multiplex1_control.ab1")

    output_path = os.path.join(temp_dir, 'multiplex_no_editing_far_guides')

    guide = "UUCCCCUUAUUUAGAUAACU,UCCCCUUAUUUAGAUAACUC,TGAGTTTTTTTGTAAGTAGC"

    job_args = (control, sample, output_path, guide)
    job_kwargs = {'verbose': True}

    results = single_sanger_analysis(*job_args, **job_kwargs)

    pp(results)

    assert results['status'] == 'succeeded'
    assert results['ice'] == 0
    assert len(results['guides']) == 3


def test_should_skip_proposal():
    fake_donor = 'ATCG'
    size_5_insertion = MagicMock(hdr_indel_size=5)
    mocked_sanger_analysis = MagicMock(_should_skip_proposal=SangerAnalysis._should_skip_proposal, donor_odn=fake_donor,
                                       HDR_OVERLAP_FILTER_CUTOFF=3, donor_alignment=size_5_insertion)

    # does pass insert above cutoff that is same size as donor
    assert mocked_sanger_analysis._should_skip_proposal(mocked_sanger_analysis, 5)

    # does not pass if there is not a donor_odn
    mocked_sanger_analysis.donor_odn = None
    assert not mocked_sanger_analysis._should_skip_proposal(mocked_sanger_analysis, 5)
    mocked_sanger_analysis.donor_odn = fake_donor  # add back fake donor for rest of tests

    # does not pass insert above cutoff that isn't same size as donor
    assert not mocked_sanger_analysis._should_skip_proposal(mocked_sanger_analysis, 6)

    # does pass deletion above cutoff
    size_7_deletion = MagicMock(hdr_indel_size=-7)
    mocked_sanger_analysis.donor_alignment = size_7_deletion
    assert mocked_sanger_analysis._should_skip_proposal(mocked_sanger_analysis, -7)

    # does not pass deletion above cutoff that isn't same size as donor
    assert not mocked_sanger_analysis._should_skip_proposal(mocked_sanger_analysis, -5)

    # does not pass insert same size as deletion
    assert not mocked_sanger_analysis._should_skip_proposal(mocked_sanger_analysis, 7)

    # does not pass insert equal to cutoff
    size_3_insertion = MagicMock(hdr_indel_size=3)
    mocked_sanger_analysis.donor_alignment = size_3_insertion
    assert not mocked_sanger_analysis._should_skip_proposal(mocked_sanger_analysis, 3)

    # does not pass deletion below cutoff
    size_2_deletion = MagicMock(hdr_indel_size=-2)
    mocked_sanger_analysis.donor_alignment = size_2_deletion
    assert not mocked_sanger_analysis._should_skip_proposal(mocked_sanger_analysis, -2)
