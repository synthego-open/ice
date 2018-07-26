import os
import shutil
import tempfile

import pytest

from ice.classes.pair_alignment import DonorAlignment, PairAlignment


@pytest.fixture
def data_dir():
    """ returns the test data directory """

    path = os.path.dirname(os.path.realpath(__file__))
    test_data = os.path.join(path, "test_data/")
    return test_data


@pytest.fixture
def example_alignment():
    """ returns an example alignment obj """
    seq1 = 'AATGTAATGATAG'
    seq2 = 'AATGTATGATAG'
    pa = PairAlignment(seq1, seq2)
    return pa


@pytest.fixture
def donor_alignment_contiguous_insert():
    """
    Returns example of DonorAlignment with size 5 contiguous insert
    Expected alignments:
    aligned control seq 'CCCCTGAAATGTA-----ATGATAGCC'
    aligned donor seq   '----TGAAATGTATTTTTATGATAGCC'
    """
    control_seq = 'CCCCTGAAATGTAATGATAGCC'
    donor_seq = 'TGAAATGTATTTTTATGATAGCC'
    return DonorAlignment(control_seq, donor_seq)


@pytest.fixture
def donor_alignment_noncontiguous_insert():
    """
    Returns example of DonorAlignment with size 9 noncontiguous insert
    Expected alignments:
    aligned control seq 'CCCCTGAAATGTA-----ATGATAGCC--ATGACT'
    aligned donor seq   '----TGAAATGTATTTTTATGATAGCCAATTGACT'
    """
    control_seq = 'CCCCTGAAATGTAATGATAGCCATGACT'
    donor_seq = 'TGAAATGTATTTTTATGATAGCCAATTGACT'
    return DonorAlignment(control_seq, donor_seq)


@pytest.fixture
def donor_alignment_insert_and_deletion():
    """
    Returns example of DonorAlignment with size 5 insert and size 2 deletion
    Expected alignments:
    aligned control seq 'CCCCTGAAATGTA-----ATGATAGCCAATTGACT'
    aligned donor seq   '----TGAAATGTATTTTTATGATAGCC--TTGACT'
    """
    control_seq = 'CCCCTGAAATGTAATGATAGCCAATTGACT'
    donor_seq = 'TGAAATGTATTTTTATGATAGCCTTGACT'
    return DonorAlignment(control_seq, donor_seq)


@pytest.fixture
def donor_alignment_deletion():
    """
    Returns example of DonorAlignment with size 4 deletion
    Expected alignments:
    aligned control seq 'CCCCTGAAATGTAATGATAGCCAATTGACT'
    aligned donor seq   '----TGAAATGTAATGATAG----TTGACT'
    """
    control_seq = 'CCCCTGAAATGTAATGATAGCCAATTGACT'
    donor_seq = 'TGAAATGTAATGATAGTTGACT'
    return DonorAlignment(control_seq, donor_seq)


@pytest.fixture
def shorter_alignment():
    """ returns an example alignment obj """
    seq1 = 'AATGTAATGATAG'
    seq3 = 'AATGTATGATAGTGGG'
    pa = PairAlignment(seq1, seq3)
    return pa


@pytest.fixture
def no_alignment():
    """ returns an example alignment obj """
    seq1 = 'AATGTAATGATAG'
    seq4 = 'GGACAGACTTAAA'
    pa = PairAlignment(seq1, seq4)
    return pa


@pytest.yield_fixture(scope='session', autouse=True)
def temp_dir():
    """ path to a temporary directory that is deleted after testing"""

    # Will be executed before the first test
    dirpath = tempfile.mkdtemp()

    yield dirpath

    # Will be executed after the last test
    shutil.rmtree(dirpath)
