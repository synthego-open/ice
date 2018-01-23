import os
import shutil
import tempfile

import pytest

from ice.classes.pair_alignment import PairAlignment


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
