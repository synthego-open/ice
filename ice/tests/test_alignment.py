import pytest

from ice.classes.pair_alignment import PairAlignment
from ice.tests.fixtures import example_alignment, shorter_alignment, no_alignment


def test_windowed_alignment():
    pa = example_alignment()
    pa.align_with_window([0,3])
    assert pa.aln_seqs[0] == "AAT---------"
    assert pa.aln_seqs[1] == "AATGTATGATAG"


def test_coord_conversion():
    '''
    GGGAATGTCCTGATAG
    ---AATGT--TGATAG
    :return:
    '''
    seq1 = 'GGGAATGTCCTGATAG'
    seq2 = 'AATGTTGATAG'
    pa = PairAlignment(seq1, seq2)
    pa.align_with_window([0,14])
    assert 2 == pa.ctrl2sample_coords(5)
    assert 0 == pa.ctrl2sample_coords(3)
    assert 9 == pa.ctrl2sample_coords(14)
    assert None == pa.ctrl2sample_coords(9)
    assert 4 == pa.ctrl2sample_coords(9, closest=True)

    with pytest.raises(Exception):
        pa.ctrl2sample_coords(None)

def test_no_aln_possible():
    seq1 = 'TTTTTTTTTT'
    seq2 = 'AAAAAAAAAAGGG'
    pa = PairAlignment(seq1, seq2)
    flag, msg = pa.align_with_window([0, 14])
    assert not flag
    assert 'No alignment found' in msg

def test_no_aln_done():
    pa = example_alignment()
    with pytest.raises(Exception):
        pa.ctrl2sample_coords(5)


def test_all_alignment():
    pa = example_alignment()
    pa.align_all()
    assert "AATGT-ATGATAG" in pa.all_aligned_clustal

    with pytest.raises(Exception):
        pa.ctrl2sample_coords(None)

    #only windowed alignment results in alingment_pairs
    with pytest.raises(Exception):
        assert 7 == pa.ctrl2sample_coords(8)


def test_aln_str():
    pa = example_alignment()
    pa.align_all()
    assert (str(pa) == 'AATGTAATGATAG\nAATGT-ATGATAG') or \
             (str(pa) == 'AATGTAATGATAG\nAATGTA-TGATAG')

def test_aln_data():
    '''
    The alignment after windowing should ignore the last base of the ref seq

    control                             AATGTAATGATA
    edited                              AATGTATGATAG

    '''
    pa = example_alignment()
    pa.align_with_window([0, 3])
    for idx, pair in enumerate(pa.alignment_pairs):
        assert pair[0] == pair[1]

def test_aln_shorter_sample():
    '''

    The alignment after windowing should have gaps
    AATGTAATGATAG---
    AATGTATGATAGTGGG

    '''
    pa = shorter_alignment()

    pa.align_with_window([0, 3])

    for idx, pair in enumerate(pa.alignment_pairs):
        if idx < 13:
            assert pair[0] == pair[1]
        else:
            assert pair[0] == None

def test_no_aln():
    '''
    This alignment should be None
    '''

    pa = no_alignment()
    pa.align_all()
    assert pa.alignment_pairs is None
