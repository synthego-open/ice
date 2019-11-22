from ice.utility import sequence
import pytest

def test_rev_comp():
    seq = "ATTCCG"

    assert sequence.reverse_complement(seq) == "CGGAAT"

def test_complement():
    seq = "ATTCCG"

    assert sequence.complement(seq) == "TAAGGC"

def test_rev_transcribe():
    seq = "GCCUuAAA"
    assert sequence.reverse_transcribe(seq) == "CGGAaTTT"


def test_transcribe():
    seq = "GCCAAATTtt"
    assert sequence.transcribe(seq) == "CGGUUUAAaa"


def test_is_nuc_acid():
    good_seq = "GCCAat"
    assert sequence.is_nuc_acid(good_seq)

    bad_seq = "AAAAXGG"
    assert not sequence.is_nuc_acid(bad_seq)

    not_str = 1234
    assert not sequence.is_nuc_acid(not_str)


def test_rna2dna():
    rna = "UGUUAACG"

    assert sequence.RNA2DNA(rna) == "TGTTAACG"

    with pytest.raises(Exception):
        assert sequence.RNA2DNA(123)

def test_dna2rna():
    dna = "TGTTAAGC"
    assert sequence.DNA2RNA(dna) == "UGUUAAGC"