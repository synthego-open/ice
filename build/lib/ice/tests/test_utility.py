import pytest

from ice.utility.misc import version
from ice.utility.sequence import reverse_transcribe, transcribe, is_nuc_acid, reverse_complement, DNA2RNA, RNA2DNA


def test_sequence_util():

    sequence = 'ACTG'
    rna_sequence = 'ACUG'

    transcribed = transcribe(sequence)
    print(transcribed)

    reverse_comp = reverse_complement(sequence)
    print(reverse_comp)

    reverse_transcribed_seq = reverse_transcribe(rna_sequence)
    print(reverse_transcribed_seq)

    assert is_nuc_acid('A')
    assert not is_nuc_acid('L')

    assert not is_nuc_acid(None)

    assert RNA2DNA('U') == 'T'
    assert DNA2RNA('T') == 'U'

    try:
        RNA2DNA('L')
    except Exception as e:
        assert 'Input sequence L is not valid nucleic acid sequence' in str(e)


def test_dna_bad_input():
    with pytest.raises(Exception):
        DNA2RNA('L')


def test_version():

    v = version()

    assert v is not None
