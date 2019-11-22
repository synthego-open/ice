from ice.classes.edit_proposal_creator import EditProposalCreator

import pytest


def test_single_cut():
    wt_basecalls = "GCACCTGTCCCCATAGAAAT"
    epc = EditProposalCreator(wt_basecalls, use_ctrl_trace=False)

    proposal = epc.single_cut_edit_proposal(10, "test", del_before=4)
    assert proposal.sequence == "GCACCTCCATAGAAAT"

    del_after_proposal = epc.single_cut_edit_proposal(10, "test", del_after=3)
    assert del_after_proposal.sequence == "GCACCTGTCCTAGAAAT"

    ins_proposal = epc.single_cut_edit_proposal(6, "test", insertion=2)
    assert ins_proposal.sequence.upper() == "GCACCTNNGTCCCCATAGAAAT"

    bad_config = epc.single_cut_edit_proposal(2, "bad", del_before=5)
    assert bad_config.sequence == "ACCTGTCCCCATAGAAAT"


def test_human_readable_sequence_single_edit_proposal():
    wt_basecalls = 'ACTGCCTGCACCTGTCCCCATAGAAATCCGCTTTACGGAGAGTCACTAGACTAGGACCCCCATAAATTTCACGACAGGGACTAGACCTTGACTGGATA'
    epc = EditProposalCreator(wt_basecalls, use_ctrl_trace=False)
    proposal = epc.single_cut_edit_proposal(27, 'g1', del_before=4)

    # test default padding
    seq = proposal.human_readable_sequence()
    assert seq == 'TGCCTGCACCTGTCCCCATAG----|CCGCTTTACGGAGAGTCACTAGACTAGGACCCCCATAAATTTCACGACAG'

    # test custom padding
    seq = proposal.human_readable_sequence(bp_before_cutsite=10, bp_after_cutsite=15)
    assert seq == 'CCATAG----|CCGCTTTACGGAGAG'

    # test passing in out of bounds padding values (should default to returning the entire proposed sequence)
    seq = proposal.human_readable_sequence(bp_before_cutsite=500, bp_after_cutsite=500)
    assert seq == 'ACTGCCTGCACCTGTCCCCATAG----|CCGCTTTACGGAGAGTCACTAGACTAGGACCCCCATAAATTTCACGACAGGGACTAGACCTTGACTGGAT'

    # test passing in inf values (should return the the entire proposed sequence)
    seq = proposal.human_readable_sequence(bp_before_cutsite=float('inf'), bp_after_cutsite=float('inf'))
    assert seq == 'ACTGCCTGCACCTGTCCCCATAG----|CCGCTTTACGGAGAGTCACTAGACTAGGACCCCCATAAATTTCACGACAGGGACTAGACCTTGACTGGAT'


def test_human_readable_sequence_multiplex_proposal():
    wt_basecalls = ('CCCATAAATCCGCTTTACAGTCACTAGACTAGGACCAATTTCACGACAGGGACTAGACCTTGCTGGATAGCACCTGTCCCCATAGAAATGGCTATGGA'
                    'AAGCCTTTGGGTTATTTGCGGCACCTGTCCCCATAGAAATGGCTATGGAAAGCCTTTGGGTTATTTGCG')
    epc = EditProposalCreator(wt_basecalls, use_ctrl_trace=False)

    # test close together dropout default padding
    close_together_dropout_proposal = epc.multiplex_proposal(30, 50, 'g1', 'g2', dropout=True)
    seq = close_together_dropout_proposal.human_readable_sequence()
    assert seq == 'AAATCCGCTTTACAGTCACTAGACT|--------------------|GACTAGACCTTGCTGGATAGCACCTGTCCC'

    # test far apart dropout default padding
    far_apart_dropout_proposal = epc.multiplex_proposal(30, 90, 'g1', 'g2', dropout=True)
    seq = far_apart_dropout_proposal.human_readable_sequence()
    assert seq == ('AAATCCGCTTTACAGTCACTAGACT|------------------------------------------------------------|GCTATGGAAAGC'
                   'CTTTGGGTTATTT')

    # test far apart insertion default padding
    far_apart_dropout_proposal = epc.multiplex_proposal(30, 90, 'g1', 'g2', dropout=False, cut1_ins=2, cut2_ins=4)
    seq = far_apart_dropout_proposal.human_readable_sequence()
    assert seq == ('ATCCGCTTTACAGTCACTAGACTnn|AGGACCAATTTCACGACAGGGACTAGACCTTGCTGGATAGCACCTGTCCCCATAGAAATGnnnn|GCTATGGA'
                   'AAGCCTTTGGGTTATTT')


def test_multiplex_cut():
    wt_basecalls = "GCACCTGTCCCCATAGAAATGGCTATGGAAAGCCTTTGGGTTATTTGCG"
    "GCACCTGTCC-CCATAGAAATGGCTATGGAAAGCCT-TTGGGTTATTTGCG"
    epc = EditProposalCreator(wt_basecalls, use_ctrl_trace=False)

    dropout_proposal = epc.multiplex_proposal(10, 35, "test1", "test2", dropout=True)
    assert dropout_proposal.sequence == "GCACCTGTCCTTGGGTTATTTGCGnnnnnnnnnnnnnnnnnnnnnnnnn"
    assert dropout_proposal.summary.startswith('-25:md')

    dropout_ins_proposal = epc.multiplex_proposal(10, 35, "test1", "test2", dropout=True, cut1_ins=2, cut2_ins=3)
    assert dropout_ins_proposal.sequence == "GCACCTGTCCnnTTGGGTTATTTGCGnnnnnnnnnnnnnnnnnnnnnnnnn"

    indep_cut = epc.multiplex_proposal(10, 35, "test1", "test2", cut1_del=(4, 0), cut2_del=(0, 3))
    assert indep_cut.sequence == "GCACCTCCATAGAAATGGCTATGGAAAGCCTGGTTATTTGCGnnnnnnn"
    assert indep_cut.summary.startswith('-7:m-4[test1]')

    indep_cut2 = epc.multiplex_proposal(10, 35, "test1", "test2", cut1_del=(-4, 0), cut2_del=(3, 0))
    assert indep_cut2.sequence == "GCACCTGTCCCCATAGAAATGGCTATGGAAAGTTGGGTTATTTGCGnnn"


def test_bad_multiplex():
    wt_basecalls = "GCACCTGTCCCCATAGAAATGGCTATGGAAAGCCTTTGGGTTATTTGCG"
    epc = EditProposalCreator(wt_basecalls, use_ctrl_trace=False)
    with pytest.raises(Exception):
        epc.multiplex_proposal(20, 10, "test1", "test2", dropout=True)


def test_wt_proposal():
    wt_basecalls = "GCACCTGTCCCCGGGGAAAT"
    epc = EditProposalCreator(wt_basecalls, use_ctrl_trace=False)

    proposal = epc.single_cut_edit_proposal(10, "test")
    assert proposal.sequence == "GCACCTGTCCCCGGGGAAAT"


def test_req_ctrl_trace():
    with pytest.raises(Exception):
        epc = EditProposalCreator("ATCG", use_ctrl_trace=True)


def test_generated_trace():
    seq = "ATNCG"

    epc = EditProposalCreator(seq, use_ctrl_trace=False)

    expected_trace = {'A': [], 'T': [], 'C': [], 'G': []}
    for base in seq:
        if base == 'N':
            for color in expected_trace.keys():
                expected_trace[color].append(0.25)
        else:
            for color in epc.base_order:
                if color == base:
                    expected_trace[color].append(0.97)
                else:
                    expected_trace[color].append(0.01)
    assert epc.wt_trace == expected_trace
