from ice.classes.edit_proposal_creator import  EditProposalCreator

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


def test_pretty_print():
    wt_basecalls = "GCACCTGTCCCCATAGAAAT"
    epc = EditProposalCreator(wt_basecalls, use_ctrl_trace=False)

    with pytest.raises(Exception) as e:
        proposal = epc.single_cut_edit_proposal(10, "test", del_before=4)
        proposal.human_readable_sequence()

    proposal = epc.single_cut_edit_proposal(10, "test", del_before=4)
    seq = proposal.human_readable_sequence(bp_before_cutsite=9, bp_after_cutsite=3)

    assert seq == "CACCT----|CCA"

def test_pretty_print2():
    wt_basecalls = "GCACCTGTCCCCATAGAAAT"
    epc = EditProposalCreator(wt_basecalls, use_ctrl_trace=False)

    with pytest.raises(Exception) as e:
        proposal = epc.single_cut_edit_proposal(10, "test", del_before=4)
        proposal.human_readable_sequence()

    proposal = epc.single_cut_edit_proposal(10, "test", del_before=4)
    seq = proposal.human_readable_sequence(bp_before_cutsite=None, bp_after_cutsite=10)

    assert seq == "GCACCT----|CCATAGAAAT"


def test_multiplex_cut():
    wt_basecalls = "GCACCTGTCCCCATAGAAATGGCTATGGAAAGCCTTTGGGTTATTTGCG"
    "GCACCTGTCC-CCATAGAAATGGCTATGGAAAGCCT-TTGGGTTATTTGCG"
    epc = EditProposalCreator(wt_basecalls, use_ctrl_trace=False)

    dropout_proposal = epc.multiplex_proposal(10, 35, "test1", "test2", dropout=True)
    assert dropout_proposal.sequence == "GCACCTGTCCTTGGGTTATTTGCG"

    dropout_ins_proposal = epc.multiplex_proposal(10, 35, "test1", "test2", dropout=True, cut1_ins=2, cut2_ins=3)
    assert dropout_ins_proposal.sequence == "GCACCTGTCCnnTTGGGTTATTTGCG"

    indep_cut = epc.multiplex_proposal(10, 35, "test1", "test2", cut1_del=(4, 0), cut2_del=(0, 3))
    assert indep_cut.sequence == "GCACCTCCATAGAAATGGCTATGGAAAGCCTGGTTATTTGCG"

    indep_cut2 = epc.multiplex_proposal(10, 35, "test1", "test2", cut1_del=(-4, 0), cut2_del=(3, 0))
    assert indep_cut2.sequence == "GCACCTGTCCCCATAGAAATGGCTATGGAAAGTTGGGTTATTTGCG"

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

    expected_trace = {'A':[], 'T':[], 'C':[], 'G':[]}
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
