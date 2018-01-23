"""
Copyright 2018 Synthego Corporation All Rights Reserved

The Synthego ICE software was developed at Synthego Corporation.

Permission to use, copy, modify and distribute any part of Synthego ICE for
educational, research and non-profit purposes, without fee, and without a
written agreement is hereby granted, provided that the above copyright notice,
this paragraph and the following paragraphs appear in all copies.

Those desiring to incorporate this Synthego ICE software into commercial
products or use for commercial purposes should contact Synthego support at Ph:
(888) 611-6883 ext:1, E-MAIL: support@synthego.com.

In no event shall Synthego Corporation be liable to any party for direct,
indirect, special, incidental, or consequential damages, including lost
profits, arising out of the use of Synthego ICE, even if Synthego Corporation
has been advised of the possibility of such damage.

The Synthego ICE tool provided herein is on an "as is" basis, and the Synthego
Corporation has no obligation to provide maintenance, support, updates,
enhancements, or modifications. The Synthego Corporation makes no
representations and extends no warranties of any kind, either implied or
express, including, but not limited to, the implied warranties of
merchantability or fitness for a particular purpose, or that the use of
Synthego ICE will not infringe any patent, trademark or other rights.
"""

from ice.classes.edit_proposal import EditProposal
from ice.classes.pair_alignment import PairAlignment
from ice.classes.proposal_base import ProposalBase
from ice.utility.sequence import reverse_complement


class EditProposalCreator:
    MIN_HOMOLOGY_ARM = 15

    def __init__(self, wt_basecalls, use_ctrl_trace=True, sanger_object=None):
        self.sanger_object = sanger_object
        self.base_order = 'ATCG'
        self.use_ctrl_trace = use_ctrl_trace
        self.wt_basecalls = wt_basecalls
        self.wt_trace = None

        if self.use_ctrl_trace:
            if self.sanger_object is None:
                raise Exception('If you wish to use a ctrl Trace, you must supply a SangerObject')
            self.base_order = self.sanger_object.base_order
            self.wt_trace = self.trace_values
        else:
            proposal_trace = {'A': [], 'T': [], 'C': [], 'G': []}
            for idx, base in enumerate(self.wt_basecalls):
                if base == 'N':
                    for base_color in self.base_order:
                        proposal_trace[base_color].append(0.25)
                else:
                    for base_color in self.base_order:
                        if base_color == base:
                            proposal_trace[base_color].append(0.97)
                        else:
                            proposal_trace[base_color].append(0.01)
            self.wt_trace = proposal_trace

    @property
    def trace_values(self):
        return self.sanger_object.get_peak_values()

    def single_cut_edit_proposal(self, cutsite, label, del_before=0, del_after=0, insertion=0):
        cutsite = cutsite - 1
        proposal_bases = []
        proposal_trace = []

        # deletion case
        if del_before > 0 or del_after > 0:
            deleted_bases = [cutsite - i for i in range(del_before)] + [cutsite + i + 1 for i in range(del_after)]
            for idx, base in enumerate(self.wt_basecalls):
                if idx in deleted_bases:
                    proposal_base = ProposalBase('-', ProposalBase.DELETION, idx)
                else:
                    proposal_base = ProposalBase(base, ProposalBase.WILD_TYPE, idx)
                    for base_index, base_color in enumerate(self.base_order):
                        proposal_trace.append(self.wt_trace[base_color][idx])
                proposal_bases.append(proposal_base)
            ep = EditProposal()
            ep.sequence_data = proposal_bases
            ep.cutsite = cutsite
            ep.bases_changed = -(del_before + del_after)
            ep.summary = "{}[{}]".format(-(del_before + del_after), label)
            ep.trace_data = proposal_trace
        # insertion case
        elif insertion > 0:
            for idx, base in enumerate(self.wt_basecalls):
                if idx == cutsite:

                    proposal_bases.append(ProposalBase(base, ProposalBase.WILD_TYPE, idx))
                    for base_index, base_color in enumerate(self.base_order):
                        proposal_trace.append(self.wt_trace[base_color][idx])
                    for i in range(insertion):

                        proposal_bases.append(ProposalBase('n', ProposalBase.INSERTION, idx))
                        for base_index, base in enumerate(self.base_order):
                            proposal_trace.append(0.25)
                else:
                    proposal_bases.append(ProposalBase(base, ProposalBase.WILD_TYPE, idx))
                    for base_index, base_color in enumerate(self.base_order):
                        proposal_trace.append(self.wt_trace[base_color][idx])
            ep = EditProposal()
            ep.sequence_data = proposal_bases
            ep.cutsite = cutsite
            ep.bases_changed = insertion
            ep.summary = "{}[{}]".format(insertion, label)
            ep.trace_data = proposal_trace
        # wildtype case
        else:
            for idx, base in enumerate(self.wt_basecalls):
                proposal_base = ProposalBase(base, ProposalBase.WILD_TYPE, idx)
                proposal_bases.append(proposal_base)
                for base_index, base_color in enumerate(self.base_order):
                    proposal_trace.append(self.wt_trace[base_color][idx])
            ep = EditProposal()
            ep.sequence_data = proposal_bases
            ep.cutsite = cutsite
            ep.bases_changed = 0
            ep.summary = "0[{}]".format(label)
            ep.trace_data = proposal_trace
        return ep

    def multiplex_proposal(self, cutsite1, cutsite2, label1, label2,
                           cut1_del=(0, 0), cut1_ins=0, cut2_del=(0, 0), cut2_ins=0,
                           dropout=False
                           ):
        '''
        this method will be used for two guides cutting at the same time
        Note that both sites MUST be edited with either a deletion or an insertion.
        Cases where one guide did not result in an edit should be handled by the single edit method
        :param cutsite1: must come before cutsite2
        :param cutsite2:
        :param cut1_del:
        :param cut1_ins:
        :param cut2_del:
        :param cut2_ins:
        :return:
        '''
        cutsite1 = cutsite1 - 1
        cutsite2 = cutsite2 - 1
        proposal_bases = []
        proposal_trace = []
        if cutsite2 <= cutsite1:
            raise Exception("cutsite1 must come before cutsite2 values are ({}, {})".format(cutsite1, cutsite2))
        summary_code = "m"
        deleted_bases = []

        if dropout:
            summary_code = "md"
            for i in range(cutsite1 + 1, cutsite2 + 1):
                deleted_bases.append(i)
        #            deleted_bases = [cutsite-i for i in range(del_before)] + [cutsite+i+1 for i in range(del_after)]

        # both cutsites results in deletions
        # todo, there are no safeguards for if the deletions for one cutsite exceed the boundaries of the other cutsite
        # that info matters a lot for the summary
        # for cut1 allow all deletions before
        # for cut2 allow all deletions after
        # for cut1 deletions after, stop if cut2 reached
        # for cut2 deletions before, stop if cut1 reached
        if cut1_del != (0, 0) and cut2_del != (0, 0):
            # todo, this should also handle the straight dropout case w/ no extra deletions
            if dropout:
                cut1_del_after = 0
                cut2_del_before = 0
            else:
                cut1_del_after = cut1_del[1]
                cut2_del_before = cut2_del[0]
            cut1 = (cutsite1, cut1_del[0], cut1_del_after)
            cut2 = (cutsite2, cut2_del_before, cut2_del[1])
            for cutsite, del_before, del_after in [cut1, cut2]:
                deleted_bases += [cutsite - i for i in range(del_before)] + [cutsite + i + 1 for i in range(del_after)]

            for idx, base in enumerate(self.wt_basecalls):
                if idx in deleted_bases:
                    proposal_base = ProposalBase('-', ProposalBase.DELETION, idx)
                else:
                    proposal_base = ProposalBase(base, ProposalBase.WILD_TYPE, idx)
                    for base_index, base_color in enumerate(self.base_order):
                        proposal_trace.append(self.wt_trace[base_color][idx])
                proposal_bases.append(proposal_base)
            ep = EditProposal()
            ep.sequence_data = proposal_bases
            ep.cutsite = cutsite1
            ep.cutsite2 = cutsite2
            total_deleted = cut1_del[0] + cut1_del_after + cut2_del_before + cut2_del[1]
            total_deleted += -(cutsite2 - cutsite1)
            ep.bases_changed = total_deleted
            ep.summary = "{}:{}-{}[{}],-{}[{}]".format(total_deleted, summary_code,
                                                       cut1_del[0] + cut1_del_after,
                                                       label1,
                                                       cut2_del_before + cut2_del[1],
                                                       label2
                                                       )
            ep.trace_data = proposal_trace

        # insertion case
        elif cut1_ins != 0 or cut2_ins != 0:
            if dropout:
                cut2_ins = 0

            cut1 = (cutsite1, cut1_ins)
            cut2 = (cutsite2, cut2_ins)

            for idx, base in enumerate(self.wt_basecalls):

                if idx in deleted_bases:
                    proposal_bases.append(ProposalBase('-', ProposalBase.DELETION, idx))
                # if base is in range where we need to do insertion
                elif idx in [cutsite1, cutsite2]:
                    for cutsite, insertion_length in [cut1, cut2]:
                        if cutsite == idx:
                            # if cutsite != cutsite2:
                            proposal_bases.append(ProposalBase(base, ProposalBase.WILD_TYPE, idx))
                            for base_index, base_color in enumerate(self.base_order):
                                proposal_trace.append(self.wt_trace[base_color][idx])
                            # elif dropout:
                            # proposal_bases.append(ProposalBase('-', ProposalBase.DELETION, idx))
                            for i in range(insertion_length):
                                proposal_bases.append(ProposalBase('n', ProposalBase.INSERTION, idx))
                                for base_index, base in enumerate(self.base_order):
                                    proposal_trace.append(0.25)
                else:
                    proposal_bases.append(ProposalBase(base, ProposalBase.WILD_TYPE, idx))
                    for base_index, base_color in enumerate(self.base_order):
                        proposal_trace.append(self.wt_trace[base_color][idx])
            ep = EditProposal()
            ep.sequence_data = proposal_bases
            ep.cutsite = cutsite1 + cut1_ins
            ep.cutsite2 = cutsite1 + cut1_ins + (cutsite2 - cutsite1) + cut2_ins
            ep.bases_changed = cut1_ins + cut2_ins
            if dropout:
                ep.bases_changed -= (cutsite2 - cutsite1)
            ep.summary = "{}:{}+{}[{}],+{}[{}]".format(
                ep.bases_changed,
                summary_code,
                cut1_ins,
                label1,
                cut2_ins,
                label2
            )
            ep.trace_data = proposal_trace
        # the intervening sequence is dropped out and no bases inserted or deleted
        else:
            if dropout:
                for idx, base in enumerate(self.wt_basecalls):
                    if idx in deleted_bases:
                        proposal_base = ProposalBase('-', ProposalBase.DELETION, idx)
                    else:
                        proposal_base = ProposalBase(base, ProposalBase.WILD_TYPE, idx)
                        for base_index, base_color in enumerate(self.base_order):
                            proposal_trace.append(self.wt_trace[base_color][idx])
                    proposal_bases.append(proposal_base)
                ep = EditProposal()
                ep.sequence_data = proposal_bases
                ep.cutsite = cutsite1
                ep.cutsite2 = cutsite2

                total_deleted = -(cutsite2 - cutsite1)
                ep.bases_changed = total_deleted
                ep.summary = "{}:{}-0[{}],-0[{}]".format(total_deleted, summary_code, label1, label2)
                ep.trace_data = proposal_trace
            else:
                # wild type case, use the single edit to model this case
                return None
        return ep

    def recombined_outcome_sequence(self, alignment, donor_sequence):
        proposal_bases = []
        proposal_trace = []
        wt_seq = list(alignment[0])
        aligned_odn = alignment[1]
        changed_bases = []
        odn_start_pos = None
        odn_idx = 0
        ctrl_seq_idx = 0
        for idx, base in enumerate(aligned_odn):
            if base is not '-':
                odn_idx += 1
            # within the region affected by the donor sequence
            if odn_idx > 0 and odn_idx <= len(donor_sequence) - 2:
                if base != wt_seq[idx]:
                    proposal_bases.append(ProposalBase(base, ProposalBase.HDR, ctrl_seq_idx))
                    changed_bases.append(ctrl_seq_idx)
                    for base_index, base_color in enumerate(self.base_order):
                        if base_color == base:
                            proposal_trace.append(0.97)
                        else:
                            proposal_trace.append(0.01)
                else:
                    # same base as wt
                    proposal_bases.append(ProposalBase(base, ProposalBase.WILD_TYPE, ctrl_seq_idx))
                    for base_index, base_color in enumerate(self.base_order):
                        proposal_trace.append(self.wt_trace[base_color][ctrl_seq_idx])
                if odn_start_pos is None:
                    odn_start_pos = idx
            # all other positions outside of the region affected by HDR
            else:
                proposal_bases.append(ProposalBase(wt_seq[idx], ProposalBase.WILD_TYPE, ctrl_seq_idx))
                for base_index, base_color in enumerate(self.base_order):
                    proposal_trace.append(self.wt_trace[base_color][ctrl_seq_idx])
            # HDR may also cause deletions or insertions, so we need to keep track of ctrl_seq_idx separately
            if wt_seq[idx] != '-':
                ctrl_seq_idx += 1
        return proposal_bases, proposal_trace, changed_bases, odn_start_pos

    def homologous_recombination_proposal(self, cutsite, donor_sequence):
        cutsite = cutsite - 1
        # check orientation and verify at least 15nt matching on both ends of ODN

        if reverse_complement(donor_sequence[:self.__class__.MIN_HOMOLOGY_ARM]) in self.wt_basecalls:
            donor_sequence = reverse_complement(donor_sequence)

        if donor_sequence[:self.__class__.MIN_HOMOLOGY_ARM] in self.wt_basecalls and \
                donor_sequence[-self.__class__.MIN_HOMOLOGY_ARM:] in self.wt_basecalls:

            # todo confirm that cutsite is inside the HDR region
            donor_alignment = PairAlignment(self.wt_basecalls, donor_sequence)
            donor_alignment.align_ssodn()

            proposal_bases, proposal_trace, changed_bases, odn_start_pos = self.recombined_outcome_sequence(
                donor_alignment.all_aligned_seqs, donor_sequence)

            ep = EditProposal()
            ep.sequence_data = proposal_bases
            ep.trace_data = proposal_trace
            ep.cutsite = cutsite
            ep.summary = ProposalBase.HDR
            ep.bases_changed = ProposalBase.HDR
            return ep, changed_bases, odn_start_pos, donor_alignment
        else:
            raise Exception(
                "Homology arms of length {} not found in control sequence".format(self.__class__.MIN_HOMOLOGY_ARM))
