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

from ice.classes.proposal_base import ProposalBase


class EditProposal:
    """
    Stores the edit information for a single proposal.
    """
    BP_MIN_DEFAULT_PADDING = 25  # min default padding around cutsite for human readable sequence

    def __init__(self, wildtype=False):
        # if proposal is wildtype
        self.wildtype = wildtype

        # list of ProposalBases
        self.sequence_data = None

        # the index in the sequence_data @ where the cutsite is
        self.cutsite = None

        # to store the second cutsite, if there is any
        self.cutsite2 = None
        self.trace_data = None

        # -1, +3, -10, etc for indels
        # 'HDR' for recombination, etc
        self.bases_changed = None

        # more comprehensive code for bases_changed
        self.summary = None

        # details of bases changed in json format
        self.summary_json = None

        # absolute and relative coefficients from the regression
        self.x_rel = None
        self.x_abs = None

    @property
    def multiplex(self):
        """
        :return: bool indicating if proposal is multiplex (more than one proposed cutsite)
        """
        return self.cutsite and self.cutsite2

    @property
    def sequence(self):
        seq = ""
        for base in self.sequence_data:
            if base.base in 'ATCGNatcgn':
                seq += base.base
        return seq

    @property
    def min_cutsite(self):
        if self.cutsite2 is None:
            return self.cutsite
        return min(self.cutsite, self.cutsite2)

    @property
    def max_cutsite(self):
        if self.cutsite2 is None:
            return self.cutsite
        return max(self.cutsite, self.cutsite2)

    def position_at_cutsite(self, pos):
        """
        :param pos: position to check if at cutsite
        :return: bool indicating if pos is at a cutsite
        """
        return (self.cutsite and self.cutsite == pos) or (self.cutsite2 and self.cutsite2 == pos)

    def _calc_default_bp_after_cutsite(self):
        """
        Calculates default basepair after cutsite padding for human_readable_sequence
        :return: number of basepairs to pad
        """
        default_bp_after_cutsite = 50

        # trim some of extra bp after based on distance between cutsites
        if self.multiplex:
            dist_between_cutsites = self.cutsite2 - self.cutsite
            return max(default_bp_after_cutsite - dist_between_cutsites, self.BP_MIN_DEFAULT_PADDING)
        return default_bp_after_cutsite

    def _calc_default_bp_before_cutsite(self):
        """
        Calculates default basepair before cutsite padding for human_readable_sequence
        :return: number of basepairs to pad
        """
        return self.BP_MIN_DEFAULT_PADDING

    def human_readable_sequence(self, bp_before_cutsite=None, bp_after_cutsite=None):
        """
        Returns sequence around proposed edit site(s) showing indels.

        :param bp_before_cutsite: padding basepair count before min cutsite position
        :param bp_after_cutsite: padding basepair count after max cutsite position
        :return: sequence string around proposed cutsites with cutsites and indels shown
        """
        if bp_before_cutsite is None:
            bp_before_cutsite = self._calc_default_bp_before_cutsite()

        if bp_after_cutsite is None:
            bp_after_cutsite = self._calc_default_bp_after_cutsite()

        start_pos = max(self.min_cutsite - bp_before_cutsite + 1, 0)
        end_pos = min(self.max_cutsite + bp_after_cutsite + 1, len(self.sequence_data) - 1)
        seq = ''
        for seq_idx, base in enumerate(self.sequence_data[start_pos:end_pos], start_pos):
            if base.base_type is ProposalBase.WILD_TYPE:
                seq += base.base
            else:
                seq += base.base.lower()

            if self.position_at_cutsite(seq_idx):
                seq += '|'
        return seq
