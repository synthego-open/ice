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
    """This class is to store the data for one edit proposal"""

    def __init__(self):
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

        # absolute and relative coefficients from the regression
        self.x_rel = None
        self.x_abs = None

    @property
    def sequence(self):
        seq = ""
        for base in self.sequence_data:
            if base.base in 'ATCGNatcgn':
                seq += base.base
        return seq

    def human_readable_sequence(self, bp_before_cutsite=25, bp_after_cutsite=50):
        """

        :return:
        """
        seq = ""
        if bp_before_cutsite is None:
            for idx, base in enumerate(self.sequence_data):
                if base.base_type is ProposalBase.WILD_TYPE:
                    seq += base.base
                else:
                    seq += base.base.lower()
                if idx == self.cutsite:
                    seq += "|"
        else:
            if bp_before_cutsite - 1 > self.cutsite:
                err = "bp_before_cutsite of {} requested exceeds sequence length".format(bp_before_cutsite)
                raise Exception(err)
            for idx, base in enumerate(self.sequence_data):
                if self.cutsite - idx > bp_before_cutsite - 1:
                    continue
                if idx - self.cutsite > bp_after_cutsite:
                    continue
                if base.base_type is ProposalBase.WILD_TYPE:
                    seq += base.base
                else:
                    seq += base.base.lower()
                if idx == self.cutsite or idx == self.cutsite2:
                    seq += "|"
        return seq
