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


class ProposalBase:
    WILD_TYPE = 'wt'
    INSERTION = 'insertion'
    DELETION = 'deletion'
    BASE_EDIT = 'base_edit'
    HDR = 'hdr'
    BASE_TYPE = [WILD_TYPE, INSERTION, DELETION, BASE_EDIT, HDR]

    def __init__(self, base, base_type, wt_coord):
        self.base = base
        assert base_type in ProposalBase.BASE_TYPE
        self.base_type = base_type
        self.wt_coord = wt_coord
