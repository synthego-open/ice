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


class ICEResult:

    def __init__(self):
        self.control_discordances = []
        self.edited_discordances = []
        self.mean_discord_before = None
        self.mean_discord_after = None
        self.aggregate = None
        self.contribs = None
        self.edit_eff = None
        self.ice_d = None
        self.r_squared = None
        self.hdr_percentage = None
        self.ko_score= None

    def to_json(self, guide_targets, warnings):

        results = {}
        if self.edit_eff is not None:
            results['ice'] = int(round(self.edit_eff * 100, 0))
        else:
            results['ice'] = None
        if self.mean_discord_before is not None:
            results['mean_discord_before'] = round(self.mean_discord_before, 3)
        else:
            results['mean_discord_before'] = None
        if self.mean_discord_after is not None:
            results['mean_discord_after'] = round(self.mean_discord_after, 3)
        else:
            results['mean_discord_after'] = None
        if self.ice_d is not None:
            results['ice_d'] = int(round(self.ice_d, 0))
        else:
            results['ice_d'] = None
        if self.r_squared is not None:
            results['rsq'] = round(self.r_squared, 2)
        else:
            results['rsq'] = None

        results['hdr_pct'] = self.hdr_percentage
        results['ko_score'] = self.ko_score
        results['notes'] = "; ".join(warnings)
        results['guides'] = [vars(g) for g in guide_targets]

        results['status'] = 'succeeded'

        return results
