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

import json

from ice.utility.misc import round_list


def generate_discordance_indel_files(sanger_analysis, ice_results, to_file=None):
    # generate dict/json
    data = {}
    discord_plot = {}
    discord_plot['bp'] = list(range(len(ice_results.control_discordances)))
    discord_plot['control_discord'] = round_list(ice_results.control_discordances, 3)
    discord_plot['edited_discord'] = round_list(ice_results.edited_discordances, 3)
    discord_plot['aln_start'] = sanger_analysis.alignment.ctrl2sample_coords(sanger_analysis.alignment_window[0])
    discord_plot['aln_end'] = sanger_analysis.alignment.ctrl2sample_coords(sanger_analysis.alignment_window[1])
    discord_plot['cut_site'] = sanger_analysis.guide_targets[0].cutsite
    discord_plot['inf_start'] = sanger_analysis.alignment.ctrl2sample_coords(sanger_analysis.inference_window[0])
    discord_plot['inf_end'] = sanger_analysis.alignment.ctrl2sample_coords(sanger_analysis.inference_window[1])
    data['discord_plot'] = discord_plot

    lowerbound = -2 * sanger_analysis.indel_max_size
    upperbound = 1 * sanger_analysis.indel_max_size
    editing_outcomes = {}
    for i in range(lowerbound, upperbound, 1):
        editing_outcomes[i] = round(sanger_analysis.results.aggregate[i] * 100, 2)
    data['editing_outcomes'] = editing_outcomes
    data['editing_eff'] = round(100 * sanger_analysis.results.edit_eff, 2)
    data['r_sq'] = round(sanger_analysis.results.r_squared, 3)

    with open(to_file, 'w') as f:
        json.dump(data, f, ensure_ascii=False)
