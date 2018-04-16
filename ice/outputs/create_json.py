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


def write_individual_contribs(sanger_analysis, to_file=None):
    # human readable
    with open(to_file, 'w') as out_file:
        print('Relative contribution of each sequence (normalized)', file=out_file)
        print("--------------------------------------------------------", file=out_file)
        for entry in sanger_analysis.results.contribs:
            if entry.x_rel > 0.0:
                print('%02.04f \t %s \t %s' % (entry.x_rel, entry.summary,
                                               entry.human_readable_sequence(bp_after_cutsite=150)),
                      file=out_file)


def write_contribs_json(sanger_analysis, to_file):
    out_list = []
    for entry in sanger_analysis.results.contribs:
        if entry.x_rel > 0.0:
            abundance = round(entry.x_rel, 3)
            out_list.append({'rel_abundance': abundance,
                             'human_readable': entry.human_readable_sequence(),
                             'indel': entry.summary_json, 'wt': entry.wildtype})

    editing_outcomes = {}
    for i in sanger_analysis.results.aggregate:
        editing_outcomes[i] = round(sanger_analysis.results.aggregate[i] * 100, 2)

    output = {'contribs_list': out_list, 'editing_outcomes': editing_outcomes, 'multiplex': sanger_analysis.multiplex}

    with open(to_file, 'w') as out_file:
        json.dump(output, out_file)


def write_all_proposals_json(sanger_analysis, to_file):
    out_list = []
    for entry in sanger_analysis.results.contribs:
        # if entry.x_rel > 0.0:
        abundance = round(entry.x_rel, 3)
        wt = False
        if entry.summary == 0:
            wt = True
        out_list.append({'rel_abundance': abundance,
                         'human_readable': entry.human_readable_sequence(),
                         'indel': entry.summary, 'wt': wt})
    with open(to_file, 'w') as out_file:
        json.dump(out_list, out_file)
