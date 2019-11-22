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
import xml.etree.ElementTree as ET

ET.register_namespace("", "http://www.w3.org/2000/svg")


def get_edited_base_px_positions(sanger_analysis, sample, changed_bases, start_px):
    changed_bases_px = []
    if sanger_analysis.donor_odn is not None:
        for changed_base in changed_bases:
            changed_base_left = sample.get_trough_position(changed_base, "left") - start_px
            changed_base_right = sample.get_trough_position(changed_base, "right") - start_px

            if changed_base_right > 0 and changed_base_right > 0:
                changed_bases_px.append((changed_base_left, changed_base_right))
    return changed_bases_px


def generate_chromatogram_data(data, start, end, cutsite, peaks, basecalls, basecoords):
    channels = [('DATA9', 'G'), ('DATA10', 'A'), ('DATA11', 'T'), ('DATA12', 'C')]
    output = {}
    output['channels'] = {}
    for (channel, base) in channels:
        # axis.plot(data[c[0]][start:end], color=colors[c[1]])
        output['channels'][base] = data[channel][start:end]
    output['cutsite'] = cutsite - start
    output['basecalls'] = []
    for idx, base in enumerate(basecalls):
        output['basecalls'].append([peaks[idx] - start, base, basecoords[0] + idx])
    return output


def reformat_to_json(data):
    list_of_dicts = []
    basecall_ptr = 0
    for idx, val in enumerate(data['channels']['G']):
        coord_vals = {}
        for base in 'ATCG':
            coord_vals[base] = data['channels'][base][idx]
        if basecall_ptr < len(data['basecalls']) and idx == data['basecalls'][basecall_ptr][0]:
            coord_vals['basecall'] = data['basecalls'][basecall_ptr][1]
            coord_vals['basecall_bp_coord'] = data['basecalls'][basecall_ptr][2]
            basecall_ptr += 1
        list_of_dicts.append(coord_vals)
    output = {}
    output['trace_data'] = list_of_dicts
    for key in ['cutsite', 'guide_start', 'guide_end', 'pam_start', 'pam_end']:
        if key in data:
            output[key] = data[key]
    return output


def generate_trace_files(sanger_analysis, to_file=None):
    ctrl_sample = sanger_analysis.control_sample
    edited_sample = sanger_analysis.edited_sample

    # convert base coords to trace coords
    upstream_in_bp = sanger_analysis.guide_targets[0].cutsite - 35
    downstream_in_bp = sanger_analysis.guide_targets[0].cutsite + 30
    ctrl_start = ctrl_sample.peak_locations[upstream_in_bp]
    ctrl_end = ctrl_sample.peak_locations[downstream_in_bp]

    ctrl_cutsite = ctrl_sample.get_trough_position(sanger_analysis.guide_targets[0].cutsite, "left")

    ctrl_peaks = ctrl_sample.peak_locations[upstream_in_bp:downstream_in_bp]
    ctrl_basecalls = ctrl_sample.primary_base_calls[upstream_in_bp:downstream_in_bp]
    ctrl_basecoords = (upstream_in_bp, downstream_in_bp)

    ctrl_data = generate_chromatogram_data(ctrl_sample.data, ctrl_start,
                                           ctrl_end, ctrl_cutsite, ctrl_peaks, ctrl_basecalls, ctrl_basecoords)
    guide_end_offset = len(sanger_analysis.gRNA_sequences[0]) - 1
    ctrl_data["guide_start"] = ctrl_sample.get_trough_position(
        sanger_analysis.guide_targets[0].cutsite - sanger_analysis.guide_targets[0].cut_offset, "left") - ctrl_start
    ctrl_data["guide_end"] = ctrl_sample.get_trough_position(sanger_analysis.guide_targets[0].cutsite + \
                                                             (guide_end_offset - sanger_analysis.guide_targets[
                                                                 0].cut_offset), "right") - ctrl_start

    if sanger_analysis.guide_targets[0].orientation == "fwd":
        pam_end = ctrl_sample.get_trough_position(sanger_analysis.guide_targets[0].cutsite + (
                guide_end_offset - sanger_analysis.guide_targets[0].cut_offset) + 3, "right")
        ctrl_data["pam_start"] = min(ctrl_data["guide_end"], pam_end - ctrl_start)
        ctrl_data["pam_end"] = max(ctrl_data["guide_end"], pam_end - ctrl_start)
    else:
        pam_start = ctrl_sample.get_trough_position(sanger_analysis.guide_targets[0].cutsite - 6, "left")
        ctrl_data["pam_start"] = min(pam_start - ctrl_start, ctrl_data["guide_start"])
        ctrl_data["pam_end"] = max(pam_start - ctrl_start, ctrl_data["guide_start"])

    edited_start_bp = sanger_analysis.alignment.ctrl2sample_coords(upstream_in_bp, closest=True)
    edited_end_bp = sanger_analysis.alignment.ctrl2sample_coords(downstream_in_bp, closest=True)
    edited_cutsite_bp = sanger_analysis.alignment.ctrl2sample_coords(sanger_analysis.guide_targets[0].cutsite,
                                                                     closest=True)
    edited_start = edited_sample.peak_locations[edited_start_bp]
    edited_end = edited_sample.peak_locations[edited_end_bp]

    edited_peaks = edited_sample.peak_locations[edited_start_bp:edited_end_bp]
    edited_basecalls = edited_sample.primary_base_calls[edited_start_bp:edited_end_bp]
    edited_basecoords = (edited_start_bp, edited_end_bp)

    if sanger_analysis.guide_targets[0].orientation == 'fwd':
        edited_cutsite = (edited_sample.peak_locations[edited_cutsite_bp]
                          + edited_sample.peak_locations[edited_cutsite_bp - 1]) / 2
    else:
        edited_cutsite = (edited_sample.peak_locations[edited_cutsite_bp]
                          + edited_sample.peak_locations[edited_cutsite_bp - 1]) / 2

    edited_data = generate_chromatogram_data(edited_sample.data, edited_start,
                                             edited_end, edited_cutsite, edited_peaks, edited_basecalls,
                                             edited_basecoords)

    ctrl_data['donor_odn_bases'] = get_edited_base_px_positions(sanger_analysis, ctrl_sample,
                                                                sanger_analysis.recombination_changed_bases,
                                                                ctrl_start)

    edited_bases = []
    for base in sanger_analysis.recombination_changed_bases:
        edited_bases.append(sanger_analysis.alignment.ctrl2sample_coords(base, closest=True))

    edited_data['donor_odn_bases'] = get_edited_base_px_positions(sanger_analysis, edited_sample,
                                                                  edited_bases,
                                                                  edited_start)

    out_json = {}

    out_json['ctrl_sample'] = reformat_to_json(ctrl_data)
    out_json['edited_sample'] = reformat_to_json(edited_data)

    with open(to_file, 'w') as f:
        json.dump(out_json, f, ensure_ascii=False)
