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

import os
import numpy as np
from Bio import SeqIO

from ice.utility.misc import running_mean
from ice.utility.misc import find_regions_greater_than


class SangerObject:
    """

    """

    PRIMARY_BASE_CALL = 'PBAS1'
    BASE_ORDER = 'FWO_1'
    PEAK_LOCATIONS = 'PLOC1'
    CHANNELS = ['DATA9', 'DATA10', 'DATA11', 'DATA12']

    def __init__(self, verbose=False):
        self.data = None
        # self.left_trim_bases = 5
        # self.right_trim_bases = 5
        # self.alignment_offset = 0
        self.external_base_calls = None
        self.path = None
        self.basename = None

    def estimate_quality_scores(self):
        '''
        a good enough way to estimate quality scores
        As ICE uses the quality scores for determining windows for alignment, we just need a rough estimate
        of base quality.  If the primary color is 90% of the signal, the phred score is 50% of the maximum.
        If the primary color is 80% of the total signal at that position, the phred score is 0% of the maximum.
        :return: list of phred scores
        '''

        MAX_PHRED = 60
        peak_values = self.get_peak_values()
        phred_scores = []
        for idx, base in enumerate(peak_values['C']):
            sum_values = 0
            max_value = 0
            for color in 'ACGT':
                sum_values += peak_values[color][idx]
                if peak_values[color][idx] > max_value:
                    max_value = peak_values[color][idx]

            if sum_values < 100 or sum_values > 5000:
                phred_score = 0
            else:
                max_value_percent = max_value / sum_values * 100
                
                # linear relationship of PHRED score 0 to 60 and max_value_percent 80 to 100
                phred_score = max(3 * max_value_percent - 240, 0)
            #check if absolute intensities are reasonable

            phred_scores.append(phred_score)
        return phred_scores

    def initialize_from_path(self, ab1_file_path, *args, **kwargs):
        '''


        :param ab1_file_path:
        :param args:
        :param kwargs:
        :return:

        'FWO_1': 'GATC' contains the base order
        DATA9: G, DATA10: A, DATA11: T, DATA12: C
        good documentation of what all the data channels are here:
        https://github.com/biopython/biopython/blob/master/Bio/SeqIO/AbiIO.py
        '''
        self.path = ab1_file_path
        self.basename = os.path.basename(ab1_file_path)
        record = SeqIO.read(ab1_file_path, 'abi')
        
        # In Biopython 1.74, some strings are converted to bytes in py3
        # Convert strings with bytes-object back to regular string
        traces_ = record.annotations['abif_raw']
        traces  = {}
        for k, v in traces_.items():
            if isinstance(v, bytes):
                v = v.decode()
            else:
                v = v
            traces[k] = v
        self.data = traces

        phred_scores = []
        for c in traces['PCON2']:
            phred_scores.append(ord(c))
        self.phred_scores = phred_scores

        override_quality_score = kwargs.get('override_quality_score', False)
        if np.count_nonzero(phred_scores) == 0 or override_quality_score:
            self.phred_scores = self.estimate_quality_scores()

    def find_alignable_window(self, window_size=30, QUAL_CUTOFF=40):
        '''
        This code finds the optimal window for aligning the edited trace to the reference sequence
        1) It looks at the phred scores and calculates a moving average
        2) It then selects the largest window where the moving average is greater than QUAL_CUTOFF
        3) This segment and other data is then returned in a dictionary
        all coordinates are absolute with respect to the sanger object
        '''

        # TODO: are there any edge cases where an edit is so good that the alignable window would extend into
        # regions where there is mismatch between sample and control?

        phreds = np.asarray(self.phred_scores)
        # find a running mean
        window_size = window_size
        windowed_phred = running_mean(phreds, window_size)

        # find regions with windowed mean > QUAL_CUTOFF
        y = find_regions_greater_than(windowed_phred, QUAL_CUTOFF)



        max_window_size = 0
        max_window = None
        for region in y:
            w_size = region[1] - region[0]
            if w_size > max_window_size:
                max_window_size = w_size
                max_window = (region[0], region[1])

        return {
            'regions': y,
            'max_window': max_window,
            'windowed_phred': windowed_phred
        }

    @property
    def primary_base_calls(self):
        return self.data[SangerObject.PRIMARY_BASE_CALL]

    @property
    def aligned_bases(self):
        return self.primary_base_calls

    @property
    def base_order(self):
        return self.data[SangerObject.BASE_ORDER]

    @property
    def peak_locations(self):
        """ Returns the peak locations for the final sequence

        :return:
        """
        return list(self.data[SangerObject.PEAK_LOCATIONS])

    def get_peak_values(self):
        output_dict = {}
        for n, channel in enumerate(SangerObject.CHANNELS):
            peak_data = np.array(self.data[channel])
            # peak_values = peak_data[self.peak_locations]
            peak_values = peak_data[list(self.data[SangerObject.PEAK_LOCATIONS])]
            output_dict[self.base_order[n]] = peak_values

        return output_dict

    def get_trough_position(self, loc1, side="left"):
        if side == "left":
            loc2 = loc1 - 1
        else:
            loc2 = loc1 + 1
        return (self.peak_locations[loc1]
                + self.peak_locations[loc2]) / 2

    @classmethod
    def trace_to_base_calls(cls, trace, base_order):
        BASE_CALL_CUTOFF = 0.75
        assert len(trace) % 4 == 0
        seq_len = int(len(trace) / 4)
        seq = ""
        for base_idx in range(seq_len):
            slice_values = trace[base_idx * 4: base_idx * 4 + 4]
            slice_max = max(slice_values)
            if slice_max > BASE_CALL_CUTOFF:
                max_idx = slice_values.index(slice_max)
                base = base_order[max_idx]
                seq += base
            else:
                seq += 'N'
        return seq
