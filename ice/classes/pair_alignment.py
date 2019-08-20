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

from io import StringIO
import json

from Bio import pairwise2, AlignIO


class PairAlignment:
    '''
    Uses a list of lists and two dictionaries to store data about a pairwise alignment
    '''

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.name1 = "control"
        self.name2 = "edited"
        self.alignment_pairs = None
        self.all_aligned_seqs = None
        self.aln_seqs = None

        self.sample_to_control = {}
        self.control_to_sample = {}
        self.all_aligned_clustal = None
        self.aligned_clustal = None
        self.last_aligned_pair_idx = None

        self.aln_clustal = None

    # def convert to clustal format
    # def write alignment to file

    def __str__(self):
        return self.all_aligned_seqs[0] + "\n" + self.all_aligned_seqs[1]

    @staticmethod
    def write_aln(aln_clustal, to_file=None):
        if to_file is not None:
            with open(to_file, 'w') as out_file:
                print(aln_clustal, file=out_file)

    @staticmethod
    def write_json(aln_seqs, to_file):
        if aln_seqs:
            out_dict = {'control': aln_seqs[0], 'edited': aln_seqs[1]}
            with open(to_file, 'w') as f:
                json.dump(out_dict, f)

    def align_list_to_clustal(self, aln, name1, name2):
        my_aln = ">{}\n".format(name1) + aln[0] + "\n" + ">{}\n".format(name2) + aln[1]
        f = StringIO(my_aln)
        aln_objs = list(AlignIO.parse(f, "fasta"))
        alignment_txt = aln_objs[0].format("clustal").split('\n', 2)[2]
        return alignment_txt

    def align_all(self):
        seq1 = self.seq1
        seq2 = self.seq2

        match_bonus = 2
        #(aseq, bseq, match, mismatch, gap_open, gap_extend)
        alignments = pairwise2.align.localms(seq1,
                                             seq2,
                                             match_bonus, -1, -3, -1)
        aln = alignments[0]
        self.all_aligned_seqs = (aln[0], aln[1])
        self.all_aligned_clustal = self.align_list_to_clustal(aln, "control", "edited")

    def ctrl2sample_coords(self, coord, closest=False):
        if self.alignment_pairs is None:
            raise Exception("alignment has not been performed and alignment_pairs is None")
        elif coord is None:
            raise Exception("ctrl coordinate to convert cannot be None")

        if self.control_to_sample[coord] is None and closest is True:
            # tries to find a sample coordinate to the left in the alignment
            pos = coord
            sample_coord = None
            while pos > -1 and sample_coord == None:
                sample_coord = self.control_to_sample[pos]
                pos -= 1
            return sample_coord
        else:
            return self.control_to_sample[coord]

    @property
    def has_alignment(self):
        return not not self.alignment_pairs

    def align_with_window(self, alignment_window):
        """
        This function uses a subsequence in seq1 to align to seq2.
        The remaining bases are then padded onto the aln

        Example:
        seq1 = AATGTAATGATAG
        seq2 = AATGTATGATAG
        window is 0, 3, so the sequence used for alignment is ATGT
        After initial alignment:
        AATG
        AATGTATGATAG
        After alignment padding:
        AATGTAATGATAG
        AATGTATGATAG

        Whereas an alignment of the full sequences would result in
        AATGTAATGATAG
        AATGTA-TGATAG

        :param alignment_window: specifies subsequence of seq1 to align to seq2
        :param seq1: control sequence
        :param seq2: edited sequence
        :return: True if completed
        """
        seq1 = self.seq1
        seq2 = self.seq2
        aw = alignment_window
        window_size = aw[1] - aw[0]

        match_bonus = 2
        # align windowed seq1 to seq2
        alignments = pairwise2.align.localms(seq1[aw[0]:aw[1]],
                                             seq2,
                                             match_bonus, -1, -2, -1)
        # a, b, c, d  scoring parameters for alignments
        # a = points added for match
        # b = points added for mismatch
        # c = gap open penalty
        # d = gap extension penalty

        # Initial Alignment
        try:
            # get the top alignment
            aln = alignments[0]
        except Exception as e:
            self.aln_clustal = "No alignment found\nCtrl {} to {}:{}\nEdited:{}\n".format(
                aw[0], aw[1], seq1, seq2)
            debug_txt = " window: %s to %s " % (aw[0], aw[1])
            return False, ('No alignment found ' + str(e) + debug_txt)

        aln_score = aln[2]
        aln_score_normalized = aln_score / (window_size * match_bonus) * 100
        self.aln_clustal = self.align_list_to_clustal(aln, "control_aln_window", "edited")

        # we'd like to have at least around 40% of bases matching)
        if 50 > aln_score_normalized:
            return False, "Poor alignment upstream of cutsite: %0.2f percent of full score" % (aln_score_normalized)

        self.aln_seqs = (aln[0], aln[1])

        src_alignment = aln[0]
        dst_alignment = aln[1]

        last_aligned_ctrl_base = None
        index = len(src_alignment) - 1

        # find the last aligned base (closest to right end), then we will fill in the alignment downstream
        while last_aligned_ctrl_base is None and index > -1:
            if src_alignment[index] == '-':
                index -= 1
            else:
                last_aligned_ctrl_base = index

        # use the alignment to generate a list of tuples that matches ref bases with sample bases
        sample_index = 0
        reference_index = aw[0]

        alignment_pairs = []
        first_alned_ctrl_base = None

        for aln_idx, sample_base in enumerate(dst_alignment):
            r = None
            s = None
            ref_base = src_alignment[aln_idx]
            if sample_base != "-":
                s = sample_index
                sample_index += 1

            if ref_base != "-":
                if first_alned_ctrl_base is None:
                    first_alned_ctrl_base = aln_idx
                r = reference_index
                reference_index += 1

            # alignment padding
            # use the aligned segment to get the downstream ctrl subsequence that matches up with the inference window
            elif aw[1] <= reference_index < len(seq1):
                r = reference_index
                reference_index += 1

            alignment_pairs.append([r, s])
            self.sample_to_control[s] = r
            self.control_to_sample[r] = s

        aln_idx = first_alned_ctrl_base
        reference_index = alignment_pairs[aln_idx][0]
        aln_idx -= 1
        reference_index -= 1
        while aln_idx > -1 and reference_index > -1:
            alignment_pairs[aln_idx][0] = reference_index
            aln_idx -= 1
            reference_index -= 1

        self.alignment_pairs = alignment_pairs

        # now we traverse backwards and find the last aligned bases
        for aln_idx, aln_pair in reversed(list(enumerate(self.alignment_pairs))):
            ctrl_idx = aln_pair[0]
            edit_idx = aln_pair[1]
            if ctrl_idx is not None and edit_idx is not None:
                self.last_aligned_pair_idx = aln_idx
                break

        return True, "Aln succeeded"


class DonorAlignment(PairAlignment):
    """
    Adds functionality that is useful in the specific case of pair alignment between control sequence and donor sequence
    """
    MATCH_BONUS = 2

    def __init__(self, control_seq, donor_seq):
        super(DonorAlignment, self).__init__(control_seq, donor_seq)
        self.align_ssodn()
        self.hdr_indel_size = self._calc_hdr_indel_size()  # how much the donor HDR would change control seq length

    @property
    def aligned_control_seq(self):
        return self.all_aligned_seqs[0]

    @property
    def aligned_donor_seq(self):
        return self.all_aligned_seqs[1]

    def _calc_hdr_indel_size(self):
        insert_total_len = self.aligned_control_seq.strip('-').count('-')
        deletion_total_len = self.aligned_donor_seq.strip('-').count('-')
        return insert_total_len - deletion_total_len

    @property
    def control_seq(self):
        return self.seq1

    @property
    def donor_seq(self):
        return self.seq2

    def align_ssodn(self):
        print('starting to align ssodn')
        alignments = pairwise2.align.localms(self.control_seq, self.donor_seq, self.MATCH_BONUS, -1, -6, -0.5)
        alignment = alignments[0]
        self.all_aligned_seqs = (alignment[0], alignment[1])
        self.all_aligned_clustal = self.align_list_to_clustal(alignment, 'control', 'donor')
