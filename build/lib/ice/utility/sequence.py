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

DNA_COMPLEMENT = {
    'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
    'a': 't', 't': 'a', 'c': 'g', 'g': 'c'
}

RT_MAP = {
    'U': 'A', 'A': 'T', 'C': 'G', 'G': 'C',
    'u': 'a', 'a': 't', 'c': 'g', 'g': 'c',
}

VALID_RNA = ['A', 'C', 'G', 'U', 'a', 'c', 'g', 'u']

VALID_NUC_ACID = ['A', 'C', 'G', 'T', 'U']

TRANS_MAP = {v: k for k, v in RT_MAP.items()}


def reverse_transcribe(input_seq):
    output = ''
    for c in input_seq:
        output += RT_MAP.get(c, 'x')
    return output


def transcribe(input_seq):
    output = ''
    for c in input_seq:
        output += TRANS_MAP.get(c, 'x')
    return output


def reverse_complement(input_seq):
    output = ''
    for c in reversed(input_seq):
        output += DNA_COMPLEMENT.get(c, 'x')
    return output


def complement(input_seq):
    output = ''
    for c in input_seq:
        output += DNA_COMPLEMENT.get(c, 'x')
    return output


def is_nuc_acid(input_seq):
    try:
        return all(nuc in VALID_NUC_ACID for nuc in input_seq.upper())
    except Exception:
        return False


def RNA2DNA(input_seq):
    if is_nuc_acid(input_seq):
        return input_seq.replace('U', 'T').replace('u', 't')
    else:
        raise Exception("Input sequence {} is not valid nucleic acid sequence".format(input_seq))


def DNA2RNA(input_seq):
    if is_nuc_acid(input_seq):
        return input_seq.replace('T', 'U').replace('t', 'u')
    else:
        raise Exception("Input sequence {} is not valid nucleic acid sequence".format(input_seq))
