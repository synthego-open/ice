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

import argparse
import os
import traceback

import pandas as pd
import datetime
import json

from ice.classes.ice_result import ICEResult
from ice.classes.sanger_analysis import SangerAnalysis
from ice.utility.misc import version
from ice.utility.sequence import is_nuc_acid

from .__version__ import __version__


def single_sanger_analysis(control_path, sample_path, base_outputname, guide, donor=None, verbose=False,
                           allprops=False):
    """

    :param control_path: path to control ab1 file
    :param sample_path: path to sample ab1 file
    :param base_outputname: path to output (eg, /path/to/out/basename) will be saved using basename
    :param guide: RNA/DNA sequence of guide
    :param verbose: verbosity True/False
    :return: results dictionary
    """

    if control_path is None or not os.path.exists(control_path):
        raise Exception('Control @ {} not found'.format(control_path))

    if not os.path.exists(sample_path):
        raise Exception('Experiment sample @ {} not found'.format(sample_path))

    # create the base directory if this location doesn't exist
    base_dir = os.path.join(*os.path.split(os.path.abspath(base_outputname))[:-1])
    if verbose:
        print('Base dir: %s' % base_dir)

    if not os.path.exists(base_dir):
        os.makedirs(base_dir, exist_ok=True)

    sa = SangerAnalysis(verbose=verbose)
    sa.debug = False
    sa.allprops = allprops
    sa.initialize_with(control_path=control_path,
                       edited_path=sample_path,
                       gRNA_sequences=guide,
                       indel_max_size=20,
                       base_outputname=base_outputname,
                       donor=donor,
                       allprops=allprops)
    try:
        sa.analyze_sample()
        return sa.results.to_json(sa.guide_targets, sa.warnings)

    except Exception as e:
        results = ICEResult()

        print('Exception Caught!')
        traceback.print_exc()
        if isinstance(e,KeyError):
            return results.to_json(sa.guide_targets, [';'.join(sa.warnings)])
        else:
            return results.to_json(sa.guide_targets, [str(e)])


def single_sanger_analysis_cli():
    """ provides command line arg parsing for single sample analysis"""

    parser = argparse.ArgumentParser(description='Analyze Sanger reads to Infer Crispr Edit outcomes')
    parser.add_argument('--control', dest='control', help='The wildtype / unedited ab1 file (REQUIRED)', required=True)
    parser.add_argument('--edited', dest='edited', help='The edited ab1 file (REQUIRED)', required=True, default=None)
    parser.add_argument('--target', dest='target', help='Target sequence(s) (17-23 bases, RNA or DNA, comma separated), (REQUIRED)',
                        required=True)
    parser.add_argument('--out', dest='out', help='Output base path (Defaults to ./results/single)', required=False,
                        default=None)
    parser.add_argument('--donor', dest='donor', help='Donor DNA sequence for HDR (Optional)', required=False,
                        default=None)
    parser.add_argument('--verbose', dest='verbose', action='store_true')
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('--allprops', dest='allprops', action='store_true', default=False, help="Output all Edit Proposals, even if they have zero contribution")

    args = parser.parse_args()

    assert os.path.isfile(args.control)
    assert os.path.isfile(args.edited)

    if args.out is None:
        out_dir = os.path.join(os.path.abspath('.'), 'results', 'single')
    else:
        out_dir = os.path.abspath(args.out)

    base_dir = os.path.join(*os.path.split(os.path.abspath(out_dir))[:-1])
    if not os.path.exists(base_dir):
        os.makedirs(base_dir, exist_ok=True)

    print('Synthego ICE (https://synthego.com)')
    print('Version: {}'.format(__version__))

    single_sanger_analysis(control_path=args.control,
                           sample_path=args.edited,
                           base_outputname=out_dir,
                           guide=args.target,
                           donor=args.donor,
                           verbose=args.verbose,
                           allprops=args.allprops)


def XLSDictReader(f, sheet_index=0):
    book    = xlrd.open_workbook(file_contents=mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ))
    sheet   = book.sheet_by_index(sheet_index)
    headers = dict( (i, sheet.cell_value(0, i) ) for i in range(sheet.ncols) )

    return ( dict( (headers[j], sheet.cell_value(i, j)) for j in headers ) for i in range(1, sheet.nrows) )

def multiple_sanger_analysis(definition_file, output_dir,
                             data_dir=None,
                             verbose=False,
                             single_line=None,
                             allprops=False):
    '''
    :param definition_file: input excel file that defines sample/control/data associations
    :param output_dir: output directory
    :return:
    '''




    input_df = pd.read_excel(definition_file)

    input_df = input_df.rename(columns={"Donor Sequence": "Donor", "Control": "Control File", "Experiment": "Experiment File"})

    results = []

    fails = []

    jobs = []
    n = 0
    for m, experiment in input_df.iterrows():

        label = experiment['Label']

        base_outputname = os.path.join(output_dir, '%s-%s' % (n, label))

        control_sequence_file = experiment['Control File']

        edit_sequence_file = experiment['Experiment File']

        guide = experiment['Guide Sequence']

        if 'Donor' in experiment and is_nuc_acid(experiment['Donor']):
            donor = experiment['Donor']
        else:
            donor = None

        print(donor)
        try:
            if pd.isnull(control_sequence_file):
                raise IOError("Control filepath not specified at line {} in definition file".format(n+1))
            if pd.isnull(edit_sequence_file):
                raise IOError("Edit filepath not specified at line {} in definition file".format(n+1))

            control_sequence_path = os.path.join(data_dir, control_sequence_file)
            edit_sequence_path = os.path.join(data_dir, edit_sequence_file)

            if single_line is not None:
                if n != single_line:
                    continue

            msg = "analyzing"
            print("-" * 50, msg, n, experiment['Label'], guide)


            job_args = (control_sequence_path, edit_sequence_path, base_outputname, guide)
            job_kwargs = {
                'verbose': verbose,
                'allprops': allprops,
                'donor': donor
            }
            result = single_sanger_analysis(*job_args, **job_kwargs)
            jobs.append((experiment, result))

        except Exception as e:
            fails.append(experiment)
            print("Single Sanger analysis failed", e)
            import traceback, sys
            traceback.print_exc(file=sys.stdout)

        n += 1

    for job in jobs:
        r = job[1]
        experiment = job[0]
        if 'Donor' in experiment and is_nuc_acid(experiment['Donor']):
            donor = experiment['Donor']
        else:
            donor = None

        if r is not None:
            tmp = [experiment['Label'], r['ice'], r['ice_d'], r['rsq'], r['hdr_pct'],r['ko_score'], r['guides'],
                   r['notes'], experiment['Experiment File'], experiment['Control File'], donor]
        else:
            tmp = [experiment['Label'], 'Failed', '', '', '', '', '', '', '', '']
        results.append(tmp)

    if results:

        input_df = pd.DataFrame(results)
        timestamp = '{:%Y-%m-%d-%H%M%S}'.format(datetime.datetime.now())
        out_file = os.path.join(output_dir, "ice.results.{}.xlsx".format(timestamp))

        header = ["sample_name", "ice", 'ice_d', "r_squared", "hdr_pct","ko_score", "guides", "notes",
                  "experiment_file", "control_file", "donor"]
        input_df.columns = header
        # to json
        out_dict = []
        for r in results:
            row = {}
            for idx, c in enumerate(header):
                row[c] = r[idx]
            out_dict.append(row)
        with open(out_file.replace('.xlsx', '.json'), 'w') as f:
            json.dump(out_dict, f, ensure_ascii=False)

        with pd.ExcelWriter(out_file) as writer:
            input_df.to_excel(writer, sheet_name="Results")

            md = {'version': __version__}
            metadata = pd.DataFrame.from_dict([md])
            metadata.to_excel(writer, sheet_name='Metadata')

        writer.save()

        return out_dict
    else:
        print("None of the samples were able to be analyzed")
        return False


def multiple_sanger_analysis_cli():
    """ provides command line arg parsing for batch analysis"""

    parser = argparse.ArgumentParser(description='Analyze Sanger reads to infer crispr edit outcomes')
    parser.add_argument('--in', dest='input', help='Input definition file in Excel xlsx format (required)',
                        required=True)
    parser.add_argument('--out', dest='out', help='Output directory path (defaults to .)', required=False, default=None)
    parser.add_argument('--data', dest='data', help='Data path, where .ab1 files are located (required)', required=True,
                        default=None)
    parser.add_argument('--verbose', dest='verbose', action='store_true', help='Display verbose output')
    parser.add_argument('--line', dest='line', default=None, type=int,
                        help="Only run specified line in the Excel xlsx definition file")
    parser.add_argument('--allprops', dest='allprops', action='store_true', default=False, help="Output all Edit Proposals, even if they have zero contribution")
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))

    # parse args
    args = parser.parse_args()
    verbose = False
    if args.verbose:
        verbose = True
    if args.allprops:
        allprops = True
    else:
        allprops = False

    assert os.path.isfile(args.input)

    if args.out is None:
        out_dir = os.path.abspath(os.path.curdir)
    else:
        out_dir = os.path.abspath(args.out)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if args.data is None:
        data_dir = os.path.abspath(os.path.curdir)
    else:
        data_dir = os.path.abspath(args.data)

    print('Synthego ICE (https://synthego.com)')
    print('Version: {}'.format(__version__))

    multiple_sanger_analysis(args.input,
                             output_dir=out_dir,
                             data_dir=data_dir,
                             verbose=verbose,
                             single_line=args.line,
                             allprops=allprops)
