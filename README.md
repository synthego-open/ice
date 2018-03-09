
# Synthego Inference of CRISPR Edits (ICE)

[![CircleCI Status](https://circleci.com/gh/synthego-open/ice.svg?style=shield&circle-token=:circle-token)](https://circleci.com/gh/synthego-open/ice)

![Synthego ICE](./doc/ice_logo.png)

ICE tool is a CRISPR editing analysis tool that infers presence of indels and other mutation.  ICE uses non-negative least squares regression to detect the presence or evidence of edits.  In contrast to TIDE<sup>[1](#ref1)</sup>, ICE can analyze insertions, deletions, HDR, multiplex edits, and base editing and is available under a open-source license for non-commercial use.  ICE can also be used for analysis of other genome engineering methods, such as TALEN and homing endonucleases.

### Audience

This document is intended for technical users who have prior experience with CRISPR editing analysis. For more detailed user documentation, please visit Synthego’s [Help Center](https://help.synthego.com/), where you can find our ICE User Guide and additional documentation.

### Citing ICE

A preprint **Inference of CRISPR Edits from Sanger Trace Data**<sup>[2](#ref2)</sup> provides an overview and empirical test of ICE on over 1,800 real world edits. We ask that you cite our paper if you use ICE in work that leads to publication.

We will continue to improve ICE, so please refer to the version number in your publication. The version number can be found by running ICE with the --version option.

### Inputs

1. Control sample Sanger ab1 file
2. Edited sample Sanger ab1 file
3. Sequence of protospacer or gRNA target

### Outputs:

Overall Editing efficiency, plots of distribution of edit types, plots of discordance (a calculation of signal agreement to control sequence), annotated Sanger traces of the region flanking the cut site, and JSON files containing the data for all of those plots.

**Trace JSON**
The Sanger sequence traces for the region around the cut site for the edited and control samples are shown.

**Discordance & indel distribution files**
The Discord json shows the agreement of both the control and edited sample to the called base. A discordance of 1 indicates that there is zero signal in the base called at that particular position, whereas a discord of 0 would mean that all non-reference bases have zero signal (all signal agrees with called base). The alignment window is used to align the control sample to the edited sample, while the inference window denotes the subsection of data used for NNLS regression.<br><br>The indel distribution json shows the distribution of indel identified, summarized by length.  Thus, two different -1 indels would be summarized to the same bin.

**Sequence Contributions**
<pre>
Relative contribution of each sequence (normalized)
--------------------------------------------------------
0.3006 	 -1[g1] 	 CCCAACACAACCAGTTGCAGGCGCC|-CATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.1996 	 0[g1] 	     CCCAACACAACCAGTTGCAGGCGCC|CCATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.1818 	 1[g1] 	     CCCAACACAACCAGTTGCAGGCGCC|nCCATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.1128 	 -2[g1] 	 CCCAACACAACCAGTTGCAGGCGC-|-CATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0541 	 2[g1] 	     CCCAACACAACCAGTTGCAGGCGCC|nnCCATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0317 	 -1[g1] 	 CCCAACACAACCAGTTGCAGGCGC-|CCATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0296 	 -4[g1] 	 CCCAACACAACCAGTTGCAGGC---|-CATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0185 	 -3[g1] 	 CCCAACACAACCAGTTGCAGGC---|CCATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0134 	 -19[g1] 	 CCCAACACAACCAGT----------|---------GCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0078 	 -16[g1] 	 CCCAACACAACCAGTTG--------|--------AGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0067 	 -18[g1] 	 CCCAACACAA---------------|---TGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0060 	 -20[g1] 	 CCCAACACAA---------------|-----GTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0053 	 -16[g1] 	 CCCAACACAACCAGTTGCAGGC---|-------------CAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0042 	 -16[g1] 	 CCCAACACAA---------------|-CATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0032 	 -15[g1] 	 CCCAACACAACC-------------|--ATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0028 	 -4[g1] 	 CCCAACACAACCAGTTGCAGGCGCC|----GGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0021 	 -20[g1] 	 CCCAACACAACCA------------|--------AGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0012 	 -17[g1] 	 CCCAACACAA---------------|--ATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
0.0005 	 -13[g1] 	 CCCAACACAACC-------------|CCATGGTGAGCATCAGCCTCTGGGTGGCCCTCCCTCTGGGCCTCGGGTATTTATGGAGCTGGATCCAAGGTCACATGCTTGTTCATGAGCTCTCAGGCA
</pre>

The first column indicates the proportion of that sequence inferred in the pool. The second column is a summarized identity indicating the size of the indel and which guide. The third column is a human-readable representation of the sequence with dashes indicating deletions and 'n' indicating insertions.

Additional files such as alignment verification are generated for each sample.

## Using ICE

A hosted free version of the ICE tool is available online at [https://ice.synthego.com](https://ice.synthego.com). The online ICE tool supports batch analysis, figure generation, and error handling & sample QC.

The source code behind the core ICE analysis is open source and free to use for non-commercial applications. Commercial use and other licensing options are available. For details, see [LICENSE](./LICENSE).

## Installation

Synthego ICE can be installed as a docker container or directly via pip. Additional developer instructions are located in [DEVELOP.md](./DEVELOP.md). All examples below use test data found in [ice/tests/test_data](./ice/tests/test_data).  The test file (./ice/tests/test_data/batch_example.xlsx) is an example of how to specify batch inputs.

### Method 1. Pip install

Install into your favorite python3 virtual environment (virtualenv, conda).

```bash
conda create --name ice_env python=3 # create a python3 virtual environment

source activate ice_env # activate the virtual environment

pip install sythego_ice # install synthego ice from pip
```

After installation, you can use Synthego ICE as a module (see  [python_example.py](./python_example.py)) or directly via command line.

#### Command line tools

`synthego_ice`

    usage: synthego_ice [-h] --control CONTROL --edited EDITED --target TARGET
                        [--out OUT] [--donor DONOR] [--verbose] [--version]

    Analyze Sanger reads to Infer Crispr Edit outcomes

    optional arguments:
      -h, --help         show this help message and exit
      --control CONTROL  The wildtype / unedited ab1 file (REQUIRED)
      --edited EDITED    The edited ab1 file (REQUIRED)
      --target TARGET    Target sequence(s) (17-23 bases, RNA or DNA, comma
                         separated), (REQUIRED)
      --out OUT          Output base path (Defaults to ./results/single)
      --donor DONOR      Donor DNA sequence for HDR (Optional)
      --verbose
      --version          show program's version number and exit

`synthego_ice_batch`

    usage: synthego_ice_batch [-h] --in INPUT [--out OUT] --data DATA [--verbose]
                              [--line LINE] [--allprops] [--version]

    Analyze Sanger reads to infer crispr edit outcomes

    optional arguments:
      -h, --help   show this help message and exit
      --in INPUT   Input definition file in Excel xlsx format (required)
      --out OUT    Output directory path (defaults to .)
      --data DATA  Data path, where .ab1 files are located (required)
      --verbose    Display verbose output
      --line LINE  Only run specified line in the Excel xlsx definition file
      --allprops   Output all Edit Proposals, even if they have zero contribution
      --version    show program's version number and exit


#### Analyzing example data via command line tools

After installing via pip, grab the example data by cloning this repository:

```bash
git clone git@github.com:synthego-open/ice.git ice

cd ice # change into the ice directory
```

**Analyzing a single sample**

```bash
synthego_ice \
	--control ./ice/tests/test_data/good_example_control.ab1 \
	--edited ./ice/tests/test_data/good_example_edited.ab1 \
	--target AACCAGTTGCAGGCGCCCCA \
	--out results/testing \
	--verbose
```

When complete, you'll find the ICE analysis outputs in the `./results` folder.

**Analyzing multiple samples (batch analysis)**

```bash
synthego_ice_batch \
	--in ./ice/tests/test_data/batch_example.xlsx \
	--out ./results/ \
	--data ./ice/tests/test_data/
	--verbose
```
When complete, you'll find the ICE analysis outputs in the `./results` folder.


### Method 2. Installing via Docker container

#### Requirements:

* Docker http://docker.com

#### Installation

From a command line, grab the latest version of Synthego ICE from Docker Hub.

```bash
docker pull synthego/ice
```

After installation, you'll be able to run Synthego ICE from the docker container.

#### Analyzing example data with docker

Grab the example data by cloning this repository:

```bash
git clone git@github.com:synthego-open/ice.git ice

cd ice # change into the ice directory
```

**Analyzing a single sample**

```bash
docker run -it -v ${PWD}:/data -w /ice -i ice:latest \
	python ice_analysis_single.py \
	--control /data/ice/tests/test_data/good_example_control.ab1 \
	--edited /data/ice/tests/test_data/good_example_edited.ab1 \
	--target AACCAGTTGCAGGCGCCCCA \
	--out /data/results/testing \
	--verbose
```

When complete, you'll find the ICE analysis outputs in the `./results` folder.

**Analyzing multiple samples (batch analysis)**

```bash
docker run -it -v ${PWD}:/data -w /ice -i ice:latest \
	python ice_analysis_batch.py \
	--in /data/ice/tests/test_data/batch_example.xlsx \
	--data /data/ice/tests/test_data/ \
	--out /data/results/ \
	--verbose
```

When complete, you'll find the ICE analysis outputs in the `./results` folder.

#### References

<a name="ref1">[1] </a>Eva K. Brinkman, Tao Chen, Mario Amendola, and Bas van Steensel. <b>Easy quantitative assessment of genome editing by sequence trace decomposition.</b>Nucleic Acids Res. 2014 Dec 16; 42(22): e168. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267669/

<a name="ref2">[2] </a>Hsiau et. al, <b>Inference of CRISPR Edits from Sanger Trace Data</b>. BioArxiv. 2018 https://www.biorxiv.org/content/early/2018/01/20/251082

### License

Copyright 2018 [Synthego Corporation](http://synthego.com) All Rights Reserved

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
