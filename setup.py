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


"""
Based on:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
long_description = '''
ICE tool is a CRISPR editing analysis tool that infers presence of indels and other mutation. ICE uses non-negative least squares regression to detect the presence or evidence of edits. In contrast to TIDE1, ICE can analyze insertions, deletions, HDR, multiplex edits, and base editing and is available open source for non-commercial use (See LICENSE FOR DETAILS). ICE can also be used for analysis of other genome engineering methods, such as TALEN and homing endonucleases.

For more information, visit https://github.com/synthego-open/ice .

Algorithm described in https://www.biorxiv.org/content/early/2018/01/20/251082

'''

with open('requirements.txt') as fp:
    install_requires = fp.read()


about = {}
with open(path.join(here, 'ice', '__version__.py'), 'r', 'utf-8') as f:
    exec(f.read(), about)

setup(

    name='synthego_ice',
    version=about['__version__'],
    description='Synthego - Inference of CRISPR Edits (ICE)',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/synthego-open/ice',
    author='Synthego',  # Optional

    # For a list of valid classifiers, see
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish
        'License :: Free for non-commercial use',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    # This field adds keywords for your project which will appear on the
    # project page. What does your project relate to?
    # Note that this is a string of words separated by whitespace, not a list.
    keywords='sanger crispr indel genomics',  # Optional

    packages=find_packages(exclude=['doc']),  # Required

    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=install_requires,

    extras_require={
        'test': ['py.test','coverage'],
    },

    entry_points={  # Optional
        'console_scripts': [
            'synthego_ice=ice.analysis:single_sanger_analysis_cli',
            'synthego_ice_batch=ice.analysis:multiple_sanger_analysis_cli',
        ],
    },
)
