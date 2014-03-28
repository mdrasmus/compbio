#!/usr/bin/env python
#
# setup for the compbio package
#
# use the following to install:
#   python setup.py install
#

from distutils.core import setup
import os

VERSION = '0.9'

scripts = [os.path.join('bin', x) for x in os.listdir('bin')]

setup(
    name='compbio',
    version=VERSION,
    description='Python libraries and utilities for computational biology',
    long_description = """
            """,
    author='Matt Rasmussen',
    author_email='matt.rasmus@gmail.edu',
    url='https://github.com/mdrasmus/compbio',
    download_url=
        'https://codeload.github.com/mdrasmus/compbio/tar.gz/v%s' % VERSION,

    packages=[
        'compbio',
        'compbio.synteny',
        'compbio.vis',
        'rasmus',
        'rasmus.ply',
        'rasmus.sexp',
        'rasmus.vis',
    ],
    scripts=scripts,
)
