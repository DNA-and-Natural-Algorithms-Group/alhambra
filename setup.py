#!/usr/bin/env python
from setuptools import setup

setup(
    name="alhambra",
    version="1.1.0",
    packages=['alhambra'],
    install_requires=[
        'numpy', 'stickydesign >= 0.8.1', 'svgwrite', 'lxml', 'shutilwhich',
        'peppercompiler >= 0.1.2', 'ruamel.yaml >= 0.15.100', 'cssutils'
    ],
    include_package_data=True,
    entry_points={'console_scripts': ['alhambra = alhambra.scripts:alhambra']},
    author="Constantine Evans",
    author_email="cevans@dna.caltech.edu",
    description="DX Tile Set Designer", )
