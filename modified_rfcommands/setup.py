#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import os

# Read version from _version.py
version_file = os.path.join(os.path.dirname(__file__), 'rfcommands', '_version.py')
version = '0.0.0'
if os.path.exists(version_file):
    with open(version_file) as f:
        version_dict = {}
        exec(f.read(), version_dict)
        version = version_dict.get('__version__', '0.0.0')

# Read long description from README if it exists
long_description = ""
readme_file = os.path.join(os.path.dirname(__file__), 'README.md')
if os.path.exists(readme_file):
    with open(readme_file, 'r', encoding='utf-8') as f:
        long_description = f.read()

setup(
    name='rfcommands-riboflow',
    version=version,
    description='Modified RiboFlow Commands (rfc) with genome alignment support',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='RiboFlow Genome',
    url='https://github.com/ribosomeprofiling/riboflow',
    packages=find_packages(),
    python_requires='>=3.7',
    install_requires=[
        'click',
        'pandas',
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'rfc=rfcommands.cli.main:cli',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)

