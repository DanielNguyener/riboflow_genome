#!/usr/bin/env python

from setuptools import setup, find_packages

# Read the version
exec(open('rfcommands/_version.py').read())

setup(
    name='rfcommands-riboflow',
    version=__version__,
    description='RiboFlow Commands with genome alignment support',
    long_description='Modified RiboFlow Commands (rfc) package with added genome alignment support, including HISAT2 log processing and enhanced statistics compilation with --label-prefix functionality.',
    author='RiboFlow Team',
    url='https://github.com/your-org/riboflow',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'rfc=rfcommands.cli.main:cli',
        ],
    },
    install_requires=[
        'click>=7.0',
        'pandas',
        'numpy'
    ],
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='bioinformatics ribosome profiling rna-seq nextflow',
)