#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Setup module for the RDKit Utilities module
#
# July 2026

from setuptools import setup
import os

# Pull in the essential run-time requirements
with open('requirements.txt') as file:
    requirements = file.read().splitlines()


# Use the README.rst as the long description.
def get_long_description():
    return open('README.rst').read()


setup(

    name='im-rdkit-utilities',
    version=os.environ.get('GITHUB_REF_SLUG', '1.0.0'),
    author='Informatics Matters',
    author_email='info@informaticsmatters.com',
    url='https://github.com/informaticsmatters/squonk2-rdkit-utilities',
    license='MIT',
    description='RDKit utilities for Squonk2 Data Manager Jobs',
    long_description=get_long_description(),
    keywords='rdkit',
    platforms=['any'],

    # Our modules to package
    packages=['rdkit_utils'],
    include_package_data=True,

    # Project classification:
    # https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Other Environment',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: 3.14',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Operating System :: POSIX :: Linux',
    ],

    python_requires='>=3.10',

    install_requires=requirements,

    zip_safe=False,

)
